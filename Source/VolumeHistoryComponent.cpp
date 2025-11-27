#include "VolumeHistoryComponent.h"
#include "PluginProcessor.h"

#include <cmath>
#include <algorithm>

//==============================================================================

VolumeHistoryComponent::VolumeHistoryComponent (LevelScopeAudioProcessor& proc)
    : processor (proc),
      visualFrameRate (proc.getLoudnessFrameRate()),
      historyLengthSeconds (1800.0),      // 30 minutes = 1800 s
      minDb (-90.0f),
      maxDb (  0.0f),
      baseDbRange (std::abs (minDb)),     // 90 dB
      historyCapacityFrames ((int) std::ceil (historyLengthSeconds * visualFrameRate)),
      historyDbMomentary  ((size_t) historyCapacityFrames, minDb),
      historyDbShortTerm  ((size_t) historyCapacityFrames, minDb)
{
    jassert (visualFrameRate > 0.0);
    jassert (historyCapacityFrames > 0);

    setOpaque (true);
    startTimerHz (60); // UI update rate
}

VolumeHistoryComponent::~VolumeHistoryComponent()
{
    stopTimer();
}

//==============================================================================

void VolumeHistoryComponent::resized()
{
    // Initialise X zoom so that, roughly, 10 seconds are visible by default.
    // visibleSeconds ≈ (width / zoomX) / visualFrameRate
    // => zoomX ≈ width / (visibleSeconds * visualFrameRate)
    if (! hasCustomZoomX)
    {
        const int w = getWidth();
        if (w > 0)
        {
            const double desiredVisibleSeconds = 10.0;
            zoomX = (double) w / (desiredVisibleSeconds * visualFrameRate);
            zoomX = juce::jlimit (minZoomX, maxZoomX, zoomX);
        }
    }
}

//==============================================================================

void VolumeHistoryComponent::timerCallback()
{
    drainProcessorFifo();
    repaint();
}

void VolumeHistoryComponent::drainProcessorFifo()
{
    // Drain all available loudness frames from the processor
    constexpr int chunkSize = 512;
    float momentaryValues [chunkSize];
    float shortTermValues [chunkSize];

    for (;;)
    {
        const int numRead = processor.readLoudnessFromFifo (momentaryValues,
                                                            shortTermValues,
                                                            chunkSize);
        if (numRead <= 0)
            break;

        pushLoudnessBatchToHistory (momentaryValues, shortTermValues, numRead);
    }
}

void VolumeHistoryComponent::pushLoudnessBatchToHistory (const float* momentaryValues,
                                                         const float* shortTermValues,
                                                         int numValues)
{
    if (numValues <= 0)
        return;

    for (int i = 0; i < numValues; ++i)
    {
        const float rmsM = momentaryValues[i];
        const float rmsS = shortTermValues[i];

        // Convert RMS (linear) to dB, clamped to minDb
        const float dbM  = juce::Decibels::gainToDecibels (rmsM, minDb);
        const float dbS  = juce::Decibels::gainToDecibels (rmsS, minDb);

        const juce::int64 writeIndex = historyFrameCount;
        const int ringIndex = (int) (writeIndex % (juce::int64) historyCapacityFrames);

        historyDbMomentary[(size_t) ringIndex] = dbM;
        historyDbShortTerm[(size_t) ringIndex] = dbS;

        ++historyFrameCount;
    }

    // viewOffsetFrames stays the same; "now" has moved,
    // but our offset behind "now" remains constant.
}

//==============================================================================

float VolumeHistoryComponent::getDbAtFrame (CurveKind kind,
                                            juce::int64 frameIndex) const noexcept
{
    if (historyFrameCount <= 0)
        return minDb;

    const juce::int64 latestIndex = historyFrameCount - 1;
    const juce::int64 earliestIndex =
        (historyFrameCount > (juce::int64) historyCapacityFrames
            ? historyFrameCount - (juce::int64) historyCapacityFrames
            : 0);

    if (frameIndex < earliestIndex || frameIndex > latestIndex)
        return minDb;

    const int ringIndex = (int) (frameIndex % (juce::int64) historyCapacityFrames);

    switch (kind)
    {
        case CurveKind::Momentary: return historyDbMomentary[(size_t) ringIndex];
        case CurveKind::ShortTerm: return historyDbShortTerm[(size_t) ringIndex];
    }

    return minDb;
}

float VolumeHistoryComponent::dbToY (float db, float height) const noexcept
{
    if (height <= 0.0f)
        return 0.0f;

    // Visible dB range shrinks as yZoom increases (zoom in)
    const float effectiveRange = (float) (baseDbRange / yZoom);
    const float topDb          = maxDb;
    const float bottomDb       = topDb - effectiveRange;

    const float clamped = juce::jlimit (bottomDb, topDb, db);

    const float norm = (clamped - bottomDb) / effectiveRange; // 0 at bottom, 1 at top
    const float y    = height * (1.0f - norm);

    return y;
}

//==============================================================================

void VolumeHistoryComponent::paint (juce::Graphics& g)
{
    auto bounds = getLocalBounds().toFloat();

    g.fillAll (juce::Colours::black);

    const int   width  = (int) bounds.getWidth();
    const float height = bounds.getHeight();

    if (width <= 1 || height <= 0.0f)
        return;

    // Draw auto-fit toggle bar on the left
    {
        juce::Rectangle<float> bar (0.0f, 0.0f, (float) autoFitBarWidth, height);
        g.setColour (autoFitEnabled ? juce::Colours::red.withAlpha (0.7f)
                                    : juce::Colours::blue.withAlpha (0.4f));
        g.fillRect (bar);
    }

    // Optional horizontal reference lines (e.g. -90, -60, -30, 0 dB)
    g.setColour (juce::Colours::darkgrey.withMultipliedAlpha (0.6f));
    const int numLines = 4;
    for (int i = 0; i <= numLines; ++i)
    {
        const float db = maxDb - (baseDbRange / (float) numLines) * (float) i;
        const float y  = dbToY (db, height);
        g.drawHorizontalLine ((int) std::round (y), 0.0f, bounds.getWidth());
    }

    if (historyFrameCount < 2)
        return;

    const juce::int64 latestIndex = historyFrameCount - 1;
    const juce::int64 earliestIndex =
        (historyFrameCount > (juce::int64) historyCapacityFrames
            ? historyFrameCount - (juce::int64) historyCapacityFrames
            : 0);

    const double maxOffsetFrames = (double) (latestIndex - earliestIndex);

    if (autoFitEnabled && width > autoFitBarWidth)
    {
        // Auto-fit: keep the entire history visible from earliest to latest
        const double widthPixels     = (double) (width - autoFitBarWidth);
        const double framesInHistory = (double) (latestIndex - earliestIndex + 1);

        if (framesInHistory > 1.0 && widthPixels > 1.0)
        {
            zoomX = widthPixels / framesInHistory;
            zoomX = juce::jmin (zoomX, maxZoomX);  // clamp only to max; allow very small zoomX

            viewOffsetFrames = 0.0;   // newest frame at the right edge
            hasCustomZoomX   = true;
        }
    }
    else
    {
        // Manual mode: clamp existing offset/zoom as before
        viewOffsetFrames = juce::jlimit (0.0, maxOffsetFrames, viewOffsetFrames);
        zoomX            = juce::jlimit (minZoomX, maxZoomX, zoomX);
    }

    // Right edge: which frame index does it show?
    const double rightIndexDouble = (double) latestIndex - viewOffsetFrames;

    // Start at the integer frame at/just before the right edge
    juce::int64 currentIndex = (juce::int64) std::floor (rightIndexDouble);

    if (currentIndex < earliestIndex || currentIndex > latestIndex)
        return;

    const int rightPixel = width - 1;

    auto drawCurve = [&] (CurveKind kind, juce::Colour colour, float strokeWidth)
    {
        juce::Path path;
        bool started = false;

        g.setColour (colour);

        int frameStep = 0;
        juce::int64 idx = currentIndex;

        for (;;)
        {
            if (idx < earliestIndex)
                break;

            const double xPos = (double) rightPixel - (double) frameStep * zoomX;

            // Optimization: stop if we are far off the left edge
            if (xPos < -50.0)
                break;

            const float db = getDbAtFrame (kind, idx);
            const float y  = dbToY (db, height);

            if (! started)
            {
                path.startNewSubPath ((float) xPos + 0.5f, y);
                started = true;
            }
            else
            {
                path.lineTo ((float) xPos + 0.5f, y);
            }

            --idx;
            ++frameStep;
        }

        if (started)
            g.strokePath (path, juce::PathStrokeType (strokeWidth));
    };

    // Draw short-term first (smoother curve), then momentary on top
    drawCurve (CurveKind::ShortTerm, juce::Colours::cyan.withMultipliedAlpha (0.8f), 2.4f);
    drawCurve (CurveKind::Momentary, juce::Colours::limegreen,                    3.0f);

    // Overlay info: curves + zooms
    g.setColour (juce::Colours::white);
    g.setFont (14.0f);

    juce::String info1 = "Momentary (400 ms)  |  Short-term (3 s)";
    juce::String info2 = "ZoomX: " + juce::String (zoomX, 2)
                       + "  ZoomY: " + juce::String (yZoom, 2);

    g.drawText (info1,  8,  8, 260, 20, juce::Justification::topLeft);
    g.drawText (info2,  8, 28, 260, 20, juce::Justification::topLeft);
}

//==============================================================================

void VolumeHistoryComponent::applyHorizontalZoom (float wheelDelta, float mouseX)
{
    if (historyFrameCount <= 1 || getWidth() <= 1 || wheelDelta == 0.0f)
        return;

    const int w = getWidth();
    const double widthPixels = (double) (w - 1);

    const double clampedMouseX = juce::jlimit (0.0, widthPixels, (double) mouseX);
    const double dxFromRight   = widthPixels - clampedMouseX; // pixels from right edge

    const juce::int64 latestIndex = historyFrameCount - 1;
    const juce::int64 earliestIndex =
        (historyFrameCount > (juce::int64) historyCapacityFrames
            ? historyFrameCount - (juce::int64) historyCapacityFrames
            : 0);

    const double maxOffsetFrames = (double) (latestIndex - earliestIndex);

    // Ensure current offset is valid
    viewOffsetFrames = juce::jlimit (0.0, maxOffsetFrames, viewOffsetFrames);

    // Current right-edge frame index
    const double rightIndex = (double) latestIndex - viewOffsetFrames;

    // How many frames from right edge to mouse (in frame steps)?
    const double k = dxFromRight / zoomX;

    // Frame index currently under the mouse cursor
    double frameAtMouse = rightIndex - k;
    frameAtMouse = juce::jlimit ((double) earliestIndex, (double) latestIndex, frameAtMouse);

    // Zoom factor: positive deltaY = zoom in (larger zoomX => more pixels per frame)
    const double zoomBase   = 1.1;
    const double zoomFactor = std::pow (zoomBase, (double) wheelDelta);

    zoomX *= zoomFactor;
    zoomX = juce::jlimit (minZoomX, maxZoomX, zoomX);

    // Recompute how many frame steps from right edge to that same frame under new zoom
    const double kNew          = dxFromRight / zoomX;
    const double rightIndexNew = frameAtMouse + kNew;

    double newOffset = (double) latestIndex - rightIndexNew;
    newOffset = juce::jlimit (0.0, maxOffsetFrames, newOffset);

    viewOffsetFrames = newOffset;
    hasCustomZoomX   = true;
}

void VolumeHistoryComponent::applyVerticalZoom (float wheelDelta)
{
    if (wheelDelta == 0.0f)
        return;

    // Positive deltaY = zoom in (smaller dB range)
    const double zoomBase   = 1.1;
    const double zoomFactor = std::pow (zoomBase, (double) wheelDelta);

    yZoom *= zoomFactor;
    yZoom = juce::jlimit (minYZoom, maxYZoom, yZoom);
}

//==============================================================================

void VolumeHistoryComponent::mouseWheelMove (const juce::MouseEvent& event,
                                             const juce::MouseWheelDetails& wheel)
    if (autoFitEnabled)
    {
        autoFitEnabled = false;
        repaint();
        return;
    }

void VolumeHistoryComponent::mouseDown (const juce::MouseEvent& event)
{
    // Toggle auto-fit when clicking inside the left bar
    if (event.position.x <= (float) autoFitBarWidth)
    {
        autoFitEnabled = ! autoFitEnabled;
        repaint();
        return;
    }

    // Otherwise, default behavior (if any) can go here
}
#include "VolumeHistoryComponent.h"
#include "PluginProcessor.h"

#include <cmath>
#include <algorithm>

// Helper: format seconds as HH:MM:SS
static juce::String formatTimeHMS (double seconds)
{
    if (seconds < 0.0)
        seconds = 0.0;

    const int total = (int) std::floor (seconds + 0.5);
    const int h = total / 3600;
    const int m = (total % 3600) / 60;
    const int s = total % 60;

    return juce::String::formatted ("%02d:%02d:%02d", h, m, s);
}

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
                                    : juce::Colours::blue.withAlpha (0.6f));
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

    // Auto-fit: force zoomX so entire history fits; manual: clamp zoomX.
    if (autoFitEnabled && width > autoFitBarWidth)
    {
        const double widthPixels     = (double) (width - autoFitBarWidth);
        const double framesInHistory = (double) (latestIndex - earliestIndex + 1);

        if (framesInHistory > 1.0 && widthPixels > 1.0)
        {
            zoomX = widthPixels / framesInHistory;
            zoomX = juce::jmin (zoomX, maxZoomX);  // clamp only to max; allow very small zoomX

            hasCustomZoomX   = true;
            minAllowedZoomX  = juce::jmax (minZoomX, zoomX);  // auto-fit defines minimum zoom-out
        }
    }
    else
    {
        zoomX = juce::jlimit (minAllowedZoomX, maxZoomX, zoomX);
    }

    // Right edge always follows "now"
    const double rightIndexDouble = (double) latestIndex;
    const juce::int64 currentIndex0 = (juce::int64) std::floor (rightIndexDouble);

    if (currentIndex0 < earliestIndex || currentIndex0 > latestIndex)
        return;

    const int rightPixel = width - 1;

    // ----------------------------------------------------------------------
    // Aggregation (decimation) only in auto-fit to avoid lag for long history.
    // In manual mode, frameStepSize = 1 (full resolution).
    // ----------------------------------------------------------------------
    int frameStepSize = 1;

    if (autoFitEnabled)
    {
        const double framesInHistory = (double) (latestIndex - earliestIndex + 1);
        const int    maxPoints       = juce::jmax (1, width - autoFitBarWidth);

        if (framesInHistory > (double) maxPoints)
            frameStepSize = (int) std::ceil (framesInHistory / (double) maxPoints);
    }

    auto drawCurve = [&] (CurveKind kind, juce::Colour colour, float strokeWidth)
    {
        juce::Path path;
        bool started = false;

        g.setColour (colour);

        int frameStep = 0;
        juce::int64 idx = currentIndex0;

        for (;;)
        {
            if (idx < earliestIndex)
                break;

            const double xPos = (double) rightPixel - (double) frameStep * zoomX;

            // Optimization: stop if we are far off the left edge
            if (xPos < -50.0)
                break;

            float dbGroup = getDbAtFrame (kind, idx);

            if (frameStepSize > 1)
            {
                float maxDb = dbGroup;

                const juce::int64 groupEnd =
                    juce::jmax<juce::int64> (earliestIndex, idx - (juce::int64) (frameStepSize - 1));

                for (juce::int64 j = idx - 1; j >= groupEnd; --j)
                {
                    const float v = getDbAtFrame (kind, j);
                    if (v > maxDb)
                        maxDb = v;
                }

                dbGroup = maxDb;
            }

            const float y = dbToY (dbGroup, height);

            if (! started)
            {
                path.startNewSubPath ((float) xPos + 0.5f, y);
                started = true;
            }
            else
            {
                path.lineTo ((float) xPos + 0.5f, y);
            }

            idx       -= (juce::int64) frameStepSize;
            frameStep += frameStepSize;
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

    //==========================================================================
    // Time ruler at bottom: HH:MM:SS, auto-adjusted spacing
    //==========================================================================

    const float rulerHeight   = 16.0f;
    const float rulerBaseY    = height - 2.0f;
    const float tickTopY      = rulerBaseY - 6.0f;
    const float textTopY      = rulerBaseY - 14.0f;
    const float leftPixel     = (float) autoFitBarWidth;

    const double widthPixelsForTime = (double) (width - autoFitBarWidth);

    if (widthPixelsForTime > 1.0 && zoomX > 0.0 && visualFrameRate > 0.0)
    {
        double visibleFramesApprox = 0.0;
        double rightIndexForTime   = (double) latestIndex;

        if (autoFitEnabled)
        {
            visibleFramesApprox = (double) (latestIndex - earliestIndex + 1);
        }
        else
        {
            visibleFramesApprox = widthPixelsForTime / zoomX;
        }

        if (visibleFramesApprox > 1.0)
        {
            const double visibleSeconds = visibleFramesApprox / visualFrameRate;
            double tRight = rightIndexForTime / visualFrameRate;
            double tLeft  = tRight - visibleSeconds;
            if (tLeft < 0.0) tLeft = 0.0;

            // Choose tick spacing so that we have at most ~20 labels.
            const double maxLabels = 20.0;
            const double tickSteps[] = {
                0.5, 1.0, 2.0, 5.0,
                10.0, 15.0, 30.0,
                60.0, 120.0, 300.0, 600.0, 1200.0, 3600.0
            };

            double tickStep = tickSteps[sizeof(tickSteps)/sizeof(tickSteps[0]) - 1];

            for (double s : tickSteps)
            {
                const double count = visibleSeconds / s;
                if (count <= maxLabels)
                {
                    tickStep = s;
                    break;
                }
            }

            const double firstTick = std::ceil (tLeft / tickStep) * tickStep;

            // Draw baseline for the ruler
            g.setColour (juce::Colours::darkgrey.withMultipliedAlpha (0.8f));
            g.drawHorizontalLine ((int) std::round (rulerBaseY), 0.0f, (float) width);

            // Draw ticks and labels
            g.setColour (juce::Colours::white);
            g.setFont (10.0f);

            for (double t = firstTick; t <= tRight + 1e-9; t += tickStep)
            {
                const double norm = (t - tLeft) / visibleSeconds;
                if (norm < 0.0 || norm > 1.0)
                    continue;

                const float x = leftPixel + (float) (norm * widthPixelsForTime);

                g.drawVerticalLine ((int) std::round (x), tickTopY, rulerBaseY);

                auto label = formatTimeHMS (t);
                const float textWidth = 52.0f;
                g.drawText (label,
                            x - textWidth * 0.5f,
                            textTopY,
                            textWidth,
                            rulerHeight,
                            juce::Justification::centred);
            }
        }
    }
}

//==============================================================================

void VolumeHistoryComponent::applyHorizontalZoom (float wheelDelta)
{
    if (historyFrameCount <= 1 || getWidth() <= 1 || wheelDelta == 0.0f)
        return;

    // Zoom factor: positive deltaY = zoom in (larger zoomX => more pixels per frame)
    const double zoomBase   = 1.1;
    const double zoomFactor = std::pow (zoomBase, (double) wheelDelta);

    zoomX *= zoomFactor;
    zoomX = juce::jlimit (minAllowedZoomX, maxZoomX, zoomX);

    hasCustomZoomX = true;
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
{
    if (wheel.deltaY == 0.0f)
        return;

    if (autoFitEnabled)
    {
        // In auto-fit:
        // - wheel down (deltaY < 0): ignore (already fully zoomed out)
        // - wheel up (deltaY > 0): exit auto-fit, keep following, and zoom in from current view
        if (wheel.deltaY < 0.0f)
            return;

        // Exit auto-fit. We do NOT restore old manual zoom here.
        autoFitEnabled = false;
    }

    // Manual mode (always following "now")
    if (event.mods.isShiftDown())
        applyVerticalZoom (wheel.deltaY);
    else
        applyHorizontalZoom (wheel.deltaY);

    repaint();
}

void VolumeHistoryComponent::mouseDown (const juce::MouseEvent& event)
{
    // Toggle auto-fit when clicking inside the left bar
    if (event.position.x <= (float) autoFitBarWidth)
    {
        if (! autoFitEnabled)
        {
            // Enter auto-fit: save current manual zoom
            savedManualZoomX    = zoomX;
            hasSavedManualZoomX = true;

            autoFitEnabled = true;
        }
        else
        {
            // Leave auto-fit: restore last manual zoom if available
            autoFitEnabled = false;

            if (hasSavedManualZoomX)
                zoomX = juce::jlimit (minAllowedZoomX, maxZoomX, savedManualZoomX);
        }

        repaint();
        return;
    }

    // Otherwise, default behavior (if any) can go here
}
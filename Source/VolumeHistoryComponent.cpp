#include "VolumeHistoryComponent.h"
#include "PluginProcessor.h"

#include <cmath>
#include <algorithm>

//==============================================================================

VolumeHistoryComponent::VolumeHistoryComponent (LevelScopeAudioProcessor& proc)
    : processor (proc),
      visualSampleRate (proc.getVisualSampleRate()),
      historyLengthSeconds (180.0),       // 180 s total history
      minDb (-90.0f),
      maxDb (  0.0f),
      baseDbRange (std::abs (minDb)),     // 90 dB
      historyCapacitySamples ((int) std::ceil (historyLengthSeconds * visualSampleRate)),
      historyDbRms  ((size_t) historyCapacitySamples, minDb),
      historyDbPeak ((size_t) historyCapacitySamples, minDb)
{
    jassert (visualSampleRate > 0.0);
    jassert (historyCapacitySamples > 0);

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
    // visibleSeconds ≈ (width / zoomX) / visualSampleRate
    // => zoomX ≈ width / (visibleSeconds * visualSampleRate)
    if (! hasCustomZoomX)
    {
        const int w = getWidth();
        if (w > 0)
        {
            const double desiredVisibleSeconds = 10.0;
            zoomX = (double) w / (desiredVisibleSeconds * visualSampleRate);
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
    // Drain all available envelope values (RMS+Peak) from the processor
    constexpr int chunkSize = 2048;
    float rmsValues [chunkSize];
    float peakValues[chunkSize];

    for (;;)
    {
        const int numRead = processor.readEnvelopeFromFifo (rmsValues, peakValues, chunkSize);
        if (numRead <= 0)
            break;

        pushEnvelopeBatchToHistory (rmsValues, peakValues, numRead);
    }
}

void VolumeHistoryComponent::pushEnvelopeBatchToHistory (const float* rmsValues,
                                                         const float* peakValues,
                                                         int numValues)
{
    if (numValues <= 0)
        return;

    for (int i = 0; i < numValues; ++i)
    {
        const float rms  = rmsValues[i];
        const float peak = peakValues[i];

        // Convert to dB, clamped to minDb
        const float dbRms  = juce::Decibels::gainToDecibels (rms,  minDb);
        const float dbPeak = juce::Decibels::gainToDecibels (peak, minDb);

        const juce::int64 writeIndex = historySampleCount;
        const int ringIndex = (int) (writeIndex % (juce::int64) historyCapacitySamples);

        historyDbRms [(size_t) ringIndex] = dbRms;
        historyDbPeak[(size_t) ringIndex] = dbPeak;

        ++historySampleCount;
    }

    // viewOffsetSamples stays the same; "now" has moved,
    // but our offset behind "now" remains constant.
}

//==============================================================================

float VolumeHistoryComponent::sampleIndexToDbNearest (juce::int64 sampleIndex) const noexcept
{
    if (historySampleCount <= 0)
        return minDb;

    const juce::int64 latestIndex = historySampleCount - 1;
    const juce::int64 earliestIndex =
        (historySampleCount > (juce::int64) historyCapacitySamples
            ? historySampleCount - (juce::int64) historyCapacitySamples
            : 0);

    if (sampleIndex < earliestIndex || sampleIndex > latestIndex)
        return minDb;

    const int ringIndex = (int) (sampleIndex % (juce::int64) historyCapacitySamples);

    if (envelopeMode == EnvelopeMode::RMS)
        return historyDbRms[(size_t) ringIndex];

    return historyDbPeak[(size_t) ringIndex];
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

    const int width    = (int) bounds.getWidth();
    const float height = bounds.getHeight();

    if (width <= 1 || height <= 0.0f)
        return;

    // Optional horizontal reference lines (e.g. -90, -60, -30, 0 dB)
    g.setColour (juce::Colours::darkgrey.withMultipliedAlpha (0.6f));
    const int numLines = 4;
    for (int i = 0; i <= numLines; ++i)
    {
        const float db = maxDb - (baseDbRange / (float) numLines) * (float) i;
        const float y  = dbToY (db, height);
        g.drawHorizontalLine ((int) std::round (y), 0.0f, bounds.getWidth());
    }

    if (historySampleCount < 2)
        return;

    // Clamp zoom and offset to valid ranges (defensive)
    const juce::int64 latestIndex = historySampleCount - 1;
    const juce::int64 earliestIndex =
        (historySampleCount > (juce::int64) historyCapacitySamples
            ? historySampleCount - (juce::int64) historyCapacitySamples
            : 0);

    const double maxOffsetSamples = (double) (latestIndex - earliestIndex);
    viewOffsetSamples = juce::jlimit (0.0, maxOffsetSamples, viewOffsetSamples);
    zoomX             = juce::jlimit (minZoomX, maxZoomX, zoomX);

    // Right edge: which sample index does it show?
    const double rightIndexDouble = (double) latestIndex - viewOffsetSamples;

    // Start at the integer sample at/just before the right edge
    juce::int64 currentIndex = (juce::int64) std::floor (rightIndexDouble);

    if (currentIndex < earliestIndex || currentIndex > latestIndex)
        return;

    const int rightPixel = width - 1;

    // Build path: newest (or offset) at the right, older samples to the left
    juce::Path path;
    bool started = false;

    // Different colours for RMS vs Peak
    const auto curveColour = (envelopeMode == EnvelopeMode::RMS
                                ? juce::Colours::limegreen
                                : juce::Colours::orange);
    g.setColour (curveColour);

    int sampleStep = 0;
    for (;;)
    {
        if (currentIndex < earliestIndex)
            break;

        const double xPos = (double) rightPixel - (double) sampleStep * zoomX;

        // Optimization: stop if we are far off the left edge
        if (xPos < -50.0)
            break;

        const float db = sampleIndexToDbNearest (currentIndex);
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

        --currentIndex;
        ++sampleStep;
    }

    if (started)
        g.strokePath (path, juce::PathStrokeType (1.5f));

    // Overlay info: mode + zooms
    g.setColour (juce::Colours::white);
    g.setFont (14.0f);

    juce::String modeText = (envelopeMode == EnvelopeMode::RMS ? "Mode: RMS" : "Mode: Peak");
    juce::String zoomInfo = "ZoomX: " + juce::String (zoomX, 2)
                          + "  ZoomY: " + juce::String (yZoom, 2);

    g.drawText (modeText,  8,  8, 160, 20, juce::Justification::topLeft);
    g.drawText (zoomInfo,  8, 28, 220, 20, juce::Justification::topLeft);
}

//==============================================================================

void VolumeHistoryComponent::applyHorizontalZoom (float wheelDelta, float mouseX)
{
    if (historySampleCount <= 1 || getWidth() <= 1 || wheelDelta == 0.0f)
        return;

    const int w = getWidth();
    const double widthPixels = (double) (w - 1);

    const double clampedMouseX = juce::jlimit (0.0, widthPixels, (double) mouseX);
    const double dxFromRight   = widthPixels - clampedMouseX; // pixels from right edge

    const juce::int64 latestIndex = historySampleCount - 1;
    const juce::int64 earliestIndex =
        (historySampleCount > (juce::int64) historyCapacitySamples
            ? historySampleCount - (juce::int64) historyCapacitySamples
            : 0);

    const double maxOffsetSamples = (double) (latestIndex - earliestIndex);

    // Ensure current offset is valid
    viewOffsetSamples = juce::jlimit (0.0, maxOffsetSamples, viewOffsetSamples);

    // Current right-edge sample index
    const double rightIndex = (double) latestIndex - viewOffsetSamples;

    // How many samples from right edge to mouse (in sample steps)?
    const double k = dxFromRight / zoomX;

    // Sample index currently under the mouse cursor
    double sampleAtMouse = rightIndex - k;
    sampleAtMouse = juce::jlimit ((double) earliestIndex, (double) latestIndex, sampleAtMouse);

    // Zoom factor: positive deltaY = zoom in (larger zoomX => more pixels per sample)
    const double zoomBase   = 1.1;
    const double zoomFactor = std::pow (zoomBase, (double) wheelDelta);

    zoomX *= zoomFactor;
    zoomX = juce::jlimit (minZoomX, maxZoomX, zoomX);

    // Recompute how many sample steps from right edge to that same sample under new zoom
    const double kNew          = dxFromRight / zoomX;
    const double rightIndexNew = sampleAtMouse + kNew;

    double newOffset = (double) latestIndex - rightIndexNew;
    newOffset = juce::jlimit (0.0, maxOffsetSamples, newOffset);

    viewOffsetSamples = newOffset;
    hasCustomZoomX    = true;
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

    if (event.mods.isShiftDown())
        applyVerticalZoom (wheel.deltaY);
    else
        applyHorizontalZoom (wheel.deltaY, event.position.x);

    repaint();
}

void VolumeHistoryComponent::mouseDown (const juce::MouseEvent& event)
{
    juce::ignoreUnused (event);

    // Toggle between RMS and Peak mode on mouse click
    envelopeMode = (envelopeMode == EnvelopeMode::RMS
                      ? EnvelopeMode::Peak
                      : EnvelopeMode::RMS);

    repaint();
}
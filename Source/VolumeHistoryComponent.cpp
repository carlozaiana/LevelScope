#include "VolumeHistoryComponent.h"
#include "PluginProcessor.h"

#include <cmath>
#include <algorithm>

//==============================================================================

VolumeHistoryComponent::VolumeHistoryComponent (LevelScopeAudioProcessor& proc)
    : processor (proc),
      visualSampleRate (proc.getVisualSampleRate()),
      historyLengthSeconds (180.0),       // spec: 180 s total history
      minDb (-90.0f),
      maxDb (  0.0f),
      baseDbRange (std::abs (minDb)),     // 90 dB
      historyCapacitySamples ((int) std::ceil (historyLengthSeconds * visualSampleRate)),
      historyDb ((size_t) historyCapacitySamples, minDb),
      visibleDurationSeconds (10.0),      // spec: default 10 s window
      minVisibleSeconds (0.1),
      maxVisibleSeconds (historyLengthSeconds)
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
    // No special layout: fills entire editor.
}

//==============================================================================

void VolumeHistoryComponent::timerCallback()
{
    drainProcessorFifo();
    repaint();
}

void VolumeHistoryComponent::drainProcessorFifo()
{
    // Drain all available RMS values from the processor into our history buffer
    constexpr int chunkSize = 2048;
    float temp[chunkSize];

    for (;;)
    {
        const int numRead = processor.readRmsFromFifo (temp, chunkSize);
        if (numRead <= 0)
            break;

        pushRmsBatchToHistory (temp, numRead);
    }
}

void VolumeHistoryComponent::pushRmsBatchToHistory (const float* values, int numValues)
{
    if (numValues <= 0)
        return;

    for (int i = 0; i < numValues; ++i)
    {
        const float rms = values[i];

        // Convert RMS (linear) to dB, clamped to minDb
        const float db = juce::Decibels::gainToDecibels (rms, minDb);

        const juce::int64 writeIndex = historySampleCount;
        const int ringIndex = (int) (writeIndex % (juce::int64) historyCapacitySamples);

        historyDb[(size_t) ringIndex] = db;
        ++historySampleCount;
    }

    // "Now" has advanced by numValues samples; keep the same offset behind now.
    // (viewOffsetSamples is in samples from now, so it stays the same)
}

//==============================================================================

float VolumeHistoryComponent::sampleIndexToDb (double sampleIndex) const noexcept
{
    if (historySampleCount <= 0)
        return minDb;

    const juce::int64 latestIndex = historySampleCount - 1;
    const juce::int64 earliestIndex =
        (historySampleCount > (juce::int64) historyCapacitySamples
            ? historySampleCount - (juce::int64) historyCapacitySamples
            : 0);

    if (sampleIndex < (double) earliestIndex || sampleIndex > (double) latestIndex)
        return minDb;

    const juce::int64 index0 = (juce::int64) std::floor (sampleIndex);
    juce::int64 index1       = (juce::int64) std::ceil (sampleIndex);

    if (index1 > latestIndex)
        index1 = latestIndex;

    const double frac = sampleIndex - (double) index0;

    const int ringIndex0 = (int) (index0 % (juce::int64) historyCapacitySamples);
    const int ringIndex1 = (int) (index1 % (juce::int64) historyCapacitySamples);

    const float v0 = historyDb[(size_t) ringIndex0];
    const float v1 = historyDb[(size_t) ringIndex1];

    return v0 + (float) ((v1 - v0) * frac);
}

float VolumeHistoryComponent::dbToY (float db, float height) const noexcept
{
    if (height <= 0.0f)
        return 0.0f;

    // Visible dB range shrinks as yZoom increases (zoom in)
    const float effectiveRange = (float) (baseDbRange / yZoom);
    const float topDb = maxDb;
    const float bottomDb = topDb - effectiveRange;

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

    const int width  = (int) bounds.getWidth();
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

    const double nowIndex      = (double) historySampleCount;
    const double viewEndSample = nowIndex - viewOffsetSamples;
    const double visibleSamples = visibleDurationSeconds * visualSampleRate;

    juce::Path path;
    bool started = false;

    g.setColour (juce::Colours::limegreen);

    for (int x = 0; x < width; ++x)
    {
        const double fracX = (width > 1 ? (double) x / (double) (width - 1) : 0.0);

        // Right edge is viewEndSample, left edge is viewEndSample - visibleSamples
        const double sampleIndex = viewEndSample - (1.0 - fracX) * visibleSamples;

        const float db = sampleIndexToDb (sampleIndex);
        const float y  = dbToY (db, height);

        if (! started)
        {
            path.startNewSubPath ((float) x + 0.5f, y);
            started = true;
        }
        else
        {
            path.lineTo ((float) x + 0.5f, y);
        }
    }

    g.strokePath (path, juce::PathStrokeType (1.5f));
}

//==============================================================================

void VolumeHistoryComponent::applyHorizontalZoom (float wheelDelta, float mouseX)
{
    if (historySampleCount <= 0 || getWidth() <= 1 || wheelDelta == 0.0f)
        return;

    const double nowIndex = (double) historySampleCount;
    const double viewEndSample = nowIndex - viewOffsetSamples;

    double visibleSamples = visibleDurationSeconds * visualSampleRate;

    const double mouseRatio =
        juce::jlimit (0.0, 1.0, (double) mouseX / (double) (getWidth() - 1));

    // Current sample under the mouse cursor
    const double sampleAtMouse = viewEndSample - (1.0 - mouseRatio) * visibleSamples;

    // Zoom factor: positive deltaY = zoom in (shorter window)
    const double zoomBase   = 1.1;
    const double zoomFactor = std::pow (zoomBase, (double) wheelDelta);

    visibleSamples /= zoomFactor;

    const double minVisibleSamples = minVisibleSeconds * visualSampleRate;
    const double maxVisibleSamples = maxVisibleSeconds * visualSampleRate;

    visibleSamples = juce::jlimit (minVisibleSamples, maxVisibleSamples, visibleSamples);
    visibleDurationSeconds = visibleSamples / visualSampleRate;

    // New end sample so that the sample at mouse stays under the cursor
    double newViewEndSample = sampleAtMouse + (1.0 - mouseRatio) * visibleSamples;

    // Enforce that we don't look into the future
    if (newViewEndSample > nowIndex)
        newViewEndSample = nowIndex;

    // Convert changed end-sample into offset (in samples) behind "now"
    double newOffset = nowIndex - newViewEndSample;
    if (newOffset < 0.0)
        newOffset = 0.0;

    // Don't allow offset to exceed the total number of samples (can't go before zero)
    if (newOffset > nowIndex)
        newOffset = nowIndex;

    viewOffsetSamples = newOffset;
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
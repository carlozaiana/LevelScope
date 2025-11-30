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
      historyLengthSeconds (1800.0),      // 30 minutes
      minDb (-90.0f),
      maxDb (  0.0f),
      baseDbRange (std::abs (minDb)),
      rawCapacityFrames ((int) std::ceil (historyLengthSeconds * visualFrameRate)),
      overviewCapacityFrames (juce::jmax (1, rawCapacityFrames / decimationFactor)),
      rawHistory ((size_t) rawCapacityFrames),
      overviewHistory ((size_t) overviewCapacityFrames)
{
    jassert (visualFrameRate > 0.0);
    jassert (rawCapacityFrames > 0);
    jassert (overviewCapacityFrames > 0);

    currentOverviewMax.momentaryDb = minDb;
    currentOverviewMax.shortTermDb = minDb;

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
        Frame f;
        f.momentaryDb = juce::Decibels::gainToDecibels (rmsM, minDb);
        f.shortTermDb = juce::Decibels::gainToDecibels (rmsS, minDb);

        // RAW history (ring buffer)
        rawHistory[(size_t) rawWriteIndex] = f;
        rawWriteIndex = (rawWriteIndex + 1) % rawCapacityFrames;
        ++totalRawFrames;

        // OVERVIEW accumulation (fixed decimationFactor)
        if (f.momentaryDb > currentOverviewMax.momentaryDb) currentOverviewMax.momentaryDb = f.momentaryDb;
        if (f.shortTermDb > currentOverviewMax.shortTermDb) currentOverviewMax.shortTermDb = f.shortTermDb;

        ++currentOverviewCount;

        if (currentOverviewCount >= decimationFactor)
        {
            overviewHistory[(size_t) overviewWriteIndex] = currentOverviewMax;
            overviewWriteIndex = (overviewWriteIndex + 1) % overviewCapacityFrames;
            ++totalOverviewFrames;

            currentOverviewMax.momentaryDb = minDb;
            currentOverviewMax.shortTermDb = minDb;
            currentOverviewCount = 0;
        }
    }
}

//==============================================================================

VolumeHistoryComponent::Frame VolumeHistoryComponent::getRawFrameAgo (int framesAgo) const noexcept
{
    Frame out;
    out.momentaryDb = minDb;
    out.shortTermDb = minDb;

    const juce::int64 available = std::min<juce::int64> (rawCapacityFrames, totalRawFrames);
    if (framesAgo < 0 || (juce::int64) framesAgo >= available)
        return out;

    const int latestIndexInRing = (rawWriteIndex - 1 + rawCapacityFrames) % rawCapacityFrames;
    int idx = latestIndexInRing - framesAgo;
    if (idx < 0) idx += rawCapacityFrames;

    return rawHistory[(size_t) idx];
}

VolumeHistoryComponent::Frame VolumeHistoryComponent::getOverviewFrameAgo (int groupsAgo) const noexcept
{
    Frame out;
    out.momentaryDb = minDb;
    out.shortTermDb = minDb;

    const juce::int64 available = std::min<juce::int64> (overviewCapacityFrames, totalOverviewFrames);
    if (groupsAgo < 0 || (juce::int64) groupsAgo >= available)
        return out;

    const int latestIndexInRing = (overviewWriteIndex - 1 + overviewCapacityFrames) % overviewCapacityFrames;
    int idx = latestIndexInRing - groupsAgo;
    if (idx < 0) idx += overviewCapacityFrames;

    return overviewHistory[(size_t) idx];
}

//==============================================================================

float VolumeHistoryComponent::dbToY (float db, float height) const noexcept
{
    if (height <= 0.0f)
        return 0.0f;

    const float effectiveRange = (float) (baseDbRange / zoomY);
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

    // Horizontal reference lines (e.g. -90, -60, -30, 0 dB)
    g.setColour (juce::Colours::darkgrey.withMultipliedAlpha (0.6f));
    const int numLines = 4;
    for (int i = 0; i <= numLines; ++i)
    {
        const float db = maxDb - (baseDbRange / (float) numLines) * (float) i;
        const float y  = dbToY (db, height);
        g.drawHorizontalLine ((int) std::round (y), 0.0f, bounds.getWidth());
    }

    const juce::int64 availableRaw = std::min<juce::int64> (rawCapacityFrames, totalRawFrames);
    if (availableRaw < 2)
        return;

    const float w = bounds.getWidth();
    const float h = bounds.getHeight();
    const float midY = h * 0.5f;

    // Decide whether to use RAW or OVERVIEW based on zoomX (pixels per raw frame)
    const bool useOverview = (zoomX < 0.05); // heuristic threshold, like SmoothScope

    juce::Path pathMomentary;
    juce::Path pathShortTerm;
    bool startedM = false;
    bool startedS = false;

    if (useOverview)
    {
        // OVERVIEW mode: decimated data, 1 point per group of decimationFactor frames
        const int groupsAvailable = (int) std::min<juce::int64> (overviewCapacityFrames, totalOverviewFrames);
        if (groupsAvailable <= 0)
            return;

        const float pointSpacing = (float) (zoomX * decimationFactor);

        for (int i = 0; i < groupsAvailable; ++i)
        {
            const float x = w - ((float) i * pointSpacing);
            if (x < -10.0f)
                break;

            const Frame f = getOverviewFrameAgo (i);

            const float yM = dbToY (f.momentaryDb, h);
            const float yS = dbToY (f.shortTermDb, h);

            if (! startedM) { pathMomentary.startNewSubPath (x, yM); startedM = true; }
            else            { pathMomentary.lineTo          (x, yM); }

            if (! startedS) { pathShortTerm.startNewSubPath (x, yS); startedS = true; }
            else            { pathShortTerm.lineTo          (x, yS); }
        }
    }
    else
    {
        // RAW mode: full resolution, iterate RAW frames with step = 1
        const int samplesToDraw = (int) std::min<juce::int64> (
            (juce::int64) std::ceil (w / zoomX) + 2, availableRaw);

        for (int i = 0; i < samplesToDraw; ++i)
        {
            const float x = w - ((float) i * (float) zoomX);
            if (x < -10.0f)
                break;

            const Frame f = getRawFrameAgo (i);

            const float yM = dbToY (f.momentaryDb, h);
            const float yS = dbToY (f.shortTermDb, h);

            if (! startedM) { pathMomentary.startNewSubPath (x, yM); startedM = true; }
            else            { pathMomentary.lineTo          (x, yM); }

            if (! startedS) { pathShortTerm.startNewSubPath (x, yS); startedS = true; }
            else            { pathShortTerm.lineTo          (x, yS); }
        }
    }

    // Draw curves
    if (startedS)
    {
        g.setColour (juce::Colours::cyan.withMultipliedAlpha (0.8f));
        g.strokePath (pathShortTerm, juce::PathStrokeType (1.5f));
    }

    if (startedM)
    {
        g.setColour (juce::Colours::limegreen);
        g.strokePath (pathMomentary, juce::PathStrokeType (2.0f));
    }

    //==========================================================================
    // Time ruler at bottom: HH:MM:SS, auto-adjusted spacing
    //==========================================================================

    const float rulerHeight   = 16.0f;
    const float rulerBaseY    = height - 2.0f;
    const float tickTopY      = rulerBaseY - 6.0f;
    const float textTopY      = rulerBaseY - 14.0f;
    const float leftPixel     = 0.0f;

    const double widthPixelsForTime = (double) width;

    if (widthPixelsForTime > 1.0 && zoomX > 0.0 && visualFrameRate > 0.0)
    {
        const double totalTimeSeconds = (double) totalRawFrames / visualFrameRate;
        const double maxHistorySeconds = (double) rawCapacityFrames / visualFrameRate;

        // Approx visible frames in RAW units
        double visibleFramesApprox = widthPixelsForTime / zoomX;
        if (visibleFramesApprox > (double) availableRaw)
            visibleFramesApprox = (double) availableRaw;

        if (visibleFramesApprox > 1.0)
        {
            const double visibleSeconds = visibleFramesApprox / visualFrameRate;

            const double tRight = totalTimeSeconds;
            const double tLeftLimit = juce::jmax (0.0, tRight - maxHistorySeconds);
            double tLeft = tRight - visibleSeconds;
            if (tLeft < tLeftLimit)
                tLeft = tLeftLimit;

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
                const double norm = (t - tLeft) / (visibleSeconds > 1e-9 ? visibleSeconds : 1.0);
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

    // Overlay info
    g.setColour (juce::Colours::white);
    g.setFont (14.0f);
    juce::String mode = useOverview ? "OVERVIEW" : "RAW";
    g.drawText (mode + " | ZoomX: " + juce::String (zoomX, 4) + " | ZoomY: " + juce::String (zoomY, 2),
                8, 8, 300, 20, juce::Justification::topLeft);
}

//==============================================================================

void VolumeHistoryComponent::applyHorizontalZoom (float wheelDelta)
{
    if (historyFrameCount <= 1 || getWidth() <= 1 || wheelDelta == 0.0f)
        return;

    const double zoomBase   = 1.1;
    const double zoomFactor = std::pow (zoomBase, (double) wheelDelta);

    zoomX *= zoomFactor;
    zoomX = juce::jlimit (minZoomX, maxZoomX, zoomX);

    hasCustomZoomX = true;
}

void VolumeHistoryComponent::applyVerticalZoom (float wheelDelta)
{
    if (wheelDelta == 0.0f)
        return;

    const double zoomBase   = 1.1;
    const double zoomFactor = std::pow (zoomBase, (double) wheelDelta);

    zoomY *= zoomFactor;
    zoomY = juce::jlimit (minZoomY, maxZoomY, zoomY);
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
        applyHorizontalZoom (wheel.deltaY);

    repaint();
}
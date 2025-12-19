#include "VolumeHistoryComponent.h"
#include "PluginProcessor.h"

#include <cmath>
#include <algorithm>
#include <limits>

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
// Constructor / Destructor
//==============================================================================

VolumeHistoryComponent::VolumeHistoryComponent (LevelScopeAudioProcessor& proc)
    : processor (proc),
      visualFrameRate (proc.getLoudnessFrameRate()),
      historyLengthSeconds (3.0 * 3600.0),    // 3 hours
      minDb (-90.0f),
      maxDb (  0.0f),
      baseDbRange (std::abs (minDb))
{
    jassert (visualFrameRate > 0.0);

    initialiseHistoryLevels();
    resetHistoryLevels();

    setOpaque (true);
    startTimerHz (60); // GUI update timer at 60 Hz
}

VolumeHistoryComponent::~VolumeHistoryComponent()
{
    stopTimer();
}

//==============================================================================
// [HISTORY-INIT]
//==============================================================================

void VolumeHistoryComponent::initialiseHistoryLevels()
{
    // RAW level capacity: number of frames to cover 3 hours at visualFrameRate.
    rawCapacityFrames = (int) std::ceil (historyLengthSeconds * visualFrameRate);
    if (rawCapacityFrames < 1)
        rawCapacityFrames = 1;

    // Level 0 (RAW)
    {
        auto& L0 = levels[0];
        L0.levelIndex     = 0;
        L0.groupsPerGroup = 1;   // not used at L0
        L0.spanFrames     = 1;   // 1 RAW frame per group
        L0.capacity       = rawCapacityFrames;
        L0.groups.assign ((size_t) rawCapacityFrames,
                          FrameGroup { minDb, minDb, minDb, minDb });
        L0.writeIndex     = 0;
        L0.totalGroups    = 0;
        L0.pendingCount   = 0;
    }

    // Derived levels (L1..LmaxLevels-1)
    int prevSpanFrames = levels[0].spanFrames;

    for (int level = 1; level < maxLevels; ++level)
    {
        auto& L = levels[(size_t) level];

        L.levelIndex     = level;
        L.groupsPerGroup = groupsPerLevel;
        L.spanFrames     = prevSpanFrames * groupsPerLevel;
        prevSpanFrames   = L.spanFrames;

        int capacity = (int) std::ceil ((double) rawCapacityFrames / (double) L.spanFrames);
        if (capacity < 1)
            capacity = 1;

        L.capacity     = capacity;
        L.groups.assign ((size_t) capacity,
                         FrameGroup { minDb, minDb, minDb, minDb });
        L.writeIndex   = 0;
        L.totalGroups  = 0;
        L.pendingCount = 0;

        L.pending.momentaryMinDb =  std::numeric_limits<float>::infinity();
        L.pending.momentaryMaxDb = -std::numeric_limits<float>::infinity();
        L.pending.shortTermMinDb =  std::numeric_limits<float>::infinity();
        L.pending.shortTermMaxDb = -std::numeric_limits<float>::infinity();
    }
}

void VolumeHistoryComponent::resetHistoryLevels()
{
    for (auto& L : levels)
    {
        L.writeIndex   = 0;
        L.totalGroups  = 0;
        L.pendingCount = 0;

        for (auto& g : L.groups)
        {
            g.momentaryMinDb = minDb;
            g.momentaryMaxDb = minDb;
            g.shortTermMinDb = minDb;
            g.shortTermMaxDb = minDb;
        }

        L.pending.momentaryMinDb =  std::numeric_limits<float>::infinity();
        L.pending.momentaryMaxDb = -std::numeric_limits<float>::infinity();
        L.pending.shortTermMinDb =  std::numeric_limits<float>::infinity();
        L.pending.shortTermMaxDb = -std::numeric_limits<float>::infinity();
    }
}

//==============================================================================
// [TIMER]
//==============================================================================

void VolumeHistoryComponent::timerCallback()
{
    drainProcessorFifo();
    repaint();
}

//==============================================================================
// [HISTORY-UPDATE]
//==============================================================================

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

        for (int i = 0; i < numRead; ++i)
            pushFrameToHistory (momentaryValues[i], shortTermValues[i]);
    }
}

void VolumeHistoryComponent::pushFrameToHistory (float momentaryRms,
                                                 float shortTermRms)
{
    const float dbM = juce::Decibels::gainToDecibels (momentaryRms, minDb);
    const float dbS = juce::Decibels::gainToDecibels (shortTermRms, minDb);

    FrameGroup fg;
    fg.momentaryMinDb = dbM;
    fg.momentaryMaxDb = dbM;
    fg.shortTermMinDb = dbS;
    fg.shortTermMaxDb = dbS;

    // Level 0 (RAW): each frame is a fully-formed group.
    writeGroupToLevel (0, fg);

    // Propagate upwards: each higher level accumulates lower-level groups.
    accumulateToHigherLevels (1, fg);
}

void VolumeHistoryComponent::writeGroupToLevel (int levelIndex, const FrameGroup& group)
{
    if (levelIndex < 0 || levelIndex >= maxLevels)
        return;

    auto& L = levels[(size_t) levelIndex];

    if (L.capacity <= 0)
        return;

    L.groups[(size_t) L.writeIndex] = group;
    L.writeIndex = (L.writeIndex + 1) % L.capacity;
    ++L.totalGroups;
}

void VolumeHistoryComponent::accumulateToHigherLevels (int levelIndex,
                                                       const FrameGroup& sourceGroup)
{
    if (levelIndex < 0 || levelIndex >= maxLevels)
        return;

    auto& L = levels[(size_t) levelIndex];

    if (L.capacity <= 0 || L.groupsPerGroup <= 0)
        return;

    if (L.pendingCount == 0)
    {
        // First group in this new aggregate.
        L.pending = sourceGroup;
    }
    else
    {
        // Expand min/max to include the incoming group's range.
        L.pending.momentaryMinDb = std::min (L.pending.momentaryMinDb, sourceGroup.momentaryMinDb);
        L.pending.momentaryMaxDb = std::max (L.pending.momentaryMaxDb, sourceGroup.momentaryMaxDb);
        L.pending.shortTermMinDb = std::min (L.pending.shortTermMinDb, sourceGroup.shortTermMinDb);
        L.pending.shortTermMaxDb = std::max (L.pending.shortTermMaxDb, sourceGroup.shortTermMaxDb);
    }

    ++L.pendingCount;

    if (L.pendingCount >= L.groupsPerGroup)
    {
        // Aggregate is complete: write it to this level, reset pending, and
        // propagate the finished group to the next higher level.
        FrameGroup finished = L.pending;

        L.pendingCount = 0;
        L.pending.momentaryMinDb =  std::numeric_limits<float>::infinity();
        L.pending.momentaryMaxDb = -std::numeric_limits<float>::infinity();
        L.pending.shortTermMinDb =  std::numeric_limits<float>::infinity();
        L.pending.shortTermMaxDb = -std::numeric_limits<float>::infinity();

        writeGroupToLevel (levelIndex, finished);

        // Recurse upwards
        accumulateToHigherLevels (levelIndex + 1, finished);
    }
}

//==============================================================================
// [HISTORY-ACCESS]
//==============================================================================

int VolumeHistoryComponent::getAvailableGroups (int levelIndex) const noexcept
{
    if (levelIndex < 0 || levelIndex >= maxLevels)
        return 0;

    const auto& L = levels[(size_t) levelIndex];
    const auto available = std::min<juce::int64> ((juce::int64) L.capacity, L.totalGroups);

    return (int) available;
}

int VolumeHistoryComponent::getPendingFramesAtLevel (int levelIndex) const noexcept
{
    if (levelIndex <= 0 || levelIndex >= maxLevels)
        return 0;

    const auto& L     = levels[(size_t) levelIndex];
    const auto& Lprev = levels[(size_t) (levelIndex - 1)];

    if (L.pendingCount <= 0)
        return 0;

    return L.pendingCount * Lprev.spanFrames;
}

VolumeHistoryComponent::FrameGroup VolumeHistoryComponent::getGroupAgo (int levelIndex,
                                                                        int groupsAgo) const noexcept
{
    FrameGroup out { minDb, minDb, minDb, minDb };

    if (levelIndex < 0 || levelIndex >= maxLevels)
        return out;

    const auto& L = levels[(size_t) levelIndex];

    const juce::int64 available = std::min<juce::int64> ((juce::int64) L.capacity, L.totalGroups);
    if (groupsAgo < 0 || (juce::int64) groupsAgo >= available || L.capacity <= 0)
        return out;

    const int latestIndexInRing = (L.writeIndex - 1 + L.capacity) % L.capacity;
    int idx = latestIndexInRing - groupsAgo;
    if (idx < 0) idx += L.capacity;

    return L.groups[(size_t) idx];
}

juce::int64 VolumeHistoryComponent::getTotalFramesL0() const noexcept
{
    // Level 0: one group per RAW loudness frame.
    const auto& L0 = levels[0];
    return L0.totalGroups;
}

//==============================================================================
// [LOD-SELECTION]
//==============================================================================

int VolumeHistoryComponent::selectBestLevelForCurrentZoom() const noexcept
{
    // framesPerPixel: how many RAW frames correspond to 1 pixel horizontally.
    const double safeZoomX       = (zoomX > 1.0e-9 ? zoomX : 1.0e-9);
    const double framesPerPixel  = 1.0 / safeZoomX;

    int    bestLevel = 0;
    double bestScore = std::numeric_limits<double>::infinity();

    for (int level = 0; level < maxLevels; ++level)
    {
        const auto& L = levels[(size_t) level];
        const int   available = getAvailableGroups (level);

        if (available < 2)
            continue; // skip empty/inadequate levels

        const double span = (double) L.spanFrames;
        if (span <= 0.0)
            continue;

        const double ratio = span / framesPerPixel; // we want ratio ~ 1
        if (ratio <= 0.0)
            continue;

        // Score: closeness of ratio to 1 in log2 space (so 0.5 and 2 both give score 1).
        const double score = std::abs (std::log2 (ratio));

        if (score < bestScore)
        {
            bestScore = score;
            bestLevel = level;
        }
    }

    return bestLevel;
}

void VolumeHistoryComponent::buildVisibleGroupsForLevel (int levelIndex,
                                                         int widthPixels,
                                                         std::vector<FrameGroup>& outGroups,
                                                         std::vector<int>& outFramesAgo) const
{
    outGroups.clear();
    outFramesAgo.clear();

    if (levelIndex < 0 || levelIndex >= maxLevels)
        return;

    const auto& L = levels[(size_t) levelIndex];

    const int availableGroups = getAvailableGroups (levelIndex);
    if (availableGroups <= 0 || widthPixels <= 0)
        return;

    const int   spanFrames     = L.spanFrames;
    const int   pendingFrames  = getPendingFramesAtLevel (levelIndex);
    const float overscanPixels = 10.0f;

    if (spanFrames <= 0 || zoomX <= 0.0)
        return;

    const double maxFramesVisible = (double) (widthPixels + overscanPixels) / zoomX;

    if (maxFramesVisible <= 0.0)
        return;

    int maxGroupsByX = 0;
    {
        const double numerator = maxFramesVisible - (double) pendingFrames;
        if (numerator >= 0.0)
        {
            const double maxGroupsFloat = numerator / (double) spanFrames;
            maxGroupsByX = (int) std::floor (maxGroupsFloat) + 1; // +1 for groupsAgo = 0
        }
    }

    if (maxGroupsByX <= 0)
        return;

    const int groupsToUse = juce::jlimit (0, availableGroups, maxGroupsByX);
    if (groupsToUse <= 0)
        return;

    outGroups.resize ((size_t) groupsToUse);
    outFramesAgo.resize ((size_t) groupsToUse);

    // We want outGroups in chronological order: index 0 = oldest, last = most recent.
    // getGroupAgo(0) = most recent, getGroupAgo(1) = just before that, etc.

    for (int groupsAgo = 0; groupsAgo < groupsToUse; ++groupsAgo)
    {
        const FrameGroup g = getGroupAgo (levelIndex, groupsAgo);
        const int        chronologicalIndex = groupsToUse - 1 - groupsAgo;

        outGroups[(size_t) chronologicalIndex]     = g;
        outFramesAgo[(size_t) chronologicalIndex] = pendingFrames + groupsAgo * spanFrames;
    }
}

//==============================================================================
// [REP-LINE]
//==============================================================================

void VolumeHistoryComponent::computeRepresentativeCurves (const std::vector<FrameGroup>& groups,
                                                          std::vector<float>& repMomentary,
                                                          std::vector<float>& repShortTerm) const
{
    const size_t n = groups.size();

    repMomentary.clear();
    repShortTerm.clear();

    if (n == 0)
        return;

    repMomentary.resize (n);
    repShortTerm.resize (n);

    const float epsilonTrend = 0.1f;  // dB threshold to avoid flicker

    float prevCenterM = 0.0f;
    float prevCenterS = 0.0f;

    bool hasPrev = false;

    for (size_t i = 0; i < n; ++i)
    {
        const auto& g = groups[i];

        const float minM    = g.momentaryMinDb;
        const float maxM    = g.momentaryMaxDb;
        const float centerM = 0.5f * (minM + maxM);

        const float minS    = g.shortTermMinDb;
        const float maxS    = g.shortTermMaxDb;
        const float centerS = 0.5f * (minS + maxS);

        float rawRepM = centerM;
        float rawRepS = centerS;

        if (! hasPrev)
        {
            // First point: use center
            hasPrev     = true;
            prevCenterM = centerM;
            prevCenterS = centerS;
        }
        else
        {
            const float trendM = centerM - prevCenterM;
            const float trendS = centerS - prevCenterS;

            // Hard version:
            //   trend > +ε  -> use max
            //   trend < -ε  -> use min
            //   otherwise   -> center
            float alphaM = 0.5f;
            if (trendM >  epsilonTrend)      alphaM = 1.0f;  // upward -> max
            else if (trendM < -epsilonTrend) alphaM = 0.0f;  // downward -> min

            float alphaS = 0.5f;
            if (trendS >  epsilonTrend)      alphaS = 1.0f;
            else if (trendS < -epsilonTrend) alphaS = 0.0f;

            // Interpolate inside [min, max]
            rawRepM = minM + alphaM * (maxM - minM);
            rawRepS = minS + alphaS * (maxS - minS);

            // Clamp to band, just to be safe
            rawRepM = juce::jlimit (minM, maxM, rawRepM);
            rawRepS = juce::jlimit (minS, maxS, rawRepS);

            prevCenterM = centerM;
            prevCenterS = centerS;
        }

        // Smoothing disabled: use rawRep directly
        repMomentary[i] = rawRepM;
        repShortTerm[i] = rawRepS;
    }
}

//==============================================================================
// [GEOMETRY HELPERS]
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
// [COMPONENT LIFECYCLE]
//==============================================================================

void VolumeHistoryComponent::resized()
{
    // Initialise X zoom so that, roughly, 10 seconds are visible by default,
    // unless the user has already chosen a custom zoom.
    if (! hasCustomZoomX)
    {
        const int w = getWidth();
        if (w > 0 && visualFrameRate > 0.0)
        {
            const double desiredVisibleSeconds = 10.0;
            zoomX = (double) w / (desiredVisibleSeconds * visualFrameRate);
            zoomX = juce::jlimit (minZoomX, maxZoomX, zoomX);
        }
    }
}

//==============================================================================
// [DRAW]
//==============================================================================

void VolumeHistoryComponent::paint (juce::Graphics& g)
{
    auto bounds = getLocalBounds().toFloat();

    // Background: darker variant of the momentary lime
    auto backgroundColour = juce::Colours::limegreen.darker (3.0f);
    g.fillAll (backgroundColour);

    const int   width  = (int) bounds.getWidth();
    const float height = bounds.getHeight();

    if (width <= 1 || height <= 0.0f)
        return;

    // Horizontal reference lines (e.g. -90, -60, -30, 0 dB)
    g.setColour (juce::Colours::darkgrey.withMultipliedAlpha (0.4f));
    const int numLines = 4;
    for (int i = 0; i <= numLines; ++i)
    {
        const float db = maxDb - (baseDbRange / (float) numLines) * (float) i;
        const float y  = dbToY (db, height);
        g.drawHorizontalLine ((int) std::round (y), 0.0f, bounds.getWidth());
    }

    const juce::int64 totalFrames   = getTotalFramesL0();
    const juce::int64 availableRaw  = std::min<juce::int64> ((juce::int64) rawCapacityFrames, totalFrames);

    if (availableRaw < 2)
        return;

    const float w = bounds.getWidth();
    const float h = bounds.getHeight();

    // Choose best level-of-detail based on current zoom
    const int selectedLevel = selectBestLevelForCurrentZoom();

    // Build visible groups and their framesAgo
    std::vector<FrameGroup> groups;
    std::vector<int>        framesAgo;
    buildVisibleGroupsForLevel (selectedLevel, width, groups, framesAgo);

        const size_t n = groups.size();
    if (n < 2)
        return;

    // Compute representative curves inside bands
    std::vector<float> repM, repS;
    computeRepresentativeCurves (groups, repM, repS);

    //==========================================================================
    // [LINE-LOD] Decide how (or if) to draw lines, based on level & group count
    //==========================================================================

    const size_t maxFullResLinePoints   = 1000;  // up to this: draw every point
    const size_t maxDecimatedLinePoints = 2000;  // above this: no lines, bands only

    bool drawLinesThisFrame = showLines; // user toggle
    int  lineDecimation     = 1;         // 1 = no decimation; 2 = every 2nd, etc.

    if (drawLinesThisFrame)
    {
        if (selectedLevel >= 4)
        {
            // Levels 4 and 5: bands only, no lines.
            drawLinesThisFrame = false;
        }
        else
        {
            if (n > maxDecimatedLinePoints)
            {
                // Too many groups -> bands only.
                drawLinesThisFrame = false;
            }
            else if (n > maxFullResLinePoints)
            {
                // Moderate overload: decimate line points.
                lineDecimation = 2; // draw every 2nd point
            }
        }
    }

    // Band thickness: at coarse levels, make bands a bit thicker.
    const float shortTermThickness = (selectedLevel >= 4 ? 1.5f : 1.0f);
    const float momentaryThickness = (selectedLevel >= 4 ? 2.0f : 1.2f);

    // Draw bands and lines
    juce::Path pathRepM, pathRepS;
    bool startedRepM = false, startedRepS = false;

        for (size_t i = 0; i < n; ++i)
    {
        const float x = w - (float) framesAgo[i] * (float) zoomX;
        if (x < -10.0f)
            continue; // off-screen to the left

        const auto& gGroup = groups[i];

        const float yMM = dbToY (gGroup.momentaryMaxDb, h);
        const float yMm = dbToY (gGroup.momentaryMinDb, h);
        const float ySM = dbToY (gGroup.shortTermMaxDb, h);
        const float ySm = dbToY (gGroup.shortTermMinDb, h);

        const float yRepM = dbToY (repM[i], h);
        const float yRepS = dbToY (repS[i], h);

        //======================================================================
        // Bands (vertical min->max lines)
        //======================================================================
        //
        // - Level 0: no bands (lines only).
        // - Level >= 1: selective bands:
        //     draw only if internal range (max-min) >= threshold.
        //======================================================================

        if (showBands && selectedLevel > 0) // L0: no bands at all
        {
            // Per-group dynamic range in dB
            const float rangeMM = gGroup.momentaryMaxDb - gGroup.momentaryMinDb;
            const float rangeSM = gGroup.shortTermMaxDb - gGroup.shortTermMinDb;

            // Default: draw bands for this group
            bool drawMomentaryBand = true;
            bool drawShortTermBand = true;

            // At all levels >= 1: only draw bands where there is real variation.
            const float bandRangeThresholdDb = 3.0f; // tweak: 3–6 dB typical

            drawMomentaryBand = (rangeMM >= bandRangeThresholdDb);
            drawShortTermBand = (rangeSM >= bandRangeThresholdDb);

            // Short-term band (cyan-ish, behind)
            if (drawShortTermBand)
            {
                g.setColour (juce::Colours::cyan.withMultipliedAlpha (0.6f));
                g.drawLine (x, ySM, x, ySm, shortTermThickness);
            }

            // Momentary band (lime-ish, on top)
            if (drawMomentaryBand)
            {
                g.setColour (juce::Colours::limegreen.withMultipliedAlpha (0.7f));
                g.drawLine (x, yMM, x, yMm, momentaryThickness);
            }
        }

        //======================================================================
        // Lines (representative curves)
        //======================================================================

        if (drawLinesThisFrame)
        {
            // Decimate lines if requested (e.g. every 2nd point).
            if ((lineDecimation <= 1) || ((i % (size_t) lineDecimation) == 0))
            {
                // Short-term representative curve
                if (! startedRepS) { pathRepS.startNewSubPath (x, yRepS); startedRepS = true; }
                else              { pathRepS.lineTo         (x, yRepS); }

                // Momentary representative curve
                if (! startedRepM) { pathRepM.startNewSubPath (x, yRepM); startedRepM = true; }
                else              { pathRepM.lineTo          (x, yRepM); }
            }
        }
    }

        if (drawLinesThisFrame)
    {
        // Short-term representative curve
        g.setColour (juce::Colours::cyan.withMultipliedAlpha (0.9f));
        g.strokePath (pathRepS, juce::PathStrokeType (1.5f));

        // Momentary representative curve
        g.setColour (juce::Colours::limegreen);
        g.strokePath (pathRepM, juce::PathStrokeType (2.0f));
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
        const double totalTimeSeconds   = (double) totalFrames / visualFrameRate;
        const double maxHistorySeconds  = (double) rawCapacityFrames / visualFrameRate;

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

            // Choose tick spacing
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

            g.setColour (juce::Colours::darkgrey.withMultipliedAlpha (0.7f));
            g.drawHorizontalLine ((int) std::round (rulerBaseY), 0.0f, (float) width);

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

    // Overlay info (debug/feedback)
    g.setColour (juce::Colours::white);
    g.setFont (14.0f);

    const auto& Lsel        = levels[(size_t) selectedLevel];
    const double spanSeconds = (double) Lsel.spanFrames / visualFrameRate;

    juce::String info = "Level: " + juce::String (selectedLevel) +
                        " (span " + juce::String (Lsel.spanFrames) + " frames, " +
                        juce::String (spanSeconds, 3) + " s)" +
                        " | ZoomX: " + juce::String (zoomX, 4) +
                        " | ZoomY: " + juce::String (zoomY, 2) +
                        " | Bands: " + juce::String (showBands ? "ON" : "OFF") +
                        " | Lines: " + juce::String (showLines ? "ON" : "OFF");

    g.drawText (info, 8, 8, (int) std::min (w - 16.0f, 520.0f), 20, juce::Justification::topLeft);
}

//==============================================================================
// [ZOOM]
//==============================================================================

void VolumeHistoryComponent::applyHorizontalZoom (float wheelDelta)
{
    if (getWidth() <= 1 || wheelDelta == 0.0f)
        return;

    const double zoomBase   = 1.1;
    const double zoomFactor = std::pow (zoomBase, (double) wheelDelta);

    zoomX *= zoomFactor;
    zoomX  = juce::jlimit (minZoomX, maxZoomX, zoomX);

    hasCustomZoomX = true;
}

void VolumeHistoryComponent::applyVerticalZoom (float wheelDelta)
{
    if (wheelDelta == 0.0f)
        return;

    const double zoomBase   = 1.1;
    const double zoomFactor = std::pow (zoomBase, (double) wheelDelta);

    zoomY *= zoomFactor;
    zoomY  = juce::jlimit (minZoomY, maxZoomY, zoomY);
}

//==============================================================================
// [MOUSE]
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

void VolumeHistoryComponent::mouseDown (const juce::MouseEvent& event)
{
    // Shift-click toggles bands; plain click toggles lines.
    if (event.mods.isShiftDown())
        showBands = ! showBands;
    else
        showLines = ! showLines;

    repaint();
}
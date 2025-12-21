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
      historyLengthSeconds (3.0 * 3600.0),
      minDb (-90.0f),
      maxDb (  0.0f),
      baseDbRange (std::abs (minDb))
{
    jassert (visualFrameRate > 0.0);

    initialiseHistoryLevels();
    resetHistoryLevels();

    setOpaque (true);
    startTimerHz (30); // Timer tick 30 Hz; we repaint only on new data.
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
    rawCapacityFrames = (int) std::ceil (historyLengthSeconds * visualFrameRate);
    if (rawCapacityFrames < 1)
        rawCapacityFrames = 1;

    // Level 0 (RAW)
    {
        auto& L0 = levels[0];
        L0.levelIndex     = 0;
        L0.groupsPerGroup = 1;
        L0.spanFrames     = 1;
        L0.capacity       = rawCapacityFrames;
        L0.groups.assign ((size_t) rawCapacityFrames,
                          FrameGroup { minDb, minDb, minDb, minDb });
        L0.writeIndex     = 0;
        L0.totalGroups    = 0;
        L0.pendingCount   = 0;
    }

    // Derived levels
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
    // [STEP1-PERF] Only repaint if we actually received new data.
    const bool gotNewData = drainProcessorFifo();
    if (gotNewData)
        repaint();
}

//==============================================================================
// [HISTORY-UPDATE]
//==============================================================================

bool VolumeHistoryComponent::drainProcessorFifo()
{
    constexpr int chunkSize = 512;
    float momentaryValues [chunkSize];
    float shortTermValues [chunkSize];

    bool readAny = false;

    for (;;)
    {
        const int numRead = processor.readLoudnessFromFifo (momentaryValues,
                                                            shortTermValues,
                                                            chunkSize);
        if (numRead <= 0)
            break;

        readAny = true;

        for (int i = 0; i < numRead; ++i)
            pushFrameToHistory (momentaryValues[i], shortTermValues[i]);
    }

    return readAny;
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

    writeGroupToLevel (0, fg);
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
        L.pending = sourceGroup;
    }
    else
    {
        L.pending.momentaryMinDb = std::min (L.pending.momentaryMinDb, sourceGroup.momentaryMinDb);
        L.pending.momentaryMaxDb = std::max (L.pending.momentaryMaxDb, sourceGroup.momentaryMaxDb);
        L.pending.shortTermMinDb = std::min (L.pending.shortTermMinDb, sourceGroup.shortTermMinDb);
        L.pending.shortTermMaxDb = std::max (L.pending.shortTermMaxDb, sourceGroup.shortTermMaxDb);
    }

    ++L.pendingCount;

    if (L.pendingCount >= L.groupsPerGroup)
    {
        FrameGroup finished = L.pending;

        L.pendingCount = 0;
        L.pending.momentaryMinDb =  std::numeric_limits<float>::infinity();
        L.pending.momentaryMaxDb = -std::numeric_limits<float>::infinity();
        L.pending.shortTermMinDb =  std::numeric_limits<float>::infinity();
        L.pending.shortTermMaxDb = -std::numeric_limits<float>::infinity();

        writeGroupToLevel (levelIndex, finished);
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
    return levels[0].totalGroups;
}

//==============================================================================
// [LOD-SELECTION]
//==============================================================================

int VolumeHistoryComponent::getMaxDrawablePoints (int widthPixels) const noexcept
{
    // [STEP2-LOD-CAP]
    // Goal: keep draw cost bounded ~ O(component width), not O(history length).
    //
    // 1 point per pixel is typically plenty for a history curve.
    // We allow a small margin so overscan and rounding don't cause oscillation.
    const int w = juce::jmax (1, widthPixels);
    const int target = (int) std::round ((double) w * 1.10); // +10% margin
    return juce::jlimit (256, 8192, target);
}

int VolumeHistoryComponent::selectBestLevelForCurrentZoom (int widthPixels) const noexcept
{
    // [STEP2-LOD-CAP]
    // Choose the *lowest* (most detailed) level where predicted drawable groups
    // stays under our cap. This prevents the "approaching next level" slowdown
    // where you can end up drawing ~2x width just before switching.
    const int maxPoints = getMaxDrawablePoints (widthPixels);

    const double safeZoomX = (zoomX > 1.0e-9 ? zoomX : 1.0e-9);
    const double overscanPixels = 10.0;
    const double maxFramesVisible = ((double) widthPixels + overscanPixels) / safeZoomX;

    int bestLevel = -1;

    for (int level = 0; level < maxLevels; ++level)
    {
        const int available = getAvailableGroups (level);
        if (available < 2)
            continue;

        const auto& L = levels[(size_t) level];
        const int spanFrames = L.spanFrames;
        if (spanFrames <= 0)
            continue;

        const int pendingFrames = getPendingFramesAtLevel (level);

        // Predicted number of groups needed to cover visible time window
        double numerator = maxFramesVisible - (double) pendingFrames;
        if (numerator < 0.0)
            numerator = 0.0;

        const int predicted = (int) std::floor (numerator / (double) spanFrames) + 1;
        const int predictedClamped = juce::jlimit (1, available, predicted);

        if (predictedClamped <= maxPoints)
        {
            bestLevel = level;
            break; // first (lowest) level that stays within budget
        }
    }

    if (bestLevel >= 0)
        return bestLevel;

    // If none fit within the cap, fall back to the highest non-empty level.
    for (int level = maxLevels - 1; level >= 0; --level)
        if (getAvailableGroups (level) >= 2)
            return level;

    return 0;
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

    const int spanFrames = L.spanFrames;
    if (spanFrames <= 0 || zoomX <= 0.0)
        return;

    const int pendingFrames = getPendingFramesAtLevel (levelIndex);

    const double overscanPixels = 10.0;
    const double maxFramesVisible = ((double) widthPixels + overscanPixels) / zoomX;
    if (maxFramesVisible <= 0.0)
        return;

    // How many groups are potentially visible (newest backwards), before any capping.
    int maxGroupsByX = 0;
    {
        double numerator = maxFramesVisible - (double) pendingFrames;
        if (numerator < 0.0)
            numerator = 0.0;

        const double maxGroupsFloat = numerator / (double) spanFrames;
        maxGroupsByX = (int) std::floor (maxGroupsFloat) + 1; // include newest
    }

    if (maxGroupsByX <= 0)
        return;

    const int groupsToUse = juce::jlimit (0, availableGroups, maxGroupsByX);
    if (groupsToUse <= 0)
        return;

    // [STEP2-LOD-CAP] Hard cap output size. If too many groups, aggregate chunks.
    const int maxDrawablePoints = getMaxDrawablePoints (widthPixels);

    const int step = (groupsToUse > maxDrawablePoints
                        ? (int) std::ceil ((double) groupsToUse / (double) maxDrawablePoints)
                        : 1);

    const int outCount = (step > 1
                            ? (int) std::ceil ((double) groupsToUse / (double) step)
                            : groupsToUse);

    if (outCount < 2)
        return;

    outGroups.resize ((size_t) outCount);
    outFramesAgo.resize ((size_t) outCount);

    // We output in chronological order: 0 = oldest, last = newest.
    //
    // We fetch source groups by "groupsAgo" where 0 = newest.
    // To produce chronological output, we walk chunks from oldest->newest.
    const int numChunks = outCount;

    for (int chunkChronoIndex = 0; chunkChronoIndex < numChunks; ++chunkChronoIndex)
    {
        // Oldest chunk corresponds to the largest groupsAgo values.
        const int baseGroupsAgo = (numChunks - 1 - chunkChronoIndex) * step;
        const int endGroupsAgo  = juce::jmin (groupsToUse - 1, baseGroupsAgo + step - 1);

        FrameGroup agg;
        agg.momentaryMinDb =  std::numeric_limits<float>::infinity();
        agg.momentaryMaxDb = -std::numeric_limits<float>::infinity();
        agg.shortTermMinDb =  std::numeric_limits<float>::infinity();
        agg.shortTermMaxDb = -std::numeric_limits<float>::infinity();

        for (int ga = baseGroupsAgo; ga <= endGroupsAgo; ++ga)
        {
            const FrameGroup g = getGroupAgo (levelIndex, ga);

            agg.momentaryMinDb = std::min (agg.momentaryMinDb, g.momentaryMinDb);
            agg.momentaryMaxDb = std::max (agg.momentaryMaxDb, g.momentaryMaxDb);
            agg.shortTermMinDb = std::min (agg.shortTermMinDb, g.shortTermMinDb);
            agg.shortTermMaxDb = std::max (agg.shortTermMaxDb, g.shortTermMaxDb);
        }

        // Place the x-position at the *center* of the aggregated time span.
        const double centerGroupsAgo = 0.5 * (double) (baseGroupsAgo + endGroupsAgo);
        const double framesAgoD = (double) pendingFrames + centerGroupsAgo * (double) spanFrames;
        const int framesAgoI = (int) std::round (framesAgoD);

        outGroups[(size_t) chunkChronoIndex] = agg;
        outFramesAgo[(size_t) chunkChronoIndex] = framesAgoI;
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

    const float epsilonTrend = 0.1f;

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
            hasPrev     = true;
            prevCenterM = centerM;
            prevCenterS = centerS;
        }
        else
        {
            const float trendM = centerM - prevCenterM;
            const float trendS = centerS - prevCenterS;

            float alphaM = 0.5f;
            if (trendM >  epsilonTrend)      alphaM = 1.0f;
            else if (trendM < -epsilonTrend) alphaM = 0.0f;

            float alphaS = 0.5f;
            if (trendS >  epsilonTrend)      alphaS = 1.0f;
            else if (trendS < -epsilonTrend) alphaS = 0.0f;

            rawRepM = minM + alphaM * (maxM - minM);
            rawRepS = minS + alphaS * (maxS - minS);

            rawRepM = juce::jlimit (minM, maxM, rawRepM);
            rawRepS = juce::jlimit (minS, maxS, rawRepS);

            prevCenterM = centerM;
            prevCenterS = centerS;
        }

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

    const float norm = (clamped - bottomDb) / effectiveRange;
    return height * (1.0f - norm);
}

//==============================================================================
// [COMPONENT LIFECYCLE]
//==============================================================================

void VolumeHistoryComponent::resized()
{
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

    // [STEP1-PERF] Reserve scratch buffers.
    {
        const int w = juce::jmax (1, getWidth());
        const size_t reserveCount = (size_t) juce::jlimit (256, 8192, (int) std::round (w * 1.25) + 512);

        scratchVisibleGroups.reserve (reserveCount);
        scratchVisibleFramesAgo.reserve (reserveCount);
        scratchRepMomentaryDb.reserve (reserveCount);
        scratchRepShortTermDb.reserve (reserveCount);
    }
}

//==============================================================================
// [DRAW]
//==============================================================================

void VolumeHistoryComponent::paint (juce::Graphics& g)
{
    auto bounds = getLocalBounds().toFloat();

    auto backgroundColour = juce::Colours::limegreen.darker (3.0f);
    g.fillAll (backgroundColour);

    const int   width  = (int) bounds.getWidth();
    const float height = bounds.getHeight();

    if (width <= 1 || height <= 0.0f)
        return;

    // Reference lines
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

    // [STEP2-LOD-CAP] improved LOD selection uses width budget
    const int selectedLevel = selectBestLevelForCurrentZoom (width);

    // Build visible groups (bounded size)
    buildVisibleGroupsForLevel (selectedLevel, width, scratchVisibleGroups, scratchVisibleFramesAgo);

    const size_t n = scratchVisibleGroups.size();
    if (n < 2)
        return;

    computeRepresentativeCurves (scratchVisibleGroups, scratchRepMomentaryDb, scratchRepShortTermDb);

    scratchPathRepM.clear();
    scratchPathRepS.clear();

    bool startedRepM = false, startedRepS = false;

    const auto colBandS = juce::Colours::cyan.withMultipliedAlpha (0.6f);
    const auto colBandM = juce::Colours::limegreen.withMultipliedAlpha (0.7f);
    const auto colLineS = juce::Colours::cyan.withMultipliedAlpha (0.9f);
    const auto colLineM = juce::Colours::limegreen;

    for (size_t i = 0; i < n; ++i)
    {
        const float x = w - (float) scratchVisibleFramesAgo[i] * (float) zoomX;
        if (x < -10.0f)
            continue;

        const auto& gGroup = scratchVisibleGroups[i];

        const float yMM = dbToY (gGroup.momentaryMaxDb, h);
        const float yMm = dbToY (gGroup.momentaryMinDb, h);
        const float ySM = dbToY (gGroup.shortTermMaxDb, h);
        const float ySm = dbToY (gGroup.shortTermMinDb, h);

        const float yRepM = dbToY (scratchRepMomentaryDb[i], h);
        const float yRepS = dbToY (scratchRepShortTermDb[i], h);

        // Bands
        if (showBands && selectedLevel > 0)
        {
            const float rangeMM = gGroup.momentaryMaxDb - gGroup.momentaryMinDb;
            const float rangeSM = gGroup.shortTermMaxDb - gGroup.shortTermMinDb;

            const float bandRangeThresholdDb = 3.0f;

            if (rangeSM >= bandRangeThresholdDb)
            {
                g.setColour (colBandS);
                g.drawLine (x, ySM, x, ySm, 1.0f);
            }

            if (rangeMM >= bandRangeThresholdDb)
            {
                g.setColour (colBandM);
                g.drawLine (x, yMM, x, yMm, 1.2f);
            }
        }

        // Lines
        if (showLines)
        {
            if (! startedRepS) { scratchPathRepS.startNewSubPath (x, yRepS); startedRepS = true; }
            else               { scratchPathRepS.lineTo         (x, yRepS); }

            if (! startedRepM) { scratchPathRepM.startNewSubPath (x, yRepM); startedRepM = true; }
            else               { scratchPathRepM.lineTo          (x, yRepM); }
        }
    }

    if (showLines)
    {
        g.setColour (colLineS);
        g.strokePath (scratchPathRepS, juce::PathStrokeType (1.5f));

        g.setColour (colLineM);
        g.strokePath (scratchPathRepM, juce::PathStrokeType (2.0f));
    }

    // Time ruler (unchanged)
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

    // Overlay info
    g.setColour (juce::Colours::white);
    g.setFont (14.0f);

    const auto& Lsel         = levels[(size_t) selectedLevel];
    const double spanSeconds = (double) Lsel.spanFrames / visualFrameRate;

    juce::String info = "Level: " + juce::String (selectedLevel) +
                        " (span " + juce::String (Lsel.spanFrames) + " frames, " +
                        juce::String (spanSeconds, 3) + " s)" +
                        " | ZoomX: " + juce::String (zoomX, 4) +
                        " | ZoomY: " + juce::String (zoomY, 2) +
                        " | Bands: " + juce::String (showBands ? "ON" : "OFF") +
                        " | Lines: " + juce::String (showLines ? "ON" : "OFF") +
                        " | CapPts: " + juce::String (getMaxDrawablePoints (width));

    g.drawText (info, 8, 8, (int) std::min (w - 16.0f, 620.0f), 20, juce::Justification::topLeft);
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
    if (event.mods.isShiftDown())
        showBands = ! showBands;
    else
        showLines = ! showLines;

    repaint();
}
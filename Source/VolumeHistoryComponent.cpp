// ============================================================================
// VolumeHistoryComponent.cpp
// ============================================================================
//
// [FIX-ACCORDION]     : Suppress pendingFrames offset for curve at coarse levels
// [FIX-RULER-JITTER]  : Suppress pendingFrames offset for ruler at coarse levels
//
// ============================================================================

#include "VolumeHistoryComponent.h"
#include "PluginProcessor.h"

#include <cmath>
#include <algorithm>
#include <limits>

//==============================================================================
// Helper: format seconds as HH:MM:SS
//==============================================================================

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

    // UI update rate
    startTimerHz (30);

    markStaticBackgroundDirty();
}

VolumeHistoryComponent::~VolumeHistoryComponent()
{
    stopTimer();
}

//==============================================================================
// History init
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

    tickStepIndex = -1;
}

//==============================================================================
// Timer
//==============================================================================

void VolumeHistoryComponent::timerCallback()
{
    const bool gotNewData = drainProcessorFifo(); // [STEP1-PERF]
    if (gotNewData)
        repaint();
}

//==============================================================================
// History update
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
// History access
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
// LOD selection  [STEP2-LOD-CAP] + [FIX-OSCILLATION]
//==============================================================================

int VolumeHistoryComponent::getMaxDrawablePoints (int widthPixels) const noexcept
{
    const int w = juce::jmax (1, widthPixels);
    const int target = (int) std::round ((double) w * 1.10);
    return juce::jlimit (256, 8192, target);
}

//==============================================================================
// [FIX-OSCILLATION]
//
// Level selection now uses totalVisibleFrames (without pendingFrames) to
// ensure stable level selection. pendingFrames oscillates 0 -> (span-1) -> 0,
// which could cause level selection to flip at boundary zoom values.
//==============================================================================

int VolumeHistoryComponent::selectBestLevelForCurrentZoom (int widthPixels) const noexcept
{
    const int maxPoints = getMaxDrawablePoints (widthPixels);

    const double safeZoomX = (zoomX > 1.0e-9 ? zoomX : 1.0e-9);
    const double overscanPixels = 10.0;
    const double totalVisibleFrames = ((double) widthPixels + overscanPixels) / safeZoomX;

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

        // [FIX-OSCILLATION]
        // Use totalVisibleFrames for prediction, NOT adjusted by pendingFrames.
        // This makes level selection stable regardless of pending frame count.
        const int predicted = (int) std::floor (totalVisibleFrames / (double) spanFrames) + 1;
        const int predictedClamped = juce::jlimit (1, available, predicted);

        if (predictedClamped <= maxPoints)
        {
            bestLevel = level;
            break;
        }
    }

    if (bestLevel >= 0)
        return bestLevel;

    for (int level = maxLevels - 1; level >= 0; --level)
        if (getAvailableGroups (level) >= 2)
            return level;

    return 0;
}

//==============================================================================
// [FIX-OSCILLATION]
//
// buildVisibleGroupsForLevel
//
// Root cause of oscillation: pendingFrames affects maxGroupsByX, which affects
// step calculation. When pendingFrames oscillates (0 -> span-1 -> 0), step
// can oscillate near threshold values, causing point spacing to oscillate.
//
// Fix: compute step based on TOTAL visible frames (zoom-dependent only).
// This makes step stable. Point positions can still use pendingFrames for
// smooth scrolling at fine levels, but the SPACING is always uniform.
//
// Additional fix: iterate from groupsAgo=0 with fixed step instead of 
// computing positions from outCount. This anchors the sampling grid to the
// right edge and prevents the entire grid from shifting when outCount changes.
//==============================================================================

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

    const int rawPendingFrames = getPendingFramesAtLevel (levelIndex);

    const double overscanPixels = 10.0;
    const double totalVisibleFrames = ((double) widthPixels + overscanPixels) / zoomX;
    if (totalVisibleFrames <= 0.0)
        return;

    // -------------------------------------------------------------------------
    // [FIX-OSCILLATION] STEP 1: Compute step based on total visible frames
    //
    // This makes step depend ONLY on zoom level, not on pendingFrames.
    // pendingFrames oscillates 0 -> (groupsPerGroup * prevSpanFrames - 1) -> 0,
    // which would cause step to oscillate if used in this calculation.
    // -------------------------------------------------------------------------
    const int maxGroupsForStep = (int) std::floor (totalVisibleFrames / (double) spanFrames) + 1;
    const int maxDrawablePoints = getMaxDrawablePoints (widthPixels);

    const int step = (maxGroupsForStep > maxDrawablePoints)
                        ? ((maxGroupsForStep + maxDrawablePoints - 1) / maxDrawablePoints)
                        : 1;

    // -------------------------------------------------------------------------
    // [FIX-OSCILLATION] STEP 2: Decide pendingFrames usage for positions
    //
    // Fine levels (< coarseLevelStartForPixelAdvance): use pendingFrames for
    // smooth scrolling. The spacing is stable (step is stable), only the
    // absolute positions shift smoothly.
    //
    // Coarse levels (>= coarseLevelStartForPixelAdvance): suppress pendingFrames.
    // Updates happen in discrete group-sized steps (all points move together).
    // -------------------------------------------------------------------------
    const bool suppressPendingForPositions = (levelIndex >= coarseLevelStartForPixelAdvance);
    const int pendingForPositions = suppressPendingForPositions ? 0 : rawPendingFrames;

    // -------------------------------------------------------------------------
    // [FIX-OSCILLATION] STEP 3: Compute how many groups to use
    //
    // For clipping, we need to account for pendingFrames to know what fits.
    // But use the position-adjusted value, not the raw value.
    // -------------------------------------------------------------------------
    const double effectiveVisibleFrames = totalVisibleFrames - (double) pendingForPositions;
    const int maxGroupsByX = (effectiveVisibleFrames > 0.0)
                                ? (int) std::floor (effectiveVisibleFrames / (double) spanFrames) + 1
                                : 0;

    const int groupsToUse = juce::jlimit (0, availableGroups, juce::jmax (0, maxGroupsByX));
    if (groupsToUse < 2)
        return;

    // -------------------------------------------------------------------------
    // [FIX-OSCILLATION] STEP 4: Iterate from groupsAgo=0 with fixed step
    //
    // This anchors the sampling grid to the right edge (most recent data).
    // When new data arrives, old points shift left but spacing stays uniform.
    // -------------------------------------------------------------------------
    const int estimatedCount = (groupsToUse / step) + 2;
    outGroups.reserve ((size_t) estimatedCount);
    outFramesAgo.reserve ((size_t) estimatedCount);

    for (int groupsAgo = 0; groupsAgo < groupsToUse; groupsAgo += step)
    {
        // Aggregate over [groupsAgo, groupsAgo + step - 1], clamped to available
        const int endGroupsAgo = juce::jmin (groupsToUse - 1, groupsAgo + step - 1);

        FrameGroup agg;
        agg.momentaryMinDb =  std::numeric_limits<float>::infinity();
        agg.momentaryMaxDb = -std::numeric_limits<float>::infinity();
        agg.shortTermMinDb =  std::numeric_limits<float>::infinity();
        agg.shortTermMaxDb = -std::numeric_limits<float>::infinity();

        for (int ga = groupsAgo; ga <= endGroupsAgo; ++ga)
        {
            const FrameGroup gg = getGroupAgo (levelIndex, ga);

            agg.momentaryMinDb = std::min (agg.momentaryMinDb, gg.momentaryMinDb);
            agg.momentaryMaxDb = std::max (agg.momentaryMaxDb, gg.momentaryMaxDb);
            agg.shortTermMinDb = std::min (agg.shortTermMinDb, gg.shortTermMinDb);
            agg.shortTermMaxDb = std::max (agg.shortTermMaxDb, gg.shortTermMaxDb);
        }

        // Center of the aggregated range for position calculation
        const double centerGroupsAgo = 0.5 * (double) (groupsAgo + endGroupsAgo);
        const double framesAgoD = (double) pendingForPositions + centerGroupsAgo * (double) spanFrames;
        const int framesAgoI = (int) std::round (framesAgoD);

        outGroups.push_back (agg);
        outFramesAgo.push_back (framesAgoI);
    }

    // -------------------------------------------------------------------------
    // Reverse to chronological order (oldest first).
    // computeRepresentativeCurves expects chronological order for trend calc.
    // -------------------------------------------------------------------------
    std::reverse (outGroups.begin(), outGroups.end());
    std::reverse (outFramesAgo.begin(), outFramesAgo.end());
}

//==============================================================================
// Representative curves
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
        const auto& grp = groups[i];

        const float minM    = grp.momentaryMinDb;
        const float maxM    = grp.momentaryMaxDb;
        const float centerM = 0.5f * (minM + maxM);

        const float minS    = grp.shortTermMinDb;
        const float maxS    = grp.shortTermMaxDb;
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
// Geometry helpers
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

float VolumeHistoryComponent::quantizeXToPixelCenter (float x) const noexcept
{
    const int xPix = (int) std::floor (x + 0.5f);
    return (float) xPix + 0.5f;
}

//==============================================================================
// Cached background  [CACHE-STATIC]
//==============================================================================

void VolumeHistoryComponent::markStaticBackgroundDirty() noexcept
{
    staticBackgroundDirty = true;
}

void VolumeHistoryComponent::rebuildStaticBackgroundIfNeeded()
{
    const int w = getWidth();
    const int h = getHeight();

    if (w <= 0 || h <= 0)
        return;

    const bool sizeChanged  = (w != cachedBgW || h != cachedBgH);
    const bool zoomYChanged = (std::abs (zoomY - cachedBgZoomY) > 1.0e-12);

    if (! staticBackgroundDirty && ! sizeChanged && ! zoomYChanged && cachedStaticBackground.isValid())
        return;

    cachedBgW = w;
    cachedBgH = h;
    cachedBgZoomY = zoomY;
    staticBackgroundDirty = false;

    cachedStaticBackground = juce::Image (juce::Image::RGB, w, h, true);
    juce::Graphics gg (cachedStaticBackground);

    auto backgroundColour = juce::Colours::limegreen.darker (3.0f);
    gg.fillAll (backgroundColour);

    gg.setColour (juce::Colours::darkgrey.withMultipliedAlpha (0.4f));
    const int numLines = 4;
    for (int i = 0; i <= numLines; ++i)
    {
        const float db = maxDb - (baseDbRange / (float) numLines) * (float) i;
        const float y  = dbToY (db, (float) h);
        gg.drawHorizontalLine ((int) std::round (y), 0.0f, (float) w);
    }

    const float rulerBaseY = (float) h - 2.0f;
    gg.setColour (juce::Colours::darkgrey.withMultipliedAlpha (0.7f));
    gg.drawHorizontalLine ((int) std::round (rulerBaseY), 0.0f, (float) w);
}

//==============================================================================
// Line quality / render mode  [LINE-QUALITY]
//==============================================================================

bool VolumeHistoryComponent::isModifierForQualityToggle (const juce::ModifierKeys& mods) const noexcept
{
   #if JUCE_MAC
    return mods.isCommandDown();
   #else
    return mods.isCtrlDown();
   #endif
}

void VolumeHistoryComponent::cycleLineRenderMode() noexcept
{
    lineRenderMode = (lineRenderMode + 1) % 3;
}

bool VolumeHistoryComponent::shouldUsePolylineForLines (int selectedLevel) const noexcept
{
    if (lineRenderMode == 2) return true;   // Force Polyline
    if (lineRenderMode == 1) return false;  // Force Stroke
    return selectedLevel >= coarseLevelStartForPolyline; // Auto
}

//==============================================================================
// Polyline drawing  [POLYLINE-PEAK] + [PIXEL-ADVANCE-FIX] + [FIX-LINE-DROPOUT]
//
// Fixed version: 
//   1. Extended clipping margins to include off-screen points for proper 
//      edge connections
//   2. Ensure we don't lose boundary points due to step-related shifts
//==============================================================================

void VolumeHistoryComponent::buildPolylinePoints (const std::vector<int>& framesAgo,
                                                  const std::vector<float>& repDb,
                                                  float width,
                                                  float height,
                                                  std::vector<juce::Point<float>>& outPoints) const
{
    outPoints.clear();

    const size_t n = juce::jmin (framesAgo.size(), repDb.size());
    if (n < 2 || zoomX <= 0.0 || width <= 1.0f)
        return;

    const int wInt = (int) std::round (width);
    outPoints.reserve ((size_t) juce::jlimit (128, 4096, wInt + 64));

    // -------------------------------------------------------------------------
    // [FIX-LINE-DROPOUT]
    // Extended margins: include points slightly off-screen on both sides.
    // This ensures:
    //   1. Lines connect properly at screen edges
    //   2. Step-related shifts don't push boundary points out of bounds
    //   3. We always have enough points to draw a continuous line
    //
    // The extra off-screen points will be clipped by JUCE's rendering,
    // but the line segments connecting them to on-screen points will draw.
    // -------------------------------------------------------------------------
    const int leftClipMargin  = 50;   // pixels off left edge to include
    const int rightClipMargin = 10;   // pixels off right edge to include

    int   currentXPix = std::numeric_limits<int>::min();
    float colYMin = 0.0f;
    float colYMax = 0.0f;
    bool  haveCol = false;

    float prevY = 0.0f;
    bool  havePrev = false;

    auto emitColumn = [&]()
    {
        if (! haveCol)
            return;

        const float span = std::abs (colYMax - colYMin);

        float yRep = 0.5f * (colYMin + colYMax);

        if (havePrev)
        {
            // Preserve peaks when the column span is meaningful in pixels.
            constexpr float spanThresholdPx = 1.0f;

            if (span >= spanThresholdPx)
            {
                const float dMin = std::abs (colYMin - prevY);
                const float dMax = std::abs (colYMax - prevY);
                yRep = (dMin >= dMax ? colYMin : colYMax);
            }
            else
            {
                // For tiny spans, stay stable to reduce shimmer.
                yRep = juce::jlimit (colYMin, colYMax, prevY);
            }
        }

        const float xAligned = (float) currentXPix + 0.5f;
        outPoints.emplace_back (xAligned, yRep);

        prevY = yRep;
        havePrev = true;
    };

    for (size_t i = 0; i < n; ++i)
    {
        const float xRaw = width - (float) framesAgo[i] * (float) zoomX;

        // [FIX-LINE-DROPOUT] Extended left margin for early exit
        if (xRaw < (float) -leftClipMargin)
            continue;

        const int xPix = (int) std::floor (xRaw + 0.5f);

        // [FIX-LINE-DROPOUT] Extended clipping bounds
        if (xPix < -leftClipMargin || xPix > wInt + rightClipMargin)
            continue;

        const float y = dbToY (repDb[i], height);

        if (! haveCol)
        {
            haveCol = true;
            currentXPix = xPix;
            colYMin = y;
            colYMax = y;
            continue;
        }

        if (xPix != currentXPix)
        {
            emitColumn();

            currentXPix = xPix;
            colYMin = y;
            colYMax = y;
        }
        else
        {
            colYMin = std::min (colYMin, y);
            colYMax = std::max (colYMax, y);
        }
    }

    emitColumn();
}

void VolumeHistoryComponent::drawPolyline (juce::Graphics& g,
                                          const std::vector<juce::Point<float>>& pts,
                                          float thickness) const
{
    if (pts.size() < 2)
        return;

    for (size_t i = 1; i < pts.size(); ++i)
    {
        const auto& a = pts[i - 1];
        const auto& b = pts[i];
        g.drawLine (a.x, a.y, b.x, b.y, thickness);
    }
}

//==============================================================================
// Ruler tickStep hysteresis  [RULER-HYST-FIX]
//==============================================================================

double VolumeHistoryComponent::getTickStepSecondsWithHysteresis (int widthPixels) noexcept
{
    const double pixelsPerSecond = zoomX * visualFrameRate;
    if (pixelsPerSecond <= 1.0e-12 || widthPixels <= 0)
        return 60.0;

    // Keep spacing in a band; add margins so we don't thrash near boundaries.
    const double minSpacingPx = (double) widthPixels / 20.0; // target <= 20 labels
    const double maxSpacingPx = (double) widthPixels / 10.0; // target >= 10 labels

    const double minSwitchPx = minSpacingPx * 0.95; // go coarser only if clearly too dense
    const double maxSwitchPx = maxSpacingPx * 1.05; // go finer only if clearly too sparse

    static constexpr double tickSteps[] = {
        0.5, 1.0, 2.0, 5.0,
        10.0, 15.0, 30.0,
        60.0, 120.0, 300.0, 600.0, 1200.0, 3600.0
    };
    static constexpr int numSteps = (int) (sizeof (tickSteps) / sizeof (tickSteps[0]));

    // Initialize once; DO NOT reset on every zoom tick, otherwise hysteresis can't work.
    if (tickStepIndex < 0 || tickStepIndex >= numSteps)
    {
        int idx = numSteps - 1;
        for (int i = 0; i < numSteps; ++i)
        {
            const double spacing = tickSteps[i] * pixelsPerSecond;
            if (spacing >= minSpacingPx)
            {
                idx = i;
                break;
            }
        }
        tickStepIndex = idx;
        return tickSteps[tickStepIndex];
    }

    // Adjust with guard (no infinite loops)
    for (int guard = 0; guard < numSteps + 2; ++guard)
    {
        const double spacing = tickSteps[tickStepIndex] * pixelsPerSecond;

        if (spacing < minSwitchPx && tickStepIndex < numSteps - 1)
        {
            ++tickStepIndex; // coarser
            continue;
        }

        if (spacing > maxSwitchPx && tickStepIndex > 0)
        {
            --tickStepIndex; // finer
            continue;
        }

        break;
    }

    tickStepIndex = juce::jlimit (0, numSteps - 1, tickStepIndex);
    return tickSteps[tickStepIndex];
}

//==============================================================================
// Component lifecycle
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

    // Reserve scratch buffers
    {
        const int w = juce::jmax (1, getWidth());
        const size_t reserveCount = (size_t) juce::jlimit (256, 8192, w + 512);

        scratchVisibleGroups.reserve (reserveCount);
        scratchVisibleFramesAgo.reserve (reserveCount);
        scratchRepMomentaryDb.reserve (reserveCount);
        scratchRepShortTermDb.reserve (reserveCount);

        scratchPolylinePtsM.reserve ((size_t) juce::jlimit (256, 4096, w + 64));
        scratchPolylinePtsS.reserve ((size_t) juce::jlimit (256, 4096, w + 64));
    }

    // Do NOT reset tickStepIndex here; hysteresis should smoothly adapt.
    markStaticBackgroundDirty();
}

//==============================================================================
// DRAW
//==============================================================================

void VolumeHistoryComponent::paint (juce::Graphics& g)
{
    auto bounds = getLocalBounds().toFloat();
    const int width  = (int) bounds.getWidth();
    const int height = (int) bounds.getHeight();

    if (width <= 1 || height <= 0)
        return;

    rebuildStaticBackgroundIfNeeded();
    if (cachedStaticBackground.isValid())
        g.drawImageAt (cachedStaticBackground, 0, 0);

    const juce::int64 totalFrames   = getTotalFramesL0();
    const juce::int64 availableRaw  = std::min<juce::int64> ((juce::int64) rawCapacityFrames, totalFrames);

    const int selectedLevel = selectBestLevelForCurrentZoom (width);

    // -------------------------------------------------------------------------
    // [FIX-RULER-JITTER]
    // For coarse levels, suppress pendingFramesOffset to prevent the ruler
    // from jumping back-and-forth as the pendingFrames value oscillates.
    // -------------------------------------------------------------------------
    const bool suppressPendingOffset = (selectedLevel >= coarseLevelStartForPixelAdvance);
    const int rawPendingFramesOffset = getPendingFramesAtLevel (selectedLevel);
    const int pendingFramesOffset = suppressPendingOffset ? 0 : rawPendingFramesOffset;

    buildVisibleGroupsForLevel (selectedLevel, width, scratchVisibleGroups, scratchVisibleFramesAgo);

    const size_t n = scratchVisibleGroups.size();
    const float w = bounds.getWidth();
    const float h = bounds.getHeight();

    const bool usePolylineLines = (showLines && shouldUsePolylineForLines (selectedLevel));
    const bool pixelAdvanceActive = (usePolylineLines && selectedLevel >= coarseLevelStartForPixelAdvance);

    if (n >= 2)
    {
        computeRepresentativeCurves (scratchVisibleGroups, scratchRepMomentaryDb, scratchRepShortTermDb);

        scratchPathBandM.clear();
        scratchPathBandS.clear();

        if (showLines && ! usePolylineLines)
        {
            scratchPathRepM.clear();
            scratchPathRepS.clear();
        }

        const float bandRangeThresholdDb = 3.0f;

        bool startedRepM = false, startedRepS = false;

        for (size_t i = 0; i < n; ++i)
        {
            float x = w - (float) scratchVisibleFramesAgo[i] * (float) zoomX;
            if (x < -10.0f)
                continue;

            // [PIXEL-ADVANCE-FIX] In coarse polyline mode, quantize X to pixel center
            if (pixelAdvanceActive)
                x = quantizeXToPixelCenter (x);

            const auto& grp = scratchVisibleGroups[i];

            const float yMM = dbToY (grp.momentaryMaxDb, h);
            const float yMm = dbToY (grp.momentaryMinDb, h);
            const float ySM = dbToY (grp.shortTermMaxDb, h);
            const float ySm = dbToY (grp.shortTermMinDb, h);

            if (showBands && selectedLevel > 0)
            {
                const float rangeMM = grp.momentaryMaxDb - grp.momentaryMinDb;
                const float rangeSM = grp.shortTermMaxDb - grp.shortTermMinDb;

                if (rangeSM >= bandRangeThresholdDb)
                {
                    scratchPathBandS.startNewSubPath (x, ySM);
                    scratchPathBandS.lineTo (x, ySm);
                }

                if (rangeMM >= bandRangeThresholdDb)
                {
                    scratchPathBandM.startNewSubPath (x, yMM);
                    scratchPathBandM.lineTo (x, yMm);
                }
            }

            if (showLines && ! usePolylineLines)
            {
                const float yRepM = dbToY (scratchRepMomentaryDb[i], h);
                const float yRepS = dbToY (scratchRepShortTermDb[i], h);

                if (! startedRepS) { scratchPathRepS.startNewSubPath (x, yRepS); startedRepS = true; }
                else               { scratchPathRepS.lineTo         (x, yRepS); }

                if (! startedRepM) { scratchPathRepM.startNewSubPath (x, yRepM); startedRepM = true; }
                else               { scratchPathRepM.lineTo          (x, yRepM); }
            }
        }

        // Bands
        if (showBands && selectedLevel > 0)
        {
            g.setColour (juce::Colours::cyan.withMultipliedAlpha (0.6f));
            g.strokePath (scratchPathBandS, juce::PathStrokeType (1.0f));

            g.setColour (juce::Colours::limegreen.withMultipliedAlpha (0.7f));
            g.strokePath (scratchPathBandM, juce::PathStrokeType (1.2f));
        }

        // Lines
        if (showLines)
        {
            if (usePolylineLines)
            {
                scratchPolylinePtsM.clear();
                scratchPolylinePtsS.clear();

                buildPolylinePoints (scratchVisibleFramesAgo, scratchRepShortTermDb, w, h, scratchPolylinePtsS);
                buildPolylinePoints (scratchVisibleFramesAgo, scratchRepMomentaryDb, w, h, scratchPolylinePtsM);

                g.setColour (juce::Colours::cyan.withMultipliedAlpha (0.9f));
                drawPolyline (g, scratchPolylinePtsS, 1.0f);

                g.setColour (juce::Colours::limegreen);
                drawPolyline (g, scratchPolylinePtsM, 1.0f);
            }
            else
            {
                g.setColour (juce::Colours::cyan.withMultipliedAlpha (0.9f));
                g.strokePath (scratchPathRepS, juce::PathStrokeType (1.5f));

                g.setColour (juce::Colours::limegreen);
                g.strokePath (scratchPathRepM, juce::PathStrokeType (2.0f));
            }
        }
    }

    //==========================================================================
    // [RULER-FRAMES] + [FIX-RULER-JITTER]
    // Ruler ticks computed in integer frames to eliminate flip-flop and jitter.
    // For coarse levels, pendingFramesOffset is suppressed (set to 0 above).
    //==========================================================================
    const float rulerHeight   = 16.0f;
    const float rulerBaseY    = (float) height - 2.0f;
    const float tickTopY      = rulerBaseY - 6.0f;
    const float textTopY      = rulerBaseY - 14.0f;

    if (zoomX > 0.0 && visualFrameRate > 0.0 && totalFrames > 0)
    {
        // Visible frames based on x-mapping (no fake stretching)
        const double framesByWidth = (double) width / (zoomX > 1.0e-12 ? zoomX : 1.0e-12);
        const double effectiveFramesByWidth = juce::jmax (0.0, framesByWidth - (double) pendingFramesOffset);

        const juce::int64 visibleFrames = (juce::int64) std::floor (juce::jmin ((double) availableRaw, effectiveFramesByWidth));
        if (visibleFrames > 1)
        {
            const juce::int64 earliestAvailableFrame = totalFrames - availableRaw;
            const juce::int64 leftFrame = juce::jmax (earliestAvailableFrame, totalFrames - visibleFrames);

            // Tick step (seconds) from zoom with hysteresis
            const double tickStepSec = getTickStepSecondsWithHysteresis (width);

            // Convert to integer frames at 60 Hz (stable)
            const juce::int64 tickStepFrames = (juce::int64) juce::jmax (1.0, std::round (tickStepSec * visualFrameRate));

            // Anchor ticks to "now" (right edge) in integer frames
            const juce::int64 lastTickFrame = (totalFrames / tickStepFrames) * tickStepFrames;

            g.setColour (juce::Colours::white);
            g.setFont (10.0f);

            for (juce::int64 tickFrame = lastTickFrame; tickFrame >= leftFrame; tickFrame -= tickStepFrames)
            {
                const juce::int64 framesAgo = totalFrames - tickFrame;

                // [FIX-RULER-JITTER] pendingFramesOffset is 0 for coarse levels
                const double framesAgoWithOffset = (double) framesAgo + (double) pendingFramesOffset;

                float x = (float) ((double) width - framesAgoWithOffset * zoomX);

                // In coarse polyline mode, quantize ticks to match curve pixel advance
                if (pixelAdvanceActive)
                    x = quantizeXToPixelCenter (x);

                if (x < -2.0f)
                    break;

                if (x > (float) width + 2.0f)
                    continue;

                g.drawLine (x, tickTopY, x, rulerBaseY, 1.0f);

                const double t = (double) tickFrame / visualFrameRate;
                const float textWidth = 52.0f;

                g.drawText (formatTimeHMS (t),
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

    const char* modeName = "Auto";
    if (lineRenderMode == 1) modeName = "ForceStroke";
    if (lineRenderMode == 2) modeName = "ForcePolyline";

    const double tickStepNow = getTickStepSecondsWithHysteresis (width);

    juce::String info = "Level: " + juce::String (selectedLevel) +
                        " (span " + juce::String (Lsel.spanFrames) + " frames, " +
                        juce::String (spanSeconds, 3) + " s)" +
                        " | ZoomX: " + juce::String (zoomX, 4) +
                        " | ZoomY: " + juce::String (zoomY, 2) +
                        " | Tick: " + juce::String (tickStepNow, 3) + "s" +
                        " | Bands: " + juce::String (showBands ? "ON" : "OFF") +
                        " | Lines: " + juce::String (showLines ? "ON" : "OFF") +
                        " | LineMode: " + juce::String (modeName);

    g.drawText (info, 8, 8, (int) std::min (bounds.getWidth() - 16.0f, 980.0f), 20, juce::Justification::topLeft);
}

//==============================================================================
// Zoom
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

    // Do NOT reset tickStepIndex here; hysteresis should do its job.
}

void VolumeHistoryComponent::applyVerticalZoom (float wheelDelta)
{
    if (wheelDelta == 0.0f)
        return;

    const double zoomBase   = 1.1;
    const double zoomFactor = std::pow (zoomBase, (double) wheelDelta);

    zoomY *= zoomFactor;
    zoomY  = juce::jlimit (minZoomY, maxZoomY, zoomY);

    markStaticBackgroundDirty();
}

//==============================================================================
// Mouse
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
    // Cmd-click (mac) / Ctrl-click (win): cycle line render mode
    if (isModifierForQualityToggle (event.mods))
    {
        cycleLineRenderMode();
        repaint();
        return;
    }

    // Shift-click toggles bands; plain click toggles lines.
    if (event.mods.isShiftDown())
        showBands = ! showBands;
    else
        showLines = ! showLines;

    repaint();
}
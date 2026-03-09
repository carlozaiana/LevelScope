#include "VolumeHistoryComponent.h"
#include "PluginProcessor.h"

//==============================================================================
// VolumeHistoryComponent_History.cpp
// - FIFO drain / timeline writes
// - overwrite-safe LOD pyramid
// - visible group building
// - persistence bootstrap from processor
//==============================================================================
//==============================================================================
// [EDIT-BLOCKS]
//   - [BEGIN VHC-HIST-TIMER-AND-DRAIN]   ... [END VHC-HIST-TIMER-AND-DRAIN]
//   - [BEGIN VHC-HIST-PUSH-FRAME]        ... [END VHC-HIST-PUSH-FRAME]
//   - [BEGIN VHC-HIST-RING-HELPERS]      ... [END VHC-HIST-RING-HELPERS]
//   - [BEGIN VHC-HIST-INIT-RESET]        ... [END VHC-HIST-INIT-RESET]
//   - [BEGIN VHC-HIST-ACCESS]            ... [END VHC-HIST-ACCESS]
//   - [BEGIN VHC-HIST-LOD]               ... [END VHC-HIST-LOD]
//   - [BEGIN VHC-HIST-VISIBLE-BUILDER]   ... [END VHC-HIST-VISIBLE-BUILDER]
//   - [BEGIN VHC-HIST-BOOTSTRAP]         ... [END VHC-HIST-BOOTSTRAP]
//   - [BEGIN VHC-HIST-HELPERS]           ... [END VHC-HIST-HELPERS]
//==============================================================================

//==============================================================================
// Timer
//==============================================================================

// [BEGIN VHC-HIST-TIMER-AND-DRAIN]
void VolumeHistoryComponent::timerCallback()
{
    const bool gotNewData = drainProcessorFifo();

    // [ROLLING-LRA] if window selection changed, start rebuild; then do small chunks per tick
    startRollingRebuildIfWindowChanged();
    if (rollingRebuildInProgress)
        rebuildRollingLraStep (200); // process up to 200 seconds per timer tick

    // [BEGIN UI-METERS-IDLE-DECAY-CALL]
    // Update cached right-strip meter display values (also decays when host stops calling audio).
    updateRightStripMetersFromTimer();
    // [END UI-METERS-IDLE-DECAY-CALL]

    // Repaint if we got new history data OR rolling is rebuilding OR meters need animation.
    // Meters need repaint even when gotNewData is false (e.g. host stopped calling processBlock()).
    const bool metersNeedRepaint = updateRightStripMetersFromTimer();
    // [BEGIN UI-METERS-IDLE-DECAY-CONDITIONAL-REPAINT]
    if (gotNewData || rollingRebuildInProgress)
    {
        // Full repaint when curves/rolling rebuild changed
        repaint();
    }
    else if (metersNeedRepaint)
    {
        // Only meters changed/decayed -> repaint only the right strip (cheap)
        repaint (getDbRulerArea());
    }
    // [END UI-METERS-IDLE-DECAY-CONDITIONAL-REPAINT]
}

bool VolumeHistoryComponent::drainProcessorFifo()
{
    // [TIMEBASE-PLAYHEAD] cache samples-per-visual-frame (60 Hz) from processor
    if (frameSamples <= 0)
    frameSamples = processor.getFrameSamples();
    constexpr int chunkSize = 512;
    float momentaryValues [chunkSize];
    float shortTermValues [chunkSize];
    float gateValues [chunkSize]; // [LRAG]
    juce::int64 frameIndex60Hz [chunkSize];
    int playingFlags [chunkSize];

    bool readAny = false;

    for (;;)
    {
        const int numRead = processor.readLoudnessFromFifoEx (momentaryValues,
                                                             shortTermValues,
                                                             gateValues,
                                                             frameIndex60Hz,
                                                             playingFlags,
                                                             chunkSize);
        if (numRead <= 0)
            break;

        readAny = true;

        for (int i = 0; i < numRead; ++i)
            pushFrameToHistory (momentaryValues[i], shortTermValues[i], gateValues[i],
                                frameIndex60Hz[i], playingFlags[i]);
    }

    return readAny;
}
// [END VHC-HIST-TIMER-AND-DRAIN]

//==============================================================================
// Timeline write path
//==============================================================================

// [BEGIN VHC-HIST-PUSH-FRAME]
void VolumeHistoryComponent::pushFrameToHistory (float momentaryVal,
                                                 float shortTermVal,
                                                 float gateLufs,
                                                 juce::int64 frameIndex60Hz,
                                                 int isPlaying)
{
    // Only write history while transport plays (your requested behavior).
    // If the host doesn't provide play state reliably in some cases, you can relax this later.
    if (isPlaying == 0)
        return;

    if (frameSamples <= 0)
        return;

    const juce::int64 frameIndex = frameIndex60Hz;

    // [TIMEBASE-PLAYHEAD] actual DAW playhead position (can go backwards)
    playheadFrameIndex = frameIndex;
    havePlayheadFrameIndex = true;

    // [TIMEBASE-PLAYHEAD] Keep a monotonic "furthest written" cursor.
    // This prevents loops/rewinds from truncating the available range and making
    // already-written data after the loop appear "lost".
    if (! haveNowFrameIndex)
    {
        nowFrameIndex = frameIndex;
        haveNowFrameIndex = true;
    }
    else
    {
        nowFrameIndex = juce::jmax<juce::int64> (nowFrameIndex, frameIndex);
    }

    // Phase 2: processor FIFO already delivers LUFS (dB), so do NOT convert.
    const float lufsM = momentaryVal;
    const float lufsS = shortTermVal;
    const float lufsG = gateLufs;

    FrameGroup fg;
    fg.momentaryMinDb = lufsM;
    fg.momentaryMaxDb = lufsM;
    fg.shortTermMinDb = lufsS;
    fg.shortTermMaxDb = lufsS;

    // [LRAG]
    fg.gateMinDb = lufsG;
    fg.gateMaxDb = lufsG;

    // L0 overwrite at absolute frame index
    writeGroupAbs (0, frameIndex, fg);

    // [ROLLING-LRA] sample once per second from the 60 Hz stream
    const juce::int64 absSecondIndex = floorDivInt64 (frameIndex, 60);

    if (absSecondIndex != lastSecondPushed)
    {
        lastSecondPushed = absSecondIndex;
        pushSecondSample (absSecondIndex, lufsS, lufsG);
    }

    // Recompute parent groups up the pyramid (overwrite-safe)
    for (int level = 1; level < maxLevels; ++level)
    {
        const auto& L = levels[(size_t) level];
        const juce::int64 absGroupIndex = frameIndex / (juce::int64) L.spanFrames;
        recomputeGroupAbsFromChildren (level, absGroupIndex);
    }
}
// [END VHC-HIST-PUSH-FRAME]

//==============================================================================
// Overwrite-safe ring helpers
//==============================================================================

// [BEGIN VHC-HIST-RING-HELPERS]
void VolumeHistoryComponent::writeGroupAbs (int levelIndex,
                                            juce::int64 absGroupIndex,
                                            const FrameGroup& group)
{
    if (levelIndex < 0 || levelIndex >= maxLevels)
        return;

    auto& L = levels[(size_t) levelIndex];
    if (L.capacity <= 0 || absGroupIndex < 0)
        return;

    juce::int64 m = absGroupIndex % (juce::int64) L.capacity;
    if (m < 0) m += (juce::int64) L.capacity;
    const int slot = (int) m;
    L.groups[(size_t) slot] = group;
    L.absGroupIndexTag[(size_t) slot] = absGroupIndex;
}

bool VolumeHistoryComponent::readGroupAbs (int levelIndex,
                                           juce::int64 absGroupIndex,
                                           FrameGroup& out) const noexcept
{
    out = FrameGroup { minDb, minDb, minDb, minDb };

    if (levelIndex < 0 || levelIndex >= maxLevels)
        return false;

    const auto& L = levels[(size_t) levelIndex];
    if (L.capacity <= 0 || absGroupIndex < 0)
        return false;

    juce::int64 m = absGroupIndex % (juce::int64) L.capacity;
    if (m < 0) m += (juce::int64) L.capacity;
    const int slot = (int) m;
    if (L.absGroupIndexTag[(size_t) slot] != absGroupIndex)
        return false;

    out = L.groups[(size_t) slot];
    return true;
}

void VolumeHistoryComponent::recomputeGroupAbsFromChildren (int levelIndex,
                                                           juce::int64 absGroupIndex)
{
    if (levelIndex <= 0 || levelIndex >= maxLevels)
        return;

    auto& L = levels[(size_t) levelIndex];
    if (L.capacity <= 0 || absGroupIndex < 0)
        return;

    FrameGroup agg;
    agg.momentaryMinDb =  std::numeric_limits<float>::infinity();
    agg.momentaryMaxDb = -std::numeric_limits<float>::infinity();
    agg.shortTermMinDb =  std::numeric_limits<float>::infinity();
    agg.shortTermMaxDb = -std::numeric_limits<float>::infinity();

    // [LRAG]
    agg.gateMinDb =  std::numeric_limits<float>::infinity();
    agg.gateMaxDb = -std::numeric_limits<float>::infinity();

    bool any = false;

    // Children live in levelIndex - 1
    const int childLevel = levelIndex - 1;
    const juce::int64 childStart = absGroupIndex * (juce::int64) groupsPerLevel;

    for (int k = 0; k < groupsPerLevel; ++k)
    {
        FrameGroup child;
        if (! readGroupAbs (childLevel, childStart + (juce::int64) k, child))
            continue;

        any = true;
        agg.momentaryMinDb = std::min (agg.momentaryMinDb, child.momentaryMinDb);
        agg.momentaryMaxDb = std::max (agg.momentaryMaxDb, child.momentaryMaxDb);
        agg.shortTermMinDb = std::min (agg.shortTermMinDb, child.shortTermMinDb);
        agg.shortTermMaxDb = std::max (agg.shortTermMaxDb, child.shortTermMaxDb);
        agg.gateMinDb = std::min (agg.gateMinDb, child.gateMinDb);
        agg.gateMaxDb = std::max (agg.gateMaxDb, child.gateMaxDb);
    }

    if (! any)
        return; // keep old value/tag (or empty)

    writeGroupAbs (levelIndex, absGroupIndex, agg);
}
// [END VHC-HIST-RING-HELPERS]

//==============================================================================
// History init/reset
//==============================================================================

// [BEGIN VHC-HIST-INIT-RESET]
void VolumeHistoryComponent::initialiseHistoryLevels()
{
    rawCapacityFrames = (int) std::ceil (historyLengthSeconds * visualFrameRate);
    if (rawCapacityFrames < 1)
        rawCapacityFrames = 1;

    // Level 0 (RAW)
    {
        auto& L0 = levels[0];
        L0.levelIndex     = 0;
        L0.spanFrames     = 1;
        L0.capacity       = rawCapacityFrames;
        L0.groups.assign ((size_t) rawCapacityFrames,
                          FrameGroup { minDb, minDb, minDb, minDb });
        L0.absGroupIndexTag.assign ((size_t) rawCapacityFrames, (juce::int64) -1);
    }

    // Derived levels
    int prevSpanFrames = levels[0].spanFrames;

    for (int level = 1; level < maxLevels; ++level)
    {
        auto& L = levels[(size_t) level];

        L.levelIndex     = level;
        L.spanFrames     = prevSpanFrames * groupsPerLevel;
        prevSpanFrames   = L.spanFrames;

        int capacity = (int) std::ceil ((double) rawCapacityFrames / (double) L.spanFrames);
        if (capacity < 1)
            capacity = 1;

        L.capacity     = capacity;
        L.groups.assign ((size_t) capacity,
                         FrameGroup { minDb, minDb, minDb, minDb });
        L.absGroupIndexTag.assign ((size_t) capacity, (juce::int64) -1);

    }
}

void VolumeHistoryComponent::resetHistoryLevels()
{
    for (auto& L : levels)
    {
        for (auto& g : L.groups)
        {
            g.momentaryMinDb = minDb;
            g.momentaryMaxDb = minDb;
            g.shortTermMinDb = minDb;
            g.shortTermMaxDb = minDb;
        }
        std::fill (L.absGroupIndexTag.begin(), L.absGroupIndexTag.end(), (juce::int64) -1);
    }
    haveNowFrameIndex = false;
    nowFrameIndex = 0;
    havePlayheadFrameIndex = false;
    playheadFrameIndex = 0;
    tickStepIndex = -1;
}
// [END VHC-HIST-INIT-RESET]

//==============================================================================
// History access
//==============================================================================

// [BEGIN VHC-HIST-ACCESS]
int VolumeHistoryComponent::getAvailableGroups (int levelIndex) const noexcept
{
    if (levelIndex < 0 || levelIndex >= maxLevels || ! haveNowFrameIndex)
        return 0;

    const auto& L = levels[(size_t) levelIndex];
    if (L.capacity <= 0)
        return 0;

    const juce::int64 nowGroupIndex = nowFrameIndex / (juce::int64) L.spanFrames;
    const juce::int64 available = juce::jmin<juce::int64> ((juce::int64) L.capacity, nowGroupIndex + 1);
    return (int) available;
}

VolumeHistoryComponent::FrameGroup
VolumeHistoryComponent::getGroupAgo (int levelIndex, int groupsAgo) const noexcept
{
    FrameGroup out { minDb, minDb, minDb, minDb };

    if (levelIndex < 0 || levelIndex >= maxLevels || ! haveNowFrameIndex)
        return out;

    const auto& L = levels[(size_t) levelIndex];
    if (L.capacity <= 0 || groupsAgo < 0)
        return out;

    const juce::int64 nowGroupIndex = nowFrameIndex / (juce::int64) L.spanFrames;
    const juce::int64 absIndex = nowGroupIndex - (juce::int64) groupsAgo;
    if (absIndex < 0)
        return out;

    const int slot = (int) (absIndex % (juce::int64) L.capacity);
    if (L.absGroupIndexTag[(size_t) slot] != absIndex)
        return out;

    return L.groups[(size_t) slot];
}

juce::int64 VolumeHistoryComponent::getTotalFramesL0() const noexcept
{
    return (haveNowFrameIndex ? nowFrameIndex : 0);
}
// [END VHC-HIST-ACCESS]

//==============================================================================
// LOD selection  [STEP2-LOD-CAP]
//==============================================================================

// [BEGIN VHC-HIST-LOD]
int VolumeHistoryComponent::getMaxDrawablePoints (int widthPixels) const noexcept
{
    const int w = juce::jmax (1, widthPixels);
    const int target = (int) std::round ((double) w * 1.10);
    return juce::jlimit (256, 8192, target);
}

int VolumeHistoryComponent::selectBestLevelForCurrentZoom (int widthPixels) const noexcept
{
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

        // Predict how many groups fit horizontally at this span.
        // +2 as a safety/overscan margin.
        const int predicted = (int) std::floor (maxFramesVisible / (double) spanFrames) + 2;
        const int predictedClamped = juce::jlimit (1, available, predicted);

        if (predictedClamped >= 2 && predictedClamped <= maxPoints)
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
// [END VHC-HIST-LOD]

//==============================================================================
// Visible groups builder  [TIMEBASE-FIX]  [DECIMATOR-GRID-ANCHOR]
//
// Key idea:
// - We still cap point count with a step.
// - But chunk boundaries are ANCHORED to an absolute grid (multiples of step),
//   not to the moving left edge. This prevents "accordion" for BOTH renderers.
// - Timestamp uses chunk END (no center/rounding drift).
//==============================================================================

// [BEGIN VHC-HIST-VISIBLE-BUILDER]
void VolumeHistoryComponent::buildVisibleGroupsForLevel (int levelIndex,
                                                         int widthPixels,
                                                         std::vector<FrameGroup>& outGroups,
                                                         std::vector<juce::int64>& outEndFrameIndex) const
{
    outGroups.clear();
    outEndFrameIndex.clear();

    if (levelIndex < 0 || levelIndex >= maxLevels)
        return;

    const auto& L = levels[(size_t) levelIndex];

    const int availableGroups = getAvailableGroups (levelIndex);
    if (availableGroups < 2 || widthPixels <= 0 || zoomX <= 0.0)
        return;

    const int spanFrames = L.spanFrames;
    if (spanFrames <= 0)
        return;

    // [TIMEBASE-PLAYHEAD] No playhead time yet -> nothing to show
    if (! haveNowFrameIndex)
        return;

    const juce::int64 totalFramesNow = (juce::int64) std::floor (viewRightFrame); // [VIEW-NAV]

    // Visible frames (overscanned a bit)
    const double overscanPixels = 10.0;
    const double maxFramesVisibleD = ((double) widthPixels + overscanPixels) / zoomX;
    const juce::int64 maxFramesVisible = (juce::int64) std::ceil (juce::jmax (0.0, maxFramesVisibleD));

    // Visible time window in absolute frames (Level 0 frame units)
    const juce::int64 leftFrame = totalFramesNow - maxFramesVisible;

    // [TIMEBASE-PLAYHEAD] Absolute group index range on the DAW timeline
    // [VIEW-NAV] Separate "data availability" from "view right edge".
    const juce::int64 latestWrittenAbsGroup = nowFrameIndex / (juce::int64) spanFrames;

    const juce::int64 viewRightFrameI = (juce::int64) std::floor (viewRightFrame);
    const juce::int64 latestViewAbsGroup = viewRightFrameI / (juce::int64) spanFrames;

    // What we will actually draw up to (cannot exceed what is written).
    const juce::int64 latestAbsGroup = juce::jmin (latestWrittenAbsGroup, latestViewAbsGroup);

    // Earliest group that can still be stored in our ring depends on latest WRITTEN, not the view.
    const juce::int64 earliestAbsGroup =
        latestWrittenAbsGroup - (juce::int64) (L.capacity - 1);

    if (latestAbsGroup < earliestAbsGroup)
        return;

    // Groups whose END frame is >= leftFrame are potentially visible.
    // Group absIndex has endFrame = (absIndex + 1) * spanFrames.
    // We need: (abs+1)*spanFrames >= leftFrame  =>  abs >= ceil(leftFrame/spanFrames) - 1
    // [NEG-TIME-FIX] make ceil division correct for negative leftFrame.
    auto floorDivInt64 = [] (juce::int64 a, juce::int64 b) -> juce::int64
    {
        // b must be > 0
        if (b <= 0) return 0;
        if (a >= 0) return a / b;
        return - ((-a + b - 1) / b); // floor division for negatives
    };

    auto ceilDivInt64 = [&] (juce::int64 a, juce::int64 b) -> juce::int64
    {
        // ceil(a/b) = -floor((-a)/b)
        return -floorDivInt64 (-a, b);
    };

    const juce::int64 span64 = (juce::int64) spanFrames;

    juce::int64 minAbsWanted = ceilDivInt64 (leftFrame, span64) - 1;

    // overscan one group to the left
    minAbsWanted -= 1;

    const juce::int64 minAbs = juce::jmax (earliestAbsGroup, minAbsWanted);
    const juce::int64 maxAbs = latestAbsGroup;

    const juce::int64 groupsToUse64 = maxAbs - minAbs + 1;
    if (groupsToUse64 < 2)
        return;

    const int groupsToUse = (int) juce::jmin<juce::int64> (groupsToUse64,
                                                          (juce::int64) std::numeric_limits<int>::max());

    //--------------------------------------------------------------------------
    // Choose a step for performance
    //--------------------------------------------------------------------------
    const int maxDrawablePoints = getMaxDrawablePoints (widthPixels);

    // [DECIMATOR-GRID-ANCHOR] zoom-based step (controls density by pixels)
    const double pxPerGroup = (double) spanFrames * zoomX;

    // Tuneable: larger => fewer points (more performance)
    const double desiredPxPerPoint = (widthPixels >= 2500 ? 1.6 : 1.25);

    const int stepByZoom = (pxPerGroup > 1.0e-12
                              ? (int) std::ceil (desiredPxPerPoint / pxPerGroup)
                              : maxDrawablePoints);

    // Safety cap so we never exceed maxDrawablePoints
    const int stepByCap = (groupsToUse > maxDrawablePoints
                             ? (int) std::ceil ((double) groupsToUse / (double) maxDrawablePoints)
                             : 1);

    const juce::int64 step = (juce::int64) juce::jmax (1, juce::jmax (stepByZoom, stepByCap));

    //--------------------------------------------------------------------------
    // [DECIMATOR-GRID-ANCHOR] Anchor chunk boundaries to a fixed global grid
    //--------------------------------------------------------------------------
    auto floorToGrid = [&] (juce::int64 v, juce::int64 grid) -> juce::int64
    {
        // [NEG-TIME-FIX] floor-to-grid for negative values too
        if (grid <= 0) return v;
        return floorDivInt64 (v, grid) * grid;
    };

    juce::int64 firstAbs = floorToGrid (minAbs, step);
    if (firstAbs > minAbs)
        firstAbs -= step;

    // overscan one chunk to the left (optional but helps continuity at edge)
    firstAbs = juce::jmax (earliestAbsGroup, firstAbs - step);

    // Reserve approx (not exact)
    const int approxOut = (int) std::ceil ((double) groupsToUse / (double) step) + 4;
    outGroups.reserve ((size_t) juce::jlimit (64, 8192, approxOut));
    outEndFrameIndex.reserve ((size_t) juce::jlimit (64, 8192, approxOut));

    // Emit chunks in chronological order: oldest -> newest
    for (juce::int64 absStart = firstAbs; absStart <= maxAbs; absStart += step)
    {
        const juce::int64 absEnd = juce::jmin (maxAbs, absStart + step - 1);

        // Skip chunks completely to the left of our needed range (except overscan)
        if (absEnd < minAbs - step)
            continue;

        FrameGroup agg;
        agg.momentaryMinDb =  std::numeric_limits<float>::infinity();
        agg.momentaryMaxDb = -std::numeric_limits<float>::infinity();
        agg.shortTermMinDb =  std::numeric_limits<float>::infinity();
        agg.shortTermMaxDb = -std::numeric_limits<float>::infinity();

        // [LRAG]
        agg.gateMinDb =  std::numeric_limits<float>::infinity();
        agg.gateMaxDb = -std::numeric_limits<float>::infinity();

        bool any = false;

        for (juce::int64 absIdx = absStart; absIdx <= absEnd; ++absIdx)
        {
            FrameGroup gg;
            if (! readGroupAbs (levelIndex, absIdx, gg))
                continue; // missing data at this timeline position (not written yet)

            any = true;

            agg.momentaryMinDb = std::min (agg.momentaryMinDb, gg.momentaryMinDb);
            agg.momentaryMaxDb = std::max (agg.momentaryMaxDb, gg.momentaryMaxDb);
            agg.shortTermMinDb = std::min (agg.shortTermMinDb, gg.shortTermMinDb);
            agg.shortTermMaxDb = std::max (agg.shortTermMaxDb, gg.shortTermMaxDb);
            agg.gateMinDb = std::min (agg.gateMinDb, gg.gateMinDb);
            agg.gateMaxDb = std::max (agg.gateMaxDb, gg.gateMaxDb);
        }

        if (! any)
            continue;

        // Timestamp: use CHUNK END (no center rounding drift)
        const juce::int64 endFrame = (absEnd + 1) * (juce::int64) spanFrames;

        outGroups.push_back (agg);
        outEndFrameIndex.push_back (endFrame);
    }

    if (outGroups.size() < 2)
    {
        outGroups.clear();
        outEndFrameIndex.clear();
    }
}
// [END VHC-HIST-VISIBLE-BUILDER]

//==============================================================================
// Bootstrap from processor (persistence)
//==============================================================================

// [BEGIN VHC-HIST-BOOTSTRAP]
void VolumeHistoryComponent::bootstrapHistoryFromProcessorIfNeeded()
{
    if (bootstrappedFromProcessor)
        return;

    const juce::int64 maxWritten = processor.getMaxWrittenFrameIndex();
    if (maxWritten == std::numeric_limits<juce::int64>::min())
        return;

    const juce::int64 startFrame = maxWritten - (juce::int64) (rawCapacityFrames - 1);

    // [ROLLING-LRA] reset so the first valid second will be pushed
    lastSecondPushed = std::numeric_limits<juce::int64>::min();

    for (juce::int64 fi = startFrame; fi <= maxWritten; ++fi)
    {
        float mVal = 0.0f, sVal = 0.0f;
        if (! processor.getDerivedRmsAtFrameIndex (fi, mVal, sVal))
            continue;

        float gate = -200.0f;
        processor.getLraGateLufsAtFrameIndex (fi, gate); // ok if missing

        // Phase 2: processor accessor returns LUFS (dB) already.
        const float lufsM = mVal;
        const float lufsS = sVal;

        FrameGroup fg;
        fg.momentaryMinDb = lufsM;
        fg.momentaryMaxDb = lufsM;
        fg.shortTermMinDb = lufsS;
        fg.shortTermMaxDb = lufsS;
        fg.gateMinDb = gate;
        fg.gateMaxDb = gate;

        writeGroupAbs (0, fi, fg);

        // [ROLLING-LRA] sample once per second from the 60 Hz timeline (fi is a 60 Hz frame index)
        const juce::int64 absSecondIndex = floorDivInt64 (fi, 60);

        if (absSecondIndex != lastSecondPushed)
        {
            lastSecondPushed = absSecondIndex;
            pushSecondSample (absSecondIndex, lufsS, gate);
        }
    }

    // Rebuild higher LOD levels over the same time range
    for (int level = 1; level < maxLevels; ++level)
    {
        const int span = levels[(size_t) level].spanFrames;
        const juce::int64 minAbs = floorDivInt64 (startFrame, (juce::int64) span);
        const juce::int64 maxAbs = floorDivInt64 (maxWritten, (juce::int64) span);

        for (juce::int64 absIdx = minAbs; absIdx <= maxAbs; ++absIdx)
            recomputeGroupAbsFromChildren (level, absIdx);
    }

    haveNowFrameIndex = true;
    nowFrameIndex = maxWritten;

    havePlayheadFrameIndex = true;
    playheadFrameIndex = maxWritten;

    bootstrappedFromProcessor = true;

    // [ROLLING-LRA] after bootstrap, rebuild rolling curve over all stored seconds (chunked in timerCallback)
    rollingRebuildInProgress = true;
    rollingRebuildMinSecond = std::numeric_limits<juce::int64>::max();
    rollingRebuildMaxSecond = std::numeric_limits<juce::int64>::min();

    for (int i = 0; i < secondsCapacity; ++i)
    {
        const auto t = secAbsIndexTag[(size_t) i];
        if (t != (juce::int64) -1)
        {
            rollingRebuildMinSecond = juce::jmin (rollingRebuildMinSecond, t);
            rollingRebuildMaxSecond = juce::jmax (rollingRebuildMaxSecond, t);
        }
    }

    if (rollingRebuildMinSecond <= rollingRebuildMaxSecond)
        rollingRebuildNextSecond = rollingRebuildMinSecond;
    else
        rollingRebuildInProgress = false;
}
// [END VHC-HIST-BOOTSTRAP]

//==============================================================================
// History update
//==============================================================================

// [BEGIN VHC-HIST-HELPERS]
juce::int64 VolumeHistoryComponent::floorDivInt64 (juce::int64 a, juce::int64 b) noexcept
{
    if (b <= 0) return 0;
    if (a >= 0) return a / b;
    return - ((-a + b - 1) / b);
}
// [END VHC-HIST-HELPERS]

//==============================================================================
// [ROLLING-LRA] Helpers
//==============================================================================

int VolumeHistoryComponent::wrapSecondSlot (juce::int64 absSecondIndex) const noexcept
{
    if (secondsCapacity <= 0)
        return 0;

    juce::int64 m = absSecondIndex % (juce::int64) secondsCapacity;
    if (m < 0) m += (juce::int64) secondsCapacity;
    return (int) m;
}

void VolumeHistoryComponent::pushSecondSample (juce::int64 absSecondIndex,
                                              float shortTermLufs,
                                              float gateLufs)
{
    if (secondsCapacity <= 0)
        return;

    const int slot = wrapSecondSlot (absSecondIndex);

    secShortTermLufs[(size_t) slot] = shortTermLufs;
    secGateLufs[(size_t) slot]      = gateLufs;

    secAbsIndexTag[(size_t) slot] = absSecondIndex;

    // Compute rolling LRA for this second immediately (cheap)
    recomputeRollingLraForSecond (absSecondIndex);
}

void VolumeHistoryComponent::recomputeRollingLraForSecond (juce::int64 absSecondIndex)
{
    if (secondsCapacity <= 0)
        return;

    const int slotNow = wrapSecondSlot (absSecondIndex);
    if (secAbsIndexTag[(size_t) slotNow] != absSecondIndex)
        return;

    const int windowSeconds = rollingWindowSecondsCached;
    const int N = juce::jlimit (1, 120, windowSeconds);

    const float gateHere = secGateLufs[(size_t) slotNow];
    const float gate = juce::jmax (-70.0f, gateHere);

    // Collect gated samples from last N seconds (including current)
    std::array<float, 120> tmp {};
    int count = 0;

    for (int k = 0; k < N; ++k)
    {
        const juce::int64 sIdx = absSecondIndex - (juce::int64) k;
        const int slot = wrapSecondSlot (sIdx);

        if (secAbsIndexTag[(size_t) slot] != sIdx)
            continue;

        const float v = secShortTermLufs[(size_t) slot];
        if (v >= gate)
            tmp[(size_t) count++] = v;
    }

    float out = 0.0f;

    if (count >= 4)
    {
        std::sort (tmp.begin(), tmp.begin() + count);

        const int rank10 = (int) std::ceil (0.10 * (double) count);
        const int rank95 = (int) std::ceil (0.95 * (double) count);

        const int i10 = juce::jlimit (0, count - 1, rank10 - 1);
        const int i95 = juce::jlimit (0, count - 1, rank95 - 1);

        out = juce::jmax (0.0f, tmp[(size_t) i95] - tmp[(size_t) i10]);
    }

    secRollingLraLu[(size_t) slotNow] = out;
}

void VolumeHistoryComponent::startRollingRebuildIfWindowChanged()
{
    const int current = processor.getRollingLraWindowSeconds();
    if (current == rollingWindowSecondsCached)
        return;

    rollingWindowSecondsCached = current;

    // Start a rebuild over all known seconds (chunked)
    rollingRebuildInProgress = true;
    rollingRebuildMinSecond = std::numeric_limits<juce::int64>::max();
    rollingRebuildMaxSecond = std::numeric_limits<juce::int64>::min();

    for (int i = 0; i < secondsCapacity; ++i)
    {
        const auto t = secAbsIndexTag[(size_t) i];
        if (t != (juce::int64) -1)
        {
            rollingRebuildMinSecond = juce::jmin (rollingRebuildMinSecond, t);
            rollingRebuildMaxSecond = juce::jmax (rollingRebuildMaxSecond, t);
        }
    }

    if (rollingRebuildMinSecond <= rollingRebuildMaxSecond)
        rollingRebuildNextSecond = rollingRebuildMinSecond;
    else
        rollingRebuildInProgress = false;
}

void VolumeHistoryComponent::rebuildRollingLraStep (int maxSecondsToProcess)
{
    if (! rollingRebuildInProgress)
        return;

    const int maxCount = juce::jmax (1, maxSecondsToProcess);
    int done = 0;

    while (done < maxCount && rollingRebuildNextSecond <= rollingRebuildMaxSecond)
    {
        // Only recompute if that second exists
        const int slot = wrapSecondSlot (rollingRebuildNextSecond);
        if (secAbsIndexTag[(size_t) slot] == rollingRebuildNextSecond)
            recomputeRollingLraForSecond (rollingRebuildNextSecond);

        ++rollingRebuildNextSecond;
        ++done;
    }

    if (rollingRebuildNextSecond > rollingRebuildMaxSecond)
        rollingRebuildInProgress = false;
}
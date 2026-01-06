#include "VolumeHistoryComponent.h"
#include "PluginProcessor.h"

#include <cmath>
#include <algorithm>
#include <limits>

// Helper: format seconds as HH:MM:SS
static juce::String formatTimeHMS (double seconds)
{
    const bool neg = (seconds < 0.0);
    seconds = std::abs (seconds);

    const int total = (int) std::floor (seconds + 0.5);
    const int h = total / 3600;
    const int m = (total % 3600) / 60;
    const int s = total % 60;

    const juce::String core = juce::String::formatted ("%02d:%02d:%02d", h, m, s);
    return neg ? "-" + core : core;
}

// [TIMECODE-USER] parse either seconds (e.g. -2.0) or HH:MM:SS (e.g. -00:00:02)
static bool parseUserOffsetSeconds (juce::String text, double& outSeconds)
{
    text = text.trim();

    if (text.isEmpty())
        return false;

    // If it contains ':' parse as HH:MM:SS
    if (text.containsChar (':'))
    {
        bool neg = false;
        if (text.startsWithChar ('-'))
        {
            neg = true;
            text = text.substring (1).trim();
        }

        auto parts = juce::StringArray::fromTokens (text, ":", "");
        parts.removeEmptyStrings();

        if (parts.size() != 3)
            return false;

        const int hh = parts[0].getIntValue();
        const int mm = parts[1].getIntValue();
        const int ss = parts[2].getIntValue();

        double sec = (double) hh * 3600.0 + (double) mm * 60.0 + (double) ss;
        if (neg) sec = -sec;

        outSeconds = sec;
        return true;
    }

    // Otherwise treat as seconds
    outSeconds = text.getDoubleValue();
    return true;
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

    // [VIEW-NAV]
    viewTopDb = (double) maxDb;     // top starts at 0 dBFS
    viewRightFrame = 0.0;           // will follow right edge once we have data
    followRightEdge = true;

    bootstrapHistoryFromProcessorIfNeeded(); // [STATE-PERSIST]

    setOpaque (true);

    // UI update rate
    startTimerHz (30);

    // [FOLLOW-BUTTON]
    followButton.setButtonText ("Follow");
    followButton.setToggleState (followRightEdge, juce::dontSendNotification);
    followButton.onClick = [this]
    {
        followRightEdge = followButton.getToggleState();

        if (followRightEdge)
        {
            const int w = getWidth();

            if (havePlayheadFrameIndex && zoomX > 1.0e-12 && w > 1)
            {
                const double visibleFrames = (double) w / zoomX;
                constexpr double playheadXFrac = 0.5; // center
                viewRightFrame = (double) playheadFrameIndex + visibleFrames * (1.0 - playheadXFrac);
            }
            else if (haveNowFrameIndex)
            {
                viewRightFrame = (double) nowFrameIndex;
            }

            clampViewRightFrame (w);
            repaint();
        }
    };

    addAndMakeVisible (followButton);

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

juce::int64 VolumeHistoryComponent::floorDivInt64 (juce::int64 a, juce::int64 b) noexcept
{
    if (b <= 0) return 0;
    if (a >= 0) return a / b;
    return - ((-a + b - 1) / b);
}

void VolumeHistoryComponent::bootstrapHistoryFromProcessorIfNeeded()
{
    if (bootstrappedFromProcessor)
        return;

    const juce::int64 maxWritten = processor.getMaxWrittenFrameIndex();
    if (maxWritten == std::numeric_limits<juce::int64>::min())
        return;

    // Fill L0 for the last rawCapacityFrames (3 hours @ 60 Hz)
    const juce::int64 startFrame = maxWritten - (juce::int64) (rawCapacityFrames - 1);

    for (juce::int64 fi = startFrame; fi <= maxWritten; ++fi)
    {
        float mRms = 0.0f, sRms = 0.0f;
        if (! processor.getDerivedRmsAtFrameIndex (fi, mRms, sRms))
            continue;

        const float dbM = juce::Decibels::gainToDecibels (mRms, minDb);
        const float dbS = juce::Decibels::gainToDecibels (sRms, minDb);

        FrameGroup fg;
        fg.momentaryMinDb = dbM;
        fg.momentaryMaxDb = dbM;
        fg.shortTermMinDb = dbS;
        fg.shortTermMaxDb = dbS;

        writeGroupAbs (0, fi, fg);
    }

    // Rebuild higher levels over the same time range
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
}

bool VolumeHistoryComponent::drainProcessorFifo()
{
    // [TIMEBASE-PLAYHEAD] cache samples-per-visual-frame (60 Hz) from processor
    if (frameSamples <= 0)
    frameSamples = processor.getFrameSamples();
    constexpr int chunkSize = 512;
    float momentaryValues [chunkSize];
    float shortTermValues [chunkSize];
    juce::int64 frameIndex60Hz [chunkSize];
    int playingFlags [chunkSize];

    bool readAny = false;

    for (;;)
    {
        const int numRead = processor.readLoudnessFromFifo (momentaryValues,
                                                           shortTermValues,
                                                           frameIndex60Hz,
                                                           playingFlags,
                                                           chunkSize);
        if (numRead <= 0)
            break;

        readAny = true;

        for (int i = 0; i < numRead; ++i)
            pushFrameToHistory (momentaryValues[i], shortTermValues[i],
                                frameIndex60Hz[i], playingFlags[i]);
    }

    return readAny;
}

void VolumeHistoryComponent::pushFrameToHistory (float momentaryRms,
                                                 float shortTermRms,
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

    const float dbM = juce::Decibels::gainToDecibels (momentaryRms, minDb);
    const float dbS = juce::Decibels::gainToDecibels (shortTermRms, minDb);

    FrameGroup fg;
    fg.momentaryMinDb = dbM;
    fg.momentaryMaxDb = dbM;
    fg.shortTermMinDb = dbS;
    fg.shortTermMaxDb = dbS;

    // L0 overwrite at absolute frame index
    writeGroupAbs (0, frameIndex, fg);

    // Recompute parent groups up the pyramid (overwrite-safe)
    for (int level = 1; level < maxLevels; ++level)
    {
        const auto& L = levels[(size_t) level];
        const juce::int64 absGroupIndex = frameIndex / (juce::int64) L.spanFrames;
        recomputeGroupAbsFromChildren (level, absGroupIndex);
    }
}

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
    }

    if (! any)
        return; // keep old value/tag (or empty)

    writeGroupAbs (levelIndex, absGroupIndex, agg);
}

//==============================================================================
// History access
//==============================================================================

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

//==============================================================================
// LOD selection  [STEP2-LOD-CAP]
//==============================================================================

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

//==============================================================================
// Visible groups builder  [TIMEBASE-FIX]  [DECIMATOR-GRID-ANCHOR]
//
// Key idea:
// - We still cap point count with a step.
// - But chunk boundaries are ANCHORED to an absolute grid (multiples of step),
//   not to the moving left edge. This prevents "accordion" for BOTH renderers.
// - Timestamp uses chunk END (no center/rounding drift).
//==============================================================================
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
    const float topDb          = (float) viewTopDb; // [VIEW-NAV]
    const float bottomDb       = topDb - effectiveRange;

    const float clamped = juce::jlimit (bottomDb, topDb, db);
    const float norm = (clamped - bottomDb) / effectiveRange;

    return height * (1.0f - norm);
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

    const bool topDbChanged = (std::abs (viewTopDb - cachedBgTopDb) > 1.0e-12); // [VIEW-NAV]

    if (! staticBackgroundDirty && ! sizeChanged && ! zoomYChanged && ! topDbChanged
        && cachedStaticBackground.isValid())
        return;

    cachedBgW = w;
    cachedBgH = h;
    cachedBgZoomY = zoomY;
    cachedBgTopDb = viewTopDb; // [VIEW-NAV]
    staticBackgroundDirty = false;

    cachedStaticBackground = juce::Image (juce::Image::RGB, w, h, true);
    juce::Graphics gg (cachedStaticBackground);

    auto backgroundColour = juce::Colour::fromRGB (16, 30, 50);
    gg.fillAll (backgroundColour);

    gg.setColour (juce::Colours::darkgrey.withMultipliedAlpha (0.4f));
    const int numLines = 4;
    const float effectiveRange = (float) (baseDbRange / zoomY);
    for (int i = 0; i <= numLines; ++i)
    {
        const float db = (float) viewTopDb - (effectiveRange / (float) numLines) * (float) i;
        const float y  = dbToY (db, (float) h);
        gg.drawHorizontalLine ((int) std::round (y), 0.0f, (float) w);
    }

    const float rulerBaseY = (float) h - 2.0f;
    gg.setColour (juce::Colours::darkgrey.withMultipliedAlpha (0.7f));
    gg.drawHorizontalLine ((int) std::round (rulerBaseY), 0.0f, (float) w);
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
        scratchVisibleEndFrameIndex.reserve (reserveCount); // [TIMEBASE-FIX]
        scratchRepMomentaryDb.reserve (reserveCount);
        scratchRepShortTermDb.reserve (reserveCount);
    }

    // Do NOT reset tickStepIndex here; hysteresis should smoothly adapt.
    markStaticBackgroundDirty();

    // [FOLLOW-BUTTON] small toggle in the top-right
    followButton.setBounds (getWidth() - 88, 6, 80, 22);
}

//==============================================================================
// DRAW
//==============================================================================

void VolumeHistoryComponent::paint (juce::Graphics& g)
{
    auto bounds = getLocalBounds().toFloat();
    const int width  = (int) bounds.getWidth();
    const int height = (int) bounds.getHeight();

    // [VIEW-NAV] Follow mode: follow playhead (so overwrite outside view jumps to that section)
    // Keep playhead at center while following.
    if (followRightEdge)
    {
        if (havePlayheadFrameIndex && zoomX > 1.0e-12)
        {
            const double visibleFrames = (double) width / zoomX;
            constexpr double playheadXFrac = 0.5; // 0.5 = center
            viewRightFrame = (double) playheadFrameIndex + visibleFrames * (1.0 - playheadXFrac);
        }
        else if (haveNowFrameIndex)
        {
            // fallback if playhead is unknown
            viewRightFrame = (double) nowFrameIndex;
        }
    }

followButton.setToggleState (followRightEdge, juce::dontSendNotification);

clampViewRightFrame (width);

    if (width <= 1 || height <= 0)
        return;

    rebuildStaticBackgroundIfNeeded();
    if (cachedStaticBackground.isValid())
        g.drawImageAt (cachedStaticBackground, 0, 0);

    const int selectedLevel = selectBestLevelForCurrentZoom (width);

    buildVisibleGroupsForLevel (selectedLevel, width, scratchVisibleGroups, scratchVisibleEndFrameIndex); // [TIMEBASE-FIX]

    const size_t n = scratchVisibleGroups.size();
    const float w = bounds.getWidth();
    const float h = bounds.getHeight();

    if (n >= 2)
    {
        computeRepresentativeCurves (scratchVisibleGroups, scratchRepMomentaryDb, scratchRepShortTermDb);

        scratchPathBandM.clear();
        scratchPathBandS.clear();

        if (showLines)
        {
            scratchPathRepM.clear();
            scratchPathRepS.clear();
        }

        const float bandRangeThresholdDb = 3.0f;

        bool startedRepM = false, startedRepS = false;

        for (size_t i = 0; i < n; ++i)
        {
            const double xD = (double) w - (viewRightFrame - (double) scratchVisibleEndFrameIndex[i]) * zoomX;
            float x = (float) xD;
            if (x < -10.0f)
                continue;

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

            if (showLines)
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
       // Draw momentary first, then short-term on top.
       if (showBands && selectedLevel > 0)
       {
           g.setColour (juce::Colour::fromRGB (95, 117, 140).withMultipliedAlpha (0.7f)); // momentary
           g.strokePath (scratchPathBandM, juce::PathStrokeType (1.2f));

           g.setColour (juce::Colour::fromRGB (0, 80, 180).withMultipliedAlpha (0.6f));   // short-term (on top)
           g.strokePath (scratchPathBandS, juce::PathStrokeType (1.0f));
       } 

        // Lines (stroked paths only)
        // Draw momentary first, then short-term on top.
        if (showLines)
        {
            g.setColour (juce::Colour::fromRGB (95, 117, 140)); // momentary
            g.strokePath (scratchPathRepM, juce::PathStrokeType (2.0f));

            g.setColour (juce::Colour::fromRGB (0, 80, 180).withMultipliedAlpha (0.95f)); // short-term (on top)
            g.strokePath (scratchPathRepS, juce::PathStrokeType (1.5f));
        }
    }

    //==============================================================================
    // [PLAYHEAD-LINE] draw DAW playhead position on top of the graph
    // Right edge is "furthest written" (nowFrameIndex). Playhead can be left of it
    // during loop/rewind, which is intended.
    //==============================================================================
    if (haveNowFrameIndex && havePlayheadFrameIndex && zoomX > 1.0e-12)
    {
        const double framesFromRight = viewRightFrame - (double) playheadFrameIndex;
        const float x = (float) ((double) width - framesFromRight * zoomX);

        if (x >= -2.0f && x <= (float) width + 2.0f)
        {
            // Thin, slightly transparent so it doesn't dominate
            g.setColour (juce::Colours::white.withMultipliedAlpha (0.55f));
            g.drawLine (x, 0.0f, x, (float) height, 1.0f);
        }
    }

    //==============================================================================
    // [DBFS-SCALE] Right-side dBFS scale (overlay)
    // Shows the currently visible dB range (affected by zoomY).
    //==============================================================================
    {
        // Keep the scale above the time ruler area
        const float rulerHeightPx = 16.0f;
        const auto scaleArea = bounds.withTrimmedBottom (rulerHeightPx);

        const float scaleH = scaleArea.getHeight();
        if (scaleH > 20.0f)
        {
            const float topDb = (float) viewTopDb;
            const float effectiveRange = (float) (baseDbRange / zoomY);
            const float bottomDb = topDb - effectiveRange;

            // Choose a tick step based on available pixel height
            static constexpr float candidates[] = { 1.0f, 2.0f, 3.0f, 6.0f, 10.0f, 12.0f, 20.0f, 30.0f };
            const float desiredTickSpacingPx = 32.0f;
            const float approxTicks = juce::jmax (1.0f, scaleH / desiredTickSpacingPx);
            const float approxStepDb = effectiveRange / approxTicks;

            float stepDb = candidates[(int) (sizeof (candidates) / sizeof (candidates[0])) - 1];
            for (float s : candidates)
            {
            if (s >= approxStepDb) { stepDb = s; break; }
            }

            const float rightX = scaleArea.getRight() - 1.0f;
            const float tickLen = 6.0f;
            const float labelWidth = 56.0f;
            const float labelHeight = 16.0f;

            g.setFont (14.0f); // [UI-FONTS]

            // Optional header
            g.setColour (juce::Colours::white.withMultipliedAlpha (0.6f));
            g.drawText ("dBFS",
                        (int) (rightX - labelWidth - 2.0f),
                        (int) (scaleArea.getY() + 2.0f),
                        (int) labelWidth,
                        12,
                        juce::Justification::right);

            // Major ticks: from first tick >= bottomDb up to topDb
            const float firstTick = std::ceil (bottomDb / stepDb) * stepDb;

            g.setColour (juce::Colours::white.withMultipliedAlpha (0.55f));

            for (float db = firstTick; db <= topDb + 0.001f; db += stepDb)
            {
                const float yLocal = dbToY (db, scaleH);
                const float y = scaleArea.getY() + yLocal;

                // Tick
                g.drawLine (rightX - tickLen, y, rightX, y, 1.0f);

                // Label (integer dB)
                const int dbInt = (int) std::round (db);
                g.drawText (juce::String (dbInt),
                            (int) (rightX - tickLen - labelWidth - 2.0f),
                            (int) (y - labelHeight * 0.5f),
                            (int) labelWidth,
                            (int) labelHeight,
                            juce::Justification::right);
            }

            // A faint vertical line to separate the scale
            g.setColour (juce::Colours::white.withMultipliedAlpha (0.18f));
            g.drawLine (rightX, scaleArea.getY(), rightX, scaleArea.getBottom(), 1.0f);
        }
    }

    //==========================================================================
    // [RULER-FRAMES]  [VIEW-NAV]
    // Ruler must use the SAME x-mapping as the curves:
    // x = width - (viewRightFrame - tickFrame) * zoomX
    //==========================================================================
    const float rulerHeight   = 22.0f;   // [UI-FONTS] bigger
    const float rulerBaseY    = (float) height - 2.0f;
    const float tickTopY      = rulerBaseY - 6.0f;
    const float textTopY      = rulerBaseY - 14.0f;

    if (zoomX > 1.0e-12 && visualFrameRate > 0.0 && haveNowFrameIndex)
    {
        const juce::int64 rightFrameI = (juce::int64) std::floor (viewRightFrame);
        const double safeZoomX = zoomX;

        const double framesByWidth = (double) width / safeZoomX;
        const juce::int64 visibleFrames = (juce::int64) std::ceil (juce::jmax (1.0, framesByWidth));

        // Earliest frame we can possibly still have in our history window
        const juce::int64 earliestStoredFrame = nowFrameIndex - (juce::int64) rawCapacityFrames + 1;

        const juce::int64 leftFrame = juce::jmax (earliestStoredFrame, rightFrameI - visibleFrames);

        const double tickStepSec = getTickStepSecondsWithHysteresis (width);
        const juce::int64 tickStepFrames = (juce::int64) juce::jmax (1.0, std::round (tickStepSec * visualFrameRate));

        // Negative-safe floor-to-multiple
        auto floorDivInt64 = [] (juce::int64 a, juce::int64 b) -> juce::int64
        {
            if (b <= 0) return 0;
            if (a >= 0) return a / b;
            return - ((-a + b - 1) / b);
        };

        const juce::int64 lastTickFrame = floorDivInt64 (rightFrameI, tickStepFrames) * tickStepFrames;

        g.setColour (juce::Colours::white);
        g.setFont (14.0f);                   // [UI-FONTS] bigger

        for (juce::int64 tickFrame = lastTickFrame; tickFrame >= leftFrame; tickFrame -= tickStepFrames)
        {
            const float x = (float) ((double) width - (viewRightFrame - (double) tickFrame) * zoomX);

            if (x < -2.0f)
                break;

            if (x > (float) width + 2.0f)
                continue;

            g.drawLine (x, tickTopY, x, rulerBaseY, 1.0f);

            const double tSec = (double) tickFrame / visualFrameRate
                              + processor.getTimecodeOffsetSeconds()
                              + processor.getUserTimecodeOffsetSeconds();
            const float textWidth = 84.0f;

            g.drawText (formatTimeHMS (tSec),
                        x - textWidth * 0.5f,
                        textTopY,
                        textWidth,
                        rulerHeight,
                        juce::Justification::centred);
        }
    }

    // Overlay info
    g.setColour (juce::Colours::white);
    g.setFont (14.0f);

    const auto& Lsel         = levels[(size_t) selectedLevel];
    const double spanSeconds = (double) Lsel.spanFrames / visualFrameRate;

    const double tickStepNow = getTickStepSecondsWithHysteresis (width);

    juce::String info = "Level: " + juce::String (selectedLevel) +
                        " (span " + juce::String (Lsel.spanFrames) + " frames, " +
                        juce::String (spanSeconds, 3) + " s)" +
                        " | ZoomX: " + juce::String (zoomX, 4) +
                        " | ZoomY: " + juce::String (zoomY, 2) +
                        " | Tick: " + juce::String (tickStepNow, 3) + "s" +
                        " | Bands: " + juce::String (showBands ? "ON" : "OFF") +
                        " | Lines: " + juce::String (showLines ? "ON" : "OFF");
                        // [DEBUG-OVERLAY] Put debug fields on a second line so they don't get clipped
                        juce::String info2;

                        info2 += "Follow:" + juce::String (followRightEdge ? "ON" : "OFF")
                              +  " viewRight:" + juce::String (viewRightFrame, 1);

                        if (haveNowFrameIndex)
                            info2 += " now:" + juce::String (nowFrameIndex);
                        if (havePlayheadFrameIndex)
                            info2 += " play:" + juce::String (playheadFrameIndex);

                        info2 += " TCoff:" + juce::String (processor.getTimecodeOffsetSeconds(), 3) + "s"
                              +  " TCu:"   + juce::String (processor.getUserTimecodeOffsetSeconds(), 3) + "s"
                              +  " hostSamp:" + juce::String (processor.hostHasTimeInSamples() ? "Y" : "N")
                              +  " hostSec:"  + juce::String (processor.hostHasTimeInSeconds() ? "Y" : "N");

                        if (processor.hostHasTimeInSamples())
                            info2 += " S:" + juce::String (processor.getLastHostTimeInSamples());
                        if (processor.hostHasTimeInSeconds())
                            info2 += " s:" + juce::String (processor.getLastHostTimeInSeconds(), 3);

                        // Draw it
                        const int overlayX = 8;
                        const int overlayW = (int) std::min (bounds.getWidth() - 16.0f, 1600.0f);

                        // Main line
                        g.setColour (juce::Colours::white);
                        g.setFont (14.0f);
                        g.drawFittedText (info, overlayX, 6, overlayW, 18, juce::Justification::topLeft, 1);

                        // Debug line
                        g.setColour (juce::Colours::white.withMultipliedAlpha (0.85f));
                        g.setFont (12.0f);
                        g.drawFittedText (info2, overlayX, 24, overlayW, 16, juce::Justification::topLeft, 1);
}

//==============================================================================
// [VIEW-NAV] clamp + pan/zoom
//==============================================================================

void VolumeHistoryComponent::clampViewRightFrame (int widthPixels) noexcept
{
    if (! haveNowFrameIndex || zoomX <= 1.0e-12 || widthPixels <= 1)
        return;

    const double visibleFrames = (double) widthPixels / zoomX;

    const double maxRight = (double) nowFrameIndex;
    const double earliest = (double) nowFrameIndex - (double) rawCapacityFrames + 1.0;

    double minRight = earliest + visibleFrames;
    if (minRight > maxRight)
        minRight = maxRight;

    viewRightFrame = juce::jlimit (minRight, maxRight, viewRightFrame);
}

void VolumeHistoryComponent::panTime (float wheelDelta)
{
    const int w = getWidth();
    if (w <= 1 || wheelDelta == 0.0f || zoomX <= 1.0e-12)
        return;

    if (followRightEdge)
        viewRightFrame = (double) nowFrameIndex;

    const double visibleFrames = (double) w / zoomX;

    // Pan by ~15% of the visible range per wheel "unit"
    const double panFrames = (double) wheelDelta * visibleFrames * 0.15;

    // Wheel up (positive) moves to earlier time (left)
    viewRightFrame -= panFrames;

    followRightEdge = false;
    clampViewRightFrame (w);
}

void VolumeHistoryComponent::panDb (float wheelDelta)
{
    const int h = getHeight();
    if (h <= 1 || wheelDelta == 0.0f)
        return;

    const double effectiveRange = (double) baseDbRange / zoomY;

    // Pan by ~10% of visible range per wheel "unit"
    const double panDbAmount = (double) wheelDelta * effectiveRange * 0.10;

    viewTopDb += panDbAmount;

    const double topMin = viewMinDbLimit + effectiveRange;
    const double topMax = viewMaxDbLimit;
    viewTopDb = juce::jlimit (topMin, topMax, viewTopDb);

    markStaticBackgroundDirty();
}

void VolumeHistoryComponent::applyHorizontalZoom (float wheelDelta, float anchorX)
{
    const int w = getWidth();
    if (w <= 1 || wheelDelta == 0.0f)
        return;

    if (! haveNowFrameIndex)
        return;

    if (followRightEdge)
        viewRightFrame = (double) nowFrameIndex;

    const double oldZoom = zoomX;

    const double zoomBase   = 1.1;
    const double zoomFactor = std::pow (zoomBase, (double) wheelDelta);

    zoomX *= zoomFactor;
    zoomX  = juce::jlimit (minZoomX, maxZoomX, zoomX);

    // Keep the timeline frame under the mouse fixed
    const double ax = juce::jlimit (0.0, (double) w, (double) anchorX);
    const double frameAtCursor = viewRightFrame - ((double) w - ax) / (oldZoom > 1.0e-12 ? oldZoom : 1.0e-12);

    viewRightFrame = frameAtCursor + ((double) w - ax) / (zoomX > 1.0e-12 ? zoomX : 1.0e-12);

    followRightEdge = false;
    hasCustomZoomX = true;

    clampViewRightFrame (w);
}

void VolumeHistoryComponent::applyVerticalZoom (float wheelDelta, float anchorY)
{
    const int h = getHeight();
    if (h <= 1 || wheelDelta == 0.0f)
        return;

    const double oldZoomY = zoomY;

    const double zoomBase   = 1.1;
    const double zoomFactor = std::pow (zoomBase, (double) wheelDelta);

    zoomY *= zoomFactor;
    zoomY  = juce::jlimit (minZoomY, maxZoomY, zoomY);

    // [VIEW-NAV-Y-LIMITS]
    // Prevent impossible states where the visible range exceeds the allowed view limits.
    // effectiveRange = baseDbRange / zoomY must be <= (viewMaxDbLimit - viewMinDbLimit)
    const double maxVisibleRange = viewMaxDbLimit - viewMinDbLimit;
    const double minZoomYByLimits = (double) baseDbRange / maxVisibleRange;
    zoomY = juce::jmax (zoomY, minZoomYByLimits);

    // Keep the dB under the mouse fixed while zooming
    const double ay = juce::jlimit (0.0, (double) h, (double) anchorY);
    const double norm = 1.0 - (ay / (double) h);

    const double effectiveOld = (double) baseDbRange / oldZoomY;
    const double bottomOld = viewTopDb - effectiveOld;
    const double dbAtCursor = bottomOld + norm * effectiveOld;

    const double effectiveNew = (double) baseDbRange / zoomY;
    const double bottomNew = dbAtCursor - norm * effectiveNew;
    double topNew = bottomNew + effectiveNew;

    const double topMin = viewMinDbLimit + effectiveNew; // ensures bottom >= viewMinDbLimit
    const double topMax = viewMaxDbLimit;                // allows positive dBFS
    topNew = juce::jlimit (topMin, topMax, topNew);
    
    viewTopDb = topNew;

    markStaticBackgroundDirty();
}

//==============================================================================
// [UI-RULERS] Areas + resets
//==============================================================================

juce::Rectangle<int> VolumeHistoryComponent::getTimeRulerArea() const
{
    // Keep consistent with paint() ruler height
    const int rulerH = 22; // [UI-FONTS] 150% bigger than before
    return { 0, getHeight() - rulerH, getWidth(), rulerH };
}

juce::Rectangle<int> VolumeHistoryComponent::getDbRulerArea() const
{
    // Right-side strip used for dB scale interaction
    const int rulerW = 72; // room for bigger labels
    const int timeRulerH = getTimeRulerArea().getHeight();

    return { getWidth() - rulerW, 0, rulerW, getHeight() - timeRulerH };
}

void VolumeHistoryComponent::resetXViewDefault()
{
    if (getWidth() <= 1 || visualFrameRate <= 0.0)
        return;

    const double desiredVisibleSeconds = 10.0;
    zoomX = (double) getWidth() / (desiredVisibleSeconds * visualFrameRate);
    zoomX = juce::jlimit (minZoomX, maxZoomX, zoomX);

    followRightEdge = true;
    if (haveNowFrameIndex)
        viewRightFrame = (double) nowFrameIndex;

    clampViewRightFrame (getWidth());
}

void VolumeHistoryComponent::fitXViewMaxZoomOut()
{
    if (getWidth() <= 1 || ! haveNowFrameIndex)
        return;

    // [FIT-FIRST-WRITTEN]
    // Find earliest written L0 frame index currently present in the GUI ring.
    // This avoids fitting the full 3h window when only a few minutes exist.
    const auto& tags = levels[0].absGroupIndexTag;

    juce::int64 minWritten = std::numeric_limits<juce::int64>::max();
    bool haveMin = false;

    for (const auto& t : tags)
    {
        if (t != (juce::int64) -1)
        {
            minWritten = juce::jmin (minWritten, t);
            haveMin = true;
        }
    }

    const juce::int64 right = nowFrameIndex;
    const juce::int64 left  = (haveMin ? minWritten : (right - (juce::int64) rawCapacityFrames + 1));

    const double visibleFrames = (double) (right - left + 1);
    zoomX = (double) getWidth() / juce::jmax (1.0, visibleFrames);
    zoomX = juce::jlimit (minZoomX, maxZoomX, zoomX);

    // Fit view to [left..right] and keep it there (do not force follow)
    followRightEdge = false;
    viewRightFrame = (double) right;
    clampViewRightFrame (getWidth());
}

void VolumeHistoryComponent::resetYViewDefault()
{
    zoomY = 1.0;
    viewTopDb = (double) maxDb; // default top at 0 dBFS
    markStaticBackgroundDirty();
}

//==============================================================================
// Mouse
//==============================================================================

void VolumeHistoryComponent::mouseWheelMove (const juce::MouseEvent& event,
                                             const juce::MouseWheelDetails& wheel)
{
    // Use deltaY when available, otherwise fall back to deltaX
    const float wheelDelta = (wheel.deltaY != 0.0f ? wheel.deltaY : wheel.deltaX);
    if (wheelDelta == 0.0f)
        return;

    // Gesture map:
    // - Shift + wheel: pan time
    // - Alt   + wheel: zoom Y at mouse
    // - Ctrl  + wheel: pan Y
    // - wheel: zoom X at mouse
    if (event.mods.isShiftDown())
        panTime (wheelDelta);
    else if (event.mods.isAltDown())
        applyVerticalZoom (wheelDelta, event.position.y);
    else if (
       #if JUCE_MAC
            event.mods.isCommandDown()
   #else
        event.mods.isCtrlDown()
   #endif
    )
        panDb (wheelDelta);
    else
        applyHorizontalZoom (wheelDelta, event.position.x);

    repaint();
}

void VolumeHistoryComponent::mouseDoubleClick (const juce::MouseEvent& event)
{
    const auto p = event.getPosition();

    // Double-click dB ruler resets Y view
    if (getDbRulerArea().contains (p))
    {
        resetYViewDefault();
        repaint();
        return;
    }

    // Double-click time ruler resets X view
    // Shift + double-click: fit/max zoom out
    if (getTimeRulerArea().contains (p))
    {
        if (event.mods.isShiftDown())
            fitXViewMaxZoomOut();
        else
            resetXViewDefault();

        followButton.setToggleState (followRightEdge, juce::dontSendNotification);
        repaint();
        return;
    }
}

void VolumeHistoryComponent::mouseDown (const juce::MouseEvent& event)
{
    const auto p = event.getPosition();

    // [TIMECODE-USER] Right-click on time ruler to set/reset user timecode offset
    if (getTimeRulerArea().contains (p) && event.mods.isPopupMenu())
    {
        juce::PopupMenu m;
        m.addItem (1, "Timecode offset: -2.0 s (Cubase common preroll)");
        m.addItem (2, "Timecode offset: 0.0 s (reset)");
        m.addSeparator();
        m.addItem (3, "Offset -1.0 s");
        m.addItem (4, "Offset +1.0 s");
        m.addItem (5, "Offset -0.1 s");
        m.addItem (6, "Offset +0.1 s");

        m.showMenuAsync (juce::PopupMenu::Options(),
                         [this] (int r)
                         {
                             auto setOff = [this] (double s)
                             {
                                 processor.setUserTimecodeOffsetSeconds (s);
                             };

                             auto addOff = [this] (double ds)
                             {
                                 processor.setUserTimecodeOffsetSeconds (processor.getUserTimecodeOffsetSeconds() + ds);
                             };

                             if (r == 1) { setOff (-2.0); repaint(); return; }
                             if (r == 2) { setOff ( 0.0); repaint(); return; }
                                     if (r == 3) { addOff (-1.0); repaint(); return; }
                             if (r == 4) { addOff (+1.0); repaint(); return; }
                             if (r == 5) { addOff (-0.1); repaint(); return; }
                             if (r == 6) { addOff (+0.1); repaint(); return; }
                         });
        return;
    }

    // [UI-RULERS] Drag time ruler to pan time
    if (getTimeRulerArea().contains (p))
    {
        dragMode = DragMode::timeRuler;
        dragStartPos = p;
        dragStartViewRightFrame = viewRightFrame;
        followRightEdge = false;
        followButton.setToggleState (false, juce::dontSendNotification);
        return;
    }

    // [UI-RULERS] Drag dB ruler to pan dB
    if (getDbRulerArea().contains (p))
    {
        dragMode = DragMode::dbRuler;
        dragStartPos = p;
        dragStartViewTopDb = viewTopDb;
        return;
    }

    // Existing behavior: Shift-click toggles bands; plain click toggles lines.
    if (event.mods.isShiftDown())
        showBands = ! showBands;
    else
        showLines = ! showLines;

    repaint();
}

void VolumeHistoryComponent::mouseDrag (const juce::MouseEvent& event)
{
    if (dragMode == DragMode::none)
        return;

    const auto p = event.getPosition();

    if (dragMode == DragMode::timeRuler)
    {
        const int dx = p.x - dragStartPos.x;

        // "Grab" behavior: drag right -> earlier time (content moves right)
        viewRightFrame = dragStartViewRightFrame - (double) dx / (zoomX > 1.0e-12 ? zoomX : 1.0e-12);
        clampViewRightFrame (getWidth());
        repaint();
        return;
    }

    if (dragMode == DragMode::dbRuler)
    {
        const int dy = p.y - dragStartPos.y;

        const double effectiveRange = (double) baseDbRange / zoomY;
        const double dbPerPixel = effectiveRange / (double) juce::jmax (1, getHeight());

        // [DB-DRAG-DIRECTION] drag down -> curves go down (topDb increases)
        viewTopDb = dragStartViewTopDb + (double) dy * dbPerPixel;

        const double topMin = viewMinDbLimit + effectiveRange;
        const double topMax = viewMaxDbLimit;
        viewTopDb = juce::jlimit (topMin, topMax, viewTopDb);

        markStaticBackgroundDirty();
        repaint();
        return;
    }
}

void VolumeHistoryComponent::mouseUp (const juce::MouseEvent&)
{
    dragMode = DragMode::none;
}
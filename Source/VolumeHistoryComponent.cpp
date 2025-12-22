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

        const int pendingFrames = getPendingFramesAtLevel (level);

        double numerator = maxFramesVisible - (double) pendingFrames;
        if (numerator < 0.0)
            numerator = 0.0;

        const int predicted = (int) std::floor (numerator / (double) spanFrames) + 1;
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

    int maxGroupsByX = 0;
    {
        double numerator = maxFramesVisible - (double) pendingFrames;
        if (numerator < 0.0)
            numerator = 0.0;

        const double maxGroupsFloat = numerator / (double) spanFrames;
        maxGroupsByX = (int) std::floor (maxGroupsFloat) + 1;
    }

    if (maxGroupsByX <= 0)
        return;

    const int groupsToUse = juce::jlimit (0, availableGroups, maxGroupsByX);
    if (groupsToUse <= 0)
        return;

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

    for (int chunkChronoIndex = 0; chunkChronoIndex < outCount; ++chunkChronoIndex)
    {
        const int baseGroupsAgo = (outCount - 1 - chunkChronoIndex) * step;
        const int endGroupsAgo  = juce::jmin (groupsToUse - 1, baseGroupsAgo + step - 1);

        FrameGroup agg;
        agg.momentaryMinDb =  std::numeric_limits<float>::infinity();
        agg.momentaryMaxDb = -std::numeric_limits<float>::infinity();
        agg.shortTermMinDb =  std::numeric_limits<float>::infinity();
        agg.shortTermMaxDb = -std::numeric_limits<float>::infinity();

        for (int ga = baseGroupsAgo; ga <= endGroupsAgo; ++ga)
        {
            const FrameGroup gg = getGroupAgo (levelIndex, ga);

            agg.momentaryMinDb = std::min (agg.momentaryMinDb, gg.momentaryMinDb);
            agg.momentaryMaxDb = std::max (agg.momentaryMaxDb, gg.momentaryMaxDb);
            agg.shortTermMinDb = std::min (agg.shortTermMinDb, gg.shortTermMinDb);
            agg.shortTermMaxDb = std::max (agg.shortTermMaxDb, gg.shortTermMaxDb);
        }

        const double centerGroupsAgo = 0.5 * (double) (baseGroupsAgo + endGroupsAgo);
        const double framesAgoD = (double) pendingFrames + centerGroupsAgo * (double) spanFrames;
        const int framesAgoI = (int) std::round (framesAgoD);

        outGroups[(size_t) chunkChronoIndex] = agg;
        outFramesAgo[(size_t) chunkChronoIndex] = framesAgoI;
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
    lineRenderMode = (lineRenderMode + 1) % 3; // 0->1->2->0
}

bool VolumeHistoryComponent::shouldUsePolylineForLines (int selectedLevel) const noexcept
{
    if (lineRenderMode == 2) return true;   // Force Polyline
    if (lineRenderMode == 1) return false;  // Force Stroke
    return selectedLevel >= coarseLevelStartForPolyline; // Auto
}

//==============================================================================
// Polyline drawing  [POLYLINE-DRAW] + [POLYLINE-PEAK] + [POLYLINE-SNAP]
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

    // Reserve ~width points (one per pixel column), bounded.
    const int wInt = (int) std::round (width);
    outPoints.reserve ((size_t) juce::jlimit (128, 4096, wInt + 64));

    // [POLYLINE-PEAK] We aggregate all samples that fall into the same pixel column
    // and keep a min/max Y for that column. Then we choose a representative Y that
    // preserves peaks instead of "latest wins".
    int   currentXPix = std::numeric_limits<int>::min();
    float colYMin = 0.0f;
    float colYMax = 0.0f;
    bool  haveCol = false;

    int   lastEmittedXPix = std::numeric_limits<int>::min();
    float prevY = 0.0f;
    bool  havePrev = false;

    auto emitColumn = [&]()
    {
        if (! haveCol)
            return;

        // Choose representative Y for this pixel column.
        float yRep = 0.5f * (colYMin + colYMax);

        if (havePrev)
        {
            // Preserve peaks: if one extreme is far from prevY, jump to it.
            const float dMin = std::abs (colYMin - prevY);
            const float dMax = std::abs (colYMax - prevY);
            const float dBest = std::max (dMin, dMax);

            // Threshold in pixels: below this, stay close to previous value to avoid
            // noisy "zig-zag" in dense views. Above it, show the peak.
            constexpr float peakThresholdPx = 1.25f;

            if (dBest >= peakThresholdPx)
                yRep = (dMin >= dMax ? colYMin : colYMax);
            else
                yRep = juce::jlimit (colYMin, colYMax, prevY);
        }

        // [POLYLINE-SNAP] Snap X only (to pixel center). Keep Y subpixel to reduce moiré.
        const float xAligned = (float) currentXPix + 0.5f;

        // Optional: break long gaps (avoid drawing misleading long diagonals).
        // If the data is sparse and we jump many pixels, we simply start a new segment
        // by duplicating the point (drawPolyline will still connect, but the visual
        // effect is much smaller with small typical gaps). We keep it simple here.
        if (! outPoints.empty())
        {
            const int dx = std::abs (currentXPix - lastEmittedXPix);
            if (dx > 50)
            {
                // Start a new segment by inserting a duplicate point at the new x.
                // (A true segment break would require multi-segment drawing logic.)
                outPoints.emplace_back (xAligned, yRep);
            }
        }

        outPoints.emplace_back (xAligned, yRep);

        lastEmittedXPix = currentXPix;
        prevY = yRep;
        havePrev = true;
    };

    for (size_t i = 0; i < n; ++i)
    {
        const float x = width - (float) framesAgo[i] * (float) zoomX;
        if (x < -10.0f)
            continue;

        const int xPix = (int) std::floor (x + 0.5f);
        if (xPix < 0 || xPix > (int) width)
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

    // If our duplicate-point "gap softener" created a zero-length segment at the end,
    // it's harmless; no need to clean up.
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

    buildVisibleGroupsForLevel (selectedLevel, width, scratchVisibleGroups, scratchVisibleFramesAgo);

    const size_t n = scratchVisibleGroups.size();
    const float w = bounds.getWidth();
    const float h = bounds.getHeight();

    // Decide line method once per paint
    const bool usePolylineLines = (showLines && shouldUsePolylineForLines (selectedLevel));

    if (n >= 2)
    {
        computeRepresentativeCurves (scratchVisibleGroups, scratchRepMomentaryDb, scratchRepShortTermDb);

        // Bands batched
        scratchPathBandM.clear();
        scratchPathBandS.clear();

        // Lines:
        if (showLines && usePolylineLines)
        {
            scratchPolylinePtsM.clear();
            scratchPolylinePtsS.clear();
        }
        else
        {
            scratchPathRepM.clear();
            scratchPathRepS.clear();
        }

        bool startedRepM = false, startedRepS = false;

        const float bandRangeThresholdDb = 3.0f;

        for (size_t i = 0; i < n; ++i)
        {
            const float x = w - (float) scratchVisibleFramesAgo[i] * (float) zoomX;
            if (x < -10.0f)
                continue;

            const auto& grp = scratchVisibleGroups[i];

            const float yMM = dbToY (grp.momentaryMaxDb, h);
            const float yMm = dbToY (grp.momentaryMinDb, h);
            const float ySM = dbToY (grp.shortTermMaxDb, h);
            const float ySm = dbToY (grp.shortTermMinDb, h);

            // [BAND-PATHS]
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

            // Lines in stroke mode only (polyline points are built in one pass later)
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

        // Draw bands
        if (showBands && selectedLevel > 0)
        {
            g.setColour (juce::Colours::cyan.withMultipliedAlpha (0.6f));
            g.strokePath (scratchPathBandS, juce::PathStrokeType (1.0f));

            g.setColour (juce::Colours::limegreen.withMultipliedAlpha (0.7f));
            g.strokePath (scratchPathBandM, juce::PathStrokeType (1.2f));
        }

        // Draw lines
        if (showLines)
        {
            if (usePolylineLines)
            {
                // [POLYLINE-PEAK] + [POLYLINE-SNAP]
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
    // [RULER-XMAP] + [RULER-STABLE]
    //
    // Fix the "move left then jump back" behavior by anchoring ticks to tRight:
    //   lastTick = floor(tRight / tickStep) * tickStep
    // and iterating leftwards. This is stable over time.
    //==========================================================================

    const float rulerHeight   = 16.0f;
    const float rulerBaseY    = (float) height - 2.0f;
    const float tickTopY      = rulerBaseY - 6.0f;
    const float textTopY      = rulerBaseY - 14.0f;

    if (zoomX > 0.0 && visualFrameRate > 0.0 && totalFrames > 0)
    {
        const double tRight = (double) totalFrames / visualFrameRate;

        // History availability limited by ring buffer
        const double availableSeconds = (double) availableRaw / visualFrameRate;
        const double tEarliest = juce::jmax (0.0, tRight - availableSeconds);

        // Match curve x-mapping offset at coarse levels
        const int pendingFramesOffset = getPendingFramesAtLevel (selectedLevel);

        // Compute visible seconds based on x-mapping (no fake stretching)
        const double framesByWidth = (double) width / (zoomX > 1.0e-12 ? zoomX : 1.0e-12);
        const double effectiveFramesByWidth = juce::jmax (0.0, framesByWidth - (double) pendingFramesOffset);

        const double visibleFrames = juce::jmin ((double) availableRaw, effectiveFramesByWidth);
        const double visibleSeconds = visibleFrames / visualFrameRate;

        if (visibleSeconds > 1.0e-6)
        {
            const double tLeft = juce::jmax (tEarliest, tRight - visibleSeconds);

            // Tick spacing (unchanged policy)
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

            // [RULER-STABLE] Anchor to the right edge ("now")
            const double eps = 1.0e-9;
            const double lastTick = std::floor ((tRight + eps) / tickStep) * tickStep;

            g.setColour (juce::Colours::white);
            g.setFont (10.0f);

            for (double t = lastTick; t >= tLeft - eps; t -= tickStep)
            {
                const double framesAgo = (tRight - t) * visualFrameRate;
                const double framesAgoWithOffset = framesAgo + (double) pendingFramesOffset;
                const float x = (float) ((double) width - framesAgoWithOffset * zoomX);

                if (x < -2.0f)
                    break; // further ticks will be even more left

                if (x > (float) width + 2.0f)
                    continue;

                // Draw tick at float x to avoid rounding jitter
                g.drawLine (x, tickTopY, x, rulerBaseY, 1.0f);

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

    const char* modeName = "Auto";
    if (lineRenderMode == 1) modeName = "ForceStroke";
    if (lineRenderMode == 2) modeName = "ForcePolyline";

    juce::String info = "Level: " + juce::String (selectedLevel) +
                        " (span " + juce::String (Lsel.spanFrames) + " frames, " +
                        juce::String (spanSeconds, 3) + " s)" +
                        " | ZoomX: " + juce::String (zoomX, 4) +
                        " | ZoomY: " + juce::String (zoomY, 2) +
                        " | CapPts: " + juce::String (getMaxDrawablePoints (width)) +
                        " | Bands: " + juce::String (showBands ? "ON" : "OFF") +
                        " | Lines: " + juce::String (showLines ? "ON" : "OFF") +
                        " | LineMode: " + juce::String (modeName);

    g.drawText (info, 8, 8, (int) std::min (bounds.getWidth() - 16.0f, 900.0f), 20, juce::Justification::topLeft);
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
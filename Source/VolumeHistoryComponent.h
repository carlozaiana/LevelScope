#pragma once

#include <JuceHeader.h>
#include <vector>
#include <array>

class LevelScopeAudioProcessor;

//==============================================================================
// VolumeHistoryComponent
//
// [SECTION TAGS]
//   - [STEP1-PERF]        : repaint only on new data + reused scratch buffers
//   - [STEP2-LOD-CAP]     : cap drawable points + improved LOD selection
//   - [CACHE-STATIC]      : cached static background (grid + ruler baseline)
//   - [LINE-QUALITY]      : render mode (stroke vs polyline)
//   - [POLYLINE-PEAK]     : peak-preserving per-pixel-column selection
//   - [PIXEL-ADVANCE]     : quantize scroll to whole pixels at coarse levels
//   - [RULER-HYST]        : tickStep hysteresis (prevents step switching jitter)
//   - [RULER-STABLE]      : ticks anchored to tRight (no periodic reseed jump)
//==============================================================================

class VolumeHistoryComponent : public juce::Component,
                               private juce::Timer
{
public:
    explicit VolumeHistoryComponent (LevelScopeAudioProcessor& processor);
    ~VolumeHistoryComponent() override;

    void paint (juce::Graphics& g) override;
    void resized() override;

    void mouseWheelMove (const juce::MouseEvent& event,
                         const juce::MouseWheelDetails& wheel) override;

    void mouseDown (const juce::MouseEvent& event) override;

private:
    //==============================================================================
    // History structures
    //==============================================================================

    struct FrameGroup
    {
        float momentaryMinDb = -90.0f;
        float momentaryMaxDb = -90.0f;
        float shortTermMinDb = -90.0f;
        float shortTermMaxDb = -90.0f;
    };

    struct HistoryLevel
    {
        int                 levelIndex      = 0;
        int                 groupsPerGroup  = 1;
        int                 spanFrames      = 1;
        int                 capacity        = 0;
        std::vector<FrameGroup> groups;
        int                 writeIndex      = 0;
        juce::int64         totalGroups     = 0;

        FrameGroup          pending;
        int                 pendingCount    = 0;
    };

    //==============================================================================
    // Timer
    //==============================================================================

    void timerCallback() override;

    //==============================================================================
    // History init/update
    //==============================================================================

    void initialiseHistoryLevels();
    void resetHistoryLevels();

    bool drainProcessorFifo(); // [STEP1-PERF]
    void pushFrameToHistory (float momentaryRms, float shortTermRms);

    void writeGroupToLevel (int levelIndex, const FrameGroup& group);
    void accumulateToHigherLevels (int levelIndex, const FrameGroup& sourceGroup);

    //==============================================================================
    // History access
    //==============================================================================

    int getAvailableGroups (int levelIndex) const noexcept;
    int getPendingFramesAtLevel (int levelIndex) const noexcept;
    FrameGroup getGroupAgo (int levelIndex, int groupsAgo) const noexcept;
    juce::int64 getTotalFramesL0() const noexcept;

    //==============================================================================
    // LOD selection
    //==============================================================================

    int getMaxDrawablePoints (int widthPixels) const noexcept;
    int selectBestLevelForCurrentZoom (int widthPixels) const noexcept;

    void buildVisibleGroupsForLevel (int levelIndex,
                                     int widthPixels,
                                     std::vector<FrameGroup>& outGroups,
                                     std::vector<int>& outFramesAgo) const;

    //==============================================================================
    // Representative curve
    //==============================================================================

    void computeRepresentativeCurves (const std::vector<FrameGroup>& groups,
                                      std::vector<float>& repMomentary,
                                      std::vector<float>& repShortTerm) const;

    //==============================================================================
    // Drawing helpers
    //==============================================================================

    float dbToY (float db, float height) const noexcept;

    //==============================================================================
    // Cached background
    //==============================================================================

    void markStaticBackgroundDirty() noexcept;
    void rebuildStaticBackgroundIfNeeded();

    //==============================================================================
    // [LINE-QUALITY]
    //==============================================================================

    bool isModifierForQualityToggle (const juce::ModifierKeys& mods) const noexcept;
    void cycleLineRenderMode() noexcept;
    bool shouldUsePolylineForLines (int selectedLevel) const noexcept;

    //==============================================================================
    // [PIXEL-ADVANCE]
    // Compute a sub-pixel phase (0..1) that turns continuous scrolling into
    // "advance by whole pixels" for coarse levels.
    //
    // The trick:
    //   x = width - framesAgo*zoomX + frac(nowFrames * zoomX)
    // This replaces "-nowFrames*zoomX" with "-floor(nowFrames*zoomX)".
    //==============================================================================

    float getPixelAdvancePhasePx (juce::int64 totalFramesL0,
                                 int pendingFramesOffset,
                                 int selectedLevel,
                                 bool enablePixelAdvance) const noexcept;

    //==============================================================================
    // [POLYLINE-PEAK] + [PIXEL-ADVANCE]
    // Build a polyline bounded to ~one point per pixel column, preserving peaks.
    //==============================================================================

    void buildPolylinePoints (const std::vector<int>& framesAgo,
                              const std::vector<float>& repDb,
                              float width,
                              float height,
                              float pixelAdvancePhasePx,
                              std::vector<juce::Point<float>>& outPoints) const;

    void drawPolyline (juce::Graphics& g,
                       const std::vector<juce::Point<float>>& pts,
                       float thickness) const;

    //==============================================================================
    // [RULER-HYST]
    // Tick step hysteresis based on pixels-per-second (zoom), not visibleSeconds.
    //==============================================================================

    double getTickStepSecondsWithHysteresis (int widthPixels) noexcept;

    //==============================================================================
    // Zoom
    //==============================================================================

    void applyHorizontalZoom (float wheelDelta);
    void applyVerticalZoom   (float wheelDelta);

    //==============================================================================
    // Members
    //==============================================================================

    LevelScopeAudioProcessor& processor;

    const double visualFrameRate;
    const double historyLengthSeconds;

    const float minDb;
    const float maxDb;
    const float baseDbRange;

    int rawCapacityFrames = 0;

    static constexpr int maxLevels      = 6;
    static constexpr int groupsPerLevel = 4;
    std::array<HistoryLevel, maxLevels> levels;

    // Zoom parameters
    double zoomX      = 5.0;
    double minZoomX   = 0.0005;
    double maxZoomX   = 1.333;
    double zoomY      = 1.0;
    double minZoomY   = 0.25;
    double maxZoomY   = 4.0;

    bool hasCustomZoomX = false;

    bool showBands = true;
    bool showLines = true;

    // [LINE-QUALITY]
    int lineRenderMode = 0; // 0=Auto, 1=Force Stroke, 2=Force Polyline
    int coarseLevelStartForPolyline = 3; // Auto: polyline from this level upward

    // [PIXEL-ADVANCE]
    // Only apply pixel-advance at/above this level (keeps L0..L2 smooth/AA).
    int coarseLevelStartForPixelAdvance = 3;

    // [RULER-HYST]
    int tickStepIndex = -1; // remembers current tickStep choice

    //==============================================================================
    // Scratch buffers
    //==============================================================================

    mutable std::vector<FrameGroup> scratchVisibleGroups;
    mutable std::vector<int>        scratchVisibleFramesAgo;

    mutable std::vector<float>      scratchRepMomentaryDb;
    mutable std::vector<float>      scratchRepShortTermDb;

    // Stroke-path mode scratch
    mutable juce::Path              scratchPathRepM;
    mutable juce::Path              scratchPathRepS;

    // Bands (batched)
    mutable juce::Path              scratchPathBandM;
    mutable juce::Path              scratchPathBandS;

    // Polyline mode scratch
    mutable std::vector<juce::Point<float>> scratchPolylinePtsM;
    mutable std::vector<juce::Point<float>> scratchPolylinePtsS;

    //==============================================================================
    // Cached background
    //==============================================================================

    juce::Image cachedStaticBackground;
    bool        staticBackgroundDirty = true;

    int         cachedBgW = 0;
    int         cachedBgH = 0;
    double      cachedBgZoomY = 1.0;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VolumeHistoryComponent)
};
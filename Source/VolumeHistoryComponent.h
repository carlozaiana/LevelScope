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
//   - [RULER-XMAP]        : ruler uses same x-mapping as curves
//   - [RULER-STABLE]      : ruler ticks anchored to tRight (no periodic jump)
//   - [CACHE-STATIC]      : cached static background (grid + ruler baseline)
//   - [BAND-PATHS]        : batch band segments into 2 paths
//   - [LINE-QUALITY]      : render mode (stroke vs polyline)
//   - [POLYLINE-DRAW]     : polyline draw for coarse levels
//   - [POLYLINE-PEAK]     : peak-preserving per-pixel-column selection
//   - [POLYLINE-SNAP]     : snap X only (reduce moiré), keep Y subpixel
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
    // History data structures
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
    // Line quality / render mode
    //==============================================================================

    bool isModifierForQualityToggle (const juce::ModifierKeys& mods) const noexcept;
    void cycleLineRenderMode() noexcept;

    bool shouldUsePolylineForLines (int selectedLevel) const noexcept;

    //==============================================================================
    // Polyline drawing (cheap)
    //==============================================================================

    // [POLYLINE-PEAK] + [POLYLINE-SNAP]
    // Build a polyline bounded to ~one point per pixel column, while preserving
    // peaks that would otherwise disappear under decimation.
    void buildPolylinePoints (const std::vector<int>& framesAgo,
                              const std::vector<float>& repDb,
                              float width,
                              float height,
                              std::vector<juce::Point<float>>& outPoints) const;

    void drawPolyline (juce::Graphics& g,
                       const std::vector<juce::Point<float>>& pts,
                       float thickness) const;

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

    // [BAND-PATHS]
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
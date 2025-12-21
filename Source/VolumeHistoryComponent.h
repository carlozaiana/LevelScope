#pragma once

#include <JuceHeader.h>
#include <vector>
#include <array>

class LevelScopeAudioProcessor;

//==============================================================================
// VolumeHistoryComponent
//
// [SECTION TAGS]
//   - [HISTORY-STRUCTS]
//   - [HISTORY-INIT]
//   - [HISTORY-UPDATE]
//   - [HISTORY-ACCESS]
//   - [LOD-SELECTION]
//   - [REP-LINE]
//   - [DRAW]
//   - [ZOOM]
//   - [MOUSE]
//   - [STEP1-PERF]        : persistent buffers + repaint only on new data
//   - [STEP2-LOD-CAP]     : cap drawable points + improved LOD selection
//   - [RULER-XMAP]        : ruler uses same x-mapping as curves
//   - [CACHE-STATIC]      : cached static background (grid + ruler baseline)
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
    // [HISTORY-STRUCTS]
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
        int                 spanFrames      = 1;   // RAW frames covered by one group at this level
        int                 capacity        = 0;
        std::vector<FrameGroup> groups;           // ring-buffer
        int                 writeIndex      = 0;
        juce::int64         totalGroups     = 0;

        FrameGroup          pending;
        int                 pendingCount    = 0;
    };

    //==============================================================================
    // [TIMER]
    //==============================================================================

    void timerCallback() override;

    //==============================================================================
    // [HISTORY-INIT]
    //==============================================================================

    void initialiseHistoryLevels();
    void resetHistoryLevels();

    //==============================================================================
    // [HISTORY-UPDATE]
    //==============================================================================

    // [STEP1-PERF] Return true if we actually read any frames this tick.
    bool drainProcessorFifo();

    void pushFrameToHistory (float momentaryRms, float shortTermRms);

    void writeGroupToLevel (int levelIndex, const FrameGroup& group);
    void accumulateToHigherLevels (int levelIndex, const FrameGroup& sourceGroup);

    //==============================================================================
    // [HISTORY-ACCESS]
    //==============================================================================

    int getAvailableGroups (int levelIndex) const noexcept;
    int getPendingFramesAtLevel (int levelIndex) const noexcept;
    FrameGroup getGroupAgo (int levelIndex, int groupsAgo) const noexcept;
    juce::int64 getTotalFramesL0() const noexcept;

    //==============================================================================
    // [LOD-SELECTION]
    //==============================================================================

    // [STEP2-LOD-CAP] Maximum number of drawable x-samples (points/groups).
    int getMaxDrawablePoints (int widthPixels) const noexcept;

    // [STEP2-LOD-CAP] Choose the lowest (most detailed) level that keeps the
    // predicted drawable count <= getMaxDrawablePoints(widthPixels).
    int selectBestLevelForCurrentZoom (int widthPixels) const noexcept;

    // Build visible groups (chronological order) for a given level.
    // [STEP2-LOD-CAP] If visible group count is too high, aggregate multiple
    // groups into one to keep output size bounded.
    void buildVisibleGroupsForLevel (int levelIndex,
                                     int widthPixels,
                                     std::vector<FrameGroup>& outGroups,
                                     std::vector<int>& outFramesAgo) const;

    //==============================================================================
    // [REP-LINE]
    //==============================================================================

    void computeRepresentativeCurves (const std::vector<FrameGroup>& groups,
                                      std::vector<float>& repMomentary,
                                      std::vector<float>& repShortTerm) const;

    //==============================================================================
    // [DRAW]
    //==============================================================================

    float dbToY (float db, float height) const noexcept;

    //==============================================================================
    // [CACHE-STATIC]
    //==============================================================================

    void markStaticBackgroundDirty() noexcept;
    void rebuildStaticBackgroundIfNeeded();

    //==============================================================================
    // [ZOOM]
    //==============================================================================

    void applyHorizontalZoom (float wheelDelta);
    void applyVerticalZoom   (float wheelDelta);

    //==============================================================================
    // Member variables
    //==============================================================================

    LevelScopeAudioProcessor& processor;

    const double visualFrameRate;
    const double historyLengthSeconds;

    const float minDb;
    const float maxDb;
    const float baseDbRange;

    int rawCapacityFrames = 0;

    static constexpr int maxLevels      = 6; // L0..L5
    static constexpr int groupsPerLevel = 4;
    std::array<HistoryLevel, maxLevels> levels;

    // Zoom parameters
    double zoomX      = 5.0;
    double minZoomX   = 0.0005;
    double maxZoomX   = 1.333;
    double zoomY      = 1.0;
    double minZoomY   = 0.25;
    double maxZoomY   = 4.0;

    bool   hasCustomZoomX = false;

    bool showBands = true;
    bool showLines = true;

    //==============================================================================
    // [STEP1-PERF] scratch buffers reused every repaint
    //==============================================================================

    mutable std::vector<FrameGroup> scratchVisibleGroups;
    mutable std::vector<int>        scratchVisibleFramesAgo;

    mutable std::vector<float>      scratchRepMomentaryDb;
    mutable std::vector<float>      scratchRepShortTermDb;

    mutable juce::Path              scratchPathRepM;
    mutable juce::Path              scratchPathRepS;

    //==============================================================================
    // [CACHE-STATIC]
    //   Cached background layer (background fill + horizontal dB grid lines +
    //   ruler baseline). Rebuilt only on resize or vertical zoom changes.
    //==============================================================================

    juce::Image cachedStaticBackground;
    bool        staticBackgroundDirty = true;

    int         cachedBgW = 0;
    int         cachedBgH = 0;
    double      cachedBgZoomY = 1.0;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VolumeHistoryComponent)
};
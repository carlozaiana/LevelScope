#pragma once

#include <JuceHeader.h>
#include <vector>
#include <array>

class LevelScopeAudioProcessor;

//==============================================================================
// VolumeHistoryComponent
//
// [OVERVIEW]
//   - Stores multi-level min/max loudness history (3h at 60 Hz).
//   - Level 0: RAW frames (1 frame per loudness-frame).
//   - Level N>0: groups of 4 previous-level groups (multi-resolution pyramid).
//   - For drawing, selects the best level based on zoom.
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

    // [STEP2-LOD-CAP] Maximum number of drawable "points" we allow for lines/bands.
    // Keep this roughly proportional to pixel width so work stays bounded.
    int getMaxDrawablePoints (int widthPixels) const noexcept;

    // [STEP2-LOD-CAP] Choose the lowest (most detailed) level that keeps the
    // predicted drawable count <= getMaxDrawablePoints(widthPixels).
    int selectBestLevelForCurrentZoom (int widthPixels) const noexcept;

    // Build visible groups (chronological order) for a given level.
    // [STEP2-LOD-CAP] If visible group count is too high, this function will
    // aggregate multiple groups into one to keep output size bounded.
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
    // [STEP1-PERF] scratch buffers
    //==============================================================================

    mutable std::vector<FrameGroup> scratchVisibleGroups;
    mutable std::vector<int>        scratchVisibleFramesAgo;

    mutable std::vector<float>      scratchRepMomentaryDb;
    mutable std::vector<float>      scratchRepShortTermDb;

    mutable juce::Path              scratchPathRepM;
    mutable juce::Path              scratchPathRepS;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VolumeHistoryComponent)
};
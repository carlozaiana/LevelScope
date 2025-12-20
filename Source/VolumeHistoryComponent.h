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
//   - For drawing, selects the best level based on zoom, so that the number
//     of groups to draw is ~ proportional to pixel width (fast).
//   - At each level, draws:
//       * vertical bar from min->max (envelope band),
//       * a representative line that stays inside the band and tends to
//         hug the top on upward trends and the bottom on downward trends.
//
// [SECTION TAGS]
//   - [HISTORY-STRUCTS]      : data structures for levels & groups
//   - [HISTORY-INIT]         : initialisation of levels
//   - [HISTORY-UPDATE]       : pushing new frames into history
//   - [HISTORY-ACCESS]       : helpers to read from levels
//   - [LOD-SELECTION]        : choose which level to draw based on zoom
//   - [REP-LINE]             : compute representative line inside min/max band
//   - [DRAW]                 : painting logic
//   - [ZOOM]                 : horizontal/vertical zoom handling
//   - [MOUSE]                : mouse wheel & clicks
//   - [STEP1-PERF]           : persistent buffers + repaint only on new data
//==============================================================================

class VolumeHistoryComponent : public juce::Component,
                               private juce::Timer
{
public:
    explicit VolumeHistoryComponent (LevelScopeAudioProcessor& processor);
    ~VolumeHistoryComponent() override;

    // juce::Component
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
        int                 levelIndex      = 0;   // 0 = RAW
        int                 groupsPerGroup  = 1;   // N previous-level groups per group at this level (>=1)
        int                 spanFrames      = 1;   // how many RAW (L0) frames per group at this level
        int                 capacity        = 0;   // number of groups stored (ring-buffer)
        std::vector<FrameGroup> groups;           // ring-buffer of groups
        int                 writeIndex      = 0;   // next write position in ring-buffer
        juce::int64         totalGroups     = 0;   // total groups ever written (monotonic)

        FrameGroup          pending;              // aggregator for current not-yet-final group
        int                 pendingCount    = 0;   // how many previous-level groups in 'pending'
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
    //   - Drain processor FIFO, convert RMS to dB, push as RAW L0 frames.
    //   - Propagate groups up the multi-resolution pyramid.
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

    int selectBestLevelForCurrentZoom() const noexcept;

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

    static constexpr int maxLevels              = 6; // L0..L5
    static constexpr int groupsPerLevel         = 4; // each level groups 4 previous-level groups
    std::array<HistoryLevel, maxLevels> levels;

    // Zoom parameters
    double zoomX      = 5.0;     // pixels per RAW frame
    double minZoomX   = 0.0005;
    double maxZoomX   = 1.333;
    double zoomY      = 1.0;
    double minZoomY   = 0.25;
    double maxZoomY   = 4.0;

    bool   hasCustomZoomX = false;

    bool showBands = true;
    bool showLines = true;

    //==============================================================================
    // [STEP1-PERF]
    //   Persistent buffers reused every repaint to avoid heap churn in paint().
    //   These are scratch buffers only (not state you need to save).
    //==============================================================================

    mutable std::vector<FrameGroup> scratchVisibleGroups;
    mutable std::vector<int>        scratchVisibleFramesAgo;

    mutable std::vector<float>      scratchRepMomentaryDb;
    mutable std::vector<float>      scratchRepShortTermDb;

    mutable juce::Path              scratchPathRepM;
    mutable juce::Path              scratchPathRepS;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VolumeHistoryComponent)
};
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
//   - This gives a consistent style across all zooms.
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

    void drainProcessorFifo();
    void pushFrameToHistory (float momentaryRms, float shortTermRms);

    void writeGroupToLevel (int levelIndex, const FrameGroup& group);
    void accumulateToHigherLevels (int levelIndex, const FrameGroup& sourceGroup);

    //==============================================================================
    // [HISTORY-ACCESS]
    //==============================================================================

    // Get the number of valid (written) groups at a given level.
    int getAvailableGroups (int levelIndex) const noexcept;

    // Get pending (partial) frames at this level (in RAW frames).
    int getPendingFramesAtLevel (int levelIndex) const noexcept;

    // Get group 'groupsAgo' back from the most recent group at a level.
    FrameGroup getGroupAgo (int levelIndex, int groupsAgo) const noexcept;

    // Convenience: total RAW frames ever written (L0 == RAW).
    juce::int64 getTotalFramesL0() const noexcept;

    //==============================================================================
    // [LOD-SELECTION]
    //==============================================================================

    // Choose the best level-of-detail for current zoomX.
    int selectBestLevelForCurrentZoom() const noexcept;

    // Build visible groups (in chronological order) for a given level,
    // along with their framesAgo (distance from "now" in RAW frames).
    void buildVisibleGroupsForLevel (int levelIndex,
                                     int widthPixels,
                                     std::vector<FrameGroup>& outGroups,
                                     std::vector<int>& outFramesAgo) const;

    //==============================================================================
    // [REP-LINE]
    //==============================================================================

    // Compute representative values (in dB) inside the [min, max] band, per group,
    // for both momentary and short-term curves.
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
    // [MOUSE]
    //==============================================================================

    // (mouse handlers declared in public section)

    //==============================================================================
    // Member variables
    //==============================================================================

    LevelScopeAudioProcessor& processor;

    // Loudness frame rate (frames per second, from processor: 60 Hz).
    const double visualFrameRate;

    // Total history length in seconds for RAW history (3 hours).
    const double historyLengthSeconds;

    // dB range
    const float minDb;
    const float maxDb;
    const float baseDbRange;

    // RAW capacity (number of L0 frames stored to cover 3 hours)
    int rawCapacityFrames = 0;

    // Multi-level history pyramid
    static constexpr int maxLevels              = 6; // L0..L5
    static constexpr int groupsPerLevel         = 4; // each level groups 4 previous-level groups
    std::array<HistoryLevel, maxLevels> levels;      // levels[0] = RAW

    // Zoom parameters
    double zoomX      = 5.0;     // pixels per RAW frame
    double minZoomX   = 0.0005;
    double maxZoomX   = 1.333;   // cap as discussed
    double zoomY      = 1.0;     // vertical zoom in dB
    double minZoomY   = 0.25;
    double maxZoomY   = 4.0;

    bool   hasCustomZoomX = false;

    // Band & line toggles:
    //  - showBands: draw vertical min->max bars
    //  - showLines: draw representative lines
    bool showBands = true;
    bool showLines = true;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VolumeHistoryComponent)
};
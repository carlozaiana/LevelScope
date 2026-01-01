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
//   - [BAND-PATHS]        : batch band segments into 2 paths
//   - [POLYLINE-PEAK]     : peak-preserving per-pixel-column selection
//   - [PIXEL-ADVANCE-FIX] : pixel advance by quantizing X (no wrap/jump-back)
//   - [RULER-FRAMES]      : ruler ticks computed in integer frames (stable)
//   - [RULER-HYST-FIX]    : tickStep hysteresis that doesn't thrash while zooming
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
        int levelIndex = 0;
        int spanFrames = 1;   // how many L0 frames one group covers at this level
        int capacity   = 0;   // number of groups stored at this level (ring)

        std::vector<FrameGroup> groups;

        // [TIMEBASE-PLAYHEAD] absolute group index tag per slot; -1 means empty
        std::vector<juce::int64> absGroupIndexTag;
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
    // [TIMEBASE-PLAYHEAD] projectSamplePos comes from the host playhead (in samples)
    void pushFrameToHistory (float momentaryRms, float shortTermRms,
                             juce::int64 projectSamplePos, int isPlaying);

    // [TIMEBASE-PLAYHEAD] overwrite-safe writing
    void writeGroupAbs (int levelIndex, juce::int64 absGroupIndex, const FrameGroup& group);
    bool readGroupAbs (int levelIndex, juce::int64 absGroupIndex, FrameGroup& out) const noexcept;
    void recomputeGroupAbsFromChildren (int levelIndex, juce::int64 absGroupIndex);

    //==============================================================================
    // History access
    //==============================================================================

    int getAvailableGroups (int levelIndex) const noexcept;
    FrameGroup getGroupAgo (int levelIndex, int groupsAgo) const noexcept;
    juce::int64 getTotalFramesL0() const noexcept;

    //==============================================================================
    // LOD selection
    //==============================================================================

    int getMaxDrawablePoints (int widthPixels) const noexcept;
    int selectBestLevelForCurrentZoom (int widthPixels) const noexcept;

    // [TIMEBASE-FIX] build visible groups + absolute frame index per point
    void buildVisibleGroupsForLevel (int levelIndex,
                                 int widthPixels,
                                 std::vector<FrameGroup>& outGroups,
                                 std::vector<juce::int64>& outEndFrameIndex) const;

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

    // [PIXEL-ADVANCE-FIX] quantize x to pixel centers in polyline mode
    float quantizeXToPixelCenter (float x) const noexcept;

    //==============================================================================
    // Cached background
    //==============================================================================

    void markStaticBackgroundDirty() noexcept;
    void rebuildStaticBackgroundIfNeeded();

    //==============================================================================
    // [RULER-HYST-FIX]
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

    // [TIMEBASE-PLAYHEAD]
    int         frameSamples = 0;        // samples per 60 Hz loudness frame (from processor)
    juce::int64 nowFrameIndex = 0;       // absolute 60 Hz frame index on the DAW timeline
    bool        haveNowFrameIndex = false;

    // [RULER-HYST-FIX]
    int tickStepIndex = -1; // remembered tickStep choice (hysteresis)

    //==============================================================================
    // Scratch buffers
    //==============================================================================

    mutable std::vector<FrameGroup>  scratchVisibleGroups;
    mutable std::vector<juce::int64> scratchVisibleEndFrameIndex; // [TIMEBASE-FIX]

    mutable std::vector<float>      scratchRepMomentaryDb;
    mutable std::vector<float>      scratchRepShortTermDb;

    // Stroke-path mode scratch
    mutable juce::Path              scratchPathRepM;
    mutable juce::Path              scratchPathRepS;

    // Bands (batched)
    mutable juce::Path              scratchPathBandM;
    mutable juce::Path              scratchPathBandS;

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
#pragma once

#include <JuceHeader.h>
#include <vector>
#include <array>
// [BEGIN MTDM-THRESH-UI-INCLUDE-ATOMIC]
#include <atomic>
// [END MTDM-THRESH-UI-INCLUDE-ATOMIC]

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

// [BEGIN MTDM-THRESH-UI-APVTS-LISTENER-ASYNCUPDATER-INHERIT]
class VolumeHistoryComponent : public juce::Component,
                               private juce::Timer,
                               private juce::AudioProcessorValueTreeState::Listener,
                               private juce::AsyncUpdater
// [END MTDM-THRESH-UI-APVTS-LISTENER-ASYNCUPDATER-INHERIT]
{
public:
    explicit VolumeHistoryComponent (LevelScopeAudioProcessor& processor);
    ~VolumeHistoryComponent() override;

    void paint (juce::Graphics& g) override;
    void resized() override;

    void mouseWheelMove (const juce::MouseEvent& event,
                         const juce::MouseWheelDetails& wheel) override;

    // [BEGIN UI4A1-RESIZE-CURSOR-MOUSEMOVE]
    void mouseMove (const juce::MouseEvent& event) override;
    void mouseExit (const juce::MouseEvent& event) override;
    // [END UI4A1-RESIZE-CURSOR-MOUSEMOVE]

    void mouseDown (const juce::MouseEvent& event) override;
    void mouseDrag (const juce::MouseEvent& event) override;
    void mouseUp (const juce::MouseEvent& event) override;
    void mouseDoubleClick (const juce::MouseEvent& event) override;
    // [BEGIN UI3C-FOLLOW-PUBLIC-API]
    void setFollowRightEdge (bool shouldFollow);
    bool getFollowRightEdge() const noexcept { return followRightEdge; }
    // [END UI3C-FOLLOW-PUBLIC-API]
    // [BEGIN UI3A-CURVE-VISIBILITY-PUBLIC-API]
    // Curve visibility (controlled by Mission Control top strip)
    void setShowMomentaryCurve (bool b);
    void setShowShortTermCurve (bool b);
    void setShowGateCurve      (bool b);
    void setShowRollingLraLane (bool b);

    bool getShowMomentaryCurve() const noexcept { return showMomentaryCurve; }
    bool getShowShortTermCurve() const noexcept { return showShortTermCurve; }
    bool getShowGateCurve() const noexcept      { return showGate; }
    bool getShowRollingLraLane() const noexcept { return showRollingLra; }
    // [END UI3A-CURVE-VISIBILITY-PUBLIC-API]

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

        // [LRAG] LRA relative gate curve (LUFS): typically IntegratedRunning - 20 LU
        float gateMinDb = -200.0f;
        float gateMaxDb = -200.0f;
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
    // [BEGIN MTDM-THRESH-UI-APVTS-LISTENER-DECL]
    void parameterChanged (const juce::String& parameterID, float newValue) override;
    // [END MTDM-THRESH-UI-APVTS-LISTENER-DECL]

    // [BEGIN MTDM-THRESH-UI-APVTS-LISTENER-ASYNCUPDATER-DECL]
    void handleAsyncUpdate() override;

    // Set from parameterChanged() (possibly audio thread), cleared on message thread.
    std::atomic<bool> threshUiNeedsRepaint { false };
    // [END MTDM-THRESH-UI-APVTS-LISTENER-ASYNCUPDATER-DECL]
    
    //==============================================================================
    // [STATE-PERSIST] bootstrap GUI history from processor timeline after project load
    void bootstrapHistoryFromProcessorIfNeeded();
    static juce::int64 floorDivInt64 (juce::int64 a, juce::int64 b) noexcept;
    bool bootstrappedFromProcessor = false;

    //==============================================================================
    // History init/update
    //==============================================================================

    void initialiseHistoryLevels();
    void resetHistoryLevels();

    bool drainProcessorFifo(); // [STEP1-PERF]
    // [TIMEBASE-PLAYHEAD] projectSamplePos comes from the host playhead (in samples)
    void pushFrameToHistory (float momentaryVal, float shortTermVal, float gateLufs,
                             juce::int64 frameIndex60Hz, int isPlaying);

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

    // [VIEW-NAV] zoom/pan tied to mouse
    void applyHorizontalZoom (float wheelDelta, float anchorX);
    void applyVerticalZoom   (float wheelDelta, float anchorY);
    void panTime             (float wheelDelta);
    void panDb               (float wheelDelta);
    void clampViewRightFrame (int widthPixels) noexcept;

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

    // [VIEW-NAV-Y-LIMITS] allow positive dBFS and extra space below -90
    static constexpr double viewMinDbLimit = -180.0; // lowest visible dB label allowed
    static constexpr double viewMaxDbLimit =  24.0; // highest visible dB label allowed

    bool hasCustomZoomX = false;

    bool showBands = true;
    bool showLines = true;
    // [BEGIN UI3B-RIGHT-STRIP-CONSTANTS]
    // LUFS scale itself is narrow; meters consume the remaining right strip width.
    // [BEGIN UI3C2-DBSCALE-WIDTH-SHRINK]
    static constexpr int dbScaleWidthPx = 40;
    // [END UI3C2-DBSCALE-WIDTH-SHRINK]

    // User-resizable total right strip width (scale + meters), anchored to the right edge.
    int rightStripWidthPxUser = dbScaleWidthPx + 120;

    static constexpr int rightStripMinWidthPx = dbScaleWidthPx + 72;
    static constexpr int rightStripMaxWidthPx = dbScaleWidthPx + 260;

    int dragStartRightStripWidthPxUser = dbScaleWidthPx + 120;
    // [END UI3B-RIGHT-STRIP-CONSTANTS]
    // [BEGIN UI3A-CURVE-VISIBILITY-FLAGS]
    bool showMomentaryCurve = true;
    bool showShortTermCurve = true;
    // [END UI3A-CURVE-VISIBILITY-FLAGS]

    //==============================================================================
    // [VIEW-NAV] View state
    //==============================================================================

    // Right edge of the visible window in 60 Hz timeline frames (can be fractional).
    double viewRightFrame = 0.0;

    // If true, view follows furthest-written (nowFrameIndex) automatically.
    bool followRightEdge = true;

    // Vertical view: top of the visible dB range (allows Y-zoom/pan)
    double viewTopDb = 0.0;

    bool autoRefollowArmed = true; // [AUTO-FOLLOW-HYST]

    //==============================================================================
    // [UI-RULERS] Ruler hit zones + interactions
    //==============================================================================

    // [BEGIN UI3B-RIGHT-STRIP-AREAS]
    juce::Rectangle<int> getTimeRulerArea() const;

    // Total right strip (LUFS scale + meters) used to reserve space from plot.
    juce::Rectangle<int> getDbRulerArea() const;

    // LUFS scale sub-area (left part of the right strip) used for Y pan/reset interactions.
    juce::Rectangle<int> getDbScaleArea() const;

    // Meter sub-area (right part of the right strip).
    juce::Rectangle<int> getRightMetersArea() const;
    // [END UI3B-RIGHT-STRIP-AREAS]

    void resetXViewDefault();
    void fitXViewMaxZoomOut();
    void resetYViewDefault();

    // Dragging rulers
    // [BEGIN ROLLING-LRA-SPLITTER-DRAGMODE]
    // [BEGIN UI3B-RIGHTSTRIP-DRAGMODE]
    enum class DragMode { none, timeRuler, dbRuler, rollingLraDivider, rightStripResizer };
    // [END UI3B-RIGHTSTRIP-DRAGMODE]
    // [END ROLLING-LRA-SPLITTER-DRAGMODE]
    DragMode dragMode = DragMode::none;

    juce::Point<int> dragStartPos;
    double dragStartViewRightFrame = 0.0;
    double dragStartViewTopDb = 0.0;

    // Follow toggle UI
    juce::ToggleButton followButton;

    juce::ToggleButton gateButton;   // [LRAG] show/hide LRA gate curve
    bool showGate = false;

    juce::ToggleButton rollingLraButton; // [ROLLING-LRA] show/hide rolling LRA curve
    bool showRollingLra = false;
    // [BEGIN UI-METERS-IDLE-DECAY-STATE]
    struct RightStripMeterDisplay
    {
        // I/O (dBFS)
        float inPeakDbCurrent  = -200.0f;
        float inPeakDbHold     = -200.0f;
        float inRmsDbCurrent   = -200.0f;

        float outPeakDbCurrent = -200.0f;
        float outPeakDbHold    = -200.0f;
        float outRmsDbCurrent  = -200.0f;

        // Stages (dB)
        float upBoostDbCurrent = 0.0f;
        float upBoostDbHold    = 0.0f;

        float downGrDbHold     = 0.0f;
        float limGrDbHold      = 0.0f;
    };

    RightStripMeterDisplay meterDisp;

    // Detect whether processor is still being called (host time-in-samples changes)
    juce::int64 lastSeenHostSamplesForMeters = 0;
    double      lastHostSamplesChangeMs      = 0.0;

    // [BEGIN UI-METERS-IDLE-DECAY-LASTMS]
    double lastMeterUiUpdateMs = 0.0;
    // [END UI-METERS-IDLE-DECAY-LASTMS]

    // UI-only decay when host stops calling processBlock()
    static constexpr float uiDecayDbPerSecLevel = 36.0f; // dBFS meters fall fairly quickly
    static constexpr float uiDecayDbPerSecStage = 18.0f; // GR/Boost falls a bit slower

    // [BEGIN UI-METERS-IDLE-DECAY-SIGNATURE]
    // Returns true if the meters changed (or are decaying) and need a repaint.
    // [BEGIN UI-METERS-IDLE-DECAY-SIGNATURE2]
    bool updateRightStripMetersFromTimer (bool gotNewHistoryData);
    // [END UI-METERS-IDLE-DECAY-SIGNATURE2]
    // [END UI-METERS-IDLE-DECAY-SIGNATURE]
    // [END UI-METERS-IDLE-DECAY-STATE]
    // [BEGIN ROLLING-LRA-SPLITTER-STATE]
    // Rolling LRA lane height (px). User can drag the divider to resize.
    int rollingLaneHeightPx = 46; // old hardcoded value
    static constexpr int rollingLaneMinHeightPx = 28;
    static constexpr int rollingLaneMaxHeightPx = 240;

    int dragStartRollingLaneHeightPx = 46;

    // Main plot area (curves + threshold overlay). Updated in paint and in hit-test prep.
    juce::Rectangle<float> mainPlotArea;
    // [END ROLLING-LRA-SPLITTER-STATE]

    // [BEGIN MTDM-THRESH-UI-DECL]
    //==============================================================================
    // [MTDM-THRESH-UI] Multi-threshold dynamics overlay (T0..T3)
    // UI-only, driven by APVTS params. O(1) draw: 4 lines + 4 handles.
    //==============================================================================

    struct ThresholdHandle
    {
        juce::Rectangle<float> hitBounds;   // expanded hit zone
        juce::Rectangle<float> drawBounds;  // visible handle rect
        int index = -1;                     // 0..3
    };

    // Cached parameter pointers (initialised lazily on UI thread)
    bool mtdmParamPtrsInitialised = false;
    std::array<std::atomic<float>*, 4>        mtdmThreshAtoms  { { nullptr, nullptr, nullptr, nullptr } };
    std::array<juce::RangedAudioParameter*, 4> mtdmThreshParams { { nullptr, nullptr, nullptr, nullptr } };

    // Drag state
    bool thresholdDragging = false;
    int  activeThresholdIndex = -1;

    // Gesture state: begin once (mouseDown / lazy for pushed), end once (mouseUp)
    std::array<bool, 4> threshGestureActive { { false, false, false, false } };

    // Geometry cache for drawing/hit-test (updated O(1) whenever needed)
    juce::Rectangle<float> mtdmOverlayArea;            // excludes right dB ruler only
    std::array<ThresholdHandle, 4> thresholdHandles;   // T0..T3

    // Constants
    static constexpr float kThreshMinGapLu = 0.1f;

    // Helpers
    void initMtdmParamPointersIfNeeded();
    bool mtdmParamsAvailable() const noexcept;

    void updateMtdmThresholdOverlayGeometry();               // O(1)
    void drawMtdmThresholdOverlay (juce::Graphics& g);       // O(1)

    float yToLufs (float y, float height) const noexcept;

    void computeOrderedThresholdsWithPush (int changedIndex,
                                          float newValueLufs,
                                          float outVals[4]) const noexcept;

    void applyThresholdValuesDuringDrag (const float newVals[4]);

    void endAllThresholdGestures();

    // Mouse handlers (called from existing mouseDown/Drag/Up)
    void handleThresholdMouseDown (const juce::MouseEvent& event);
    void handleThresholdMouseDrag (const juce::MouseEvent& event);
    void handleThresholdMouseUp   (const juce::MouseEvent& event);
    // [END MTDM-THRESH-UI-DECL]

    //==============================================================================
    // [ROLLING-LRA] Per-second short-term + gate history and computed rolling LRA
    // We compute short-term samples at 1 Hz by sampling the incoming 60 Hz stream.
    // Then compute rolling LRA over the last N seconds (30/60/120), gated by:
    // gate = max(-70 LUFS, gateLufsAtCurrentSecond) where gateLufsAtCurrentSecond is the
    // yellow LRAG curve (IntegratedRunning - 20).
    //==============================================================================

    int secondsCapacity = 0; // ~10800 for 3 hours

    std::vector<float> secShortTermLufs;   // per-second sampled short-term LUFS
    std::vector<float> secGateLufs;        // per-second sampled LRA gate (LUFS)
    std::vector<float> secRollingLraLu;    // computed rolling LRA (LU)
    std::vector<juce::int64> secAbsIndexTag; // -1 = empty, else absSecondIndex

    juce::int64 lastSecondPushed = std::numeric_limits<juce::int64>::min();

    // Rolling window selection (read from processor)
    int rollingWindowSecondsCached = 60;

    // Rebuild state (when user changes 30/60/120, we rebuild rolling curve over stored seconds)
    bool        rollingRebuildInProgress = false;
    juce::int64 rollingRebuildMinSecond = 0;
    juce::int64 rollingRebuildMaxSecond = -1;
    juce::int64 rollingRebuildNextSecond = 0;

    // [ROLLING-LRA] ring helpers
    int wrapSecondSlot (juce::int64 absSecondIndex) const noexcept;

    // [ROLLING-LRA] write one per-second sample (overwrite-safe via tag)
    void pushSecondSample (juce::int64 absSecondIndex,
                           float shortTermLufs,
                           float gateLufs);

    // [ROLLING-LRA] recompute rolling LRA for one second index (if present)
    void recomputeRollingLraForSecond (juce::int64 absSecondIndex);

    // [ROLLING-LRA] detect selector change and (re)start rebuild
    void startRollingRebuildIfWindowChanged();

    // [ROLLING-LRA] chunked rebuild so we don’t block UI
    void rebuildRollingLraStep (int maxSecondsToProcess);

    // [TIMEBASE-PLAYHEAD]
    int         frameSamples = 0;            // samples per 60 Hz loudness frame (from processor)

    // Furthest-written timeline frame index (monotonic while running).
    // Used as "latest" for what data exists, so loops/rewinds don't truncate history.
    juce::int64 nowFrameIndex = 0;
    bool        haveNowFrameIndex = false;

    // Actual DAW playhead position on the same 60 Hz frame grid (can move backwards on loop/seek).
    juce::int64 playheadFrameIndex = 0;
    bool        havePlayheadFrameIndex = false;

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

    mutable juce::Path scratchPathGate; // [LRAG]

    mutable juce::Path scratchPathRollingLra; // [ROLLING-LRA]

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
    double      cachedBgTopDb = 0.0; // [VIEW-NAV]

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VolumeHistoryComponent)
};
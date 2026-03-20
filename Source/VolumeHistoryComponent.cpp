#include "VolumeHistoryComponent.h"
#include "PluginProcessor.h"
// [BEGIN MTDM-THRESH-UI-APVTS-LISTENER-INCLUDE]
#include "Core/Processing/Modules/MultiThresholdDynamicsParamIDs.h"
// [END MTDM-THRESH-UI-APVTS-LISTENER-INCLUDE]

#include <cmath>
#include <algorithm>
#include <limits>

//==============================================================================
// Constructor / Destructor
//==============================================================================

VolumeHistoryComponent::VolumeHistoryComponent (LevelScopeAudioProcessor& proc)
    : processor (proc),
      visualFrameRate (proc.getLoudnessFrameRate()),
      historyLengthSeconds (3.0 * 3600.0),
      minDb (-70.0f),
      maxDb (-10.0f),
      baseDbRange (60.0f)
{
    jassert (visualFrameRate > 0.0);

    initialiseHistoryLevels();
    resetHistoryLevels();

    // [ROLLING-LRA] allocate per-second buffers (3 hours)
    secondsCapacity = (int) std::ceil (historyLengthSeconds);
    secondsCapacity = juce::jmax (1, secondsCapacity);

    secShortTermLufs.assign ((size_t) secondsCapacity, -200.0f);
    secGateLufs.assign      ((size_t) secondsCapacity, -200.0f);
    secRollingLraLu.assign  ((size_t) secondsCapacity, 0.0f);
    secAbsIndexTag.assign   ((size_t) secondsCapacity, (juce::int64) -1);

    // Cache current rolling window selection from processor (default is 60)
    rollingWindowSecondsCached = processor.getRollingLraWindowSeconds();

    // [VIEW-NAV]
    viewTopDb = (double) maxDb;     // top starts at 0 dBFS
    viewRightFrame = 0.0;           // will follow right edge once we have data
    followRightEdge = true;

    // [BEGIN LS-UIST-APPLY-HISTORY-CTOR]
    applyPersistedUiStateFromProcessor();
    // [END LS-UIST-APPLY-HISTORY-CTOR]

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
        if (! followRightEdge)
           autoRefollowArmed = false; // [AUTO-FOLLOW-HYST] prevent instant re-follow

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

    // [LRAG] Gate toggle
    gateButton.setButtonText ("Gate");
    gateButton.setToggleState (showGate, juce::dontSendNotification);
    gateButton.onClick = [this]
    {
        showGate = gateButton.getToggleState();
        repaint();
    };
    addAndMakeVisible (gateButton);

    // [ROLLING-LRA] Rolling LRA toggle
    rollingLraButton.setButtonText ("rLRA");
    rollingLraButton.setToggleState (showRollingLra, juce::dontSendNotification);
    rollingLraButton.onClick = [this]
    {
        showRollingLra = rollingLraButton.getToggleState();
        repaint();
    };
    addAndMakeVisible (rollingLraButton);

    // [BEGIN UI3C-HIDE-HISTORY-TOPRIGHT-TOGGLES]
    // These controls are now owned by the Mission Control strip.
    followButton.setVisible (false);
    gateButton.setVisible (false);
    rollingLraButton.setVisible (false);
    // [END UI3C-HIDE-HISTORY-TOPRIGHT-TOGGLES]

    markStaticBackgroundDirty();
    // [BEGIN MTDM-THRESH-UI-APVTS-LISTENER-REGISTER]
    {
        auto& apvts = processor.getAPVTS();
        using namespace levelscope::mtdm::ParamIDs;

        apvts.addParameterListener (t0Lufs, this);
        apvts.addParameterListener (t1Lufs, this);
        apvts.addParameterListener (t2Lufs, this);
        apvts.addParameterListener (t3Lufs, this);
    }
    // [END MTDM-THRESH-UI-APVTS-LISTENER-REGISTER]
}

// [BEGIN MTDM-THRESH-UI-APVTS-LISTENER-UNREGISTER]
VolumeHistoryComponent::~VolumeHistoryComponent()
{
    {
        auto& apvts = processor.getAPVTS();
        using namespace levelscope::mtdm::ParamIDs;

        apvts.removeParameterListener (t0Lufs, this);
        apvts.removeParameterListener (t1Lufs, this);
        apvts.removeParameterListener (t2Lufs, this);
        apvts.removeParameterListener (t3Lufs, this);
    }

    stopTimer();
}
// [END MTDM-THRESH-UI-APVTS-LISTENER-UNREGISTER]

// [BEGIN UI3A-HISTORY-RELOAD-FROM-PROCESSOR-IMPL]
void VolumeHistoryComponent::applyPersistedUiStateFromProcessor()
{
    const auto ui = processor.getPersistedUIStateSnapshot();

    rightStripWidthPxUser = juce::jlimit (rightStripMinWidthPx, rightStripMaxWidthPx, ui.rightStripWidthPxUser);
    rollingLaneHeightPx   = juce::jlimit (rollingLaneMinHeightPx, rollingLaneMaxHeightPx, ui.rollingLaneHeightPx);

    dragStartRightStripWidthPxUser = rightStripWidthPxUser;
    dragStartRollingLaneHeightPx   = rollingLaneHeightPx;

    showMomentaryCurve = ui.showMomentaryCurve;
    showShortTermCurve = ui.showShortTermCurve;
    showGate           = ui.showGateCurve;
    showRollingLra     = ui.showRollingLra;
    followRightEdge    = ui.followRightEdge;
    showLines          = (showMomentaryCurve || showShortTermCurve);

    followButton.setToggleState (followRightEdge, juce::dontSendNotification);
    gateButton.setToggleState   (showGate, juce::dontSendNotification);
    rollingLraButton.setToggleState (showRollingLra, juce::dontSendNotification);

    processor.setRollingLraWindowSeconds (ui.rollingWindowSeconds);
    rollingWindowSecondsCached = processor.getRollingLraWindowSeconds();

    if (ui.historyViewStateValid)
    {
        zoomX = juce::jlimit (minZoomX, maxZoomX, ui.historyZoomX);
        zoomY = juce::jlimit (minZoomY, maxZoomY, ui.historyZoomY);

        const double maxVisibleRange = viewMaxDbLimit - viewMinDbLimit;
        const double minZoomYByLimits = (double) baseDbRange / maxVisibleRange;
        zoomY = juce::jmax (zoomY, minZoomYByLimits);

        viewTopDb = ui.historyViewTopDb;
        {
            const double effectiveRange = (double) baseDbRange / zoomY;
            const double topMin = viewMinDbLimit + effectiveRange;
            const double topMax = viewMaxDbLimit;
            viewTopDb = juce::jlimit (topMin, topMax, viewTopDb);
        }

        viewRightFrame = ui.historyViewRightFrame;
        hasCustomZoomX = ui.historyHasCustomZoomX;
    }

    markStaticBackgroundDirty();
}

void VolumeHistoryComponent::syncPersistedUiStateToProcessor() const
{
    processor.setUiHistoryToggleState (showMomentaryCurve,
                                       showShortTermCurve,
                                       showGate,
                                       showRollingLra,
                                       followRightEdge);

    processor.setUiHistoryViewState (zoomX,
                                     zoomY,
                                     viewRightFrame,
                                     viewTopDb,
                                     hasCustomZoomX,
                                     true);
}

void VolumeHistoryComponent::reloadFromProcessorState()
{
    applyPersistedUiStateFromProcessor();

    resetHistoryLevels();

    std::fill (secShortTermLufs.begin(),  secShortTermLufs.end(),  -200.0f);
    std::fill (secGateLufs.begin(),       secGateLufs.end(),       -200.0f);
    std::fill (secRollingLraLu.begin(),   secRollingLraLu.end(),   0.0f);
    std::fill (secAbsIndexTag.begin(),    secAbsIndexTag.end(),    (juce::int64) -1);

    lastSecondPushed = std::numeric_limits<juce::int64>::min();

    rollingWindowSecondsCached = processor.getRollingLraWindowSeconds();
    rollingRebuildInProgress = false;
    rollingRebuildMinSecond = 0;
    rollingRebuildMaxSecond = -1;
    rollingRebuildNextSecond = 0;

    nowFrameIndex = 0;
    haveNowFrameIndex = false;

    playheadFrameIndex = 0;
    havePlayheadFrameIndex = false;

    bootstrappedFromProcessor = false;
    bootstrapHistoryFromProcessorIfNeeded();

    if (followRightEdge && haveNowFrameIndex)
    {
        viewRightFrame = (double) nowFrameIndex;
        clampViewRightFrame (getWidth());
    }

    resized();
    markStaticBackgroundDirty();
    repaint();
}
// [END UI3A-HISTORY-RELOAD-FROM-PROCESSOR-IMPL]
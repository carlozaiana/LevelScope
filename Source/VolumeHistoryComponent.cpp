#include "VolumeHistoryComponent.h"
#include "PluginProcessor.h"

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

    markStaticBackgroundDirty();
}

VolumeHistoryComponent::~VolumeHistoryComponent()
{
    stopTimer();
}
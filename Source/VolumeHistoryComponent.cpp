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
      minDb (-110.0f),
      maxDb (  0.0f),
      baseDbRange (std::abs (minDb))
{
    jassert (visualFrameRate > 0.0);

    initialiseHistoryLevels();
    resetHistoryLevels();

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

    markStaticBackgroundDirty();
}

VolumeHistoryComponent::~VolumeHistoryComponent()
{
    stopTimer();
}
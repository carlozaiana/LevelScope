#pragma once

#include "StandaloneAudition.h"
#include "StandaloneSessionModel.h"
#include <juce_gui_basics/juce_gui_basics.h>

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-AUDIO-PANEL]
class StandaloneAudioPanel final : public juce::Component
{
public:
    StandaloneAudioPanel (StandaloneSessionModel&, StandaloneAudition&, std::function<void()> changed);
    void refresh();
    void sourceChanged();
    void synchronize();
    void resized() override;
private:
    StandaloneSessionModel& session;
    StandaloneAudition& audio;
    std::function<void()> onChanged;
    juce::ComboBox typeBox, deviceBox, domainBox, layoutPreset;
    juce::TextButton scanButton { "Scan devices" }, openButton { "Open device" }, closeButton { "Close device" };
    juce::TextButton prepareButton { "Prepare audition" };
    juce::Label auditionLabel;
    juce::TextButton playButton { "Play selected" }, stopButton { "Stop" }, confirmButton { "Confirm roles" };
    juce::TextButton restartButton { "Restart" }, seekButton { "Seek" };
    juce::TextEditor seekPosition;
    juce::Slider gain;
    juce::Label statusLabel, heading;
    juce::Viewport viewport;
    juce::Component rows;
    std::array<juce::Label, maxChannels> channelLabels;
    std::array<juce::ComboBox, maxChannels> roleBoxes, outputBoxes;
    juce::String actionError;
    void scan();
    void rebuildOutputs();
};
// [END LS-STANDALONE-AUDIO-PANEL]
}

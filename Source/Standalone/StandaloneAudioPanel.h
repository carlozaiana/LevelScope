#pragma once

#include "StandaloneAudition.h"
#include "RememberedChannelRoles.h"
#include "StandaloneSessionModel.h"
#include <juce_gui_basics/juce_gui_basics.h>

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-AUDIO-PANEL]
class StandaloneAudioPanel final : public juce::Component
{
public:
    StandaloneAudioPanel (StandaloneSessionModel&, StandaloneAudition&, std::function<void()> changed,juce::File roleDirectory={});
    void refresh();
    void sourceChanged(bool freshImport=false);
    std::function<void()> onAnalyse;
    void synchronize();
    void resized() override;
private:
    StandaloneSessionModel& session;
    StandaloneAudition& audio;
    std::function<void()> onChanged;
    RememberedChannelRoles roleMemory;
    bool reviewing=false,editingOrder=false;
    juce::Label suggestionLabel;
    juce::TextButton editOrder{"Edit order"};
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

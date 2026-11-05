#pragma once

#include "StandaloneAudition.h"
#include "RememberedChannelRoles.h"
#include "OutputDeviceSettings.h"
#include "StandaloneSessionModel.h"
#include <juce_gui_basics/juce_gui_basics.h>

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-AUDIO-PANEL]
class StandaloneAudioPanel final : public juce::Component
{
public:
    StandaloneAudioPanel (StandaloneSessionModel&, StandaloneAudition&, std::function<void()> changed,juce::File roleDirectory={},
                          std::shared_ptr<OutputDeviceCatalog> = {}, OutputDeviceSettings = OutputDeviceSettings());
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
    OutputDeviceConnection connection;
    uint64_t connectionRevision = 0;
    bool reviewing=false,editingOrder=false;
    juce::Label suggestionLabel;
    juce::TextButton editOrder{"Edit order"};
    juce::ComboBox deviceBox, domainBox, layoutPreset;
    juce::ToggleButton automaticReconnect { "Automatically reconnect my last output device" };
    juce::TextButton scanButton { "Scan devices" }, openButton { "Connect" }, closeButton { "Disconnect" };
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
    void refreshDeviceChoices();
    void connectionChanged();
    void rebuildOutputs();
};
// [END LS-STANDALONE-AUDIO-PANEL]
}

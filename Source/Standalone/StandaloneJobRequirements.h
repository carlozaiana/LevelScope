#pragma once
#include "StandaloneSessionModel.h"
#include <juce_gui_basics/juce_gui_basics.h>

namespace levelscope::standalone
{
std::optional<proposal::JobRequirements> loadRequirementsPreset(const juce::File&,juce::String&);
bool saveRequirementsPreset(const proposal::JobRequirements&,const juce::File&,juce::String&,
                            const std::function<bool(const char*)>& fail = {});
// Pinned measured comparison and the complete editable target definition.
class JobRequirementsPanel final : public juce::Component
{
public:
    JobRequirementsPanel(StandaloneSessionModel&,std::function<void()>,juce::File presetFolder = {});
    void resized() override;
    void paint(juce::Graphics&) override;
    void refresh();
    void refreshAvailability();
    juce::Label& informationLine(){return status;}
    bool commitFields();
private:
    StandaloneSessionModel& session;std::function<void()> changed;
    std::array<juce::TextEditor,5> fields;
    std::array<juce::Label,5> headers;
    std::array<juce::Label,10> measured;
    std::array<juce::Label,3> rowNames;
    juce::ComboBox profiles;
    juce::TextButton save{"Save as preset..."};
    juce::Label status;
    juce::String runtime,activeDigest;
    int selectedIndex=-1;
    bool editing=false, updating=false;
    juce::File presetDirectory;
    std::vector<proposal::JobRequirements> presets;
    std::unique_ptr<juce::FileChooser> chooser;
    std::unique_ptr<juce::AlertWindow> saveDialog;
    void message(const juce::String&);
    void rebuildProfiles();
    void selectProfile();
    void savePreset();
    void loadPreset();
    int rowHeight() const;
};
}

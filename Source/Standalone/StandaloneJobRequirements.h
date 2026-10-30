#pragma once
#include "StandaloneSessionModel.h"
#include <juce_gui_basics/juce_gui_basics.h>

namespace levelscope::standalone
{
std::optional<proposal::JobRequirements> loadRequirementsPreset(const juce::File&,juce::String&);
bool saveRequirementsPreset(const proposal::JobRequirements&,const juce::File&,juce::String&,
                            const std::function<bool(const char*)>& fail = {});
class JobRequirementsPanel final : public juce::Component
{
public:
    JobRequirementsPanel(StandaloneSessionModel&,std::function<void()>);
    void resized() override;
    void refresh();
private:
    StandaloneSessionModel& session;std::function<void()> changed;
    juce::ToggleButton enabled{"Job requires a maximum LRA"};
    juce::TextEditor ceiling,name,provenance;
    juce::Label units,nameLabel,provenanceLabel,status;
    juce::TextButton use{"Use EBU requirements"},save{"Save preset..."},load{"Load preset..."};
    juce::String runtime,activeDigest;
    std::unique_ptr<juce::FileChooser> chooser;
    void message(const juce::String&);
    void confirm();
    void choose(bool saving);
};
}

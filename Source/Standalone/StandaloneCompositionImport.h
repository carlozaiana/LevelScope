#pragma once
#include "ProgrammeComposition.h"
#include "StandaloneProgrammeImport.h"

namespace levelscope::standalone
{
class CompositionImportPanel final : public juce::Component, private juce::Timer
{
public:
    using Completion=std::function<void(std::optional<ProgrammeComposition>)>;
    CompositionImportPanel(std::optional<ProgrammeComposition>,Completion);
    ~CompositionImportPanel() override;
    void resized() override;
    void cancel();
private:
    ProgrammeComposition recipe;
    Completion complete;
    juce::ComboBox preset;
    juce::Label explanation,notice;
    juce::TextButton accept{"Commit confirmed composition"},cancelButton{"Cancel"};
    juce::ToggleButton confirmation{"I confirm these are complementary, correctly balanced stems intended to be combined, with exact sample-zero alignment and the displayed roles."};
    std::array<juce::TextButton,2> choose,groups,canonical;
    std::array<juce::Label,2> descriptions;
    std::array<std::array<juce::ComboBox,6>,2> roles;
    std::array<std::array<juce::Label,6>,2> channelNames;
    std::array<juce::TextEditor,2> languages,versions;
    std::unique_ptr<juce::FileChooser> chooser;
    std::unique_ptr<ProgrammeImportPanel> groupPanel;
    std::unique_ptr<juce::DocumentWindow> groupWindow;
    std::shared_ptr<JobControl> control;
    std::future<MediaValidationResult> pending;
    juce::File pendingFile;
    int pendingComponent=0;
    bool closing=false;
    juce::String operationNotice;
    ProgrammeComponent& part(int i){return i==0?recipe.me:recipe.dialogue;}
    std::vector<ChannelRole> expectedRoles(int i) const;
    void changed();
    void refresh();
    void pickFile(int);
    void pickGroup(int);
    void timerCallback() override;
};
}

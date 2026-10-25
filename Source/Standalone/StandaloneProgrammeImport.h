#pragma once
#include "ProgrammeSource.h"
#include <juce_gui_basics/juce_gui_basics.h>
#include <future>

namespace levelscope::standalone
{
class ProgrammeImportPanel final : public juce::Component, private juce::Timer
{
public:
    using Completion=std::function<void(std::optional<ProgrammeAssembly>)>;
    explicit ProgrammeImportPanel(Completion);
    ~ProgrammeImportPanel() override;
    void resized() override;
    void cancel();
private:
    Completion complete;
    juce::ComboBox layout;
    juce::TextButton choose{"Choose mono files..."}, accept{"Import complete programme"}, cancelButton{"Cancel"};
    juce::ToggleButton alignment{"I confirm these exact files form the same aligned mix; sample zero is shared."};
    juce::Label notice, explanation;
    std::array<juce::Label,6> labels;
    std::array<juce::ComboBox,6> roles;
    std::unique_ptr<juce::FileChooser> chooser;
    std::shared_ptr<JobControl> control;
    std::future<ProgrammeInspection> pending;
    std::optional<ProgrammeAssembly> group;
    bool closing=false;
    void refresh();
    void timerCallback() override;
};
}

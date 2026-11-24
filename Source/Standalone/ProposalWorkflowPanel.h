#pragma once
#include "StandaloneSessionModel.h"
#include "ProposalOperation.h"
#include <juce_gui_basics/juce_gui_basics.h>
namespace levelscope::standalone
{
class ProposalWorkflowPanel final : public juce::Component
{
public:
    static constexpr int height=184;
    explicit ProposalWorkflowPanel(StandaloneSessionModel&);
    std::function<void()> onPropose,onReview,onHistory,onPreferences,onRegions,onChanged;
    std::function<void()> onRecover,onBoundaries;
    std::function<void(size_t)> onAlternative;
    void setOperation(const ProposalOperation* value){operation=value;}
    void setDiagnosticInspection(bool value) {diagnosticInspection=value;}
    void refresh(float progress=-1);
    juce::String outcomeText() const {return headline.getText()+"\n"+reason.getText();}
    juce::String actionText() const {return availability.getText();}
    void resized() override;
private:
    StandaloneSessionModel& session;
    juce::TextButton propose{"Propose"},review{"Review proposal"},apply{"Apply proposal"},discard{"Discard proposal"},history{"Last applied explanation"},preferences{"Proposal limits..."};
    juce::TextButton regions{"Method / region planning..."};
    juce::TextButton recover{"Restore before Propose"},boundaries{"Inspect region boundaries"};
    juce::ComboBox alternatives;
    const ProposalOperation* operation=nullptr;
    juce::Label method,headline,availability,startHint;
    bool diagnosticInspection=false;
    juce::TextEditor reason;
};
}

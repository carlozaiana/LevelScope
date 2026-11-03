#pragma once
#include "StandaloneSessionModel.h"
#include <juce_gui_basics/juce_gui_basics.h>
namespace levelscope::standalone
{
class ProposalWorkflowPanel final : public juce::Component
{
public:
    static constexpr int height=150;
    explicit ProposalWorkflowPanel(StandaloneSessionModel&);
    std::function<void()> onPropose,onReview,onHistory,onPreferences,onChanged;
    void refresh(float progress=-1);
    juce::String outcomeText() const {return headline.getText()+"\n"+reason.getText();}
    juce::String actionText() const {return availability.getText();}
    void resized() override;
private:
    StandaloneSessionModel& session;
    juce::TextButton propose{"Propose"},review{"Review proposal"},apply{"Apply proposal"},discard{"Discard proposal"},history{"Last applied explanation"},preferences{"Proposal limits..."};
    juce::Label headline,availability;
    juce::TextEditor reason;
};
}

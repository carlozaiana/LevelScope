#pragma once
#include "StandaloneRegionProposal.h"
#include "StandaloneProposalPreferences.h"

namespace levelscope::standalone
{
class RegionProposalReviewPanel final : public juce::Component
{
public:
    RegionProposalReviewPanel(std::shared_ptr<const SectionProposalResult>,std::function<juce::String()>);
    std::function<void(size_t,int64_t,int64_t)> onInspect;
    void resized() override;
    void paint(juce::Graphics& g) override {g.fillAll(juce::Colour(0xff17212a));}
private:
    std::shared_ptr<const SectionProposalResult> result;
    std::function<juce::String()> statusText;
    juce::ComboBox attempts,regions;
    juce::Label label,authority;
    juce::TextButton inspect{"Inspect in timeline"};
    std::unique_ptr<ProposalInspectionPanel> report;
    void chooseAttempt();
    void inspectSelection();
};
}

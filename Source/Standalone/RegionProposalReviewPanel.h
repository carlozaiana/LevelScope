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
    // At most one complete explanation plus three immutable measured attempts.
    // Keep each bounded page view, including its current page, when switching.
    std::vector<std::unique_ptr<ProposalInspectionPanel>> reports;
    ProposalInspectionPanel* report=nullptr;
    juce::String attemptText(size_t) const;
    void chooseAttempt();
    void inspectSelection();
};
}

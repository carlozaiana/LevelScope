#pragma once
#include "TimelineAppearance.h"
#include "StandaloneRegionProposal.h"
#include "StandaloneProposalPreferences.h"

namespace levelscope::standalone
{
class RegionProposalReviewPanel final : public juce::Component, private juce::Timer
{
public:
    static constexpr size_t pageBytes=1500; // dense numeric histories on actual High Sierra
    static constexpr int sharedEvidenceId=1001,proxyComparisonsId=1002;
    RegionProposalReviewPanel(std::shared_ptr<const SectionProposalResult>,std::function<juce::String()>);
    std::function<void(size_t,int64_t,int64_t)> onInspect;
    using Selection=std::function<std::shared_ptr<const SectionProposalResult>(const juce::String&,const juce::String&)>;
    void configureMeasuredSelection(bool qualified,std::function<juce::String()> blocked,Selection select);
    void refreshSelection();
    void resized() override;
    void paint(juce::Graphics& g) override {g.fillAll(palette::panel());}
private:
    std::shared_ptr<const SectionProposalResult> result;
    std::function<juce::String()> statusText;
    juce::ComboBox attempts,regions;
    juce::Label label,authority;
    juce::TextButton inspect{"Inspect in timeline"};
    juce::TextButton useAlternative{"Use this measured alternative"};
    juce::TextEditor coverage;
    bool selectionQualified=false;
    std::function<juce::String()> selectionBlocked;
    Selection selectMeasured;
    juce::String selectionNotice;
    // Complete explanation, three immutable attempts and two documentary views.
    // Keep each bounded page view, including its current page, when switching.
    std::vector<std::unique_ptr<ProposalInspectionPanel>> reports;
    ProposalInspectionPanel* report=nullptr;
    juce::String attemptText(size_t) const;
    juce::String sharedEvidenceText() const,proxyComparisonsText() const;
    void chooseAttempt();
    void inspectSelection();
    void useMeasuredAlternative();
    void timerCallback() override {if(isShowing())refreshSelection();}
};
}

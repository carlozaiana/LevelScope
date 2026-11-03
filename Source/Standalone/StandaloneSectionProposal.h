#pragma once
#include "StandaloneProposal.h"
#include "VerifiedSource.h"
#include "Core/Proposal/Methods/SectionTrimMethod.h"

namespace levelscope::standalone
{
struct SectionProposalRequest
{
    ServiceContext sourceContext, currentContext;
    proposal::ProposalPolicy policy;
    std::shared_ptr<const AnalysisSnapshot> source, current;
    proposal::Applicability applicability=proposal::Applicability::unconfirmed;
    proposal::SectionMap annotations, appliedAnnotations;
    bool pendingManualDraft=false;
    std::shared_ptr<const proposal::SectionPeakMap> attribution;
    proposal::TrimRequest scalar() const;
    juce::String binding() const;
};
enum class SectionProposalStatus { unavailable, noChange, globalOnly, noSolution, ready };
enum class SectionProposalIssue { none, inputs, unsupportedRequirements, preferences, candidateRejected };
// Only the service constructs accepted results. Callers get immutable access and
// must still compare the request binding before Apply. Never session-cache this.
class SectionProposalResult
{
public:
    SectionProposalStatus status() const { return state; }
    const juce::String& identity() const { return bound; }
    const juce::String& text() const { return explanation; }
    const juce::String& narrativeText() const { return narrative; }
    const juce::String& outcomeTitle() const { return title; }
    const juce::String& outcomeReason() const { return reason; }
    SectionProposalIssue issue() const { return problem; }
    const proposal::SectionMap& originalAnnotations() const { return annotations; }
    std::optional<ProcessingPlan> previewPlan() const;
    bool previewHasFinalGlobal() const { return state==SectionProposalStatus::ready||state==SectionProposalStatus::globalOnly||bool(candidate); }
    const ProcessingPlan& plan() const { return proposed; }
    const proposal::SectionTrimSeed& selection() const { return seed; }
    const proposal::TrimProposal& global() const { return scalarResult; }
    const std::shared_ptr<const AnalysisSnapshot>& assessment() const { return candidate; }
    const std::shared_ptr<const proposal::SectionPeakMap>& peaks() const { return attribution; }
    int passes() const { return passCount; }
private:
    SectionProposalResult()=default;
    SectionProposalStatus state=SectionProposalStatus::unavailable;
    juce::String bound,explanation,narrative,title,reason;
    SectionProposalIssue problem=SectionProposalIssue::inputs;
    proposal::SectionMap annotations;
    ProcessingPlan proposed;
    proposal::SectionTrimSeed seed;
    proposal::TrimProposal scalarResult;
    std::shared_ptr<const AnalysisSnapshot> candidate;
    std::shared_ptr<const proposal::SectionPeakMap> attribution;
    int passCount=0;
    friend ServiceResult<std::shared_ptr<const SectionProposalResult>> generateSectionProposal(
        const SectionProposalRequest&,JobControl&,const SnapshotStorage&);
};
using SectionProposalOutcome=ServiceResult<std::shared_ptr<const SectionProposalResult>>;
SectionProposalOutcome generateSectionProposal(const SectionProposalRequest&,JobControl&,const SnapshotStorage& = {});
}

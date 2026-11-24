#pragma once
#include "StandaloneProposal.h"
#include "VerifiedSource.h"
#include "Core/Proposal/Methods/SectionTrimMethod.h"
#include "Core/Proposal/LraPassages.h"
#include "Core/Proposal/RegionPlanning.h"

namespace levelscope::standalone
{
struct RegionProposalReview;
class RegionAssessmentReceipt;
struct SectionProposalRequest
{
    ServiceContext sourceContext, currentContext;
    proposal::ProposalPolicy policy;
    std::shared_ptr<const AnalysisSnapshot> source, current;
    proposal::Applicability applicability=proposal::Applicability::unconfirmed;
    proposal::SectionMap annotations, appliedAnnotations;
    uint64_t sectionPlanningGeneration=0; // Runtime-only authority for editorial geometry edits.
    bool pendingManualDraft=false;
    std::shared_ptr<const proposal::SectionPeakMap> attribution;
    std::optional<proposal::RegionPlanningSettings> regionSettings;
    proposal::TrimRequest scalar() const;
    juce::String binding() const;
};
enum class SectionProposalStatus { unavailable, noChange, globalOnly, noSolution, ready, reviewRequired };
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
    const juce::String& methodVersion() const { return method; }
    const std::optional<proposal::LraTrimSelection>& lraSelection() const { return lraSeed; }
    const std::shared_ptr<const RegionProposalReview>& regionReview()const {return regions;}
    SectionProposalIssue issue() const { return problem; }
    const proposal::SectionMap& originalAnnotations() const { return annotations; }
    std::optional<ProcessingPlan> previewPlan() const;
    bool previewHasFinalGlobal() const { return state==SectionProposalStatus::ready||state==SectionProposalStatus::globalOnly||bool(candidate)||bool(eventReceipt); }
    const ProcessingPlan& plan() const { return proposed; }
    const proposal::SectionTrimSeed& selection() const { return seed; }
    const proposal::TrimProposal& global() const { return scalarResult; }
    const std::shared_ptr<const AnalysisSnapshot>& assessment() const { return candidate; }
    const std::shared_ptr<const RegionAssessmentReceipt>& eventAssessment()const {return eventReceipt;}
    const juce::String& selectedAuthority()const {return selectionAuthority;}
    bool selectedFrom(const std::shared_ptr<const SectionProposalResult>& source) const {return selectionSource==source;}
    bool preservesAppliedAnnotations()const {return scalarMeasured;}
    bool explicitlySelectedAlternative()const {return userSelectedAlternative;}
    const std::shared_ptr<const proposal::SectionPeakMap>& peaks() const { return attribution; }
    int passes() const { return passCount; }
private:
    SectionProposalResult()=default;
    SectionProposalStatus state=SectionProposalStatus::unavailable;
    juce::String bound,explanation,narrative,title,reason;
    juce::String method=proposal::sectionTrimMethodVersion;
    std::optional<proposal::LraTrimSelection> lraSeed;
    std::shared_ptr<const RegionProposalReview> regions;
    SectionProposalIssue problem=SectionProposalIssue::inputs;
    proposal::SectionMap annotations;
    ProcessingPlan proposed;
    proposal::SectionTrimSeed seed;
    proposal::TrimProposal scalarResult;
    std::shared_ptr<const AnalysisSnapshot> candidate;
    std::shared_ptr<const RegionAssessmentReceipt> eventReceipt;
    juce::String selectionAuthority;
    bool userSelectedAlternative=false,scalarMeasured=false;
    std::shared_ptr<const SectionProposalResult> selectionSource;
    std::shared_ptr<const proposal::SectionPeakMap> attribution;
    int passCount=0;
    friend ServiceResult<std::shared_ptr<const SectionProposalResult>> generateSectionProposal(
        const SectionProposalRequest&,JobControl&,const SnapshotStorage&);
    friend ServiceResult<std::shared_ptr<const SectionProposalResult>> generateLraProposal(
        const SectionProposalRequest&,JobControl&,const SnapshotStorage&);
    friend ServiceResult<std::shared_ptr<const SectionProposalResult>> generateRegionProposal(
        const SectionProposalRequest&,JobControl&,const SnapshotStorage&);
    friend ServiceResult<std::shared_ptr<const SectionProposalResult>> generateEventRegionProposal(
        const SectionProposalRequest&,JobControl&,const SnapshotStorage&);
    friend ServiceResult<std::shared_ptr<const SectionProposalResult>> measureGlobalProposal(
        const SectionProposalRequest&,const std::shared_ptr<const SectionProposalResult>&,JobControl&,const SnapshotStorage&);
    friend ServiceResult<std::shared_ptr<const SectionProposalResult>> selectWorkflowProposal(
        const SectionProposalRequest&,const std::shared_ptr<const SectionProposalResult>&,const juce::String&,const juce::String&,bool);
    friend ServiceResult<std::shared_ptr<const SectionProposalResult>> selectMeasuredEventAlternative(
        const SectionProposalRequest&,const std::shared_ptr<const SectionProposalResult>&,const juce::String&,const juce::String&);
};
using SectionProposalOutcome=ServiceResult<std::shared_ptr<const SectionProposalResult>>;
// Workflow plumbing: no new search, ranking or correction strategy.
SectionProposalOutcome measureGlobalProposal(const SectionProposalRequest&,const std::shared_ptr<const SectionProposalResult>&,JobControl&,const SnapshotStorage& = {});
SectionProposalOutcome selectWorkflowProposal(const SectionProposalRequest&,const std::shared_ptr<const SectionProposalResult>&,
    const juce::String& resultId,const juce::String& attemptId,bool selectedByPropose);
juce::String measuredProposalApplyError(const SectionProposalRequest&,const SectionProposalResult&);
SectionProposalOutcome generateSectionProposal(const SectionProposalRequest&,JobControl&,const SnapshotStorage& = {});
SectionProposalOutcome generateLraProposal(const SectionProposalRequest&,JobControl&,const SnapshotStorage& = {});
}

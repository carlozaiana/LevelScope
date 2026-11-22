#pragma once
#include "StandaloneSectionProposal.h"
#include "Core/Proposal/RegionPlanning.h"
#include "Core/Proposal/EventPlanning.h"
#include <limits>

namespace levelscope::standalone
{
// Runtime-only measured authority. No public constructor or deserialization path.
// The worker's complete window populations are validated, then released; only the
// finite exact metric/peak provenance needed by selection and verification remains.
class RegionAssessmentReceipt final
{
public:
    const measurement::Result& measured()const {return detail;}
    const proposal::MeasuredProgramme& input()const {return programme;}
    const ServiceContext& context()const {return measuredContext;}
    const juce::String& snapshotIdentity()const {return snapshot;}
    const juce::String& identity()const {return id;}
    const juce::String& method()const {return measurementMethod;}
    const juce::String& content()const {return measuredContent;}
    juce::String validationError(const SectionProposalRequest&,const ProcessingPlan&)const;
    proposal::Evaluation evaluation(const SectionProposalRequest&,const ProcessingPlan&)const;
    size_t storageBytes()const;
private:
    RegionAssessmentReceipt(const AnalysisSnapshot&,const SectionProposalRequest&);
    RegionAssessmentReceipt(const RegionAssessmentReceipt&)=delete;
    RegionAssessmentReceipt& operator=(const RegionAssessmentReceipt&)=delete;
    measurement::Result detail; // window vectors MUST remain empty
    proposal::MeasuredProgramme programme;
    ServiceContext measuredContext;
    juce::String snapshot,id,binding,planIdentity,profileDigest,measurementMethod,measuredContent;
    proposal::Applicability applicability=proposal::Applicability::unconfirmed;
    friend SectionProposalOutcome generateRegionProposal(const SectionProposalRequest&,JobControl&,const SnapshotStorage&);
    friend SectionProposalOutcome generateEventRegionProposal(const SectionProposalRequest&,JobControl&,const SnapshotStorage&);
};
struct RegionAttempt
{
    proposal::RegionPlan plan;
    proposal::RegionPrediction prediction;
    proposal::Evaluation evaluation;
    // Scalar documentary values only; no second copy of full window populations.
    std::array<measurement::Metric,6> measured; // IL, TP, LRA, M, S, sample peak
    juce::String snapshotIdentity;
    double seconds=0;
    std::array<std::optional<double>,5> signedPredictionError; // IL/LRA/P10/P95/sample M
    proposal::EventCoverage coverage;
    std::shared_ptr<const RegionAssessmentReceipt> receipt;
    juce::String attemptIdentity,family;
};
struct RegionProposalReview
{
    proposal::RegionConstraints constraints;
    proposal::RegionFormation formation;
    std::vector<RegionAttempt> attempts;
    std::vector<juce::String> restrictions;
    size_t proxyEvaluations=0,duplicatePlans=0,rejectedPlans=0,retainedAttempt=std::numeric_limits<size_t>::max();
    size_t sidecarBytes=0,planningBytes=0;
    double preparationSeconds=0,sidecarSeconds=0,planningSeconds=0;
    bool searchTruncated=false;
    bool eventMode=false;
    juce::String resultIdentity;
    std::optional<size_t> recommendedAttempt,selectedAttempt;
    std::shared_ptr<const proposal::EventEvidence> eventEvidence;
    size_t featureBytes=0,featureHighWaterBytes=0;
    juce::String text()const;
};
SectionProposalOutcome generateRegionProposal(const SectionProposalRequest&,JobControl&,const SnapshotStorage& = {});
SectionProposalOutcome generateEventRegionProposal(const SectionProposalRequest&,JobControl&,const SnapshotStorage& = {});
SectionProposalOutcome selectMeasuredEventAlternative(const SectionProposalRequest&,
    const std::shared_ptr<const SectionProposalResult>&,const juce::String& resultIdentity,const juce::String& attemptIdentity);
juce::String validateEventRegionProposalForApply(const SectionProposalRequest&,const SectionProposalResult&);
juce::String validateRegionProposalForApply(const SectionProposalRequest&,const SectionProposalResult&);
}

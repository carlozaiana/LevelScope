#pragma once
#include "StandaloneSectionProposal.h"
#include "Core/Proposal/RegionPlanning.h"
#include <limits>

namespace levelscope::standalone
{
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
    juce::String text()const;
};
SectionProposalOutcome generateRegionProposal(const SectionProposalRequest&,JobControl&,const SnapshotStorage& = {});
juce::String validateRegionProposalForApply(const SectionProposalRequest&,const SectionProposalResult&);
}

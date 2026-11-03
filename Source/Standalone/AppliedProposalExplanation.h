#pragma once
#include "StandaloneSectionProposal.h"
#include "Core/Proposal/Model/JobRequirements.h"

namespace levelscope::standalone
{
// Documentary scalar results only. These cannot be installed as AnalysisSnapshots.
struct ProposalMeasurementRecord
{
    juce::String identity, content, method, build, platform, originSession, identicalSource;
    uint64_t originJob=0, originRevision=0;
    MeasurementDomain domain=MeasurementDomain::source;
    std::array<measurement::Metric,6> metrics; // IL, TP, LRA, max M, max S, sample peak
    juce::String evaluation;
};
struct AppliedProposalExplanation
{
    static constexpr size_t maximumBytes=8*1024*1024;
    juce::String programmeId, inputIdentity, binding, method, narrative;
    juce::String beforeIdentity, afterIdentity;
    ProcessingPlan before, after;
    ChannelLayout roles;int sampleRate=0;int64_t frames=0;
    proposal::TargetProfile target;
    proposal::JobRequirements requirements;
    proposal::ProposalPolicy preferences;
    proposal::Applicability applicability=proposal::Applicability::unconfirmed;
    proposal::SectionMap annotations;
    proposal::SectionTrimSeed selection;
    proposal::TrimProposal global;
    ProposalMeasurementRecord source, baseline;
    std::optional<ProposalMeasurementRecord> candidate;
    int passes=0;
    juce::String text() const;
};
juce::String proposalTargetText(const proposal::TargetProfile&);
ProposalMeasurementRecord proposalMeasurementRecord(const AnalysisSnapshot&,const proposal::TargetProfile&,proposal::Applicability);
juce::String proposalPlanIdentity(const ProcessingPlan&);
juce::String proposalMeasurementText(const ProposalMeasurementRecord&,const juce::String& label);
AppliedProposalExplanation captureProposalExplanation(const SectionProposalRequest&,const proposal::JobRequirements&,
    const juce::String& programmeId,const proposal::TrimProposal&,const SectionProposalResult* = nullptr);
// Independent, versioned session attachment. Bounds and semantic checks apply on
// both write and read. Failure never promotes history into verification authority.
std::vector<uint8_t> encodeAppliedExplanation(const AppliedProposalExplanation&);
std::shared_ptr<const AppliedProposalExplanation> decodeAppliedExplanation(const std::vector<uint8_t>&,JobControl&);
}

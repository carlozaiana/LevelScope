#pragma once
#include "RuleBasedMethod.h"
#include "Core/Proposal/Model/ProposalPolicy.h"
#include "Core/Proposal/SectionPeakMap.h"
#include "Core/Processing/SectionGainPlan.h"

namespace levelscope::proposal
{
inline constexpr const char* sectionTrimMethodVersion="section-trim-proposal-v2";
struct SectionCorrection
{
    uint64_t id=0;
    int64_t begin=0,end=0,peakEighthFrame=0;
    int channel=0;
    double originalDb=0,resultingDb=0,peakDbtp=0,deficitDb=0,requestedDb=0,reductionDb=0,reserveDb=0,uncoveredDb=0;
    bool padding=false;
};
struct SectionTrimSeed
{
    bool available=false;
    ProposalPolicy policy;
    size_t processingRegions=0;
    processing::SectionGainPlan plan;
    std::vector<SectionCorrection> corrections;
    int64_t affectedFrames=-1; // -1: extent not assessed after an earlier rejection
    juce::String reason;
};
// Immutable per-request policy. No DSP or session mutation.
SectionTrimSeed selectSectionTrims(const TrimRequest&, const SectionMap&,
    const std::optional<processing::SectionGainPlan>& existing, const SectionPeakMap&, const ProposalPolicy& = {});
std::vector<std::pair<int64_t,int64_t>> changedEnvelopeSupport(const processing::SectionGainPlan&, const processing::SectionGainPlan&);
int64_t changedEnvelopeExtent(const processing::SectionGainPlan& before,
                             const processing::SectionGainPlan& after);
juce::String sectionSeedText(const SectionTrimSeed&);
}

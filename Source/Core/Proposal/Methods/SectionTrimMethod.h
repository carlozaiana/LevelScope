#pragma once
#include "RuleBasedMethod.h"
#include "Core/Proposal/SectionPeakMap.h"
#include "Core/Processing/SectionGainPlan.h"

namespace levelscope::proposal
{
inline constexpr const char* sectionTrimMethodVersion="section-trim-proposal-v1";
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
    processing::SectionGainPlan plan;
    std::vector<SectionCorrection> corrections;
    int64_t affectedFrames=-1; // -1: extent not assessed after an earlier rejection
    juce::String reason;
};
// Fixed M13 policy. No mutable preference controls, DSP or session mutation.
SectionTrimSeed selectSectionTrims(const TrimRequest&, const SectionMap&,
    const std::optional<processing::SectionGainPlan>& existing, const SectionPeakMap&);
int64_t changedEnvelopeExtent(const processing::SectionGainPlan& before,
                             const processing::SectionGainPlan& after);
juce::String sectionSeedText(const SectionTrimSeed&);
}

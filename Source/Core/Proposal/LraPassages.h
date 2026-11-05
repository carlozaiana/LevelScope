#pragma once
#include "AnalysisEngine.h"
#include "Methods/SectionTrimMethod.h"

namespace levelscope::proposal
{
inline constexpr const char* lraPassagesVersion="lra-passages-v1";
inline constexpr const char* lraTrimMethodVersion="lra-selective-trim-v1";
struct LraPassage
{
    int64_t begin=0,end=0;
    size_t windows=0,paddedWindows=0;
    double score=0,p95=0; // LU * included windows; not an additive LRA contribution
};
struct LraPassages
{
    std::shared_ptr<const ProgramAnalysis> analysis;
    std::optional<double> ceiling,upperLevel;
    juce::String identity,notice;
    size_t included=0,excluded=0,unavailable=0;
    std::vector<LraPassage> low,high; // deterministic ranked order
    size_t storageBytes() const;
};
std::shared_ptr<const LraPassages> findLraPassages(std::shared_ptr<const ProgramAnalysis>,
    std::optional<double> ceiling,const std::function<bool()>& cancelled=[] {return false;});
struct LraRegionGain
{
    LraPassage passage;
    double requested=0,reduction=0,reserve=0,uncovered=0;
};
struct LraTrimSelection
{
    processing::SectionGainPlan before,plan;
    std::vector<LraRegionGain> selected;
    std::vector<juce::String> restrictions;
    int64_t affectedFrames=0;
    size_t processingRegions=0;
    double maximumAdditional=0;
    bool available=false;
    juce::String text() const;
};
LraTrimSelection selectLraTrims(const LraPassages&,const processing::SectionGainPlan& before,
    const ProposalPolicy&,const std::function<bool()>& cancelled=[] {return false;});
LraTrimSelection refineLraTrims(const LraTrimSelection&,double measuredDeficit,
    const ProposalPolicy&,const std::function<bool()>& cancelled=[] {return false;});
// Exact linear-segment endpoints establish full-envelope constraints, not a
// decimated plot or only plateau values. Returns a policy/method explanation.
juce::String validateLraEnvelope(const processing::SectionGainPlan& before,
    const processing::SectionGainPlan& after,const ProposalPolicy&,int64_t& affected,
    size_t& regions,double& maximumAdditional);
}

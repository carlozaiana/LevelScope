#pragma once
#include "LraPassages.h"
#include "EventFeatureDefinition.h"
#include <array>
#include <chrono>

namespace levelscope::proposal
{
inline constexpr const char* regionConstraintsVersion = "region-constraints-v1";
inline constexpr const char* regionFormationVersion = "region-formation-v1";
inline constexpr const char* regionSearchVersion = "region-gain-search-v1";
inline constexpr const char* regionProxyVersion = "window-energy-proxy-v1";
inline constexpr const char* regionRankVersion = "candidate-review-rank-v2";
struct RegionSpan
{
    int64_t begin=0,end=0; // half-open source frames, never video coordinates
    bool operator==(const RegionSpan& b)const {return begin==b.begin&&end==b.end;}
};
struct RegionEdgeRange
{
    int64_t minimum=0,maximum=0; // closed; min==max is an exact pin
    bool contains(int64_t f)const {return minimum<=f&&f<=maximum;}
};
enum class RegionMode { automatic, selectedSection, selection };
struct RegionPlanningSettings
{
    int schema=1;
    juce::String inputIdentity;
    RegionMode mode=RegionMode::automatic;
    uint64_t sectionId=0;
    RegionSpan selection;
    bool allowSplitting=true,allowRefinement=true,fillShortGap=false,allowGlobal=true,compareAlternatives=false;
    bool preferCompleteEvent=false,allowAutomaticMargins=false;
    EventFeatureDefinition eventFeatures;
    int maximumSplits=2,leftMarginMs=0,rightMarginMs=0,minimumIntervalMs=3000,separationMs=1000;
    int needThresholdTenths=10,persistenceMs=1000;
    std::optional<RegionEdgeRange> start,end;
    // Explicit internal nominal-cut constraints, ordered and nonoverlapping.
    // Pins may survive as non-audible documentary dividers after canonicalization.
    std::vector<RegionEdgeRange> internalCuts;
    juce::String protectionMethod="local-envelope-v1";
    std::vector<RegionSpan> protectedLocal;
    static RegionPlanningSettings defaults(RegionMode,const juce::String& input);
    static RegionPlanningSettings eventDefaults(RegionMode,const juce::String& input);
    const char* methodVersion()const {return schema==2?eventSearchVersion:regionSearchVersion;}
    juce::String validationError()const;
    juce::var encode()const;
    static RegionPlanningSettings decode(const juce::var&); // throws; never widens unknown permissions
    juce::String digest()const;
};
struct PlanningRoot
{
    juce::String id;
    RegionSpan hint,allowed;
    RegionEdgeRange start,end;
    std::vector<RegionEdgeRange> internalCuts;
    bool userSelected=false;
    double score=0;
    RegionSpan requestedAllowed,ownership;
    juce::String permissionNotice;
};
struct RegionConstraints
{
    RegionPlanningSettings settings;
    int sampleRate=0;
    int64_t frames=0;
    std::vector<PlanningRoot> roots;
    size_t omittedRoots=0;
    juce::String identity,error;
};
RegionConstraints resolveRegionConstraints(const RegionPlanningSettings&,const SectionMap&,
    const LraPassages&,const processing::SectionGainPlan&);
struct RegionCut {int64_t frame=0;double strength=0;juce::String reason;};
struct RegionGeometry
{
    size_t root=0;
    RegionSpan support; // requested OUTER changed-support edges; plateau cuts differ
    std::vector<int64_t> cuts; // nominal internal plateau cuts
    juce::String reason;
    double strength=0;
};
struct RegionFormation
{
    std::vector<std::vector<RegionGeometry>> alternatives;
    std::vector<std::vector<RegionCut>> nominations;
    size_t constructions=0,omittedCuts=0,omittedGeometries=0;
    bool cancelled=false;
    std::vector<juce::String> reasons;
};
RegionFormation formRegions(const RegionConstraints&,const ProgramAnalysis&,double ceiling,
    const std::function<bool()>& cancelled=[] {return false;},
    const std::vector<std::vector<RegionGeometry>>* preferred=nullptr);
struct RegionAssignment
{
    RegionGeometry geometry;
    RegionSpan nominal; // solved outer plateau cuts, before/after fixed ramps
    std::vector<double> reductions; // cuts+1; internal zero parts retained
    std::vector<double> requested,reserve,uncovered;
};
struct RegionActivity
{
    std::vector<RegionSpan> support;
    std::vector<size_t> splits;
    int64_t affectedFrames=0,longestRun=0,marginFrames=0;
    size_t intervalCount=0,transitionCount=0;
    double maximumAdditional=0,maximumCumulative=0,absoluteDbFrames=0,globalDelta=0;
};
struct RegionPlan
{
    processing::SectionGainPlan envelope;
    double globalDb=0;
    std::vector<RegionAssignment> assignments;
    RegionActivity activity;
    juce::String identity,error;
};
// Constructs with the unchanged processor, solves representable outer support edges,
// then validates the complete canonical before/after envelope. No ramp truncation.
RegionPlan constructRegionPlan(const processing::SectionGainPlan&,double baselineGlobal,
    const RegionConstraints&,const ProposalPolicy&,std::vector<RegionAssignment>,double global);
juce::String validateRegionPlan(const processing::SectionGainPlan&,double baselineGlobal,
    const RegionConstraints&,const ProposalPolicy&,RegionPlan&);
std::vector<double> regionGainChoices(const processing::SectionGainPlan&,const ProgramAnalysis&,
    RegionSpan,const ProposalPolicy&,double ceiling,bool includeLowContext=false);
struct RegionGainNeed {double headroom=0,deficit=0,requested=0,reserveRequested=0;std::vector<double> choices;};
RegionGainNeed regionGainNeed(const processing::SectionGainPlan&,const ProgramAnalysis&,
    RegionSpan,const ProposalPolicy&,double ceiling,bool includeLowContext=false);
juce::String regionPlanIdentity(const processing::SectionGainPlan&,double globalDb);

struct PredictionBin
{
    int64_t begin=0,end=0;
    double weightedEnergy=0;
    std::array<double,measurement::channelCapacity> sampleMax{};
    bool drain=false;
};
struct PredictionSidecar
{
    juce::String binding,baselineSnapshot,baselinePlan;
    int sampleRate=0,channels=0;
    int64_t frames=0;
    std::vector<measurement::Role> roles;
    std::vector<PredictionBin> bins;
    size_t storageBytes()const;
};
// Auxiliary analysis only. The caller supplies exact, already processed baseline PCM.
// Does not call or change the authoritative meter or retain programme PCM.
class PredictionCollector
{
public:
    PredictionCollector(int rate,std::vector<measurement::Role>,int64_t frames,
        juce::String binding,juce::String snapshot,juce::String plan);
    void push(const float* const*,int count);
    std::shared_ptr<const PredictionSidecar> finish();
private:
    struct Filter {double b0=1,b1=0,b2=0,a1=0,a2=0,z1=0,z2=0;double process(double);};
    struct Channel {Filter shelf,highpass;double weight=1;};
    std::vector<Channel> filters;
    std::shared_ptr<PredictionSidecar> data;
    int64_t position=0;
    bool finished=false;
    void append(double,const float* const*,int,bool);
};
struct RegionPrediction
{
    std::optional<double> integrated,lra,p10,p95,momentarySampled,shortTermSampled;
    std::optional<double> truePeakUpperBound; // unknown by default; never an actual failure
    std::optional<double> integratedGate,lraGate;
    size_t integratedEntered=0,integratedLeft=0,lraEntered=0,lraLeft=0;
    bool eofUncertain=false,gateSensitive=false,identical=false,cancelled=false;
    juce::String limitations;
};
RegionPrediction predictRegionPlan(const ProgramAnalysis&,const PredictionSidecar&,
    const processing::SectionGainPlan& before,double beforeGlobal,const RegionPlan&,
    const std::function<bool()>& cancelled=[] {return false;});
struct RegionReviewScore
{
    bool complete=false,allPass=false;
    double lraExcess=0,otherMaximum=0,otherSum=0;
    size_t failedRules=0;
    RegionActivity activity;
    juce::String id;
};
// Anchored shortlist over the WHOLE set, not a non-transitive pairwise comparator.
std::optional<size_t> selectRegionReview(const std::vector<RegionReviewScore>&,double band=.01);
RegionReviewScore regionMeasuredScore(const Evaluation&,const RegionActivity&,const juce::String&);
RegionReviewScore regionPredictedScore(const RegionPrediction&,const TargetProfile&,const RegionActivity&,const juce::String&);
}

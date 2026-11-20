#pragma once
#include "RegionPlanning.h"
#include <memory>

namespace levelscope::proposal
{
// Auxiliary descriptive evidence only. None of these values can satisfy a rule.
struct EventCell
{
    int64_t begin=0,end=0;
    double broadband=0,low=0,current=0,lfe=0;
    float rawPeak=0,fastContrast=0,lowContrast=0;
    uint8_t fastWitness=0,lowWitness=0,flags=0,reserved=0;
};
static_assert(sizeof(EventCell)==64,"Event cell resource contract changed");
struct EventGroup
{
    juce::String id;
    size_t root=0;
    RegionSpan onset,lastAttack,decay,inDecay;
    int64_t onsetPeakFrame=-1,firstNeedWindow=-1,lastNeedWindow=-1;
    size_t needWindowCount=0;int onsetPeakChannel=-1;double onsetPeak=0;
    double background=0,lowBackground=0;
    uint16_t witnesses=0;
    bool corroborated=false,needAssociated=false,endingSupported=false,incomplete=false;
    juce::String reason;
    int outerGroup=-1; // Index of the retained continuing interpretation.
    int64_t entryDecision=0, observedEnd=0;
    double entryStepDb=0,capturedExcess=0;
    bool innerNominee=false;
};
struct AttackWitness
{
    size_t root=0;
    int64_t decision=0,confirmed=0,earliest=0,predecessor1=0,predecessor2=0,innerEntry=0;
    double fastDb=0,lowDb=0,stepDb=0;
    std::optional<double> beforeMomentary,afterMomentary;
    enum State { ordinary, noPredecessors, noContext, noStep, localised, continuation, conflict, roleOnly, clipped } state=ordinary;
};
struct EventEvidence
{
    juce::String binding,sourceSnapshot,currentSnapshot,baselinePlan,constraintsIdentity,definitionDigest,identity;
    int sampleRate=0;
    int64_t frames=0;
    std::vector<measurement::Role> roles;
    std::vector<EventCell> cells;
    std::vector<RegionSpan> retainedContexts;
    std::vector<EventGroup> groups;
    std::vector<AttackWitness> attacks;
    std::array<size_t,8> summaryCounts{};
    size_t omittedWitnesses=0,conflictingEntrances=0;
    bool reserveSummary(size_t root);
    juce::String innerDescription()const,innerDigest()const;
    std::vector<std::vector<RegionCut>> needCuts;
    std::vector<bool> incompleteRoots;
    size_t reservedBytes=0,highWaterBytes=0,omittedGroups=0,omittedContexts=0;
    bool finished=false;
    size_t storageBytes()const;
};
class EventFeatureCollector
{
public:
    EventFeatureCollector(const RegionConstraints&,const RegionFormation&,const ProgramAnalysis& current,
        const std::vector<measurement::Role>&,double lraCeiling,juce::String binding,
        juce::String sourceSnapshot,juce::String baselinePlan,const measurement::Result* measuredSource=nullptr);
    ~EventFeatureCollector();
    EventFeatureCollector(const EventFeatureCollector&)=delete;
    EventFeatureCollector& operator=(const EventFeatureCollector&)=delete;
    // Source is verified pre-processing programme PCM; Current is its exact frozen plan.
    void push(const float* const* source,const float* const* current,int count,
        const std::function<bool()>& cancelled=[] {return false;});
    std::shared_ptr<const EventEvidence> finish(const PredictionSidecar* baseline=nullptr,const std::function<bool()>& cancelled=[] {return false;});
    size_t reservedBytes()const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};
enum class EventCoverageClass { notApplicable, supportedComplete, knownPartial, uncertain, partialAndUncertain };
struct GroupTreatment
{
    juce::String groupId;
    RegionSpan extent;
    bool needBearing=false,knownPartial=false,incomplete=false;
    juce::String reason;
};
struct EventCoverage
{
    EventCoverageClass category=EventCoverageClass::notApplicable;
    std::vector<GroupTreatment> groups;
    std::vector<int64_t> signature;
    size_t unsupportedTransitions=0,newInternalTransitions=0;
    juce::String evidenceIdentity,identity,notice;
    EventCoverage forRanking()const;
    bool complete()const {return category==EventCoverageClass::supportedComplete;}
    juce::String text()const;
};
const char* eventCoverageName(EventCoverageClass);
EventCoverage assessEventCoverage(const processing::SectionGainPlan& before,const RegionPlan&,
    const RegionConstraints&,const EventEvidence&,const ProgramAnalysis&,double ceiling);
// Additional nominees share the original formation caps. v1 path is untouched.
RegionFormation formEventRegions(const RegionConstraints&,const ProgramAnalysis&,double ceiling,
    const EventEvidence&,const std::function<bool()>& cancelled=[] {return false;});
struct EventAssessmentChoice
{
    RegionReviewScore score;EventCoverage coverage;
    bool measured=false,original=false,screening=false,localised=false;
};
std::optional<size_t> chooseEventAssessment(const std::vector<EventAssessmentChoice>&,bool first,bool completePass);
std::optional<size_t> selectCompleteEventPass(const std::vector<RegionReviewScore>&,
    const std::vector<EventCoverage>&);
}

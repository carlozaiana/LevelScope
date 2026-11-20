#pragma once
#include "EventPlanning.h"

namespace levelscope::proposal
{
// One shared lookahead ring. At the minimum hop 21 cells give the entire
// inclusive post-attack reference, including fast/low evidence on every role.
class InnerEventTracker
{
public:
    using Bands=std::array<double,measurement::channelCapacity>;
    InnerEventTracker(const RegionConstraints&,const ProgramAnalysis&,const measurement::Result*,double,EventEvidence&);
    void begin(const EventCell&,int64_t peakFrame,int peakChannel,const Bands& fast,const Bands& low,const Bands& background,const Bands& lowBackground);
    void observe(size_t root,int outer,bool attack,bool activeRise,bool corroborated,int64_t earliest);
    void advance();
    void finish(const PredictionSidecar*,const std::function<bool()>&);
private:
    struct Decision {int outer=-1;bool inContext=false,attack=false,rise=false,confirmed=false;int64_t earliest=0;};
    struct Sample {EventCell cell;int64_t peakFrame=-1;int peakChannel=-1;Bands fast{},low{},background{},lowBackground{};std::array<Decision,8> decisions;};
    struct Reference {int64_t frame=0;double fast=0,low=0;Bands channelFast{},channelLow{};};
    struct Root
    {
        std::array<Reference,2> preceding{};
        int inner=-1,band=0,role=-1;
        double opening=0;
        std::array<double,3> recent{};
        size_t members=0;
        Bands background{},lowBackground{},lastFast{};
        int64_t referenceEnd=0,quiet=-1,fall=-1;
    };
    const RegionConstraints& constraints;
    const ProgramAnalysis& current;
    const measurement::Result* source;
    const double ceiling;
    EventEvidence& evidence;
    std::array<Root,8> roots{};
    std::array<Sample,22> ring{};
    size_t written=0,next=0;
    int rate=0,hop=0,referenceCells=0;
    void process(size_t index);
    void close(size_t root,int64_t end,bool sealed);
    void finaliseCandidates(const PredictionSidecar*,const std::function<bool()>&);
    Reference reference(size_t index)const;
    void member(Root&,EventGroup&,const Sample&,const Reference&,double);
};
}

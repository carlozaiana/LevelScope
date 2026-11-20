#include "InnerEventTracker.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace levelscope::proposal
{
namespace
{
constexpr double floorEnergy=1.e-8;
double level(double energy){return 10*std::log10(std::max(energy,floorEnergy));}

template<size_t N>
double upperMedian(std::array<double,N> values,size_t count)
{
    auto middle=values.begin()+std::ptrdiff_t(count/2);
    std::nth_element(values.begin(),middle,values.begin()+std::ptrdiff_t(count));
    return *middle;
}

template<class W,class End,class Level>
std::optional<double> contextualMedian(const std::vector<W>& windows,int rate,int64_t begin,int64_t end,End frame,Level value)
{
    if(begin<0)return {};
    const int64_t step=rate/10, width=rate*4/10;
    int64_t expected=((begin+width+step-1)/step)*step;
    auto it=std::lower_bound(windows.begin(),windows.end(),expected,[&](const W& w,int64_t f){return frame(w)<f;});
    std::array<double,8> levels{};
    size_t count=0;
    for(;expected<=end;expected+=step,++it)
    {
        if(it==windows.end()||frame(*it)!=expected||count>=levels.size())return {};
        auto v=value(*it);
        if(!v||!std::isfinite(*v))return {};
        levels[count++]=*v;
    }
    if(count==0)return {};
    std::sort(levels.begin(),levels.begin()+std::ptrdiff_t(count));
    return (levels[(count-1)/2]+levels[count/2])*.5;
}
const char* reason(AttackWitness::State state)
{
    switch(state)
    {
    case AttackWitness::ordinary:return "no corroborated inner entrance";
    case AttackWitness::noPredecessors:return "two preceding attacks within 1 s unavailable";
    case AttackWitness::noContext:return "missing or insufficient sustained level change";
    case AttackWitness::noStep:return "reference step below threshold";
    case AttackWitness::localised:return "inner entrance";
    case AttackWitness::continuation:return "comparable member";
    case AttackWitness::conflict:return "distinct entrance while inner slot occupied";
    case AttackWitness::roleOnly:return "role-specific tentative entrance";
    case AttackWitness::clipped:return "post-attack reference clipped";
    }
    return "";
}
}
bool EventEvidence::reserveSummary(size_t root)
{
    if(summaryCounts[root]>=EventFeatureDefinition::maximumGroupsPerRoot)
    {incompleteRoots[root]=true;return false;}
    ++summaryCounts[root];return true;
}
namespace
{
juce::String witnessLine(const AttackWitness& w)
{
        juce::String text="Source witness root "+juce::String(w.root+1)+" frame "+juce::String(w.decision)
            +" onset support "+juce::String(w.earliest)+" confirmation "+juce::String(w.confirmed)
            +" prior decisions "+juce::String(w.predecessor1)+"/"+juce::String(w.predecessor2)
            +" fast/low reference dB "+juce::String(w.fastDb,3)+"/"+juce::String(w.lowDb,3)
            +" entry step "+juce::String(w.stepDb,3)+" dB";
        if(w.beforeMomentary&&w.afterMomentary)text+=" surrounding M "+juce::String(*w.beforeMomentary,2)+"/"+juce::String(*w.afterMomentary,2)+" LUFS";
        text+="; "+juce::String(reason(w.state))+"; inner decision "+juce::String(w.innerEntry)+"\n";
    return text;
}
juce::String groupLine(const EventGroup& g,const std::vector<EventGroup>& groups)
{
    return "Inner interpretation "+g.id+" of "+groups[size_t(g.outerGroup)].id
            +"; last comparable ["+juce::String(g.lastAttack.begin)+", "+juce::String(g.lastAttack.end)+")"
            +"; in-decay ["+juce::String(g.inDecay.begin)+", "+juce::String(g.inDecay.end)+")"
            +"; after-return ["+juce::String(g.decay.begin)+", "+juce::String(g.decay.end)+")"
            +"; captured window excess "+juce::String(g.capturedExcess,5)
            +(g.innerNominee?" (retained pair)":" (unselected or incomplete pair)")+"\n";
}
}
juce::String EventEvidence::innerDigest()const
{
    juce::String state;
    const auto absorb=[&](const juce::String& row){const auto bytes=state+row;state=juce::SHA256(bytes.toRawUTF8(),bytes.getNumBytesAsUTF8()).toHexString();};
    for(const auto& w:attacks)absorb(witnessLine(w));
    for(const auto& g:groups)if(g.outerGroup>=0)absorb(groupLine(g,groups));
    absorb(juce::String(omittedWitnesses)+"/"+juce::String(conflictingEntrances));return state;
}
juce::String EventEvidence::innerDescription()const
{
    juce::String text;for(const auto& w:attacks)text+=witnessLine(w);
    for(const auto& g:groups)if(g.outerGroup>=0)text+=groupLine(g,groups);
    return text+"Omitted witnesses "+juce::String(omittedWitnesses)+", occupied-slot conflicts "+juce::String(conflictingEntrances)+"\n";
}
InnerEventTracker::InnerEventTracker(const RegionConstraints& c,const ProgramAnalysis& a,const measurement::Result* s,double target,EventEvidence& e)
    :constraints(c),current(a),source(s),ceiling(target),evidence(e),rate(c.sampleRate),hop(c.settings.eventFeatures.hopMs),
     referenceCells(EventFeatureDefinition::lastAttackMs/hop)
{
    evidence.attacks.reserve(c.roots.size()*EventFeatureDefinition::maximumGroupsPerRoot);
}
void InnerEventTracker::begin(const EventCell& cell,int64_t peakFrame,int peakChannel,const Bands& fast,const Bands& low,const Bands& background,const Bands& lowBackground)
{
    auto& s=ring[written%ring.size()];s={};s.cell=cell;s.peakFrame=peakFrame;s.peakChannel=peakChannel;s.fast=fast;s.low=low;s.background=background;s.lowBackground=lowBackground;
}
void InnerEventTracker::observe(size_t r,int outer,bool attack,bool rise,bool confirmed,int64_t earliest)
{
    ring[written%ring.size()].decisions[r]={outer,true,attack,rise,confirmed,earliest};
}
void InnerEventTracker::advance()
{
    ++written;
    if(written-next>size_t(referenceCells))process(next++);
}
InnerEventTracker::Reference InnerEventTracker::reference(size_t index)const
{
    Reference ref;ref.frame=ring[index%ring.size()].cell.end;
    const size_t count=size_t(referenceCells)+1;
    std::array<double,21> fast{},low{};
    for(size_t n=0;n<count;++n){const auto& c=ring[(index+n)%ring.size()].cell;fast[n]=c.broadband;low[n]=c.low;}
    ref.fast=upperMedian(fast,count);ref.low=upperMedian(low,count);
    for(size_t ch=0;ch<evidence.roles.size();++ch)
    {
        for(size_t n=0;n<count;++n){const auto& c=ring[(index+n)%ring.size()];fast[n]=c.fast[ch];low[n]=c.low[ch];}
        ref.channelFast[ch]=upperMedian(fast,count);ref.channelLow[ch]=upperMedian(low,count);
    }
    return ref;
}
void InnerEventTracker::member(Root& state,EventGroup& g,const Sample& s,const Reference& ref,double strength)
{
    state.recent[state.members++%state.recent.size()]=strength;state.referenceEnd=s.cell.end+rate/10;
    state.lastFast=ref.channelFast;state.quiet=state.fall=-1;
    g.lastAttack={std::max<int64_t>(0,s.cell.end-rate/50),s.cell.end};g.decay={};g.inDecay={};g.endingSupported=false;
}
void InnerEventTracker::close(size_t r,int64_t end,bool sealed)
{
    auto& state=roots[r];if(state.inner<0)return;
    auto& g=evidence.groups[size_t(state.inner)];g.observedEnd=end;
    if(!sealed){g.endingSupported=false;g.incomplete=true;g.reason+="; return unsealed at context limit";}
    state.inner=-1;
}
void InnerEventTracker::process(size_t index)
{
    const auto& s=ring[index%ring.size()];const auto t=s.cell.end;
    const bool full=index+size_t(referenceCells)<written;
    std::optional<Reference> ref;
    const auto measuredM=[&](int64_t a,int64_t b)
    {
        return contextualMedian(current.momentary,rate,a,b,[](const auto& w){return w.end;},
            [](const DiagnosticWindow& w){return w.padding?std::optional<double>{}:w.lufs;});
    };
    for(size_t r=0;r<constraints.roots.size();++r)
    {
        auto& state=roots[r];const auto& d=s.decisions[r];
        if(!d.inContext){close(r,s.cell.begin,false);state.preceding={};continue;}
        if(d.attack)
        {
            AttackWitness w;w.root=r;w.decision=t;w.predecessor1=state.preceding[0].frame;w.predecessor2=state.preceding[1].frame;
            for(size_t j=index;j<written;++j)
            {
                const auto& look=ring[j%ring.size()];const auto& decision=look.decisions[r];
                if(!decision.inContext||!decision.rise||look.cell.end-t>rate*80/1000)break;
                if(decision.confirmed){w.confirmed=look.cell.end;w.earliest=decision.earliest;break;}
            }
            if(!full||!ring[(index+size_t(referenceCells))%ring.size()].decisions[r].inContext)
            {w.state=AttackWitness::clipped;evidence.incompleteRoots[r]=true;state.preceding={};}
            else
            {
                if(!ref)ref=reference(index);
                w.fastDb=level(ref->fast);w.lowDb=level(ref->low);
                const auto& p=state.preceding;
                int band=0,role=-1;
                bool distinct=false,sourceStep=false;
                if(w.confirmed&&d.outer>=0)
                {
                    if(!p[0].frame||t-p[0].frame>rate)w.state=AttackWitness::noPredecessors;
                    else
                    {
                        const double fastStep=w.fastDb-level(std::max(p[0].fast,p[1].fast));
                        const double lowStep=w.lowDb-level(std::max(p[0].low,p[1].low));
                        band=fastStep>=EventFeatureDefinition::referenceStepDb?0:1;
                        w.stepDb=band==0?fastStep:lowStep;
                        if(std::max(fastStep,lowStep)<EventFeatureDefinition::referenceStepDb)
                        {
                            w.stepDb=std::max(fastStep,lowStep);
                            for(size_t ch=0;ch<evidence.roles.size();++ch)if(evidence.roles[ch]!=measurement::Role::lfe)
                            {
                                for(int b=0;b<2;++b)
                                {
                                    const auto& x=b?ref->channelLow:ref->channelFast;
                                    const auto& y=b?p[0].channelLow:p[0].channelFast;
                                    const auto& z=b?p[1].channelLow:p[1].channelFast;
                                    const double step=level(x[ch])-level(std::max(y[ch],z[ch]));
                                    if(step>w.stepDb){w.stepDb=step;role=int(ch);band=b;}
                                }
                            }
                        }
                        w.state=AttackWitness::noStep;
                        if(w.stepDb>=EventFeatureDefinition::referenceStepDb)
                        {
                            w.beforeMomentary=measuredM(t-rate,t);w.afterMomentary=measuredM(t,t+rate);
                            distinct=w.beforeMomentary&&w.afterMomentary&&*w.afterMomentary-*w.beforeMomentary>=EventFeatureDefinition::momentaryStepLu;
                            w.state=AttackWitness::noContext;
                            if(distinct&&source&&source->sampleRate==rate&&source->frames==current.frames)
                            {
                                const auto src=[&](int64_t a,int64_t b){return contextualMedian(source->momentary,rate,a,b,
                                    [](const auto& v){return v.endFrame;},[](const auto& v)->std::optional<double>{
                                        if(v.energy<=0)return {};return measurement::OfflineMeasurement::loudness(v.energy);});};
                                const auto left=src(t-rate,t),right=src(t,t+rate);
                                sourceStep=left&&right&&*right-*left>=EventFeatureDefinition::momentaryStepLu;
                            }
                        }
                    }
                }
                const auto strength=[&](int b,int ch){return ch<0?(b?ref->low:ref->fast):(b?ref->channelLow[size_t(ch)]:ref->channelFast[size_t(ch)]);};
                if(state.inner<0&&distinct&&evidence.reserveSummary(r))
                {
                    EventGroup g;g.root=r;g.outerGroup=d.outer;g.onset={w.earliest,w.confirmed};g.entryDecision=t;g.entryStepDb=w.stepDb;
                    g.corroborated=true;g.incomplete=role>=0||!sourceStep;g.witnesses=uint16_t((1u<<s.cell.fastWitness)|(1u<<s.cell.lowWitness));
                    g.onsetPeakFrame=s.peakFrame;g.onsetPeakChannel=s.peakChannel;g.onsetPeak=s.cell.rawPeak;
                    g.reason=role>=0?"Tentative role-specific inner rise":"Corroborated inner Source rise within continuing activity";
                    if(!sourceStep)g.reason+="; independent Source momentary step unavailable or conflicting";
                    auto key=(evidence.binding+";inner;"+juce::String(r)+";"+juce::String(t)).toStdString();
                    g.id=juce::SHA256(key.data(),key.size()).toHexString();
                    evidence.groups.push_back(g);state.inner=int(evidence.groups.size()-1);state.role=role;state.band=band;
                    state.opening=strength(band,role);state.members=0;state.background=s.background;state.lowBackground=s.lowBackground;
                    member(state,evidence.groups.back(),s,*ref,state.opening);w.innerEntry=t;
                    w.state=role>=0?AttackWitness::roleOnly:AttackWitness::localised;
                }
                else if(state.inner>=0)
                {
                    auto& g=evidence.groups[size_t(state.inner)];
                    const auto cutoff=std::max(state.opening,upperMedian(state.recent,std::min(state.members,state.recent.size())));
                    const auto value=strength(state.band,state.role);
                    if(distinct)
                    {
                        ++evidence.conflictingEntrances;g.incomplete=true;w.state=AttackWitness::conflict;
                        if(!g.reason.contains("occupied"))g.reason+="; another distinct rise while the inner slot was occupied";
                    }
                    if(level(value)>=level(cutoff)-EventFeatureDefinition::comparableDb)
                    {member(state,g,s,*ref,value);if(!distinct)w.state=AttackWitness::continuation;w.innerEntry=g.entryDecision;}
                }
                state.preceding[0]=state.preceding[1];state.preceding[1]=*ref;
            }
            if(evidence.reserveSummary(r))evidence.attacks.push_back(w);else ++evidence.omittedWitnesses;
        }
        if(state.inner<0||t<state.referenceEnd)continue;
        auto& g=evidence.groups[size_t(state.inner)];
        bool quiet=true,fall=true;
        const auto near=[](double x,double bg){return bg<=floorEnergy?x<=floorEnergy:x<=bg*std::pow(10.,.3);};
        for(size_t ch=0;ch<evidence.roles.size();++ch)
        {
            quiet&=near(s.fast[ch],state.background[ch])&&near(s.low[ch],state.lowBackground[ch]);
            if(state.lastFast[ch]>floorEnergy)fall&=s.fast[ch]<=state.lastFast[ch]*std::pow(10.,-.6);
        }
        if(quiet){if(state.quiet<0)state.quiet=s.cell.begin;}
        else{state.quiet=-1;g.decay={};g.endingSupported=false;}
        if(fall){if(state.fall<0)state.fall=s.cell.begin;}else state.fall=-1;
        if(state.fall>=0&&!g.inDecay.end&&t-state.fall>=rate/5)g.inDecay={state.fall,t};
        if(state.quiet>=0&&t-state.quiet>=rate*3/10)
        {
            if(!g.decay.end)g.decay={state.quiet,t};
            if(t-g.lastAttack.end>=rate){g.endingSupported=true;close(r,t,true);}
        }
    }
}
void InnerEventTracker::finaliseCandidates(const PredictionSidecar* bins,const std::function<bool()>& cancelled)
{
    std::array<EventGroup*,8> best{};
    const bool valid=bins&&bins->binding==evidence.binding&&bins->baselineSnapshot==evidence.currentSnapshot
        &&bins->baselinePlan==evidence.baselinePlan&&bins->sampleRate==rate&&bins->frames==evidence.frames&&bins->roles==evidence.roles;
    const auto detail=[&](RegionSpan span){return span.end>span.begin&&std::any_of(evidence.retainedContexts.begin(),evidence.retainedContexts.end(),
        [&](RegionSpan s){return s.begin<=span.begin&&span.end<=s.end;});};
    for(auto& g:evidence.groups)if(g.outerGroup>=0)
    {
        const auto& hint=constraints.roots[g.root].hint;
        const int64_t end=g.endingSupported?g.decay.end:g.observedEnd;
        const auto start=std::max(hint.begin,g.onset.begin),stop=std::min(hint.end,end);
        if(!detail(g.onset)||!detail(g.decay)){g.incomplete=true;g.reason+="; onset or ending detail incomplete";}
        auto w=std::lower_bound(current.lraWindows.begin(),current.lraWindows.end(),start,
            [](const auto& x,int64_t f){return x.begin+(x.end-x.begin)/2<f;});
        for(;w!=current.lraWindows.end()&&w->begin+(w->end-w->begin)/2<stop;++w)
        {
            if(cancelled())throw std::runtime_error("Inner association cancelled");
            if(!w->included||!w->lufs||*w->lufs<=*current.p10+ceiling)continue;
            const auto idx=int64_t(w-current.lraWindows.begin());if(!g.needWindowCount)g.firstNeedWindow=idx;
            g.lastNeedWindow=idx;++g.needWindowCount;g.needAssociated=true;
            if(!valid||w->energy<=0)continue;
            const auto first=std::max(g.onset.begin,w->begin),last=std::min(end,w->end);
            double captured=0;
            auto b=std::lower_bound(bins->bins.begin(),bins->bins.end(),first,[](const auto& bin,int64_t f){return bin.end<=f;});
            for(;b!=bins->bins.end()&&b->begin<last&&!b->drain;++b)
            {
                const auto overlap=std::min(last,b->end)-std::max(first,b->begin);
                captured+=b->weightedEnergy*double(overlap)/double(b->end-b->begin);
            }
            g.capturedExcess+=(*w->lufs-*current.p10-ceiling)*std::clamp(captured/(w->energy*double(w->end-w->begin)),0.,1.);
        }
        if(!valid){g.incomplete=true;g.reason+="; matching Current attribution unavailable";}
        if(!g.needAssociated||g.capturedExcess<=0||(!g.endingSupported&&!g.inDecay.end))continue;
        auto& incumbent=best[g.root];
        if(!incumbent||std::make_tuple(g.capturedExcess,g.entryStepDb,-g.onset.begin)>
                       std::make_tuple(incumbent->capturedExcess,incumbent->entryStepDb,-incumbent->onset.begin))incumbent=&g;
    }
    for(auto* g:best)if(g)g->innerNominee=true;
}
void InnerEventTracker::finish(const PredictionSidecar* bins,const std::function<bool()>& cancelled)
{
    while(next<written)process(next++);
    for(size_t r=0;r<constraints.roots.size();++r)close(r,evidence.frames,false);
    finaliseCandidates(bins,cancelled);
}
}

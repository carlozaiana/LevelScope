#include "EventPlanning.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace levelscope::proposal
{
namespace
{
void require(bool v,const char* message){if(!v)throw std::invalid_argument(message);}
juce::String digest(const juce::String& s){auto b=s.toStdString();return juce::SHA256(b.data(),b.size()).toHexString();}
bool overlap(RegionSpan a,RegionSpan b){return a.begin<b.end&&b.begin<a.end;}
constexpr size_t capacity=measurement::channelCapacity;
constexpr double floorEnergy=1.e-8;
double ratio(double a,double b){return a>floorEnergy&&b>floorEnergy?10*std::log10(a/b):0;}
double median(std::array<double,64> v,size_t count)
{if(!count)return 0;std::sort(v.begin(),v.begin()+std::ptrdiff_t(count));return v[count/2];}
}
size_t EventEvidence::storageBytes()const
{
    size_t n=sizeof(*this)+roles.capacity()*sizeof(measurement::Role)+cells.capacity()*sizeof(EventCell)
        +retainedContexts.capacity()*sizeof(RegionSpan)+groups.capacity()*sizeof(EventGroup)+(incompleteRoots.capacity()+7)/8;
    for(const auto& g:groups)n+=(g.id.getNumBytesAsUTF8()+g.reason.getNumBytesAsUTF8())*4;
    n+=needCuts.capacity()*sizeof(std::vector<RegionCut>);for(const auto& row:needCuts){n+=row.capacity()*sizeof(RegionCut);for(const auto& cut:row)n+=cut.reason.getNumBytesAsUTF8()*4;}
    for(const auto* s:{&binding,&sourceSnapshot,&currentSnapshot,&baselinePlan,&constraintsIdentity,&definitionDigest,&identity})n+=s->getNumBytesAsUTF8()*4;
    return n;
}
struct EventFeatureCollector::Impl
{
    struct Filter
    {
        double x=0,hp=0,lp=0,cx=0,chp=0;
        void process(double source,double current,double ah,double al)
        {hp=ah*(hp+source-x);x=source;lp=(1-al)*hp+al*lp;chp=ah*(chp+current-cx);cx=current;}
    };
    struct Hop {int64_t begin=0,end=0;std::array<double,capacity> fast{},low{};double current=0,lfe=0,peak=0;int64_t peakFrame=-1;int peakChannel=-1;};
    struct Feature {int64_t begin=0,end=0;std::array<double,capacity> fast{},low{};double current=0;};
    struct RootState
    {
        std::optional<size_t> active;
        int64_t fastStart=-1,lowStart=-1,quietStart=-1,decayStart=-1,referenceEnd=0;
        bool latched=false;
        size_t referenceCount=0,groups=0;
        std::array<std::array<double,64>,capacity> reference{};
        std::array<double,capacity> background{},lowBackground{},lastReference{};
    };
    const RegionConstraints& constraints;
    const ProgramAnalysis& currentAnalysis;
    double ceiling=0,ah=0,al=0;
    std::shared_ptr<EventEvidence> data=std::make_shared<EventEvidence>();
    std::array<Filter,capacity> filters;
    std::array<Hop,16> hops;
    std::array<Feature,64> backgrounds;
    std::array<RootState,8> roots;
    size_t hopCount=0,backgroundCount=0,nonLfe=0;
    Hop accumulating;
    int64_t position=0,cellIndex=0,nextEnd=0;
    bool finished=false;
    Impl(const RegionConstraints& c,const RegionFormation& formation,const ProgramAnalysis& a,
        const std::vector<measurement::Role>& roles,double target,juce::String binding,
        juce::String sourceSnapshot,juce::String before):constraints(c),currentAnalysis(a),ceiling(target)
    {
        require(c.error.isEmpty()&&c.settings.schema==2&&c.roots.size()<=8,"Invalid event feature constraints");
        require(a.sampleRate==c.sampleRate&&a.frames==c.frames&&roles.size()<=capacity&&roles.size()>0,"Event feature programme mismatch");
        require(std::isfinite(target)&&a.p10.has_value(),"Event feature need unavailable");
        data->binding=std::move(binding);data->sourceSnapshot=std::move(sourceSnapshot);data->currentSnapshot=a.snapshotIdentity;
        data->baselinePlan=std::move(before);data->sampleRate=c.sampleRate;data->frames=c.frames;data->roles=roles;
        data->constraintsIdentity=c.identity;data->definitionDigest=c.settings.eventFeatures.digest();data->incompleteRoots.resize(c.roots.size(),false);
        data->needCuts=formation.nominations;
        for(auto r:roles)if(r!=measurement::Role::lfe)++nonLfe;require(nonLfe>0,"No non-LFE event context");
        const int64_t halo=int64_t(c.sampleRate)*3;
        std::vector<RegionSpan> contexts;contexts.reserve(c.roots.size()*14);
        const auto retain=[&](int64_t centre){contexts.push_back({std::max<int64_t>(0,centre-halo),std::min(c.frames,centre+halo)});};
        for(size_t n=0;n<c.roots.size();++n)
        {
            retain(c.roots[n].hint.begin);retain(c.roots[n].hint.end);
            if(n<formation.nominations.size())for(size_t k=0;k<std::min<size_t>(12,formation.nominations[n].size());++k)retain(formation.nominations[n][k].frame);
        }
        std::sort(contexts.begin(),contexts.end(),[](auto x,auto y){return x.begin<y.begin;});
        data->retainedContexts.reserve(contexts.size());
        for(auto span:contexts)if(span.begin<span.end)
        {if(!data->retainedContexts.empty()&&span.begin<=data->retainedContexts.back().end)data->retainedContexts.back().end=std::max(span.end,data->retainedContexts.back().end);else data->retainedContexts.push_back(span);}
        size_t cells=0;const auto hop=c.settings.eventFeatures.hopMs;
        for(auto span:data->retainedContexts)cells+=size_t((span.end-span.begin)*1000/(int64_t(c.sampleRate)*hop)+3);
        const size_t groups=c.roots.size()*EventFeatureDefinition::maximumGroupsPerRoot;
        // Account worst-case bounded reason/ID storage before any large reservation.
        const size_t planned=sizeof(Impl)+sizeof(EventEvidence)+cells*sizeof(EventCell)+groups*(sizeof(EventGroup)+2048)+16384;
        require(planned<=EventFeatureDefinition::maximumBytes,"Event feature reservations exceed 10 MiB sub-budget");
        data->cells.reserve(cells);data->groups.reserve(groups);data->reservedBytes=planned;data->highWaterBytes=planned;
        constexpr double pi=3.141592653589793238462643383279502884;
        ah=std::exp(-2*pi*EventFeatureDefinition::hpHz/c.sampleRate);al=std::exp(-2*pi*EventFeatureDefinition::lowHz/c.sampleRate);
        nextEnd=std::min(c.frames,int64_t(c.sampleRate)*hop/1000);
    }
    bool hasDetail(RegionSpan span)const
    {return std::any_of(data->retainedContexts.begin(),data->retainedContexts.end(),[&](auto r){return r.begin<=span.begin&&r.end>=span.end;});}
    void finalizeGroup(size_t r)
    {
        auto& state=roots[r];if(!state.active)return;auto& g=data->groups[*state.active];
        const auto end=g.endingSupported?g.decay.end:std::min(data->frames,position);
        auto first=std::lower_bound(currentAnalysis.lraWindows.begin(),currentAnalysis.lraWindows.end(),std::max(g.onset.begin,constraints.roots[r].hint.begin),[](const auto& w,int64_t f){return w.begin+(w.end-w.begin)/2<f;});
        for(auto i=first;i!=currentAnalysis.lraWindows.end();++i)
        {
            const auto& w=*i;
            const auto mid=w.begin+(w.end-w.begin)/2;if(mid>=end)break;
            if(mid>=std::max(g.onset.begin,constraints.roots[r].hint.begin)&&mid<constraints.roots[r].hint.end&&w.included&&w.lufs&&*w.lufs>*currentAnalysis.p10+ceiling){g.needAssociated=true;const auto index=int64_t(i-currentAnalysis.lraWindows.begin());if(g.firstNeedWindow<0)g.firstNeedWindow=index;g.lastNeedWindow=index;++g.needWindowCount;}
        }
        if(!hasDetail(g.onset)||!g.endingSupported||!hasDetail(g.decay))
        {g.incomplete=true;if(g.reason.isNotEmpty())g.reason+="; ";g.reason+="Onset/ending detail or confirmed ending unavailable";}
        g.id=digest(data->binding+";event;"+juce::String(r)+";"+juce::String(g.onset.begin)+";"+juce::String(g.lastAttack.end));state.active.reset();
    }
    void emit()
    {
        const int hopMs=constraints.settings.eventFeatures.hopMs;
        hops[hopCount%hops.size()]=accumulating;++hopCount;
        Feature feature;feature.end=position;feature.begin=position;
        std::array<double,capacity> bg{},bgLow{};double current=0;int64_t fastCount=0,lowCount=0;
        for(size_t k=0;k<std::min(hopCount,size_t(EventFeatureDefinition::lowMs/hopMs));++k)
        {
            const auto& h=hops[(hopCount-1-k)%hops.size()];const auto count=h.end-h.begin;
            for(size_t ch=0;ch<data->roles.size();++ch)feature.low[ch]+=h.low[ch];lowCount+=count;
            if(k<size_t(EventFeatureDefinition::fastMs/hopMs))
            {feature.begin=h.begin;fastCount+=count;current+=h.current;for(size_t ch=0;ch<data->roles.size();++ch)feature.fast[ch]+=h.fast[ch];}
        }
        const size_t count=std::min(backgroundCount,size_t(EventFeatureDefinition::backgroundMs/hopMs));
        for(size_t ch=0;ch<data->roles.size();++ch)
        {
            feature.fast[ch]/=double(fastCount);feature.low[ch]/=double(lowCount);
            std::array<double,64> f{},l{};for(size_t k=0;k<count;++k){const auto& x=backgrounds[(backgroundCount-1-k)%backgrounds.size()];f[k]=x.fast[ch];l[k]=x.low[ch];}
            bg[ch]=median(f,count);bgLow[ch]=median(l,count);
        }
        EventCell cell;cell.begin=accumulating.begin;cell.end=position;cell.current=current/double(fastCount*int64_t(nonLfe));cell.rawPeak=float(accumulating.peak);
        bool fast=false,low=false;double bestFast=-1e10,bestLow=-1e10;
        for(size_t ch=0;ch<data->roles.size();++ch)
        {
            if(data->roles[ch]==measurement::Role::lfe){cell.lfe+=feature.low[ch];continue;}
            cell.broadband+=feature.fast[ch]/double(nonLfe);cell.low+=feature.low[ch]/double(nonLfe);
            const double f=ratio(feature.fast[ch],bg[ch]),l=ratio(feature.low[ch],bgLow[ch]);
            const bool fromSilence=bg[ch]<=floorEnergy&&feature.fast[ch]>floorEnergy;
            fast|=feature.fast[ch]>floorEnergy&&(fromSilence||f>=EventFeatureDefinition::fastRiseDb);
            low|=feature.low[ch]>floorEnergy&&(bgLow[ch]<=floorEnergy||l>=EventFeatureDefinition::lowRiseDb);
            if(feature.fast[ch]>floorEnergy&&f>bestFast){bestFast=f;cell.fastWitness=uint8_t(ch);}if(feature.low[ch]>floorEnergy&&l>bestLow){bestLow=l;cell.lowWitness=uint8_t(ch);}
            if(fromSilence)cell.flags|=1;
        }
        feature.current=cell.current;std::array<double,64> currentBackground{};
        for(size_t k=0;k<count;++k)currentBackground[k]=backgrounds[(backgroundCount-1-k)%backgrounds.size()].current;
        const double cb=median(currentBackground,count);
        if(!fast&&!low&&cell.current>floorEnergy&&(cb<=floorEnergy||ratio(cell.current,cb)>=EventFeatureDefinition::fastRiseDb))cell.flags|=2; // Current-only rise: processing-associated context, never Source event authority
        cell.fastContrast=float(bestFast==-1e10?0:bestFast);cell.lowContrast=float(bestLow==-1e10?0:bestLow);
        require(std::isfinite(cell.broadband)&&std::isfinite(cell.low)&&std::isfinite(cell.current),"Non-finite event feature accumulation");
        const bool detail=std::any_of(data->retainedContexts.begin(),data->retainedContexts.end(),[&](auto s){return overlap(s,{cell.begin,cell.end});});
        if(detail){require(data->cells.size()<data->cells.capacity(),"Event detail reservation exhausted");data->cells.push_back(cell);}
        for(size_t r=0;r<constraints.roots.size();++r)
        {
            const auto& root=constraints.roots[r];const auto halo=int64_t(data->sampleRate)*3;auto& s=roots[r];
            if(cell.end<=std::max<int64_t>(0,root.hint.begin-halo))continue;
            if(cell.begin>=std::min(data->frames,root.hint.end+halo)){finalizeGroup(r);continue;}
            if(fast){if(s.fastStart<0)s.fastStart=cell.begin;}else s.fastStart=-1;
            if(low){if(s.lowStart<0)s.lowStart=cell.begin;}else s.lowStart=-1;
            const bool f=fast&&s.fastStart>=0&&position-s.fastStart>=int64_t(data->sampleRate)*20/1000;
            const bool l=low&&s.lowStart>=0&&position-s.lowStart>=int64_t(data->sampleRate)*40/1000;
            const bool dual=f&&l&&std::abs(s.fastStart-s.lowStart)<=int64_t(data->sampleRate)*80/1000;
            const bool sustained=f&&position-s.fastStart>=int64_t(data->sampleRate)*80/1000;
            const bool corroborated=(dual||sustained)&&count==size_t(EventFeatureDefinition::backgroundMs/hopMs);
            const bool attack=(f||l)&&!s.latched;
            if(!fast&&!low)s.latched=false;
            if(attack)
            {
                s.latched=true;
                if(!s.active)
                {
                    if(s.groups>=EventFeatureDefinition::maximumGroupsPerRoot){data->incompleteRoots[r]=true;++data->omittedGroups;continue;}
                    ++s.groups;EventGroup g;g.root=r;
                    const auto firstFast=s.fastStart<0?position:std::max<int64_t>(0,s.fastStart-int64_t(data->sampleRate)*(20-hopMs)/1000);
                    const auto firstLow=s.lowStart<0?position:std::max<int64_t>(0,s.lowStart-int64_t(data->sampleRate)*(80-hopMs)/1000);
                    g.onset={std::min(firstFast,firstLow),position};g.onsetPeak=accumulating.peak;g.onsetPeakFrame=accumulating.peakFrame;g.onsetPeakChannel=accumulating.peakChannel;
                    g.background=std::accumulate(bg.begin(),bg.begin()+std::ptrdiff_t(data->roles.size()),0.)/double(data->roles.size());
                    g.lowBackground=std::accumulate(bgLow.begin(),bgLow.begin()+std::ptrdiff_t(data->roles.size()),0.)/double(data->roles.size());
                    g.corroborated=corroborated;g.incomplete=!corroborated;g.reason=corroborated?"Broadband and low-band rise with persistent Source evidence":"Tentative Source rise; corroboration incomplete";
                    data->groups.push_back(std::move(g));s.active=data->groups.size()-1;s.background=bg;s.lowBackground=bgLow;
                }
                auto& g=data->groups[*s.active];g.lastAttack={feature.begin,position};g.endingSupported=false;g.decay={};g.inDecay={};
                g.witnesses|=uint16_t((1u<<cell.fastWitness)|(1u<<cell.lowWitness));
                s.referenceEnd=position+int64_t(data->sampleRate)*100/1000;s.referenceCount=0;s.quietStart=s.decayStart=-1;
            }
            if(s.active)
            {
                auto& g=data->groups[*s.active];
                // Later corroboration of the SAME entrance may complete its evidence.
                if(corroborated&&position-g.onset.end<=int64_t(data->sampleRate)*80/1000&&!g.corroborated)
                {g.corroborated=true;g.incomplete=false;g.reason=dual?"Broadband and low-band persistent Source rise":"Source rise corroborated by sustained fast-energy elevation";}
                if(position<=s.referenceEnd&&s.referenceCount<64)
                {for(size_t ch=0;ch<data->roles.size();++ch)s.reference[ch][s.referenceCount]=feature.fast[ch];++s.referenceCount;}
                if(position>=s.referenceEnd&&s.referenceCount)
                {
                    bool quiet=true,decaying=true;
                    for(size_t ch=0;ch<data->roles.size();++ch)
                    {
                        s.lastReference[ch]=median(s.reference[ch],s.referenceCount);
                        const auto near=[](double x,double bg){return bg<=floorEnergy?x<=floorEnergy:x<=bg*std::pow(10.,.3);};
                        quiet&=near(feature.fast[ch],s.background[ch])&&near(feature.low[ch],s.lowBackground[ch]);
                        if(s.lastReference[ch]>floorEnergy)decaying&=feature.fast[ch]<=s.lastReference[ch]*std::pow(10.,-.6);
                    }
                    if(quiet){if(s.quietStart<0)s.quietStart=cell.begin;}else {s.quietStart=-1;g.endingSupported=false;g.decay={};}
                    if(decaying){if(s.decayStart<0)s.decayStart=cell.begin;}else s.decayStart=-1;
                    if(s.decayStart>=0&&position-s.decayStart>=int64_t(data->sampleRate)*200/1000&&g.inDecay.end==0)g.inDecay={s.decayStart,position};
                    if(s.quietStart>=0&&position-s.quietStart>=int64_t(data->sampleRate)*300/1000)
                    {if(!g.endingSupported){g.decay={s.quietStart,position};g.endingSupported=true;}if(position-g.lastAttack.end>=data->sampleRate)finalizeGroup(r);}
                }
            }
        }
        backgrounds[backgroundCount%backgrounds.size()]=feature;++backgroundCount;
        accumulating={};accumulating.begin=position;++cellIndex;
        nextEnd=std::min(data->frames,(cellIndex+1)*int64_t(data->sampleRate)*hopMs/1000);
    }
};
EventFeatureCollector::EventFeatureCollector(const RegionConstraints& c,const RegionFormation& f,const ProgramAnalysis& a,
    const std::vector<measurement::Role>& roles,double ceiling,juce::String binding,juce::String source,juce::String before)
:impl(std::make_unique<Impl>(c,f,a,roles,ceiling,std::move(binding),std::move(source),std::move(before))){}
EventFeatureCollector::~EventFeatureCollector()=default;
void EventFeatureCollector::push(const float* const* source,const float* const* current,int count,const std::function<bool()>& cancelled)
{
    auto& x=*impl;require(!x.finished&&source&&current&&count>=0&&count<=x.data->frames-x.position,"Invalid event feature chunk");
    for(size_t ch=0;ch<x.data->roles.size();++ch)require(source[ch]&&current[ch],"Missing event feature channel");
    for(int i=0;i<count;++i)
    {
        if((i%256)==0&&cancelled())throw std::runtime_error("Event collection cancelled");
        for(size_t ch=0;ch<x.data->roles.size();++ch)
        {
            const double s=source[ch][i],c=current[ch][i];require(std::isfinite(s)&&std::isfinite(c),"Non-finite event PCM");
            auto& f=x.filters[ch];f.process(s,c,x.ah,x.al);x.accumulating.fast[ch]+=f.hp*f.hp;x.accumulating.low[ch]+=f.lp*f.lp;
            if(std::abs(s)>x.accumulating.peak){x.accumulating.peak=std::abs(s);x.accumulating.peakFrame=x.position;x.accumulating.peakChannel=int(ch);}
            if(x.data->roles[ch]!=measurement::Role::lfe)x.accumulating.current+=f.chp*f.chp;
        }
        ++x.position;x.accumulating.end=x.position;if(x.position==x.nextEnd)x.emit();
    }
}
std::shared_ptr<const EventEvidence> EventFeatureCollector::finish()
{
    auto& x=*impl;if(x.finished)return x.data;require(x.position==x.data->frames,"Incomplete Source event pass");
    for(size_t r=0;r<x.constraints.roots.size();++r)x.finalizeGroup(r);
    juce::String summaries;for(const auto& g:x.data->groups)summaries+=g.id+":"+juce::String(g.onset.begin)+":"+juce::String(g.decay.end)+":"+juce::String(int(g.corroborated))+":"+juce::String(int(g.incomplete))+":"+juce::String(int(g.needAssociated))+";";
    const auto cellHash=juce::SHA256(x.data->cells.data(),x.data->cells.size()*sizeof(EventCell)).toHexString();
    x.data->identity=digest(x.data->binding+";"+x.data->sourceSnapshot+";"+x.data->currentSnapshot+";"+x.data->baselinePlan+";"+x.data->constraintsIdentity+";"+x.data->definitionDigest+";"+cellHash+";"+summaries);
    x.data->highWaterBytes=std::max(x.data->reservedBytes,x.data->storageBytes()+sizeof(Impl));
    require(x.data->highWaterBytes<=EventFeatureDefinition::maximumBytes,"Event high-water exceeds 10 MiB");
    x.finished=x.data->finished=true;return x.data;
}
size_t EventFeatureCollector::reservedBytes()const{return impl->data->reservedBytes;}
}

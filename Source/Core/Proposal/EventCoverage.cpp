#include "EventPlanning.h"
#include <algorithm>
#include <set>
#include <tuple>

namespace levelscope::proposal
{
namespace
{
bool overlap(RegionSpan a,RegionSpan b){return a.begin<b.end&&b.begin<a.end;}
juce::String hash(const juce::String& s){auto b=s.toStdString();return juce::SHA256(b.data(),b.size()).toHexString();}
using Plan=processing::SectionGainPlan;
void points(const Plan& p,std::vector<int64_t>& out)
{for(size_t n=1;n<p.intervals.size();++n){const auto f=p.intervals[n].begin,r=p.radius(n);for(auto x:{f-r,f-r+1,f+r-1,f+r})if(x>=0&&x<p.frames)out.push_back(x);}}
bool needIn(const ProgramAnalysis& a,RegionSpan span,double ceiling)
{
    if(!a.p10)return false;
    auto i=std::lower_bound(a.lraWindows.begin(),a.lraWindows.end(),span.begin,[](const auto& w,int64_t x){return w.begin+(w.end-w.begin)/2<x;});
    for(;i!=a.lraWindows.end();++i){const auto f=i->begin+(i->end-i->begin)/2;if(f>=span.end)break;if(i->included&&i->lufs&&*i->lufs>*a.p10+ceiling)return true;}return false;
}
}
const char* eventCoverageName(EventCoverageClass c)
{
    switch(c){case EventCoverageClass::supportedComplete:return "Supported complete";case EventCoverageClass::knownPartial:return "Known partial";
        case EventCoverageClass::uncertain:return "Uncertain";case EventCoverageClass::partialAndUncertain:return "Known partial + uncertain";default:return "Not an event-treatment candidate";}
}
juce::String EventCoverage::text()const
{
    juce::String s=juce::String("Event coverage: ")+eventCoverageName(category)+". Editorial evidence, not a measured delivery rule.\n"+notice+"\n";
    for(const auto& g:groups)s+=g.groupId+" ["+juce::String(g.extent.begin)+", "+juce::String(g.extent.end)+"): "+g.reason+"\n";
    return s;
}
EventCoverage assessEventCoverage(const Plan& before,const RegionPlan& plan,const RegionConstraints& c,
    const EventEvidence& evidence,const ProgramAnalysis& analysis,double ceiling)
{
    EventCoverage out;out.evidenceIdentity=evidence.identity;
    if(!evidence.finished||evidence.constraintsIdentity!=c.identity||evidence.definitionDigest!=c.settings.eventFeatures.digest()||evidence.currentSnapshot!=analysis.snapshotIdentity||evidence.frames!=before.frames||evidence.sampleRate!=before.sampleRate||plan.error.isNotEmpty())
    {out.category=EventCoverageClass::uncertain;out.notice="Missing/stale evidence or invalid plan";return out;}
    const auto delta=[&](int64_t f){return before.gainAt(f)-plan.envelope.gainAt(f);};
    bool need=false,partial=false,incomplete=false;std::set<size_t> changedRoots;
    for(const auto& a:plan.assignments)
    {
        if(a.geometry.root>=c.roots.size())continue;
        for(auto d:plan.activity.support)if(overlap(d,a.geometry.support)&&needIn(analysis,{std::max(d.begin,c.roots[a.geometry.root].hint.begin),std::min(d.end,c.roots[a.geometry.root].hint.end)},ceiling))need=true;
        if(std::any_of(a.reductions.begin(),a.reductions.end(),[](double g){return g>0;}))changedRoots.insert(a.geometry.root);
        for(size_t n=1;n<a.reductions.size();++n)if(a.reductions[n]!=a.reductions[n-1])
        {
            ++out.newInternalTransitions;const auto cut=a.geometry.cuts[n-1];
            const bool supported=a.geometry.root<evidence.needCuts.size()&&std::any_of(evidence.needCuts[a.geometry.root].begin(),evidence.needCuts[a.geometry.root].end(),[&](auto x){return x.frame==cut;});
            if(!supported)++out.unsupportedTransitions;
        }
    }
    if(!need||changedRoots.empty()){out.notice="No nonzero local correction associated with measured need; no complete-event credit";return out;}
    for(const auto& g:evidence.groups)
    {
        if(g.root>=c.roots.size())continue;
        const RegionSpan extent{g.onset.begin,g.endingSupported?g.decay.end:std::min(before.frames,c.roots[g.root].hint.end+int64_t(c.sampleRate)*3)};
        if(!std::any_of(plan.activity.support.begin(),plan.activity.support.end(),[&](auto d){return overlap(d,extent);}))continue;
        GroupTreatment t;t.groupId=g.id;t.extent=extent;t.incomplete=g.incomplete||!g.corroborated||!g.endingSupported;
        t.needBearing=g.needAssociated&&std::any_of(plan.activity.support.begin(),plan.activity.support.end(),[&](auto d){return needIn(analysis,{std::max({d.begin,extent.begin,c.roots[g.root].hint.begin}),std::min({d.end,extent.end,c.roots[g.root].hint.end})},ceiling);});
        if(t.needBearing){out.signature.push_back(extent.begin);out.signature.push_back(extent.end);}
        const auto assignment=std::find_if(plan.assignments.begin(),plan.assignments.end(),[&](const auto& a){return a.geometry.root==g.root;});
        double first=0,last=0;
        if(assignment==plan.assignments.end()){t.incomplete=true;t.reason+="Incidental event has no independently established intended schedule; ";}
        if(assignment!=plan.assignments.end())
        {
            for(double d:assignment->reductions)if(d>0){if(first==0)first=d;last=d;}
            // Match the processor's representable plateau difference, not a tolerance.
            const auto representable=[&](double requested,int64_t f){const double b=before.gainAt(std::clamp<int64_t>(f,0,before.frames-1));double after=b-requested;if(b-after>requested)after=std::nextafter(after,b);return b-after;};
            first=representable(first,assignment->nominal.begin+std::llround(c.sampleRate*.05));
            last=representable(last,assignment->nominal.end-1-std::llround(c.sampleRate*.05));
        }
        if(g.corroborated&&first>0&&delta(std::clamp<int64_t>(g.onset.begin,0,before.frames-1))<first)
        {t.knownPartial=true;t.reason+="Intended correction not reached by earliest plausible onset; ";}
        if(g.endingSupported&&last>0&&delta(std::clamp<int64_t>(g.decay.end-1,0,before.frames-1))<last)
        {t.knownPartial=true;t.reason+="Return precedes supported decay end; ";}
        std::vector<int64_t> frames{extent.begin,std::min(before.frames-1,extent.end-1)};points(before,frames);points(plan.envelope,frames);
        for(auto f:frames)if(f>=extent.begin&&f<extent.end&&delta(f)<=0)
        {if(g.corroborated&&g.endingSupported){t.knownPartial=true;t.reason+="Additional correction absent inside supported event; ";}break;}
        if(assignment!=plan.assignments.end()&&g.corroborated&&g.endingSupported)
        {
            // Intended additional schedule is frozen by the chosen reductions.
            // Extend its outer gains through the event; internal changes retain
            // the qualified fixed transition convention, never an epsilon floor.
            Plan intended=before;intended.intervals.clear();
            for(size_t n=0;n<assignment->reductions.size();++n)
                intended.intervals.push_back({n?assignment->geometry.cuts[n-1]:0,n<assignment->geometry.cuts.size()?assignment->geometry.cuts[n]:before.frames,-assignment->reductions[n]});
            intended=intended.canonical();points(intended,frames);
            for(auto f:frames)if(f>=extent.begin&&f<extent.end)
            {
                const double requested=-intended.gainAt(f),b=before.gainAt(f);double target=b-requested;if(b-target>requested)target=std::nextafter(target,b);
                if(requested<=0||delta(f)<b-target){t.knownPartial=true;t.reason+="Intended internal correction schedule not maintained; ";break;}
            }
        }
        if(t.incomplete)t.reason+="Event evidence incomplete/uncertain: "+g.reason+"; ";
        if(!t.incomplete&&!t.knownPartial)t.reason="Intended nonzero correction covers supported onset through decay";
        partial|=t.knownPartial;incomplete|=t.incomplete;out.groups.push_back(std::move(t));
    }
    for(auto root:changedRoots)
    {
        const bool missing=root>=evidence.incompleteRoots.size()||evidence.incompleteRoots[root]
            ||std::none_of(evidence.groups.begin(),evidence.groups.end(),[&](const auto& g){return g.root==root&&g.needAssociated;});
        if(missing){incomplete=true;out.notice+="Root "+juce::String(root+1)+": missing/capped event association; ";}
    }
    if(out.signature.empty()){incomplete=true;out.notice+="No supported event signature for changed measured need; ";}
    out.category=partial?(incomplete?EventCoverageClass::partialAndUncertain:EventCoverageClass::knownPartial):(incomplete?EventCoverageClass::uncertain:EventCoverageClass::supportedComplete);
    std::sort(out.signature.begin(),out.signature.end());
    juce::String key=evidence.identity+";"+plan.identity+";"+eventCoverageVersion+";"+juce::String(int(out.category));
    for(const auto& g:out.groups)key+=";"+g.groupId+":"+juce::String(int(g.knownPartial))+":"+juce::String(int(g.incomplete));out.identity=hash(key);return out;
}
RegionFormation formEventRegions(const RegionConstraints& c,const ProgramAnalysis& a,double ceiling,
    const EventEvidence& evidence,const std::function<bool()>& cancelled)
{
    std::vector<std::vector<RegionGeometry>> preferred(c.roots.size());
    const auto ramp=int64_t(std::llround(c.sampleRate*.1));
    for(size_t root=0;root<c.roots.size();++root)
    {
        std::optional<RegionSpan> complete;int64_t decayReturn=0;
        for(const auto& g:evidence.groups)if(g.root==root&&g.needAssociated&&g.corroborated&&g.endingSupported&&!g.incomplete)
        {
            RegionSpan span{std::max<int64_t>(0,g.onset.begin-ramp),std::min(c.frames,g.decay.end+ramp)};
            if(!complete)complete=span;else{complete->begin=std::min(complete->begin,span.begin);complete->end=std::max(complete->end,span.end);}
            if(g.inDecay.end)decayReturn=std::max(decayReturn,g.inDecay.end+ramp);
        }
        if(complete&&c.settings.allowRefinement)
        {
            RegionGeometry g;g.root=root;g.support=*complete;g.reason="Supported complete-event boundary nominee; exact coverage still checked";g.strength=1e6;preferred[root].push_back(g);
            if(decayReturn>complete->begin&&decayReturn<complete->end){g.support.end=decayReturn;g.reason="Explicit in-decay return alternative; not equivalent complete coverage";preferred[root].push_back(g);}
        }
    }
    auto out=formRegions(c,a,ceiling,cancelled,&preferred);
    for(size_t r=0;r<preferred.size();++r)for(const auto& g:preferred[r])
        if(!c.roots[r].start.contains(g.support.begin)||!c.roots[r].end.contains(g.support.end))out.reasons.push_back("Root "+juce::String(r+1)+": complete-event/decay alternative conflicts with explicit permission, ownership or minimum length");
    return out;
}
std::optional<size_t> selectCompleteEventPass(const std::vector<RegionReviewScore>& scores,const std::vector<EventCoverage>& coverage)
{
    std::optional<size_t> best;
    const auto key=[&](size_t n){const auto& a=scores[n].activity;const auto& c=coverage[n];return std::tuple{c.unsupportedTransitions,a.affectedFrames,a.absoluteDbFrames,a.maximumAdditional,c.newInternalTransitions,std::abs(a.globalDelta),c.signature,scores[n].id};};
    for(size_t n=0;n<scores.size()&&n<coverage.size();++n)if(scores[n].complete&&scores[n].allPass&&coverage[n].complete()&&!coverage[n].signature.empty())
        if(!best||key(n)<key(*best))best=n;
    return best;
}
std::optional<size_t> chooseEventAssessment(const std::vector<EventAssessmentChoice>& choices,bool first,bool completePass)
{
    const auto best=[&](int family,bool screeningOnly)->std::optional<size_t>
    {
        std::vector<RegionReviewScore> scores;std::vector<EventCoverage> coverage;
        std::optional<size_t> fallback;
        for(size_t n=0;n<choices.size();++n)
        {
            const auto& c=choices[n];auto s=c.score;
            const bool eligible=!c.measured&&(family==0?c.coverage.complete():family==1?c.original:family==2?!c.original&&!c.coverage.complete():true)&&(!screeningOnly||c.screening);
            if(!eligible)s.complete=false;
            else
            {
                if(!fallback||s.id<choices[*fallback].score.id)fallback=n;
                if(c.screening)s.allPass=true; // forecast ranking ONLY; never measured authority
            }
            scores.push_back(s);coverage.push_back(c.coverage);
        }
        if(family==0)if(auto n=selectCompleteEventPass(scores,coverage))return n;
        if(auto n=selectRegionReview(scores))return n;
        return fallback;
    };
    if(first){if(auto n=best(0,true))return n;if(auto n=best(1,false))return n;return best(3,false);}
    if(!completePass)if(auto n=best(0,false))return n;
    const bool originalMeasured=std::any_of(choices.begin(),choices.end(),[](const auto& c){return c.measured&&c.original;});
    if(!originalMeasured)if(auto n=best(1,false))return n;
    if(auto n=best(2,false))return n;
    return best(3,false);
}

}

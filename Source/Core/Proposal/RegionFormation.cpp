#include "RegionPlanning.h"
#include <algorithm>
#include <set>
#include <numeric>

namespace levelscope::proposal
{
namespace
{
struct Need {int64_t centre;double q;};
double median(std::vector<double> values)
{std::sort(values.begin(),values.end());return values[values.size()/2];}
std::vector<RegionCut> nominate(const RegionConstraints& c,const PlanningRoot& root,
    const ProgramAnalysis& a,double ceiling,const std::function<bool()>& cancel)
{
    std::vector<Need> windows;
    if(!a.p10)return {};
    const int64_t halo=int64_t(c.sampleRate)*3,step=c.sampleRate/10,persistence=int64_t(c.settings.persistenceMs)*c.sampleRate/1000;
    for(const auto& w:a.lraWindows)
    {
        const auto centre=w.begin+(w.end-w.begin)/2;
        if(w.lufs&&centre>=root.allowed.begin-halo&&centre<root.allowed.end+halo)
            windows.push_back({centre,std::max(0.,*w.lufs-*a.p10-ceiling)});
    }
    std::vector<RegionCut> cuts;
    struct InteriorGap {int64_t begin,end;double contrast;};std::vector<InteriorGap> gaps;
    const double threshold=c.settings.needThresholdTenths/10.;
    const auto at=[&](int64_t f){return std::lower_bound(windows.begin(),windows.end(),f,[](auto w,int64_t x){return w.centre<x;});};
    size_t iterations=0;
    for(int64_t f=((root.allowed.begin+step-1)/step)*step;f<root.allowed.end;f+=step)
    {
        if((iterations++%64)==0&&cancel())return {};
        std::vector<double> left,right;
        for(auto it=at(f-halo);it!=windows.end()&&it->centre<f;++it)left.push_back(it->q);
        for(auto it=at(f);it!=windows.end()&&it->centre<f+halo;++it)right.push_back(it->q);
        if(left.empty()||right.empty())continue;
        const double l=median(left),r=median(right),mid=(l+r)/2;
        if(std::abs(l-r)<threshold)continue;
        bool stable=true;size_t nl=0,nr=0;
        for(auto it=at(f-persistence);it!=windows.end()&&it->centre<f;++it){++nl;if(l<r?it->q>mid:it->q<mid)stable=false;}
        for(auto it=at(f);it!=windows.end()&&it->centre<f+persistence;++it){++nr;if(l<r?it->q<mid:it->q>mid)stable=false;}
        // Missing samples are unavailable context, not invented zero need.
        if(stable&&nl>=size_t(persistence/step)&&nr>=size_t(persistence/step))cuts.push_back({f,std::abs(l-r),"Persistent median correction-need change; centre hypothesis, not event onset",l<r?EdgeDirection::rising:EdgeDirection::falling,f});
    }
    for(size_t n=0;n<windows.size();)
    {
        if(cancel())return {};if(windows[n].q>threshold/2){++n;continue;}
        const size_t first=n++;while(n<windows.size()&&windows[n].q<=threshold/2&&windows[n].centre-windows[n-1].centre==step)++n;
        if(windows[n-1].centre-windows[first].centre+step<persistence)continue;
        std::vector<double> previous,following;
        for(auto it=at(windows[first].centre-halo);it!=windows.end()&&it->centre<windows[first].centre;++it)previous.push_back(it->q);
        for(auto it=at(windows[n-1].centre+step);it!=windows.end()&&it->centre<windows[n-1].centre+step+halo;++it)following.push_back(it->q);
        const double before=previous.empty()?-1:median(previous),after=following.empty()?-1:median(following);
        if(before<=threshold/2&&after<=threshold/2)continue;
        if(before>threshold/2&&after>threshold/2&&windows[first].centre>root.allowed.begin&&windows[n-1].centre+step<root.allowed.end)
            gaps.push_back({windows[first].centre,windows[n-1].centre+step,std::min(before,after)});
        for(auto f:{windows[first].centre,windows[n-1].centre+step})
            if(f>root.allowed.begin&&f<root.allowed.end)cuts.push_back({f,std::max(before,after),"Sustained low-need run boundary; preserve satisfactory material if complete measurement permits",f==windows[first].centre?EdgeDirection::falling:EdgeDirection::rising,f});
    }
    std::sort(cuts.begin(),cuts.end(),[](const auto& a,const auto& b){if(a.strength!=b.strength)return a.strength>b.strength;return a.frame<b.frame;});
    const bool directed=c.settings.schema==3;
    std::vector<std::pair<int64_t,EdgeDirection>> origins;
    std::vector<RegionCut> expanded;std::set<int64_t> seen;
    const auto expand=[&](const RegionCut& cut)
    {
        if(directed)
        {
            // Nearby median and run edges can be two observations of the same transition.
            if(std::any_of(origins.begin(),origins.end(),[&](auto origin){
                return origin.second==cut.direction&&std::abs(origin.first-cut.frame)<=step;}))return;
            origins.emplace_back(cut.frame,cut.direction);
        }
        for(auto delta:{int64_t(0),-int64_t(c.sampleRate)*3/2,int64_t(c.sampleRate)*3/2})
        {const auto f=cut.frame+delta;if(f<=root.allowed.begin||f>=root.allowed.end||!seen.insert(f).second)continue;expanded.push_back({f,cut.strength,cut.reason+(delta?"; +/-1.5 s window-support alternative":""),cut.direction,cut.frame});}
    };
    // Reserve a paired interior-gap hypothesis before the twelve-location cap.
    // Otherwise near-equal peaks at one edge can consume all slots and make the
    // promised two-cut spared-middle alternative impossible to even nominate.
    std::sort(gaps.begin(),gaps.end(),[](auto a,auto b){if(a.contrast!=b.contrast)return a.contrast>b.contrast;return a.begin<b.begin;});
    if(!gaps.empty())for(auto edge:{gaps.front().begin,gaps.front().end})expand({edge,gaps.front().contrast,"Reserved paired interior low-need boundary",edge==gaps.front().begin?EdgeDirection::falling:EdgeDirection::rising,edge});
    for(const auto& cut:cuts)expand(cut);
    return expanded;
}
juce::String geometryKey(const RegionGeometry& g)
{juce::String s=juce::String(g.support.begin)+":"+juce::String(g.support.end);for(auto f:g.cuts)s+=";"+juce::String(f);return s;}
}
RegionFormation formRegions(const RegionConstraints& c,const ProgramAnalysis& a,double ceiling,const std::function<bool()>& cancel,const std::vector<std::vector<RegionGeometry>>* preferred)
{
    RegionFormation out;out.alternatives.resize(c.roots.size());out.nominations.resize(c.roots.size());
    if(c.error.isNotEmpty()){out.reasons.push_back(c.error);return out;}
    for(size_t rootIndex=0;rootIndex<c.roots.size();++rootIndex)
    {
        if(cancel()){out.cancelled=true;return out;}const auto& root=c.roots[rootIndex];
        auto nominations=nominate(c,root,a,ceiling,cancel);if(cancel()){out.cancelled=true;return out;}
        if(preferred&&rootIndex<preferred->size())
        {
            std::vector<RegionCut> eventCuts;
            for(const auto& g:(*preferred)[rootIndex])for(auto f:{g.support.begin,g.support.end})
                if(eventCuts.size()<(c.settings.schema==3?6u:4u)&&std::none_of(eventCuts.begin(),eventCuts.end(),[&](auto x){return x.frame==f;}))eventCuts.push_back({f,1e5,c.settings.schema==3?"Event pair edge; permission and full-ramp validation required":"Supported event edge; permission and full-ramp validation still required",EdgeDirection::paired,f});
            const auto keep=12-eventCuts.size();if(nominations.size()>keep){out.omittedCuts+=nominations.size()-keep;nominations.resize(keep);}
            nominations.insert(nominations.begin(),eventCuts.begin(),eventCuts.end());
        }
        if(nominations.size()>12){out.omittedCuts+=nominations.size()-12;nominations.resize(12);}out.nominations[rootIndex]=nominations;
        const int maximum=c.settings.allowSplitting?c.settings.maximumSplits:0;
        std::vector<RegionCut> locations=nominations;
        for(auto mandatory:root.internalCuts)
        {
            locations.push_back({mandatory.minimum,1e6,"Explicit internal boundary range/pin"});
            if(mandatory.maximum!=mandatory.minimum)locations.push_back({mandatory.maximum,1e6,"Explicit internal boundary range endpoint"});
        }
        std::sort(locations.begin(),locations.end(),[](const auto& x,const auto& y){if(x.frame!=y.frame)return x.frame<y.frame;return x.strength>y.strength;});
        locations.erase(std::unique(locations.begin(),locations.end(),[](auto x,auto y){return x.frame==y.frame;}),locations.end());
        struct Cuts {std::vector<int64_t> frames;double score=0;};
        std::array<std::vector<Cuts>,5> combinations;
        for(int k=0;k<=maximum;++k)
        {
            std::vector<int64_t> frames;
            std::function<void(size_t,double)> enumerate=[&](size_t first,double score)
            {
                if(cancel())return;
                if(frames.size()==size_t(k))
                {
                    for(auto pin:root.internalCuts)if(std::none_of(frames.begin(),frames.end(),[&](auto f){return pin.contains(f);}))return;
                    combinations[size_t(k)].push_back({frames,score});return;
                }
                for(size_t n=first;n<locations.size();++n){frames.push_back(locations[n].frame);enumerate(n+1,score+locations[n].strength);frames.pop_back();}
            };
            enumerate(0,0);auto& list=combinations[size_t(k)];
            std::sort(list.begin(),list.end(),[](const auto& x,const auto& y){if(x.score!=y.score)return x.score>y.score;return x.frames<y.frames;});
            // This is a nomination ordering, not full candidate construction/search.
            if(list.size()>256){out.omittedGeometries+=list.size()-256;list.resize(256);}
        }
        std::vector<std::pair<int64_t,double>> starts{{root.start.minimum,0},{root.start.maximum,0}},ends{{root.end.minimum,0},{root.end.maximum,0}};
        if(c.settings.schema==3)
        {
            starts.clear();ends.clear();
            if(c.settings.start){starts.emplace_back(root.start.minimum,0);starts.emplace_back(root.start.maximum,0);}
            if(c.settings.end){ends.emplace_back(root.end.minimum,0);ends.emplace_back(root.end.maximum,0);}
        }
        if(root.start.contains(root.hint.begin))starts.push_back({root.hint.begin,0});if(root.end.contains(root.hint.end))ends.push_back({root.hint.end,0});
        if(c.settings.allowRefinement)for(auto cut:nominations){if((c.settings.schema<3||cut.direction==EdgeDirection::rising)&&root.start.contains(cut.frame))starts.push_back({cut.frame,cut.strength});if((c.settings.schema<3||cut.direction==EdgeDirection::falling)&&root.end.contains(cut.frame))ends.push_back({cut.frame,cut.strength});}
        struct Edges {RegionSpan span;double score=0;int64_t movement=0;};std::vector<Edges> edges;
        for(auto x:starts)for(auto y:ends)if(x.first<y.first)edges.push_back({{x.first,y.first},x.second+y.second,std::abs(x.first-root.hint.begin)+std::abs(y.first-root.hint.end)});
        std::sort(edges.begin(),edges.end(),[](auto x,auto y){if(x.score!=y.score)return x.score>y.score;if(x.movement!=y.movement)return x.movement<y.movement;if(x.span.begin!=y.span.begin)return x.span.begin<y.span.begin;return x.span.end<y.span.end;});
        std::set<juce::String> keys;size_t attempts=0;auto& retained=out.alternatives[rootIndex];
        const auto add=[&](RegionSpan bounds,const Cuts& cuts,const juce::String& reason,BoundaryFamily family=BoundaryFamily::coarse)
        {
            if(attempts>=256||retained.size()>=16)return;
            ++attempts;++out.constructions;
            if(!root.start.contains(bounds.begin)||!root.end.contains(bounds.end)||bounds.begin>=bounds.end)return;
            if(!cuts.frames.empty()&&(cuts.frames.front()<=bounds.begin||cuts.frames.back()>=bounds.end))return;
            // Retained topologies must have room for a nonzero outer decision and
            // either a legal nonzero interior decision or an unchanged gap. Exact
            // support legality still comes from the unchanged processor below.
            const int64_t minimum=int64_t(c.settings.minimumIntervalMs)*c.sampleRate/1000;
            const int64_t radius=std::llround(c.sampleRate*.05);
            std::vector<int64_t> nominal{bounds.begin+(bounds.begin?radius:0)};
            nominal.insert(nominal.end(),cuts.frames.begin(),cuts.frames.end());nominal.push_back(bounds.end-(bounds.end<c.frames?radius:0));
            for(size_t n=1;n<nominal.size();++n)
            {
                const bool outer=n==1||n+1==nominal.size();
                const auto least=outer?minimum:std::min(minimum,int64_t(c.settings.separationMs)*c.sampleRate/1000+2*radius);
                if(nominal[n]-nominal[n-1]<least)return;
            }
            RegionGeometry g;g.root=rootIndex;g.support=bounds;g.cuts=cuts.frames;g.strength=cuts.score;g.reason=reason;g.family=family;
            if(keys.insert(geometryKey(g)).second)retained.push_back(std::move(g));
        };
        if(root.internalCuts.empty())add(root.hint,{},"Unsplit original root, fixed-ramp support edges");
        if(preferred&&rootIndex<preferred->size())for(const auto& g:(*preferred)[rootIndex])
            if(g.cuts.size()<=size_t(maximum))add(g.support,{g.cuts,g.strength},g.reason,g.family);
        if(root.internalCuts.empty()&&!edges.empty())add(edges.front().span,{},"Unsplit outer refinement with strongest combined boundary nominations");
        // Reserve a two-cut topology before single-cut alternatives fill the list.
        if(maximum>=2)
        {
            std::optional<Cuts> spared;int64_t longest=0;
            for(const auto& cuts:combinations[2])
            {
                const int64_t begin=cuts.frames.front(),end=cuts.frames.back();double maximumNeed=0;size_t count=0;
                const int64_t minimum=int64_t(c.settings.minimumIntervalMs)*c.sampleRate/1000+std::llround(c.sampleRate*.05);
                if(begin-root.hint.begin<minimum||root.hint.end-end<minimum)continue;
                for(const auto& w:a.lraWindows){const auto mid=w.begin+(w.end-w.begin)/2;if(mid>=begin&&mid<end&&w.lufs&&a.p10){maximumNeed=std::max(maximumNeed,std::max(0.,*w.lufs-*a.p10-ceiling));++count;}}
                if(count&&maximumNeed<=c.settings.needThresholdTenths/20.&&end-begin>longest)
                {spared=cuts;longest=end-begin;}
            }
            if(spared)
            {
                add(root.hint,*spared,"Reserved two-cut unchanged-middle alternative, longest nominated low-need gap");
                auto inward=*spared;inward.frames[0]+=int64_t(c.sampleRate)*3/2;inward.frames[1]-=int64_t(c.sampleRate)*3/2;
                if(inward.frames[0]<inward.frames[1]&&std::all_of(inward.frames.begin(),inward.frames.end(),[&](auto f){return std::any_of(locations.begin(),locations.end(),[&](auto x){return x.frame==f;});}))
                    add(root.hint,inward,"Reserved two-cut unchanged-middle alternative with inward window-support hypotheses");
            }
        }
        std::array<size_t,5> cursors{};
        while(attempts<256&&retained.size()<16)
        {
            bool any=false;
            for(int k=1;k<=maximum&&attempts<256&&retained.size()<16;++k)
            {
                auto& cursor=cursors[size_t(k)];const auto& combos=combinations[size_t(k)];
                if(cursor>=combos.size())continue;any=true;const auto& cuts=combos[cursor++];
                add(root.hint,cuts,"Bounded split topology ordered by need-change strength and source frame");
                if(!edges.empty())add(edges.front().span,cuts,"Split topology with refined outer boundaries");
                if(cancel()){out.cancelled=true;return out;}
            }
            if(!any)break;
        }
        if(root.internalCuts.empty())for(const auto& e:edges){if(retained.size()>=16||attempts>=256)break;add(e.span,{},"Alternative exact outer processing support");}
        for(int k=1;k<=maximum;++k)out.omittedGeometries+=combinations[size_t(k)].size()-cursors[size_t(k)];
        if(retained.empty())out.reasons.push_back(root.id+": no nominated geometry fits the declared ranges; this is a bounded method limitation");
    }
    return out;
}
}

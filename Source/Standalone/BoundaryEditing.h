#pragma once
#include "ManualEditHistory.h"

namespace levelscope::standalone
{
struct BoundaryEditBatch
{
    ManualEditState before,next;
    juce::String runtime;
    uint64_t revision=0;
};
struct AppliedBoundaryDrag
{
    BoundaryEditBatch edit;
    size_t right=0;
    uint64_t token=0;
};
// Applied regions share a nominal cut. Move that cut only; retain both plateau
// gains and let the existing processor determine the centred transition radius.
inline std::optional<ManualEditState> moveAppliedBoundary(const ManualEditState& before,size_t right,int64_t frame)
{
    if(!before.gains||right==0||right>=before.gains->intervals.size())return {};
    auto next=before;auto& p=*next.gains;auto& left=p.intervals[right-1];auto& r=p.intervals[right];
    if(frame<left.begin+2||frame>r.end-2)return {};
    left.end=frame;r.begin=frame;
    if(p.validationError().isNotEmpty())return {};
    return next;
}
// Rebuild once from original sections and final edges. Retained interiors keep
// their gain changes at their original samples; extensions inherit edge plateaus.
// No intermediate mouse position is written into the applied processing plan.
inline std::optional<ManualEditState> moveBoundaryGains(const ManualEditState& before,const proposal::SectionMap& after)
{
    if(!after.valid()||after.frames!=before.map.frames||after.sampleRate!=before.map.sampleRate
        ||after.sections.size()!=before.map.sections.size()||after.inputIdentity!=before.map.inputIdentity)return {};
    auto next=before;next.map=after;if(!before.gains)return next;
    auto plan=*before.gains;plan.intervals.clear();const auto& old=before.gains->intervals;
    const auto plateau=[&](int64_t frame){auto i=std::upper_bound(old.begin(),old.end(),frame,[](auto f,const auto& r){return f<r.end;});return i==old.end()?0.:i->db;};
    const auto add=[&](int64_t a,int64_t b,double db){if(b>a)plan.intervals.push_back({a,b,db});};
    for(size_t n=0;n<after.sections.size();++n)
    {
        const auto& a=before.map.sections[n];const auto& b=after.sections[n];if(a.id!=b.id)return {};
        add(b.begin,std::min(b.end,a.begin),plateau(a.begin));
        auto i=std::upper_bound(old.begin(),old.end(),std::max(a.begin,b.begin),[](auto f,const auto& r){return f<r.end;});
        for(;i!=old.end()&&i->begin<std::min(a.end,b.end);++i)add(std::max({a.begin,b.begin,i->begin}),std::min({a.end,b.end,i->end}),i->db);
        add(std::max(b.begin,a.end),b.end,plateau(a.end-1));
    }
    plan=plan.canonical();if(plan.validationError().isNotEmpty())return {};next.gains=plan.isUnity()?std::nullopt:std::optional<processing::SectionGainPlan>(std::move(plan));return next;
}
}

#pragma once
#include "Core/Processing/SectionGainPlan.h"
#include "Core/Proposal/Model/ProgrammeSections.h"

namespace levelscope::standalone::section_edit
{
// Edit nominal plateaus only. The unchanged processor supplies boundary smoothing.
inline void split(processing::SectionGainPlan& plan,int64_t frame)
{
    for(size_t i=0;i<plan.intervals.size();++i)if(plan.intervals[i].begin<frame&&frame<plan.intervals[i].end)
    {auto right=plan.intervals[i];plan.intervals[i].end=frame;right.begin=frame;plan.intervals.insert(plan.intervals.begin()+ptrdiff_t(i+1),right);return;}
}
inline void setRange(processing::SectionGainPlan& plan,int64_t begin,int64_t end,double db)
{
    split(plan,begin);split(plan,end);
    for(auto& r:plan.intervals)if(r.begin>=begin&&r.end<=end)r.db=db;
}
inline std::optional<double> uniform(const processing::SectionGainPlan* plan,int64_t begin,int64_t end)
{
    if(!plan)return 0.;
    std::optional<double> value;
    for(const auto& r:plan->intervals)if(r.begin<end&&r.end>begin)
    {if(value&&*value!=r.db)return {};value=r.db;}
    return value;
}
inline proposal::SectionMap provenance(const processing::SectionGainPlan& plan,const proposal::SectionMap& map)
{
    auto result=map;result.sections.clear();result.nextId=1;
    for(const auto& r:plan.canonical().intervals)result.sections.push_back({result.nextId++,r.begin,r.end,"Manual gain region"});
    return result;
}
}

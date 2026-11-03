#pragma once
#include "Core/Processing/SectionGainPlan.h"
#include <array>

namespace levelscope::standalone
{
// Values at the same ramp endpoints as the processor. At display resolution
// retain the first, last, lowest and highest knot of each column in time order.
// Narrow sections survive even when several fit into a single pixel.
struct GainKnot {int64_t frame;double db;};
inline std::vector<GainKnot> gainOutline(const processing::SectionGainPlan* p,int64_t frames,
                                        double begin,double end,int pixels)
{
    std::vector<GainKnot> result;
    if(frames<=0||pixels<=0||end<=begin)return result;
    const int columns=std::min(pixels,4096);
    const auto first=std::clamp<int64_t>(int64_t(std::floor(begin)),0,frames-1);
    const auto last=std::clamp<int64_t>(int64_t(std::ceil(end)),first,frames-1);
    std::array<GainKnot,4> bucket{};int column=-1;
    const auto flush=[&]
    {
        if(column<0)return;
        std::sort(bucket.begin(),bucket.end(),[](auto a,auto b){return a.frame<b.frame;});
        for(auto v:bucket)if(result.empty()||result.back().frame<v.frame)result.push_back(v);
    };
    const auto add=[&](int64_t frame,double db)
    {
        const GainKnot knot{frame,db};
        const int index=std::clamp(int((frame-begin)*columns/(end-begin)),0,columns-1);
        if(index!=column){flush();bucket.fill(knot);column=index;}
        else {bucket[1]=knot;if(db<bucket[2].db)bucket[2]=knot;if(db>bucket[3].db)bucket[3]=knot;}
    };
    add(first,p?p->gainAt(first):0);
    if(p)
    {
        const auto here=std::upper_bound(p->intervals.begin(),p->intervals.end(),first,[](int64_t f,const auto& r){return f<r.end;});
        for(size_t i=std::max<size_t>(1,size_t(here-p->intervals.begin()));i<p->intervals.size();++i)
        {
            const auto b=p->intervals[i].begin,r=p->radius(i);
            if(b-r>last)break;
            if(b-r>first)add(b-r,p->intervals[i-1].db);
            if(b+r-1>first&&b+r-1<last)add(b+r-1,p->intervals[i].db);
        }
    }
    add(last,p?p->gainAt(last):0);
    if(last==frames-1)add(frames,p?p->intervals.back().db:0);
    flush();return result;
}
}

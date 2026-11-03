#include "WaveformOverview.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace levelscope::waveform
{
Extrema merge(const Extrema& a,const Extrema& b) noexcept
{
    auto result=a;
    if(b.minimum<a.minimum || (b.minimum==a.minimum && b.minimumFrame<a.minimumFrame))
    {result.minimum=b.minimum;result.minimumFrame=b.minimumFrame;}
    if(b.maximum>a.maximum || (b.maximum==a.maximum && b.maximumFrame<a.maximumFrame))
    {result.maximum=b.maximum;result.maximumFrame=b.maximumFrame;}
    return result;
}
Budget::Budget(size_t limit):limitBytes(limit)
{if(limit<displayReserveBytes || limit>combinedBudgetBytes)throw std::invalid_argument("Invalid combined waveform budget");}
bool Budget::reserve(size_t amount) noexcept
{
    auto current=bytes.load();
    do {if(amount>limitBytes-current)return false;}
    while(!bytes.compare_exchange_weak(current,current+amount));
    return true;
}
std::shared_ptr<Overview> Overview::create(Spec spec,std::shared_ptr<Budget> budget)
{
    if(!budget || spec.frames<=0 || spec.sampleRate<=0 || spec.channels<1 || spec.channels>channelCapacity)
        throw std::invalid_argument("Invalid waveform programme");
    Layout plan;
    const size_t allowance=std::min(maximumOverviewBytes,budget->remaining());
    for(;;)
    {
        const auto count=uint64_t((spec.frames-1)/plan.width)+1;
        const auto maximumNodes=allowance/(sizeof(Extrema)*size_t(spec.channels));
        if(count<=maximumNodes)
        {
            plan.nodes=0;plan.levels=0;
            for(size_t n=size_t(count);;n=(n+1)/2)
            {
                if(plan.levels==int(plan.level.size()) || n>maximumNodes-plan.nodes){plan.nodes=maximumNodes+1;break;}
                plan.level[size_t(plan.levels++)]={n,plan.nodes};plan.nodes+=n;
                if(n==1)break;
            }
            const size_t overhead=sizeof(Overview)+512; // shared ownership + allocation bookkeeping
            if(plan.nodes<=maximumNodes && allowance>=overhead
                && plan.nodes<=(allowance-overhead)/(sizeof(Extrema)*size_t(spec.channels)))
            {plan.bytes=overhead+plan.nodes*sizeof(Extrema)*size_t(spec.channels);break;}
        }
        if(plan.width>=spec.frames || plan.width>std::numeric_limits<int64_t>::max()/2)return {};
        plan.width*=2;
    }
    if(!budget->reserve(plan.bytes))return {};
    Overview* raw=nullptr;
    try {raw=new Overview(spec,budget,plan);}
    catch(...) {budget->release(plan.bytes);throw;}
    // shared_ptr deletes raw if its control-block allocation fails.
    return std::shared_ptr<Overview>(raw);
}
Overview::Overview(Spec s,std::shared_ptr<Budget> b,const Layout& p)
    :spec(s),binFrames(p.width),storageBytes(p.bytes),budget(std::move(b)),layout(p),
     bins(new Extrema[p.nodes*size_t(s.channels)]) {}
Overview::~Overview(){budget->release(storageBytes);}
int64_t Overview::frameAt(size_t bin) const noexcept
{
    if(bin>=layout.level[0].bins)return spec.frames;
    return int64_t(bin)*binFrames;
}
void Overview::seal() noexcept
{
    for(int c=0;c<spec.channels;++c)bins[sealed*size_t(spec.channels)+size_t(c)]=pending[size_t(c)];
    size_t index=sealed++;
    for(int level=1;level<layout.levels;++level)
    {
        const auto& previous=layout.level[size_t(level-1)];
        // Seal complete sibling pairs; EOF also seals the shortened right tail.
        if((index&1)==0 && index+1<previous.bins)break;
        const size_t parent=index/2, first=parent*2;
        for(int c=0;c<spec.channels;++c)
        {
            const auto a=bins[(previous.offset+first)*size_t(spec.channels)+size_t(c)];
            bins[(layout.level[size_t(level)].offset+parent)*size_t(spec.channels)+size_t(c)]=
                first+1<previous.bins?merge(a,bins[(previous.offset+first+1)*size_t(spec.channels)+size_t(c)]):a;
        }
        index=parent;
    }
}
void Overview::publish(bool final) noexcept
{
    const auto now=std::chrono::steady_clock::now();
    if(final || now-lastPublication>=std::chrono::milliseconds(200))
    {published.store((uint64_t(sealed)<<1)|(final?1u:0u),std::memory_order_release);lastPublication=now;}
}
bool Overview::append(const float* const* audio,int channels,int count,int64_t first) noexcept
{
    if(failed || finished || !audio || channels!=spec.channels || count<0 || first!=consumed || count>spec.frames-consumed)
    {failed=true;return false;}
    for(int c=0;c<channels;++c)if(!audio[c]){failed=true;return false;}
    for(int n=0;n<count;)
    {
        const auto offset=consumed%binFrames;
        const int take=int(std::min<int64_t>(count-n,binFrames-offset));
        for(int c=0;c<channels;++c)
        {
            auto& extrema=pending[size_t(c)];
            for(int i=0;i<take;++i)
            {
                const float x=audio[c][n+i];const int64_t frame=consumed+i;
                if(!std::isfinite(x)){failed=true;return false;}
                if(offset==0 && i==0)extrema={x,x,frame,frame};
                else {if(x<extrema.minimum){extrema.minimum=x;extrema.minimumFrame=frame;}
                      if(x>extrema.maximum){extrema.maximum=x;extrema.maximumFrame=frame;}}
            }
        }
        consumed+=take;n+=take;
        if(consumed%binFrames==0 || consumed==spec.frames)seal();
    }
    publish(false);return true;
}
bool Overview::finish() noexcept
{
    if(failed || finished || consumed!=spec.frames)return false;
    finished=true;publish(true);return true;
}
int64_t Overview::coveredFrames() const noexcept {return coveredFrames(publication());}
int64_t Overview::coveredFrames(uint64_t snapshot) const noexcept {return frameAt(size_t(std::min(snapshot,publication())>>1));}
Extrema Overview::query(size_t begin,size_t end,int channel) const noexcept
{
    Extrema result{};bool first=true;
    while(begin<end)
    {
        size_t span=1;int level=0;
        while(level+1<layout.levels && begin%(span*2)==0 && span*2<=end-begin){span*=2;++level;}
        const auto& value=bins[(layout.level[size_t(level)].offset+begin/span)*size_t(spec.channels)+size_t(channel)];
        result=first?value:merge(result,value);first=false;begin+=span;
    }
    return result;
}
std::vector<Column> Overview::reduce(int64_t begin,int64_t end,int columns,uint64_t snapshot) const
{
    std::vector<Column> result;
    begin=std::clamp<int64_t>(begin,0,spec.frames);end=std::clamp<int64_t>(end,0,spec.frames);
    if(end<=begin || columns<1)return result;
    const size_t first=size_t(begin/binFrames),last=size_t((end-1)/binFrames)+1,span=last-first;
    const size_t count=std::min<size_t>(std::min(columns,maximumColumns),span);
    result.resize(count);
    const size_t ready=size_t(std::min(snapshot,publication())>>1);
    for(size_t n=0;n<count;++n)
    {
        const size_t a=first+span*n/count,b=first+span*(n+1)/count,covered=std::min(b,ready);
        auto& column=result[n];column.begin=frameAt(a);column.end=frameAt(b);column.coveredEnd=frameAt(std::max(a,covered));
        if(covered>a)for(int c=0;c<spec.channels;++c)column.channels[size_t(c)]=query(a,covered,c);
    }
    return result;
}
std::optional<Column> Overview::inspect(int64_t frame) const
{
    if(frame<0 || frame>=spec.frames)return {};
    auto columns=reduce(frame,frame+1,1);
    return columns.empty()?std::nullopt:std::optional<Column>(columns.front());
}
}

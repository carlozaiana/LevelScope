#include "StandaloneLoudnessPreview.h"
#include <cmath>
#include <stdexcept>

namespace levelscope::standalone
{
LoudnessPreview::LoudnessPreview(int64_t f,int rate):frames(f),sampleRate(rate)
{
    if(frames<=0 || (rate!=44100&&rate!=48000&&rate!=96000&&rate!=192000))throw std::invalid_argument("Invalid loudness preview dimensions");
    capacity={count(false,frames),count(true,frames)};
    if(capacity[0]+capacity[1]>maximumBytes/sizeof(measurement::Window))throw std::length_error("Loudness preview exceeds bounded storage");
    for(size_t n=0;n<2;++n)windows[n]=std::make_unique<measurement::Window[]>(capacity[n]);
}
size_t LoudnessPreview::count(bool shortTerm,int64_t covered) const
{return size_t(std::max<int64_t>(0,covered/(sampleRate/10)-(shortTerm?30:4)+1));}
void LoudnessPreview::window(bool shortTerm,const measurement::Window& value) noexcept
{
    const size_t n=shortTerm?1:0;
    if(failed)return;
    if(written[n]>=capacity[n] || value.endFrame!=int64_t(written[n]+(shortTerm?30:4))*(sampleRate/10)
        || !std::isfinite(value.energy) || value.energy<0){failed=true;return;}
    windows[n][written[n]++]=value;
}
void LoudnessPreview::publish(int64_t covered,bool force) noexcept
{
    if(failed || covered<coverage() || covered>frames || count(false,covered)>written[0] || count(true,covered)>written[1])return;
    const auto now=std::chrono::steady_clock::now();
    if(force || now-lastPublication>=std::chrono::milliseconds(200))
    {published.store(covered,std::memory_order_release);lastPublication=now;}
}
proposal::TimelinePlot LoudnessPreview::reduce(double begin,double end,int columns,int64_t covered) const
{
    proposal::TimelinePlot plot;plot.begin=std::clamp(begin,0.,double(frames));plot.end=std::clamp(end,plot.begin,double(frames));
    if(!std::isfinite(begin)||!std::isfinite(end)||plot.end<=plot.begin||columns<1)return plot;
    columns=std::min(4096,columns);plot.momentary.resize(size_t(columns));plot.shortTerm.resize(size_t(columns));
    covered=std::min(covered,coverage());
    for(bool s:{false,true})
    {
        const size_t length=count(s,covered);if(length==0)continue;
        const auto* first=windows[s?1:0].get();auto it=std::lower_bound(first,first+length,plot.begin,[](const auto& w,double frame){return w.endFrame<frame;});
        auto& buckets=s?plot.shortTerm:plot.momentary;
        for(;it!=first+length && it->endFrame<=plot.end;++it)
        {
            ++plot.fullResolutionPoints;if(it->energy==0)continue; // actual silence is a gap, never 0 LUFS
            const auto index=size_t(std::min(columns-1,int((it->endFrame-plot.begin)*columns/(plot.end-plot.begin))));
            buckets[index].ungated.add(measurement::OfflineMeasurement::loudness(it->energy),it->endFrame);
        }
    }
    return plot;
}
std::optional<measurement::Window> LoudnessPreview::inspect(bool s,double frame,int64_t covered) const
{
    covered=std::min(covered,coverage());const auto length=count(s,covered);
    if(!std::isfinite(frame)||frame<0||frame>covered||length==0)return {};
    const auto* first=windows[s?1:0].get();const auto* it=std::lower_bound(first,first+length,frame,[](const auto& w,double f){return w.endFrame<f;});
    if(it==first+length)return first[length-1];
    if(it!=first&&frame-(it-1)->endFrame<=it->endFrame-frame)--it;
    return *it;
}
std::shared_ptr<LoudnessPreview> LoudnessProgress::begin(int64_t frames,int rate) noexcept
{
    try
    {
        auto next=std::make_shared<LoudnessPreview>(frames,rate);std::lock_guard<std::mutex> lock(mutex);
        if(abandoned)return {};data=next;return next;
    }
    catch(...){return {};}
}
std::shared_ptr<const LoudnessPreview> LoudnessProgress::view(const WaveformBinding& expected) const
{std::lock_guard<std::mutex> lock(mutex);return !abandoned&&binding.sameJob(expected)?data:nullptr;}
void LoudnessProgress::abandon(){std::lock_guard<std::mutex> lock(mutex);abandoned=true;data.reset();}
}

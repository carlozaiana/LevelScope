#include "SectionPeakMap.h"
#include <cmath>
#include <stdexcept>

namespace levelscope::proposal
{
namespace
{
LocatedPeak maximum(LocatedPeak a, LocatedPeak b)
{
    return b.magnitude > a.magnitude || (b.magnitude == a.magnitude && b.eighthFrame < a.eighthFrame) ? b : a;
}
}
juce::String sectionGeometryKey(const SectionMap& map)
{
    juce::String key=map.inputIdentity+":"+juce::String(map.sampleRate)+":"+juce::String(map.frames);
    for(const auto& row:map.sections) key += ";"+juce::String(row.id)+":"+juce::String(row.begin)+":"+juce::String(row.end);
    return key;
}
size_t SectionPeakMap::storageBytes() const
{
    size_t bytes=sizeof(*this)+geometry.getNumBytesAsUTF8()+snapshotIdentity.getNumBytesAsUTF8()+sections.capacity()*sizeof(SectionPeaks);
    for(const auto& s:sections) bytes+=(s.sample.capacity()+s.truePeak.capacity())*sizeof(LocatedPeak);
    return bytes;
}
SectionPeakCollector::SectionPeakCollector(const SectionMap& m,int channels,const juce::String& snapshot):map(m)
{
    if(!map.valid()||channels<1||channels>measurement::channelCapacity||snapshot.isEmpty())
        throw std::invalid_argument("Invalid continuous section-peak request");
    result.geometry=sectionGeometryKey(map);result.snapshotIdentity=snapshot;result.channels=channels;
    const int64_t margin=int64_t(std::llround(map.sampleRate*.050))+64;
    edges={0,map.frames*8};
    for(const auto& s:map.sections)
    {
        const auto lo=std::max<int64_t>(0,s.begin-margin)*8;
        const auto hi=std::min(map.frames,s.end+margin)*8;
        neighbourhoods.emplace_back(lo,hi);edges.push_back(lo);edges.push_back(hi);
        result.sections.push_back({std::vector<LocatedPeak>(size_t(channels)),std::vector<LocatedPeak>(size_t(channels))});
    }
    std::sort(edges.begin(),edges.end());edges.erase(std::unique(edges.begin(),edges.end()),edges.end());
    buckets.resize((edges.size()-1)*size_t(channels));
}
void SectionPeakCollector::sample(size_t channel,double magnitude,int64_t frame)
{
    if(finished||channel>=size_t(result.channels)||!std::isfinite(magnitude)||magnitude<0||frame<0||frame>=map.frames)
        throw std::invalid_argument("Invalid section sample observation");
    auto& index=sampleCursor[channel];
    if(frame<map.sections[index].begin||frame>=map.sections[index].end)
        index=size_t(std::upper_bound(map.sections.begin(),map.sections.end(),frame,
            [](auto f,const auto& s){return f<s.end;})-map.sections.begin());
    auto& p=result.sections[index].sample[channel];p=maximum(p,{magnitude,frame*8,false});
}
void SectionPeakCollector::truePeak(size_t channel,double magnitude,int64_t frame,bool padding)
{
    if(finished||channel>=size_t(result.channels)||!std::isfinite(magnitude)||magnitude<0)
        throw std::invalid_argument("Invalid section true-peak observation");
    // Clamp only the attribution lookup. The reported location is never clamped.
    const auto coordinate=std::clamp<int64_t>(frame,0,map.frames*8-1);
    auto& index=peakCursor[channel];
    if(coordinate<edges[index]||coordinate>=edges[index+1])
        index=size_t(std::upper_bound(edges.begin(),edges.end(),coordinate)-edges.begin()-1);
    auto& p=buckets[index*size_t(result.channels)+channel];
    p=maximum(p,{magnitude,frame,padding||frame<0||frame>=map.frames*8});
}
std::shared_ptr<const SectionPeakMap> SectionPeakCollector::finish()
{
    if(finished)throw std::logic_error("Section peak collector already finished");
    finished=true;size_t base=1;while(base<edges.size()-1)base*=2;
    // One channel at a time: bounded storage and O(n log n), including 10,000
    // very short sections with heavily overlapping attribution neighbourhoods.
    for(size_t ch=0;ch<size_t(result.channels);++ch)
    {
        std::vector<LocatedPeak> tree(base*2);
        for(size_t n=0;n<edges.size()-1;++n)tree[base+n]=buckets[n*size_t(result.channels)+ch];
        for(size_t n=base-1;n>0;--n)tree[n]=maximum(tree[n*2],tree[n*2+1]);
        for(size_t n=0;n<map.sections.size();++n)
        {
            size_t lo=base+size_t(std::lower_bound(edges.begin(),edges.end(),neighbourhoods[n].first)-edges.begin());
            size_t hi=base+size_t(std::lower_bound(edges.begin(),edges.end(),neighbourhoods[n].second)-edges.begin());
            LocatedPeak peak;
            while(lo<hi){if(lo&1)peak=maximum(peak,tree[lo++]);if(hi&1)peak=maximum(peak,tree[--hi]);lo/=2;hi/=2;}
            result.sections[n].truePeak[ch]=peak;
        }
    }
    return std::make_shared<const SectionPeakMap>(std::move(result));
}
}

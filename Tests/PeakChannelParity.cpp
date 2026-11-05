// Exact serial oracle for channel-parallel peak processing, independent of JUCE.
#include "Core/Measurement/OfflineMeasurement.h"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <thread>

using namespace levelscope::measurement;
namespace
{
void check (bool ok, const char* what) { if (!ok) throw std::runtime_error (what); }
bool equal (double a, double b) { return std::memcmp (&a, &b, sizeof a) == 0; }
bool equal (const std::optional<double>& a, const std::optional<double>& b)
{ return bool(a) == bool(b) && (!a || equal (*a,*b)); }
void compare (const Metric& a, const Metric& b)
{ check (a.status == b.status && equal (a.value,b.value), "metric value/status"); }
void compare (const std::vector<Window>& a, const std::vector<Window>& b)
{
    check (a.size()==b.size(),"window count");
    for (size_t i=0; i<a.size(); ++i)
        check (a[i].endFrame==b[i].endFrame && equal(a[i].energy,b[i].energy),"window frame/energy");
}
void compare (const Result& a, const Result& b)
{
    check(a.sampleRate==b.sampleRate && a.frames==b.frames && a.integratedIncluded==b.integratedIncluded
          && a.lraIncluded==b.lraIncluded && a.lraStable==b.lraStable,"coverage/gate counts");
    compare(a.integrated,b.integrated); compare(a.lra,b.lra); compare(a.truePeak,b.truePeak);
    compare(a.samplePeak,b.samplePeak); compare(a.momentaryMax,b.momentaryMax); compare(a.shortTermMax,b.shortTermMax);
    check(equal(a.integratedRelativeGate,b.integratedRelativeGate) && equal(a.lraRelativeGate,b.lraRelativeGate)
          && equal(a.p10,b.p10) && equal(a.p95,b.p95),"gates/percentiles");
    compare(a.momentary,b.momentary); compare(a.shortTerm,b.shortTerm); compare(a.lraWindows,b.lraWindows);
    check(a.samplePeakFrames==b.samplePeakFrames && a.truePeakEighthFrames==b.truePeakEighthFrames,"peak positions");
    for(size_t c=0;c<a.channelTruePeaks.size();++c)
        check(equal(a.channelSamplePeaks.at(c),b.channelSamplePeaks.at(c)) &&
              equal(a.channelTruePeaks.at(c),b.channelTruePeaks.at(c)),"channel peak bits");
}
struct Windows final : WindowObserver
{
    std::vector<Window> m,s;
    void window(bool shortTerm, const Window& w) noexcept override { (shortTerm?s:m).push_back(w); }
};
struct Event
{
    bool sample, padding;
    int64_t at;
    double magnitude;
};
struct Events final : PeakObserver
{
    explicit Events(bool allowed):parallel(allowed),caller(std::this_thread::get_id()) {}
    bool acceptsConcurrentChannels() const noexcept override { return parallel; }
    void sample(size_t ch,double v,int64_t f) override { append(ch,{true,false,f,v}); }
    void truePeak(size_t ch,double v,int64_t f,bool pad) override { append(ch,{false,pad,f,v}); }
    void append(size_t ch,Event e)
    {
        if(std::this_thread::get_id()!=caller) offCaller[ch]=true;
        perChannel[ch].push_back(e);
        if(!parallel) sequence.push_back(ch);
    }
    bool parallel;
    std::thread::id caller;
    std::array<bool,channelCapacity> offCaller{};
    std::array<std::vector<Event>,channelCapacity> perChannel;
    std::vector<size_t> sequence;
};
using Audio=std::vector<std::vector<float>>;
Audio programme(size_t channels,int frames,bool silence)
{
    Audio a(channels,std::vector<float>(frames));
    if(silence) return a;
    for(size_t ch=0;ch<channels;++ch)
    {
        for(int n=0;n<frames;++n)
            a[ch][n]=float((.24+.055*ch)*std::sin(.23*n+.17*ch)*(.65+.35*std::cos(.000027*n)));
        for(int n:{0,255,4095,4096,frames-1})
            if(n>=0 && n<frames) a[ch][n]=(ch%2 ? -1.f : 1.f);
    }
    return a;
}
void pair(int rate, const std::vector<Role>& roles, int frames, bool silence,
          bool observe=false, bool allowConcurrent=true)
{
    auto a=programme(roles.size(),frames,silence);
    Windows leftWindows,rightWindows; Events left(allowConcurrent),right(allowConcurrent);
    OfflineMeasurement serial(rate,roles,observe?&left:nullptr,&leftWindows,PeakExecution::serial);
    OfflineMeasurement parallel(rate,roles,observe?&right:nullptr,&rightWindows,PeakExecution::boundedChannels);
    const std::array<int,7> sizes{1,255,256,4096,17,1023,4096};
    std::array<const float*,channelCapacity> p{};
    int pos=0;
    for(size_t step=0;pos<frames;++step)
    {
        const int count=std::min(frames-pos,sizes[step%sizes.size()]);
        for(size_t ch=0;ch<roles.size();++ch) p[ch]=a[ch].data()+pos;
        serial.push(p.data(),count); parallel.push(p.data(),count);
        pos+=count;
    }
    compare(serial.finish(),parallel.finish());
    compare(serial.finish(),parallel.finish()); // idempotent drain
    compare(leftWindows.m,rightWindows.m); compare(leftWindows.s,rightWindows.s);
    if(observe)
    {
        for(size_t ch=0;ch<roles.size();++ch)
        {
            auto& x=left.perChannel[ch]; auto& y=right.perChannel[ch];
            check(x.size()==y.size(),"channel event count");
            for(size_t n=0;n<x.size();++n)
                check(x[n].at==y[n].at && x[n].padding==y[n].padding && x[n].sample==y[n].sample
                      && equal(x[n].magnitude,y[n].magnitude),"event including interpolation phase and EOF padding");
            if(!silence) check(parallel.finish().samplePeakFrames[ch]==0,"earliest equal-magnitude sample tie");
        }
        if(!allowConcurrent)
            check(left.sequence==right.sequence && std::none_of(right.offCaller.begin(),right.offCaller.end(),[](bool x){return x;}),
                  "legacy observer callback order and thread");
        if(roles.size()!=2 && roles.size()!=6)
            check(std::none_of(right.offCaller.begin(),right.offCaller.end(),[](bool x){return x;}),
                  "mono and other layouts retain serial execution");
    }
}
void failures()
{
    auto a=programme(2,4096,false);
    const float* p[]{a[0].data(),a[1].data()};
    struct Throwing : PeakObserver
    {
        bool acceptsConcurrentChannels() const noexcept override { return true; }
        void sample(size_t ch,double,int64_t frame) override
        { if(ch==1 && frame==137) throw std::runtime_error("observer error"); }
        void truePeak(size_t,double,int64_t,bool) override {}
    } observer;
    OfflineMeasurement monitored(48000,{Role::left,Role::right},&observer,nullptr,PeakExecution::boundedChannels);
    bool failed=false;
    try { monitored.push(p,4096); } catch(const std::runtime_error& e) { failed=std::string(e.what())=="observer error"; }
    check(failed,"worker error delivered to job owner");
    for(auto mode:{PeakExecution::serial,PeakExecution::boundedChannels})
    {
        a[1][1024]=std::numeric_limits<float>::infinity();
        OfflineMeasurement engine(48000,{Role::left,Role::right},nullptr,nullptr,mode);
        failed=false;
        try { engine.push(p,4096); } catch(const std::invalid_argument&) { failed=true; }
        check(failed,"non-finite input rejected");
    }
}
}
int main(int argc,char** argv)
{
    try
    {
        const std::vector<Role> stereo{Role::right,Role::left};
        const std::vector<Role> surround{Role::lfe,Role::left,Role::centre,Role::surroundRight,Role::right,Role::surroundLeft};
        if(argc==2 && std::string(argv[1])=="--thread-audit")
        {
            pair(48000,surround,8209,false,true);
            pair(48000,stereo,8209,false,true,false);
            failures();
            std::cout << "PASS: peak callback, filter history, fallback and failure paths.\n";
            return 0;
        }
        for(int rate:{44100,48000,96000,192000})
            for(const auto& roles:{stereo,surround})
            {
                pair(rate,roles,rate*3+rate/2+137,false);
                pair(rate,roles,8209,false,true);
                pair(rate,roles,1024,true,true);
            }
        pair(48000,stereo,8209,false,true,false);
        pair(48000,surround,0,true,true);
        pair(48000,stereo,1,false,true);
        pair(48000,{Role::centre},48017,false,true);
        pair(48000,{Role::left,Role::right,Role::centre,Role::surroundLeft,Role::surroundRight},48017,false,true);
        auto eight=surround; eight.push_back(Role::rearLeft); eight.push_back(Role::rearRight);
        pair(48000,eight,48017,false,true);
        failures();
        const auto original=std::fegetround();
        for(auto rounding:{FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO,FE_TONEAREST})
        {
            check(std::fesetround(rounding)==0,"set rounding mode");
            pair(48000,surround,8209,false,true);
            check(std::fegetround()==rounding,"worker changed caller rounding mode");
        }
        std::fesetround(original);
        std::cout << "PASS: all result bits, windows, channel peaks and sample/eighth-frame locations; observer events, fallback and errors.\n";
        return 0;
    }
    catch(const std::exception& e) { std::cerr << "FAIL: " << e.what() << '\n'; return 1; }
}

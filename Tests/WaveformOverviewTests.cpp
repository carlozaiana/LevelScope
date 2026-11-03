#include "Core/Waveform/WaveformOverview.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <random>
#include <stdexcept>
#include <thread>
using namespace levelscope::waveform;
namespace
{
int checks=0;
void expect(bool ok,const char* text){++checks;if(!ok)throw std::runtime_error(text);}
bool same(const Extrema& a,const Extrema& b)
{return a.minimum==b.minimum&&a.maximum==b.maximum&&a.minimumFrame==b.minimumFrame&&a.maximumFrame==b.maximumFrame;}
Extrema exact(const std::vector<float>& x,int64_t begin,int64_t end)
{
    Extrema e{x[size_t(begin)],x[size_t(begin)],begin,begin};
    for(auto f=begin+1;f<end;++f)e=merge(e,{x[size_t(f)],x[size_t(f)],f,f});
    return e;
}
void oracle()
{
    for(int channels:{1,2,6,8,12})for(int count:{1,63,64,65,129,4097,14013})
    {
        std::vector<std::vector<float>> samples(size_t(channels),std::vector<float>(size_t(count),0));
        std::mt19937 rng(729);
        for(int c=0;c<channels;++c)for(int f=0;f<count;++f)samples[size_t(c)][size_t(f)]=float(int(rng()%31)-15)/4;
        for(int block:{1,31,64,4096})
        {
            auto budget=std::make_shared<Budget>();auto overview=Overview::create({count,48000,channels},budget);
            expect(bool(overview),"allocate bounded overview");
            std::array<const float*,channelCapacity> pointers{};
            for(int f=0;f<count;f+=block)
            {for(int c=0;c<channels;++c)pointers[size_t(c)]=samples[size_t(c)].data()+f;expect(overview->append(pointers.data(),channels,std::min(block,count-f),f),"append actual samples");}
            expect(overview->finish()&&overview->complete()&&overview->coveredFrames()==count,"exact odd EOF coverage");
            for(int columns:{1,2,7,31,4096,100000})for(auto bounds:{std::pair<int64_t,int64_t>{0,count},{count/3,count}})
            {
                const auto reduced=overview->reduce(bounds.first,bounds.second,columns);expect(reduced.size()<=maximumColumns,"bounded display columns");
                for(const auto& col:reduced)for(int c=0;c<channels;++c)
                    expect(col.available()&&col.end<=count&&same(col.channels[size_t(c)],exact(samples[size_t(c)],col.begin,col.coveredEnd)),"every pyramid/display reduction preserves extrema and exact positions");
            }
            const auto tail=overview->inspect(count-1);expect(tail&&tail->end==count&&tail->begin==(count-1)/overview->binFrames*overview->binFrames,"finest interval shortened at EOF");
            expect(!overview->inspect(count)&&!overview->inspect(-1),"no invented samples outside programme");
        }
    }
}
void memoryAndThreads()
{
    auto budget=std::make_shared<Budget>();
    auto original=Overview::create({3LL*3600*192000,192000,12},budget);
    auto current=Overview::create({3LL*3600*192000,192000,12},budget);
    expect(original&&current&&original->binFrames>64,"long twelve-channel overview adapts finest resolution");
    expect(original->storageBytes<=maximumOverviewBytes&&current->storageBytes<=maximumOverviewBytes&&budget->used()<=combinedBudgetBytes,"all levels plus Source/Current and display reservation within 128 MiB");
    auto shared=original;const auto used=budget->used();expect(used==budget->used(),"sharing immutable overview uses no second slab");
    original.reset();expect(budget->used()==used,"retained reference stays charged");shared.reset();current.reset();expect(budget->used()==displayReserveBytes,"retired storage releases budget");
    auto noRoom=std::make_shared<Budget>(displayReserveBytes);expect(!Overview::create({128,48000,2},noRoom),"budget refusal is explicit");
    auto wave=Overview::create({300003,48000,2},budget);std::atomic<bool> done{false};std::atomic<int> failures{0};
    std::thread reader([&]{while(!done.load())for(const auto& c:wave->reduce(0,wave->spec.frames,73))if(c.available())
        if(c.channels[0].minimumFrame!=c.begin||c.channels[0].maximumFrame!=c.begin||c.channels[0].minimum!=.75f
            ||c.channels[1].minimumFrame!=c.begin||c.channels[1].minimum!=-.75f)++failures;});
    float a[137],b[137];std::fill_n(a,137,.75f);std::fill_n(b,137,-.75f);const float* data[]={a,b};
    for(int64_t f=0;f<wave->spec.frames;f+=137)expect(wave->append(data,2,int(std::min<int64_t>(137,wave->spec.frames-f)),f),"sequential progressive append");
    expect(wave->finish(),"finish concurrent publication");done=true;reader.join();expect(failures==0,"acquired published prefixes are immutable and keep anti-phase channels separate");
    auto col=wave->reduce(0,wave->spec.frames,1).front();expect(col.channels[0].minimumFrame==0&&col.channels[1].maximumFrame==0,"equal extrema keep earliest frame through all reductions");
    auto partial=Overview::create({128,48000,1},budget);const float* mono[]={a};expect(partial->append(mono,1,64,0),"one complete bin published");
    auto missing=partial->inspect(100);expect(missing&&!missing->available()&&!partial->complete(),"unprocessed passage stays pending, not zero");
    expect(!partial->append(mono,1,2,65)&&!partial->finish(),"discontinuous input cannot finish");
}
void publicationAndInvalidInput()
{
    auto budget=std::make_shared<Budget>();auto wave=Overview::create({256,48000,1},budget);
    float silence[64]{};const float* samples[]={silence};
    expect(wave->append(samples,1,64,0),"first sealed silence bin");const auto first=wave->publication();
    expect(wave->append(samples,1,64,64)&&wave->publication()==first,"intermediate publications coalesce within 200 ms");
    std::this_thread::sleep_for(std::chrono::milliseconds(210));
    expect(wave->append(samples,1,64,128)&&wave->coveredFrames()==192,"next publication includes every newly completed bin");
    expect(wave->append(samples,1,64,192)&&wave->finish()&&wave->coveredFrames()==256,"final publication is immediate even within throttle interval");
    const auto silent=wave->inspect(0);expect(silent&&silent->available()&&silent->channels[0].minimum==0&&silent->channels[0].maximum==0,"actual digital silence is covered PCM");
    auto bad=Overview::create({64,48000,1},budget);silence[13]=std::numeric_limits<float>::infinity();
    expect(!bad->append(samples,1,64,0)&&!bad->finish()&&bad->coveredFrames()==0,"nonfinite sample never publishes partial bin or completes");
}
}
int main()
{
    try{oracle();memoryAndThreads();publicationAndInvalidInput();std::cout<<"PASS "<<checks<<" waveform core checks\n";return 0;}
    catch(const std::exception& e){std::cerr<<"FAIL "<<checks<<" "<<e.what()<<'\n';return 1;}
}

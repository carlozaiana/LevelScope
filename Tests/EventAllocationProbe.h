// Qualification executable only. Observe Darwin's allocation callbacks without
// allocating in the callback. Charge allocator block sizes, including JUCE,
// vectors and transient reallocations, across every thread during a scope.
#pragma once
#if defined(__APPLE__)
#include <malloc/malloc.h>
#include <atomic>
#include <array>
#include <cstdint>
#include <stdexcept>
namespace event_allocation
{
using Logger=void(uint32_t,uintptr_t,uintptr_t,uintptr_t,uintptr_t,uint32_t);
extern "C" Logger* malloc_logger;
struct Block {uintptr_t address=0;size_t bytes=0;};
inline std::array<Block,65536> blocks;
inline std::atomic_flag busy=ATOMIC_FLAG_INIT;
inline bool overflow=false,active=false;
inline size_t live=0,peak=0,occupied=0;
inline constexpr size_t mask=65535;
inline size_t bucket(uintptr_t ptr){return ((ptr>>4)*11400714819323198485ull)>>48;}
inline void remove(uintptr_t ptr)
{
    if(!ptr)return;
    size_t hole=bucket(ptr);
    while(blocks[hole].address&&blocks[hole].address!=ptr)hole=(hole+1)&mask;
    if(!blocks[hole].address)return;
    live-=blocks[hole].bytes;--occupied;
    for(size_t pos=(hole+1)&mask;blocks[pos].address;pos=(pos+1)&mask)
    {
        const size_t home=bucket(blocks[pos].address);
        if(((pos-home)&mask)>((hole-home)&mask)){blocks[hole]=blocks[pos];hole=pos;}
    }
    blocks[hole]={};
}
inline void record(uint32_t type,uintptr_t,uintptr_t arg2,uintptr_t arg3,uintptr_t result,uint32_t)
{
    if(!(type&8)||!(type&6))return;
    while(busy.test_and_set(std::memory_order_acquire)){}
    if(active&&!overflow)
    {
        const bool allocate=(type&2)!=0,deallocate=(type&4)!=0;
        if(deallocate&&(!allocate||result||arg3==0))remove(arg2);
        if(allocate&&result&&occupied>=mask)overflow=true;
        if(allocate&&result&&!overflow)
        {
            size_t pos=bucket(result),start=pos;
            while(blocks[pos].address){pos=(pos+1)&mask;if(pos==start){overflow=true;break;}}
            if(!overflow)
            {
                const auto bytes=malloc_size(reinterpret_cast<const void*>(result));
                blocks[pos]={result,bytes};++occupied;live+=bytes;peak=std::max(peak,live);
            }
        }
    }
    busy.clear(std::memory_order_release);
}
class Scope
{
public:
    static constexpr bool supported=true;
    Scope()
    {
        if(malloc_logger||active)throw std::runtime_error("Qualification malloc callback already in use");
        blocks.fill({});live=peak=occupied=0;overflow=false;active=true;malloc_logger=record;
    }
    Scope(const Scope&)=delete;
    ~Scope(){if(!ended)stop();}
    size_t stop()
    {
        if(ended)return maximum;
        while(busy.test_and_set(std::memory_order_acquire)){}
        active=false;maximum=peak;ended=true;
        busy.clear(std::memory_order_release);malloc_logger=nullptr;return maximum;
    }
    bool valid()const{return !overflow&&ended&&maximum>0;}
private:
    size_t maximum=0;bool ended=false;
};
}
#else
namespace event_allocation
{
struct Scope
{
    static constexpr bool supported=false;
    size_t stop()const{return 0;}
    bool valid()const{return false;}
};
}
#endif

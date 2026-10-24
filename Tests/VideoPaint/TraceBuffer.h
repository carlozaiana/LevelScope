#pragma once
#include <atomic>
#include <cstdint>
#include <memory>

namespace levelscope::paintdiagnostic
{
enum class Kind : uint32_t { timerEnter=1, timerExit, selection, repaintRequest,
    paintEnter, paintExit, mainVblank, displayTick, windowState, stopping };
struct Event
{
    uint64_t hostNs=0;
    Kind kind=Kind::timerEnter;
    uint32_t flags=0;
    int64_t timeNumerator=0, timeDenominator=1;
    uint64_t generation=0, epoch=0;
    int64_t desired=-1, picture=-1, queueFirst=-1, queueLast=-1;
    uint64_t queueCount=0, display=0, auxiliary=0;
};
// One producer per buffer. Read only after that producer has stopped. Never
// overwrite earlier events: overflow invalidates the diagnostic, visibly.
template<size_t Capacity> class Buffer
{
public:
    Buffer() : events(new Event[Capacity]) {}
    void push(const Event& e) noexcept
    {
        if (written == Capacity) { ++dropped; return; }
        events[written++]=e;
    }
    size_t size() const noexcept { return written; }
    uint64_t lost() const noexcept { return dropped; }
    const Event& at(size_t i) const noexcept { return events[i]; }
    static constexpr size_t bytes=sizeof(Event)*Capacity;
private:
    std::unique_ptr<Event[]> events;
    size_t written=0;
    uint64_t dropped=0;
};
}

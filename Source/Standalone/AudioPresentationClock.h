#pragma once
#include <array>
#include <atomic>
#include <cstdint>
#include <algorithm>

namespace levelscope::standalone
{
// A bounded history of submitted source spans mapped to predicted DAC times.
// Writer: audio callback. Reader: message/test thread. Atomic fields avoid the
// C++ data race in a conventional non-atomic seqlock. No allocation or RT locks.
class AudioPresentationClock
{
public:
    struct Span { uint64_t generation = 0, hostNs = 0; int64_t first = 0; int frames = 0; };
    static constexpr size_t capacity = 4096;
    void reset() noexcept { written.store (0); } // callback must be quiescent
    void submit (Span span) noexcept
    {
        const auto index = written.load (std::memory_order_relaxed);
        auto& slot = slots[index % capacity];
        slot.serial.store (index * 2 + 1, std::memory_order_seq_cst);
        slot.generation.store (span.generation, std::memory_order_seq_cst);
        slot.host.store (span.hostNs, std::memory_order_seq_cst);
        slot.first.store (span.first, std::memory_order_seq_cst);
        slot.frames.store (span.frames, std::memory_order_seq_cst);
        slot.serial.store (index * 2 + 2, std::memory_order_seq_cst);
        written.store (index + 1, std::memory_order_release);
    }
    int64_t position (uint64_t nowNs, int64_t rate, uint64_t generation, int64_t fallback) const noexcept
    {
        const auto count = written.load (std::memory_order_acquire);
        const auto begin = count > capacity ? count - capacity : 0;
        for (auto i = count; i > begin; --i)
        {
            const auto& slot = slots[(i - 1) % capacity];
            const auto serial = slot.serial.load (std::memory_order_seq_cst);
            if (serial != (i - 1) * 2 + 2) continue;
            Span s { slot.generation.load (std::memory_order_seq_cst), slot.host.load (std::memory_order_seq_cst),
                     slot.first.load (std::memory_order_seq_cst), slot.frames.load (std::memory_order_seq_cst) };
            if (slot.serial.load (std::memory_order_seq_cst) != serial || s.generation != generation) continue;
            if (nowNs < s.hostNs) { fallback = s.first; continue; }
            const auto elapsed = static_cast<__int128> (nowNs - s.hostNs) * rate / 1000000000;
            return s.first + static_cast<int64_t> (std::min<__int128> (s.frames, elapsed));
        }
        return fallback;
    }
private:
    static_assert (std::atomic<int64_t>::is_always_lock_free && std::atomic<uint64_t>::is_always_lock_free
                   && std::atomic<int>::is_always_lock_free, "Audio clock requires lock-free atomics");
    struct Slot
    {
        std::atomic<uint64_t> serial { 0 }, generation { 0 }, host { 0 };
        std::atomic<int64_t> first { 0 };
        std::atomic<int> frames { 0 };
    };
    std::array<Slot, capacity> slots;
    std::atomic<uint64_t> written { 0 };
};
}

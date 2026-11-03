#pragma once
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

namespace levelscope::waveform
{
constexpr size_t combinedBudgetBytes = 128u * 1024u * 1024u;
// Covers the two bounded display buffers during replacement, channel labels,
// mailbox/control metadata and allocator overhead outside overview slabs.
constexpr size_t displayReserveBytes = 8u * 1024u * 1024u;
constexpr size_t maximumOverviewBytes = (combinedBudgetBytes-displayReserveBytes)/2;
constexpr int channelCapacity = 12, maximumColumns = 4096;
struct Spec { int64_t frames=0; int sampleRate=0, channels=0; };
struct Extrema
{
    float minimum=0, maximum=0;
    int64_t minimumFrame=0, maximumFrame=0;
};
Extrema merge(const Extrema&,const Extrema&) noexcept;
struct Column
{
    int64_t begin=0,end=0,coveredEnd=0;
    std::array<Extrema,channelCapacity> channels{};
    bool available() const noexcept {return coveredEnd>begin;}
};
class Budget
{
public:
    explicit Budget(size_t limit=combinedBudgetBytes);
    size_t used() const noexcept {return bytes.load();}
    size_t remaining() const noexcept {return limitBytes-used();}
    size_t limit() const noexcept {return limitBytes;}
private:
    friend class Overview;
    const size_t limitBytes;
    std::atomic<size_t> bytes{displayReserveBytes};
    bool reserve(size_t) noexcept;
    void release(size_t n) noexcept {bytes.fetch_sub(n);}
};
// One worker appends sequentially. Readers acquire a sealed immutable prefix.
// No sample storage, growing queue, per-publication copy, or audio callback use.
class Overview
{
public:
    static std::shared_ptr<Overview> create(Spec,std::shared_ptr<Budget>);
    ~Overview();
    Overview(const Overview&)=delete;
    Overview& operator=(const Overview&)=delete;
    const Spec spec;
    const int64_t binFrames;
    const size_t storageBytes;
    bool append(const float* const*,int channels,int count,int64_t firstFrame) noexcept;
    bool finish() noexcept;
    uint64_t publication() const noexcept {return published.load(std::memory_order_acquire);}
    int64_t coveredFrames() const noexcept;
    int64_t coveredFrames(uint64_t snapshot) const noexcept;
    bool complete() const noexcept {return (publication()&1)!=0;}
    // Returned buffers are bounded to 4096 columns. The UI retains at most two
    // during replacement, covered by the shared display reservation above.
    std::vector<Column> reduce(int64_t begin,int64_t end,int columns,uint64_t snapshot=UINT64_MAX) const;
    std::optional<Column> inspect(int64_t frame) const;
private:
    struct Level {size_t bins=0,offset=0;};
    struct Layout {int64_t width=64;size_t nodes=0,bytes=0;int levels=0;std::array<Level,64> level{};};
    Overview(Spec,std::shared_ptr<Budget>,const Layout&);
    std::shared_ptr<Budget> budget;
    const Layout layout;
    std::unique_ptr<Extrema[]> bins;
    std::array<Extrema,channelCapacity> pending{};
    int64_t consumed=0;
    size_t sealed=0;
    bool failed=false,finished=false;
    std::atomic<uint64_t> published{0}; // sealed finest bins << 1, plus EOF flag
    std::chrono::steady_clock::time_point lastPublication{};
    void seal() noexcept;
    void publish(bool final) noexcept;
    Extrema query(size_t begin,size_t end,int channel) const noexcept;
    int64_t frameAt(size_t bin) const noexcept;
};
}

#pragma once
#include "StandaloneWaveform.h"
#include "Core/Proposal/Model/ProgramAnalysis.h"
#include <chrono>

namespace levelscope::standalone
{
// One analysis worker writes each window once, readers acquire a sealed prefix.
// Fixed arrays, no growing queue or complete-history copies on publication.
class LoudnessPreview final : public measurement::WindowObserver
{
public:
    static constexpr size_t maximumBytes=4u*1024u*1024u;
    LoudnessPreview(int64_t frames,int rate);
    const int64_t frames;
    const int sampleRate;
    void window(bool shortTerm,const measurement::Window&) noexcept override;
    void publish(int64_t coveredFrames,bool force=false) noexcept;
    int64_t coverage() const noexcept {return published.load(std::memory_order_acquire);}
    proposal::TimelinePlot reduce(double begin,double end,int columns,int64_t coverage) const;
    std::optional<measurement::Window> inspect(bool shortTerm,double frame,int64_t coverage) const;
    size_t storageBytes() const {return (capacity[0]+capacity[1])*sizeof(measurement::Window);}
private:
    std::array<std::unique_ptr<measurement::Window[]>,2> windows;
    std::array<size_t,2> capacity{},written{};
    std::atomic<int64_t> published{0};
    std::chrono::steady_clock::time_point lastPublication{};
    bool failed=false;
    size_t count(bool shortTerm,int64_t covered) const;
};
class LoudnessProgress
{
public:
    explicit LoudnessProgress(WaveformBinding b):binding(std::move(b)){}
    const WaveformBinding binding;
    std::shared_ptr<LoudnessPreview> begin(int64_t frames,int rate) noexcept;
    std::shared_ptr<const LoudnessPreview> view(const WaveformBinding&) const;
    void abandon();
private:
    mutable std::mutex mutex;
    std::shared_ptr<LoudnessPreview> data;
    bool abandoned=false;
};
}

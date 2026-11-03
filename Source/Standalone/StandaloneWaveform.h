#pragma once
#include "StandaloneServices.h"
#include "Core/Waveform/WaveformOverview.h"
#include <mutex>

namespace levelscope::standalone
{
// Worker-only publication boundary. Allocation failures affect the replaceable
// waveform, never measurement arithmetic or acceptance. This is not an RT sink.
class WaveformProgress
{
public:
    enum class State { waiting, collecting, complete, unavailable, abandoned };
    explicit WaveformProgress(std::shared_ptr<waveform::Budget>);
    std::shared_ptr<waveform::Overview> begin(waveform::Spec) noexcept;
    std::shared_ptr<const waveform::Overview> overview() const;
    void finish(const std::shared_ptr<waveform::Overview>&) noexcept;
    void unavailable() noexcept;
    void abandon() noexcept;
    State state() const noexcept {return phase.load();}
private:
    std::shared_ptr<waveform::Budget> budget;
    mutable std::mutex mutex;
    std::shared_ptr<waveform::Overview> data;
    std::atomic<State> phase{State::waiting};
};
struct WaveformBinding
{
    juce::String runtime,identity;
    JobTicket ticket;
    MeasurementDomain domain=MeasurementDomain::source;
    static WaveformBinding from(const ServiceContext&,MeasurementDomain);
    bool sameAudio(const WaveformBinding&) const;
    bool sameJob(const WaveformBinding&) const;
};
struct WaveformView
{
    std::shared_ptr<const waveform::Overview> data;
    bool accepted=false;
    WaveformProgress::State state=WaveformProgress::State::waiting;
};
// Runtime-only Source/Current cache. A single pool accounts for live, obsolete
// worker-owned and shared overviews alike, including temporary display storage.
class WaveformCache
{
public:
    std::shared_ptr<WaveformProgress> begin(const WaveformBinding&);
    bool accept(const WaveformBinding& completed,const WaveformBinding& expected);
    void reject(const WaveformBinding&);
    void synchronize(const WaveformBinding& source,const WaveformBinding& current,bool sourceValid,bool currentInitialized);
    bool shareInitialCurrent(const AnalysisSnapshot& source,const AnalysisSnapshot& current);
    WaveformView view(const WaveformBinding&) const;
    size_t retainedBytes() const {return budget->used();}
private:
    struct Slot {WaveformBinding binding;std::shared_ptr<const waveform::Overview> data;};
    std::shared_ptr<waveform::Budget> budget=std::make_shared<waveform::Budget>();
    std::array<Slot,2> slots{};
    WaveformBinding active;
    std::shared_ptr<WaveformProgress> progress;
};
}

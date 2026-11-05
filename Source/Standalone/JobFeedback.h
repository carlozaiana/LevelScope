#pragma once
#include <atomic>
#include <cmath>
#include <cstdint>
#include <juce_core/juce_core.h>

namespace levelscope::standalone
{
// Worker metadata only: never used by measurement decisions or audio callbacks.
// A coherent, bounded snapshot; no retained event history or per-block allocation.
enum class JobPhase : uint32_t { starting, metadata, samples, checksum, sourceHash, copy, copyHash, measuring, finalising };
struct JobPhaseSnapshot { JobPhase phase = JobPhase::starting; double fraction = -1; };
class JobFeedback
{
public:
    virtual ~JobFeedback() = default;
    virtual void publish(JobPhase phase, double fraction)
    {
        const uint32_t value = std::isfinite(fraction) && fraction >= 0
            ? uint32_t(juce::jlimit(0.0, 1.0, fraction) * 1000000) : UINT32_MAX;
        state.store((uint64_t(phase) << 32) | value, std::memory_order_relaxed);
    }
    JobPhaseSnapshot snapshot() const
    {
        const auto value = state.load(std::memory_order_relaxed);
        const auto progress = uint32_t(value);
        return {JobPhase(value >> 32), progress == UINT32_MAX ? -1.0 : progress / 1000000.0};
    }
private:
    std::atomic<uint64_t> state { UINT32_MAX };
};
struct JobFeedbackText { juce::String title, detail; double fraction = -1; };
inline JobFeedbackText describeJobFeedback(JobPhaseSnapshot stage, bool analysis, bool current, bool cancelling)
{
    if(cancelling) return {"Cancelling - waiting for the worker to stop", "Incomplete results will not become completed measurements.", -1};
    JobFeedbackText text;
    const juce::String domain = current ? "Current analysis" : "Original analysis";
    const bool measuring = stage.phase == JobPhase::measuring || stage.phase == JobPhase::finalising;
    text.title = analysis ? domain + (measuring ? ": measurement" : ": preparation / verification") : "Import: source validation";
    switch(stage.phase)
    {
        case JobPhase::starting: text.detail = "Starting; no new completed measurement yet."; break;
        case JobPhase::metadata: text.detail = "Checking file format and channel metadata"; break;
        case JobPhase::samples: text.detail = "Checking decoded samples"; break;
        case JobPhase::checksum: text.detail = "Verifying FLAC audio checksum"; break;
        case JobPhase::sourceHash: text.detail = "Verifying file content identity"; break;
        case JobPhase::copy: text.detail = "Creating protected source copy"; break;
        case JobPhase::copyHash: text.detail = "Verifying source copy identity"; break;
        case JobPhase::measuring: text.detail = "Measuring programme; waveform and loudness previews are incomplete"; break;
        case JobPhase::finalising: text.detail = "Finalising measurement; waiting for completed results"; break;
    }
    text.fraction = stage.fraction;
    if(stage.fraction >= 0)
        text.detail += " - " + juce::String(juce::roundToInt(stage.fraction * 100)) + "% of " + (measuring ? "programme" : "current file / pass");
    else if(stage.phase != JobPhase::starting) text.detail += "...";
    return text;
}
}

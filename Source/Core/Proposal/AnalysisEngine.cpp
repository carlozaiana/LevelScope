#include "AnalysisEngine.h"
#include "AttributionEngine.h"
#include "Core/Measurement/GatePopulation.h"
#include <juce_cryptography/juce_cryptography.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace levelscope::proposal
{
namespace
{
void require (bool valid, const char* reason = "Inconsistent or unsupported measured diagnostics input") { if (! valid) throw std::invalid_argument (reason); }
void checkMetric (const measurement::Metric& m)
{
    require ((m.status == measurement::Status::available) == m.value.has_value());
    if (m.value) require (std::isfinite (*m.value));
}
}
size_t ProgramAnalysis::storageBytes() const
{
    return sizeof (*this) + (momentary.capacity() + shortTerm.capacity() + lraWindows.capacity()) * sizeof (DiagnosticWindow)
        + (integratedRank.capacity() + excludedIntegrated.capacity() + lraOrder.capacity() + lowerTail.capacity() + upperTail.capacity()) * sizeof (size_t)
        + histogram.capacity() * sizeof (DistributionBin) + peaks.capacity() * sizeof (PeakLocation);
}

std::shared_ptr<const ProgramAnalysis> analyzeProgramme (const measurement::Result& r, const juce::String& snapshot,
                                                        const std::function<bool()>& cancelled)
{
    if (cancelled()) return {};
    require (snapshot.isNotEmpty() && (r.sampleRate == 44100 || r.sampleRate == 48000 || r.sampleRate == 96000 || r.sampleRate == 192000)
             && r.frames > 0 && r.frames <= static_cast<int64_t> (r.sampleRate) * 24 * 3600);
    for (const auto* metric : { &r.integrated, &r.lra, &r.momentaryMax, &r.shortTermMax, &r.samplePeak, &r.truePeak }) checkMetric (*metric);
    for (const auto& value : { r.integratedRelativeGate, r.lraRelativeGate, r.p10, r.p95 }) if (value) require (std::isfinite (*value));
    auto a = std::make_shared<ProgramAnalysis>(); a->snapshotIdentity = snapshot;
    const auto binding = juce::String (diagnosticsVersion) + ";" + snapshot;
    a->identity = juce::SHA256 (binding.toRawUTF8(), binding.getNumBytesAsUTF8()).toHexString();
    a->sampleRate = r.sampleRate; a->frames = r.frames;
    a->integrated = r.integrated; a->lra = r.lra; a->momentaryMax = r.momentaryMax; a->shortTermMax = r.shortTermMax;
    a->samplePeak = r.samplePeak; a->truePeak = r.truePeak;
    a->integratedGate = r.integratedRelativeGate; a->lraGate = r.lraRelativeGate;
    a->p10 = r.p10; a->p95 = r.p95; a->lraStable = r.lraStable;
    const auto ig = measurement::relativeGateEnergy (r.momentary, false, cancelled);
    const auto lg = measurement::relativeGateEnergy (r.lraWindows, true, cancelled);
    if (cancelled()) return {};
    require (r.integratedRelativeGate == (ig ? std::optional<double> (measurement::OfflineMeasurement::loudness (*ig)) : std::nullopt), "Integrated gate does not match measured population");
    const auto eligibleLra = r.frames >= static_cast<int64_t> (r.sampleRate) * 3 && lg.has_value();
    require (r.lraRelativeGate == (eligibleLra ? std::optional<double> (measurement::OfflineMeasurement::loudness (*lg)) : std::nullopt), "LRA gate does not match measured population");
    const auto fill = [&] (const std::vector<measurement::Window>& in, std::vector<DiagnosticWindow>& out,
                          int64_t length, bool lra, bool ungated)
    {
        const auto step = r.sampleRate / 10;
        const auto limit = r.frames + (lra ? r.sampleRate * 3 / 2 : 0);
        const auto expectedCount = limit < length ? 0 : (limit - length) / step + 1;
        require (in.size() == static_cast<size_t> (expectedCount));
        out.reserve (in.size());
        for (size_t n = 0; n < in.size(); ++n)
        {
            if ((n % 256) == 0 && cancelled()) return false;
            const auto& w = in[n];
            require (w.endFrame == length + static_cast<int64_t> (n) * step && std::isfinite (w.energy) && w.energy >= 0);
            DiagnosticWindow d; d.begin = w.endFrame - length; d.end = w.endFrame; d.energy = w.energy;
            if (w.energy > 0) d.lufs = measurement::OfflineMeasurement::loudness (w.energy);
            d.padding = w.endFrame > r.frames;
            d.included = ! ungated && (lra ? eligibleLra && measurement::includedInPopulation (w.energy, *lg, true)
                                                       : ig && measurement::includedInPopulation (w.energy, *ig, false));
            out.push_back (d);
        }
        return true;
    };
    if (! fill (r.momentary, a->momentary, r.sampleRate * 4 / 10, false, false)
        || ! fill (r.shortTerm, a->shortTerm, r.sampleRate * 3, false, true)
        || ! fill (r.lraWindows, a->lraWindows, r.sampleRate * 3, true, false)) return {};
    if (! buildAttribution (*a, cancelled)) return {};
    require (a->integratedRank.size() == r.integratedIncluded && a->lraOrder.size() == r.lraIncluded, "Included window counts differ from measured result");
    double sum = 0;
    for (const auto& w : a->momentary) if (w.included) sum += w.energy;
    if (! a->integratedRank.empty()) require (r.integrated.value == measurement::OfflineMeasurement::loudness (sum / static_cast<double> (a->integratedRank.size())), "Integrated value differs from measured population");
    else require (! r.integrated.value);
    if (! a->lraOrder.empty())
    {
        const auto low = a->lraWindows[a->lraOrder[measurement::percentileIndex (a->lraOrder.size(), .10)]].lufs;
        const auto high = a->lraWindows[a->lraOrder[measurement::percentileIndex (a->lraOrder.size(), .95)]].lufs;
        require (low == r.p10 && high == r.p95 && r.lra.value == *high - *low, "LRA percentiles differ from measured population");
    }
    else require (! r.p10 && ! r.p95 && ! r.lra.value);
    const auto channels = r.channelSamplePeaks.size();
    require (channels > 0 && channels <= measurement::channelCapacity && r.channelTruePeaks.size() == channels
             && r.samplePeakFrames.size() == channels && r.truePeakEighthFrames.size() == channels);
    for (size_t n = 0; n < channels; ++n)
    {
        require (std::isfinite (r.channelSamplePeaks[n]) && r.channelSamplePeaks[n] >= 0
                 && std::isfinite (r.channelTruePeaks[n]) && r.channelTruePeaks[n] >= r.channelSamplePeaks[n]);
        PeakLocation peak; peak.channel = static_cast<int> (n); peak.sampleFrame = r.samplePeakFrames[n]; peak.truePeakEighthFrame = r.truePeakEighthFrames[n];
        if (r.channelSamplePeaks[n] > 0) { require (peak.sampleFrame >= 0 && peak.sampleFrame < r.frames); peak.sampleDbfs = 20 * std::log10 (r.channelSamplePeaks[n]); }
        if (r.channelTruePeaks[n] > 0)
        {
            require (peak.truePeakEighthFrame >= -32 * 8 && peak.truePeakEighthFrame <= (r.frames + 64) * 8);
            peak.trueDbtp = 20 * std::log10 (r.channelTruePeaks[n]);
        }
        a->peaks.push_back (peak);
    }
    return cancelled() ? nullptr : a;
}

TimelinePlot prepareTimeline (const ProgramAnalysis& a, double begin, double end, int columns)
{
    require (std::isfinite (begin) && std::isfinite (end) && end > begin && columns > 0);
    TimelinePlot plot; plot.begin = std::clamp (begin, 0., static_cast<double> (a.frames));
    plot.end = std::clamp (end, plot.begin, static_cast<double> (a.frames));
    require (plot.end > plot.begin);
    columns = std::min (4096, columns);
    plot.momentary.resize (static_cast<size_t> (columns)); plot.shortTerm.resize (static_cast<size_t> (columns));
    auto fill = [&] (const auto& windows, auto& bins, bool ungated)
    {
        auto it = std::lower_bound (windows.begin(), windows.end(), plot.begin,
                                    [] (const auto& w, double f) { return w.end < f; });
        for (; it != windows.end() && it->end <= plot.end; ++it)
        {
            ++plot.fullResolutionPoints;
            if (! it->lufs) continue;
            auto& b = bins[static_cast<size_t> (std::min (columns - 1, static_cast<int> ((it->end - plot.begin) * columns / (plot.end - plot.begin))))];
            (ungated ? b.ungated : it->included ? b.included : b.excluded).add (*it->lufs);
        }
    };
    fill (a.momentary, plot.momentary, false); fill (a.shortTerm, plot.shortTerm, true);
    return plot;
}

std::optional<size_t> inspectWindow (const std::vector<DiagnosticWindow>& windows, double frame)
{
    if (windows.empty() || ! std::isfinite (frame)) return {};
    auto it = std::lower_bound (windows.begin(), windows.end(), frame, [] (const auto& w, double f) { return w.end < f; });
    if (it == windows.end()) return windows.size() - 1;
    if (it != windows.begin() && frame - (it - 1)->end <= it->end - frame) --it;
    return static_cast<size_t> (it - windows.begin());
}
}

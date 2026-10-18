#pragma once
#include "Core/Measurement/OfflineMeasurement.h"
#include <juce_core/juce_core.h>
#include <memory>

namespace levelscope::proposal
{
inline constexpr const char* diagnosticsVersion = "measured-diagnostics-m4.1";
struct DiagnosticWindow
{
    int64_t begin = 0, end = 0; // half-open audio-source frames; plot at end
    double energy = 0, energyShare = 0;
    std::optional<double> lufs;
    bool included = false, padding = false;
};
struct PeakLocation
{
    int channel = 0; // confirmed file order, never device order
    int64_t sampleFrame = 0, truePeakEighthFrame = 0;
    std::optional<double> sampleDbfs, trueDbtp;
};
struct DistributionBin
{
    double lower = 0; // [lower, lower+1) LUFS; descriptive display bins only
    size_t included = 0, excluded = 0, cumulativeIncluded = 0;
};
struct ProgramAnalysis
{
    juce::String identity, snapshotIdentity;
    int sampleRate = 0;
    int64_t frames = 0;
    measurement::Metric integrated, lra, momentaryMax, shortTermMax, samplePeak, truePeak;
    std::optional<double> integratedGate, lraGate, p10, p95;
    bool lraStable = false;
    size_t missingLraWindows = 0;
    std::vector<DiagnosticWindow> momentary, shortTerm, lraWindows;
    std::vector<size_t> integratedRank, excludedIntegrated, lraOrder, lowerTail, upperTail;
    std::vector<DistributionBin> histogram;
    std::vector<PeakLocation> peaks;
    size_t storageBytes() const;
};
struct ValueExtent
{
    size_t count = 0;
    double minimum = 0, maximum = 0;
    void add (double v) { if (count++ == 0) minimum = maximum = v; else { minimum = std::min (minimum, v); maximum = std::max (maximum, v); } }
};
struct PlotBucket
{
    // Separate ranges preserve mixed gate membership in a single display column.
    ValueExtent included, excluded, ungated;
};
struct TimelinePlot
{
    double begin = 0, end = 0; // source frames, clamped to [0, programme frames]
    std::vector<PlotBucket> momentary, shortTerm;
    size_t fullResolutionPoints = 0;
};
}

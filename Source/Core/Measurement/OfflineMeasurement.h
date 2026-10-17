#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <vector>

namespace levelscope::measurement
{
// Independent of plugin detectors/history, GUI and device state. Worker-thread only.
inline constexpr const char* version = "bs1770-5-ebu2023-m1.1";
enum class Role { left, right, centre, surroundLeft, surroundRight, lfe };
enum class Status { available, silence, insufficientDuration, belowGate };
struct Metric
{
    std::optional<double> value;
    Status status = Status::insufficientDuration;
};
struct Window
{
    int64_t endFrame = 0; // half-open, source sample clock; LRA drain may exceed EOF
    double energy = 0;
};
struct Result
{
    int sampleRate = 0;
    int64_t frames = 0;
    Metric integrated, lra, momentaryMax, shortTermMax, samplePeak, truePeak;
    std::optional<double> integratedRelativeGate, lraRelativeGate, p10, p95;
    size_t integratedIncluded = 0, lraIncluded = 0;
    bool lraStable = false;
    std::vector<Window> momentary, shortTerm, lraWindows;
    std::vector<double> channelSamplePeaks, channelTruePeaks; // linear, silence = 0
    std::vector<int64_t> samplePeakFrames, truePeakEighthFrames;
};

// Invalid configuration/samples throw std::invalid_argument; accumulation overflow
// throws std::overflow_error. finish() is idempotent; push after finish is invalid.
class OfflineMeasurement
{
public:
    OfflineMeasurement (int sampleRate, const std::vector<Role>& fileRoles);
    void push (const float* const* channels, int count);
    const Result& finish();
    static bool supports (int sampleRate, const std::vector<Role>& roles);
    static double loudness (double energy);
private:
    struct Filter
    {
        double b0 = 1, b1 = 0, b2 = 0, a1 = 0, a2 = 0, z1 = 0, z2 = 0;
        double process (double x);
    };
    struct Channel
    {
        Filter shelf, highpass;
        double weight = 1;
        std::array<double, 128> history {};
        int position = 0;
    };
    void addEnergy (double energy, bool drain);
    void addPeak (size_t channel, double sample, int64_t frame, bool drain);
    Result result;
    std::vector<Channel> channels;
    std::array<std::array<double, 64>, 8> interpolation {};
    std::vector<double> energyRing;
    int64_t clock = 0;
    int mLength, sLength, step;
    double mSum = 0, sSum = 0, mMaximum = 0, sMaximum = 0;
    bool finished = false;
};
}

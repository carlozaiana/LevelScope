#include "OfflineMeasurement.h"
#include "GatePopulation.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace levelscope::measurement
{
namespace
{
constexpr double pi = 3.1415926535897932384626433832795;
Metric fromEnergy (double energy, bool enough)
{
    if (! enough) return {};
    if (energy <= 0) return { {}, Status::silence };
    return { OfflineMeasurement::loudness (energy), Status::available };
}
Metric fromPeak (double peak)
{
    return peak > 0 ? Metric { 20 * std::log10 (peak), Status::available }
                    : Metric { {}, Status::silence };
}
}

double OfflineMeasurement::loudness (double energy) { return -0.691 + 10 * std::log10 (energy); }
double OfflineMeasurement::Filter::process (double x)
{
    const double y = b0 * x + z1;
    z1 = b1 * x - a1 * y + z2;
    z2 = b2 * x - a2 * y;
    return y;
}

bool OfflineMeasurement::supports (int rate, const std::vector<Role>& roles)
{
    if (rate != 44100 && rate != 48000 && rate != 96000 && rate != 192000) return false;
    auto sorted = roles;
    std::sort (sorted.begin(), sorted.end());
    if (std::adjacent_find (sorted.begin(), sorted.end()) != sorted.end()) return false;
    if (roles.size() == 1) return roles[0] == Role::left || roles[0] == Role::right || roles[0] == Role::centre;
    if (roles.size() == 2) return sorted == std::vector<Role> { Role::left, Role::right };
    return sorted == std::vector<Role> { Role::left, Role::right, Role::centre, Role::surroundLeft, Role::surroundRight }
        || sorted == std::vector<Role> { Role::left, Role::right, Role::centre, Role::surroundLeft, Role::surroundRight, Role::lfe }
        || sorted == std::vector<Role> { Role::left, Role::right, Role::centre, Role::surroundLeft, Role::surroundRight, Role::lfe, Role::rearLeft, Role::rearRight };
}

OfflineMeasurement::OfflineMeasurement (int rate, const std::vector<Role>& roles, PeakObserver* observer)
    : mLength (0), sLength (0), step (0), peakObserver (observer)
{
    if (! supports (rate, roles)) throw std::invalid_argument ("Unqualified measurement sample rate or layout");
    mLength = rate * 4 / 10; sLength = rate * 3; step = rate / 10;
    result.sampleRate = rate;
    energyRing.resize (static_cast<size_t> (sLength));
    channels.resize (roles.size());
    result.channelSamplePeaks.resize (roles.size()); result.channelTruePeaks.resize (roles.size());
    result.samplePeakFrames.resize (roles.size(), -1); result.truePeakEighthFrames.resize (roles.size(), -1);
    for (size_t i = 0; i < roles.size(); ++i)
    {
        auto& c = channels[i];
        c.weight = roles[i] == Role::lfe ? 0 : (roles[i] == Role::surroundLeft || roles[i] == Role::surroundRight ? 1.41 : 1);
        if (rate == 48000)
        {
            c.shelf = { 1.53512485958697, -2.69169618940638, 1.19839281085285, -1.69065929318241, 0.73248077421585 };
            c.highpass = { 1, -2, 1, -1.99004745483398, 0.99007225036621 };
        }
        else
        {
            // Bilinear design matching the published 48 kHz response (De Man form).
            double k = std::tan (pi * 1681.974450955533 / rate), q = 0.7071752369554196;
            const double vh = std::pow (10.0, 3.999843853973347 / 20), vb = std::pow (vh, 0.4996667741545416);
            double a0 = 1 + k / q + k * k;
            c.shelf = { (vh + vb*k/q + k*k)/a0, 2*(k*k-vh)/a0,
                        (vh-vb*k/q+k*k)/a0, 2*(k*k-1)/a0, (1-k/q+k*k)/a0 };
            k = std::tan (pi * 38.13547087602444 / rate); q = 0.5003270373238773;
            a0 = 1 + k/q + k*k;
            c.highpass = { 1, -2, 1, 2*(k*k-1)/a0, (1-k/q+k*k)/a0 };
        }
    }
    // 8x, 64-tap/phase Blackman-windowed sinc, unity DC gain per phase.
    // Fractional phases use a common continuous window; delay = 31 + phase/8.
    for (int phase = 0; phase < 8; ++phase)
    {
        double sum = 0;
        for (int tap = 0; tap < 64; ++tap)
        {
            const double x = tap - 31.0 - phase / 8.0;
            const double sinc = std::abs (x) < 1e-12 ? 1 : std::sin (pi*x)/(pi*x);
            const double window = std::abs (x) <= 32 ? 0.42 + 0.5*std::cos (pi*x/32) + 0.08*std::cos (2*pi*x/32) : 0;
            interpolation[phase][tap] = sinc * window;
            sum += interpolation[phase][tap];
        }
        for (auto& value : interpolation[phase]) value /= sum;
    }
}

void OfflineMeasurement::addPeak (size_t index, double sample, int64_t frame, bool drain)
{
    auto& c = channels[index];
    const double magnitude = std::abs (sample);
    if (peakObserver)
    {
        if (!drain) peakObserver->sample(index, magnitude, frame);
        peakObserver->truePeak(index, magnitude, frame * 8, drain);
    }
    if (! drain && magnitude > result.channelSamplePeaks[index])
    {
        result.channelSamplePeaks[index] = magnitude;
        result.samplePeakFrames[index] = frame;
    }
    if (magnitude > result.channelTruePeaks[index])
    {
        result.channelTruePeaks[index] = magnitude;
        result.truePeakEighthFrames[index] = frame * 8;
    }
    c.position = (c.position + 63) % 64;
    c.history[c.position] = c.history[c.position + 64] = sample;
    for (int phase = 0; phase < 8; ++phase)
    {
        double interpolated = 0;
        for (int tap = 0; tap < 64; ++tap) interpolated += c.history[c.position + tap] * interpolation[phase][tap];
        if (peakObserver) peakObserver->truePeak(index, std::abs(interpolated), (frame - 31) * 8 - phase, drain);
        if (std::abs (interpolated) > result.channelTruePeaks[index])
        {
            result.channelTruePeaks[index] = std::abs (interpolated);
            result.truePeakEighthFrames[index] = (frame - 31) * 8 - phase;
        }
    }
}

void OfflineMeasurement::addEnergy (double energy, bool drain)
{
    const auto position = static_cast<size_t> (clock % sLength);
    const double oldM = clock >= mLength ? energyRing[static_cast<size_t> ((clock - mLength) % sLength)] : 0;
    mSum += energy - oldM;
    sSum += energy - energyRing[position];
    energyRing[position] = energy;
    ++clock;
    if (! std::isfinite (sSum) || ! std::isfinite (mSum)) throw std::overflow_error ("Measurement energy overflow");
    const double mEnergy = std::max (0.0, mSum / mLength), sEnergy = std::max (0.0, sSum / sLength);
    if (! drain)
    {
        if (clock >= mLength) mMaximum = std::max (mMaximum, mEnergy);
        if (clock >= sLength) sMaximum = std::max (sMaximum, sEnergy);
    }
    if (clock % step == 0)
    {
        if (! drain && clock >= mLength) result.momentary.push_back ({ clock, mEnergy });
        if (clock >= sLength)
        {
            result.lraWindows.push_back ({ clock, sEnergy });
            if (! drain) result.shortTerm.push_back ({ clock, sEnergy });
        }
    }
}

void OfflineMeasurement::push (const float* const* input, int count)
{
    if (finished || count < 0 || input == nullptr) throw std::invalid_argument ("Invalid measurement push");
    for (size_t c = 0; c < channels.size(); ++c)
        if (input[c] == nullptr) throw std::invalid_argument ("Null channel");
    for (int frame = 0; frame < count; ++frame)
    {
        double energy = 0;
        for (size_t i = 0; i < channels.size(); ++i)
        {
            auto& c = channels[i];
            const double x = input[i][frame];
            if (! std::isfinite (x)) throw std::invalid_argument ("Non-finite audio");
            addPeak (i, x, result.frames, false);
            const double y = c.highpass.process (c.shelf.process (x));
            energy += c.weight * y * y;
        }
        addEnergy (energy, false);
        ++result.frames;
    }
}

const Result& OfflineMeasurement::finish()
{
    if (finished) return result;
    result.momentaryMax = fromEnergy (mMaximum, result.frames >= mLength);
    result.shortTermMax = fromEnergy (sMaximum, result.frames >= sLength);
    // Separate LRA drain. It never enters I/Mmax/Smax/sample peak.
    for (int n = 0; n < result.sampleRate * 3 / 2; ++n)
    {
        double energy = 0;
        for (size_t i = 0; i < channels.size(); ++i)
        {
            auto& c = channels[i];
            const double y = c.highpass.process (c.shelf.process (0));
            energy += c.weight * y*y;
            if (n < 64) addPeak (i, 0, result.frames + n, true);
        }
        addEnergy (energy, true);
    }
    const double sp = *std::max_element (result.channelSamplePeaks.begin(), result.channelSamplePeaks.end());
    const double tp = *std::max_element (result.channelTruePeaks.begin(), result.channelTruePeaks.end());
    result.samplePeak = fromPeak (sp); result.truePeak = fromPeak (tp);
    double sum = 0;
    const auto integratedGate = relativeGateEnergy (result.momentary, false);
    if (integratedGate)
    {
        const double relative = *integratedGate;
        result.integratedRelativeGate = loudness (relative);
        sum = 0;
        for (auto w : result.momentary)
            if (includedInPopulation (w.energy, relative, false)) { sum += w.energy; ++result.integratedIncluded; }
        result.integrated = fromEnergy (sum / static_cast<double> (result.integratedIncluded), true);
    }
    else result.integrated = { {}, result.frames < mLength ? Status::insufficientDuration : sp == 0 ? Status::silence : Status::belowGate };
    const auto lraGate = relativeGateEnergy (result.lraWindows, true);
    // A source shorter than 3 s does not acquire an LRA merely through drain padding.
    if (lraGate && result.frames >= sLength)
    {
        const double relative = *lraGate;
        result.lraRelativeGate = loudness (relative);
        std::vector<double> population;
        for (auto w : result.lraWindows)
            if (includedInPopulation (w.energy, relative, true)) population.push_back (loudness (w.energy));
        std::sort (population.begin(), population.end());
        result.lraIncluded = population.size();
        auto percentile = [&] (double p) { return population[percentileIndex (population.size(), p)]; };
        result.p10 = percentile (0.10); result.p95 = percentile (0.95);
        result.lra = { *result.p95 - *result.p10, Status::available };
        result.lraStable = result.frames >= static_cast<int64_t> (result.sampleRate) * 60;
    }
    else result.lra = { {}, result.frames < sLength ? Status::insufficientDuration : sp == 0 ? Status::silence : Status::belowGate };
    finished = true;
    return result;
}
}

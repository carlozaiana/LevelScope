#pragma once

#include "../../BS1770KWeighting.h"

#include <array>
#include <cstdint>
#include <cmath>
#include <algorithm>

// Momentary loudness detector (BS.1770-style):
// - K-weighting per channel
// - energy accumulation into 60 Hz frames
// - momentary window = 0.4 s = 24 frames @ 60 Hz
//
// Notes:
// - This is not a full EBU R128 implementation (no gating, no integrated).
// - Intended as a stable control signal for dynamics stages (Upward).
// - RT-safe: no allocations in process path (allocations only in prepare() via BS1770KWeighting).
namespace levelscope::dsp
{
class BS1770MomentaryLufsDetector
{
public:
    static constexpr int kMaxChannels = 16;
    static constexpr double kFrameRateHz = 60.0;
    static constexpr int kMomentaryFrames = 24; // 0.4s @ 60Hz

    void prepare (double sampleRate, int numChannels)
    {
        fs = (sampleRate > 0.0 ? sampleRate : 48000.0);
        channels = std::clamp (numChannels, 1, kMaxChannels);

        frameSamples = std::max (1, (int) std::lround (fs / kFrameRateHz));
        samplesUntilNextFrame = frameSamples;

        kWeight.prepare (fs, channels);
        reset();
    }

    void reset() noexcept
    {
        kWeight.reset();

        frameEnergyAccum.fill (0.0);

        for (auto& chHist : frameEnergyHistory)
            chHist.fill (0.0f);

        momentaryEnergySums.fill (0.0f);

        writeIndex = 0;
        framesWritten = 0;
        samplesUntilNextFrame = frameSamples;
    }

    // Process one channel sample (returns K-weighted sample, useful for instantaneous guards).
    float processSample (int ch, float x) noexcept
    {
        if (ch < 0 || ch >= channels)
            return x;

        const float y = kWeight.processSample (ch, x);
        frameEnergyAccum[(size_t) ch] += (double) y * (double) y;
        return y;
    }

    // Call once per audio sample (after processing all channels for that sample).
    // Returns true if a 60 Hz detector frame was pushed this sample.
    bool advanceSample() noexcept
    {
        if (--samplesUntilNextFrame > 0)
            return false;

        samplesUntilNextFrame += frameSamples;
        pushFrameNoAlloc();
        return true;
    }

    float getMomentaryLufsForChannel (int ch) const noexcept
    {
        if (ch < 0 || ch >= channels)
            return -200.0f;

        const int count = std::min (framesWritten, kMomentaryFrames);
        if (count <= 0)
            return -200.0f;

        const float e = momentaryEnergySums[(size_t) ch] / (float) count;
        return energyToLufsSafe (e);
    }

    float getMomentaryLufsForMask (uint16_t maskBits, int numChannelsForMask) const noexcept
    {
        const int n = std::clamp (numChannelsForMask, 0, channels);
        const int count = std::min (framesWritten, kMomentaryFrames);

        if (n <= 0 || count <= 0 || maskBits == 0)
            return -200.0f;

        float totalEnergy = 0.0f;

        for (int ch = 0; ch < n; ++ch)
        {
            const uint16_t bit = (uint16_t) (1u << (unsigned) ch);
            if ((maskBits & bit) == 0)
                continue;

            totalEnergy += (momentaryEnergySums[(size_t) ch] / (float) count);
        }

        return energyToLufsSafe (totalEnergy);
    }

    int getNumChannels() const noexcept { return channels; }

private:
    static float energyToLufsSafe (float energy) noexcept
    {
        if (! std::isfinite (energy) || energy <= 0.0f)
            return -200.0f;

        return -0.691f + 10.0f * std::log10 (energy);
    }

    void pushFrameNoAlloc() noexcept
    {
        for (int ch = 0; ch < channels; ++ch)
        {
            const float newEnergy =
                (frameSamples > 0 ? (float) std::max (0.0, frameEnergyAccum[(size_t) ch] / (double) frameSamples)
                                  : 0.0f);

            const float old =
                (framesWritten >= kMomentaryFrames
                    ? frameEnergyHistory[(size_t) ch][(size_t) writeIndex]
                    : 0.0f);

            momentaryEnergySums[(size_t) ch] += newEnergy - old;
            frameEnergyHistory[(size_t) ch][(size_t) writeIndex] = newEnergy;

            frameEnergyAccum[(size_t) ch] = 0.0;
        }

        writeIndex = (writeIndex + 1) % kMomentaryFrames;
        ++framesWritten;
    }

    double fs = 48000.0;
    int channels = 2;

    int frameSamples = 800;
    int samplesUntilNextFrame = 800;

    BS1770KWeighting kWeight;

    std::array<double, kMaxChannels> frameEnergyAccum {};
    std::array<std::array<float, kMomentaryFrames>, kMaxChannels> frameEnergyHistory {};
    std::array<float, kMaxChannels> momentaryEnergySums {};

    int writeIndex = 0;
    int framesWritten = 0;
};
} // namespace levelscope::dsp
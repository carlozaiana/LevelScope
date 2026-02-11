#include "LookaheadLimiter.h"

namespace levelscope::dsp
{
void LookaheadLimiter::prepare (double sampleRate,
                                int numChannels,
                                const juce::AudioChannelSet& channelSet,
                                int maxBlockSize)
{
    fs = (sampleRate > 0.0 ? sampleRate : 48000.0);
    preparedNumChannels = std::max (1, numChannels);
    preparedChannelSet = channelSet;
    preparedMaxBlockSize = std::max (0, maxBlockSize);

    const int maxLookaheadSamples = (int) std::ceil ((double) fs * (double) kMaxLookaheadMs * 0.001);
    delayCapacity = std::max (1, maxLookaheadSamples + 1);

    delay.clear();
    delay.resize ((size_t) preparedNumChannels);

    for (auto& ch : delay)
        ch.buf.assign ((size_t) delayCapacity, 0.0f);

    gainDelay.assign ((size_t) delayCapacity, 1.0f);
    inputScratch.assign ((size_t) preparedNumChannels, 0.0f);

    writePos = 0;

    lastLookaheadMs = -999.0f;
    lastReleaseMs = -999.0f;
    lookaheadSamples = 0;

    updateLookaheadIfNeeded();
    updateReleaseCoeffIfNeeded();

    reset();
}

void LookaheadLimiter::reset() noexcept
{
    for (auto& ch : delay)
        std::fill (ch.buf.begin(), ch.buf.end(), 0.0f);

    std::fill (gainDelay.begin(), gainDelay.end(), 1.0f);

    writePos = 0;
    gainZ = 1.0f;
    pendingReset = false;
}

void LookaheadLimiter::setParametersAudioThread (const Parameters& p) noexcept
{
    params = p;

    // Changes to lookahead imply delay topology change -> reset safely (rare user action).
    updateLookaheadIfNeeded();
    updateReleaseCoeffIfNeeded();

    if (pendingReset)
        reset();
}

void LookaheadLimiter::updateLookaheadIfNeeded() noexcept
{
    const float ms = juce::jlimit (0.0f, kMaxLookaheadMs, params.lookaheadMs);

    if (ms == lastLookaheadMs)
        return;

    lastLookaheadMs = ms;

    const int newSamples = (int) std::lround ((double) fs * (double) ms * 0.001);

    if (newSamples != lookaheadSamples)
    {
        lookaheadSamples = juce::jlimit (0, std::max (0, delayCapacity - 1), newSamples);
        pendingReset = true;
    }
}

void LookaheadLimiter::updateReleaseCoeffIfNeeded() noexcept
{
    const float rel = std::max (5.0f, params.releaseMs);

    if (rel == lastReleaseMs)
        return;

    lastReleaseMs = rel;

    const float relS = std::max (1.0e-4f, rel * 0.001f);
    const float sr = (float) std::max (1.0, fs);

    // per-sample one-pole release coefficient
    aRelease = std::exp (-1.0f / (relS * sr));
}

void LookaheadLimiter::process (juce::AudioBuffer<float>& buffer) noexcept
{
    juce::ScopedNoDenormals noDenormals;

    if (! params.enabled)
        return;

    const int numSamples = buffer.getNumSamples();
    const int numChInBuf = buffer.getNumChannels();

    if (numSamples <= 0 || numChInBuf <= 0 || preparedNumChannels <= 0)
        return;

    const int chToProcess = std::min (preparedNumChannels, numChInBuf);

    const float ceilingLin = dbToLin (juce::jmin (0.0f, params.ceilingDb));

    float* const* chans = buffer.getArrayOfWritePointers();

    if (lookaheadSamples <= 0)
    {
        // No lookahead: apply immediate gain (still useful as a simple limiter).
        for (int i = 0; i < numSamples; ++i)
        {
            float peak = 0.0f;
            for (int ch = 0; ch < chToProcess; ++ch)
            {
                const float x = chans[ch][i];
                peak = std::max (peak, std::abs (x));
            }

            const float target = (peak > ceilingLin ? (ceilingLin / (peak + 1.0e-12f)) : 1.0f);

            // Limiter attack is instantaneous (drop immediately), release is smoothed upward
            if (target < gainZ)
                gainZ = target;
            else
                gainZ = aRelease * gainZ + (1.0f - aRelease) * target;

            for (int ch = 0; ch < chToProcess; ++ch)
                chans[ch][i] *= gainZ;
        }

        return;
    }

    // Lookahead: schedule gain alongside audio delay line.
    for (int i = 0; i < numSamples; ++i)
    {
        // Capture input before we overwrite the buffer with delayed output
        float peak = 0.0f;
        for (int ch = 0; ch < chToProcess; ++ch)
        {
            const float x = chans[ch][i];
            inputScratch[(size_t) ch] = x;
            peak = std::max (peak, std::abs (x));
        }

        const float target = (peak > ceilingLin ? (ceilingLin / (peak + 1.0e-12f)) : 1.0f);

        if (target < gainZ)
            gainZ = target;
        else
            gainZ = aRelease * gainZ + (1.0f - aRelease) * target;

        const int readPos = (writePos - lookaheadSamples + delayCapacity) % delayCapacity;

        const float gOut = gainDelay[(size_t) readPos];

        // Output delayed samples with delayed gain
        for (int ch = 0; ch < chToProcess; ++ch)
        {
            const float y = delay[(size_t) ch].buf[(size_t) readPos] * gOut;
            chans[ch][i] = y;
        }

        // Write current input and gain into delay lines
        for (int ch = 0; ch < chToProcess; ++ch)
            delay[(size_t) ch].buf[(size_t) writePos] = inputScratch[(size_t) ch];

        gainDelay[(size_t) writePos] = gainZ;

        writePos = (writePos + 1) % delayCapacity;
    }
}
} // namespace levelscope::dsp
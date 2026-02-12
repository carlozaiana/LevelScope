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
    // [BEGIN LS-LIM-PREPARE-PREVDRIVEN]
    prevDriven.assign ((size_t) preparedNumChannels, 0.0f);
    // [END LS-LIM-PREPARE-PREVDRIVEN]

    writePos = 0;

    lastLookaheadMs = -999.0f;
    lastReleaseMs = -999.0f;
    lookaheadSamples = 0;

    updateLookaheadIfNeeded();
    updateReleaseCoeffIfNeeded();
    // [BEGIN LS-LIM-PREPARE-TP]
    updateAttackIfNeeded();
    updateDriveIfNeeded();
    // [END LS-LIM-PREPARE-TP]

    reset();
}

void LookaheadLimiter::reset() noexcept
{
    for (auto& ch : delay)
        std::fill (ch.buf.begin(), ch.buf.end(), 0.0f);

    std::fill (gainDelay.begin(), gainDelay.end(), 1.0f);

    writePos = 0;
    gainZ = 1.0f;
    // [BEGIN LS-LIM-RESET-PREVDRIVEN]
    std::fill (prevDriven.begin(), prevDriven.end(), 0.0f);
    // [END LS-LIM-RESET-PREVDRIVEN]
    pendingReset = false;
}

void LookaheadLimiter::setParametersAudioThread (const Parameters& p) noexcept
{
    params = p;

    // Changes to lookahead imply delay topology change -> reset safely (rare user action).
    updateLookaheadIfNeeded();
    updateReleaseCoeffIfNeeded();
    // [BEGIN LS-LIM-SETPARAMS-TP]
    updateAttackIfNeeded();
    updateDriveIfNeeded();
    // [END LS-LIM-SETPARAMS-TP]

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

// [BEGIN LS-LIM-HELPERS-TP-IMPL]
void LookaheadLimiter::updateAttackIfNeeded() noexcept
{
    const float atk = juce::jlimit (0.0f, 5.0f, params.attackMs);

    if (atk == lastAttackMs)
        return;

    lastAttackMs = atk;

    const double sr = std::max (1.0, fs);
    const int atkSamp = (int) std::lround (sr * (double) atk * 0.001);

    // Only meaningful if we have lookahead; clamp to lookahead window.
    attackSamples = (lookaheadSamples > 0 ? juce::jlimit (0, lookaheadSamples, atkSamp) : 0);

    // We use a simple linear ramp in the ring-buffer; no one-pole needed here.
    aAttack = 0.0f;
}

void LookaheadLimiter::updateDriveIfNeeded() noexcept
{
    const float d = params.driveDb;

    if (d == lastDriveDb)
        return;

    lastDriveDb = d;
    driveLin = dbToLin (d);
}

// Detector peak with optional 2x/4x linear-interp oversampling between prev and current driven samples.
float LookaheadLimiter::computeLinkedPeakDetector (const float* const* chans,
                                                   int chToProcess,
                                                   int sampleIndex) noexcept
{
    const int os = juce::jlimit (0, 2, params.oversamplingChoice); // 0=off,1=2x,2=4x

    float peak = 0.0f;

    for (int ch = 0; ch < chToProcess; ++ch)
    {
        const float x0 = prevDriven[(size_t) ch];
        const float x1 = chans[ch][sampleIndex] * driveLin;

        // always consider current sample
        peak = std::max (peak, std::abs (x1));

        if (os == 1) // 2x
        {
            const float mid = 0.5f * (x0 + x1);
            peak = std::max (peak, std::abs (mid));
        }
        else if (os == 2) // 4x
        {
            const float d = (x1 - x0);
            const float q1 = x0 + 0.25f * d;
            const float q2 = x0 + 0.50f * d;
            const float q3 = x0 + 0.75f * d;
            peak = std::max (peak, std::abs (q1));
            peak = std::max (peak, std::abs (q2));
            peak = std::max (peak, std::abs (q3));
        }

        prevDriven[(size_t) ch] = x1; // update prev for next sample
    }

    return peak;
}
// [END LS-LIM-HELPERS-TP-IMPL]

// [BEGIN LS-LIM-PROCESS-TP]
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

    // Ensure derived values are up to date (safe: no allocations)
    updateLookaheadIfNeeded();
    updateReleaseCoeffIfNeeded();
    updateAttackIfNeeded();
    updateDriveIfNeeded();

    if (lookaheadSamples <= 0)
    {
        // No lookahead: hard limiter behavior (attack smoothing cannot be guaranteed).
        for (int i = 0; i < numSamples; ++i)
        {
            float peak = 0.0f;
            for (int ch = 0; ch < chToProcess; ++ch)
            {
                const float x = chans[ch][i] * driveLin;
                peak = std::max (peak, std::abs (x));
            }

            const float required = (peak > ceilingLin ? (ceilingLin / (peak + 1.0e-12f)) : 1.0f);

            // Instant drop, smoothed release
            if (required < gainZ)
                gainZ = required;
            else
                gainZ = aRelease * gainZ + (1.0f - aRelease) * 1.0f;

            const float gOut = gainZ;

            for (int ch = 0; ch < chToProcess; ++ch)
                chans[ch][i] = (chans[ch][i] * driveLin) * gOut;
        }

        return;
    }

    // Lookahead limiter with attack ramp:
    // - schedule per-sample gain in ring buffer
    // - when a lower gain is required at writePos, ramp the scheduled gains down over attackSamples
    //   leading up to writePos (within the lookahead window).
    for (int i = 0; i < numSamples; ++i)
    {
        // Detector on driven signal (optionally oversampled)
        const float peak = computeLinkedPeakDetector ((const float* const*) chans, chToProcess, i);

        const float required = (peak > ceilingLin ? (ceilingLin / (peak + 1.0e-12f)) : 1.0f);

        // Gain state: instant reduction to required, exponential recovery toward 1
        if (required < gainZ)
            gainZ = required;
        else
            gainZ = aRelease * gainZ + (1.0f - aRelease) * 1.0f;

        const float scheduled = std::min (gainZ, required);

        // Store driven audio at writePos (for later output)
        for (int ch = 0; ch < chToProcess; ++ch)
            delay[(size_t) ch].buf[(size_t) writePos] = chans[ch][i] * driveLin;

        // Write scheduled gain for this sample
        gainDelay[(size_t) writePos] = std::min (gainDelay[(size_t) writePos], scheduled);

        // Attack ramp: adjust previous scheduled gains in the ring (not yet output)
        if (attackSamples > 0 && scheduled < 0.999999f)
        {
            for (int k = 1; k <= attackSamples; ++k)
            {
                const int pos = (writePos - k + delayCapacity) % delayCapacity;

                const float t = 1.0f - ((float) k / (float) attackSamples); // 1..0
                const float rampG = 1.0f + (scheduled - 1.0f) * t;

                gainDelay[(size_t) pos] = std::min (gainDelay[(size_t) pos], rampG);
            }
        }

        const int readPos = (writePos - lookaheadSamples + delayCapacity) % delayCapacity;
        const float gOut = gainDelay[(size_t) readPos];

        for (int ch = 0; ch < chToProcess; ++ch)
            chans[ch][i] = delay[(size_t) ch].buf[(size_t) readPos] * gOut;

        // After reading, clear the gain slot so future scheduling can restore toward unity
        gainDelay[(size_t) readPos] = 1.0f;

        writePos = (writePos + 1) % delayCapacity;
    }
}
// [END LS-LIM-PROCESS-TP]
} // namespace levelscope::dsp
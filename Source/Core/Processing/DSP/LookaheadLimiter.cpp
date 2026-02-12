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
    // [BEGIN LS-LIM-PREPARE-FIR-OVERSAMPLERS]
    // FIR oversamplers for true-peak-ish detection (detector path).
    // 2x = 1 stage, 4x = 2 stages.
    os2 = std::make_unique<juce::dsp::Oversampling<float>> (
        (size_t) preparedNumChannels, 1,
        juce::dsp::Oversampling<float>::filterHalfBandFIREquiripple, true);

    os4 = std::make_unique<juce::dsp::Oversampling<float>> (
        (size_t) preparedNumChannels, 2,
        juce::dsp::Oversampling<float>::filterHalfBandFIREquiripple, true);

    os2->initProcessing ((size_t) juce::jmax (1, preparedMaxBlockSize));
    os4->initProcessing ((size_t) juce::jmax (1, preparedMaxBlockSize));

    detectorBuffer.setSize (preparedNumChannels, juce::jmax (1, preparedMaxBlockSize), false, false, true);

    activeOs = nullptr;
    osFactor = 1;
    detectorDelaySamples = 0;
    // [END LS-LIM-PREPARE-FIR-OVERSAMPLERS]

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
    // [BEGIN LS-LIM-RESET-FIR-OVERSAMPLERS]
    if (os2) os2->reset();
    if (os4) os4->reset();
    // [END LS-LIM-RESET-FIR-OVERSAMPLERS]
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

// [BEGIN LS-LIM-SELECT-ACTIVE-OS]
static int osChoiceToFactor (int choice) noexcept
{
    if (choice == 1) return 2;
    if (choice == 2) return 4;
    return 1;
}
// [END LS-LIM-SELECT-ACTIVE-OS]

// [BEGIN LS-LIM-SAMPLEPEAK-DETECTOR]
float LookaheadLimiter::computeLinkedPeakSamplePeak (float* const* chans,
                                                     int chToProcess,
                                                     int sampleIndex) noexcept
{
    float peak = 0.0f;
    for (int ch = 0; ch < chToProcess; ++ch)
    {
        const float x = chans[ch][sampleIndex] * driveLin;
        peak = std::max (peak, std::abs (x));
    }
    return peak;
}
// [END LS-LIM-SAMPLEPEAK-DETECTOR]

// [BEGIN LS-LIM-PROCESS-FIR-TRUEPEAK]
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

    float* const* chans = buffer.getArrayOfWritePointers();

    // Keep derived values fresh (no allocations)
    updateLookaheadIfNeeded();
    updateReleaseCoeffIfNeeded();
    updateAttackIfNeeded();
    updateDriveIfNeeded();

    // Select FIR oversampler for detector (Off/2x/4x)
    const int osChoice = juce::jlimit (0, 2, params.oversamplingChoice);
    const int factor = osChoiceToFactor (osChoice);

    if (factor == 2)
        activeOs = (os2 ? os2.get() : nullptr);
    else if (factor == 4)
        activeOs = (os4 ? os4.get() : nullptr);
    else
        activeOs = nullptr;

    osFactor = factor;
    detectorDelaySamples = (activeOs != nullptr ? activeOs->getLatencyInSamples() : 0);

    // Effective limiter latency = lookahead + detector FIR delay
    const int effectiveDelay = lookaheadSamples + detectorDelaySamples;

    const float ceilingLin = dbToLin (juce::jmin (0.0f, params.ceilingDb));

    // If lookahead is zero, we fall back to immediate limiting (still uses drive and optional sample-peak).
    if (effectiveDelay <= 0)
    {
        for (int i = 0; i < numSamples; ++i)
        {
            const float peak = computeLinkedPeakSamplePeak (chans, chToProcess, i);
            const float required = (peak > ceilingLin ? (ceilingLin / (peak + 1.0e-12f)) : 1.0f);

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

    // FIR detector path: build driven detectorBuffer (before we overwrite output with delayed samples)
    juce::dsp::AudioBlock<float> osBlock;

    if (activeOs != nullptr)
    {
        for (int ch = 0; ch < chToProcess; ++ch)
        {
            const float* in = buffer.getReadPointer (ch);
            float* dst = detectorBuffer.getWritePointer (ch);

            for (int i = 0; i < numSamples; ++i)
                dst[i] = in[i] * driveLin;
        }

        juce::dsp::AudioBlock<float> inBlock (detectorBuffer);
        auto slice = inBlock.getSubBlock (0, (size_t) numSamples);
        osBlock = activeOs->processSamplesUp (slice);
    }

    // Per-sample processing: schedule gain for the sample that will be output after effectiveDelay.
    for (int i = 0; i < numSamples; ++i)
    {
        float peak = 0.0f;

        if (activeOs != nullptr)
        {
            // Oversampled detector: peak over this base-rate slice across channels.
            const int osStart = i * osFactor;
            const int osEnd   = osStart + osFactor;

            for (int ch = 0; ch < chToProcess; ++ch)
            {
                const float* osCh = osBlock.getChannelPointer ((size_t) ch);
                for (int j = osStart; j < osEnd; ++j)
                    peak = std::max (peak, std::abs (osCh[(size_t) j]));
            }
        }
        else
        {
            // Sample peak detector
            peak = computeLinkedPeakSamplePeak (chans, chToProcess, i);
        }

        const float required = (peak > ceilingLin ? (ceilingLin / (peak + 1.0e-12f)) : 1.0f);

        // Gain state: instantaneous down, smoothed up
        if (required < gainZ)
            gainZ = required;
        else
            gainZ = aRelease * gainZ + (1.0f - aRelease) * 1.0f;

        const float scheduled = std::min (gainZ, required);

        // Store driven audio at writePos
        for (int ch = 0; ch < chToProcess; ++ch)
            delay[(size_t) ch].buf[(size_t) writePos] = buffer.getReadPointer (ch)[i] * driveLin;

        // [BEGIN LS-LIM-FIR-SCHEDULE-WITH-DETECTOR-DELAY]
        // IMPORTANT: FIR oversampling has detector latency. The "required" gain we computed at time i
        // corresponds to an earlier sample by detectorDelaySamples. Therefore schedule gain at:
        //   gainWritePos = writePos - detectorDelaySamples
        const int gainWritePos = (writePos - detectorDelaySamples + delayCapacity) % delayCapacity;

        gainDelay[(size_t) gainWritePos] = std::min (gainDelay[(size_t) gainWritePos], scheduled);

        // Rounded attack: ramp scheduled gain backwards over attackSamples leading up to gainWritePos.
        if (attackSamples > 0 && scheduled < 0.999999f)
        {
            for (int k = 1; k <= attackSamples; ++k)
            {
                const int pos = (gainWritePos - k + delayCapacity) % delayCapacity;

                const float t = 1.0f - ((float) k / (float) attackSamples); // 1..0
                const float rampG = 1.0f + (scheduled - 1.0f) * t;

                gainDelay[(size_t) pos] = std::min (gainDelay[(size_t) pos], rampG);
            }
        }
        // [END LS-LIM-FIR-SCHEDULE-WITH-DETECTOR-DELAY]

        const int readPos = (writePos - effectiveDelay + delayCapacity) % delayCapacity;
        const float gOut = gainDelay[(size_t) readPos];

        for (int ch = 0; ch < chToProcess; ++ch)
            chans[ch][i] = delay[(size_t) ch].buf[(size_t) readPos] * gOut;

        // Clear after use (so future scheduling can recover to unity)
        gainDelay[(size_t) readPos] = 1.0f;

        writePos = (writePos + 1) % delayCapacity;
    }
}
// [END LS-LIM-PROCESS-FIR-TRUEPEAK]
} // namespace levelscope::dsp
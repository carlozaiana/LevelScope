#include "SpectralUpwardCompressor.h"
// [BEGIN LS-SUC-UPWARD-GAINLAW-INCLUDE]
#include "../DSP/UpwardGainLaw.h"
// [END LS-SUC-UPWARD-GAINLAW-INCLUDE]

namespace levelscope::dsp
{
float SpectralUpwardCompressor::softKnee01 (float levelDb, float threshold, float kneeWidthDb) noexcept
{
    kneeWidthDb = std::max (1.0e-4f, kneeWidthDb);
    const float d = levelDb - threshold;

    if (d < -kneeWidthDb) return 0.0f;
    if (d >  kneeWidthDb) return 1.0f;

    const float t = (d + kneeWidthDb) / (2.0f * kneeWidthDb);
    return t * t * (3.0f - 2.0f * t);
}

void SpectralUpwardCompressor::prepare (double sampleRate,
                                       int numChannels,
                                       const juce::AudioChannelSet& channelSet,
                                       int maxBlockSize)
{
    fs = (sampleRate > 0.0 ? sampleRate : 48000.0);
    preparedNumChannels = std::max (1, numChannels);
    preparedChannelSet = channelSet;
    preparedMaxBlockSize = std::max (0, maxBlockSize);

    // Prebuild FFT profiles for all choices (so we can switch at runtime without allocating).
    for (int pi = 0; pi < kNumFftChoices; ++pi)
    {
        auto& pr = profiles[(size_t) pi];

        pr.fftSize = kFftSizes[pi];
        pr.hopSize = pr.fftSize / 4;
        pr.overlapCount = (pr.hopSize > 0 ? pr.fftSize / pr.hopSize : 4);

        pr.fft = std::make_unique<juce::dsp::FFT> ((int) std::log2 ((double) pr.fftSize));

        // sqrt-Hann window
        pr.window.resize ((size_t) pr.fftSize);
        for (int n = 0; n < pr.fftSize; ++n)
        {
            const double hann =
                0.5 - 0.5 * std::cos (2.0 * juce::MathConstants<double>::pi * (double) n / (double) (pr.fftSize - 1));
            pr.window[(size_t) n] = (float) std::sqrt (std::max (0.0, hann));
        }

        // coherent gain (sum window / N)
        pr.coherentGain = 0.0;
        for (float w : pr.window)
            pr.coherentGain += (double) w;
        pr.coherentGain /= (double) pr.fftSize;

        // [BEGIN LS-SUC-HOP-NORM-PREPARE]
                // COLA normalization for emitted hop samples:
                // When we output ola[0..hop-1], those samples will (in steady state) contain
                // contributions from window indices i + m*hop. We normalize at output time.
                pr.hopNorm.resize ((size_t) pr.hopSize);

                for (int i = 0; i < pr.hopSize; ++i)
                {
                    double s = 0.0;
                    for (int m = 0; m < pr.overlapCount; ++m)
                    {
                        const int idx = i + m * pr.hopSize;
                        if (idx >= 0 && idx < pr.fftSize)
                        {
                            const double w = (double) pr.window[(size_t) idx];
                            s += w * w; // window^2
                        }
                    }
                    pr.hopNorm[(size_t) i] = (float) std::max (1.0e-8, s);
                }
        // [END LS-SUC-HOP-NORM-PREPARE]
    }

    maxFftSize = kFftSizes[kNumFftChoices - 1];

    // Allocate per-channel buffers sized to max FFT size (no allocations in process()).
    ch.clear();
    ch.resize ((size_t) preparedNumChannels);

    for (auto& st : ch)
    {
        st.input.assign ((size_t) maxFftSize, 0.0f);
        st.fftBuf.assign ((size_t) (2 * maxFftSize), 0.0f);
        st.ola.assign ((size_t) maxFftSize, 0.0f);

        // FIFO capacity: 2*maxFft is plenty for startup latency and any transient mismatch.
        st.fifo.assign ((size_t) (2 * maxFftSize), 0.0f);
        st.fifoRead = st.fifoWrite = st.fifoCount = 0;
    }

    // [BEGIN LS-SUC-STAGE-E-PREPARE-BUILD-PREPARED-MASKBITS]
    preparedAllMaskBits    = 0;
    preparedNonLfeMaskBits = 0;
    preparedLfeMaskBits    = 0;

    const int nMaskCh = juce::jlimit (0, kMaxMaskChannels, preparedNumChannels);

    for (int chIdx = 0; chIdx < nMaskCh; ++chIdx)
    {
        const uint16_t b = (uint16_t) (1u << (unsigned) chIdx);
        preparedAllMaskBits |= b;

        const bool isLfe = (preparedChannelSet.getTypeOfChannel (chIdx) == juce::AudioChannelSet::LFE);
        if (isLfe) preparedLfeMaskBits |= b;
        else       preparedNonLfeMaskBits |= b;
    }

    // Safety fallback
    if (preparedNonLfeMaskBits == 0)
        preparedNonLfeMaskBits = preparedAllMaskBits;

    // Default: no external override
    externalMasksActive = false;
    externalDetectMaskBits = 0;
    externalApplyMaskBits  = 0;
    externalUnlinked       = false;

    effectiveDetectMaskBitsCached = 0;
    effectiveApplyMaskBitsCached  = 0;
    detectCount = applyCount = 0;
    // [END LS-SUC-STAGE-E-PREPARE-BUILD-PREPARED-MASKBITS]

    // Apply initial selection + derived values
    activeFftChoice = clampChoice (params.fftSizeChoice, kNumFftChoices);
    activeFftSize   = profiles[(size_t) activeFftChoice].fftSize;
    activeHopSize   = profiles[(size_t) activeFftChoice].hopSize;

    inputWritePos = 0;

    // Build bands + smoothers
    lastBandsPerOctChoice = params.bandsPerOctChoice;
    lastMinFreqHz = params.minFreqHz;
    lastMaxFreqHz = params.maxFreqHz;

    rebuildBandsNoAlloc();

    lastAttackMs  = params.attackMs;
    lastReleaseMs = params.releaseMs;
    updateSmoothersNoAlloc();

    // Offset smoother (slow)
    // [BEGIN LS-SUC-PREPARE-GLOBAL-ZONE]
    offsetSmoother.prepare (fs, activeHopSize, 4.0 /*seconds*/);
    smoothedOffsetDb = 0.0;

    // Global zone smoothing (Option A). Slow enough to avoid audible pumping.
    globalZoneSmoother.prepare (fs, activeHopSize, 0.5 /*attackSeconds*/, 0.05 /*releaseSeconds*/);
    smoothedGlobalZoneAmount01 = 1.0f;

    // [BEGIN LS-SUC-STAGE-E-PREPARE-UNLINKED-SMOOTHERS]
    for (int chIdx = 0; chIdx < kMaxMaskChannels; ++chIdx)
    {
        offsetSmootherUnlinked[(size_t) chIdx].prepare (fs, activeHopSize, 4.0);
        smoothedOffsetDbUnlinked[(size_t) chIdx] = 0.0;

        globalZoneSmootherUnlinked[(size_t) chIdx].prepare (fs, activeHopSize, 0.5, 0.05);
        smoothedGlobalZoneAmount01Unlinked[(size_t) chIdx] = 1.0f;
    }
    // [END LS-SUC-STAGE-E-PREPARE-UNLINKED-SMOOTHERS]

    reset();
    // [END LS-SUC-PREPARE-GLOBAL-ZONE]
}

void SpectralUpwardCompressor::reset() noexcept
{
    for (auto& st : ch)
    {
        std::fill (st.input.begin(), st.input.end(), 0.0f);
        std::fill (st.fftBuf.begin(), st.fftBuf.end(), 0.0f);
        std::fill (st.ola.begin(), st.ola.end(), 0.0f);
        std::fill (st.fifo.begin(), st.fifo.end(), 0.0f);

        st.fifoRead = st.fifoWrite = st.fifoCount = 0;
    }

    inputWritePos = 0;

    for (int i = 0; i < numBands; ++i)
        bandSmoothers[(size_t) i].reset();

    // [BEGIN LS-SUC-STAGE-E-RESET-UNLINKED]
    for (int chIdx = 0; chIdx < kMaxMaskChannels; ++chIdx)
    {
        offsetSmootherUnlinked[(size_t) chIdx].reset();
        smoothedOffsetDbUnlinked[(size_t) chIdx] = 0.0;

        globalZoneSmootherUnlinked[(size_t) chIdx].reset();
        smoothedGlobalZoneAmount01Unlinked[(size_t) chIdx] = 1.0f;

        for (int bi = 0; bi < numBands; ++bi)
            bandSmoothersUnlinked[(size_t) chIdx][(size_t) bi].reset();
    }
    // [END LS-SUC-STAGE-E-RESET-UNLINKED]

    // [BEGIN LS-SUC-RESET-GLOBAL-ZONE]
    offsetSmoother.reset();
    smoothedOffsetDb = 0.0;

    globalZoneSmoother.reset();
    smoothedGlobalZoneAmount01 = 1.0f;

    pendingHardReset = false;
    // [END LS-SUC-RESET-GLOBAL-ZONE]
}

void SpectralUpwardCompressor::setParametersAudioThread (const Parameters& p) noexcept
{
    params = p;

    // Switching FFT size choice changes internal latency and state. We handle it without allocations by:
    // - switching to a prebuilt profile
    // - doing a hard reset at the next safe point
    const int newChoice = clampChoice (params.fftSizeChoice, kNumFftChoices);
    if (newChoice != activeFftChoice)
    {
        activeFftChoice = newChoice;
        activeFftSize   = profiles[(size_t) activeFftChoice].fftSize;
        activeHopSize   = profiles[(size_t) activeFftChoice].hopSize;

        // [BEGIN LS-SUC-PARAMS-REPREPARE-SMOOTHERS]
        offsetSmoother.prepare (fs, activeHopSize, 4.0);
        globalZoneSmoother.prepare (fs, activeHopSize, 0.5, 0.05);

        for (int chIdx = 0; chIdx < kMaxMaskChannels; ++chIdx)
        {
            offsetSmootherUnlinked[(size_t) chIdx].prepare (fs, activeHopSize, 4.0);
            globalZoneSmootherUnlinked[(size_t) chIdx].prepare (fs, activeHopSize, 0.5, 0.05);
        }

        pendingHardReset = true;
        // [END LS-SUC-PARAMS-REPREPARE-SMOOTHERS]
    }

    // Band mapping changes: recompute into fixed storage, no allocations.
    if (params.bandsPerOctChoice != lastBandsPerOctChoice
        || params.minFreqHz != lastMinFreqHz
        || params.maxFreqHz != lastMaxFreqHz)
    {
        lastBandsPerOctChoice = params.bandsPerOctChoice;
        lastMinFreqHz = params.minFreqHz;
        lastMaxFreqHz = params.maxFreqHz;

        rebuildBandsNoAlloc();
        pendingHardReset = true; // conservative: avoids mismatched OLA when bands change
    }

    // Smoother coefficient updates (no allocations)
    if (params.attackMs != lastAttackMs || params.releaseMs != lastReleaseMs)
    {
        lastAttackMs  = params.attackMs;
        lastReleaseMs = params.releaseMs;
        updateSmoothersNoAlloc();
    }
}

// [BEGIN LS-SUC-STAGE-E-MASK-IMPL]
void SpectralUpwardCompressor::setChannelMasksAudioThread (uint16_t detectBits,
                                                           uint16_t applyBits,
                                                           bool unlinked) noexcept
{
    if (detectBits == 0 && applyBits == 0)
    {
        externalMasksActive = false;
        externalDetectMaskBits = 0;
        externalApplyMaskBits  = 0;
        externalUnlinked       = false;
        return;
    }

    externalMasksActive = true;
    externalDetectMaskBits = detectBits;
    externalApplyMaskBits  = applyBits;
    externalUnlinked       = unlinked;
}

void SpectralUpwardCompressor::rebuildIndexListsIfNeededNoAlloc (uint16_t detBits,
                                                                 uint16_t appBits) noexcept
{
    if (detBits == effectiveDetectMaskBitsCached && appBits == effectiveApplyMaskBitsCached)
        return;

    effectiveDetectMaskBitsCached = detBits;
    effectiveApplyMaskBitsCached  = appBits;

    const int nMaskCh = juce::jlimit (0, kMaxMaskChannels, preparedNumChannels);

    detectCount = 0;
    applyCount  = 0;

    for (int ch = 0; ch < nMaskCh; ++ch)
    {
        const uint16_t b = (uint16_t) (1u << (unsigned) ch);

        if ((detBits & b) != 0 && detectCount < kMaxMaskChannels)
            detectIdx[(size_t) detectCount++] = (uint8_t) ch;

        if ((appBits & b) != 0 && applyCount < kMaxMaskChannels)
            applyIdx[(size_t) applyCount++] = (uint8_t) ch;
    }

    // Safety fallbacks
    if (detectCount <= 0)
    {
        for (int ch = 0; ch < nMaskCh && detectCount < kMaxMaskChannels; ++ch)
            detectIdx[(size_t) detectCount++] = (uint8_t) ch;
    }

    if (applyCount <= 0)
    {
        for (int ch = 0; ch < nMaskCh && applyCount < kMaxMaskChannels; ++ch)
            applyIdx[(size_t) applyCount++] = (uint8_t) ch;
    }
}
// [END LS-SUC-STAGE-E-MASK-IMPL]

void SpectralUpwardCompressor::rebuildBandsNoAlloc() noexcept
{
    numBands = 0;

    const int fftSize = activeFftSize;
    const int maxBin = fftSize / 2;

    const double minHz = (double) std::max (1.0f, params.minFreqHz);
    const double maxHz = (double) std::max (minHz + 1.0, (double) params.maxFreqHz);

    const int bpo = kBandsPerOct[clampChoice (params.bandsPerOctChoice, kNumBandsPerOctChoices)];
    const double step = std::pow (2.0, 1.0 / (double) bpo);
    const double edge = std::pow (2.0, 1.0 / (2.0 * (double) bpo));

    const double df = fs / (double) fftSize;

    double fc = minHz;
    int lastEnd = 0;

    while (fc < maxHz && numBands < kMaxBands)
    {
        const double f0 = fc / edge;
        const double f1 = fc * edge;

        int b0 = (int) std::ceil (f0 / df);
        int b1 = (int) std::floor (f1 / df);

        // restrict to 1..(N/2-1) to avoid DC/Nyquist complications
        b0 = juce::jlimit (1, maxBin - 1, b0);
        b1 = juce::jlimit (1, maxBin - 1, b1);

        // enforce non-overlap to avoid double scaling bins
        if (b0 <= lastEnd)
            b0 = lastEnd + 1;

        if (b1 >= b0)
        {
            bands[(size_t) numBands] = Band { b0, b1 };
            lastEnd = b1;
            ++numBands;
        }

        fc *= step;
    }
}

void SpectralUpwardCompressor::updateSmoothersNoAlloc() noexcept
{
    // [BEGIN LS-SUC-STAGE-E-PREPARE-BANDSMOOTHERS-LINKED-AND-UNLINKED]
    for (int i = 0; i < numBands; ++i)
    {
        bandSmoothers[(size_t) i].prepare (fs, activeHopSize, params.attackMs, params.releaseMs);

        for (int chIdx = 0; chIdx < kMaxMaskChannels; ++chIdx)
            bandSmoothersUnlinked[(size_t) chIdx][(size_t) i].prepare (fs, activeHopSize, params.attackMs, params.releaseMs);
    }
    // [END LS-SUC-STAGE-E-PREPARE-BANDSMOOTHERS-LINKED-AND-UNLINKED]
}

float SpectralUpwardCompressor::computeTargetGainLinear (float bandLevelDb,
                                                         float t0SpectralDb,
                                                         float t1SpectralDb) noexcept
{
    UpwardGainLaw::Params law;
    law.t0Db       = t0SpectralDb;
    law.t1Db       = t1SpectralDb;
    law.lowKneeDb  = params.lowKneeDb;
    law.highKneeDb = params.highKneeDb;
    law.maxBoostDb = std::max (0.0f, params.maxBoostDb);
    law.curve01    = juce::jlimit (0.0f, 1.0f, params.curve);
    law.curveType  = (params.curveType == CurveType::bell ? UpwardGainLaw::CurveType::bell
                                                          : UpwardGainLaw::CurveType::monotonic);

    // NOTE: This returns the *base* gain (no amount fade). Caller applies amount in linear domain.
    return UpwardGainLaw::computeBaseGainLin (bandLevelDb, law);
}

// [BEGIN LS-SUC-STAGE-E-PROCESSFRAME-REPLACE]
void SpectralUpwardCompressor::processFrameAllChannels() noexcept
{
    if (pendingHardReset)
        reset();

    // [BEGIN LS-SUC-UPWARD-METERING-FRAME-INIT]
    lastFrameMaxLinearGain  = 1.0f;
    lastFrameLastLinearGain = 1.0f;
    // [END LS-SUC-UPWARD-METERING-FRAME-INIT]

    const auto& pr = profiles[(size_t) activeFftChoice];
    const int fftSize = pr.fftSize;
    const int hopSize = pr.hopSize;
    const int maxBin  = fftSize / 2;

    // --- Determine effective masks (external override if set; otherwise legacy LFE policy)
    const uint16_t legacyDetectBits = (params.lfeInDetector ? preparedAllMaskBits : preparedNonLfeMaskBits);
    const uint16_t legacyApplyBits  = (params.lfeInApply    ? preparedAllMaskBits : preparedNonLfeMaskBits);

    uint16_t effDetectBits = (externalMasksActive ? externalDetectMaskBits : legacyDetectBits);
    uint16_t effApplyBits  = (externalMasksActive ? externalApplyMaskBits  : legacyApplyBits);

    // Clamp to prepared channel count / 16ch range
    effDetectBits = (uint16_t) (effDetectBits & preparedAllMaskBits);
    effApplyBits  = (uint16_t) (effApplyBits  & preparedAllMaskBits);

    if (effDetectBits == 0) effDetectBits = preparedAllMaskBits;
    if (effApplyBits  == 0) effApplyBits  = preparedAllMaskBits;

    rebuildIndexListsIfNeededNoAlloc (effDetectBits, effApplyBits);

    const bool doUnlinked = (externalMasksActive && externalUnlinked);

    // --- 1) Analysis window into per-channel fft buffers.
    // Also compute broadband proxy energies.
    std::array<double, kMaxMaskChannels> sumSqWinPerCh {};
    double sumSqLinked = 0.0;

    for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
    {
        auto& st = ch[(size_t) chIdx];

        double sumSqThis = 0.0;

        for (int i = 0; i < fftSize; ++i)
        {
            const float x = st.input[(size_t) i] * pr.window[(size_t) i];
            st.fftBuf[(size_t) i] = x;
            sumSqThis += (double) x * (double) x;
        }

        std::fill (st.fftBuf.begin() + (size_t) fftSize,
                   st.fftBuf.begin() + (size_t) (2 * fftSize),
                   0.0f);

        if (chIdx < kMaxMaskChannels)
            sumSqWinPerCh[(size_t) chIdx] = sumSqThis;
    }

    // --- 2) Forward FFT per channel
    for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
    {
        auto& st = ch[(size_t) chIdx];
        pr.fft->performRealOnlyForwardTransform (st.fftBuf.data());
    }

    // --- Audition bypass: preserve delay pipeline, no gains
    if (params.auditionBypass)
    {
        // [BEGIN LS-SUC-UPWARD-METERING-BYPASS-ZERO]
        // No boost applied while bypassed
        lastFrameMaxLinearGain  = 1.0f;
        lastFrameLastLinearGain = 1.0f;
        // [END LS-SUC-UPWARD-METERING-BYPASS-ZERO]
        for (int bi = 0; bi < numBands; ++bi)
            bandSmoothers[(size_t) bi].reset();

        for (int chIdx = 0; chIdx < kMaxMaskChannels; ++chIdx)
            for (int bi = 0; bi < numBands; ++bi)
                bandSmoothersUnlinked[(size_t) chIdx][(size_t) bi].reset();

        // Inverse FFT + synthesis (no spectral modifications)
        for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
        {
            auto& st = ch[(size_t) chIdx];

            pr.fft->performRealOnlyInverseTransform (st.fftBuf.data());

            for (int i = 0; i < fftSize; ++i)
                st.ola[(size_t) i] += (st.fftBuf[(size_t) i] * pr.window[(size_t) i]);

            for (int i = 0; i < hopSize; ++i)
                pushFifo (st, st.ola[(size_t) i] / pr.hopNorm[(size_t) i]);

            const int keep2 = fftSize - hopSize;
            std::memmove (st.ola.data(), st.ola.data() + hopSize, (size_t) keep2 * sizeof (float));
            std::fill (st.ola.begin() + keep2, st.ola.begin() + fftSize, 0.0f);
        }

        // shift input
        const int keep = fftSize - hopSize;
        for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
        {
            auto& st = ch[(size_t) chIdx];
            std::memmove (st.input.data(), st.input.data() + hopSize, (size_t) keep * sizeof (float));
            std::fill (st.input.begin() + keep, st.input.begin() + fftSize, 0.0f);
        }

        inputWritePos = keep;
        return;
    }

    const double ref = (pr.coherentGain * (double) fftSize * 0.5);
    const double refPower = std::max (1.0e-18, ref * ref);
    const int numBins = std::max (1, (maxBin - 1));

    const float userAmount = juce::jlimit (0.0f, 1.0f, params.amount01);

    // --- Linked vs Unlinked processing
    if (! doUnlinked)
    {
        // Linked broadband proxy over detector channels
        for (int di = 0; di < detectCount; ++di)
        {
            const int chDet = (int) detectIdx[(size_t) di];
            if (chDet >= 0 && chDet < kMaxMaskChannels)
                sumSqLinked += sumSqWinPerCh[(size_t) chDet];
        }

        const int linkedCh = std::max (1, detectCount);
        const double meanSq = sumSqLinked / (double) (fftSize * linkedCh);
        const double broadbandDb = -0.691 + 10.0 * std::log10 (meanSq + 1.0e-12);

        // Spectral proxy over detector channels
        double pAll = 0.0;
        for (int di = 0; di < detectCount; ++di)
        {
            const int chDet = (int) detectIdx[(size_t) di];
            if (chDet < 0 || chDet >= preparedNumChannels)
                continue;

            auto& st = ch[(size_t) chDet];
            for (int bin = 1; bin <= maxBin - 1; ++bin)
            {
                float re = 0.0f, im = 0.0f;
                getBin (st.fftBuf, fftSize, bin, re, im);
                pAll += (double) re * (double) re + (double) im * (double) im;
            }
        }

        const double meanBinPowerAll = pAll / (double) (numBins * linkedCh);
        const double spectralDb = 10.0 * std::log10 (meanBinPowerAll / refPower + 1.0e-18);

        // Adaptive offset (linked)
        if (meanSq > 1.0e-12)
        {
            const double instantOffset = broadbandDb - spectralDb;
            const double clamped = juce::jlimit (-30.0, 30.0, instantOffset);
            smoothedOffsetDb = offsetSmoother.process (clamped);
        }

        const double effectiveOffsetDb = smoothedOffsetDb + (double) params.calibrationTrimDb;
        const float t0SpectralDb = (float) ((double) params.t0Lufs - effectiveOffsetDb);
        const float t1SpectralDb = (float) ((double) params.t1Lufs - effectiveOffsetDb);

        // Global zone amount (linked)
        // [BEGIN LS-SUC-GLOBAL-ZONE-USING-UPWARDGAINLAW-LINKED]
        float zoneTarget01Linked = 1.0f;

        {
            UpwardGainLaw::Params zoneLaw;
            zoneLaw.t0Db       = params.t0Lufs;
            zoneLaw.t1Db       = params.t1Lufs;
            zoneLaw.lowKneeDb  = params.lowKneeDb;
            zoneLaw.highKneeDb = params.highKneeDb;
            zoneLaw.maxBoostDb = 0.0f;
            zoneLaw.curve01    = 0.0f;
            zoneLaw.curveType  = UpwardGainLaw::CurveType::monotonic;

            const float L = (float) broadbandDb;

            // 0..1 instantaneous target (0 at/above T1)
            zoneTarget01Linked = UpwardGainLaw::computeActiveZone01 (L, zoneLaw);

            // Smoothed applied zone (fast release) for "smooth return to unity"
            smoothedGlobalZoneAmount01 = globalZoneSmoother.process (zoneTarget01Linked);
        }

        // If we're at/above T1, steer band smoothers toward unity (prevents hidden gain buildup)
        const bool zoneOffNow = (zoneTarget01Linked <= 1.0e-6f);
        // [END LS-SUC-GLOBAL-ZONE-USING-UPWARDGAINLAW-LINKED]

        // Per-band linked gains (measure over detector channels, apply over apply channels)
        for (int bi = 0; bi < numBands; ++bi)
        {
            const int b0 = juce::jlimit (1, maxBin - 1, bands[(size_t) bi].startBin);
            const int b1 = juce::jlimit (1, maxBin - 1, bands[(size_t) bi].endBin);
            if (b1 < b0)
            {
                bandTargetGain[(size_t) bi] = 1.0f;
                continue;
            }

            double pBand = 0.0;
            const int nBins = (b1 - b0 + 1);

            for (int di = 0; di < detectCount; ++di)
            {
                const int chDet = (int) detectIdx[(size_t) di];
                if (chDet < 0 || chDet >= preparedNumChannels)
                    continue;

                auto& st = ch[(size_t) chDet];
                for (int bin = b0; bin <= b1; ++bin)
                {
                    float re = 0.0f, im = 0.0f;
                    getBin (st.fftBuf, fftSize, bin, re, im);
                    pBand += (double) re * (double) re + (double) im * (double) im;
                }
            }

            const double meanBinPower = pBand / (double) (nBins * std::max (1, detectCount));
            const float bandDb = (float) (10.0 * std::log10 (meanBinPower / refPower + 1.0e-18));

            // [BEGIN LS-SUC-ZONE-LAST-LINKED-TARGET]
            const float baseG = computeTargetGainLinear (bandDb, t0SpectralDb, t1SpectralDb); // >= 1, no amount applied
            float targetG = 1.0f + (baseG - 1.0f) * userAmount;

            // Above T1: steer smoother targets to unity so no hidden boost accumulates
            if (zoneOffNow)
                targetG = 1.0f;
            // [END LS-SUC-ZONE-LAST-LINKED-TARGET]

            bandTargetGain[(size_t) bi] = targetG;
        }

        // freq smooth
        if (numBands > 0)
        {
            for (int bi = 0; bi < numBands; ++bi)
            {
                const float g0 = bandTargetGain[(size_t) bi];
                const float gL = (bi > 0 ? bandTargetGain[(size_t) (bi - 1)] : g0);
                const float gR = (bi + 1 < numBands ? bandTargetGain[(size_t) (bi + 1)] : g0);
                bandTargetGainFreqSmoothed[(size_t) bi] = 0.25f * gL + 0.5f * g0 + 0.25f * gR;
            }
        }

        // time smooth + apply
        for (int bi = 0; bi < numBands; ++bi)
        {
            const int b0 = juce::jlimit (1, maxBin - 1, bands[(size_t) bi].startBin);
            const int b1 = juce::jlimit (1, maxBin - 1, bands[(size_t) bi].endBin);
            if (b1 < b0) continue;

            // [BEGIN LS-SUC-ZONE-LAST-LINKED-APPLY]
            const float gSmoothed = bandSmoothers[(size_t) bi].process (bandTargetGainFreqSmoothed[(size_t) bi]);

            // Zone last: forces unity when zone -> 0, regardless of smoother memory
            const float zone01 = juce::jlimit (0.0f, 1.0f, smoothedGlobalZoneAmount01);
            const float g = 1.0f + (gSmoothed - 1.0f) * zone01;
            // [END LS-SUC-ZONE-LAST-LINKED-APPLY]

            // [BEGIN LS-SUC-UPWARD-METERING-LINKED-UPDATE]
            lastFrameMaxLinearGain  = std::max (lastFrameMaxLinearGain,  g);
            lastFrameLastLinearGain = g;
            // [END LS-SUC-UPWARD-METERING-LINKED-UPDATE]

            for (int ai = 0; ai < applyCount; ++ai)
            {
                const int chAp = (int) applyIdx[(size_t) ai];
                if (chAp < 0 || chAp >= preparedNumChannels)
                    continue;

                auto& st = ch[(size_t) chAp];
                for (int bin = b0; bin <= b1; ++bin)
                {
                    scaleBin (st.fftBuf, fftSize, bin, g);
                    const int mirror = fftSize - bin;
                    if (mirror != bin)
                        scaleBin (st.fftBuf, fftSize, mirror, g);
                }
            }
        }
        // [BEGIN LS-SUC-ZONE-LAST-LINKED-RESET-WHEN-FULLY-OFF]
        // Once the *smoothed* zone is fully off, clear any remaining band smoother state so re-entry can't "reveal" stored gain.
        if (zoneOffNow && smoothedGlobalZoneAmount01 <= 1.0e-4f)
        {
            for (int bi = 0; bi < numBands; ++bi)
                bandSmoothers[(size_t) bi].reset();
        }
        // [END LS-SUC-ZONE-LAST-LINKED-RESET-WHEN-FULLY-OFF]
    }
    else
    {
        // --- Unlinked: per-channel detector and per-channel smoothing
        for (int ai = 0; ai < applyCount; ++ai)
        {
            const int chAp = (int) applyIdx[(size_t) ai];
            if (chAp < 0 || chAp >= preparedNumChannels || chAp >= kMaxMaskChannels)
                continue;

            const double meanSqCh = sumSqWinPerCh[(size_t) chAp] / (double) fftSize;
            const double broadbandDb = -0.691 + 10.0 * std::log10 (meanSqCh + 1.0e-12);

            // spectral proxy for this channel
            double pAll = 0.0;
            {
                auto& st = ch[(size_t) chAp];
                for (int bin = 1; bin <= maxBin - 1; ++bin)
                {
                    float re = 0.0f, im = 0.0f;
                    getBin (st.fftBuf, fftSize, bin, re, im);
                    pAll += (double) re * (double) re + (double) im * (double) im;
                }
            }

            const double meanBinPowerAll = pAll / (double) numBins;
            const double spectralDb = 10.0 * std::log10 (meanBinPowerAll / refPower + 1.0e-18);

            if (meanSqCh > 1.0e-12)
            {
                const double instantOffset = broadbandDb - spectralDb;
                const double clamped = juce::jlimit (-30.0, 30.0, instantOffset);
                smoothedOffsetDbUnlinked[(size_t) chAp] =
                    offsetSmootherUnlinked[(size_t) chAp].process (clamped);
            }

            const double effectiveOffsetDb = smoothedOffsetDbUnlinked[(size_t) chAp] + (double) params.calibrationTrimDb;
            const float t0SpectralDb = (float) ((double) params.t0Lufs - effectiveOffsetDb);
            const float t1SpectralDb = (float) ((double) params.t1Lufs - effectiveOffsetDb);

            // [BEGIN LS-SUC-GLOBAL-ZONE-USING-UPWARDGAINLAW-UNLINKED]
            const float L = (float) broadbandDb;

            UpwardGainLaw::Params zoneLaw;
            zoneLaw.t0Db       = params.t0Lufs;
            zoneLaw.t1Db       = params.t1Lufs;
            zoneLaw.lowKneeDb  = params.lowKneeDb;
            zoneLaw.highKneeDb = params.highKneeDb;
            zoneLaw.maxBoostDb = 0.0f;
            zoneLaw.curve01    = 0.0f;
            zoneLaw.curveType  = UpwardGainLaw::CurveType::monotonic;

            const float zoneTarget01 = UpwardGainLaw::computeActiveZone01 (L, zoneLaw);
            const bool zoneOffNow = (zoneTarget01 <= 1.0e-6f);

            smoothedGlobalZoneAmount01Unlinked[(size_t) chAp] =
                globalZoneSmootherUnlinked[(size_t) chAp].process (zoneTarget01);
            // [END LS-SUC-GLOBAL-ZONE-USING-UPWARDGAINLAW-UNLINKED]

            // per-band targets for this channel (reuse scratch arrays)
            for (int bi = 0; bi < numBands; ++bi)
            {
                const int b0 = juce::jlimit (1, maxBin - 1, bands[(size_t) bi].startBin);
                const int b1 = juce::jlimit (1, maxBin - 1, bands[(size_t) bi].endBin);
                if (b1 < b0)
                {
                    bandTargetGain[(size_t) bi] = 1.0f;
                    continue;
                }

                double pBand = 0.0;
                const int nBins = (b1 - b0 + 1);

                auto& st = ch[(size_t) chAp];
                for (int bin = b0; bin <= b1; ++bin)
                {
                    float re = 0.0f, im = 0.0f;
                    getBin (st.fftBuf, fftSize, bin, re, im);
                    pBand += (double) re * (double) re + (double) im * (double) im;
                }

                const double meanBinPower = pBand / (double) nBins;
                const float bandDb = (float) (10.0 * std::log10 (meanBinPower / refPower + 1.0e-18));

                // [BEGIN LS-SUC-ZONE-LAST-UNLINKED-TARGET]
                const float baseG = computeTargetGainLinear (bandDb, t0SpectralDb, t1SpectralDb);
                float targetG = 1.0f + (baseG - 1.0f) * userAmount;

                if (zoneOffNow)
                    targetG = 1.0f;
                // [END LS-SUC-ZONE-LAST-UNLINKED-TARGET]
                bandTargetGain[(size_t) bi] = targetG;
            }

            // freq smooth
            if (numBands > 0)
            {
                for (int bi = 0; bi < numBands; ++bi)
                {
                    const float g0 = bandTargetGain[(size_t) bi];
                    const float gL = (bi > 0 ? bandTargetGain[(size_t) (bi - 1)] : g0);
                    const float gR = (bi + 1 < numBands ? bandTargetGain[(size_t) (bi + 1)] : g0);
                    bandTargetGainFreqSmoothed[(size_t) bi] = 0.25f * gL + 0.5f * g0 + 0.25f * gR;
                }
            }

            // time smooth + apply bins (this channel only)
            for (int bi = 0; bi < numBands; ++bi)
            {
                const int b0 = juce::jlimit (1, maxBin - 1, bands[(size_t) bi].startBin);
                const int b1 = juce::jlimit (1, maxBin - 1, bands[(size_t) bi].endBin);
                if (b1 < b0) continue;

                // [BEGIN LS-SUC-ZONE-LAST-UNLINKED-APPLY]
                const float gSmoothed =
                    bandSmoothersUnlinked[(size_t) chAp][(size_t) bi].process (bandTargetGainFreqSmoothed[(size_t) bi]);

                const float zone01 = juce::jlimit (0.0f, 1.0f, smoothedGlobalZoneAmount01Unlinked[(size_t) chAp]);
                const float g = 1.0f + (gSmoothed - 1.0f) * zone01;
                // [END LS-SUC-ZONE-LAST-UNLINKED-APPLY]

                // [BEGIN LS-SUC-UPWARD-METERING-UNLINKED-UPDATE]
                lastFrameMaxLinearGain  = std::max (lastFrameMaxLinearGain,  g);
                lastFrameLastLinearGain = g;
                // [END LS-SUC-UPWARD-METERING-UNLINKED-UPDATE]

                auto& st = ch[(size_t) chAp];

                for (int bin = b0; bin <= b1; ++bin)
                {
                    scaleBin (st.fftBuf, fftSize, bin, g);
                    const int mirror = fftSize - bin;
                    if (mirror != bin)
                        scaleBin (st.fftBuf, fftSize, mirror, g);
                }
            }
            // [BEGIN LS-SUC-ZONE-LAST-UNLINKED-RESET-WHEN-FULLY-OFF]
            if (zoneOffNow && smoothedGlobalZoneAmount01Unlinked[(size_t) chAp] <= 1.0e-4f)
            {
                for (int bi = 0; bi < numBands; ++bi)
                    bandSmoothersUnlinked[(size_t) chAp][(size_t) bi].reset();
            }
            // [END LS-SUC-ZONE-LAST-UNLINKED-RESET-WHEN-FULLY-OFF]
        }
    }

    // --- 6) Inverse FFT + synthesis
    for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
    {
        auto& st = ch[(size_t) chIdx];

        pr.fft->performRealOnlyInverseTransform (st.fftBuf.data());

        for (int i = 0; i < fftSize; ++i)
            st.ola[(size_t) i] += (st.fftBuf[(size_t) i] * pr.window[(size_t) i]);

        for (int i = 0; i < hopSize; ++i)
            pushFifo (st, st.ola[(size_t) i] / pr.hopNorm[(size_t) i]);

        const int keep = fftSize - hopSize;
        std::memmove (st.ola.data(), st.ola.data() + hopSize, (size_t) keep * sizeof (float));
        std::fill (st.ola.begin() + keep, st.ola.begin() + fftSize, 0.0f);
    }

    // --- 7) Shift input buffers
    const int keep = fftSize - hopSize;
    for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
    {
        auto& st = ch[(size_t) chIdx];
        std::memmove (st.input.data(), st.input.data() + hopSize, (size_t) keep * sizeof (float));
        std::fill (st.input.begin() + keep, st.input.begin() + fftSize, 0.0f);
    }

    inputWritePos = keep;
}
// [END LS-SUC-STAGE-E-PROCESSFRAME-REPLACE]

void SpectralUpwardCompressor::process (juce::AudioBuffer<float>& buffer) noexcept
{
    juce::ScopedNoDenormals noDenormals;

    const int numSamples = buffer.getNumSamples();
    const int numChInBuf = buffer.getNumChannels();

    if (numSamples <= 0 || numChInBuf <= 0 || preparedNumChannels <= 0)
        return;

    const int chToProcess = std::min (preparedNumChannels, numChInBuf);

    // [BEGIN LS-SUC-UPWARD-METERING-BLOCK-INIT]
    bool hadFrameThisCall = false;
    float blockMaxG  = 1.0f;
    float blockLastG = lastBlockLastLinearGain;
    // [END LS-SUC-UPWARD-METERING-BLOCK-INIT]

    // Sample-by-sample framing keeps all channels aligned and avoids extra ring buffers.
    for (int i = 0; i < numSamples; ++i)
    {
        // Push input into per-channel frame buffers
        for (int chIdx = 0; chIdx < chToProcess; ++chIdx)
        {
            auto& st = ch[(size_t) chIdx];
            st.input[(size_t) inputWritePos] = buffer.getReadPointer (chIdx)[i];
        }

        ++inputWritePos;

        // [BEGIN LS-SUC-UPWARD-METERING-CAPTURE-FRAME]
        if (inputWritePos >= activeFftSize)
        {
            processFrameAllChannels();

            hadFrameThisCall = true;
            blockMaxG  = std::max (blockMaxG, lastFrameMaxLinearGain);
            blockLastG = lastFrameLastLinearGain;
        }
        // [END LS-SUC-UPWARD-METERING-CAPTURE-FRAME]

        // Pop output
        for (int chIdx = 0; chIdx < chToProcess; ++chIdx)
        {
            auto& st = ch[(size_t) chIdx];
            buffer.getWritePointer (chIdx)[i] = (st.fifoCount > 0 ? popFifo (st) : 0.0f);
        }
    }
    // [BEGIN LS-SUC-UPWARD-METERING-BLOCK-STORE]
    if (hadFrameThisCall)
    {
        lastBlockMaxLinearGain  = std::max (1.0f, blockMaxG);
        lastBlockLastLinearGain = std::max (1.0f, blockLastG);
    }
    // [END LS-SUC-UPWARD-METERING-BLOCK-STORE]
}
} // namespace levelscope::dsp
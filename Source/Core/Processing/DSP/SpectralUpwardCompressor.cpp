#include "SpectralUpwardCompressor.h"

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

    // [BEGIN LS-SUC-PREPARE-BUILD-ALLCHANNELS]
    allChannels.clear();
    allChannels.reserve ((size_t) preparedNumChannels);
    for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
        allChannels.push_back (chIdx);
    // [END LS-SUC-PREPARE-BUILD-ALLCHANNELS]

    // [BEGIN LS-SUC-PREPARE-BUILD-MASKS]
    detectChannels.clear();
    applyChannels.clear();
    detectChannels.reserve ((size_t) preparedNumChannels);
    applyChannels.reserve ((size_t) preparedNumChannels);

    for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
    {
        const bool isLfe = (preparedChannelSet.getTypeOfChannel (chIdx) == juce::AudioChannelSet::LFE);

        if (! isLfe)
        {
            detectChannels.push_back (chIdx);
            applyChannels.push_back (chIdx);
        }
    }

    if (detectChannels.empty() || applyChannels.empty())
    {
        detectChannels.clear();
        applyChannels.clear();
        for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
        {
            detectChannels.push_back (chIdx);
            applyChannels.push_back (chIdx);
        }
    }
    // [END LS-SUC-PREPARE-BUILD-MASKS]

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
    globalZoneSmoother.prepare (fs, activeHopSize, 0.5 /*seconds*/);
    smoothedGlobalZoneAmount01 = 1.0f;

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
        globalZoneSmoother.prepare (fs, activeHopSize, 0.5);
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
    for (int i = 0; i < numBands; ++i)
        bandSmoothers[(size_t) i].prepare (fs, activeHopSize, params.attackMs, params.releaseMs);
}

float SpectralUpwardCompressor::computeTargetGainLinear (float bandLevelDb,
                                                        float t0SpectralDb,
                                                        float t1SpectralDb) noexcept
{
    // Active zone fades (soft knees)
    const float lowFade  = softKnee01 (bandLevelDb, t0SpectralDb, params.lowKneeDb);
    const float highFade = 1.0f - softKnee01 (bandLevelDb, t1SpectralDb, params.highKneeDb);
    const float active = lowFade * highFade;

    if (active < 0.001f)
        return 1.0f;

    const float range = std::max (1.0f, t1SpectralDb - t0SpectralDb);
    float pos = (bandLevelDb - t0SpectralDb) / range; // 0..1 across zone
    pos = juce::jlimit (0.0f, 1.0f, pos);

    const float expo = 1.0f + juce::jlimit (0.0f, 1.0f, params.curve) * 3.0f;

    float boostDb = 0.0f;

    if (params.curveType == CurveType::monotonic)
    {
        // quiet end (near T0) => max boost, approaching T1 => 0 boost
        const float shaped = std::pow (std::max (0.0f, 1.0f - pos), expo);
        boostDb = params.maxBoostDb * shaped;
    }
    else // CurveType::bell
    {
        const float d = std::abs (pos - 0.5f) * 2.0f; // 0..1
        const float shaped = std::pow (std::max (0.0f, 1.0f - d), expo);
        boostDb = params.maxBoostDb * shaped;
    }

    boostDb *= active;
    boostDb = juce::jlimit (0.0f, params.maxBoostDb, boostDb);

    const float g = std::pow (10.0f, boostDb / 20.0f);
    return g;
}

void SpectralUpwardCompressor::processFrameAllChannels() noexcept
{
    if (pendingHardReset)
        reset();

    const auto& pr = profiles[(size_t) activeFftChoice];
    const int fftSize = pr.fftSize;
    const int hopSize = pr.hopSize;
    const int maxBin = fftSize / 2;

    // [BEGIN LS-SUC-SELECT-DETECT-APPLY-LISTS]
    const auto& detectList = (params.lfeInDetector ? allChannels : detectChannels);
    const auto& applyList  = (params.lfeInApply    ? allChannels : applyChannels);

    const bool haveDetect = (! detectList.empty());
    const bool haveApply  = (! applyList.empty());
    // [END LS-SUC-SELECT-DETECT-APPLY-LISTS]

    // [BEGIN LS-SUC-LFE-EXCLUDE-BROADBAND-PROXY]
    // 1) Analysis window into per-channel fft buffers.
    // We compute the broadband proxy energy from DETECTOR channels only (default excludes LFE).
    double sumSq = 0.0;

    for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
    {
        auto& st = ch[(size_t) chIdx];

        for (int i = 0; i < fftSize; ++i)
        {
            const float x = st.input[(size_t) i] * pr.window[(size_t) i];
            st.fftBuf[(size_t) i] = x;
        }

        std::fill (st.fftBuf.begin() + (size_t) fftSize,
                   st.fftBuf.begin() + (size_t) (2 * fftSize),
                   0.0f);
    }

    // Sum energy only over detector channels (fallback: all channels)
    const int linkedCh = (haveDetect ? (int) detectList.size() : std::max (1, preparedNumChannels));

    if (haveDetect)
    {
        for (int di = 0; di < (int) detectList.size(); ++di)
        {
            const int chIdx = detectList[(size_t) di];
            if (chIdx < 0 || chIdx >= preparedNumChannels)
            continue;

            const auto& st = ch[(size_t) chIdx];
            for (int i = 0; i < fftSize; ++i)
            {
                const double x = (double) st.fftBuf[(size_t) i];
                sumSq += x * x;
            }
        }
    }
    else
    {
        for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
        {
            const auto& st = ch[(size_t) chIdx];
            for (int i = 0; i < fftSize; ++i)
            {
                const double x = (double) st.fftBuf[(size_t) i];
                sumSq += x * x;
            }
        }
    }

    const double meanSq = sumSq / (double) (fftSize * std::max (1, linkedCh));
    // [END LS-SUC-LFE-EXCLUDE-BROADBAND-PROXY]

    // cheap LUFS-ish proxy (no K-weighting, no gating)
    const double broadbandDb = -0.691 + 10.0 * std::log10 (meanSq + 1.0e-12);

    // [BEGIN LS-SUC-GLOBAL-ZONE-AMOUNT]
    // Option A global zone scaler:
    // - Fade IN around T0 (below T0 => less processing)
    // - Fade OUT around T1 (above T1 => less processing)
    // Uses the same knee widths as the per-band zone for now (simple + parameter-light).
        {
            const float L = (float) broadbandDb;

            const float inAroundT0  = softKnee01 (L, params.t0Lufs, params.lowKneeDb);          // 0..1
            const float outAroundT1 = 1.0f - softKnee01 (L, params.t1Lufs, params.highKneeDb); // 1..0

            const float zoneTarget01 = juce::jlimit (0.0f, 1.0f, inAroundT0 * outAroundT1);
            smoothedGlobalZoneAmount01 = globalZoneSmoother.process (zoneTarget01);
        }
    // [END LS-SUC-GLOBAL-ZONE-AMOUNT]

    // [BEGIN LS-SUC-FORWARD-NONNEGATIVE]
        // 2) Forward FFT per channel
        // IMPORTANT: We rely on JUCE's packed "non-negative frequencies only" real-only format:
        //   bin 0    -> data[0] (real), imag=0
        //   bin N/2  -> data[1] (real), imag=0
        //   bins 1..N/2-1 -> data[2*bin] (real), data[2*bin+1] (imag)
        // So we must request only non-negative frequencies here.
        for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
        {
            auto& st = ch[(size_t) chIdx];
            pr.fft->performRealOnlyForwardTransform (st.fftBuf.data());
        }
    // [END LS-SUC-FORWARD-NONNEGATIVE]

    // [BEGIN LS-SUC-AUDITION-BYPASS-EARLYRETURN]
    // Audition bypass: preserve STFT delay pipeline, but do not measure/apply any gains.
    if (params.auditionBypass)
    {
        // Ensure band smoothers don't carry non-unity state while bypassed.
        for (int bi = 0; bi < numBands; ++bi)
            bandSmoothers[(size_t) bi].reset();

        // Inverse FFT + synthesis (same as normal path)
        for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
        {
            auto& st = ch[(size_t) chIdx];

            pr.fft->performRealOnlyInverseTransform (st.fftBuf.data());

            for (int i = 0; i < fftSize; ++i)
            {
                const float x = st.fftBuf[(size_t) i];
                const float y = (x * pr.window[(size_t) i]);
                st.ola[(size_t) i] += y;
            }

            for (int i = 0; i < hopSize; ++i)
                pushFifo (st, st.ola[(size_t) i] / pr.hopNorm[(size_t) i]);

            const int keep2 = fftSize - hopSize;
            std::memmove (st.ola.data(),
                          st.ola.data() + hopSize,
                          (size_t) keep2 * sizeof (float));
            std::fill (st.ola.begin() + keep2, st.ola.begin() + fftSize, 0.0f);
        }

        // Shift input buffers left by hop (common framing)
        const int keep = fftSize - hopSize;
        for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
        {
            auto& st = ch[(size_t) chIdx];

            std::memmove (st.input.data(),
                          st.input.data() + hopSize,
                          (size_t) keep * sizeof (float));

            std::fill (st.input.begin() + keep, st.input.begin() + fftSize, 0.0f);
        }

        inputWritePos = keep;
        return;
    }
    // [END LS-SUC-AUDITION-BYPASS-EARLYRETURN]

    // 3) Spectral proxy level for adaptive offset
    // Sum power across bins 1..N/2-1 across linked channels
    double pAll = 0.0;
    const int numBins = std::max (1, (maxBin - 1));

    // [BEGIN LS-SUC-LFE-EXCLUDE-SPECTRAL-PROXY]
    if (haveDetect)
    {
        for (int di = 0; di < (int) detectList.size(); ++di)
        {
            const int chIdx = detectList[(size_t) di];
            if (chIdx < 0 || chIdx >= preparedNumChannels)
                continue;

            auto& st = ch[(size_t) chIdx];

            for (int bin = 1; bin <= maxBin - 1; ++bin)
            {
                float re = 0.0f, im = 0.0f;
                getBin (st.fftBuf, fftSize, bin, re, im);
                pAll += (double) re * (double) re + (double) im * (double) im;
            }
        }
    }
    else
    {
        for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
        {
            auto& st = ch[(size_t) chIdx];

            for (int bin = 1; bin <= maxBin - 1; ++bin)
            {
                float re = 0.0f, im = 0.0f;
                getBin (st.fftBuf, fftSize, bin, re, im);
                pAll += (double) re * (double) re + (double) im * (double) im;
            }
        }
    }
    // [END LS-SUC-LFE-EXCLUDE-SPECTRAL-PROXY]

    const double ref = (pr.coherentGain * (double) fftSize * 0.5);
    const double refPower = std::max (1.0e-18, ref * ref);

    const double meanBinPowerAll = pAll / (double) (numBins * linkedCh);
    const double spectralDb = 10.0 * std::log10 (meanBinPowerAll / refPower + 1.0e-18);

    // 4) Adaptive offset: broadband - spectral, smoothed slowly and clamped.
    // Freeze update near silence to avoid wild offsets.
    if (meanSq > 1.0e-12)
    {
        const double instantOffset = broadbandDb - spectralDb;
        const double clamped = juce::jlimit (-30.0, 30.0, instantOffset);
        smoothedOffsetDb = offsetSmoother.process (clamped);
    }

    // [BEGIN LS-SUC-APPLY-CAL-TRIM]
        const double effectiveOffsetDb = smoothedOffsetDb + (double) params.calibrationTrimDb;

        const float t0SpectralDb = (float) ((double) params.t0Lufs - effectiveOffsetDb);
        const float t1SpectralDb = (float) ((double) params.t1Lufs - effectiveOffsetDb);
    // [END LS-SUC-APPLY-CAL-TRIM]

    // [BEGIN LS-SUC-EFFECTIVE-AMOUNT]
        const float userAmount = juce::jlimit (0.0f, 1.0f, params.amount01);
        const float effectiveAmount = juce::jlimit (0.0f, 1.0f, userAmount * smoothedGlobalZoneAmount01);
    // [END LS-SUC-EFFECTIVE-AMOUNT]

    // [BEGIN LS-SUC-PERBAND-GAINS-WITH-FREQ-SMOOTH]
    // 5) Per-band gains computed from LINKED power.
    // First pass: measure band levels and compute per-band target gains.
    // Second pass: cross-band smooth target gains, then time-smooth, then apply to bins.

    // ---- Pass 1: compute target gains per band (linear)
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

        // [BEGIN LS-SUC-LFE-EXCLUDE-PERBAND-MEASURE]
        if (haveDetect)
        {
            for (int di = 0; di < (int) detectList.size(); ++di)
            {
                const int chIdx = detectList[(size_t) di];
                if (chIdx < 0 || chIdx >= preparedNumChannels)
                    continue;

                auto& st = ch[(size_t) chIdx];

                for (int bin = b0; bin <= b1; ++bin)
                {
                    float re = 0.0f, im = 0.0f;
                    getBin (st.fftBuf, fftSize, bin, re, im);
                    pBand += (double) re * (double) re + (double) im * (double) im;
                }
            }
        }
        else
        {
            for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
            {
                auto& st = ch[(size_t) chIdx];

                for (int bin = b0; bin <= b1; ++bin)
                {
                    float re = 0.0f, im = 0.0f;
                    getBin (st.fftBuf, fftSize, bin, re, im);
                    pBand += (double) re * (double) re + (double) im * (double) im;
                }
            }
        }
        // [END LS-SUC-LFE-EXCLUDE-PERBAND-MEASURE]

        const double meanBinPower = pBand / (double) (nBins * linkedCh);
        const float bandDb = (float) (10.0 * std::log10 (meanBinPower / refPower + 1.0e-18));

        float targetG = computeTargetGainLinear (bandDb, t0SpectralDb, t1SpectralDb);

        // Global strength: user amount * global zone amount
        targetG = 1.0f + (targetG - 1.0f) * effectiveAmount;

        bandTargetGain[(size_t) bi] = targetG;
    }

    // ---- Pass 2a: cross-band smoothing (3-tap)
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

    // ---- Pass 2b: time smoothing + apply gains to bins (and mirror bins)
    for (int bi = 0; bi < numBands; ++bi)
    {
        const int b0 = juce::jlimit (1, maxBin - 1, bands[(size_t) bi].startBin);
        const int b1 = juce::jlimit (1, maxBin - 1, bands[(size_t) bi].endBin);
        if (b1 < b0)
            continue;

        const float g = bandSmoothers[(size_t) bi].process (bandTargetGainFreqSmoothed[(size_t) bi]);

        // [BEGIN LS-SUC-LFE-EXCLUDE-PERBAND-APPLY]

        if (haveApply)
        {
            for (int ai = 0; ai < (int) applyList.size(); ++ai)
            {
                const int chIdx = applyList[(size_t) ai];
                if (chIdx < 0 || chIdx >= preparedNumChannels)
                    continue;

                auto& st = ch[(size_t) chIdx];

                for (int bin = b0; bin <= b1; ++bin)
                {
                    scaleBin (st.fftBuf, fftSize, bin, g);

                    const int mirror = fftSize - bin;
                    if (mirror != bin)
                        scaleBin (st.fftBuf, fftSize, mirror, g);
                }
            }
        }
        else
        {
            for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
            {
                auto& st = ch[(size_t) chIdx];

                for (int bin = b0; bin <= b1; ++bin)
                {
                    scaleBin (st.fftBuf, fftSize, bin, g);

                    const int mirror = fftSize - bin;
                    if (mirror != bin)
                        scaleBin (st.fftBuf, fftSize, mirror, g);
                }
            }
        }
        // [END LS-SUC-LFE-EXCLUDE-PERBAND-APPLY]
    }
    // [END LS-SUC-PERBAND-GAINS-WITH-FREQ-SMOOTH]

    // [BEGIN LS-SUC-REMOVE-EXTRA-INV-N]
    // 6) Inverse FFT + synthesis (window + norm + OLA), then emit hop samples
    for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
    {
        auto& st = ch[(size_t) chIdx];

        pr.fft->performRealOnlyInverseTransform (st.fftBuf.data());

        // [BEGIN LS-SUC-SYNTH-ACCUM-THEN-NORMALIZE-ON-EMIT]
            // Accumulate un-normalized WOLA synthesis (normalize at emit time).
            for (int i = 0; i < fftSize; ++i)
            {
                const float x = st.fftBuf[(size_t) i];
                const float y = (x * pr.window[(size_t) i]);
                st.ola[(size_t) i] += y;
            }

            // Emit next hop and normalize using COLA sum for hop positions.
            for (int i = 0; i < hopSize; ++i)
                pushFifo (st, st.ola[(size_t) i] / pr.hopNorm[(size_t) i]);
        // [END LS-SUC-SYNTH-ACCUM-THEN-NORMALIZE-ON-EMIT]

        const int keep = fftSize - hopSize;
        std::memmove (st.ola.data(),
                      st.ola.data() + hopSize,
                      (size_t) keep * sizeof (float));
        std::fill (st.ola.begin() + keep, st.ola.begin() + fftSize, 0.0f);
    }
    // [END LS-SUC-REMOVE-EXTRA-INV-N]

    // 7) Shift input buffers left by hop (common framing)
    const int keep = fftSize - hopSize;
    for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
    {
        auto& st = ch[(size_t) chIdx];

        std::memmove (st.input.data(),
                      st.input.data() + hopSize,
                      (size_t) keep * sizeof (float));

        std::fill (st.input.begin() + keep, st.input.begin() + fftSize, 0.0f);
    }

    inputWritePos = keep;
}

void SpectralUpwardCompressor::process (juce::AudioBuffer<float>& buffer) noexcept
{
    juce::ScopedNoDenormals noDenormals;

    const int numSamples = buffer.getNumSamples();
    const int numChInBuf = buffer.getNumChannels();

    if (numSamples <= 0 || numChInBuf <= 0 || preparedNumChannels <= 0)
        return;

    const int chToProcess = std::min (preparedNumChannels, numChInBuf);

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

        if (inputWritePos >= activeFftSize)
            processFrameAllChannels();

        // Pop output
        for (int chIdx = 0; chIdx < chToProcess; ++chIdx)
        {
            auto& st = ch[(size_t) chIdx];
            buffer.getWritePointer (chIdx)[i] = (st.fifoCount > 0 ? popFifo (st) : 0.0f);
        }
    }
}
} // namespace levelscope::dsp
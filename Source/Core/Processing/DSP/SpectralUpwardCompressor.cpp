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

        // OLA normalization (periodic with hop)
        pr.olaNorm.resize ((size_t) pr.fftSize);
        for (int i = 0; i < pr.fftSize; ++i)
        {
            double s = 0.0;
            for (int m = 0; m < pr.overlapCount; ++m)
            {
                const int idx = i + m * pr.hopSize;
                if (idx >= 0 && idx < pr.fftSize)
                {
                    const double w = (double) pr.window[(size_t) idx];
                    s += w * w;
                }
            }
            pr.olaNorm[(size_t) i] = (float) std::max (1.0e-8, s);
        }
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
    offsetSmoother.prepare (fs, activeHopSize, 4.0 /*seconds*/);
    smoothedOffsetDb = 0.0;

    reset();
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

    offsetSmoother.reset();
    smoothedOffsetDb = 0.0;

    pendingHardReset = false;
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

        offsetSmoother.prepare (fs, activeHopSize, 4.0);
        pendingHardReset = true;
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

    // 1) Analysis window into per-channel fft buffers, and compute broadband proxy energy
    double sumSq = 0.0;

    for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
    {
        auto& st = ch[(size_t) chIdx];

        for (int i = 0; i < fftSize; ++i)
        {
            const float x = st.input[(size_t) i] * pr.window[(size_t) i];
            st.fftBuf[(size_t) i] = x;
            sumSq += (double) x * (double) x;
        }

        std::fill (st.fftBuf.begin() + (size_t) fftSize,
                   st.fftBuf.begin() + (size_t) (2 * fftSize),
                   0.0f);
    }

    const int linkedCh = std::max (1, preparedNumChannels);
    const double meanSq = sumSq / (double) (fftSize * linkedCh);

    // cheap LUFS-ish proxy (no K-weighting, no gating)
    const double broadbandDb = -0.691 + 10.0 * std::log10 (meanSq + 1.0e-12);

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

    // 3) Spectral proxy level for adaptive offset
    // Sum power across bins 1..N/2-1 across linked channels
    double pAll = 0.0;
    const int numBins = std::max (1, (maxBin - 1));

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

    const float t0SpectralDb = (float) (params.t0Lufs - (float) smoothedOffsetDb);
    const float t1SpectralDb = (float) (params.t1Lufs - (float) smoothedOffsetDb);

    const float amount = juce::jlimit (0.0f, 1.0f, params.amount01);

    // 5) Per-band gains computed from LINKED power, applied to all channels
    for (int bi = 0; bi < numBands; ++bi)
    {
        const int b0 = juce::jlimit (1, maxBin - 1, bands[(size_t) bi].startBin);
        const int b1 = juce::jlimit (1, maxBin - 1, bands[(size_t) bi].endBin);
        if (b1 < b0) continue;

        double pBand = 0.0;
        const int nBins = (b1 - b0 + 1);

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

        const double meanBinPower = pBand / (double) (nBins * linkedCh);
        const float bandDb = (float) (10.0 * std::log10 (meanBinPower / refPower + 1.0e-18));

        float targetG = computeTargetGainLinear (bandDb, t0SpectralDb, t1SpectralDb);
        targetG = 1.0f + (targetG - 1.0f) * amount;

        const float g = bandSmoothers[(size_t) bi].process (targetG);

        // Apply gain to each channel's bins in this band
        // [BEGIN LS-SUC-SCALE-MIRROR-BINS]
                for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
                {
                    auto& st = ch[(size_t) chIdx];

                    for (int bin = b0; bin <= b1; ++bin)
                    {
                        scaleBin (st.fftBuf, fftSize, bin, g);

                        const int mirror = fftSize - bin; // conjugate bin for real signals
                        if (mirror != bin)
                            scaleBin (st.fftBuf, fftSize, mirror, g);
                    }
                }
        // [END LS-SUC-SCALE-MIRROR-BINS]
    }

    // [BEGIN LS-SUC-REMOVE-EXTRA-INV-N]
    // 6) Inverse FFT + synthesis (window + norm + OLA), then emit hop samples
    for (int chIdx = 0; chIdx < preparedNumChannels; ++chIdx)
    {
        auto& st = ch[(size_t) chIdx];

        pr.fft->performRealOnlyInverseTransform (st.fftBuf.data());

        for (int i = 0; i < fftSize; ++i)
        {
            // IMPORTANT:
            // Do NOT apply an extra 1/N scale here. In practice this caused output attenuation
            // proportional to FFT size (e.g. ~-60 dB at N=1024, ~-72 dB at N=4096).
            const float x = st.fftBuf[(size_t) i];
            const float y = (x * pr.window[(size_t) i]) / pr.olaNorm[(size_t) i];
            st.ola[(size_t) i] += y;
        }

        for (int i = 0; i < hopSize; ++i)
            pushFifo (st, st.ola[(size_t) i]);

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
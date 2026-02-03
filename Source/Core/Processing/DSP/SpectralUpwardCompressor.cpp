// [BEGIN SUC-CPP]
#include "SpectralUpwardCompressor.h"

namespace levelscope::dsp
{
namespace
{
    inline float dbToLinear (float db) noexcept
    {
        return std::pow (10.0f, db / 20.0f);
    }

    inline float linearToDb (float x) noexcept
    {
        const float v = juce::jmax (1.0e-20f, x);
        return 20.0f * std::log10 (v);
    }
}

void SpectralUpwardCompressor::GainSmoother::prepare (float sampleRate,
                                                     float hopSamples,
                                                     float attackMs,
                                                     float releaseMs) noexcept
{
    const float attackS  = juce::jmax (1.0e-4f, attackMs * 0.001f);
    const float releaseS = juce::jmax (1.0e-4f, releaseMs * 0.001f);

    const float dt = hopSamples / juce::jmax (1.0f, sampleRate);

    aA = std::exp (-dt / attackS);
    aR = std::exp (-dt / releaseS);

    z = 1.0f;
}

float SpectralUpwardCompressor::GainSmoother::process (float target) noexcept
{
    target = juce::jlimit (0.01f, 100.0f, target);
    const float a = (target > z ? aA : aR);
    z = a * z + (1.0f - a) * target;
    return z;
}

//==============================================================================

bool SpectralUpwardCompressor::isSupportedFftSize (int n) noexcept
{
    for (int s : kFftSizes)
        if (n == s) return true;
    return false;
}

int SpectralUpwardCompressor::clampToSupportedFftSize (int n) noexcept
{
    if (isSupportedFftSize (n))
        return n;

    // Pick nearest supported size (simple heuristic).
    int best = kFftSizes[0];
    int bestDist = std::abs (n - best);

    for (int s : kFftSizes)
    {
        const int d = std::abs (n - s);
        if (d < bestDist) { bestDist = d; best = s; }
    }
    return best;
}

float SpectralUpwardCompressor::energyToLufs (float meanSquare) noexcept
{
    const double e = (double) meanSquare;
    if (e <= 0.0)
        return -200.0f;

    const double lufs = -0.691 + 10.0 * std::log10 (e);
    return (float) lufs;
}

void SpectralUpwardCompressor::prepare (double sampleRate,
                                       int numChannels,
                                       const juce::AudioChannelSet& channelSet)
{
    fs = (sampleRate > 0.0 ? sampleRate : 48000.0);
    preparedChannelSet = channelSet;
    numChPrepared = juce::jlimit (1, kMaxChannels, numChannels);

    // Default: detect/apply all channels (Stage D1a is stereo only anyway).
    detectorMask = (numChPrepared >= 32 ? 0xFFFFFFFFu : ((1u << (uint32_t) numChPrepared) - 1u));
    applyMask    = detectorMask;

    // K-weighting filter + channel weights (LFE excluded by default)
    kWeight.prepare (fs, numChPrepared);

    bs1770ChannelWeights.assign ((size_t) numChPrepared, 1.0f);
    for (int ch = 0; ch < numChPrepared; ++ch)
        if (preparedChannelSet.getTypeOfChannel (ch) == juce::AudioChannelSet::LFE)
            bs1770ChannelWeights[(size_t) ch] = 0.0f;

    buildProfiles();

    // Allocate per-channel state (max sized; no allocations in process)
    chState.clear();
    chState.resize ((size_t) numChPrepared);

    for (auto& st : chState)
    {
        st.input.assign  ((size_t) kMaxFftSize, 0.0f);
        st.fftBuf.assign ((size_t) (2 * kMaxFftSize), 0.0f);
        st.ola.assign    ((size_t) kMaxFftSize, 0.0f);
        st.outFifo.assign((size_t) (4 * kMaxFftSize), 0.0f);

        st.fifoRead = st.fifoWrite = st.fifoCount = 0;
    }

    // Apply initial requested FFT size (safe clamp)
    setActiveFftSize (clampToSupportedFftSize (params.fftSize));

    // Build bands + prepare smoothers
    rebuildBands();
    prepareBandSmoothers();

    reset();
}

void SpectralUpwardCompressor::buildProfiles()
{
    profiles.clear();
    profiles.reserve (4);

    for (int size : kFftSizes)
    {
        FftProfile p;
        p.fftSize = size;
        p.hopSize = size / 4;
        p.overlapCount = size / p.hopSize; // for N/4 => 4

        p.fft = std::make_unique<juce::dsp::FFT> ((int) std::log2 ((double) size));

        // sqrt-Hann window (analysis/synthesis)
        p.window.resize ((size_t) size);
        for (int n = 0; n < size; ++n)
        {
            const double hann = 0.5 - 0.5 * std::cos (2.0 * juce::MathConstants<double>::pi
                                                      * (double) n / (double) (size - 1));
            p.window[(size_t) n] = (float) std::sqrt (juce::jmax (0.0, hann));
        }

        // Coherent gain for bin-centered sine reference
        p.coherentGain = 0.0;
        for (float w : p.window)
            p.coherentGain += (double) w;
        p.coherentGain /= (double) size;

        // Overlap-add normalization profile for shift-left OLA
        p.olaNorm.resize ((size_t) size);
        for (int i = 0; i < size; ++i)
        {
            double s = 0.0;
            for (int m = 0; m < p.overlapCount; ++m)
            {
                const int idx = i + m * p.hopSize;
                if (idx >= 0 && idx < size)
                {
                    const double w = (double) p.window[(size_t) idx];
                    s += w * w;
                }
            }
            p.olaNorm[(size_t) i] = (float) juce::jmax (1.0e-8, s);
        }

        profiles.push_back (std::move (p));
    }
}

void SpectralUpwardCompressor::setActiveFftSize (int fftSize) noexcept
{
    activeProfile = nullptr;
    for (auto& p : profiles)
        if (p.fftSize == fftSize)
            activeProfile = &p;

    if (activeProfile == nullptr && ! profiles.empty())
        activeProfile = &profiles.front();

    activeFftSize = (activeProfile != nullptr ? activeProfile->fftSize : 0);
    activeHopSize = (activeProfile != nullptr ? activeProfile->hopSize : 0);

    // Adaptive offset smoothing coefficient (slow)
    const float tauSeconds = 5.0f;
    const float dt = (activeHopSize > 0 ? (float) activeHopSize / (float) juce::jmax (1.0, fs) : 0.02f);
    offsetAlpha = std::exp (-dt / juce::jmax (0.25f, tauSeconds));
}

void SpectralUpwardCompressor::reset() noexcept
{
    clearAllState();

    offsetDb = 0.0f;
    broadbandEnergyAccum = 0.0;
    broadbandSamplesAccum = 0;
    lastBroadbandLufs = -200.0f;

    for (auto& b : bands)
        b.smoother.reset();

    kWeight.reset();
}

void SpectralUpwardCompressor::clearAllState() noexcept
{
    inputWritePos = 0;

    for (auto& st : chState)
    {
        std::fill (st.input.begin(), st.input.end(), 0.0f);
        std::fill (st.fftBuf.begin(), st.fftBuf.end(), 0.0f);
        std::fill (st.ola.begin(), st.ola.end(), 0.0f);
        std::fill (st.outFifo.begin(), st.outFifo.end(), 0.0f);
        st.fifoRead = st.fifoWrite = st.fifoCount = 0;
    }
}

void SpectralUpwardCompressor::setParameters (const Parameters& p) noexcept
{
    // AUDIO-THREAD-ONLY
    pendingParams = p;
    pendingDirty = true;
}

void SpectralUpwardCompressor::setChannelMasks (uint32_t newDetectorMask, uint32_t newApplyMask) noexcept
{
    // AUDIO-THREAD-ONLY
    detectorMask = newDetectorMask;
    applyMask    = newApplyMask;
}

void SpectralUpwardCompressor::rebuildBandsIfNeeded() noexcept
{
    // Rebuild bands + smoother coeffs if mapping-related params changed.
    const int curBpo = juce::jlimit (1, 12, params.bandsPerOctave);
    const float curMin = params.minFreqHz;
    const float curMax = params.maxFreqHz;

    if (curBpo != lastBandsPerOct || curMin != lastMinHz || curMax != lastMaxHz)
    {
        rebuildBands();
        prepareBandSmoothers();
    }

    if (params.attackMs != lastAttackMs || params.releaseMs != lastReleaseMs)
        prepareBandSmoothers();
}

void SpectralUpwardCompressor::rebuildBands() noexcept
{
    numBands = 0;

    const int N = activeFftSize;
    if (N <= 0) return;

    const int maxBin = N / 2;
    const double df = fs / (double) N;

    const int bandsPerOct = juce::jlimit (1, 12, params.bandsPerOctave);

    float fMin = juce::jlimit (10.0f, 24000.0f, params.minFreqHz);
    float fMax = juce::jlimit (10.0f, 24000.0f, params.maxFreqHz);
    if (fMax < fMin) std::swap (fMin, fMax);

    const double step = std::pow (2.0, 1.0 / (double) bandsPerOct);
    const double edge = std::pow (2.0, 0.5 / (double) bandsPerOct);

    double fc = (double) fMin;

    while (fc <= (double) fMax && numBands < kMaxBands)
    {
        const double f0 = fc / edge;
        const double f1 = fc * edge;

        int b0 = (int) std::ceil  (f0 / df);
        int b1 = (int) std::floor (f1 / df);

        b0 = juce::jlimit (1, maxBin - 1, b0);
        b1 = juce::jlimit (1, maxBin - 1, b1);

        if (b1 >= b0)
        {
            bands[(size_t) numBands].startBin = b0;
            bands[(size_t) numBands].endBin   = b1;
            ++numBands;
        }

        fc *= step;
    }

    lastBandsPerOct = bandsPerOct;
    lastMinHz = fMin;
    lastMaxHz = fMax;
}

void SpectralUpwardCompressor::prepareBandSmoothers() noexcept
{
    const float hop = (float) juce::jmax (1, activeHopSize);
    const float sr  = (float) juce::jmax (1.0, fs);

    for (int i = 0; i < numBands; ++i)
        bands[(size_t) i].smoother.prepare (sr, hop, params.attackMs, params.releaseMs);

    lastAttackMs  = params.attackMs;
    lastReleaseMs = params.releaseMs;
}

inline void SpectralUpwardCompressor::pushFifo (ChannelState& st, float s) noexcept
{
    const int cap = (int) st.outFifo.size();
    if (cap <= 0) return;

    if (st.fifoCount >= cap)
    {
        // overflow: drop oldest
        st.fifoRead = (st.fifoRead + 1) % cap;
        --st.fifoCount;
    }

    st.outFifo[(size_t) st.fifoWrite] = s;
    st.fifoWrite = (st.fifoWrite + 1) % cap;
    ++st.fifoCount;
}

inline float SpectralUpwardCompressor::popFifo (ChannelState& st) noexcept
{
    const int cap = (int) st.outFifo.size();
    if (cap <= 0 || st.fifoCount <= 0)
        return 0.0f;

    const float s = st.outFifo[(size_t) st.fifoRead];
    st.fifoRead = (st.fifoRead + 1) % cap;
    --st.fifoCount;
    return s;
}

inline void SpectralUpwardCompressor::scaleBin (float* fftData, int fftSize, int bin, float g) noexcept
{
    // bin 0    -> re=fftData[0], im=0
    // bin N/2  -> re=fftData[1], im=0
    if (bin == 0)          { fftData[0] *= g; return; }
    if (bin == fftSize/2)  { fftData[1] *= g; return; }

    const int i = 2 * bin;
    fftData[(size_t) i]     *= g;
    fftData[(size_t) (i+1)] *= g;
}

float SpectralUpwardCompressor::softKnee01 (float x, float threshold, float kneeWidth) noexcept
{
    kneeWidth = juce::jmax (1.0e-3f, kneeWidth);
    const float d = x - threshold;

    if (d < -kneeWidth) return 0.0f;
    if (d >  kneeWidth) return 1.0f;

    const float t = (d + kneeWidth) / (2.0f * kneeWidth);
    return t * t * (3.0f - 2.0f * t);
}

float SpectralUpwardCompressor::calculateBandGainLinear (float bandLevelDb,
                                                        float t0BandDb,
                                                        float t1BandDb) const noexcept
{
    // Outside zone -> unity (with knees)
    const float lowFade  = softKnee01 (bandLevelDb, t0BandDb, params.lowKneeDb);
    const float highFade = 1.0f - softKnee01 (bandLevelDb, t1BandDb, params.highKneeDb);
    const float activeZone = lowFade * highFade;

    if (activeZone < 0.001f)
        return 1.0f;

    const float range = juce::jmax (1.0f, t1BandDb - t0BandDb);
    float pos = (bandLevelDb - t0BandDb) / range; // 0..1
    pos = juce::jlimit (0.0f, 1.0f, pos);

    // Monotonic: max boost at 

#pragma once

// [BEGIN LS-SUC-HEADER]
// SpectralUpwardCompressor (Stage D1a building block)
// - STFT-based spectral upward processing with per-band floor (T0) and upper bound (T1)
// - Supports LUFS-oriented UI thresholds by translating them into spectral-domain thresholds
//   via an adaptive + smoothed offset (cheap proxy; not full BS.1770).
//
// Threading contract (IMPORTANT):
// - All setters are AUDIO-THREAD-ONLY unless explicitly documented otherwise.
// - prepare()/reset() are non-audio-thread (or at least not concurrently with process()).
// - process() performs no allocations, no locks.
//
// This is a DSP building block; MultiThresholdDynamicsModule decides when/how to enable it.
// [END LS-SUC-HEADER]

#include <juce_audio_basics/juce_audio_basics.h>
#include <juce_dsp/juce_dsp.h>

#include <array>
#include <vector>
#include <memory>
#include <cstring>
#include <cmath>
#include <algorithm>
// [BEGIN LS-SUC-INCLUDE-CSTDINT]
#include <cstdint>
// [END LS-SUC-INCLUDE-CSTDINT]

namespace levelscope::dsp
{
class SpectralUpwardCompressor
{
public:
    enum class CurveType : int
    {
        monotonic = 0, // quieter (above T0) => more boost, approaching T1 => less boost
        bell      = 1  // mid-zone boosted most (legacy option; useful later)
    };

    struct Parameters
    {
        // Thresholds shown on LUFS history curve (UI orientation).
        // Internally translated to spectral domain using adaptive offset.
        float t0Lufs = -45.0f;
        float t1Lufs = -30.0f;

        float amount01   = 1.0f;  // 0..1 fade in/out
        float maxBoostDb = 8.0f;  // clamp

        // shaping
        CurveType curveType = CurveType::monotonic;
        float curve = 0.5f;       // 0..1, shapes the curve exponent

        // knees for the "active zone" (soft transitions near T0 and T1)
        float lowKneeDb  = 3.0f;
        float highKneeDb = 3.0f;

        // smoothing of band gains (per-hop)
        float attackMs  = 5.0f;
        float releaseMs = 50.0f;

        // advanced
        int   fftSizeChoice     = 2;      // 0=1024,1=2048,2=4096,3=8192
        int   bandsPerOctChoice = 2;      // 0=1,1=2,2=3,3=6
        float minFreqHz = 25.0f;
        float maxFreqHz = 20000.0f;

        // [BEGIN LS-SUC-CAL-TRIM-PARAM]
                // User trim (dB) applied on top of adaptive smoothed offset.
                // Purpose: make LUFS-threshold handles line up intuitively with SUC behavior.
                float calibrationTrimDb = 0.0f;
        // [END LS-SUC-CAL-TRIM-PARAM]

        // [BEGIN LS-SUC-LFE-MASK-PARAMS]
        // LFE routing policies:
        // - lfeInDetector: include LFE when computing control signals
        // - lfeInApply: apply computed gain to LFE
        bool lfeInDetector = false;
        bool lfeInApply    = false;
        // [END LS-SUC-LFE-MASK-PARAMS]

        // [BEGIN LS-SUC-AUDITION-BYPASS-PARAM]
        bool auditionBypass = false; // if true: preserve STFT delay but apply no spectral gain changes
        // [END LS-SUC-AUDITION-BYPASS-PARAM]
    };

    SpectralUpwardCompressor() = default;
    ~SpectralUpwardCompressor() = default;

    SpectralUpwardCompressor (const SpectralUpwardCompressor&) = delete;
    SpectralUpwardCompressor& operator= (const SpectralUpwardCompressor&) = delete;

    void prepare (double sampleRate,
                  int numChannels,
                  const juce::AudioChannelSet& channelSet,
                  int maxBlockSize);

    void reset() noexcept;

    // AUDIO-THREAD-ONLY (see header comment)
    void setParametersAudioThread (const Parameters& p) noexcept;
    // [BEGIN LS-SUC-STAGE-E-MASK-API]
    // AUDIO-THREAD-ONLY:
    // If detect/apply bits are both 0, clears override and reverts to legacy LFE policy in params.
    // Bits are channel indices (bit 0 => channel 0). Only first 16 channels are addressable.
    void setChannelMasksAudioThread (uint16_t detectMaskBits,
                                     uint16_t applyMaskBits,
                                     bool unlinked) noexcept;
    // [END LS-SUC-STAGE-E-MASK-API]

    void process (juce::AudioBuffer<float>& buffer) noexcept;

    int getLatencySamples() const noexcept { return activeFftSize; }

    // [BEGIN LS-SUC-UPWARD-METERING-GETTERS]
    // AUDIO THREAD writes, UI thread may read via MTDM snapshot on message thread.
    // Linear gain >= 1.0. (1.0 means no boost)
    float getLastBlockMaxLinearGain()  const noexcept { return lastBlockMaxLinearGain; }
    float getLastBlockLastLinearGain() const noexcept { return lastBlockLastLinearGain; }
    // [END LS-SUC-UPWARD-METERING-GETTERS]

private:
    //==============================================================================
    // Constants / limits
    //==============================================================================
    static constexpr int kNumFftChoices = 4;
    static constexpr int kFftSizes[kNumFftChoices] = { 1024, 2048, 4096, 8192 };

    static constexpr int kNumBandsPerOctChoices = 4;
    static constexpr int kBandsPerOct[kNumBandsPerOctChoices] = { 1, 2, 3, 6 };

    static constexpr int kMaxBands = 128;

    //==============================================================================
    struct GainSmoother
    {
        void prepare (double sampleRate, int hopSamples, float attackMs, float releaseMs) noexcept
        {
            const float attackS  = std::max (1.0e-4f, attackMs  * 0.001f);
            const float releaseS = std::max (1.0e-4f, releaseMs * 0.001f);

            const float dt = (float) hopSamples / std::max (1.0f, (float) sampleRate);

            aA = std::exp (-dt / attackS);
            aR = std::exp (-dt / releaseS);
            z  = 1.0f;
        }

        void reset() noexcept { z = 1.0f; }

        float process (float target) noexcept
        {
            target = juce::jlimit (0.01f, 100.0f, target);
            const float a = (target > z ? aA : aR);
            z = a * z + (1.0f - a) * target;
            return z;
        }

        float aA = 0.9f, aR = 0.99f;
        float z  = 1.0f;
    };

    struct OffsetSmoother
    {
        void prepare (double sampleRate, int hopSamples, double timeSeconds) noexcept
        {
            const double dt = (double) hopSamples / std::max (1.0, sampleRate);
            const double tau = std::max (0.05, timeSeconds);
            a = std::exp (-dt / tau);
            z = 0.0;
            hasValue = false;
        }

        void reset() noexcept
        {
            z = 0.0;
            hasValue = false;
        }

        double process (double x) noexcept
        {
            if (! hasValue)
            {
                z = x;
                hasValue = true;
                return z;
            }

            z = a * z + (1.0 - a) * x;
            return z;
        }

        double a = 0.99;
        double z = 0.0;
        bool hasValue = false;
    };

    // [BEGIN LS-SUC-GLOBAL-AMOUNT-SMOOTHER]
    struct OnePoleSmoother
    {
        void prepare (double sampleRate, int hopSamples, double timeSeconds) noexcept
        {
            const double dt  = (double) hopSamples / std::max (1.0, sampleRate);
            const double tau = std::max (0.01, timeSeconds);
            a = std::exp (-dt / tau);
            z = 0.0f;
            hasValue = false;
        }

        void reset() noexcept
        {
            z = 0.0f;
            hasValue = false;
        }

        float process (float x) noexcept
        {
            x = juce::jlimit (0.0f, 1.0f, x);

            if (! hasValue)
            {
                z = x;
                hasValue = true;
                return z;
            }

            z = (float) (a * (double) z + (1.0 - a) * (double) x);
            return z;
        }

        double a = 0.99;
        float  z = 0.0f;
        bool   hasValue = false;
    };
    // [END LS-SUC-GLOBAL-AMOUNT-SMOOTHER]

    struct Band
    {
        int startBin = 1;
        int endBin   = 1;
    };

    struct ChannelState
    {
        std::vector<float> input;   // maxFftSize
        std::vector<float> fftBuf;  // 2*maxFftSize (JUCE real-only packed)
        std::vector<float> ola;     // maxFftSize

        std::vector<float> fifo;    // ring for output samples
        int fifoRead = 0, fifoWrite = 0, fifoCount = 0;
    };

    //==============================================================================
    // FFT profiles (prebuilt for each allowed size to allow switching without allocations)
    //==============================================================================
    struct FftProfile
    {
        int fftSize = 0;
        int hopSize = 0;
        int overlapCount = 4;

        std::unique_ptr<juce::dsp::FFT> fft;
        // [BEGIN LS-SUC-HOP-NORM-DECL]
            std::vector<float> window;   // sqrt-Hann
            std::vector<float> hopNorm;  // size = hopSize; COLA norm for emitted hop samples
            double coherentGain = 0.5;
        // [END LS-SUC-HOP-NORM-DECL]
    };

    //==============================================================================
    // Internal helpers
    //==============================================================================
    static int clampChoice (int v, int numChoices) noexcept
    {
        return juce::jlimit (0, std::max (0, numChoices - 1), v);
    }

    static float softKnee01 (float levelDb, float threshold, float kneeWidthDb) noexcept;

    void rebuildBandsNoAlloc() noexcept;
    void updateSmoothersNoAlloc() noexcept;

    void processFrameAllChannels() noexcept;

    // [BEGIN LS-SUC-FFT-FULL-COMPLEX-FORMAT]
        // JUCE real-only FFT (full spectrum format):
        // After performRealOnlyForwardTransform (default mode), the buffer contains N complex bins:
        //   bin k: re = data[2*k], im = data[2*k+1], for k = 0..N-1.
        //
        // IMPORTANT: To keep the inverse real-valued, any magnitude scaling applied to bin k in
        // 1..N/2-1 must also be applied to its conjugate mirror bin (N - k).
        inline void getBin (const std::vector<float>& fftBuf, int /*fftSize*/, int bin, float& re, float& im) const noexcept
        {
            const int i = 2 * bin;
            re = fftBuf[(size_t) i];
            im = fftBuf[(size_t) (i + 1)];
        }

        inline void scaleBin (std::vector<float>& fftBuf, int /*fftSize*/, int bin, float g) noexcept
        {
            const int i = 2 * bin;
            fftBuf[(size_t) i]     *= g;
            fftBuf[(size_t) (i+1)] *= g;
        }
    // [END LS-SUC-FFT-FULL-COMPLEX-FORMAT]

    inline void pushFifo (ChannelState& st, float s) noexcept
    {
        const int cap = (int) st.fifo.size();
        if (cap <= 0) return;

        if (st.fifoCount >= cap)
        {
            st.fifoRead = (st.fifoRead + 1) % cap;
            --st.fifoCount;
        }

        st.fifo[(size_t) st.fifoWrite] = s;
        st.fifoWrite = (st.fifoWrite + 1) % cap;
        ++st.fifoCount;
    }

    inline float popFifo (ChannelState& st) noexcept
    {
        const int cap = (int) st.fifo.size();
        if (cap <= 0 || st.fifoCount <= 0) return 0.0f;

        const float s = st.fifo[(size_t) st.fifoRead];
        st.fifoRead = (st.fifoRead + 1) % cap;
        --st.fifoCount;
        return s;
    }

    float computeTargetGainLinear (float bandLevelDb,
                                   float t0SpectralDb,
                                   float t1SpectralDb) noexcept;

    //==============================================================================
    // State
    //==============================================================================
    Parameters params;

    double fs = 48000.0;
    int preparedNumChannels = 0;
    juce::AudioChannelSet preparedChannelSet;
    int preparedMaxBlockSize = 0;

    // common write position (channels are synchronous)
    int inputWritePos = 0;

    // active FFT selection
    int activeFftChoice = 2;
    int activeFftSize = 4096;
    int activeHopSize = 1024;

    int maxFftSize = 8192;

    std::array<FftProfile, kNumFftChoices> profiles;

    // bands (fixed storage, variable count)
    std::array<Band, kMaxBands> bands;
    int numBands = 0;

    // band smoothing (linked for now; later we can add per-channel smoothers)
    std::array<GainSmoother, kMaxBands> bandSmoothers;

    // smoothed LUFS<->spectral translation
    OffsetSmoother offsetSmoother;
    double smoothedOffsetDb = 0.0;

    // [BEGIN LS-SUC-GLOBAL-ZONE-AND-FREQ-SMOOTH-MEMBERS]
    // Global zone scaler (Option A): fades SUC in around T0 and fades out around T1 based on broadband loudness proxy.
    OnePoleSmoother globalZoneSmoother;
    float smoothedGlobalZoneAmount01 = 1.0f; // 0..1

    // Scratch arrays (no allocations in process): per-band target gains before/after cross-band smoothing
    std::array<float, kMaxBands> bandTargetGain {};
    std::array<float, kMaxBands> bandTargetGainFreqSmoothed {};
    // [END LS-SUC-GLOBAL-ZONE-AND-FREQ-SMOOTH-MEMBERS]

    // [BEGIN LS-SUC-UPWARD-METERING-MEMBERS]
    // Upward boost metering (linear gain >= 1). Updated on audio thread.
    float lastBlockMaxLinearGain  = 1.0f;
    float lastBlockLastLinearGain = 1.0f;

    // Updated by processFrameAllChannels(), consumed by process() when a frame is processed.
    float lastFrameMaxLinearGain  = 1.0f;
    float lastFrameLastLinearGain = 1.0f;
    // [END LS-SUC-UPWARD-METERING-MEMBERS]

    // change tracking (no allocations; triggers reset/recompute safely)
    int lastBandsPerOctChoice = 2;
    float lastMinFreqHz = 25.0f;
    float lastMaxFreqHz = 20000.0f;

    float lastAttackMs  = 5.0f;
    float lastReleaseMs = 50.0f;

    bool pendingHardReset = false;

    std::vector<ChannelState> ch;

    // [BEGIN LS-SUC-STAGE-E-MASK-MEMBERS]
    static constexpr int kMaxMaskChannels = 16;

    // Prepared masks (computed in prepare() from channel roles)
    uint16_t preparedAllMaskBits    = 0;
    uint16_t preparedNonLfeMaskBits = 0;
    uint16_t preparedLfeMaskBits    = 0;

    // External override masks (set by MTDM Stage E). If both == 0 => not set.
    uint16_t externalDetectMaskBits = 0;
    uint16_t externalApplyMaskBits  = 0;
    bool     externalUnlinked       = false;
    bool     externalMasksActive    = false;

    // Cached effective masks + fixed index lists (rebuilt when masks change, no alloc)
    uint16_t effectiveDetectMaskBitsCached = 0;
    uint16_t effectiveApplyMaskBitsCached  = 0;

    std::array<uint8_t, kMaxMaskChannels> detectIdx {};
    std::array<uint8_t, kMaxMaskChannels> applyIdx  {};
    int detectCount = 0;
    int applyCount  = 0;

    void rebuildIndexListsIfNeededNoAlloc (uint16_t effectiveDetectBits,
                                          uint16_t effectiveApplyBits) noexcept;

    // Unlinked state: per-channel smoothing/translation
    std::array<OffsetSmoother, kMaxMaskChannels> offsetSmootherUnlinked;
    std::array<double,        kMaxMaskChannels> smoothedOffsetDbUnlinked {};
    std::array<OnePoleSmoother, kMaxMaskChannels> globalZoneSmootherUnlinked;
    std::array<float,          kMaxMaskChannels> smoothedGlobalZoneAmount01Unlinked {};

    std::array<std::array<GainSmoother, kMaxBands>, kMaxMaskChannels> bandSmoothersUnlinked;
    // [END LS-SUC-STAGE-E-MASK-MEMBERS]
};
} // namespace levelscope::dsp
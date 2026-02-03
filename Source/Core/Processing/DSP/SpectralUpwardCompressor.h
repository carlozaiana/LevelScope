#pragma once

// [BEGIN SUC-HEADER]
#include <juce_dsp/juce_dsp.h>
#include <juce_audio_basics/juce_audio_basics.h>

#include "Core/BS1770KWeighting.h" // for a cheap LUFS-ish proxy aligned with existing analysis constant

#include <vector>
#include <array>
#include <memory>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>

//==============================================================================
// SpectralUpwardCompressor (Stage D1a DSP building block)
//
// Threading contract (IMPORTANT):
// - All setters (setParameters/setChannelMasks) are AUDIO-THREAD-ONLY.
// - No locks, no allocations in processBlock().
// - All heap allocations happen in prepare().
//
// Behavior summary:
// - STFT (hop = N/4, 75% overlap), sqrt-Hann analysis/synthesis.
// - Per-band upward gain with a per-band floor (T0) to avoid lifting noise.
// - User thresholds are provided in "LUFS" UI scale.
// - Internally, we derive a per-frame adaptive + smoothed offset that maps
//   broadband LUFS-ish (time-domain proxy) to spectral dB (FFT proxy).
//   Then we translate:
//       T0_bandDb = T0_userLufs - offsetDb
//       T1_bandDb = T1_userLufs - offsetDb
//
// Notes:
// - This is not a full BS.1770 measurement engine (no gating/integration spec).
//   It uses the same LUFS conversion constant (-0.691) and optional K-weighting
//   filter to better align UI thresholds with loudness curves.
// - Designed to be embedded inside MultiThresholdDynamicsModule.
//==============================================================================

namespace levelscope::dsp
{
class SpectralUpwardCompressor
{
public:
    // Supported FFT sizes for Stage D1a (must be powers of two).
    static constexpr int kFftSizes[] = { 1024, 2048, 4096, 8192 };

    struct Parameters
    {
        // User-facing thresholds (on LUFS axis in UI)
        float t0UserLufs = -45.0f;
        float t1UserLufs = -30.0f;

        // Effect
        float amount01   = 1.0f;  // global wet amount 0..1
        float maxBoostDb = 8.0f;

        // Curve shaping (monotonic upward curve; higher -> more concentrated toward low end)
        float curve = 0.5f; // 0..1

        // Soft knees around T0 and T1 (in the *band domain* after translation)
        float lowKneeDb  = 3.0f;
        float highKneeDb = 3.0f;

        // Per-band gain smoothing (per-hop update)
        float attackMs  = 5.0f;
        float releaseMs = 50.0f;

        // Advanced quality / mapping
        int   fftSize        = 4096;   // must be one of kFftSizes
        int   bandsPerOctave = 3;      // typical: 1,2,3,6
        float minFreqHz      = 25.0f;
        float maxFreqHz      = 20000.0f;
    };

    SpectralUpwardCompressor() = default;

    void prepare (double sampleRate,
                  int numChannels,
                  const juce::AudioChannelSet& channelSet);

    void reset() noexcept;

    // AUDIO-THREAD-ONLY setters
    void setParameters (const Parameters& p) noexcept;
    void setChannelMasks (uint32_t detectorMask, uint32_t applyMask) noexcept;

    // In-place processing. Writes output for all channels.
    void processBlock (juce::AudioBuffer<float>& buffer) noexcept;

    int getLatencySamples() const noexcept { return activeFftSize; }

private:
    //==============================================================================
    // Internal helpers
    //==============================================================================

    struct GainSmoother
    {
        void prepare (float sampleRate, float hopSamples, float attackMs, float releaseMs) noexcept;
        void reset() noexcept { z = 1.0f; }

        float process (float target) noexcept;

        float aA = 0.9f, aR = 0.99f;
        float z  = 1.0f;
    };

    struct Band
    {
        int startBin = 1;
        int endBin   = 1;
        GainSmoother smoother;
    };

    struct FftProfile
    {
        int fftSize = 0;
        int hopSize = 0;
        int overlapCount = 4;

        std::unique_ptr<juce::dsp::FFT> fft;

        std::vector<float> window;   // sqrt-Hann
        std::vector<float> olaNorm;  // overlap-add normalization profile
        double coherentGain = 0.5;
    };

    struct ChannelState
    {
        std::vector<float> input;     // size = maxFftSize
        std::vector<float> fftBuf;    // size = 2*maxFftSize (JUCE real-only packing)
        std::vector<float> ola;       // size = maxFftSize
        std::vector<float> outFifo;   // size = 4*maxFftSize

        int fifoRead = 0, fifoWrite = 0, fifoCount = 0;
    };

    static constexpr int kMaxFftSize   = 8192;
    static constexpr int kMaxBands     = 128;
    static constexpr int kMaxChannels  = 16; // future-proof up to 7.1.4 (12) + headroom

    static bool isSupportedFftSize (int n) noexcept;
    static int  clampToSupportedFftSize (int n) noexcept;

    void buildProfiles();
    void setActiveFftSize (int fftSize) noexcept;

    void rebuildBandsIfNeeded() noexcept;
    void rebuildBands() noexcept;
    void prepareBandSmoothers() noexcept;

    void clearAllState() noexcept;

    // Packed FFT helpers (JUCE real-only format)
    inline void scaleBin (float* fftData, int fftSize, int bin, float g) noexcept;

    // Main STFT frame processing (runs once per hop after initial fill)
    void processFrame() noexcept;

    // Gain curve (monotonic upward curve with soft knees)
    float calculateBandGainLinear (float bandLevelDb,
                                  float t0BandDb,
                                  float t1BandDb) const noexcept;

    static float softKnee01 (float x, float threshold, float kneeWidth) noexcept;

    static float energyToLufs (float meanSquare) noexcept;

    // FIFO helpers
    inline void pushFifo (ChannelState& st, float s) noexcept;
    inline float popFifo (ChannelState& st) noexcept;

    //==============================================================================
    // State
    //==============================================================================

    double fs = 48000.0;
    int numChPrepared = 0;
    juce::AudioChannelSet preparedChannelSet;

    // Channel masks (bit i = channel i). Stage D1a: stereo only, but kept for future.
    uint32_t detectorMask = 0xFFFFFFFFu;
    uint32_t applyMask    = 0xFFFFFFFFu;

    // Profiles for supported FFT sizes
    std::vector<FftProfile> profiles;
    FftProfile* activeProfile = nullptr;

    int activeFftSize = 0;
    int activeHopSize = 0;

    // Shared band list + smoothers (linked gains)
    std::array<Band, kMaxBands> bands {};
    int numBands = 0;

    // Cached band-mapping parameters to detect changes
    int   lastBandsPerOct = 0;
    float lastMinHz = 0.0f;
    float lastMaxHz = 0.0f;
    float lastAttackMs = 0.0f;
    float lastReleaseMs = 0.0f;

    // Input write head (shared across channels; channels stay frame-aligned)
    int inputWritePos = 0;

    // Per-channel STFT state
    std::vector<ChannelState> chState;

    // -------------------------------------------------------------------------
    // Parameters (audio-thread-owned via setParameters())
    Parameters params;
    Parameters pendingParams;
    bool pendingDirty = false;

    // -------------------------------------------------------------------------
    // Adaptive offset: broadband LUFS-ish vs spectral dB
    //
    // offsetDb is smoothed slowly (seconds) to provide stable UI-to-spectral mapping.
    float offsetDb = 0.0f;
    float offsetAlpha = 0.999f; // computed from tau + hop dt

    // Broadband LUFS-ish proxy using K-weighting filter + channel weights
    BS1770KWeighting kWeight;
    std::vector<float> bs1770ChannelWeights;

    double broadbandEnergyAccum = 0.0;
    int    broadbandSamplesAccum = 0;

    // Cached last broadband LUFS for debugging/inspection (not currently exposed)
    float lastBroadbandLufs = -200.0f;
//==============================================================================
};
// [END SUC-HEADER]
} // namespace levelscope::dsp

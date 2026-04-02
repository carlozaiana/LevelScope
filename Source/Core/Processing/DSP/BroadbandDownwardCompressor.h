#pragma once

// [BEGIN LS-BDC-HEADER]
// BroadbandDownwardCompressor (Stage D2a building block)
// - Time-domain (broadband) downward compressor controlled by a T2–T3 engagement zone:
//     * below T2: no compression (unity)
//     * between T2 and T3: crossfade from unity -> full compression
//     * above T3: full compression
// - Linked detection across a detector channel list, apply to an apply channel list.
// - Default policy (built in prepare): exclude LFE from detector/apply unless opted in.
// - RT safe: no allocations / no locks in process().
//
// Threading contract:
// - setParametersAudioThread() is AUDIO-THREAD-ONLY.
// - setChannelListsAudioThread() is AUDIO-THREAD-ONLY (but lists are prebuilt in prepare()).
// - prepare()/reset() are non-audio-thread (or not concurrent with process()).
// [END LS-BDC-HEADER]

#include <juce_audio_basics/juce_audio_basics.h>
#include <vector>
#include <cmath>
#include <algorithm>

// [BEGIN LS-BDC-INCLUDE-CSTDINT]
#include <cstdint>
// [END LS-BDC-INCLUDE-CSTDINT]

namespace levelscope::dsp
{
class BroadbandDownwardCompressor
{
public:
    struct Parameters
    {
        bool  enabled = false;

        float t2Lufs = -18.0f;   // zone start (no compression below)
        float t3Lufs =  -6.0f;   // zone full engagement

        float ratio   = 2.0f;    // >= 1.0
        float kneeDb  = 6.0f;    // soft knee around threshold (T2)

        float attackMs  = 10.0f; // detector + gain smoothing (minimal set)
        float releaseMs = 100.0f;

        float makeupDb = 0.0f;   // optional; default 0
    };

    void prepare (double sampleRate,
                  int numChannels,
                  const juce::AudioChannelSet& channelSet,
                  int maxBlockSize);

    void reset() noexcept;

    // AUDIO-THREAD-ONLY
    void setParametersAudioThread (const Parameters& p) noexcept;

    // [BEGIN LS-BDC-STAGE-E-MASK-API]
    // AUDIO-THREAD-ONLY:
    // If detect/apply bits are both 0, clears override and uses legacy (LFE policy / lists).
    void setChannelMasksAudioThread (uint16_t detectMaskBits,
                                     uint16_t applyMaskBits,
                                     bool unlinked) noexcept;
    // [END LS-BDC-STAGE-E-MASK-API]

    // AUDIO-THREAD-ONLY: selects among prebuilt lists (no allocations)
    void setChannelListsAudioThread (const std::vector<int>* detectList,
                                     const std::vector<int>* applyList) noexcept;

    void process (juce::AudioBuffer<float>& buffer) noexcept;

    // [BEGIN LS-BDC-METERING-GETTERS]
    // Non-RT/UI safe to read; written by audio thread once per process() call.
    // These gains are the compressor gain (zone crossfade + smoothing), EXCLUDING makeup gain.
    float getLastBlockMinCompGain() const noexcept   { return lastBlockMinCompGain; }
    float getLastBlockLastCompGain() const noexcept  { return lastBlockLastCompGain; }

    // Current detector loudness/proxy (LUFS-ish) used for downward engagement.
    // Linked/dialog-mask: shared detector at end of block.
    // Unlinked: max current detector loudness across active channels at end of block.
    float getLastBlockDetectorLufs() const noexcept  { return lastBlockDetectorLufs; }
    // [END LS-BDC-METERING-GETTERS]

    // [BEGIN LS-BDC-LFE-POLICY-API]
    // AUDIO-THREAD-ONLY: selects internal prebuilt lists based on LFE policy flags.
    void setLfePolicyAudioThread (bool lfeInDetector, bool lfeInApply) noexcept;
    // [END LS-BDC-LFE-POLICY-API]

private:
    static float dbToLin (float db) noexcept
    {
        return std::pow (10.0f, db / 20.0f);
    }

    static float smoothstep01 (float x) noexcept
    {
        x = juce::jlimit (0.0f, 1.0f, x);
        return x * x * (3.0f - 2.0f * x);
    }

    void updateCoefficientsIfNeeded() noexcept;

    // Returns linear gain <= 1 based on L (LUFS-ish), threshold, ratio, knee.
    float computeCompressionGain (float levelLufs) const noexcept;

    Parameters params;

    double fs = 48000.0;
    int preparedNumChannels = 0;
    juce::AudioChannelSet preparedChannelSet;
    int preparedMaxBlockSize = 0;

    // [BEGIN LS-BDC-STAGE-E-MASK-MEMBERS]
    static constexpr int kMaxMaskChannels = 16;

    uint16_t preparedAllMaskBits    = 0;
    uint16_t preparedNonLfeMaskBits = 0;

    // Legacy selection (what setLfePolicyAudioThread / setChannelListsAudioThread modifies)
    uint16_t legacyDetectMaskBits = 0;
    uint16_t legacyApplyMaskBits  = 0;

    // External override selection (set by MTDM Stage E)
    uint16_t externalDetectMaskBits = 0;
    uint16_t externalApplyMaskBits  = 0;
    bool     externalUnlinked       = false;
    bool     externalMasksActive    = false;

    // Cached effective masks + fixed index lists
    uint16_t effectiveDetectMaskBitsCached = 0;
    uint16_t effectiveApplyMaskBitsCached  = 0;

    std::array<uint8_t, kMaxMaskChannels> detectIdx {};
    std::array<uint8_t, kMaxMaskChannels> applyIdx  {};
    int detectCount = 0;
    int applyCount  = 0;

    void rebuildIndexListsIfNeededNoAlloc (uint16_t effectiveDetectBits,
                                          uint16_t effectiveApplyBits) noexcept;

    // Unlinked per-channel detector/gain state
    std::array<float, kMaxMaskChannels> envMSUnlinked {};
    std::array<float, kMaxMaskChannels> gainZUnlinked {};
    // [END LS-BDC-STAGE-E-MASK-MEMBERS]

    // detector state (mean-square)
    float envMS = 0.0f;

    // gain smoothing state
    float gainZ = 1.0f;

    // smoothing coefficients (per-sample one-pole)
    float aA = 0.99f, aR = 0.999f;

    float lastAttackMs  = -1.0f;
    float lastReleaseMs = -1.0f;

    // [BEGIN LS-BDC-METERING-MEMBERS]
    float lastBlockMinCompGain   = 1.0f;    // minimum compressor gain (<=1) in last process() call
    float lastBlockLastCompGain  = 1.0f;    // last compressor gain in last process() call
    float lastBlockDetectorLufs  = -200.0f; // current detector loudness/proxy at end of last process() call
    // [END LS-BDC-METERING-MEMBERS]
};
} // namespace levelscope::dsp
#pragma once

// [BEGIN LS-BUC-HEADER]
// BroadbandUpwardCompressor (Stage D1c building block)
// - Time-domain (broadband) upward processor with T0–T1 zone behavior (Option A):
//     * fades in around T0
//     * fades out around T1
// - CurveType: monotonic or bell
// - RT safe: no allocations or locks in process()
//
// Threading contract (IMPORTANT):
// - setParametersAudioThread() is AUDIO-THREAD-ONLY (called from MTDM::process()).
// - prepare()/reset() are non-audio-thread (or at least not concurrent with process()).
// [END LS-BUC-HEADER]

#include <juce_audio_basics/juce_audio_basics.h>
#include <cmath>
#include <algorithm>
// [BEGIN LS-BUC-INCLUDE-CSTDINT]
#include <cstdint>
// [END LS-BUC-INCLUDE-CSTDINT]
// [BEGIN LS-BUC-INCLUDE-VECTOR]
#include <vector>
// [END LS-BUC-INCLUDE-VECTOR]

namespace levelscope::dsp
{
class BroadbandUpwardCompressor
{
public:
    enum class CurveType : int
    {
        monotonic = 0,
        bell      = 1
    };

    struct Parameters
    {
        float t0Lufs = -45.0f;
        float t1Lufs = -30.0f;

        float amount01   = 1.0f;  // 0..1
        float maxBoostDb = 8.0f;  // cap

        CurveType curveType = CurveType::monotonic;
        float curve = 0.5f;       // 0..1 (shaping exponent)

        float lowKneeDb  = 3.0f;
        float highKneeDb = 3.0f;

        // Used both for detector smoothing and gain smoothing (minimal parameter set).
        float attackMs  = 10.0f;
        float releaseMs = 100.0f;

        // [BEGIN LS-BUC-LFE-MASK-PARAMS]
        bool lfeInDetector = false;
        bool lfeInApply    = false;
        // [END LS-BUC-LFE-MASK-PARAMS]
    };

    void prepare (double sampleRate,
                  int numChannels,
                  const juce::AudioChannelSet& channelSet,
                  int maxBlockSize);

    void reset() noexcept;

    // AUDIO-THREAD-ONLY
    void setParametersAudioThread (const Parameters& p) noexcept;

    // [BEGIN LS-BUC-STAGE-E-MASK-API]
    // AUDIO-THREAD-ONLY:
    // If detect/apply bits are both 0, clears override and reverts to legacy LFE policy in params.
    // Bits are channel indices (bit 0 => channel 0), up to 16 channels.
    void setChannelMasksAudioThread (uint16_t detectMaskBits,
                                     uint16_t applyMaskBits,
                                     bool unlinked) noexcept;
    // [END LS-BUC-STAGE-E-MASK-API]

    void process (juce::AudioBuffer<float>& buffer) noexcept;

    int getLatencySamples() const noexcept { return 0; }

    // [BEGIN LS-BUC-UPWARD-METERING-GETTERS]
    float getLastBlockMaxLinearGain()  const noexcept { return lastBlockMaxLinearGain; }
    float getLastBlockLastLinearGain() const noexcept { return lastBlockLastLinearGain; }
    // [END LS-BUC-UPWARD-METERING-GETTERS]

private:
    static float softKnee01 (float levelDb, float thresholdDb, float kneeWidthDb) noexcept
    {
        kneeWidthDb = std::max (1.0e-4f, kneeWidthDb);
        const float d = levelDb - thresholdDb;

        if (d < -kneeWidthDb) return 0.0f;
        if (d >  kneeWidthDb) return 1.0f;

        const float t = (d + kneeWidthDb) / (2.0f * kneeWidthDb);
        return t * t * (3.0f - 2.0f * t);
    }

    static float dbToLin (float db) noexcept
    {
        return std::pow (10.0f, db / 20.0f);
    }

    void updateCoefficientsIfNeeded() noexcept;

    Parameters params;

    double fs = 48000.0;
    int preparedNumChannels = 0;
    juce::AudioChannelSet preparedChannelSet;
    int preparedMaxBlockSize = 0;

    // detector (mean-square) state
    float envMS = 0.0f;

    // gain smoothing state
    float gainZ = 1.0f;

    // smoothing coefficients (per-sample one-pole)
    float aDetA = 0.99f, aDetR = 0.999f;
    float aGainA = 0.99f, aGainR = 0.999f;

    float lastAttackMs  = -1.0f;
    float lastReleaseMs = -1.0f;

    // [BEGIN LS-BUC-UPWARD-METERING-MEMBERS]
    float lastBlockMaxLinearGain  = 1.0f;
    float lastBlockLastLinearGain = 1.0f;
    // [END LS-BUC-UPWARD-METERING-MEMBERS]

    // [BEGIN LS-BUC-STAGE-E-MASK-MEMBERS]
    static constexpr int kMaxMaskChannels = 16;

    uint16_t preparedAllMaskBits    = 0;
    uint16_t preparedNonLfeMaskBits = 0;

    uint16_t externalDetectMaskBits = 0;
    uint16_t externalApplyMaskBits  = 0;
    bool     externalUnlinked       = false;
    bool     externalMasksActive    = false;

    uint16_t effectiveDetectMaskBitsCached = 0;
    uint16_t effectiveApplyMaskBitsCached  = 0;

    std::array<uint8_t, kMaxMaskChannels> detectIdx {};
    std::array<uint8_t, kMaxMaskChannels> applyIdx  {};
    int detectCount = 0;
    int applyCount  = 0;

    void rebuildIndexListsIfNeededNoAlloc (uint16_t effectiveDetectBits,
                                          uint16_t effectiveApplyBits) noexcept;

    // Unlinked per-channel state (detector + gain)
    std::array<float, kMaxMaskChannels> envMSUnlinked {};
    std::array<float, kMaxMaskChannels> gainZUnlinked {};
    // [END LS-BUC-STAGE-E-MASK-MEMBERS]

};
} // namespace levelscope::dsp
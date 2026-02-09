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
    };

    void prepare (double sampleRate,
                  int numChannels,
                  const juce::AudioChannelSet& channelSet,
                  int maxBlockSize);

    void reset() noexcept;

    // AUDIO-THREAD-ONLY
    void setParametersAudioThread (const Parameters& p) noexcept;

    void process (juce::AudioBuffer<float>& buffer) noexcept;

    int getLatencySamples() const noexcept { return 0; }

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

    // [BEGIN LS-BUC-LFE-MASK-MEMBERS]
    // Channel masks (prepared once; RT-safe usage in process()).
    // Default policy: exclude LFE from detector and gain application.
    std::vector<int> detectChannels;
    std::vector<int> applyChannels;
    // [END LS-BUC-LFE-MASK-MEMBERS]

    // [BEGIN LS-BUC-LFE-MASK-PARAMS]
    bool lfeInDetector = false;
    bool lfeInApply    = false;
    // [END LS-BUC-LFE-MASK-PARAMS]
};
} // namespace levelscope::dsp
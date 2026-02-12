#pragma once

// [BEGIN LS-LIM-HEADER]
// LookaheadLimiter v0 (Stage D3a building block)
// - Sample-peak lookahead limiter (true-peak later)
// - Adds algorithmic latency = lookaheadSamples when enabled and lookaheadMs > 0
// - RT safe: no allocations/locks in process()
// - Intended as "output safety ceiling" stage after upward+downward.
//
// Threading contract:
// - setParametersAudioThread() is AUDIO-THREAD-ONLY (called from MTDM::process()).
// - prepare()/reset() are non-audio-thread (or not concurrent with process()).
// [END LS-LIM-HEADER]

#include <juce_audio_basics/juce_audio_basics.h>
#include <vector>
#include <cmath>
#include <algorithm>

namespace levelscope::dsp
{
class LookaheadLimiter
{
public:
    struct Parameters
    {
        bool  enabled = false;

        float ceilingDb = -1.0f;     // dBFS, <= 0 typically
        float lookaheadMs = 5.0f;    // 0..50ms typical
        // [BEGIN LS-LIM-PARAMS-TP]
        float attackMs = 0.5f;        // ramp-down length (requires lookahead to be effective)
        float driveDb  = 0.0f;        // pre-limiter gain
        int   oversamplingChoice = 0; // 0=Off, 1=2x, 2=4x (detector only)
        // [END LS-LIM-PARAMS-TP]
        float releaseMs = 100.0f;    // release smoothing
    };

    void prepare (double sampleRate,
                  int numChannels,
                  const juce::AudioChannelSet& channelSet,
                  int maxBlockSize);

    void reset() noexcept;

    // AUDIO-THREAD-ONLY
    void setParametersAudioThread (const Parameters& p) noexcept;

    void process (juce::AudioBuffer<float>& buffer) noexcept;

    int getLatencySamples() const noexcept { return (params.enabled ? lookaheadSamples : 0); }

private:
    static float dbToLin (float db) noexcept
    {
        return std::pow (10.0f, db / 20.0f);
    }

    void updateLookaheadIfNeeded() noexcept;
    void updateReleaseCoeffIfNeeded() noexcept;

    // [BEGIN LS-LIM-HELPERS-TP]
    void updateAttackIfNeeded() noexcept;
    void updateDriveIfNeeded() noexcept;

    float computeLinkedPeakDetector (const float* const* chans, int chToProcess, int sampleIndex) noexcept;
    // [END LS-LIM-HELPERS-TP]

    Parameters params;

    double fs = 48000.0;
    int preparedNumChannels = 0;
    juce::AudioChannelSet preparedChannelSet;
    int preparedMaxBlockSize = 0;

    // Preallocated delay lines
    struct ChannelDelay
    {
        std::vector<float> buf;
    };

    std::vector<ChannelDelay> delay;
    std::vector<float> gainDelay;      // parallel ring for gain scheduling
    std::vector<float> inputScratch;   // per-sample input capture (size = channels)

    int delayCapacity = 0;
    int writePos = 0;

    // current lookahead
    int lookaheadSamples = 0;
    float lastLookaheadMs = -999.0f;

    // gain state + release smoothing
    float gainZ = 1.0f;
    float aRelease = 0.999f;
    float lastReleaseMs = -999.0f;
    // [BEGIN LS-LIM-STATE-TP]
    float aAttack = 0.0f;
    float lastAttackMs = -999.0f;
    int attackSamples = 0;

    float lastDriveDb = -999.0f;
    float driveLin = 1.0f;

    int lastOversamplingChoice = -999;

    std::vector<float> prevDriven; // per-channel previous driven sample (for 2x/4x detector)
    // [END LS-LIM-STATE-TP]

    bool pendingReset = false;

    static constexpr float kMaxLookaheadMs = 50.0f;
};
} // namespace levelscope::dsp
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
// [BEGIN LS-LIM-INCLUDE-DSP]
#include <juce_dsp/juce_dsp.h>
// [END LS-LIM-INCLUDE-DSP]
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

    // [BEGIN LS-LIM-METERING-GETTERS]
    // Non-RT / UI safe to read; written by audio thread once per process() call.
    float getLastBlockMinGain() const noexcept { return lastBlockMinGain; }
    float getLastBlockLastGain() const noexcept { return lastBlockLastGain; }
    // [END LS-LIM-METERING-GETTERS]

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

    // [BEGIN LS-LIM-DETECTOR-API-FIR]
    float computeLinkedPeakSamplePeak (float* const* chans, int chToProcess, int sampleIndex) noexcept;
    // [END LS-LIM-DETECTOR-API-FIR]
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

    // [BEGIN LS-LIM-FIR-OVERSAMPLING-MEMBERS]
    std::unique_ptr<juce::dsp::Oversampling<float>> os2; // 2x FIR
    std::unique_ptr<juce::dsp::Oversampling<float>> os4; // 4x FIR
    juce::dsp::Oversampling<float>* activeOs = nullptr;

    int osFactor = 1;
    int detectorDelaySamples = 0; // FIR group delay at base rate samples

    // [BEGIN LS-LIM-METERING-MEMBERS]
    float lastBlockMinGain  = 1.0f; // minimum gOut used during last process() call
    float lastBlockLastGain = 1.0f; // gOut of last sample in last process() call
    // [END LS-LIM-METERING-MEMBERS]

    juce::AudioBuffer<float> detectorBuffer; // driven input copy for oversampled detector
    // [END LS-LIM-FIR-OVERSAMPLING-MEMBERS]

    bool pendingReset = false;

    static constexpr float kMaxLookaheadMs = 50.0f;
};
} // namespace levelscope::dsp
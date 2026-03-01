#pragma once

#include <JuceHeader.h>
#include <vector>
#include <atomic>
// [BEGIN LS-UI-METERING-INCLUDES]
#include <memory>
// [END LS-UI-METERING-INCLUDES]

#include "Core/LevelScopeHistoryModel.h"
#include "Core/BS1770KWeighting.h"
#include "Core/RunningLoudnessStats.h"
// [BEGIN LS-PROCESSORCORE-INCLUDES]
#include "Core/Processing/ProcessContext.h"
#include "Core/Processing/ProcessorCore.h"
// [END LS-PROCESSORCORE-INCLUDES]

//==============================================================================
// Main audio processor for LevelScope
// Phase 1 refactor: timeline truth + FIFO + persistence moved to Core model.
// Timebase/discontinuity detection remains here for now.
//==============================================================================

// [BEGIN LS-UI-METERING-FWDDECL]
namespace levelscope { class MultiThresholdDynamicsModule; }
// [END LS-UI-METERING-FWDDECL]

// [BEGIN LS-LATENCY-LISTENER-CLASS]
class LevelScopeAudioProcessor  : public juce::AudioProcessor,
                                  private juce::AudioProcessorValueTreeState::Listener,
                                  private juce::Timer
// [END LS-LATENCY-LISTENER-CLASS]
{
public:
    static constexpr double loudnessFrameRate       = 60.0;  // 60 Hz
    static constexpr double momentaryWindowSeconds  = 0.4;   // 400 ms
    static constexpr double shortTermWindowSeconds  = 3.0;   // 3 s

    LevelScopeAudioProcessor();
    ~LevelScopeAudioProcessor() override;

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlockExpected) override;
    void releaseResources() override;

   #if JUCE_PLUGIN_ENABLE_ARA
    bool isARACompatible() const override { return false; }
   #endif

    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
    void processBlock (juce::AudioBuffer<float>&, juce::MidiBuffer&) override;

    //==============================================================================
    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override { return true; }

    //==============================================================================
    const juce::String getName() const override                     { return JucePlugin_Name; }
    bool acceptsMidi() const override                              { return false; }
    bool producesMidi() const override                             { return false; }
    bool isMidiEffect() const override                             { return false; }
    double getTailLengthSeconds() const override                   { return 0.0; }

    //==============================================================================
    int getNumPrograms() override                                  { return 1; }
    int getCurrentProgram() override                               { return 0; }
    void setCurrentProgram (int) override                          {}
    const juce::String getProgramName (int) override               { return {}; }
    void changeProgramName (int, const juce::String&) override     {}

    //==============================================================================
    void getStateInformation (juce::MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;
    // [BEGIN MTDM-APVTS-ACCESSOR]
        juce::AudioProcessorValueTreeState&       getAPVTS() noexcept       { return apvts; }
        const juce::AudioProcessorValueTreeState& getAPVTS() const noexcept { return apvts; }
    // [END MTDM-APVTS-ACCESSOR]

    // [BEGIN LS-LIMITER-METERING-SNAPSHOT-API]
    struct LimiterMeteringSnapshot
    {
        float grDbCurrent   = 0.0f;
        float grDbBlockPeak = 0.0f;
        float grDbHold      = 0.0f;
    };

    LimiterMeteringSnapshot getLimiterMeteringSnapshot() const noexcept;
    // [END LS-LIMITER-METERING-SNAPSHOT-API]

    // [BEGIN LS-UPWARD-METERING-SNAPSHOT-API]
    struct UpwardMeteringSnapshot
    {
        float boostDbCurrent   = 0.0f;
        float boostDbBlockPeak = 0.0f;
        float boostDbHold      = 0.0f;
    };

    UpwardMeteringSnapshot getUpwardMeteringSnapshot() const noexcept;
    // [END LS-UPWARD-METERING-SNAPSHOT-API]

    // [BEGIN LS-DOWNWARD-METERING-SNAPSHOT-API]
    struct DownwardMeteringSnapshot
    {
        float grDbCurrent   = 0.0f;
        float grDbBlockPeak = 0.0f;
        float grDbHold      = 0.0f;
    };

    DownwardMeteringSnapshot getDownwardMeteringSnapshot() const noexcept;
    // [END LS-DOWNWARD-METERING-SNAPSHOT-API]

    // [BEGIN LS-IO-METERING-SNAPSHOT-API]
    struct IOMeteringSnapshot
    {
        float inPeakDbCurrent  = -200.0f;
        float inPeakDbHold     = -200.0f;
        float inRmsDbCurrent   = -200.0f;

        float outPeakDbCurrent = -200.0f;
        float outPeakDbHold    = -200.0f;
        float outRmsDbCurrent  = -200.0f;
    };

    IOMeteringSnapshot getIOMeteringSnapshot() const noexcept;
    // [END LS-IO-METERING-SNAPSHOT-API]

    //==============================================================================
    // GUI access helpers (UI API unchanged; forwarded to core model)

    int readLoudnessFromFifo (float* momentaryDest,
                              float* shortTermDest,
                              juce::int64* frameIndexDest,
                              int* isPlayingDest,
                              int maxNumToRead) noexcept;

    // Extended FIFO read: also returns LRA gate (LUFS) per frame (can be nullptr in the core call).
    
    int readLoudnessFromFifoEx (float* momentaryDest,
                                float* shortTermDest,
                                float* lraGateDest,
                                juce::int64* frameIndexDest,
                                int* isPlayingDest,
                                int maxNumToRead) noexcept;

    double getTimecodeOffsetSeconds() const noexcept { return timecodeOffsetSeconds.load (std::memory_order_relaxed); }

    double getUserTimecodeOffsetSeconds() const noexcept { return historyModel.getUserTimecodeOffsetSeconds(); }
    void   setUserTimecodeOffsetSeconds (double s) noexcept { historyModel.setUserTimecodeOffsetSeconds (s); }

    juce::int64 getLastHostTimeInSamples() const noexcept { return lastHostTimeSamples.load (std::memory_order_relaxed); }
    double      getLastHostTimeInSeconds() const noexcept { return lastHostTimeSeconds.load (std::memory_order_relaxed); }
    bool        hostHasTimeInSamples() const noexcept { return haveHostTimeSamples.load (std::memory_order_relaxed) != 0; }
    bool        hostHasTimeInSeconds() const noexcept { return haveHostTimeSeconds.load (std::memory_order_relaxed) != 0; }

    double getLoudnessFrameRate() const noexcept { return loudnessFrameRate; }
    int getFrameSamples() const noexcept { return frameSamples; }

    bool getDerivedRmsAtFrameIndex (juce::int64 frameIndex,
                                    float& momentaryRms,
                                    float& shortTermRms) const noexcept
    {
        return historyModel.getDerivedLufsAtFrameIndex (frameIndex, momentaryRms, shortTermRms);
    }

    bool getLraGateLufsAtFrameIndex (juce::int64 frameIndex, float& gateLufs) const noexcept
    {
        return historyModel.getLraGateLufsAtFrameIndex (frameIndex, gateLufs);
    }

    juce::int64 getMaxWrittenFrameIndex() const noexcept
    {
        return historyModel.getMaxWrittenFrameIndex();
    }

    float getRunningIntegratedLufs() const noexcept { return runningStats.getIntegratedLufs(); }
    float getRunningLraLu() const noexcept          { return runningStats.getLraLu(); }
    float getRunningLraGateLufs() const noexcept    { return runningStats.getLraGateLufs(); }

    float getRollingLraLu() const noexcept { return runningStats.getRollingLraLu(); }

    int  getRollingLraWindowSeconds() const noexcept { return runningStats.getRollingWindowSeconds(); }
    void setRollingLraWindowSeconds (int s) noexcept { runningStats.setRollingWindowSeconds (s); }

private:
    //==============================================================================
    // Analysis state (still owned by plugin wrapper for Phase 1)

    double currentSampleRate     = 44100.0;
    int    frameSamples          = 0;   // samples between 60 Hz frames
    int    samplesUntilNextFrame = 0;

    // timeline truth (Core)
    LevelScopeHistoryModel historyModel;

    // Per-frame energy accumulator (audio thread)
    double frameEnergyAccum = 0.0;

    // [BS1770] K-weighting filter state (per channel)
    BS1770KWeighting kWeight;

    // [BS1770] Per-channel weights (LFE=0, all others=1). For now mono/stereo only.
    std::vector<float> bs1770ChannelWeights;
    
    // Member RunningLoudnessStats
    RunningLoudnessStats runningStats;

    // [BEGIN MTDM-APVTS-DECL]
        static juce::AudioProcessorValueTreeState::ParameterLayout createParameterLayout();
        juce::AudioProcessorValueTreeState apvts;
    // [END MTDM-APVTS-DECL]

    // [BEGIN LS-PROCESSORCORE-MEMBER]
    // DSP module-chain host (Stage A: empty graph => no-op, no audible change)
    levelscope::ProcessorCore processorCore;
    // [END LS-PROCESSORCORE-MEMBER]

    // [BEGIN LS-MTDM-UI-HANDLE]
    // Non-RT handle for UI/meter polling. Stored atomically on graph rebuild.
    std::shared_ptr<levelscope::MultiThresholdDynamicsModule> mtdmForUI;
    // [END LS-MTDM-UI-HANDLE]

    // [FIX-START-RAMP] guard against host start ramp overwriting existing truth
    int transportStartOverwriteGuardFrames = 0;

    // [TIMECODE-OFFSET] display-only offset so our ruler can match DAW timecode
    std::atomic<double> timecodeOffsetSeconds { 0.0 };

    // [TIMECODE-DEBUG] last raw host values
    std::atomic<juce::int64> lastHostTimeSamples { 0 };
    std::atomic<double>      lastHostTimeSeconds { 0.0 };
    std::atomic<int>         haveHostTimeSamples { 0 };
    std::atomic<int>         haveHostTimeSeconds { 0 };

    // [BEGIN LS-IO-METERING-ATOMICS]
    // Input/Output metering (max across channels). Updated on audio thread, read on UI thread.
    std::atomic<float> ioInPeakDbCurrent  { -200.0f };
    std::atomic<float> ioInPeakDbHold     { -200.0f };
    std::atomic<float> ioInRmsDbCurrent   { -200.0f };

    std::atomic<float> ioOutPeakDbCurrent { -200.0f };
    std::atomic<float> ioOutPeakDbHold    { -200.0f };
    std::atomic<float> ioOutRmsDbCurrent  { -200.0f };

    // Peak-hold decay rate (dB per second)
    float ioPeakHoldDecayDbPerSec = 12.0f;
    // [END LS-IO-METERING-ATOMICS]

    // [FIX-RESTART-PARTIAL-FRAME]
    juce::int64 lastBlockEndProjectSample = 0;
    int         lastBlockIsPlaying = 0;
    bool        haveLastBlockEnd = false;

    bool        skipNextPartialFrameWrite = false;

    juce::int64 currentBlockStartProjectSample = 0;
    int         currentBlockIsPlaying         = 1;
    int         currentBlockSampleIndex       = 0;

    // Helpers (negative-safe)
    static juce::int64 floorDivInt64 (juce::int64 a, juce::int64 b) noexcept;

    void resetLoudnessState() noexcept;

    void processSampleForLoudness (const float* const* channelData,
                                   int numChannels,
                                   int sampleIndex) noexcept;

    // [BEGIN LS-C2-MODGRAPH-DECL]
        void rebuildModuleGraphFromState (const juce::MemoryBlock* modgChunkData);

        bool processorCorePrepared = false;
        int  lastMaxBlockSizeForSpec = 0;
    // [END LS-C2-MODGRAPH-DECL]

    // [BEGIN LS-LATENCY-HELPER-DECL]
        void updateLatencyFromAPVTS_NonRT();
    // [END LS-LATENCY-HELPER-DECL]

    // [BEGIN LS-LATENCY-LISTENER-DECL]
    void parameterChanged (const juce::String& parameterID, float newValue) override;

    void timerCallback() override;

    std::atomic<bool> latencyDirty { false };

    void registerLatencyParamListeners();
    void unregisterLatencyParamListeners();
    // [END LS-LATENCY-LISTENER-DECL]

    // [BEGIN LS-LATENCY-CACHED-PARAMS]
    // Cached raw parameter pointers for latency computation (no APVTS lookup in getLatencySamples()).
    std::atomic<float>* pLatMtdmEnabled01          = nullptr;
    std::atomic<float>* pLatUpwardModeChoice       = nullptr;
    std::atomic<float>* pLatSucFftSizeChoice       = nullptr;

    std::atomic<float>* pLatLimEnabled01           = nullptr;
    std::atomic<float>* pLatLimLookaheadMs         = nullptr;
    std::atomic<float>* pLatLimOversamplingChoice  = nullptr;

    int computeTotalLatencySamplesFromCachedParams() const noexcept;
    void cacheLatencyParamPointers();
// [END LS-LATENCY-CACHED-PARAMS]

    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (LevelScopeAudioProcessor)
};
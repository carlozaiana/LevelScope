#pragma once

#include <JuceHeader.h>
#include <vector>
#include <atomic>
#include <memory>

//==============================================================================
// Main audio processor for LevelScope (prototype loudness display)
//==============================================================================

class LevelScopeAudioProcessor : public juce::AudioProcessor
{
public:
    // Loudness "frame" rate for visualization (frames per second)
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

    //==============================================================================
    // GUI access helpers

    // Read up to maxNumToRead loudness frames from the FIFO into dest arrays.
    // Returns the number actually read (non-blocking).
    // Each frame has both momentary and short-term RMS values (linear).
    // Read up to maxNumToRead loudness frames from the FIFO into dest arrays.
    // Returns the number actually read (non-blocking).
    int readLoudnessFromFifo (float* momentaryDest,
                              float* shortTermDest,
                              juce::int64* frameIndexDest,
                              int* isPlayingDest,
                              int maxNumToRead) noexcept;

    // Visual loudness frame rate accessor (used by GUI)
    double getLoudnessFrameRate() const noexcept { return loudnessFrameRate; }
    int getFrameSamples() const noexcept { return frameSamples; } // samples between 60 Hz loudness frames

    juce::int64 getMaxWrittenFrameIndex() const noexcept { return maxWrittenFrameIndex.load(); }

private:
    //==============================================================================
    // Loudness / envelope aggregation (momentary + short-term)

    double currentSampleRate      = 44100.0;
    int    momentaryWindowSamples = 0;   // samples in 400 ms window
    int    shortTermWindowSamples = 0;   // samples in 3 s window

    int    frameSamples           = 0;   // samples between loudness frames (for 60 Hz)
    int    samplesUntilNextFrame  = 0;   // countdown to next frame

    // Sliding energy windows for momentary and short-term (per-sample energy)
    std::vector<double> momentaryEnergyBuffer;
    std::vector<double> shortTermEnergyBuffer;

    int    momentaryIndex   = 0;
    int    shortTermIndex   = 0;
    double momentarySum     = 0.0;
    double shortTermSum     = 0.0;

    juce::int64 totalSamplesProcessed = 0; // for startup warm-up

    //==============================================================================
    // [TIMELINE-ENERGY] timeline-truth storage at 60 Hz
    // Store mean-square energy per 60 Hz frame, keyed by absolute frameIndex.
    // Derived momentary/short-term are computed from this stored energy.
    //==============================================================================

    static constexpr double historyLengthSeconds = 3.0 * 3600.0;
    static constexpr int historyCapacityFrames =
        (int) (historyLengthSeconds * loudnessFrameRate + 0.5); // ~648000

    // Published-by-tag ring (single writer: audio thread; readers: GUI thread for snapshots later)
    std::vector<float> energyMeanSquare;   // base measure
    std::vector<float> momentaryRmsHist;   // derived
    std::vector<float> shortTermRmsHist;   // derived
    std::unique_ptr<std::atomic<juce::int64>[]> frameIndexTag; // -1 = empty, else abs frameIndex

    std::atomic<juce::int64> maxWrittenFrameIndex { std::numeric_limits<juce::int64>::min() };

    // Per-frame energy accumulator (audio thread)
    double frameEnergyAccum = 0.0;

    // Helpers
    static juce::int64 floorDivInt64 (juce::int64 a, juce::int64 b) noexcept; // b>0
    static int wrapSlot (juce::int64 absIndex, int capacity) noexcept;

    bool readEnergyAbs (juce::int64 absFrameIndex, float& outEnergy) const noexcept;
    void writeFrameAbs (juce::int64 absFrameIndex, float energyMS, float mRms, float sRms) noexcept;

    float computeWindowRmsFromEnergy (juce::int64 endFrameIndex, int windowFrames) const noexcept;

    juce::int64 lastFrameProjectSample = 0; // [TIMEBASE-PLAYHEAD]
    int         lastFrameIsPlaying     = 1; // [TIMEBASE-PLAYHEAD]

    juce::int64 currentBlockStartProjectSample = 0; // [TIMEBASE-PLAYHEAD]
    int         currentBlockIsPlaying         = 1; // [TIMEBASE-PLAYHEAD]
    int         currentBlockSampleIndex       = 0; // [TIMEBASE-PLAYHEAD] updated per-sample

    struct LoudnessFrame
    {
        float momentaryRms = 0.0f;
        float shortTermRms = 0.0f;

        // [TIMELINE-ENERGY] absolute 60 Hz timeline frame index (can be negative)
        juce::int64 frameIndex = 0;

        // [TIMEBASE-PLAYHEAD]
        int isPlaying = 1;
    };

    static constexpr int loudnessFifoSize = 4096;
    juce::AbstractFifo              loudnessFifo;
    std::vector<LoudnessFrame>      loudnessBuffer;

    void resetLoudnessState() noexcept;

    void processSampleForLoudness (const float* const* channelData,
                                   int numChannels,
                                   int sampleIndex) noexcept;

    void pushLoudnessFrame (float momentaryRms,
                            float shortTermRms) noexcept;

    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (LevelScopeAudioProcessor)
};
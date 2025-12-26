#pragma once

#include <JuceHeader.h>
#include <vector>

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
    int readLoudnessFromFifo (float* momentaryDest,
                              float* shortTermDest,
                              int maxNumToRead) noexcept;

    // Visual loudness frame rate accessor (used by GUI)
    double getLoudnessFrameRate() const noexcept { return loudnessFrameRate; }

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

    struct LoudnessFrame
    {
        float momentaryRms = 0.0f;
        float shortTermRms = 0.0f;
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
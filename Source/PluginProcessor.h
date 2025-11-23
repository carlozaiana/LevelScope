#pragma once

#include <JuceHeader.h>
#include <vector>

//==============================================================================
// Main audio processor for LevelScope
//==============================================================================

class LevelScopeAudioProcessor : public juce::AudioProcessor
{
public:
    // Visual envelope sampling rate (envelope "frames" per second)
    static constexpr double visualSampleRate = 2000.0; // Hz

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

    // Read up to maxNumToRead envelope samples from the FIFO into dest arrays.
    // Returns the number actually read (non-blocking).
    // Each sample has both RMS and Peak values.
    int readEnvelopeFromFifo (float* rmsDest,
                              float* peakDest,
                              int maxNumToRead) noexcept;

    // Visual sample rate accessor (used by GUI)
    double getVisualSampleRate() const noexcept { return visualSampleRate; }

private:
    //==============================================================================
    // Envelope aggregation for visualization (RMS + Peak per window)

    double currentSampleRate = 44100.0;
    int    windowSamples     = 0;     // how many audio frames per envelope chunk (RMS+Peak)

    double rmsSumSquares        = 0.0; // accumulated mean-square energy
    int    framesAccumulated    = 0;   // frames so far in current window
    float  peakAccumulator      = 0.0f; // max abs across all channels and samples in window

    struct EnvelopeSample
    {
        float rms  = 0.0f;
        float peak = 0.0f;
    };

    static constexpr int envelopeFifoSize = 16384; // plenty for UI latency
    juce::AbstractFifo              envelopeFifo;
    std::vector<EnvelopeSample>     envelopeBuffer;

    void resetEnvelopeCollector() noexcept;
    void pushEnvelopeIntoFifo (float rms, float peak) noexcept;

    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (LevelScopeAudioProcessor)
};
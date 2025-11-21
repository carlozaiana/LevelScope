#pragma once

#include <JuceHeader.h>
#include <vector>

//==============================================================================
// Main audio processor for LevelScope
//==============================================================================

class LevelScopeAudioProcessor : public juce::AudioProcessor
{
public:
    // Visual RMS sampling rate (RMS "frames" per second)
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
    // Read up to maxNumToRead RMS values from the FIFO into dest.
    // Returns the number actually read (non-blocking).
    int readRmsFromFifo (float* dest, int maxNumToRead) noexcept;

    // Visual sample rate accessor (used by GUI)
    double getVisualSampleRate() const noexcept { return visualSampleRate; }

private:
    //==============================================================================
    // RMS aggregation for visualization
    double currentSampleRate = 44100.0;
    int rmsWindowSamples = 0;          // how many audio frames per RMS chunk
    double rmsSumSquares = 0.0;        // accumulated mean-square energy
    int rmsFramesAccumulated = 0;      // frames so far in current RMS window

    static constexpr int rmsFifoSize = 16384; // plenty for UI latency
    juce::AbstractFifo rmsFifo;
    std::vector<float> rmsBuffer;

    void resetRmsCollector() noexcept;
    void pushRmsIntoFifo (float rms) noexcept;

    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (LevelScopeAudioProcessor)
};
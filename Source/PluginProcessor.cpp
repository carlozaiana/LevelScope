#include "PluginProcessor.h"
#include "PluginEditor.h"

#include <cmath>
#include <algorithm>

//==============================================================================

LevelScopeAudioProcessor::LevelScopeAudioProcessor()
    : AudioProcessor (BusesProperties()
                        .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
                        .withOutput ("Output", juce::AudioChannelSet::stereo(), true)),
      rmsFifo (rmsFifoSize),
      rmsBuffer ((size_t) rmsFifoSize, 0.0f)
{
}

LevelScopeAudioProcessor::~LevelScopeAudioProcessor() = default;

//==============================================================================

void LevelScopeAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlockExpected)
{
    juce::ignoreUnused (samplesPerBlockExpected);

    currentSampleRate = (sampleRate > 0.0 ? sampleRate : 44100.0);

    // Number of audio frames per visual RMS sample
    rmsWindowSamples = juce::jmax (1, (int) std::round (currentSampleRate / visualSampleRate));

    resetRmsCollector();
}

void LevelScopeAudioProcessor::releaseResources()
{
    resetRmsCollector();
}

void LevelScopeAudioProcessor::resetRmsCollector() noexcept
{
    rmsSumSquares = 0.0;
    rmsFramesAccumulated = 0;
    rmsFifo.reset();
}

//==============================================================================

bool LevelScopeAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
    // Only support mono or stereo on main bus
    const auto& mainIn  = layouts.getMainInputChannelSet();
    const auto& mainOut = layouts.getMainOutputChannelSet();

    if (mainIn.isDisabled() || mainOut.isDisabled())
        return false;

    if (mainIn.size() > 2 || mainOut.size() > 2)
        return false;

    // Require same channel count in/out (typical for metering effect)
    if (mainIn != mainOut)
        return false;

    return true;
}

//==============================================================================

void LevelScopeAudioProcessor::processBlock (juce::AudioBuffer<float>& buffer,
                                             juce::MidiBuffer& midiMessages)
{
    juce::ignoreUnused (midiMessages);

    const int numChannels = getTotalNumInputChannels();
    const int numSamples  = buffer.getNumSamples();

    // Ensure any output channels beyond inputs are silent (JUCE convention)
    for (int ch = numChannels; ch < getTotalNumOutputChannels(); ++ch)
        buffer.clear (ch, 0, numSamples);

    if (numChannels <= 0 || numSamples <= 0)
        return;

    // Prepare read pointers for faster inner loop
    juce::HeapBlock<const float*> channelData ((size_t) numChannels);
    for (int ch = 0; ch < numChannels; ++ch)
        channelData[ch] = buffer.getReadPointer (ch);

    for (int i = 0; i < numSamples; ++i)
    {
        // Average energy across all input channels for this frame
        double energy = 0.0;
        for (int ch = 0; ch < numChannels; ++ch)
        {
            const float s = channelData[ch][i];
            energy += (double) s * (double) s;
        }

        energy /= (double) numChannels;

        rmsSumSquares += energy;
        ++rmsFramesAccumulated;

        if (rmsFramesAccumulated >= rmsWindowSamples)
        {
            const double meanSq = rmsSumSquares / (double) rmsFramesAccumulated;
            const float rms     = (float) std::sqrt (meanSq);

            pushRmsIntoFifo (rms);

            rmsSumSquares       = 0.0;
            rmsFramesAccumulated = 0;
        }
    }

    // Audio is passed through unchanged (we don't touch buffer data)
}

//==============================================================================

void LevelScopeAudioProcessor::pushRmsIntoFifo (float rms) noexcept
{
    int start1 = 0, size1 = 0, start2 = 0, size2 = 0;
    rmsFifo.prepareToWrite (1, start1, size1, start2, size2);

    const int writable = size1 + size2;
    if (writable > 0)
    {
        if (size1 > 0)
            rmsBuffer[(size_t) start1] = rms;
        else if (size2 > 0)
            rmsBuffer[(size_t) start2] = rms;

        rmsFifo.finishedWrite (1);
    }
    else
    {
        // No space: drop sample (rare with large FIFO and low visual rate)
        rmsFifo.finishedWrite (0);
    }
}

int LevelScopeAudioProcessor::readRmsFromFifo (float* dest, int maxNumToRead) noexcept
{
    if (dest == nullptr || maxNumToRead <= 0)
        return 0;

    int start1 = 0, size1 = 0, start2 = 0, size2 = 0;
    rmsFifo.prepareToRead (maxNumToRead, start1, size1, start2, size2);

    const int totalToRead = size1 + size2;

    if (size1 > 0)
        std::copy (rmsBuffer.data() + start1,
                   rmsBuffer.data() + start1 + size1,
                   dest);

    if (size2 > 0)
        std::copy (rmsBuffer.data() + start2,
                   rmsBuffer.data() + start2 + size2,
                   dest + size1);

    rmsFifo.finishedRead (totalToRead);

    return totalToRead;
}

//==============================================================================

void LevelScopeAudioProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    juce::ignoreUnused (destData);
    // No parameters yet; could store zoom state here later if desired.
}

void LevelScopeAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    juce::ignoreUnused (data, sizeInBytes);
    // No parameters yet.
}

//==============================================================================

juce::AudioProcessorEditor* LevelScopeAudioProcessor::createEditor()
{
    return new LevelScopeAudioProcessorEditor (*this);
}

//==============================================================================

// Factory function required by JUCE
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new LevelScopeAudioProcessor();
}
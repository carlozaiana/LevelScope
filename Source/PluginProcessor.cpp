#include "PluginProcessor.h"
#include "PluginEditor.h"

#include <cmath>
#include <algorithm>

//==============================================================================

LevelScopeAudioProcessor::LevelScopeAudioProcessor()
    : AudioProcessor (BusesProperties()
                        .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
                        .withOutput ("Output", juce::AudioChannelSet::stereo(), true)),
      envelopeFifo (envelopeFifoSize),
      envelopeBuffer ((size_t) envelopeFifoSize)
{
}

LevelScopeAudioProcessor::~LevelScopeAudioProcessor() = default;

//==============================================================================

void LevelScopeAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlockExpected)
{
    juce::ignoreUnused (samplesPerBlockExpected);

    currentSampleRate = (sampleRate > 0.0 ? sampleRate : 44100.0);

    // Number of audio frames per visual envelope sample
    windowSamples = juce::jmax (1, (int) std::round (currentSampleRate / visualSampleRate));

    resetEnvelopeCollector();
}

void LevelScopeAudioProcessor::releaseResources()
{
    resetEnvelopeCollector();
}

void LevelScopeAudioProcessor::resetEnvelopeCollector() noexcept
{
    rmsSumSquares     = 0.0;
    framesAccumulated = 0;
    peakAccumulator   = 0.0f;

    envelopeFifo.reset();
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
        // Average energy across all input channels for this frame (for RMS)
        double energy = 0.0;

        // For peak: track max absolute across channels at this frame
        float framePeakAbs = 0.0f;

        for (int ch = 0; ch < numChannels; ++ch)
        {
            const float s = channelData[ch][i];
            energy += (double) s * (double) s;

            const float absS = std::abs (s);
            if (absS > framePeakAbs)
                framePeakAbs = absS;
        }

        energy /= (double) numChannels;

        rmsSumSquares += energy;
        ++framesAccumulated;

        if (framePeakAbs > peakAccumulator)
            peakAccumulator = framePeakAbs;

        if (framesAccumulated >= windowSamples)
        {
            const double meanSq = rmsSumSquares / (double) framesAccumulated;
            const float  rms    = (float) std::sqrt (meanSq);
            const float  peak   = peakAccumulator;

            pushEnvelopeIntoFifo (rms, peak);

            rmsSumSquares     = 0.0;
            framesAccumulated = 0;
            peakAccumulator   = 0.0f;
        }
    }

    // Audio is passed through unchanged (we don't touch buffer data)
}

//==============================================================================

void LevelScopeAudioProcessor::pushEnvelopeIntoFifo (float rms, float peak) noexcept
{
    int start1 = 0, size1 = 0, start2 = 0, size2 = 0;
    envelopeFifo.prepareToWrite (1, start1, size1, start2, size2);

    const int writable = size1 + size2;
    if (writable > 0)
    {
        // We only ever write in one contiguous block for count=1
        const int index = (size1 > 0 ? start1 : start2);

        envelopeBuffer[(size_t) index].rms  = rms;
        envelopeBuffer[(size_t) index].peak = peak;

        envelopeFifo.finishedWrite (1);
    }
    else
    {
        // No space: drop sample (rare with large FIFO and low visual rate)
        envelopeFifo.finishedWrite (0);
    }
}

int LevelScopeAudioProcessor::readEnvelopeFromFifo (float* rmsDest,
                                                    float* peakDest,
                                                    int maxNumToRead) noexcept
{
    if (rmsDest == nullptr || peakDest == nullptr || maxNumToRead <= 0)
        return 0;

    int start1 = 0, size1 = 0, start2 = 0, size2 = 0;
    envelopeFifo.prepareToRead (maxNumToRead, start1, size1, start2, size2);

    const int totalToRead = size1 + size2;

    int destIndex = 0;

    if (size1 > 0)
    {
        for (int i = 0; i < size1; ++i)
        {
            const auto& s = envelopeBuffer[(size_t) (start1 + i)];
            rmsDest [destIndex] = s.rms;
            peakDest[destIndex] = s.peak;
            ++destIndex;
        }
    }

    if (size2 > 0)
    {
        for (int i = 0; i < size2; ++i)
        {
            const auto& s = envelopeBuffer[(size_t) (start2 + i)];
            rmsDest [destIndex] = s.rms;
            peakDest[destIndex] = s.peak;
            ++destIndex;
        }
    }

    envelopeFifo.finishedRead (totalToRead);

    return totalToRead;
}

//==============================================================================

void LevelScopeAudioProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    juce::ignoreUnused (destData);
    // No parameters yet; could store zoom/mode state here later if desired.
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
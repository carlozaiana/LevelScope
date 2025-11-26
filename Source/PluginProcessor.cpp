#include "PluginProcessor.h"
#include "PluginEditor.h"

#include <cmath>
#include <algorithm>

//==============================================================================

LevelScopeAudioProcessor::LevelScopeAudioProcessor()
    : AudioProcessor (BusesProperties()
                        .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
                        .withOutput ("Output", juce::AudioChannelSet::stereo(), true)),
      loudnessFifo (loudnessFifoSize),
      loudnessBuffer ((size_t) loudnessFifoSize)
{
}

LevelScopeAudioProcessor::~LevelScopeAudioProcessor() = default;

//==============================================================================

void LevelScopeAudioProcessor::prepareToPlay (double sampleRate,
                                              int /*samplesPerBlockExpected*/)
{
    currentSampleRate = (sampleRate > 0.0 ? sampleRate : 44100.0);

    momentaryWindowSamples = juce::jmax (1,
        (int) std::round (momentaryWindowSeconds * currentSampleRate));

    shortTermWindowSamples = juce::jmax (1,
        (int) std::round (shortTermWindowSeconds * currentSampleRate));

    frameSamples = juce::jmax (1,
        (int) std::round (currentSampleRate / loudnessFrameRate));

    momentaryEnergyBuffer.assign ((size_t) momentaryWindowSamples, 0.0);
    shortTermEnergyBuffer.assign ((size_t) shortTermWindowSamples, 0.0);

    resetLoudnessState();
}

void LevelScopeAudioProcessor::releaseResources()
{
    resetLoudnessState();
}

void LevelScopeAudioProcessor::resetLoudnessState() noexcept
{
    momentaryIndex   = 0;
    shortTermIndex   = 0;
    momentarySum     = 0.0;
    shortTermSum     = 0.0;
    totalSamplesProcessed = 0;

    samplesUntilNextFrame = frameSamples;

    std::fill (momentaryEnergyBuffer.begin(), momentaryEnergyBuffer.end(), 0.0);
    std::fill (shortTermEnergyBuffer.begin(), shortTermEnergyBuffer.end(), 0.0);

    loudnessFifo.reset();
}

//==============================================================================

bool LevelScopeAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
    // Only support mono or stereo on main bus for this prototype
    const auto& mainIn  = layouts.getMainInputChannelSet();
    const auto& mainOut = layouts.getMainOutputChannelSet();

    if (mainIn.isDisabled() || mainOut.isDisabled())
        return false;

    if (mainIn.size() > 2 || mainOut.size() > 2)
        return false;

    if (mainIn != mainOut)
        return false;

    return true;
}

//==============================================================================

void LevelScopeAudioProcessor::processSampleForLoudness (const float* const* channelData,
                                                         int numChannels,
                                                         int sampleIndex) noexcept
{
    // Simple energy across channels for now (no K-weighting yet).
    double energy = 0.0;

    for (int ch = 0; ch < numChannels; ++ch)
    {
        const float s = channelData[ch][sampleIndex];
        energy += (double) s * (double) s;
    }

    if (numChannels > 0)
        energy /= (double) numChannels;

    // Update momentary sliding window
    if (momentaryWindowSamples > 0)
    {
        const double old = momentaryEnergyBuffer[(size_t) momentaryIndex];
        momentarySum     -= old;
        momentaryEnergyBuffer[(size_t) momentaryIndex] = energy;
        momentarySum     += energy;

        if (++momentaryIndex >= momentaryWindowSamples)
            momentaryIndex = 0;
    }

    // Update short-term sliding window
    if (shortTermWindowSamples > 0)
    {
        const double old = shortTermEnergyBuffer[(size_t) shortTermIndex];
        shortTermSum     -= old;
        shortTermEnergyBuffer[(size_t) shortTermIndex] = energy;
        shortTermSum     += energy;

        if (++shortTermIndex >= shortTermWindowSamples)
            shortTermIndex = 0;
    }

    ++totalSamplesProcessed;

    // Frame scheduling (10 Hz)
    if (--samplesUntilNextFrame <= 0)
    {
        samplesUntilNextFrame += frameSamples;

        // Handle startup before windows are fully "warmed up"
        const double momentaryDenom = (totalSamplesProcessed < (juce::int64) momentaryWindowSamples
                                         ? (double) juce::jmax<juce::int64> (1, totalSamplesProcessed)
                                         : (double) momentaryWindowSamples);

        const double shortTermDenom = (totalSamplesProcessed < (juce::int64) shortTermWindowSamples
                                         ? (double) juce::jmax<juce::int64> (1, totalSamplesProcessed)
                                         : (double) shortTermWindowSamples);

        const double meanMomentary = momentarySum / momentaryDenom;
        const double meanShortTerm = shortTermSum / shortTermDenom;

        const float rmsMomentary = (float) std::sqrt (juce::jmax (0.0, meanMomentary));
        const float rmsShortTerm = (float) std::sqrt (juce::jmax (0.0, meanShortTerm));

        pushLoudnessFrame (rmsMomentary, rmsShortTerm);
    }
}

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
        processSampleForLoudness (channelData.getData(), numChannels, i);

    // Audio is passed through unchanged.
}

//==============================================================================

void LevelScopeAudioProcessor::pushLoudnessFrame (float momentaryRms,
                                                  float shortTermRms) noexcept
{
    int start1 = 0, size1 = 0, start2 = 0, size2 = 0;
    loudnessFifo.prepareToWrite (1, start1, size1, start2, size2);

    const int writable = size1 + size2;
    if (writable > 0)
    {
        const int index = (size1 > 0 ? start1 : start2);

        loudnessBuffer[(size_t) index].momentaryRms = momentaryRms;
        loudnessBuffer[(size_t) index].shortTermRms = shortTermRms;

        loudnessFifo.finishedWrite (1);
    }
    else
    {
        loudnessFifo.finishedWrite (0); // FIFO full, drop frame.
    }
}

int LevelScopeAudioProcessor::readLoudnessFromFifo (float* momentaryDest,
                                                    float* shortTermDest,
                                                    int maxNumToRead) noexcept
{
    if (momentaryDest == nullptr || shortTermDest == nullptr || maxNumToRead <= 0)
        return 0;

    int start1 = 0, size1 = 0, start2 = 0, size2 = 0;
    loudnessFifo.prepareToRead (maxNumToRead, start1, size1, start2, size2);

    const int totalToRead = size1 + size2;

    int destIndex = 0;

    if (size1 > 0)
    {
        for (int i = 0; i < size1; ++i)
        {
            const auto& f = loudnessBuffer[(size_t) (start1 + i)];
            momentaryDest[destIndex] = f.momentaryRms;
            shortTermDest[destIndex] = f.shortTermRms;
            ++destIndex;
        }
    }

    if (size2 > 0)
    {
        for (int i = 0; i < size2; ++i)
        {
            const auto& f = loudnessBuffer[(size_t) (start2 + i)];
            momentaryDest[destIndex] = f.momentaryRms;
            shortTermDest[destIndex] = f.shortTermRms;
            ++destIndex;
        }
    }

    loudnessFifo.finishedRead (totalToRead);

    return totalToRead;
}

//==============================================================================

void LevelScopeAudioProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    juce::ignoreUnused (destData);
    // No parameters yet.
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

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new LevelScopeAudioProcessor();
}
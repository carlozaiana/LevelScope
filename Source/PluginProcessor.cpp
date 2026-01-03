#include "PluginProcessor.h"
#include "PluginEditor.h"

#include <cmath>
#include <algorithm>

juce::int64 LevelScopeAudioProcessor::floorDivInt64 (juce::int64 a, juce::int64 b) noexcept
{
    // b must be > 0
    if (b <= 0) return 0;

    if (a >= 0)
        return a / b;

    // floor division for negatives
    return - ( ( -a + b - 1 ) / b );
}

int LevelScopeAudioProcessor::wrapSlot (juce::int64 absIndex, int capacity) noexcept
{
    if (capacity <= 0) return 0;
    juce::int64 m = absIndex % (juce::int64) capacity;
    if (m < 0) m += (juce::int64) capacity;
    return (int) m;
}

bool LevelScopeAudioProcessor::readEnergyAbs (juce::int64 absFrameIndex, float& outEnergy) const noexcept
{
    outEnergy = 0.0f;
    if (absFrameIndex == (juce::int64) -1)
        return false;

    const int slot = wrapSlot (absFrameIndex, historyCapacityFrames);

    const juce::int64 tag1 = frameIndexTag[(size_t) slot].load (std::memory_order_acquire);
    if (tag1 != absFrameIndex)
        return false;

    const float e = energyMeanSquare[(size_t) slot];

    const juce::int64 tag2 = frameIndexTag[(size_t) slot].load (std::memory_order_acquire);
    if (tag2 != absFrameIndex)
        return false;

    outEnergy = e;
    return true;
}

void LevelScopeAudioProcessor::writeFrameAbs (juce::int64 absFrameIndex, float energyMS, float mRms, float sRms) noexcept
{
    const int slot = wrapSlot (absFrameIndex, historyCapacityFrames);

    // Write data first, publish tag last
    energyMeanSquare[(size_t) slot] = energyMS;
    momentaryRmsHist[(size_t) slot] = mRms;
    shortTermRmsHist[(size_t) slot] = sRms;

    frameIndexTag[(size_t) slot].store (absFrameIndex, std::memory_order_release);

    // Maintain monotonic max written (do not shrink on loops)
    juce::int64 cur = maxWrittenFrameIndex.load (std::memory_order_relaxed);
    while (absFrameIndex > cur && ! maxWrittenFrameIndex.compare_exchange_weak (cur, absFrameIndex))
    {}
}

float LevelScopeAudioProcessor::computeWindowRmsFromEnergy (juce::int64 endFrameIndex, int windowFrames) const noexcept
{
    if (windowFrames <= 0)
        return 0.0f;

    double sum = 0.0;
    int count = 0;

    const juce::int64 start = endFrameIndex - (juce::int64) (windowFrames - 1);

    for (juce::int64 fi = start; fi <= endFrameIndex; ++fi)
    {
        float e = 0.0f;
        if (! readEnergyAbs (fi, e))
            continue;

        sum += (double) e;
        ++count;
    }

    if (count <= 0)
        return 0.0f;

    const double mean = sum / (double) count;
    return (float) std::sqrt (juce::jmax (0.0, mean));
}

//==============================================================================

LevelScopeAudioProcessor::LevelScopeAudioProcessor()
    : AudioProcessor (BusesProperties()
                        .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
                        .withOutput ("Output", juce::AudioChannelSet::stereo(), true)),
      loudnessFifo (loudnessFifoSize),
      loudnessBuffer ((size_t) loudnessFifoSize)
{
    // [TIMELINE-ENERGY] allocate history storage (done once, not on audio thread)
    energyMeanSquare.assign ((size_t) historyCapacityFrames, 0.0f);
    momentaryRmsHist.assign ((size_t) historyCapacityFrames, 0.0f);
    shortTermRmsHist.assign ((size_t) historyCapacityFrames, 0.0f);
    frameIndexTag.resize ((size_t) historyCapacityFrames);
    for (auto& t : frameIndexTag)
        t.store ((juce::int64) -1, std::memory_order_relaxed);
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
    frameEnergyAccum = 0.0;
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
    // Mean-square energy across channels
    double e = 0.0;
    for (int ch = 0; ch < numChannels; ++ch)
    {
        const float s = channelData[ch][sampleIndex];
        e += (double) s * (double) s;
    }
    if (numChannels > 0)
        e /= (double) numChannels;

    frameEnergyAccum += e;

    ++totalSamplesProcessed;

    // 60 Hz frame scheduling
    if (--samplesUntilNextFrame <= 0)
    {
        samplesUntilNextFrame += frameSamples;

        const double meanSquare = frameEnergyAccum / (double) juce::jmax (1, frameSamples);
        frameEnergyAccum = 0.0;

        // Timestamp this loudness frame on the DAW timeline
        lastFrameProjectSample = currentBlockStartProjectSample + (juce::int64) currentBlockSampleIndex;
        lastFrameIsPlaying     = currentBlockIsPlaying;

        // Quantize to our 60 Hz frame grid using floor division (supports negative time)
        const juce::int64 frameIndex = floorDivInt64 (lastFrameProjectSample, (juce::int64) frameSamples);

        // Store base energy first (publish via tag later in writeFrameAbs)
        const float energyMS = (float) juce::jmax (0.0, meanSquare);

        // Fixed windows in 60 Hz frames
        const int momentaryFrames = (int) std::round (momentaryWindowSeconds * loudnessFrameRate); // 0.4s -> 24
        const int shortTermFrames = (int) std::round (shortTermWindowSeconds * loudnessFrameRate); // 3s -> 180

        // Make current frame's energy visible before computing windows
        // (we store it with provisional RMS=0, then overwrite with final RMS)
        writeFrameAbs (frameIndex, energyMS, 0.0f, 0.0f);

        const float rmsMomentary = computeWindowRmsFromEnergy (frameIndex, juce::jmax (1, momentaryFrames));
        const float rmsShortTerm = computeWindowRmsFromEnergy (frameIndex, juce::jmax (1, shortTermFrames));

        // Store final derived values (same abs index, overwrite-safe)
        writeFrameAbs (frameIndex, energyMS, rmsMomentary, rmsShortTerm);

        // Push to GUI FIFO
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

    // [TIMEBASE-PLAYHEAD] Capture DAW playhead position once per block.
    // If not available, fall back to our internal sample counter.
    juce::int64 blockStartProjectSample = totalSamplesProcessed;
    int blockIsPlaying = 1;

    if (auto* ph = getPlayHead())
    {
        juce::AudioPlayHead::CurrentPositionInfo pos;
        if (ph->getCurrentPosition (pos))
        {
            blockStartProjectSample = pos.timeInSamples;
            blockIsPlaying = (pos.isPlaying ? 1 : 0);
        }
    }

    // [FIX-STOP-DIP]
    // When the transport is stopped, many hosts keep calling processBlock with silence.
    // If we keep updating our loudness windows with zeros, the internal state decays
    // and causes a sharp loudness drop on the next restart.
    // Freeze loudness analysis while not playing.
    if (blockIsPlaying == 0)
        return;

    currentBlockStartProjectSample = blockStartProjectSample;
    currentBlockIsPlaying = blockIsPlaying;

    // Prepare read pointers for faster inner loop
    juce::HeapBlock<const float*> channelData ((size_t) numChannels);
    for (int ch = 0; ch < numChannels; ++ch)
        channelData[ch] = buffer.getReadPointer (ch);

    for (int i = 0; i < numSamples; ++i)
    {
        // [TIMEBASE-PLAYHEAD] allow processSampleForLoudness() to timestamp frames
        currentBlockSampleIndex = i;

        processSampleForLoudness (channelData.getData(), numChannels, i);
    }

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

        loudnessBuffer[(size_t) index].frameIndex = floorDivInt64 (lastFrameProjectSample, (juce::int64) frameSamples);
        loudnessBuffer[(size_t) index].isPlaying  = lastFrameIsPlaying;

        loudnessBuffer[(size_t) index].projectSample = lastFrameProjectSample;
        loudnessBuffer[(size_t) index].isPlaying     = lastFrameIsPlaying;

        loudnessFifo.finishedWrite (1);
    }
    else
    {
        loudnessFifo.finishedWrite (0); // FIFO full, drop frame.
    }
}

int LevelScopeAudioProcessor::readLoudnessFromFifo (float* momentaryDest,
                                                    float* shortTermDest,
                                                    juce::int64* frameIndexDest,
                                                    int* isPlayingDest,
                                                    int maxNumToRead) noexcept
{
    if (momentaryDest == nullptr || shortTermDest == nullptr
        || frameIndexDest == nullptr || isPlayingDest == nullptr
        || maxNumToRead <= 0)
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
            momentaryDest[destIndex]  = f.momentaryRms;
            shortTermDest[destIndex]  = f.shortTermRms;
            frameIndexDest[destIndex] = f.frameIndex;
            isPlayingDest[destIndex]  = f.isPlaying;
            ++destIndex;
        }
    }

    if (size2 > 0)
    {
        for (int i = 0; i < size2; ++i)
        {
            const auto& f = loudnessBuffer[(size_t) (start2 + i)];
            momentaryDest[destIndex]  = f.momentaryRms;
            shortTermDest[destIndex]  = f.shortTermRms;
            frameIndexDest[destIndex] = f.frameIndex;
            isPlayingDest[destIndex]  = f.isPlaying;
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
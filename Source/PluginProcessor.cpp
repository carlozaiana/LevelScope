#include "PluginProcessor.h"
#include "PluginEditor.h"

#include <cmath>
#include <algorithm>
#include <type_traits>

namespace
{
    // Treat arithmetic types (int, int64, bool, double) as always "present".
    template <typename T, std::enable_if_t<std::is_arithmetic_v<std::decay_t<T>>, int> = 0>
    bool optHasValue (T) noexcept { return true; }

    template <typename T, std::enable_if_t<std::is_arithmetic_v<std::decay_t<T>>, int> = 0>
    T optValue (T v) noexcept { return v; }

    // Treat optional-like types (std::optional, juce::Optional, etc.) as "present" if bool(o) is true.
    template <typename Opt, std::enable_if_t<!std::is_arithmetic_v<std::decay_t<Opt>>, int> = 0>
    bool optHasValue (const Opt& o) noexcept { return (bool) o; }

    template <typename Opt, std::enable_if_t<!std::is_arithmetic_v<std::decay_t<Opt>>, int> = 0>
    auto optValue (const Opt& o) noexcept -> decltype(*o) { return *o; }
}

juce::int64 LevelScopeAudioProcessor::floorDivInt64 (juce::int64 a, juce::int64 b) noexcept
{
    if (b <= 0) return 0;
    if (a >= 0) return a / b;
    return - ( ( -a + b - 1 ) / b );
}

//==============================================================================

LevelScopeAudioProcessor::LevelScopeAudioProcessor()
    : AudioProcessor (BusesProperties()
                        .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
                        .withOutput ("Output", juce::AudioChannelSet::stereo(), true))
{
    // Phase 2A: core publishes LUFS (UI will be updated in Phase 2B)
    historyModel.setOutputMode (LevelScopeHistoryModel::OutputMode::lufs);
}

LevelScopeAudioProcessor::~LevelScopeAudioProcessor() = default;

//==============================================================================

void LevelScopeAudioProcessor::prepareToPlay (double sampleRate,
                                             int /*samplesPerBlockExpected*/)
{
    currentSampleRate = (sampleRate > 0.0 ? sampleRate : 44100.0);

    frameSamples = juce::jmax (1,
        (int) std::round (currentSampleRate / loudnessFrameRate));

    historyModel.setFrameSamplesForMetadata (frameSamples);

    // [BS1770] prepare K-weighting for current channel count
    const int numCh = juce::jmax (1, getTotalNumInputChannels());
    kWeight.prepare (currentSampleRate, numCh);

    // [BS1770] channel weights (LFE = 0, everything else = 1)
    bs1770ChannelWeights.assign ((size_t) numCh, 1.0f);
    const auto layout = getBusesLayout().getMainInputChannelSet();
    for (int ch = 0; ch < numCh; ++ch)
        if (layout.getTypeOfChannel (ch) == juce::AudioChannelSet::LFE)
            bs1770ChannelWeights[(size_t) ch] = 0.0f;

    resetLoudnessState();
}

void LevelScopeAudioProcessor::releaseResources()
{
    resetLoudnessState();
}

void LevelScopeAudioProcessor::resetLoudnessState() noexcept
{
    samplesUntilNextFrame = frameSamples;
    frameEnergyAccum = 0.0;
    historyModel.resetRealtimeFifo();
    kWeight.reset();
    runningStats.reset();
}

//==============================================================================

bool LevelScopeAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
    const auto& mainIn  = layouts.getMainInputChannelSet();
    const auto& mainOut = layouts.getMainOutputChannelSet();

    if (mainIn.isDisabled() || mainOut.isDisabled())
        return false;

    // Prototype: mono or stereo only
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
    // [BS1770] K-weighted mean-square energy (sum over channels; LFE weight=0)
    double e = 0.0;

    for (int ch = 0; ch < numChannels; ++ch)
    {
        const float w = (ch < (int) bs1770ChannelWeights.size() ? bs1770ChannelWeights[(size_t) ch] : 1.0f);
        if (w <= 0.0f)
            continue;

        const float x = channelData[ch][sampleIndex];
        const float y = kWeight.processSample (ch, x);

        e += (double) w * (double) y * (double) y;
    }

    frameEnergyAccum += e;

    // 60 Hz frame scheduling
    if (--samplesUntilNextFrame <= 0)
    {
        samplesUntilNextFrame += frameSamples;

        const double meanSquare = frameEnergyAccum / (double) juce::jmax (1, frameSamples);
        frameEnergyAccum = 0.0;

        if (skipNextPartialFrameWrite)
        {
            skipNextPartialFrameWrite = false;
            return;
        }

        const juce::int64 frameProjectSample =
            currentBlockStartProjectSample + (juce::int64) currentBlockSampleIndex;

        const juce::int64 frameIndex =
            floorDivInt64 (frameProjectSample, (juce::int64) frameSamples);

        // [FIX-START-RAMP] do not overwrite existing timeline truth during guard period
        if (transportStartOverwriteGuardFrames > 0)
        {
            const bool alreadyExists = historyModel.frameExists (frameIndex);
            --transportStartOverwriteGuardFrames;

            if (alreadyExists)
                return;
        }

        const float energyMS = (float) juce::jmax (0.0, meanSquare);

        const int momentaryFrames = (int) std::round (momentaryWindowSeconds * loudnessFrameRate); // 24
        const int shortTermFrames = (int) std::round (shortTermWindowSeconds * loudnessFrameRate); // 180

        // [LOUDNESS-STATS] update running integrated/LRA
        runningStats.pushFrameEnergy (energyMS);

        historyModel.pushEnergyFrame (frameIndex,
                                     energyMS,
                                     momentaryFrames,
                                     shortTermFrames,
                                     currentBlockIsPlaying);
    }
}

void LevelScopeAudioProcessor::processBlock (juce::AudioBuffer<float>& buffer,
                                            juce::MidiBuffer& midiMessages)
{
    juce::ScopedNoDenormals noDenormals;
    juce::ignoreUnused (midiMessages);

    const int numChannels = getTotalNumInputChannels();
    const int numSamples  = buffer.getNumSamples();

    for (int ch = numChannels; ch < getTotalNumOutputChannels(); ++ch)
        buffer.clear (ch, 0, numSamples);

    if (numChannels <= 0 || numSamples <= 0)
        return;

    // [TIMEBASE-PLAYHEAD]
    juce::int64 blockStartProjectSample = 0;
    int blockIsPlaying = 0;
    bool haveTimeInSamples = false;
    bool haveTimeInSeconds = false;

    juce::int64 hostSamplesValue = 0;
    double hostSecondsValue = 0.0;

    if (auto* ph = getPlayHead())
    {
        if (auto pos = ph->getPosition())
        {
            const auto tSamp = pos->getTimeInSamples();
            if (optHasValue (tSamp))
            {
                hostSamplesValue = (juce::int64) optValue (tSamp);
                haveTimeInSamples = true;
                blockStartProjectSample = hostSamplesValue;
            }

            const auto tSec = pos->getTimeInSeconds();
            if (optHasValue (tSec))
            {
                hostSecondsValue = (double) optValue (tSec);
                haveTimeInSeconds = true;
            }

            const auto pl = pos->getIsPlaying();
            if (optHasValue (pl))
                blockIsPlaying = (optValue (pl) ? 1 : 0);

            haveHostTimeSamples.store (haveTimeInSamples ? 1 : 0, std::memory_order_relaxed);
            haveHostTimeSeconds.store (haveTimeInSeconds ? 1 : 0, std::memory_order_relaxed);
            if (haveTimeInSamples)
                lastHostTimeSamples.store (hostSamplesValue, std::memory_order_relaxed);
            if (haveTimeInSeconds)
                lastHostTimeSeconds.store (hostSecondsValue, std::memory_order_relaxed);

            if (haveTimeInSamples && haveTimeInSeconds && currentSampleRate > 1.0e-12)
            {
                const double sampleSeconds = (double) hostSamplesValue / currentSampleRate;
                timecodeOffsetSeconds.store (hostSecondsValue - sampleSeconds, std::memory_order_relaxed);
            }
        }
    }

    const bool shouldAnalyse = (blockIsPlaying == 1 && haveTimeInSamples);

    // [FIX-RESTART-PARTIAL-FRAME]
    const bool discontinuity =
        (! haveLastBlockEnd) ||
        (lastBlockIsPlaying == 0) ||
        (blockStartProjectSample != lastBlockEndProjectSample);

    if (shouldAnalyse && discontinuity && frameSamples > 0)
    {
        frameEnergyAccum = 0.0;
        kWeight.reset();
        runningStats.reset();

        juce::int64 mod = blockStartProjectSample % (juce::int64) frameSamples;
        if (mod < 0) mod += (juce::int64) frameSamples;

        samplesUntilNextFrame = frameSamples - (int) mod;
        if (samplesUntilNextFrame <= 0)
            samplesUntilNextFrame += frameSamples;

        skipNextPartialFrameWrite = (mod != 0);
    }

    const bool transportStart = (blockIsPlaying == 1 && lastBlockIsPlaying == 0);
    if (transportStart)
        transportStartOverwriteGuardFrames = 6;

    haveLastBlockEnd = true;
    lastBlockEndProjectSample = blockStartProjectSample + (juce::int64) numSamples;
    lastBlockIsPlaying = blockIsPlaying;

    if (! shouldAnalyse)
        return;

    currentBlockStartProjectSample = blockStartProjectSample;
    currentBlockIsPlaying = blockIsPlaying;

    juce::HeapBlock<const float*> channelData ((size_t) numChannels);
    for (int ch = 0; ch < numChannels; ++ch)
        channelData[ch] = buffer.getReadPointer (ch);

    for (int i = 0; i < numSamples; ++i)
    {
        currentBlockSampleIndex = i;
        processSampleForLoudness (channelData.getData(), numChannels, i);
    }

    // Audio passthrough unchanged
}

//==============================================================================
// [CORE->UI] forwarding
//==============================================================================

int LevelScopeAudioProcessor::readLoudnessFromFifo (float* momentaryDest,
                                                   float* shortTermDest,
                                                   juce::int64* frameIndexDest,
                                                   int* isPlayingDest,
                                                   int maxNumToRead) noexcept
{
    return historyModel.readLoudnessFromFifo (momentaryDest,
                                             shortTermDest,
                                             frameIndexDest,
                                             isPlayingDest,
                                             maxNumToRead);
}

//==============================================================================
// [STATE-PERSIST] forward to core (format unchanged)
//==============================================================================

void LevelScopeAudioProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    historyModel.saveState (destData);
}

void LevelScopeAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    historyModel.loadState (data, sizeInBytes);
}

//==============================================================================

juce::AudioProcessorEditor* LevelScopeAudioProcessor::createEditor()
{
    return new LevelScopeAudioProcessorEditor (*this);
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new LevelScopeAudioProcessor();
}
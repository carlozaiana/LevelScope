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

    // [BEGIN LS-PROCESSORCORE-CONSTRUCTOR-EMPTY-GRAPH]
    // Stage A: install an empty/default module graph (no modules => no audible change).
    // Non-RT thread (constructor), allocation OK here.
    {
        auto emptyGraph = std::make_shared<levelscope::ModuleGraph>();
        emptyGraph->revision = 1;
        processorCore.setActiveGraph (std::move (emptyGraph));
    }
    // [END LS-PROCESSORCORE-CONSTRUCTOR-EMPTY-GRAPH]
}

LevelScopeAudioProcessor::~LevelScopeAudioProcessor() = default;

//==============================================================================

// [BEGIN LS-PROCESSORCORE-PREPARETOPLAY]
void LevelScopeAudioProcessor::prepareToPlay (double sampleRate,
                                             int samplesPerBlockExpected)
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

    // Stage A: prepare ProcessorCore (empty graph => no-op).
    // We avoid juce::dsp::ProcessSpec (juce_dsp not enabled).
    // Also, we don't assume ModulePrepareSpec field names here.
    juce::ignoreUnused (samplesPerBlockExpected);
    processorCore.prepare (levelscope::ModulePrepareSpec{});

    resetLoudnessState();
}
// [END LS-PROCESSORCORE-PREPARETOPLAY]

void LevelScopeAudioProcessor::releaseResources()
{
    resetLoudnessState();
}

// [BEGIN LS-PROCESSORCORE-RESETLOUDNESSSTATE]
void LevelScopeAudioProcessor::resetLoudnessState() noexcept
{
    samplesUntilNextFrame = frameSamples;
    frameEnergyAccum = 0.0;
    historyModel.resetRealtimeFifo();
    kWeight.reset();
    runningStats.reset();

    // RT-safe reset; Stage A graph is empty => no audible change
    processorCore.reset();
}
// [END LS-PROCESSORCORE-RESETLOUDNESSSTATE]

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

        // Update running stats first (integrated + LRA, etc.)
        runningStats.pushFrameEnergy (energyMS);

        // LRA "relative gate" curve value (integratedRunning - 20 LU)
        float gateLufs = runningStats.getLraGateLufs();

        // [LRAG-WARMUP-OVERWRITE-FIX]
        // After a stop/start, integrated is not valid immediately (warm-up).
        // If we're overwriting an already-existing timeline frame during warm-up,
        // do NOT overwrite its previously stored gate curve value.
        const float I = runningStats.getIntegratedLufs();
        const bool integratedValid = (I > -199.0f);

        if (! integratedValid && historyModel.frameExists (frameIndex))
        {
            float existingGate = -200.0f;
            if (historyModel.getLraGateLufsAtFrameIndex (frameIndex, existingGate))
                gateLufs = existingGate;
        }

        historyModel.pushEnergyFrame (frameIndex,
                                     energyMS,
                                     momentaryFrames,
                                     shortTermFrames,
                                     currentBlockIsPlaying,
                                     gateLufs);
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

    // [BEGIN LS-PROCESSORCORE-SKIP-ANALYSIS-STILL-PROCESS]
        if (! shouldAnalyse)
        {
            // Stage A: still run ProcessorCore (empty graph => no-op, buffer unchanged).
            levelscope::ProcessContext ctx { buffer,
                                            &midiMessages,
                                            currentSampleRate,
                                            numSamples,
                                            getBusesLayout().getMainInputChannelSet() };

            ctx.isRealtime        = ! isNonRealtime();
            ctx.isPlaying         = (blockIsPlaying == 1);
            ctx.isDiscontinuity   = discontinuity;
            ctx.freezeAnalysis    = true; // matches existing "stop-time silence freeze" behavior intent
            ctx.absoluteSampleIndex = (int64_t) blockStartProjectSample;

            if (frameSamples > 0 && haveTimeInSamples)
            {
                ctx.hasFrameIndex60Hz = true;
                ctx.absoluteFrameIndex60Hz = (int64_t) floorDivInt64 (blockStartProjectSample,
                                                                     (juce::int64) frameSamples);
            }

            processorCore.process (ctx);
            return;
        }
    // [END LS-PROCESSORCORE-SKIP-ANALYSIS-STILL-PROCESS]

    // [BEGIN LS-PROCESSORCORE-PROCESSBLOCK]
    currentBlockStartProjectSample = blockStartProjectSample;
        currentBlockIsPlaying = blockIsPlaying;

        // RT: avoid per-block heap allocation (was juce::HeapBlock)
        const float* const* channelData = buffer.getArrayOfReadPointers();

        for (int i = 0; i < numSamples; ++i)
        {
            currentBlockSampleIndex = i;
            processSampleForLoudness (channelData, numChannels, i);
        }

        // Stage A: run ProcessorCore with an empty/default graph (no-op => no audible change).
        {
            levelscope::ProcessContext ctx { buffer,
                                            &midiMessages,
                                            currentSampleRate,
                                            numSamples,
                                            getBusesLayout().getMainInputChannelSet() };

            ctx.isRealtime         = ! isNonRealtime();
            ctx.isPlaying          = (blockIsPlaying == 1);
            ctx.isDiscontinuity    = discontinuity;
            ctx.freezeAnalysis     = false;
            ctx.absoluteSampleIndex = (int64_t) blockStartProjectSample;

            if (frameSamples > 0 && haveTimeInSamples)
            {
                ctx.hasFrameIndex60Hz = true;
                ctx.absoluteFrameIndex60Hz = (int64_t) floorDivInt64 (blockStartProjectSample,
                                                                     (juce::int64) frameSamples);
            }

            processorCore.process (ctx);
        }

        // Audio passthrough unchanged (empty graph => buffer is unchanged)
    // [END LS-PROCESSORCORE-PROCESSBLOCK]
}

//==============================================================================
// [CORE->UI] forwarding
//==============================================================================

int LevelScopeAudioProcessor::readLoudnessFromFifoEx (float* momentaryDest,
                                                     float* shortTermDest,
                                                     float* lraGateDest,
                                                     juce::int64* frameIndexDest,
                                                     int* isPlayingDest,
                                                     int maxNumToRead) noexcept
{
    return historyModel.readLoudnessFromFifo (momentaryDest,
                                             shortTermDest,
                                             lraGateDest,
                                             frameIndexDest,
                                             isPlayingDest,
                                             maxNumToRead);
}

int LevelScopeAudioProcessor::readLoudnessFromFifo (float* momentaryDest,
                                                   float* shortTermDest,
                                                   juce::int64* frameIndexDest,
                                                   int* isPlayingDest,
                                                   int maxNumToRead) noexcept
{
    // Backward-compatible wrapper (no gate)
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
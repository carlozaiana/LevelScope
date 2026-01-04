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
    frameIndexTag.reset (new std::atomic<juce::int64>[historyCapacityFrames]);
    for (int i = 0; i < historyCapacityFrames; ++i)
        frameIndexTag[(size_t) i].store ((juce::int64) -1, std::memory_order_relaxed);
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

bool LevelScopeAudioProcessor::getDerivedRmsAtFrameIndex (juce::int64 frameIndex,
                                                         float& momentaryRms,
                                                         float& shortTermRms) const noexcept
{
    momentaryRms = 0.0f;
    shortTermRms = 0.0f;

    const int slot = wrapSlot (frameIndex, historyCapacityFrames);

    const juce::int64 tag1 = frameIndexTag[(size_t) slot].load (std::memory_order_acquire);
    if (tag1 != frameIndex)
        return false;

    const float m = momentaryRmsHist[(size_t) slot];
    const float s = shortTermRmsHist[(size_t) slot];

    const juce::int64 tag2 = frameIndexTag[(size_t) slot].load (std::memory_order_acquire);
    if (tag2 != frameIndex)
        return false;

    momentaryRms = m;
    shortTermRms = s;
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

        // [FIX-RESTART-PARTIAL-FRAME]
        // If we started playback mid 60 Hz frame, discard this first partial frame and
        // keep the previously stored timeline value (prevents short spikes/lines).
        if (skipNextPartialFrameWrite)
        {
            skipNextPartialFrameWrite = false;
            return;
        }

        // Timestamp this loudness frame on the DAW timeline
        lastFrameProjectSample = currentBlockStartProjectSample + (juce::int64) currentBlockSampleIndex;
        lastFrameIsPlaying     = currentBlockIsPlaying;

        // Quantize to our 60 Hz frame grid using floor division (supports negative time)
        const juce::int64 frameIndex = floorDivInt64 (lastFrameProjectSample, (juce::int64) frameSamples);

        // [FIX-START-RAMP]
        // If this frame already exists in the timeline, and we're within the short
        // transport-start guard period, do NOT overwrite it (prevents 2–4 dB edges).
        // But if it does NOT exist yet, we still write it to avoid gaps.
        if (transportStartOverwriteGuardFrames > 0)
        {
            float existingEnergy = 0.0f;
            const bool alreadyExists = readEnergyAbs (frameIndex, existingEnergy);

            --transportStartOverwriteGuardFrames;

            if (alreadyExists)
                return; // keep prior timeline truth, skip FIFO push too
        }

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

    // [TIMEBASE-PLAYHEAD]
    juce::int64 blockStartProjectSample = 0;
    int blockIsPlaying = 0;
    bool haveTimeInSamples = false;

    if (auto* ph = getPlayHead())
    {
        if (auto pos = ph->getPosition())
        {
            if (auto t = pos->getTimeInSamples())
            {
                blockStartProjectSample = *t;
                haveTimeInSamples = true;
            }

            blockIsPlaying = (pos->getIsPlaying() ? 1 : 0);
        }
    }

    // [TIMEBASE-PLAYHEAD] only analyze when host provides a real timeline position
    const bool shouldAnalyse = (blockIsPlaying == 1 && haveTimeInSamples);

    // [FIX-RESTART-PARTIAL-FRAME]
    // Detect discontinuities (stop->play, seek, loop jump, host gaps).
    const bool discontinuity =
        (! haveLastBlockEnd) ||
        (lastBlockIsPlaying == 0) ||
        (blockStartProjectSample != lastBlockEndProjectSample);

    if (shouldAnalyse && discontinuity && frameSamples > 0)
    {
        frameEnergyAccum = 0.0;

        // Align frame countdown to the global 60 Hz grid anchored at sample 0.
        juce::int64 mod = blockStartProjectSample % (juce::int64) frameSamples;
        if (mod < 0) mod += (juce::int64) frameSamples; // keep it positive

        samplesUntilNextFrame = frameSamples - (int) mod;
        if (samplesUntilNextFrame <= 0)
            samplesUntilNextFrame += frameSamples;

        // If we start mid-frame, the first computed frame would be partial -> don't overwrite it.
        skipNextPartialFrameWrite = (mod != 0);
    }

    // [FIX-START-RAMP] Detect a transport start (stop -> play)
    const bool transportStart = (blockIsPlaying == 1 && lastBlockIsPlaying == 0);

    // Guard duration in 60 Hz frames (0.1s ~= 6 frames)
    if (transportStart)
        transportStartOverwriteGuardFrames = 6;

    // Update "last block" tracking even if we return early.
    haveLastBlockEnd = true;
    lastBlockEndProjectSample = blockStartProjectSample + (juce::int64) numSamples;
    lastBlockIsPlaying = blockIsPlaying;

    // [FIX-STOP-DIP]
    // When the transport is stopped, many hosts keep calling processBlock with silence.
    // If we keep updating our loudness windows with zeros, the internal state decays
    // and causes a sharp loudness drop on the next restart.
    // Freeze loudness analysis while not playing.
    if (! shouldAnalyse)
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

        loudnessBuffer[(size_t) index].frameIndex =
            floorDivInt64 (lastFrameProjectSample, (juce::int64) frameSamples);
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
// [STATE-PERSIST] Binary, versioned, chunked state
// - Stores base timeline energy + valid mask (compressed)
// - Rebuilds derived momentary/short-term on load
// - Designed to add more chunks later (params, markers, regions, etc.)
//==============================================================================

static constexpr juce::uint32 fourcc (char a, char b, char c, char d) noexcept
{
    return (juce::uint32) (juce::uint8) a
        | ((juce::uint32) (juce::uint8) b << 8)
        | ((juce::uint32) (juce::uint8) c << 16)
        | ((juce::uint32) (juce::uint8) d << 24);
}

void LevelScopeAudioProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    juce::MemoryOutputStream out (destData, true);

    const juce::uint32 kMagic   = fourcc ('L','S','C','P');
    const juce::uint32 kVersion = 1;

    out.writeInt ((int) kMagic);
    out.writeInt ((int) kVersion);

    // chunk count (we write exactly 1 chunk for now)
    out.writeInt (1);

    //--------------------------------------------------------------------------
    // Chunk: HIST
    //--------------------------------------------------------------------------
    juce::MemoryOutputStream chunk (4096);

    // Chunk version
    chunk.writeInt (1);

    // Metadata needed to interpret the history
    chunk.writeInt (historyCapacityFrames);
    chunk.writeInt (frameSamples);

    const int momentaryFrames = (int) std::round (momentaryWindowSeconds * loudnessFrameRate); // 24
    const int shortTermFrames = (int) std::round (shortTermWindowSeconds * loudnessFrameRate); // 180
    chunk.writeInt (momentaryFrames);
    chunk.writeInt (shortTermFrames);

    const juce::int64 maxWritten = maxWrittenFrameIndex.load (std::memory_order_relaxed);
    chunk.writeInt64 (maxWritten);

    // If we have no history, store numFrames = 0
    if (maxWritten == std::numeric_limits<juce::int64>::min())
    {
        chunk.writeInt64 (0);   // startFrameIndex (unused)
        chunk.writeInt (0);     // numFrames
        chunk.writeInt (0);     // compressedBytes
    }
    else
    {
        const juce::int64 startFrameIndex = maxWritten - (juce::int64) (historyCapacityFrames - 1);
        const int numFrames = historyCapacityFrames;

        chunk.writeInt64 (startFrameIndex);
        chunk.writeInt (numFrames);

        // Build (validMask + energy) for this contiguous window
        juce::MemoryOutputStream uncompressed (numFrames * (1 + 4));

        for (int i = 0; i < numFrames; ++i)
        {
            const juce::int64 fi = startFrameIndex + (juce::int64) i;
            const int slot = wrapSlot (fi, historyCapacityFrames);

            const juce::int64 tag = frameIndexTag[(size_t) slot].load (std::memory_order_acquire);
            if (tag == fi)
            {
                const float e = energyMeanSquare[(size_t) slot];
                uncompressed.writeByte ((char) 1);
                uncompressed.writeFloat (e);
            }
            else
            {
                uncompressed.writeByte ((char) 0);
                uncompressed.writeFloat (0.0f);
            }
        }

        // Compress the payload
        juce::MemoryOutputStream compressed;
        {
            juce::GZIPCompressorOutputStream gz (compressed);
            gz.write (uncompressed.getData(), uncompressed.getDataSize());
            gz.flush();
        }

        chunk.writeInt ((int) compressed.getDataSize());
        chunk.write (compressed.getData(), compressed.getDataSize());
    }

    // Write chunk header to main stream
    const juce::uint32 chunkId = fourcc ('H','I','S','T');
    out.writeInt ((int) chunkId);
    out.writeInt ((int) chunk.getDataSize());
    out.write (chunk.getData(), chunk.getDataSize());
}

void LevelScopeAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    if (data == nullptr || sizeInBytes <= 0)
        return;

    juce::MemoryInputStream in (data, (size_t) sizeInBytes, false);

    const juce::uint32 kMagicExpected = fourcc ('L','S','C','P');
    const int magic = in.readInt();
    const int version = in.readInt();

    if ((juce::uint32) magic != kMagicExpected || version != 1)
        return;

    const int numChunks = in.readInt();
    if (numChunks <= 0)
        return;

    // Reset tags (do NOT clear stored arrays necessarily, but tags define validity)
    for (int i = 0; i < historyCapacityFrames; ++i)
        frameIndexTag[(size_t) i].store ((juce::int64) -1, std::memory_order_relaxed);

    maxWrittenFrameIndex.store (std::numeric_limits<juce::int64>::min(), std::memory_order_relaxed);

    for (int ci = 0; ci < numChunks; ++ci)
    {
        const juce::uint32 chunkId = (juce::uint32) in.readInt();
        const int chunkBytes = in.readInt();
        if (chunkBytes <= 0 || (juce::int64) in.getNumBytesRemaining() < (juce::int64) chunkBytes)
        {
            in.skipNextBytes (juce::jmax (0, chunkBytes));
            continue;
        }

        juce::MemoryBlock chunkData ((size_t) chunkBytes);
        in.read (chunkData.getData(), (size_t) chunkBytes);

        if (chunkId != fourcc ('H','I','S','T'))
            continue;

        juce::MemoryInputStream ch (chunkData.getData(), chunkData.getSize(), false);

        const int chunkVer = ch.readInt();
        if (chunkVer != 1)
            continue;

        const int capStored = ch.readInt();
        const int frameSamplesStored = ch.readInt();
        const int momentaryFramesStored = ch.readInt();
        const int shortTermFramesStored = ch.readInt();

        const juce::int64 maxWrittenStored = ch.readInt64();
        const juce::int64 startFrameIndex = ch.readInt64();
        const int numFrames = ch.readInt();
        const int compressedBytes = ch.readInt();

        juce::ignoreUnused (capStored, frameSamplesStored);

        if (maxWrittenStored == std::numeric_limits<juce::int64>::min() || numFrames <= 0 || compressedBytes <= 0)
        {
            maxWrittenFrameIndex.store (maxWrittenStored, std::memory_order_relaxed);
            continue;
        }

        if ((int) ch.getNumBytesRemaining() < compressedBytes)
            continue;

        juce::MemoryBlock compressed ((size_t) compressedBytes);
        ch.read (compressed.getData(), (size_t) compressedBytes);

        juce::MemoryInputStream compIn (compressed.getData(), compressed.getSize(), false);
        juce::GZIPDecompressorInputStream gzIn (compIn);

        // Rolling window sums with "valid count" to match timeline-truth semantics
        const int mWin = juce::jmax (1, momentaryFramesStored);
        const int sWin = juce::jmax (1, shortTermFramesStored);

        std::vector<float> mBuf ((size_t) mWin, 0.0f);
        std::vector<juce::uint8> mValid ((size_t) mWin, 0);
        int mPos = 0;
        double mSum = 0.0;
        int mCount = 0;

        std::vector<float> sBuf ((size_t) sWin, 0.0f);
        std::vector<juce::uint8> sValid ((size_t) sWin, 0);
        int sPos = 0;
        double sSum = 0.0;
        int sCount = 0;

        for (int i = 0; i < numFrames; ++i)
        {
            const juce::int64 fi = startFrameIndex + (juce::int64) i;

            const int v = (int) gzIn.readByte();
            const float e = gzIn.readFloat();

            // Update momentary window
            if (mValid[(size_t) mPos])
            {
                mSum -= (double) mBuf[(size_t) mPos];
                --mCount;
            }
            if (v != 0)
            {
                mBuf[(size_t) mPos] = e;
                mValid[(size_t) mPos] = 1;
                mSum += (double) e;
                ++mCount;
            }
            else
            {
                mBuf[(size_t) mPos] = 0.0f;
                mValid[(size_t) mPos] = 0;
            }
            mPos = (mPos + 1) % mWin;

            // Update short-term window
            if (sValid[(size_t) sPos])
            {
                sSum -= (double) sBuf[(size_t) sPos];
                --sCount;
            }
            if (v != 0)
            {
                sBuf[(size_t) sPos] = e;
                sValid[(size_t) sPos] = 1;
                sSum += (double) e;
                ++sCount;
            }
            else
            {
                sBuf[(size_t) sPos] = 0.0f;
                sValid[(size_t) sPos] = 0;
            }
            sPos = (sPos + 1) % sWin;

            if (v == 0)
                continue;

            const float mRms = (mCount > 0 ? (float) std::sqrt (juce::jmax (0.0, mSum / (double) mCount)) : 0.0f);
            const float sRms = (sCount > 0 ? (float) std::sqrt (juce::jmax (0.0, sSum / (double) sCount)) : 0.0f);

            writeFrameAbs (fi, e, mRms, sRms);
        }

        maxWrittenFrameIndex.store (maxWrittenStored, std::memory_order_relaxed);
        loudnessFifo.reset(); // avoid stale GUI frames from before load
    }
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
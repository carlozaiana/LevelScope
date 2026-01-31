#include "LevelScopeHistoryModel.h"

#include <cmath>
#include <algorithm>

//==============================================================================
// Helpers
//==============================================================================

juce::int64 LevelScopeHistoryModel::floorDivInt64 (juce::int64 a, juce::int64 b) noexcept
{
    if (b <= 0) return 0;
    if (a >= 0) return a / b;
    return - ( ( -a + b - 1 ) / b );
}

int LevelScopeHistoryModel::wrapSlot (juce::int64 absIndex, int capacity) noexcept
{
    if (capacity <= 0) return 0;
    juce::int64 m = absIndex % (juce::int64) capacity;
    if (m < 0) m += (juce::int64) capacity;
    return (int) m;
}

float LevelScopeHistoryModel::energyToLufs (float meanSquare) noexcept
{
    const double e = (double) meanSquare;
    if (e <= 0.0)
        return -200.0f;

    const double lufs = -0.691 + 10.0 * std::log10 (e);
    return (float) lufs;
}

//==============================================================================

LevelScopeHistoryModel::LevelScopeHistoryModel()
    : loudnessFifo (loudnessFifoSize),
      loudnessBuffer ((size_t) loudnessFifoSize)
{
    energyMeanSquare.assign   ((size_t) historyCapacityFrames, 0.0f);
    momentaryValueHist.assign ((size_t) historyCapacityFrames, 0.0f);
    shortTermValueHist.assign ((size_t) historyCapacityFrames, 0.0f);
    lraGateLufsHist.assign    ((size_t) historyCapacityFrames, -200.0f);

    frameIndexTag.reset (new std::atomic<juce::int64>[historyCapacityFrames]);
    for (int i = 0; i < historyCapacityFrames; ++i)
        frameIndexTag[(size_t) i].store ((juce::int64) -1, std::memory_order_relaxed);
}

//==============================================================================
// Timeline read/write
//==============================================================================

bool LevelScopeHistoryModel::readEnergyAbs (juce::int64 absFrameIndex, float& outEnergy) const noexcept
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

float LevelScopeHistoryModel::computeWindowMeanEnergy (juce::int64 endFrameIndex, int windowFrames) const noexcept
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

    return (float) (sum / (double) count);
}

void LevelScopeHistoryModel::writeFrameAbs (juce::int64 absFrameIndex,
                                           float energyMS,
                                           float momentaryValue,
                                           float shortTermValue,
                                           float lraGateLufs) noexcept
{
    const int slot = wrapSlot (absFrameIndex, historyCapacityFrames);

    energyMeanSquare[(size_t) slot]   = energyMS;
    momentaryValueHist[(size_t) slot] = momentaryValue;
    shortTermValueHist[(size_t) slot] = shortTermValue;
    lraGateLufsHist[(size_t) slot]    = lraGateLufs;

    frameIndexTag[(size_t) slot].store (absFrameIndex, std::memory_order_release);

    juce::int64 cur = maxWrittenFrameIndex.load (std::memory_order_relaxed);
    while (absFrameIndex > cur && ! maxWrittenFrameIndex.compare_exchange_weak (cur, absFrameIndex))
    {}
}

bool LevelScopeHistoryModel::frameExists (juce::int64 absFrameIndex) const noexcept
{
    const int slot = wrapSlot (absFrameIndex, historyCapacityFrames);
    return frameIndexTag[(size_t) slot].load (std::memory_order_acquire) == absFrameIndex;
}

//==============================================================================
// Queries for GUI bootstrap
//==============================================================================

bool LevelScopeHistoryModel::getDerivedLufsAtFrameIndex (juce::int64 frameIndex,
                                                        float& momentaryValue,
                                                        float& shortTermValue) const noexcept
{
    momentaryValue = 0.0f;
    shortTermValue = 0.0f;

    const int slot = wrapSlot (frameIndex, historyCapacityFrames);

    const juce::int64 tag1 = frameIndexTag[(size_t) slot].load (std::memory_order_acquire);
    if (tag1 != frameIndex)
        return false;

    const float m = momentaryValueHist[(size_t) slot];
    const float s = shortTermValueHist[(size_t) slot];

    const juce::int64 tag2 = frameIndexTag[(size_t) slot].load (std::memory_order_acquire);
    if (tag2 != frameIndex)
        return false;

    momentaryValue = m;
    shortTermValue = s;
    return true;
}

bool LevelScopeHistoryModel::getLraGateLufsAtFrameIndex (juce::int64 frameIndex,
                                                        float& lraGateLufsOut) const noexcept
{
    lraGateLufsOut = -200.0f;

    const int slot = wrapSlot (frameIndex, historyCapacityFrames);

    const juce::int64 tag1 = frameIndexTag[(size_t) slot].load (std::memory_order_acquire);
    if (tag1 != frameIndex)
        return false;

    const float g = lraGateLufsHist[(size_t) slot];

    const juce::int64 tag2 = frameIndexTag[(size_t) slot].load (std::memory_order_acquire);
    if (tag2 != frameIndex)
        return false;

    lraGateLufsOut = g;
    return true;
}

//==============================================================================
// FIFO
//==============================================================================

void LevelScopeHistoryModel::resetRealtimeFifo() noexcept
{
    loudnessFifo.reset();
}

void LevelScopeHistoryModel::pushToFifo (juce::int64 absFrameIndex,
                                        float momentaryValue,
                                        float shortTermValue,
                                        float lraGateLufs,
                                        int isPlaying) noexcept
{
    int start1 = 0, size1 = 0, start2 = 0, size2 = 0;
    loudnessFifo.prepareToWrite (1, start1, size1, start2, size2);

    const int writable = size1 + size2;
    if (writable <= 0)
    {
        loudnessFifo.finishedWrite (0);
        return;
    }

    const int index = (size1 > 0 ? start1 : start2);

    auto& f = loudnessBuffer[(size_t) index];
    f.momentaryValue = momentaryValue;
    f.shortTermValue = shortTermValue;
    f.lraGateLufs    = lraGateLufs;
    f.frameIndex     = absFrameIndex;
    f.isPlaying      = isPlaying;

    loudnessFifo.finishedWrite (1);
}

int LevelScopeHistoryModel::readLoudnessFromFifo (float* momentaryDest,
                                                 float* shortTermDest,
                                                 float* lraGateDest,
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
    int di = 0;

    auto copyOne = [&] (const LoudnessFrame& f)
    {
        momentaryDest[di]  = f.momentaryValue;
        shortTermDest[di]  = f.shortTermValue;
        if (lraGateDest != nullptr)
            lraGateDest[di] = f.lraGateLufs;
        frameIndexDest[di] = f.frameIndex;
        isPlayingDest[di]  = f.isPlaying;
        ++di;
    };

    for (int i = 0; i < size1; ++i) copyOne (loudnessBuffer[(size_t) (start1 + i)]);
    for (int i = 0; i < size2; ++i) copyOne (loudnessBuffer[(size_t) (start2 + i)]);

    loudnessFifo.finishedRead (totalToRead);
    return totalToRead;
}

//==============================================================================
// Realtime push
//==============================================================================

void LevelScopeHistoryModel::pushEnergyFrame (juce::int64 absFrameIndex,
                                             float energyMS,
                                             int momentaryFrames,
                                             int shortTermFrames,
                                             int isPlaying,
                                             float lraGateLufs) noexcept
{
    energyMS = (float) juce::jmax (0.0, (double) energyMS);

    // Placeholder write so energy is visible for the window computations
    writeFrameAbs (absFrameIndex, energyMS, -200.0f, -200.0f, lraGateLufs);

    const float meanE_M = computeWindowMeanEnergy (absFrameIndex, juce::jmax (1, momentaryFrames));
    const float meanE_S = computeWindowMeanEnergy (absFrameIndex, juce::jmax (1, shortTermFrames));

    float outM = 0.0f;
    float outS = 0.0f;

    if (getOutputMode() == OutputMode::lufs)
    {
        outM = energyToLufs (meanE_M);
        outS = energyToLufs (meanE_S);
    }
    else
    {
        outM = (float) std::sqrt (juce::jmax (0.0f, meanE_M));
        outS = (float) std::sqrt (juce::jmax (0.0f, meanE_S));
    }

    writeFrameAbs (absFrameIndex, energyMS, outM, outS, lraGateLufs);
    pushToFifo (absFrameIndex, outM, outS, lraGateLufs, isPlaying);
}

//==============================================================================
// State (HIST + LRAG + TCOF)
//==============================================================================

// [BEGIN LS-STATE-SAVESTATE-WRAPPER]
void LevelScopeHistoryModel::saveState (juce::MemoryBlock& destData) const
{
    saveState (destData, nullptr, nullptr);
}
// [END LS-STATE-SAVESTATE-WRAPPER]

// [BEGIN LS-STATE-SAVESTATE-ADDITIVE]
void LevelScopeHistoryModel::saveState (juce::MemoryBlock& destData,
                                       const juce::MemoryBlock* apvsChunkData,
                                       const juce::MemoryBlock* modgChunkData) const
{
    juce::MemoryOutputStream out (destData, true);

    const juce::uint32 kMagic   = fourcc ('L','S','C','P');
    const juce::uint32 kVersion = 1;

    out.writeInt ((int) kMagic);
    out.writeInt ((int) kVersion);

    const auto hasExtraChunk = [] (const juce::MemoryBlock* mb) noexcept
    {
        return (mb != nullptr && mb->getSize() > 0);
    };

    const int extraChunks =
        (hasExtraChunk (apvsChunkData) ? 1 : 0) +
        (hasExtraChunk (modgChunkData) ? 1 : 0);

    // Baseline chunk count is 3 (HIST + LRAG + TCOF). We append extras after.
    out.writeInt (3 + extraChunks);

    const juce::int64 maxWritten = maxWrittenFrameIndex.load (std::memory_order_relaxed);

    //--------------------------------------------------------------------------
    // Chunk: HIST (energy + valid mask)  [BASELINE - unchanged]
    //--------------------------------------------------------------------------
    {
        juce::MemoryOutputStream chunk (4096);

        chunk.writeInt (1); // chunk version
        chunk.writeInt (historyCapacityFrames);
        chunk.writeInt (frameSamplesForMetadata);

        const int momentaryFrames = (int) std::round (0.4 * loudnessFrameRate); // 24
        const int shortTermFrames = (int) std::round (3.0 * loudnessFrameRate); // 180
        chunk.writeInt (momentaryFrames);
        chunk.writeInt (shortTermFrames);

        chunk.writeInt64 (maxWritten);

        if (maxWritten == std::numeric_limits<juce::int64>::min())
        {
            chunk.writeInt64 (0);
            chunk.writeInt (0);
            chunk.writeInt (0);
        }
        else
        {
            const juce::int64 startFrameIndex = maxWritten - (juce::int64) (historyCapacityFrames - 1);
            const int numFrames = historyCapacityFrames;

            chunk.writeInt64 (startFrameIndex);
            chunk.writeInt (numFrames);

            juce::MemoryOutputStream uncompressed (numFrames * (1 + 4));

            for (int i = 0; i < numFrames; ++i)
            {
                const juce::int64 fi = startFrameIndex + (juce::int64) i;
                const int slot = wrapSlot (fi, historyCapacityFrames);

                const juce::int64 tag = frameIndexTag[(size_t) slot].load (std::memory_order_acquire);
                if (tag == fi)
                {
                    uncompressed.writeByte ((char) 1);
                    uncompressed.writeFloat (energyMeanSquare[(size_t) slot]);
                }
                else
                {
                    uncompressed.writeByte ((char) 0);
                    uncompressed.writeFloat (0.0f);
                }
            }

            juce::MemoryOutputStream compressed;
            {
                juce::GZIPCompressorOutputStream gz (compressed);
                gz.write (uncompressed.getData(), uncompressed.getDataSize());
                gz.flush();
            }

            chunk.writeInt ((int) compressed.getDataSize());
            chunk.write (compressed.getData(), compressed.getDataSize());
        }

        const juce::uint32 chunkId = fourcc ('H','I','S','T');
        out.writeInt ((int) chunkId);
        out.writeInt ((int) chunk.getDataSize());
        out.write (chunk.getData(), chunk.getDataSize());
    }

    //--------------------------------------------------------------------------
    // Chunk: LRAG (gate curve + valid mask)  [BASELINE - unchanged]
    //--------------------------------------------------------------------------
    {
        juce::MemoryOutputStream chunk (4096);

        chunk.writeInt (1); // chunk version
        chunk.writeInt (historyCapacityFrames);
        chunk.writeInt64 (maxWritten);

        if (maxWritten == std::numeric_limits<juce::int64>::min())
        {
            chunk.writeInt64 (0);
            chunk.writeInt (0);
            chunk.writeInt (0);
        }
        else
        {
            const juce::int64 startFrameIndex = maxWritten - (juce::int64) (historyCapacityFrames - 1);
            const int numFrames = historyCapacityFrames;

            chunk.writeInt64 (startFrameIndex);
            chunk.writeInt (numFrames);

            juce::MemoryOutputStream uncompressed (numFrames * (1 + 4));

            for (int i = 0; i < numFrames; ++i)
            {
                const juce::int64 fi = startFrameIndex + (juce::int64) i;
                const int slot = wrapSlot (fi, historyCapacityFrames);

                const juce::int64 tag = frameIndexTag[(size_t) slot].load (std::memory_order_acquire);
                if (tag == fi)
                {
                    uncompressed.writeByte ((char) 1);
                    uncompressed.writeFloat (lraGateLufsHist[(size_t) slot]);
                }
                else
                {
                    uncompressed.writeByte ((char) 0);
                    uncompressed.writeFloat (0.0f);
                }
            }

            juce::MemoryOutputStream compressed;
            {
                juce::GZIPCompressorOutputStream gz (compressed);
                gz.write (uncompressed.getData(), uncompressed.getDataSize());
                gz.flush();
            }

            chunk.writeInt ((int) compressed.getDataSize());
            chunk.write (compressed.getData(), compressed.getDataSize());
        }

        const juce::uint32 chunkId = fourcc ('L','R','A','G');
        out.writeInt ((int) chunkId);
        out.writeInt ((int) chunk.getDataSize());
        out.write (chunk.getData(), chunk.getDataSize());
    }

    //--------------------------------------------------------------------------
    // Chunk: TCOF  [BASELINE - unchanged]
    //--------------------------------------------------------------------------
    {
        juce::MemoryOutputStream tc (64);
        tc.writeInt (1);
        tc.writeDouble (userTimecodeOffsetSeconds.load (std::memory_order_relaxed));

        const juce::uint32 tcId = fourcc ('T','C','O','F');
        out.writeInt ((int) tcId);
        out.writeInt ((int) tc.getDataSize());
        out.write (tc.getData(), tc.getDataSize());
    }

    //--------------------------------------------------------------------------
    // Chunk: APVS (APVTS state)  [Stage C2 - additive]
    //--------------------------------------------------------------------------
    if (hasExtraChunk (apvsChunkData))
    {
        const juce::uint32 id = fourcc ('A','P','V','S');
        out.writeInt ((int) id);
        out.writeInt ((int) apvsChunkData->getSize());
        out.write (apvsChunkData->getData(), apvsChunkData->getSize());
    }

    //--------------------------------------------------------------------------
    // Chunk: MODG (module graph state)  [Stage C2 - additive]
    //--------------------------------------------------------------------------
    if (hasExtraChunk (modgChunkData))
    {
        const juce::uint32 id = fourcc ('M','O','D','G');
        out.writeInt ((int) id);
        out.writeInt ((int) modgChunkData->getSize());
        out.write (modgChunkData->getData(), modgChunkData->getSize());
    }
}
// [END LS-STATE-SAVESTATE-ADDITIVE]

// [BEGIN LS-STATE-LOADSTATE-WRAPPER]
void LevelScopeHistoryModel::loadState (const void* data, int sizeInBytes)
{
    loadState (data, sizeInBytes, nullptr, nullptr);
}
// [END LS-STATE-LOADSTATE-WRAPPER]

// [BEGIN LS-STATE-LOADSTATE-ADDITIVE]
void LevelScopeHistoryModel::loadState (const void* data, int sizeInBytes,
                                       juce::MemoryBlock* apvsChunkOut,
                                       juce::MemoryBlock* modgChunkOut)
{
    if (apvsChunkOut != nullptr) apvsChunkOut->reset();
    if (modgChunkOut != nullptr) modgChunkOut->reset();

    if (data == nullptr || sizeInBytes <= 0)
        return;

    juce::MemoryInputStream in (data, (size_t) sizeInBytes, false);

    const juce::uint32 kMagicExpected = fourcc ('L','S','C','P');
    const int magic = in.readInt();
    const int version = in.readInt();

    // Container version unchanged (v1). New chunks are optional.
    if ((juce::uint32) magic != kMagicExpected || version != 1)
        return;

    const int numChunks = in.readInt();
    if (numChunks <= 0)
        return;

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

        // Stage C2 additive chunks: capture raw payload and continue
        if (chunkId == fourcc ('A','P','V','S'))
        {
            if (apvsChunkOut != nullptr)
                *apvsChunkOut = std::move (chunkData);
            continue;
        }

        if (chunkId == fourcc ('M','O','D','G'))
        {
            if (modgChunkOut != nullptr)
                *modgChunkOut = std::move (chunkData);
            continue;
        }

        if (chunkId == fourcc ('T','C','O','F'))
        {
            juce::MemoryInputStream tc (chunkData.getData(), chunkData.getSize(), false);
            if (tc.readInt() == 1)
                userTimecodeOffsetSeconds.store (tc.readDouble(), std::memory_order_relaxed);
            continue;
        }

        //--------------------------------------------------------------------------
        // HIST
        //--------------------------------------------------------------------------
        if (chunkId == fourcc ('H','I','S','T'))
        {
            juce::MemoryInputStream ch (chunkData.getData(), chunkData.getSize(), false);
            if (ch.readInt() != 1)
                continue;

            const int capStored = ch.readInt();
            const int frameSamplesStored = ch.readInt();
            const int momentaryFramesStored = ch.readInt();
            const int shortTermFramesStored = ch.readInt();
            juce::ignoreUnused (capStored, frameSamplesStored);

            const juce::int64 maxWrittenStored = ch.readInt64();
            const juce::int64 startFrameIndex  = ch.readInt64();
            const int numFrames                = ch.readInt();
            const int compressedBytes          = ch.readInt();

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

            // Rolling windows (valid-count semantics)
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

                // momentary window
                if (mValid[(size_t) mPos]) { mSum -= (double) mBuf[(size_t) mPos]; --mCount; }
                if (v != 0) { mBuf[(size_t) mPos] = e; mValid[(size_t) mPos] = 1; mSum += (double) e; ++mCount; }
                else        { mBuf[(size_t) mPos] = 0.0f; mValid[(size_t) mPos] = 0; }
                mPos = (mPos + 1) % mWin;

                // short-term window
                if (sValid[(size_t) sPos]) { sSum -= (double) sBuf[(size_t) sPos]; --sCount; }
                if (v != 0) { sBuf[(size_t) sPos] = e; sValid[(size_t) sPos] = 1; sSum += (double) e; ++sCount; }
                else        { sBuf[(size_t) sPos] = 0.0f; sValid[(size_t) sPos] = 0; }
                sPos = (sPos + 1) % sWin;

                if (v == 0)
                    continue;

                const float meanE_M = (mCount > 0 ? (float) (mSum / (double) mCount) : 0.0f);
                const float meanE_S = (sCount > 0 ? (float) (sSum / (double) sCount) : 0.0f);

                float outM = 0.0f;
                float outS = 0.0f;

                if (getOutputMode() == OutputMode::lufs)
                {
                    outM = energyToLufs (meanE_M);
                    outS = energyToLufs (meanE_S);
                }
                else
                {
                    outM = (float) std::sqrt (juce::jmax (0.0f, meanE_M));
                    outS = (float) std::sqrt (juce::jmax (0.0f, meanE_S));
                }

                // Preserve any gate value that may have been loaded earlier (LRAG before HIST)
                float gate = -200.0f;
                {
                    const int slot = wrapSlot (fi, historyCapacityFrames);
                    const juce::int64 tag = frameIndexTag[(size_t) slot].load (std::memory_order_acquire);
                    if (tag == fi)
                        gate = lraGateLufsHist[(size_t) slot];
                }

                writeFrameAbs (fi, e, outM, outS, gate);
            }

            maxWrittenFrameIndex.store (maxWrittenStored, std::memory_order_relaxed);
            loudnessFifo.reset();
            continue;
        }

        //--------------------------------------------------------------------------
        // LRAG
        //--------------------------------------------------------------------------
        if (chunkId == fourcc ('L','R','A','G'))
        {
            juce::MemoryInputStream ch (chunkData.getData(), chunkData.getSize(), false);
            if (ch.readInt() != 1)
                continue;

            const int capStored = ch.readInt();
            juce::ignoreUnused (capStored);

            const juce::int64 maxWrittenStored = ch.readInt64();
            const juce::int64 startFrameIndex  = ch.readInt64();
            const int numFrames                = ch.readInt();
            const int compressedBytes          = ch.readInt();

            if (maxWrittenStored == std::numeric_limits<juce::int64>::min() || numFrames <= 0 || compressedBytes <= 0)
                continue;

            if ((int) ch.getNumBytesRemaining() < compressedBytes)
                continue;

            juce::MemoryBlock compressed ((size_t) compressedBytes);
            ch.read (compressed.getData(), (size_t) compressedBytes);

            juce::MemoryInputStream compIn (compressed.getData(), compressed.getSize(), false);
            juce::GZIPDecompressorInputStream gzIn (compIn);

            for (int i = 0; i < numFrames; ++i)
            {
                const juce::int64 fi = startFrameIndex + (juce::int64) i;

                const int v = (int) gzIn.readByte();
                const float g = gzIn.readFloat();

                if (v == 0)
                    continue;

                const int slot = wrapSlot (fi, historyCapacityFrames);
                lraGateLufsHist[(size_t) slot] = g;

                // Re-publish tag so readers see consistent data
                frameIndexTag[(size_t) slot].store (fi, std::memory_order_release);
            }

            loudnessFifo.reset();
            continue;
        }

        // Unknown chunk: ignore
    }
}
// [END LS-STATE-LOADSTATE-ADDITIVE]
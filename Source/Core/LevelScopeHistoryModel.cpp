#include "LevelScopeHistoryModel.h"

#include <cmath>
#include <algorithm>

//==============================================================================
// [CORE] Helpers
//==============================================================================

float LevelScopeHistoryModel::rmsToLufs (float rms) noexcept
{
    const double r = (double) rms;
    if (r <= 0.0)
        return -200.0f; // effectively -inf for display purposes

    // BS.1770 loudness: LUFS = -0.691 + 10*log10(meanSquare)
    // meanSquare = rms^2  =>  10*log10(rms^2) = 20*log10(rms)
    const double lufs = -0.691 + 20.0 * std::log10 (r);
    return (float) lufs;
}

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

float LevelScopeHistoryModel::energyToLufs (float meanSquareEnergy) noexcept
{
    // ITU-R BS.1770: LKFS/LUFS = -0.691 + 10*log10( mean_square_energy )
    // (Momentary/Short-term are ungated.)
    constexpr float kOffset = -0.691f;
    constexpr float kEps    = 1.0e-12f; // avoid -inf

    const float e = juce::jmax (kEps, meanSquareEnergy);
    return kOffset + 10.0f * (float) std::log10 ((double) e);
}

//==============================================================================

LevelScopeHistoryModel::LevelScopeHistoryModel()
    : loudnessFifo (loudnessFifoSize),
      loudnessBuffer ((size_t) loudnessFifoSize)
{
    energyMeanSquare.assign ((size_t) historyCapacityFrames, 0.0f);
    momentaryLufsHist.assign ((size_t) historyCapacityFrames, -120.0f);
    shortTermLufsHist.assign ((size_t) historyCapacityFrames, -120.0f);

    frameIndexTag.reset (new std::atomic<juce::int64>[historyCapacityFrames]);
    for (int i = 0; i < historyCapacityFrames; ++i)
        frameIndexTag[(size_t) i].store ((juce::int64) -1, std::memory_order_relaxed);
}

//==============================================================================
// [CORE-TIMELINE] read/write
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

void LevelScopeHistoryModel::writeFrameAbs (juce::int64 absFrameIndex, float energyMS, float mLufs, float sLufs) noexcept
{
    const int slot = wrapSlot (absFrameIndex, historyCapacityFrames);

    energyMeanSquare[(size_t) slot]  = energyMS;
    momentaryLufsHist[(size_t) slot] = mLufs;
    shortTermLufsHist[(size_t) slot] = sLufs;

    frameIndexTag[(size_t) slot].store (absFrameIndex, std::memory_order_release);

    juce::int64 cur = maxWrittenFrameIndex.load (std::memory_order_relaxed);
    while (absFrameIndex > cur && ! maxWrittenFrameIndex.compare_exchange_weak (cur, absFrameIndex))
    {}
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

bool LevelScopeHistoryModel::frameExists (juce::int64 absFrameIndex) const noexcept
{
    const int slot = wrapSlot (absFrameIndex, historyCapacityFrames);
    const juce::int64 tag = frameIndexTag[(size_t) slot].load (std::memory_order_acquire);
    return (tag == absFrameIndex);
}

//==============================================================================
// [CORE->UI] Derived LUFS query
//==============================================================================

bool LevelScopeHistoryModel::getDerivedLufsAtFrameIndex (juce::int64 frameIndex,
                                                        float& momentaryLufs,
                                                        float& shortTermLufs) const noexcept
{
    momentaryLufs = -120.0f;
    shortTermLufs = -120.0f;

    const int slot = wrapSlot (frameIndex, historyCapacityFrames);

    const juce::int64 tag1 = frameIndexTag[(size_t) slot].load (std::memory_order_acquire);
    if (tag1 != frameIndex)
        return false;

    const float m = momentaryLufsHist[(size_t) slot];
    const float s = shortTermLufsHist[(size_t) slot];

    const juce::int64 tag2 = frameIndexTag[(size_t) slot].load (std::memory_order_acquire);
    if (tag2 != frameIndex)
        return false;

    momentaryLufs = m;
    shortTermLufs = s;
    return true;
}

//==============================================================================
// [CORE-REALTIME] FIFO push/pop
//==============================================================================

void LevelScopeHistoryModel::resetRealtimeFifo() noexcept
{
    loudnessFifo.reset();
}

void LevelScopeHistoryModel::pushToFifo (juce::int64 absFrameIndex,
                                        float momentaryLufs,
                                        float shortTermLufs,
                                        int isPlaying) noexcept
{
    int start1 = 0, size1 = 0, start2 = 0, size2 = 0;
    loudnessFifo.prepareToWrite (1, start1, size1, start2, size2);

    const int writable = size1 + size2;
    if (writable > 0)
    {
        const int index = (size1 > 0 ? start1 : start2);

        auto& f = loudnessBuffer[(size_t) index];
        f.momentaryLufs = momentaryLufs;
        f.shortTermLufs = shortTermLufs;
        f.frameIndex    = absFrameIndex;
        f.isPlaying     = isPlaying;

        loudnessFifo.finishedWrite (1);
    }
    else
    {
        loudnessFifo.finishedWrite (0);
    }
}

int LevelScopeHistoryModel::readLoudnessFromFifo (float* momentaryDest,
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
            momentaryDest[destIndex]  = f.momentaryLufs;
            shortTermDest[destIndex]  = f.shortTermLufs;
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
            momentaryDest[destIndex]  = f.momentaryLufs;
            shortTermDest[destIndex]  = f.shortTermLufs;
            frameIndexDest[destIndex] = f.frameIndex;
            isPlayingDest[destIndex]  = f.isPlaying;
            ++destIndex;
        }
    }

    loudnessFifo.finishedRead (totalToRead);
    return totalToRead;
}

//==============================================================================
// [CORE-REALTIME] per-frame entry point
//==============================================================================

void LevelScopeHistoryModel::pushEnergyFrame (juce::int64 absFrameIndex,
                                             float energyMS,
                                             int momentaryFrames,
                                             int shortTermFrames,
                                             int isPlaying) noexcept
{
    energyMS = (float) juce::jmax (0.0, (double) energyMS);

    const auto mode = getOutputMode();

    // Make energy visible first (derived values are placeholders)
    // (Derived placeholders do not affect the energy-window computation.)
    if (mode == OutputMode::lufs)
        writeFrameAbs (absFrameIndex, energyMS, -120.0f, -120.0f);
    else
        writeFrameAbs (absFrameIndex, energyMS, 0.0f, 0.0f);

    const float meanE_M = computeWindowMeanEnergy (absFrameIndex, juce::jmax (1, momentaryFrames));
    const float meanE_S = computeWindowMeanEnergy (absFrameIndex, juce::jmax (1, shortTermFrames));

    float outM = 0.0f;
    float outS = 0.0f;

    if (mode == OutputMode::lufs)
    {
        outM = energyToLufs (meanE_M);
        outS = energyToLufs (meanE_S);
    }
    else
    {
        // Legacy mode: publish RMS (linear). UI will convert to dB.
        outM = (float) std::sqrt (juce::jmax (0.0f, meanE_M));
        outS = (float) std::sqrt (juce::jmax (0.0f, meanE_S));
    }

    writeFrameAbs (absFrameIndex, energyMS, outM, outS);
    pushToFifo (absFrameIndex, outM, outS, isPlaying);
}

//==============================================================================
// [CORE-STATE] save/load
//==============================================================================

void LevelScopeHistoryModel::saveState (juce::MemoryBlock& destData) const
{
    juce::MemoryOutputStream out (destData, true);

    const juce::uint32 kMagic   = fourcc ('L','S','C','P');
    const juce::uint32 kVersion = 1;

    out.writeInt ((int) kMagic);
    out.writeInt ((int) kVersion);

    out.writeInt (2); // chunk count (HIST + TCOF)

    //--------------------------------------------------------------------------
    // Chunk: HIST
    //--------------------------------------------------------------------------
    juce::MemoryOutputStream chunk (4096);

    chunk.writeInt (2); // chunk version (was 1)
    chunk.writeInt ((int) storedEnergyKind); // new in v2

    chunk.writeInt (historyCapacityFrames);
    chunk.writeInt (frameSamplesForMetadata);

    const int momentaryFrames = (int) std::round (0.4 * loudnessFrameRate);
    const int shortTermFrames = (int) std::round (3.0 * loudnessFrameRate);
    chunk.writeInt (momentaryFrames);
    chunk.writeInt (shortTermFrames);

    const juce::int64 maxWritten = maxWrittenFrameIndex.load (std::memory_order_relaxed);
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

    //--------------------------------------------------------------------------
    // Chunk: TCOF
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
}

void LevelScopeHistoryModel::loadState (const void* data, int sizeInBytes)
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

        if (chunkId == fourcc ('T','C','O','F'))
        {
            juce::MemoryInputStream tc (chunkData.getData(), chunkData.getSize(), false);
            if (tc.readInt() == 1)
                userTimecodeOffsetSeconds.store (tc.readDouble(), std::memory_order_relaxed);
            continue;
        }

        if (chunkId != fourcc ('H','I','S','T'))
            continue;

        juce::MemoryInputStream ch (chunkData.getData(), chunkData.getSize(), false);

        const int chunkVer = ch.readInt();

        int energyKindInt = (int) EnergyKind::unknownOrLegacy;

        if (chunkVer == 2)
            energyKindInt = ch.readInt();
        else if (chunkVer != 1)
            continue;

        storedEnergyKind = (EnergyKind) energyKindInt;

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

        // Rolling sums for mean energy
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

            writeFrameAbs (fi, e, outM, outS);
        }

        maxWrittenFrameIndex.store (maxWrittenStored, std::memory_order_relaxed);
        loudnessFifo.reset();
    }
}
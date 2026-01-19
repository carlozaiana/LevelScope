#pragma once

#include <juce_core/juce_core.h>
#include <atomic>
#include <vector>
#include <memory>
#include <limits>

//==============================================================================
// LevelScopeHistoryModel
// Owns:
// - timeline truth ring keyed by absolute 60 Hz frame index (supports negative)
// - GUI FIFO with newest frames
// - plugin state persistence (chunked, gzip)
//
// Phase 2 additions:
// - publishes momentary/short-term as LUFS (dB) (OutputMode::lufs)
// - also stores/publishes an LRA-relative-gate curve value per frame (LRAG)
//==============================================================================

class LevelScopeHistoryModel
{
public:
    static constexpr double loudnessFrameRate = 60.0;
    static constexpr double historyLengthSeconds = 3.0 * 3600.0;

    static constexpr int historyCapacityFrames =
        (int) (historyLengthSeconds * loudnessFrameRate + 0.5); // ~648000

    LevelScopeHistoryModel();
    ~LevelScopeHistoryModel() = default;

    //==============================================================================
    // [CORE-OUTPUT-MODE]
    // - rms  : derived values are RMS (linear)
    // - lufs : derived values are LUFS (dB)
    enum class OutputMode { rms = 0, lufs = 1 };

    void setOutputMode (OutputMode m) noexcept { outputMode.store ((int) m, std::memory_order_relaxed); }
    OutputMode getOutputMode() const noexcept  { return (OutputMode) outputMode.load (std::memory_order_relaxed); }

    //==============================================================================
    // [CORE-REALTIME] Push one 60 Hz frame (audio thread)
    //
    // energyMeanSquare: BS.1770 K-weighted mean-square energy for that 60 Hz frame.
    // momentaryFrames/shortTermFrames: window lengths in 60 Hz frames (e.g. 24 and 180).
    // lraGateLufs: current LRA gate value (typically IntegratedRunning - 20 LU).
    void pushEnergyFrame (juce::int64 absFrameIndex,
                          float energyMeanSquare,
                          int momentaryFrames,
                          int shortTermFrames,
                          int isPlaying,
                          float lraGateLufs) noexcept;

    bool frameExists (juce::int64 absFrameIndex) const noexcept;

    //==============================================================================
    // [CORE->UI] FIFO read
    //
    // lraGateDest can be nullptr if the caller doesn't need it.
    int readLoudnessFromFifo (float* momentaryDest,
                              float* shortTermDest,
                              float* lraGateDest,
                              juce::int64* frameIndexDest,
                              int* isPlayingDest,
                              int maxNumToRead) noexcept;

    // Backward-compatible convenience wrapper (no gate)
    int readLoudnessFromFifo (float* momentaryDest,
                              float* shortTermDest,
                              juce::int64* frameIndexDest,
                              int* isPlayingDest,
                              int maxNumToRead) noexcept
    {
        return readLoudnessFromFifo (momentaryDest, shortTermDest, nullptr, frameIndexDest, isPlayingDest, maxNumToRead);
    }

    //==============================================================================
    // [CORE->UI] Timeline truth queries (used for GUI bootstrap after state load)
    bool getDerivedLufsAtFrameIndex (juce::int64 frameIndex,
                                    float& momentaryValue,
                                    float& shortTermValue) const noexcept;

    bool getLraGateLufsAtFrameIndex (juce::int64 frameIndex,
                                     float& lraGateLufsOut) const noexcept;

    juce::int64 getMaxWrittenFrameIndex() const noexcept
    {
        return maxWrittenFrameIndex.load (std::memory_order_relaxed);
    }

    void resetRealtimeFifo() noexcept;

    //==============================================================================
    // [CORE-STATE] Binary, chunked, gzip-compressed
    // Chunks:
    // - HIST (base energy + valid mask)
    // - LRAG (gate curve + valid mask)
    // - TCOF (user timecode offset)
    void saveState (juce::MemoryBlock& destData) const;
    void loadState (const void* data, int sizeInBytes);

    //==============================================================================
    // [CORE-TIMECODE-USER]
    double getUserTimecodeOffsetSeconds() const noexcept
    {
        return userTimecodeOffsetSeconds.load (std::memory_order_relaxed);
    }

    void setUserTimecodeOffsetSeconds (double s) noexcept
    {
        userTimecodeOffsetSeconds.store (s, std::memory_order_relaxed);
    }

    void setFrameSamplesForMetadata (int fs) noexcept { frameSamplesForMetadata = fs; }

private:
    //==============================================================================
    // [CORE-TIMELINE] overwrite-safe ring
    //==============================================================================

    std::vector<float> energyMeanSquare;      // base
    std::vector<float> momentaryValueHist;    // derived (LUFS or RMS)
    std::vector<float> shortTermValueHist;    // derived (LUFS or RMS)
    std::vector<float> lraGateLufsHist;       // derived/control curve (LUFS)

    std::unique_ptr<std::atomic<juce::int64>[]> frameIndexTag; // -1 = empty, else abs frameIndex
    std::atomic<juce::int64> maxWrittenFrameIndex { std::numeric_limits<juce::int64>::min() };

    int frameSamplesForMetadata = 0;

    std::atomic<int> outputMode { (int) OutputMode::lufs };

    // Helpers (negative-time safe)
    static juce::int64 floorDivInt64 (juce::int64 a, juce::int64 b) noexcept; // b>0
    static int wrapSlot (juce::int64 absIndex, int capacity) noexcept;

    bool readEnergyAbs (juce::int64 absFrameIndex, float& outEnergy) const noexcept;

    // Mean energy over last N frames (skips missing frames)
    float computeWindowMeanEnergy (juce::int64 endFrameIndex, int windowFrames) const noexcept;

    // BS.1770 conversion helpers (meanSquare -> LUFS)
    static float energyToLufs (float meanSquare) noexcept;

    void writeFrameAbs (juce::int64 absFrameIndex,
                        float energyMS,
                        float momentaryValue,
                        float shortTermValue,
                        float lraGateLufs) noexcept;

    //==============================================================================
    // [CORE-FIFO]
    //==============================================================================

    struct LoudnessFrame
    {
        float momentaryValue = 0.0f;
        float shortTermValue = 0.0f;
        float lraGateLufs    = -200.0f;

        juce::int64 frameIndex = 0;
        int isPlaying = 1;
    };

    static constexpr int loudnessFifoSize = 4096;

    juce::AbstractFifo         loudnessFifo;
    std::vector<LoudnessFrame> loudnessBuffer;

    void pushToFifo (juce::int64 absFrameIndex,
                     float momentaryValue,
                     float shortTermValue,
                     float lraGateLufs,
                     int isPlaying) noexcept;

    //==============================================================================
    // [CORE-STATE]
    //==============================================================================

    static constexpr juce::uint32 fourcc (char a, char b, char c, char d) noexcept
    {
        return (juce::uint32) (juce::uint8) a
            | ((juce::uint32) (juce::uint8) b << 8)
            | ((juce::uint32) (juce::uint8) c << 16)
            | ((juce::uint32) (juce::uint8) d << 24);
    }

    std::atomic<double> userTimecodeOffsetSeconds { 0.0 };
};
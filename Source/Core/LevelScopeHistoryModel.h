#pragma once

#include <juce_core/juce_core.h>
#include <atomic>
#include <vector>
#include <memory>
#include <limits>

//==============================================================================
// LevelScopeHistoryModel
//
// [CORE] Owns "timeline truth" storage + lock-free GUI FIFO + state persistence.
//
// Phase 2:
// - Stored base measure = K-weighted mean-square energy per 60 Hz frame.
// - Derived values = Momentary + Short-term LOUDNESS in LUFS (float, dB-like).
//==============================================================================

class LevelScopeHistoryModel
{
public:
    static constexpr double loudnessFrameRate = 60.0;
    static constexpr double historyLengthSeconds = 3.0 * 3600.0;

    static constexpr int historyCapacityFrames =
        (int) (historyLengthSeconds * loudnessFrameRate + 0.5);

    enum class EnergyKind : int
    {
        unknownOrLegacy = 0, // old sessions (pre K-weight energy)
        kWeighted1770   = 1  // current
    };

    LevelScopeHistoryModel();
    ~LevelScopeHistoryModel() = default;

    //==============================================================================
    // [CORE-OUTPUT-MODE]
    // Controls what the model stores/publishes as "momentary" and "short-term".
    // Default is RMS (legacy behavior).
    // Phase 2A part 2 will switch to LUFS.
    enum class OutputMode
    {
        rms = 0,   // legacy: derived values are RMS (linear)
        lufs = 1   // derived values are LUFS (dB)
    };

    void setOutputMode (OutputMode m) noexcept { outputMode.store ((int) m, std::memory_order_relaxed); }
    OutputMode getOutputMode() const noexcept  { return (OutputMode) outputMode.load (std::memory_order_relaxed); }

    //==============================================================================
    // [CORE-REALTIME] The plugin calls this once per 60 Hz frame.
    // energyMeanSquare must be the BS.1770 K-weighted energy sum across channels.
    void pushEnergyFrame (juce::int64 absFrameIndex,
                          float energyMeanSquare,
                          int momentaryFrames,
                          int shortTermFrames,
                          int isPlaying) noexcept;

    bool frameExists (juce::int64 absFrameIndex) const noexcept;

    //==============================================================================
    // [CORE->UI] FIFO read (values are LUFS, not RMS)
    int readLoudnessFromFifo (float* momentaryDest,
                              float* shortTermDest,
                              juce::int64* frameIndexDest,
                              int* isPlayingDest,
                              int maxNumToRead) noexcept;

    // [CORE->UI] Timeline truth (LUFS) used by GUI bootstrap
    bool getDerivedLufsAtFrameIndex (juce::int64 frameIndex,
                                    float& momentaryLufs,
                                    float& shortTermLufs) const noexcept;

    juce::int64 getMaxWrittenFrameIndex() const noexcept
    {
        return maxWrittenFrameIndex.load (std::memory_order_relaxed);
    }

    void resetRealtimeFifo() noexcept;

    //==============================================================================
    // [CORE-STATE] same outer format: magic 'LSCP', version 1, chunked.
    // HIST chunk version is bumped to 2 to include EnergyKind.
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

    // [CORE-META]
    void setFrameSamplesForMetadata (int fs) noexcept { frameSamplesForMetadata = fs; }

    EnergyKind getStoredEnergyKind() const noexcept { return storedEnergyKind; }
    void       setStoredEnergyKind (EnergyKind k) noexcept { storedEnergyKind = k; }

private:
    //==============================================================================
    // [CORE-TIMELINE] overwrite-safe ring buffer keyed by absolute frameIndex
    //==============================================================================

    std::vector<float> energyMeanSquare;   // base energy per frame
    std::vector<float> momentaryLufsHist;  // derived LUFS
    std::vector<float> shortTermLufsHist;  // derived LUFS

    std::unique_ptr<std::atomic<juce::int64>[]> frameIndexTag; // -1 empty, else abs frame index
    std::atomic<juce::int64> maxWrittenFrameIndex { std::numeric_limits<juce::int64>::min() };

    int frameSamplesForMetadata = 0;
    EnergyKind storedEnergyKind = EnergyKind::kWeighted1770;

    // Helpers (negative-safe)
    static juce::int64 floorDivInt64 (juce::int64 a, juce::int64 b) noexcept;
    static int wrapSlot (juce::int64 absIndex, int capacity) noexcept;

    bool readEnergyAbs (juce::int64 absFrameIndex, float& outEnergy) const noexcept;
    void writeFrameAbs (juce::int64 absFrameIndex, float energyMS, float mLufs, float sLufs) noexcept;

    float computeWindowMeanEnergy (juce::int64 endFrameIndex, int windowFrames) const noexcept;

    static float energyToLufs (float meanSquareEnergy) noexcept;

    //==============================================================================
    // [CORE-REALTIME-FIFO]
    //==============================================================================

    struct LoudnessFrame
    {
        float momentaryLufs = -120.0f;
        float shortTermLufs = -120.0f;
        juce::int64 frameIndex = 0;
        int isPlaying = 1;
    };

    static constexpr int loudnessFifoSize = 4096;

    juce::AbstractFifo         loudnessFifo;
    std::vector<LoudnessFrame> loudnessBuffer;

    void pushToFifo (juce::int64 absFrameIndex,
                     float momentaryLufs,
                     float shortTermLufs,
                     int isPlaying) noexcept;

    //==============================================================================
    // [CORE-STATE] helpers
    //==============================================================================

        std::atomic<int> outputMode { (int) OutputMode::rms };

    static constexpr juce::uint32 fourcc (char a, char b, char c, char d) noexcept
    {
        return (juce::uint32) (juce::uint8) a
            | ((juce::uint32) (juce::uint8) b << 8)
            | ((juce::uint32) (juce::uint8) c << 16)
            | ((juce::uint32) (juce::uint8) d << 24);
    }

    std::atomic<double> userTimecodeOffsetSeconds { 0.0 };
};
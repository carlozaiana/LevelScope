#pragma once

#include <juce_core/juce_core.h>
#include <atomic>
#include <vector>
#include <memory>
#include <limits>

//==============================================================================
// LevelScopeHistoryModel
//
// [CORE] Owns the "timeline truth" storage + lock-free GUI FIFO + state persistence.
// This is extracted from PluginProcessor to become reusable by:
// - future LevelFlow modules
// - standalone app (offline analysis)
// - sibling products (scopes/meters)
//
// IMPORTANT: This class does NOT know about DAW playhead logic/discontinuities yet.
// The plugin wrapper decides *when* to write frames.
//==============================================================================

class LevelScopeHistoryModel
{
public:
    // Must match existing project behavior/state format assumptions.
    static constexpr double loudnessFrameRate = 60.0;          // 60 Hz timeline grid
    static constexpr double historyLengthSeconds = 3.0 * 3600.0; // 3 hours

    static constexpr int historyCapacityFrames =
        (int) (historyLengthSeconds * loudnessFrameRate + 0.5); // ~648000

    LevelScopeHistoryModel();
    ~LevelScopeHistoryModel() = default;

    //==============================================================================
    // [CORE-REALTIME] The plugin calls this once per 60 Hz frame.
    //
    // - absFrameIndex: absolute frame index on the 60 Hz grid (can be negative).
    // - energyMeanSquare: mean-square energy for that frame (linear).
    // - momentaryFrames / shortTermFrames: window lengths in 60 Hz frames
    //   (e.g. 24 and 180).
    // - isPlaying: passed through into the GUI FIFO.
    //
    // This function:
    //  1) writes energy into the ring at absFrameIndex,
    //  2) computes derived momentary/short-term RMS from stored energy,
    //  3) stores derived RMS back into the ring,
    //  4) pushes a frame into the GUI FIFO.
    void pushEnergyFrame (juce::int64 absFrameIndex,
                          float energyMeanSquare,
                          int momentaryFrames,
                          int shortTermFrames,
                          int isPlaying) noexcept;

    // Used by the plugin to implement "start ramp overwrite guard"
    bool frameExists (juce::int64 absFrameIndex) const noexcept;

    //==============================================================================
    // [CORE->UI] FIFO read API (matches your existing UI calls)
    int readLoudnessFromFifo (float* momentaryDest,
                              float* shortTermDest,
                              juce::int64* frameIndexDest,
                              int* isPlayingDest,
                              int maxNumToRead) noexcept;

    // [CORE->UI] Timeline truth queries (used for GUI bootstrap after state load)
    bool getDerivedRmsAtFrameIndex (juce::int64 frameIndex,
                                   float& momentaryRms,
                                   float& shortTermRms) const noexcept;

    juce::int64 getMaxWrittenFrameIndex() const noexcept
    {
        return maxWrittenFrameIndex.load (std::memory_order_relaxed);
    }

    void resetRealtimeFifo() noexcept;

    //==============================================================================
    // [CORE-STATE] The plugin forwards getState/setState to this model.
    // Format is intentionally kept identical to the previous PluginProcessor version:
    // magic 'LSCP', version 1, chunked, gzip-compressed HIST + TCOF.
    void saveState (juce::MemoryBlock& destData) const;
    void loadState (const void* data, int sizeInBytes);

    //==============================================================================
    // [CORE-TIMECODE-USER] persisted user display offset (seconds)
    double getUserTimecodeOffsetSeconds() const noexcept
    {
        return userTimecodeOffsetSeconds.load (std::memory_order_relaxed);
    }

    void setUserTimecodeOffsetSeconds (double s) noexcept
    {
        userTimecodeOffsetSeconds.store (s, std::memory_order_relaxed);
    }

    // [CORE-META] Stored for state metadata only (debug/backwards compatibility)
    void setFrameSamplesForMetadata (int fs) noexcept { frameSamplesForMetadata = fs; }

private:
    //==============================================================================
    // [CORE-TIMELINE-ENERGY] overwrite-safe ring buffer keyed by absolute frameIndex
    //==============================================================================

    std::vector<float> energyMeanSquare;   // base measure (per 60 Hz frame)
    std::vector<float> momentaryRmsHist;   // derived RMS from energy
    std::vector<float> shortTermRmsHist;   // derived RMS from energy

    // Publish-by-tag ring: -1 = empty, else absolute frame index present
    std::unique_ptr<std::atomic<juce::int64>[]> frameIndexTag;

    std::atomic<juce::int64> maxWrittenFrameIndex { std::numeric_limits<juce::int64>::min() };

    int frameSamplesForMetadata = 0;

    // Helpers (negative-time safe)
    static juce::int64 floorDivInt64 (juce::int64 a, juce::int64 b) noexcept; // b>0
    static int wrapSlot (juce::int64 absIndex, int capacity) noexcept;

    bool readEnergyAbs (juce::int64 absFrameIndex, float& outEnergy) const noexcept;
    void writeFrameAbs (juce::int64 absFrameIndex, float energyMS, float mRms, float sRms) noexcept;

    float computeWindowRmsFromEnergy (juce::int64 endFrameIndex, int windowFrames) const noexcept;

    //==============================================================================
    // [CORE-REALTIME-FIFO]
    //==============================================================================

    struct LoudnessFrame
    {
        float momentaryRms = 0.0f;
        float shortTermRms = 0.0f;
        juce::int64 frameIndex = 0;
        int isPlaying = 1;
    };

    static constexpr int loudnessFifoSize = 4096;

    juce::AbstractFifo         loudnessFifo;
    std::vector<LoudnessFrame> loudnessBuffer;

    void pushToFifo (juce::int64 absFrameIndex,
                     float momentaryRms,
                     float shortTermRms,
                     int isPlaying) noexcept;

    //==============================================================================
    // [CORE-STATE] helpers
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

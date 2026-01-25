#pragma once

#include <juce_core/juce_core.h>
#include <array>
#include <atomic>

//==============================================================================
// [LOUDNESS-STATS] Running loudness statistics (realtime-friendly)
//
// Computes (running):
// - Integrated loudness (EBU-gated), updated at 10 Hz
// - LRA (Tech 3342-ish), updated at 1 Hz
// - Rolling LRA (windowed), updated at 1 Hz, selectable window: 30/60/120 s
//
// Input:
// - pushFrameEnergy() called once per 60 Hz frame with BS.1770 K-weighted mean-square energy.
//
// Threading:
// - Audio thread calls pushFrameEnergy().
// - UI thread reads atomics and may set rolling window seconds via atomic.
//==============================================================================

class RunningLoudnessStats
{
public:
    RunningLoudnessStats();
    void reset() noexcept;

    void pushFrameEnergy (float frameMeanSquareEnergy) noexcept;

    float getIntegratedLufs() const noexcept { return integratedLufs.load (std::memory_order_relaxed); }
    float getLraLu() const noexcept          { return lraLu.load (std::memory_order_relaxed); }
    float getLraGateLufs() const noexcept    { return lraGateLufs.load (std::memory_order_relaxed); }

    // Rolling LRA (windowed)
    void setRollingWindowSeconds (int seconds) noexcept;
    int  getRollingWindowSeconds() const noexcept { return rollingWindowSeconds.load (std::memory_order_relaxed); }
    float getRollingLraLu() const noexcept        { return rollingLraLu.load (std::memory_order_relaxed); }

private:
    //==============================================================================
    // Conversions (BS.1770)
    //==============================================================================

    static float energyToLufs (float meanSquare) noexcept;

    // Histogram parameters (0.1 LU bins)
    static constexpr float histMinLufs = -120.0f;
    static constexpr float histMaxLufs =  20.0f;
    static constexpr float histStep    =   0.1f;
    static constexpr int   histBins    = (int) ((histMaxLufs - histMinLufs) / histStep + 1.0f); // 1401

    static int   lufsToBin (float lufs) noexcept;
    static float binToLufs (int bin) noexcept;

    // Gates
    static constexpr float absGateLufs = -70.0f;

    //==============================================================================
    // Integrated (10 Hz, 400ms blocks with 100ms step)
    //==============================================================================

    static constexpr int framesPerSecond = 60;

    static constexpr int mWindowFrames  = 24; // 400ms @ 60 Hz
    static constexpr int mStepFrames    = 6;  // 100ms @ 60 Hz

    std::array<float, mWindowFrames> winM {};
    int   winMPos = 0;
    int   winMFilled = 0;
    double winMSum = 0.0;

    int frameCounter = 0;

    // Integrated histogram (absolute-gated blocks only)
    std::array<int,    histBins> intCount {};
    std::array<double, histBins> intEnergySum {};

    int    intTotalCount = 0;
    double intTotalEnergySum = 0.0;

    //==============================================================================
    // LRA (program-style) (1 Hz, short-term from 3s window)
    //==============================================================================

    static constexpr int sWindowFrames = 180; // 3s @ 60 Hz
    static constexpr int sStepFrames   = 60;  // compute once per second

    std::array<float, sWindowFrames> winS {};
    int   winSPos = 0;
    int   winSFilled = 0;
    double winSSum = 0.0;

    std::array<int, histBins> stCount {};
    int stTotalCount = 0;

    //==============================================================================
    // Rolling LRA (windowed, 30/60/120 seconds, computed from short-term @1Hz)
    //==============================================================================

    static constexpr int maxRollingSeconds = 120;

    std::array<float, maxRollingSeconds> rollingStLufs {};
    int rollingPos = 0;
    int rollingFilled = 0;

    std::atomic<int>   rollingWindowSeconds { 60 };
    std::atomic<float> rollingLraLu { 0.0f };

    float computeRollingLra (float integratedForGate) const noexcept;

    //==============================================================================
    // Published results (atomics)
    //==============================================================================

    std::atomic<float> integratedLufs { -200.0f };
    std::atomic<float> lraLu          { 0.0f };
    std::atomic<float> lraGateLufs    { -200.0f };

    float computeIntegratedGatedLufs() const noexcept;
    float computeLraFromHistogram (float integratedForGate) const noexcept;
};
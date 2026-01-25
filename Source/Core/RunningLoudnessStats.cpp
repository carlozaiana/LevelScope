#include "RunningLoudnessStats.h"

#include <cmath>
#include <algorithm>

//==============================================================================

RunningLoudnessStats::RunningLoudnessStats()
{
    reset();
}

void RunningLoudnessStats::reset() noexcept
{
    winMPos = 0;
    winMFilled = 0;
    winMSum = 0.0;

    winSPos = 0;
    winSFilled = 0;
    winSSum = 0.0;

    frameCounter = 0;

    intTotalCount = 0;
    intTotalEnergySum = 0.0;
    intCount.fill (0);
    intEnergySum.fill (0.0);

    stTotalCount = 0;
    stCount.fill (0);

    rollingPos = 0;
    rollingFilled = 0;
    rollingStLufs.fill (-200.0f);

    integratedLufs.store (-200.0f, std::memory_order_relaxed);
    lraLu.store (0.0f, std::memory_order_relaxed);
    lraGateLufs.store (-200.0f, std::memory_order_relaxed);
    rollingLraLu.store (0.0f, std::memory_order_relaxed);
}

//==============================================================================

float RunningLoudnessStats::energyToLufs (float meanSquare) noexcept
{
    const double e = (double) meanSquare;
    if (e <= 0.0)
        return -200.0f;

    const double lufs = -0.691 + 10.0 * std::log10 (e);
    return (float) lufs;
}

int RunningLoudnessStats::lufsToBin (float lufs) noexcept
{
    const float clamped = juce::jlimit (histMinLufs, histMaxLufs, lufs);
    const int bin = (int) std::floor ((clamped - histMinLufs) / histStep + 0.5f);
    return juce::jlimit (0, histBins - 1, bin);
}

float RunningLoudnessStats::binToLufs (int bin) noexcept
{
    bin = juce::jlimit (0, histBins - 1, bin);
    return histMinLufs + (float) bin * histStep;
}

//==============================================================================
// Rolling window control
//==============================================================================

void RunningLoudnessStats::setRollingWindowSeconds (int seconds) noexcept
{
    // Only allow 30/60/120 for now.
    int s = seconds;

    if (s <= 30) s = 30;
    else if (s <= 60) s = 60;
    else s = 120;

    rollingWindowSeconds.store (s, std::memory_order_relaxed);
}

//==============================================================================
// Main update
//==============================================================================

void RunningLoudnessStats::pushFrameEnergy (float frameMeanSquareEnergy) noexcept
{
    frameMeanSquareEnergy = juce::jmax (0.0f, frameMeanSquareEnergy);

    //--------------------------------------------------------------------------
    // Update 400ms window (for integrated blocks @ 10Hz)
    //--------------------------------------------------------------------------
    if (winMFilled < mWindowFrames)
    {
        winM[(size_t) winMPos] = frameMeanSquareEnergy;
        winMSum += (double) frameMeanSquareEnergy;
        ++winMFilled;
    }
    else
    {
        winMSum -= (double) winM[(size_t) winMPos];
        winM[(size_t) winMPos] = frameMeanSquareEnergy;
        winMSum += (double) frameMeanSquareEnergy;
    }

    winMPos = (winMPos + 1) % mWindowFrames;

    //--------------------------------------------------------------------------
    // Update 3s window (for short-term @ 1Hz)
    //--------------------------------------------------------------------------
    if (winSFilled < sWindowFrames)
    {
        winS[(size_t) winSPos] = frameMeanSquareEnergy;
        winSSum += (double) frameMeanSquareEnergy;
        ++winSFilled;
    }
    else
    {
        winSSum -= (double) winS[(size_t) winSPos];
        winS[(size_t) winSPos] = frameMeanSquareEnergy;
        winSSum += (double) frameMeanSquareEnergy;
    }

    winSPos = (winSPos + 1) % sWindowFrames;

    //--------------------------------------------------------------------------
    ++frameCounter;

    // 10 Hz integrated update
    if ((frameCounter % mStepFrames) == 0 && winMFilled >= mWindowFrames)
    {
        const float blockMeanE = (float) (winMSum / (double) mWindowFrames);
        const float blockLufs  = energyToLufs (blockMeanE);

        // Absolute gate first (EBU)
        if (blockLufs >= absGateLufs)
        {
            const int bin = lufsToBin (blockLufs);

            ++intCount[(size_t) bin];
            intEnergySum[(size_t) bin] += (double) blockMeanE;

            ++intTotalCount;
            intTotalEnergySum += (double) blockMeanE;
        }

        const float I = computeIntegratedGatedLufs();
        integratedLufs.store (I, std::memory_order_relaxed);
        lraGateLufs.store (I - 20.0f, std::memory_order_relaxed);
    }

    // 1 Hz short-term + LRA update
    if ((frameCounter % sStepFrames) == 0 && winSFilled >= sWindowFrames)
    {
        const float stMeanE = (float) (winSSum / (double) sWindowFrames);
        const float stLufs  = energyToLufs (stMeanE);

        // Store for rolling LRA (store even if very low; gating happens later)
        rollingStLufs[(size_t) rollingPos] = stLufs;
        rollingPos = (rollingPos + 1) % maxRollingSeconds;
        rollingFilled = juce::jmin (rollingFilled + 1, maxRollingSeconds);

        // Program-style LRA histogram uses abs gate
        if (stLufs >= absGateLufs)
        {
            const int bin = lufsToBin (stLufs);
            ++stCount[(size_t) bin];
            ++stTotalCount;
        }

        const float I = integratedLufs.load (std::memory_order_relaxed);

        const float LRA = computeLraFromHistogram (I);
        lraLu.store (LRA, std::memory_order_relaxed);

        const float rolling = computeRollingLra (I);
        rollingLraLu.store (rolling, std::memory_order_relaxed);
    }
}

//==============================================================================
// Integrated (EBU-style gating; running)
//==============================================================================

float RunningLoudnessStats::computeIntegratedGatedLufs() const noexcept
{
    if (intTotalCount <= 0 || intTotalEnergySum <= 0.0)
        return -200.0f;

    // "Ungated" here means: after absolute gate only (standard approach).
    const float meanE_abs = (float) (intTotalEnergySum / (double) intTotalCount);
    const float I_abs = energyToLufs (meanE_abs);

    const float relGate = I_abs - 10.0f;
    const int gateBin = lufsToBin (relGate);

    double sumE = 0.0;
    int count = 0;

    for (int b = gateBin; b < histBins; ++b)
    {
        count += intCount[(size_t) b];
        sumE  += intEnergySum[(size_t) b];
    }

    if (count <= 0 || sumE <= 0.0)
        return I_abs;

    const float meanE_gated = (float) (sumE / (double) count);
    return energyToLufs (meanE_gated);
}

//==============================================================================
// LRA (running; histogram percentile based)
//==============================================================================

float RunningLoudnessStats::computeLraFromHistogram (float integratedForGate) const noexcept
{
    if (stTotalCount <= 0)
        return 0.0f;

    const float gate = juce::jmax (absGateLufs, integratedForGate - 20.0f);
    const int gateBin = lufsToBin (gate);

    int gatedCount = 0;
    for (int b = gateBin; b < histBins; ++b)
        gatedCount += stCount[(size_t) b];

    if (gatedCount <= 0)
        return 0.0f;

    const int rank10 = (int) std::ceil (0.10 * (double) gatedCount);
    const int rank95 = (int) std::ceil (0.95 * (double) gatedCount);

    int cum = 0;
    bool haveP10 = false;
    float p10 = binToLufs (gateBin);
    float p95 = binToLufs (gateBin);

    for (int b = gateBin; b < histBins; ++b)
    {
        cum += stCount[(size_t) b];

        if (! haveP10 && cum >= rank10)
        {
            p10 = binToLufs (b);
            haveP10 = true;
        }

        if (cum >= rank95)
        {
            p95 = binToLufs (b);
            break;
        }
    }

    return juce::jmax (0.0f, p95 - p10);
}

//==============================================================================
// Rolling LRA (windowed) from short-term samples
//==============================================================================

float RunningLoudnessStats::computeRollingLra (float integratedForGate) const noexcept
{
    const int wantSeconds = rollingWindowSeconds.load (std::memory_order_relaxed);
    const int N = juce::jlimit (1, maxRollingSeconds, wantSeconds);

    const int available = rollingFilled;
    const int useN = juce::jmin (N, available);

    if (useN < 4)
        return 0.0f;

    const float gate = juce::jmax (absGateLufs, integratedForGate - 20.0f);

    // Gather gated values into a small temp buffer
    std::array<float, maxRollingSeconds> tmp {};
    int count = 0;

    // Start from newest and go backwards useN samples
    for (int k = 0; k < useN; ++k)
    {
        const int idx = (rollingPos - 1 - k + maxRollingSeconds) % maxRollingSeconds;
        const float v = rollingStLufs[(size_t) idx];

        if (v >= gate)
            tmp[(size_t) count++] = v;
    }

    if (count < 4)
        return 0.0f;

    std::sort (tmp.begin(), tmp.begin() + count);

    // Percentile ranks (match our histogram rank behavior)
    const int rank10 = (int) std::ceil (0.10 * (double) count);
    const int rank95 = (int) std::ceil (0.95 * (double) count);

    const int i10 = juce::jlimit (0, count - 1, rank10 - 1);
    const int i95 = juce::jlimit (0, count - 1, rank95 - 1);

    const float p10 = tmp[(size_t) i10];
    const float p95 = tmp[(size_t) i95];

    return juce::jmax (0.0f, p95 - p10);
}
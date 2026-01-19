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

    integratedLufs.store (-200.0f, std::memory_order_relaxed);
    lraLu.store (0.0f, std::memory_order_relaxed);
    lraGateLufs.store (-200.0f, std::memory_order_relaxed);
}

//==============================================================================
// Conversions (BS.1770)
// LUFS = -0.691 + 10*log10(meanSquare)
//==============================================================================

float RunningLoudnessStats::energyToLufs (float meanSquare) noexcept
{
    const double e = (double) meanSquare;
    if (e <= 0.0)
        return -200.0f;

    const double lufs = -0.691 + 10.0 * std::log10 (e);
    return (float) lufs;
}

float RunningLoudnessStats::lufsToEnergy (float lufs) noexcept
{
    // meanSquare = 10^((LUFS + 0.691)/10)
    const double e = std::pow (10.0, ((double) lufs + 0.691) / 10.0);
    return (float) e;
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
    // Counters + scheduled computations
    //--------------------------------------------------------------------------
    ++frameCounter;

    // 10 Hz: create 400ms block loudness every 100ms
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

    // 1 Hz: compute short-term loudness distribution
    if ((frameCounter % sStepFrames) == 0 && winSFilled >= sWindowFrames)
    {
        const float stMeanE = (float) (winSSum / (double) sWindowFrames);
        const float stLufs  = energyToLufs (stMeanE);

        // Absolute gate for LRA processing is also typically -70 LUFS
        if (stLufs >= absGateLufs)
        {
            const int bin = lufsToBin (stLufs);
            ++stCount[(size_t) bin];
            ++stTotalCount;
        }

        const float I = integratedLufs.load (std::memory_order_relaxed);
        const float LRA = computeLraFromHistogram (I);
        lraLu.store (LRA, std::memory_order_relaxed);
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

    // Compute energy/count above relative gate using histogram
    double sumE = 0.0;
    int count = 0;

    const int gateBin = lufsToBin (relGate);

    for (int b = gateBin; b < histBins; ++b)
    {
        count += intCount[(size_t) b];
        sumE  += intEnergySum[(size_t) b];
    }

    if (count <= 0 || sumE <= 0.0)
        return I_abs; // fallback

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

    // Count how many values survive the gate
    int gatedCount = 0;
    for (int b = gateBin; b < histBins; ++b)
        gatedCount += stCount[(size_t) b];

    if (gatedCount <= 0)
        return 0.0f;

    const int rank10 = (int) std::ceil (0.10 * (double) gatedCount);
    const int rank95 = (int) std::ceil (0.95 * (double) gatedCount);

    int cum = 0;
    float p10 = binToLufs (gateBin);
    float p95 = binToLufs (gateBin);

    for (int b = gateBin; b < histBins; ++b)
    {
        cum += stCount[(size_t) b];

        if (cum >= rank10 && p10 == binToLufs (gateBin))
            p10 = binToLufs (b);

        if (cum >= rank95)
        {
            p95 = binToLufs (b);
            break;
        }
    }

    return juce::jmax (0.0f, p95 - p10);
}
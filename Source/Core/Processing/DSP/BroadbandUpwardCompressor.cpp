#include "BroadbandUpwardCompressor.h"

namespace levelscope::dsp
{
void BroadbandUpwardCompressor::prepare (double sampleRate,
                                        int numChannels,
                                        const juce::AudioChannelSet& channelSet,
                                        int maxBlockSize)
{
    fs = (sampleRate > 0.0 ? sampleRate : 48000.0);
    preparedNumChannels = std::max (1, numChannels);
    preparedChannelSet = channelSet;
    preparedMaxBlockSize = std::max (0, maxBlockSize);

    // [BEGIN LS-BUC-PREPARE-BUILD-MASKS]
    detectChannels.clear();
    applyChannels.clear();
    detectChannels.reserve ((size_t) preparedNumChannels);
    applyChannels.reserve ((size_t) preparedNumChannels);

    for (int ch = 0; ch < preparedNumChannels; ++ch)
    {
        const bool isLfe = (preparedChannelSet.getTypeOfChannel (ch) == juce::AudioChannelSet::LFE);

        if (! isLfe)
        {
            detectChannels.push_back (ch);
            applyChannels.push_back (ch);
        }
    }

    // Safety fallback: if layout is weird and we excluded everything, include all channels.
    if (detectChannels.empty() || applyChannels.empty())
    {
        detectChannels.clear();
        applyChannels.clear();
        for (int ch = 0; ch < preparedNumChannels; ++ch)
        {
            detectChannels.push_back (ch);
            applyChannels.push_back (ch);
        }
    }
    // [END LS-BUC-PREPARE-BUILD-MASKS]

    // [BEGIN LS-BUC-PREPARE-BUILD-ALLCHANNELS]
    allChannels.clear();
    allChannels.reserve ((size_t) preparedNumChannels);
    for (int ch = 0; ch < preparedNumChannels; ++ch)
        allChannels.push_back (ch);
    // [END LS-BUC-PREPARE-BUILD-ALLCHANNELS]

    updateCoefficientsIfNeeded();
    reset();
}

void BroadbandUpwardCompressor::reset() noexcept
{
    envMS = 0.0f;
    gainZ = 1.0f;
}

void BroadbandUpwardCompressor::setParametersAudioThread (const Parameters& p) noexcept
{
    params = p;
    updateCoefficientsIfNeeded();
}

void BroadbandUpwardCompressor::updateCoefficientsIfNeeded() noexcept
{
    // Only recompute if attack/release changed meaningfully.
    if (params.attackMs == lastAttackMs && params.releaseMs == lastReleaseMs)
        return;

    lastAttackMs  = params.attackMs;
    lastReleaseMs = params.releaseMs;

    const float attackS  = std::max (1.0e-4f, params.attackMs  * 0.001f);
    const float releaseS = std::max (1.0e-4f, params.releaseMs * 0.001f);

    const float sr = (float) std::max (1.0, fs);

    // Standard per-sample one-pole coefficients
    aDetA  = std::exp (-1.0f / (attackS  * sr));
    aDetR  = std::exp (-1.0f / (releaseS * sr));

    // Use the same time constants for gain smoothing in this minimal version
    aGainA = aDetA;
    aGainR = aDetR;
}

void BroadbandUpwardCompressor::process (juce::AudioBuffer<float>& buffer) noexcept
{
    juce::ScopedNoDenormals noDenormals;

    const int numSamples = buffer.getNumSamples();
    const int numChInBuf = buffer.getNumChannels();
    if (numSamples <= 0 || numChInBuf <= 0)
        return;

    // [BEGIN LS-BUC-PROCESS-EXCLUDE-LFE]
    const int chToProcess = std::min (preparedNumChannels > 0 ? preparedNumChannels : numChInBuf, numChInBuf);

    const float t0 = std::min (params.t0Lufs, params.t1Lufs);
    const float t1 = std::max (params.t0Lufs, params.t1Lufs);

    const float userAmount = juce::jlimit (0.0f, 1.0f, params.amount01);
    const float maxBoostDb = std::max (0.0f, params.maxBoostDb);
    const float curve01    = juce::jlimit (0.0f, 1.0f, params.curve);

    const float expo = 1.0f + curve01 * 3.0f;
    const float range = std::max (1.0f, t1 - t0);

    float* const* chans = buffer.getArrayOfWritePointers();

    // [BEGIN LS-BUC-SELECT-DETECT-APPLY-LISTS]
    const auto& detectList = (params.lfeInDetector ? allChannels : detectChannels);
    const auto& applyList  = (params.lfeInApply    ? allChannels : applyChannels);

    const int numDetect = (int) detectList.size();
    const int numApply  = (int) applyList.size();
    // [END LS-BUC-SELECT-DETECT-APPLY-LISTS]

    for (int i = 0; i < numSamples; ++i)
    {
        // Linked broadband mean-square across DETECTOR channels (default excludes LFE).
        double sumSq = 0.0;

        for (int di = 0; di < numDetect; ++di)
        {
            const int ch = detectList[(size_t) di];
            if (ch >= 0 && ch < chToProcess)
            {
                const float x = chans[ch][i];
                sumSq += (double) x * (double) x;
            }
        }

        const float e = (float) (sumSq / (double) std::max (1, numDetect));

        // Detector smoothing (mean-square)
        const float aDet = (e > envMS ? aDetA : aDetR);
        envMS = aDet * envMS + (1.0f - aDet) * e;

        // Loudness proxy in LUFS-ish units (same constant as your other LUFS conversion)
        const float L = (float) (-0.691 + 10.0 * std::log10 ((double) envMS + 1.0e-12));

        // [BEGIN LS-BUC-T1T2-UNTOUCHED-ZONE]
        // Option A zone window with ONE-SIDED knees:
        // Fade OUT reaches 0 at/above T1 (no bleed into T1–T2).
        auto kneeUpToThreshold01 = [] (float levelDb, float threshold, float kneeWidthDb) noexcept
        {
            kneeWidthDb = juce::jmax (1.0e-4f, kneeWidthDb);

            const float start = threshold - kneeWidthDb;

            if (levelDb <= start)     return 0.0f;
            if (levelDb >= threshold) return 1.0f;

            const float tt = (levelDb - start) / kneeWidthDb; // 0..1
            return tt * tt * (3.0f - 2.0f * tt);
        };

        const float inAroundT0  = kneeUpToThreshold01 (L, t0, params.lowKneeDb);
        const float outAroundT1 = 1.0f - kneeUpToThreshold01 (L, t1, params.highKneeDb);

        const float zone01 = juce::jlimit (0.0f, 1.0f, inAroundT0 * outAroundT1);

        // Hard reset above T1 to prevent lingering gain smoothing into untouched zone.
        if (L >= t1)
            gainZ = 1.0f;
        // [END LS-BUC-T1T2-UNTOUCHED-ZONE]

        // Position inside zone (0..1), used for curve shaping
        float pos = (L - t0) / range;
        pos = juce::jlimit (0.0f, 1.0f, pos);

        float shaped = 0.0f;

        if (params.curveType == CurveType::monotonic)
        {
            shaped = std::pow (std::max (0.0f, 1.0f - pos), expo);
        }
        else // bell
        {
            const float d = std::abs (pos - 0.5f) * 2.0f; // 0..1
            shaped = std::pow (std::max (0.0f, 1.0f - d), expo);
        }

        // Boost in dB, bounded and zone-weighted
        float boostDb = maxBoostDb * shaped * zone01;
        boostDb = juce::jlimit (0.0f, maxBoostDb, boostDb);

        const float baseGain = dbToLin (boostDb);

        // Apply user amount (fade toward unity)
        const float gainTarget = 1.0f + (baseGain - 1.0f) * userAmount;

        // Gain smoothing
        const float aG = (gainTarget > gainZ ? aGainA : aGainR);
        gainZ = aG * gainZ + (1.0f - aG) * gainTarget;

        // Apply to APPLY channels (default excludes LFE)
        for (int ai = 0; ai < numApply; ++ai)
        {
            const int ch = applyList[(size_t) ai];
            if (ch >= 0 && ch < chToProcess)
                chans[ch][i] *= gainZ;
        }
    }
// [END LS-BUC-PROCESS-EXCLUDE-LFE]
}
} // namespace levelscope::dsp
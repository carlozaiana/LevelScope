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

    // [BEGIN LS-BUC-STAGE-E-PREPARE-MASKBITS]
    preparedAllMaskBits    = 0;
    preparedNonLfeMaskBits = 0;

    const int nMaskCh = juce::jlimit (0, kMaxMaskChannels, preparedNumChannels);

    for (int ch = 0; ch < nMaskCh; ++ch)
    {
        const uint16_t b = (uint16_t) (1u << (unsigned) ch);
        preparedAllMaskBits |= b;

        const bool isLfe = (preparedChannelSet.getTypeOfChannel (ch) == juce::AudioChannelSet::LFE);
        if (! isLfe)
            preparedNonLfeMaskBits |= b;
    }

    if (preparedNonLfeMaskBits == 0)
        preparedNonLfeMaskBits = preparedAllMaskBits;

    externalMasksActive = false;
    externalDetectMaskBits = 0;
    externalApplyMaskBits  = 0;
    externalUnlinked       = false;

    effectiveDetectMaskBitsCached = 0;
    effectiveApplyMaskBitsCached  = 0;
    detectCount = applyCount = 0;
    // [END LS-BUC-STAGE-E-PREPARE-MASKBITS]

    updateCoefficientsIfNeeded();
    reset();
}

void BroadbandUpwardCompressor::reset() noexcept
{
    envMS = 0.0f;
    gainZ = 1.0f;
    // [BEGIN LS-BUC-STAGE-E-RESET-UNLINKED]
    envMSUnlinked.fill (0.0f);
    gainZUnlinked.fill (1.0f);
    // [END LS-BUC-STAGE-E-RESET-UNLINKED]
}

void BroadbandUpwardCompressor::setParametersAudioThread (const Parameters& p) noexcept
{
    params = p;
    updateCoefficientsIfNeeded();
}

// [BEGIN LS-BUC-STAGE-E-MASK-IMPL]
void BroadbandUpwardCompressor::setChannelMasksAudioThread (uint16_t detectBits,
                                                            uint16_t applyBits,
                                                            bool unlinked) noexcept
{
    if (detectBits == 0 && applyBits == 0)
    {
        externalMasksActive = false;
        externalDetectMaskBits = 0;
        externalApplyMaskBits  = 0;
        externalUnlinked       = false;
        return;
    }

    externalMasksActive = true;
    externalDetectMaskBits = detectBits;
    externalApplyMaskBits  = applyBits;
    externalUnlinked       = unlinked;
}

void BroadbandUpwardCompressor::rebuildIndexListsIfNeededNoAlloc (uint16_t detBits,
                                                                  uint16_t appBits) noexcept
{
    if (detBits == effectiveDetectMaskBitsCached && appBits == effectiveApplyMaskBitsCached)
        return;

    effectiveDetectMaskBitsCached = detBits;
    effectiveApplyMaskBitsCached  = appBits;

    const int nMaskCh = juce::jlimit (0, kMaxMaskChannels, preparedNumChannels);

    detectCount = 0;
    applyCount  = 0;

    for (int ch = 0; ch < nMaskCh; ++ch)
    {
        const uint16_t b = (uint16_t) (1u << (unsigned) ch);

        if ((detBits & b) != 0 && detectCount < kMaxMaskChannels)
            detectIdx[(size_t) detectCount++] = (uint8_t) ch;

        if ((appBits & b) != 0 && applyCount < kMaxMaskChannels)
            applyIdx[(size_t) applyCount++] = (uint8_t) ch;
    }

    if (detectCount <= 0)
        for (int ch = 0; ch < nMaskCh && detectCount < kMaxMaskChannels; ++ch)
            detectIdx[(size_t) detectCount++] = (uint8_t) ch;

    if (applyCount <= 0)
        for (int ch = 0; ch < nMaskCh && applyCount < kMaxMaskChannels; ++ch)
            applyIdx[(size_t) applyCount++] = (uint8_t) ch;
}
// [END LS-BUC-STAGE-E-MASK-IMPL]

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

    // [BEGIN LS-BUC-STAGE-E-PROCESS-MASKS-AND-UNLINKED]
    const int chToProcess = std::min (preparedNumChannels > 0 ? preparedNumChannels : numChInBuf, numChInBuf);

    const float t0 = std::min (params.t0Lufs, params.t1Lufs);
    const float t1 = std::max (params.t0Lufs, params.t1Lufs);

    const float userAmount = juce::jlimit (0.0f, 1.0f, params.amount01);
    const float maxBoostDb = std::max (0.0f, params.maxBoostDb);
    const float curve01    = juce::jlimit (0.0f, 1.0f, params.curve);

    const float expo  = 1.0f + curve01 * 3.0f;
    const float range = std::max (1.0f, t1 - t0);

    float* const* chans = buffer.getArrayOfWritePointers();

    const uint16_t legacyDetectBits = (params.lfeInDetector ? preparedAllMaskBits : preparedNonLfeMaskBits);
    const uint16_t legacyApplyBits  = (params.lfeInApply    ? preparedAllMaskBits : preparedNonLfeMaskBits);

    uint16_t effDetectBits = (externalMasksActive ? externalDetectMaskBits : legacyDetectBits);
    uint16_t effApplyBits  = (externalMasksActive ? externalApplyMaskBits  : legacyApplyBits);

    effDetectBits = (uint16_t) (effDetectBits & preparedAllMaskBits);
    effApplyBits  = (uint16_t) (effApplyBits  & preparedAllMaskBits);
    if (effDetectBits == 0) effDetectBits = preparedAllMaskBits;
    if (effApplyBits  == 0) effApplyBits  = preparedAllMaskBits;

    rebuildIndexListsIfNeededNoAlloc (effDetectBits, effApplyBits);

    const bool doUnlinked = (externalMasksActive && externalUnlinked);

    auto kneeUpToThreshold01 = [] (float levelDb, float threshold, float kneeWidthDb) noexcept
    {
        kneeWidthDb = juce::jmax (1.0e-4f, kneeWidthDb);
        const float start = threshold - kneeWidthDb;
        if (levelDb <= start)     return 0.0f;
        if (levelDb >= threshold) return 1.0f;
        const float tt = (levelDb - start) / kneeWidthDb;
        return tt * tt * (3.0f - 2.0f * tt);
    };

    if (! doUnlinked)
    {
        for (int i = 0; i < numSamples; ++i)
        {
            // Linked detector across detectIdx
            double sumSq = 0.0;
            for (int di = 0; di < detectCount; ++di)
            {
                const int ch = (int) detectIdx[(size_t) di];
                if (ch >= 0 && ch < chToProcess)
                {
                    const float x = chans[ch][i];
                    sumSq += (double) x * (double) x;
                }
            }

            const float e = (float) (sumSq / (double) std::max (1, detectCount));

            const float aDet = (e > envMS ? aDetA : aDetR);
            envMS = aDet * envMS + (1.0f - aDet) * e;

            const float L = (float) (-0.691 + 10.0 * std::log10 ((double) envMS + 1.0e-12));

            const float inAroundT0  = kneeUpToThreshold01 (L, t0, params.lowKneeDb);
            const float outAroundT1 = 1.0f - kneeUpToThreshold01 (L, t1, params.highKneeDb);
            const float zone01      = juce::jlimit (0.0f, 1.0f, inAroundT0 * outAroundT1);

            if (L >= t1)
                gainZ = 1.0f;

            float pos = (L - t0) / range;
            pos = juce::jlimit (0.0f, 1.0f, pos);

            float shaped = 0.0f;
            if (params.curveType == CurveType::monotonic)
                shaped = std::pow (std::max (0.0f, 1.0f - pos), expo);
            else
            {
                const float d = std::abs (pos - 0.5f) * 2.0f;
                shaped = std::pow (std::max (0.0f, 1.0f - d), expo);
            }

            float boostDb = maxBoostDb * shaped * zone01;
            boostDb = juce::jlimit (0.0f, maxBoostDb, boostDb);

            const float baseGain = dbToLin (boostDb);
            const float gainTarget = 1.0f + (baseGain - 1.0f) * userAmount;

            const float aG = (gainTarget > gainZ ? aGainA : aGainR);
            gainZ = aG * gainZ + (1.0f - aG) * gainTarget;

            for (int ai = 0; ai < applyCount; ++ai)
            {
                const int ch = (int) applyIdx[(size_t) ai];
                if (ch >= 0 && ch < chToProcess)
                    chans[ch][i] *= gainZ;
            }
        }
    }
    else
    {
        // Unlinked: per-channel detectors/gains (apply channels only)
        for (int i = 0; i < numSamples; ++i)
        {
            for (int ai = 0; ai < applyCount; ++ai)
            {
                const int ch = (int) applyIdx[(size_t) ai];
                if (ch < 0 || ch >= chToProcess || ch >= kMaxMaskChannels)
                    continue;

                const float x = chans[ch][i];
                const float e = x * x;

                float& env = envMSUnlinked[(size_t) ch];
                float& gz  = gainZUnlinked[(size_t) ch];

                const float aDet = (e > env ? aDetA : aDetR);
                env = aDet * env + (1.0f - aDet) * e;

                const float L = (float) (-0.691 + 10.0 * std::log10 ((double) env + 1.0e-12));

                const float inAroundT0  = kneeUpToThreshold01 (L, t0, params.lowKneeDb);
                const float outAroundT1 = 1.0f - kneeUpToThreshold01 (L, t1, params.highKneeDb);
                const float zone01      = juce::jlimit (0.0f, 1.0f, inAroundT0 * outAroundT1);

                if (L >= t1)
                    gz = 1.0f;

                float pos = (L - t0) / range;
                pos = juce::jlimit (0.0f, 1.0f, pos);

                float shaped = 0.0f;
                if (params.curveType == CurveType::monotonic)
                    shaped = std::pow (std::max (0.0f, 1.0f - pos), expo);
                else
                {
                    const float d = std::abs (pos - 0.5f) * 2.0f;
                    shaped = std::pow (std::max (0.0f, 1.0f - d), expo);
                }

                float boostDb = maxBoostDb * shaped * zone01;
                boostDb = juce::jlimit (0.0f, maxBoostDb, boostDb);

                const float baseGain = dbToLin (boostDb);
                const float gainTarget = 1.0f + (baseGain - 1.0f) * userAmount;

                const float aG = (gainTarget > gz ? aGainA : aGainR);
                gz = aG * gz + (1.0f - aG) * gainTarget;

                chans[ch][i] *= gz;
            }
        }
    }
    // [END LS-BUC-STAGE-E-PROCESS-MASKS-AND-UNLINKED]
}
} // namespace levelscope::dsp
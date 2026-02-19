#include "BroadbandDownwardCompressor.h"

namespace levelscope::dsp
{
void BroadbandDownwardCompressor::prepare (double sampleRate,
                                          int numChannels,
                                          const juce::AudioChannelSet& channelSet,
                                          int maxBlockSize)
{
    fs = (sampleRate > 0.0 ? sampleRate : 48000.0);
    preparedNumChannels = std::max (1, numChannels);
    preparedChannelSet = channelSet;
    preparedMaxBlockSize = std::max (0, maxBlockSize);

    // [BEGIN LS-BDC-STAGE-E-PREPARE-MASKBITS]
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

    legacyDetectMaskBits = preparedNonLfeMaskBits;
    legacyApplyMaskBits  = preparedNonLfeMaskBits;

    externalMasksActive = false;
    externalDetectMaskBits = 0;
    externalApplyMaskBits  = 0;
    externalUnlinked       = false;

    effectiveDetectMaskBitsCached = 0;
    effectiveApplyMaskBitsCached  = 0;
    detectCount = applyCount = 0;
    // [END LS-BDC-STAGE-E-PREPARE-MASKBITS]

    updateCoefficientsIfNeeded();
    reset();
}

void BroadbandDownwardCompressor::reset() noexcept
{
    envMS = 0.0f;
    gainZ = 1.0f;
    // [BEGIN LS-BDC-STAGE-E-RESET-UNLINKED]
    envMSUnlinked.fill (0.0f);
    gainZUnlinked.fill (1.0f);
    // [END LS-BDC-STAGE-E-RESET-UNLINKED]
}

void BroadbandDownwardCompressor::setParametersAudioThread (const Parameters& p) noexcept
{
    params = p;
    updateCoefficientsIfNeeded();
}

// [BEGIN LS-BDC-STAGE-E-MASK-IMPL]
void BroadbandDownwardCompressor::setChannelMasksAudioThread (uint16_t detectBits,
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

void BroadbandDownwardCompressor::rebuildIndexListsIfNeededNoAlloc (uint16_t detBits,
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
// [END LS-BDC-STAGE-E-MASK-IMPL]

// [BEGIN LS-BDC-STAGE-E-SETCHANNELLISTS-LEGACYBITS]
void BroadbandDownwardCompressor::setChannelListsAudioThread (const std::vector<int>* detectList,
                                                             const std::vector<int>* applyList) noexcept
{
    auto listToBits = [this] (const std::vector<int>* list) noexcept -> uint16_t
    {
        if (list == nullptr || list->empty())
            return preparedNonLfeMaskBits;

        uint16_t bits = 0;
        for (int idx : *list)
        {
            if (idx >= 0 && idx < kMaxMaskChannels)
                bits |= (uint16_t) (1u << (unsigned) idx);
        }

        bits = (uint16_t) (bits & preparedAllMaskBits);
        if (bits == 0) bits = preparedAllMaskBits;
        return bits;
    };

    legacyDetectMaskBits = listToBits (detectList);
    legacyApplyMaskBits  = listToBits (applyList);
}
// [END LS-BDC-STAGE-E-SETCHANNELLISTS-LEGACYBITS]

// [BEGIN LS-BDC-STAGE-E-LFE-POLICY-LEGACYBITS]
void BroadbandDownwardCompressor::setLfePolicyAudioThread (bool lfeInDetector, bool lfeInApply) noexcept
{
    legacyDetectMaskBits = (lfeInDetector ? preparedAllMaskBits : preparedNonLfeMaskBits);
    legacyApplyMaskBits  = (lfeInApply    ? preparedAllMaskBits : preparedNonLfeMaskBits);
}
// [END LS-BDC-STAGE-E-LFE-POLICY-LEGACYBITS]

void BroadbandDownwardCompressor::updateCoefficientsIfNeeded() noexcept
{
    if (params.attackMs == lastAttackMs && params.releaseMs == lastReleaseMs)
        return;

    lastAttackMs  = params.attackMs;
    lastReleaseMs = params.releaseMs;

    const float attackS  = std::max (1.0e-4f, params.attackMs  * 0.001f);
    const float releaseS = std::max (1.0e-4f, params.releaseMs * 0.001f);

    const float sr = (float) std::max (1.0, fs);

    aA = std::exp (-1.0f / (attackS  * sr));
    aR = std::exp (-1.0f / (releaseS * sr));
}

float BroadbandDownwardCompressor::computeCompressionGain (float levelLufs) const noexcept
{
    const float ratio = std::max (1.0f, params.ratio);
    const float thr   = params.t2Lufs;
    const float knee  = std::max (0.0f, params.kneeDb);

    if (ratio <= 1.0001f)
        return 1.0f;

    const float x = levelLufs - thr;

    // Soft knee around threshold (width = kneeDb, symmetric)
    float overDb = 0.0f;

    if (knee <= 1.0e-3f)
    {
        overDb = std::max (0.0f, x);
    }
    else
    {
        const float half = 0.5f * knee;

        if (x <= -half)       overDb = 0.0f;
        else if (x >= half)   overDb = x;
        else
        {
            // quadratic ease-in within knee region
            const float t = (x + half) / knee; // 0..1
            overDb = t * t * half;             // 0..half (smooth)
        }
    }

    if (overDb <= 0.0f)
        return 1.0f;

    // Gain reduction dB = overDb * (1 - 1/ratio)
    const float reductionDb = overDb * (1.0f - 1.0f / ratio);
    const float g = dbToLin (-reductionDb);
    return juce::jlimit (0.0f, 1.0f, g);
}

void BroadbandDownwardCompressor::process (juce::AudioBuffer<float>& buffer) noexcept
{
    juce::ScopedNoDenormals noDenormals;

    if (! params.enabled)
        return;

    const int numSamples = buffer.getNumSamples();
    const int numChInBuf = buffer.getNumChannels();
    if (numSamples <= 0 || numChInBuf <= 0)
        return;

    float* const* chans = buffer.getArrayOfWritePointers();

    // [BEGIN LS-BDC-STAGE-E-PROCESS-MASKS]
    const int chToProcess = std::min (preparedNumChannels > 0 ? preparedNumChannels : numChInBuf, numChInBuf);

    uint16_t effDetectBits = (externalMasksActive ? externalDetectMaskBits : legacyDetectMaskBits);
    uint16_t effApplyBits  = (externalMasksActive ? externalApplyMaskBits  : legacyApplyMaskBits);

    effDetectBits = (uint16_t) (effDetectBits & preparedAllMaskBits);
    effApplyBits  = (uint16_t) (effApplyBits  & preparedAllMaskBits);

    if (effDetectBits == 0) effDetectBits = preparedAllMaskBits;
    if (effApplyBits  == 0) effApplyBits  = preparedAllMaskBits;

    rebuildIndexListsIfNeededNoAlloc (effDetectBits, effApplyBits);

    const bool doUnlinked = (externalMasksActive && externalUnlinked);
    // [END LS-BDC-STAGE-E-PROCESS-MASKS]

    const float t2 = std::min (params.t2Lufs, params.t3Lufs);
    const float t3 = std::max (params.t2Lufs, params.t3Lufs);
    const float zoneWidth = std::max (1.0f, t3 - t2);

    const float makeupLin = dbToLin (params.makeupDb);

    // [BEGIN LS-BDC-METERING-INIT]
    float blockMinG  = 1.0f;
    float blockLastG = 1.0f;
    // [END LS-BDC-METERING-INIT]

    // [BEGIN LS-BDC-STAGE-E-LINKED-UNLINKED-LOOPS]
    if (! doUnlinked)
    {
        for (int i = 0; i < numSamples; ++i)
        {
            double sumSq = 0.0;
            for (int di = 0; di < detectCount; ++di)
            {
                const int chDet = (int) detectIdx[(size_t) di];
                if (chDet >= 0 && chDet < chToProcess)
                {
                    const float x = chans[chDet][i];
                    sumSq += (double) x * (double) x;
                }
            }

            const float e = (float) (sumSq / (double) std::max (1, detectCount));

            const float a = (e > envMS ? aA : aR);
            envMS = a * envMS + (1.0f - a) * e;

            const float L = (float) (-0.691 + 10.0 * std::log10 ((double) envMS + 1.0e-12));

            float pos = (L - t2) / zoneWidth;
            const float zone01 = smoothstep01 (pos);

            const float gComp = computeCompressionGain (L);
            const float gTarget = 1.0f + (gComp - 1.0f) * zone01;

            const float aG = (gTarget < gainZ ? aA : aR);
            gainZ = aG * gainZ + (1.0f - aG) * gTarget;

            const float gOut = gainZ * makeupLin;

            blockMinG  = std::min (blockMinG, gainZ);
            blockLastG = gainZ;

            for (int ai = 0; ai < applyCount; ++ai)
            {
                const int chAp = (int) applyIdx[(size_t) ai];
                if (chAp >= 0 && chAp < chToProcess)
                    chans[chAp][i] *= gOut;
            }
        }
    }
    else
    {
        // Unlinked: per-channel detectors/gains, applied per channel (apply channels)
        float lastMinAtEnd = 1.0f;

        for (int i = 0; i < numSamples; ++i)
        {
            float minGainThisSample = 1.0f;

            for (int ai = 0; ai < applyCount; ++ai)
            {
                const int chAp = (int) applyIdx[(size_t) ai];
                if (chAp < 0 || chAp >= chToProcess || chAp >= kMaxMaskChannels)
                    continue;

                const float x = chans[chAp][i];
                const float e = x * x;

                float& env = envMSUnlinked[(size_t) chAp];
                float& gz  = gainZUnlinked[(size_t) chAp];

                const float a = (e > env ? aA : aR);
                env = a * env + (1.0f - a) * e;

                const float L = (float) (-0.691 + 10.0 * std::log10 ((double) env + 1.0e-12));

                float pos = (L - t2) / zoneWidth;
                const float zone01 = smoothstep01 (pos);

                const float gComp = computeCompressionGain (L);
                const float gTarget = 1.0f + (gComp - 1.0f) * zone01;

                const float aG = (gTarget < gz ? aA : aR);
                gz = aG * gz + (1.0f - aG) * gTarget;

                const float gOut = gz * makeupLin;

                chans[chAp][i] *= gOut;

                minGainThisSample = std::min (minGainThisSample, gz);
                blockMinG = std::min (blockMinG, gz);
            }

            lastMinAtEnd = minGainThisSample;
        }

        blockLastG = lastMinAtEnd;
    }
    // [END LS-BDC-STAGE-E-LINKED-UNLINKED-LOOPS]
    // [BEGIN LS-BDC-METERING-STORE]
    lastBlockMinCompGain  = blockMinG;
    lastBlockLastCompGain = blockLastG;
    // [END LS-BDC-METERING-STORE]
}
} // namespace levelscope::dsp
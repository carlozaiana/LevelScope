#include "BroadbandUpwardCompressor.h"
// [BEGIN LS-BUC-UPWARD-GAINLAW-INCLUDE]
#include "../DSP/UpwardGainLaw.h"
// [END LS-BUC-UPWARD-GAINLAW-INCLUDE]

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
    // [BEGIN LS-BUC-MOMENTARY-DETECTOR-PREPARE]
    momentaryDetector.prepare (fs, preparedNumChannels);
    cachedMomentaryLufsLinked = -200.0f;
    cachedMomentaryLufsByChannel.fill (-200.0f);
    // [END LS-BUC-MOMENTARY-DETECTOR-PREPARE]
    reset();
}

void BroadbandUpwardCompressor::reset() noexcept
{
    envMS = 0.0f;
    gainZ = 1.0f;
    // [BEGIN LS-BUC-MOMENTARY-DETECTOR-RESET]
    momentaryDetector.reset();
    cachedMomentaryLufsLinked = -200.0f;
    cachedMomentaryLufsByChannel.fill (-200.0f);
    // [END LS-BUC-MOMENTARY-DETECTOR-RESET]
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

    // [BEGIN LS-BUC-UPWARD-METERING-BLOCK-INIT]
    float blockMaxG  = 1.0f;
    float blockLastG = 1.0f;
    // [END LS-BUC-UPWARD-METERING-BLOCK-INIT]

    const int chToProcess = std::min (preparedNumChannels > 0 ? preparedNumChannels : numChInBuf, numChInBuf);

    const float t0 = std::min (params.t0Lufs, params.t1Lufs);
    const float t1 = std::max (params.t0Lufs, params.t1Lufs);

    const float userAmount = juce::jlimit (0.0f, 1.0f, params.amount01);
    const float maxBoostDb = std::max (0.0f, params.maxBoostDb);
    const float curve01    = juce::jlimit (0.0f, 1.0f, params.curve);

    const float expo  = 1.0f + curve01 * 3.0f;
    // [BEGIN LS-BUC-RELATIVE-GUARD-CONSTANTS]
    // Relative guard: prevents upward boost on silence / very low noise floor.
    // Effective guard floor is (T0 - kGuardBelowT0Db), with a fade-in over kGuardFadeDb.
    static constexpr float kGuardBelowT0Db = 24.0f;
    static constexpr float kGuardFadeDb    = 6.0f;
    // [END LS-BUC-RELATIVE-GUARD-CONSTANTS]
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

    if (! doUnlinked)
    {
        // [BEGIN LS-BUC-MOMENTARY-LINKED]
        // Initialize cached control loudness at block start (will update at 60 Hz frames).
        cachedMomentaryLufsLinked = momentaryDetector.getMomentaryLufsForMask (effDetectBits, chToProcess);

        UpwardGainLaw::Params law;
        law.t0Db       = t0;
        law.t1Db       = t1;
        law.lowKneeDb  = params.lowKneeDb;
        law.highKneeDb = params.highKneeDb;
        law.maxBoostDb = maxBoostDb;
        law.curve01    = curve01;
        law.curveType  = (params.curveType == CurveType::bell ? UpwardGainLaw::CurveType::bell
                                                              : UpwardGainLaw::CurveType::monotonic);

        const float guardFloor = t0 - kGuardBelowT0Db;

        for (int i = 0; i < numSamples; ++i)
        {
            // 1) Update K-weighted detector state (all channels)
            for (int ch = 0; ch < chToProcess && ch < kMaxMaskChannels; ++ch)
                kWeightedSampleScratch[(size_t) ch] = momentaryDetector.processSample (ch, chans[ch][i]);

            // Advance the detector clock once per sample (60 Hz frames)
            if (momentaryDetector.advanceSample())
                cachedMomentaryLufsLinked = momentaryDetector.getMomentaryLufsForMask (effDetectBits, chToProcess);

            // 2) Instantaneous guard level from K-weighted instantaneous energy (detect channels)
            double sumSqInst = 0.0;
            for (int di = 0; di < detectCount; ++di)
            {
                const int ch = (int) detectIdx[(size_t) di];
                if (ch >= 0 && ch < chToProcess && ch < kMaxMaskChannels)
                {
                    const float y = kWeightedSampleScratch[(size_t) ch];
                    sumSqInst += (double) y * (double) y;
                }
            }

            const float eInst = (float) (sumSqInst / (double) std::max (1, detectCount));
            const float Linst = (float) (-0.691 + 10.0 * std::log10 ((double) eInst + 1.0e-12));

            const float guard01 =
                UpwardGainLaw::kneeUpToThreshold01 (Linst, guardFloor + kGuardFadeDb, kGuardFadeDb);

            // 3) Control loudness for gain-law (Momentary LUFS)
            const float L = cachedMomentaryLufsLinked;

            const float baseGain = UpwardGainLaw::computeBaseGainLin (L, law);

            const float effectiveAmount = userAmount * guard01;
            const float gainTarget = UpwardGainLaw::applyAmountLin (baseGain, effectiveAmount);

            // Smooth return to unity above T1 based on control loudness
            const bool aboveT1 = (L >= t1);

            float aG = 0.0f;
            if (gainTarget > gainZ)
                aG = aGainA;
            else
                aG = (aboveT1 ? aGainA : aGainR);

            gainZ = aG * gainZ + (1.0f - aG) * gainTarget;

            blockMaxG  = std::max (blockMaxG, gainZ);
            blockLastG = gainZ;

            for (int ai = 0; ai < applyCount; ++ai)
            {
                const int ch = (int) applyIdx[(size_t) ai];
                if (ch >= 0 && ch < chToProcess)
                    chans[ch][i] *= gainZ;
            }
        }
        // [END LS-BUC-MOMENTARY-LINKED]
    }
    else
    {
        // [BEGIN LS-BUC-MOMENTARY-UNLINKED]
        // Initialize per-channel cached control loudness
        for (int ch = 0; ch < chToProcess && ch < kMaxMaskChannels; ++ch)
            cachedMomentaryLufsByChannel[(size_t) ch] = momentaryDetector.getMomentaryLufsForChannel (ch);

        UpwardGainLaw::Params law;
        law.t0Db       = t0;
        law.t1Db       = t1;
        law.lowKneeDb  = params.lowKneeDb;
        law.highKneeDb = params.highKneeDb;
        law.maxBoostDb = maxBoostDb;
        law.curve01    = curve01;
        law.curveType  = (params.curveType == CurveType::bell ? UpwardGainLaw::CurveType::bell
                                                              : UpwardGainLaw::CurveType::monotonic);

        const float guardFloor = t0 - kGuardBelowT0Db;

        for (int i = 0; i < numSamples; ++i)
        {
            float sampleMaxG = 1.0f;

            // Update detector state for all channels
            for (int ch = 0; ch < chToProcess && ch < kMaxMaskChannels; ++ch)
                kWeightedSampleScratch[(size_t) ch] = momentaryDetector.processSample (ch, chans[ch][i]);

            if (momentaryDetector.advanceSample())
            {
                for (int ch = 0; ch < chToProcess && ch < kMaxMaskChannels; ++ch)
                    cachedMomentaryLufsByChannel[(size_t) ch] = momentaryDetector.getMomentaryLufsForChannel (ch);
            }

            for (int ai = 0; ai < applyCount; ++ai)
            {
                const int ch = (int) applyIdx[(size_t) ai];
                if (ch < 0 || ch >= chToProcess || ch >= kMaxMaskChannels)
                    continue;

                const float y = kWeightedSampleScratch[(size_t) ch];
                const float eInst = y * y;
                const float Linst = (float) (-0.691 + 10.0 * std::log10 ((double) eInst + 1.0e-12));

                const float guard01 =
                    UpwardGainLaw::kneeUpToThreshold01 (Linst, guardFloor + kGuardFadeDb, kGuardFadeDb);

                const float L = cachedMomentaryLufsByChannel[(size_t) ch];

                const float baseGain = UpwardGainLaw::computeBaseGainLin (L, law);

                const float effectiveAmount = userAmount * guard01;
                const float gainTarget = UpwardGainLaw::applyAmountLin (baseGain, effectiveAmount);

                const bool aboveT1 = (L >= t1);

                float& gz = gainZUnlinked[(size_t) ch];

                float aG = 0.0f;
                if (gainTarget > gz)
                    aG = aGainA;
                else
                    aG = (aboveT1 ? aGainA : aGainR);

                gz = aG * gz + (1.0f - aG) * gainTarget;

                sampleMaxG = std::max (sampleMaxG, gz);
                chans[ch][i] *= gz;
            }

            blockMaxG  = std::max (blockMaxG, sampleMaxG);
            blockLastG = sampleMaxG;
        }
        // [END LS-BUC-MOMENTARY-UNLINKED]
    }
    // [BEGIN LS-BUC-UPWARD-METERING-BLOCK-STORE]
    lastBlockMaxLinearGain  = std::max (1.0f, blockMaxG);
    lastBlockLastLinearGain = std::max (1.0f, blockLastG);
    // [END LS-BUC-UPWARD-METERING-BLOCK-STORE]
    // [END LS-BUC-STAGE-E-PROCESS-MASKS-AND-UNLINKED]
}
} // namespace levelscope::dsp
#include "LevelerModule.h"

#include <cmath>
#include <algorithm>

namespace
{
    static float energyToLufsSafe (float energy) noexcept
    {
        if (! std::isfinite (energy) || energy <= 0.0f)
            return -200.0f;

        // BS.1770-oriented loudness conversion from K-weighted mean-square energy.
        return -0.691f + 10.0f * std::log10 (energy);
    }

    static float clampFinite (float v, float lo, float hi, float fallback) noexcept
    {
        if (! std::isfinite (v))
            return fallback;

        return juce::jlimit (lo, hi, v);
    }
}

namespace levelscope
{
    // [BEGIN LVLR-MODULE-IMPL]
    LevelerModule::LevelerModule()
    {
        resetDetectorState();
        resetGainState();
    }

    // [BEGIN LVLR-SAVESTATE-SIGNATURE-FIX]
    // [BEGIN LVLR-VALUETREE-PERSISTENCE-IMPL]
    juce::ValueTree LevelerModule::saveState() const
    {
        juce::ValueTree t (getModuleID());
        t.setProperty ("bypassed", isBypassed(), nullptr);
        return t;
    }

    void LevelerModule::loadState (const juce::ValueTree& state)
    {
        const auto b = state.getProperty ("bypassed", false);
        setBypassed ((bool) b);
    }
    // [END LVLR-VALUETREE-PERSISTENCE-IMPL]

    void LevelerModule::bindParameters (std::atomic<float>* enabledParam,
                                        std::atomic<float>* targetLufsParam,
                                        std::atomic<float>* maxBoostDbParam,
                                        std::atomic<float>* maxCutDbParam,
                                        std::atomic<float>* measChoiceParam,
                                        std::atomic<float>* modeChoiceParam,
                                        std::atomic<float>* learn01Param,
                                        std::atomic<float>* rateUpDbPerSecParam,
                                        std::atomic<float>* rateDownDbPerSecParam,
                                        std::atomic<float>* controlModeChoiceParam,
                                        std::atomic<float>* hostGainDbParam,
                                        std::atomic<float>* mcPolicyChoiceParam,
                                        std::atomic<float>* dialogDetectorChoiceParam,
                                        std::atomic<float>* dialogApplyChoiceParam,
                                        std::atomic<float>* lfeInDetectorParam,
                                        std::atomic<float>* lfeInApplyParam) noexcept
    {
        pEnabled01            = enabledParam;
        pTargetLufs           = targetLufsParam;
        pMaxBoostDb           = maxBoostDbParam;
        pMaxCutDb             = maxCutDbParam;
        pMeasChoice           = measChoiceParam;
        pModeChoice           = modeChoiceParam;
        pLearn01              = learn01Param;
        pRateUpDbPerSec       = rateUpDbPerSecParam;
        pRateDownDbPerSec     = rateDownDbPerSecParam;
        pControlModeChoice    = controlModeChoiceParam;
        pHostGainDb           = hostGainDbParam;

        pMcPolicyChoice       = mcPolicyChoiceParam;
        pDialogDetectorChoice = dialogDetectorChoiceParam;
        pDialogApplyChoice    = dialogApplyChoiceParam;
        pLfeInDetector        = lfeInDetectorParam;
        pLfeInApply           = lfeInApplyParam;

        routingMasksDirty = true;
    }

    void LevelerModule::prepare (const ModulePrepareSpec& spec)
    {
        sampleRate = (spec.sampleRate > 0.0 ? spec.sampleRate : 48000.0);
        preparedChannelSet = spec.channelSet;
        preparedNumChannels = juce::jlimit (0, maxSupportedChannels, (int) spec.channelSet.size());

        frameSamples = juce::jmax (1, (int) std::lround (sampleRate / detectorFrameRateHz));
        samplesUntilNextFrame = frameSamples;

        kWeight.prepare (sampleRate, juce::jmax (1, preparedNumChannels));

        routingMasksDirty = true;
        cachedMcPolicyChoice       = std::numeric_limits<int>::min();
        cachedDialogDetectorChoice = std::numeric_limits<int>::min();
        cachedDialogApplyChoice    = std::numeric_limits<int>::min();
        cachedLfeInDetector01      = std::numeric_limits<int>::min();
        cachedLfeInApply01         = std::numeric_limits<int>::min();

        resetDetectorState();
        resetGainState();
    }

    void LevelerModule::reset()
    {
        resetDetectorState();
        resetGainState();
    }

    float LevelerModule::readParamOrDefault (const std::atomic<float>* p, float fallback) noexcept
    {
        return (p != nullptr ? p->load (std::memory_order_relaxed) : fallback);
    }

    int LevelerModule::readChoiceOrDefault (const std::atomic<float>* p, int fallback) noexcept
    {
        if (p == nullptr)
            return fallback;

        return (int) std::lround (p->load (std::memory_order_relaxed));
    }

    uint16_t LevelerModule::bitForChannel (int ch) noexcept
    {
        return (ch >= 0 && ch < maxSupportedChannels ? (uint16_t) (1u << ch) : (uint16_t) 0u);
    }

    void LevelerModule::resetDetectorState() noexcept
    {
        kWeight.reset();

        for (auto& v : frameEnergyAccumByChannel)      v = 0.0;
        for (auto& v : momentaryEnergySumsByChannel)   v = 0.0f;
        for (auto& v : shortTermEnergySumsByChannel)   v = 0.0f;

        for (auto& chHist : frameEnergyHistoryByChannel)
            for (auto& e : chHist)
                e = 0.0f;

        detectorHistoryWriteIndex = 0;
        detectorFramesWritten = 0;

        latestMomentaryLufs = -200.0f;
        latestShortTermLufs = -200.0f;
        momentaryValid = false;
        shortTermValid = false;

        samplesUntilNextFrame = frameSamples;
    }

    void LevelerModule::resetGainState() noexcept
    {
        for (auto& v : currentGainDbByChannel)         v = 0.0f;
        for (auto& v : desiredGainDbByChannel)         v = 0.0f;
        for (auto& v : learnedCandidateGainDbByChannel) v = 0.0f;
        for (auto& v : heldGainDbByChannel)            v = 0.0f;

        wasPlayingLastBlock = false;
        learnWasOnLastBlock = false;
        haveCommittedHoldGain = false;

        metering.gainDbCurrent.store   (0.0f, std::memory_order_relaxed);
        metering.gainDbBlockPeak.store (0.0f, std::memory_order_relaxed);
        metering.gainDbHold.store      (0.0f, std::memory_order_relaxed);
    }

    void LevelerModule::updateRoutingMasksIfNeeded (const juce::AudioChannelSet& channelSet,
                                                    int numChannels) noexcept
    {
        if (numChannels != preparedNumChannels || channelSet != preparedChannelSet)
        {
            preparedNumChannels = numChannels;
            preparedChannelSet = channelSet;
            routingMasksDirty = true;
        }

        const int mcPolicyChoice       = readChoiceOrDefault (pMcPolicyChoice, 0);
        const int dialogDetectorChoice = readChoiceOrDefault (pDialogDetectorChoice, 0);
        const int dialogApplyChoice    = readChoiceOrDefault (pDialogApplyChoice, 0);
        const int lfeInDetector01      = (readParamOrDefault (pLfeInDetector, 0.0f) >= 0.5f ? 1 : 0);
        const int lfeInApply01         = (readParamOrDefault (pLfeInApply,    0.0f) >= 0.5f ? 1 : 0);

        if (! routingMasksDirty
            && mcPolicyChoice       == cachedMcPolicyChoice
            && dialogDetectorChoice == cachedDialogDetectorChoice
            && dialogApplyChoice    == cachedDialogApplyChoice
            && lfeInDetector01      == cachedLfeInDetector01
            && lfeInApply01         == cachedLfeInApply01)
            return;

        rebuildRoutingMasks (channelSet,
                             numChannels,
                             mcPolicyChoice,
                             dialogDetectorChoice,
                             dialogApplyChoice,
                             lfeInDetector01 != 0,
                             lfeInApply01 != 0);

        cachedMcPolicyChoice       = mcPolicyChoice;
        cachedDialogDetectorChoice = dialogDetectorChoice;
        cachedDialogApplyChoice    = dialogApplyChoice;
        cachedLfeInDetector01      = lfeInDetector01;
        cachedLfeInApply01         = lfeInApply01;
        routingMasksDirty          = false;
    }

    void LevelerModule::rebuildRoutingMasks (const juce::AudioChannelSet& channelSet,
                                             int numChannels,
                                             int mcPolicyChoice,
                                             int dialogDetectorChoice,
                                             int dialogApplyChoice,
                                             bool lfeInDetector,
                                             bool lfeInApply) noexcept
    {
        numChannels = juce::jlimit (0, maxSupportedChannels, numChannels);

        uint16_t allMask    = 0;
        uint16_t nonLfeMask = 0;
        uint16_t lfeMask    = 0;
        uint16_t cMask      = 0;
        uint16_t lcrMask    = 0;

        for (int ch = 0; ch < numChannels; ++ch)
        {
            const auto bit  = bitForChannel (ch);
            const auto type = channelSet.getTypeOfChannel (ch);

            allMask |= bit;

            if (type == juce::AudioChannelSet::LFE)
            {
                lfeMask |= bit;
                continue;
            }

            nonLfeMask |= bit;

            if (type == juce::AudioChannelSet::centre)
                cMask |= bit;

            if (type == juce::AudioChannelSet::left
                || type == juce::AudioChannelSet::right
                || type == juce::AudioChannelSet::centre)
                lcrMask |= bit;
        }

        if (nonLfeMask == 0)
            nonLfeMask = allMask;

        if (mcPolicyChoice == 1) // Dialog-mask
        {
            detectMaskBits = (dialogDetectorChoice == 0 ? cMask : lcrMask);
            applyMaskBits  = (dialogApplyChoice    == 0 ? cMask : lcrMask);

            if (detectMaskBits == 0) detectMaskBits = nonLfeMask;
            if (applyMaskBits  == 0) applyMaskBits  = nonLfeMask;
        }
        else if (mcPolicyChoice == 2) // Unlinked
        {
            detectMaskBits = nonLfeMask;
            applyMaskBits  = nonLfeMask;
        }
        else // Linked
        {
            detectMaskBits = nonLfeMask;
            applyMaskBits  = nonLfeMask;
        }

        if (lfeInDetector)
            detectMaskBits |= lfeMask;

        if (lfeInApply)
            applyMaskBits |= lfeMask;

        if (detectMaskBits == 0) detectMaskBits = allMask;
        if (applyMaskBits  == 0) applyMaskBits  = allMask;
    }

    void LevelerModule::pushDetectorFrameFromAccumulatedEnergy (int numChannels) noexcept
    {
        numChannels = juce::jlimit (0, maxSupportedChannels, numChannels);

        for (int ch = 0; ch < numChannels; ++ch)
        {
            const float newEnergy =
                (frameSamples > 0 ? (float) juce::jmax (0.0, frameEnergyAccumByChannel[(size_t) ch] / (double) frameSamples)
                                  : 0.0f);

            const float oldShort =
                (detectorFramesWritten >= shortTermWindowFrames
                    ? frameEnergyHistoryByChannel[(size_t) ch][(size_t) detectorHistoryWriteIndex]
                    : 0.0f);

            int removeMomentaryIndex = detectorHistoryWriteIndex - momentaryWindowFrames;
            if (removeMomentaryIndex < 0)
                removeMomentaryIndex += shortTermWindowFrames;

            const float oldMomentary =
                (detectorFramesWritten >= momentaryWindowFrames
                    ? frameEnergyHistoryByChannel[(size_t) ch][(size_t) removeMomentaryIndex]
                    : 0.0f);

            shortTermEnergySumsByChannel[(size_t) ch] += newEnergy - oldShort;
            momentaryEnergySumsByChannel[(size_t) ch] += newEnergy - oldMomentary;

            frameEnergyHistoryByChannel[(size_t) ch][(size_t) detectorHistoryWriteIndex] = newEnergy;
            frameEnergyAccumByChannel[(size_t) ch] = 0.0;
        }

        detectorHistoryWriteIndex = (detectorHistoryWriteIndex + 1) % shortTermWindowFrames;
        ++detectorFramesWritten;

        momentaryValid = (detectorFramesWritten >= momentaryWindowFrames);
        shortTermValid = (detectorFramesWritten >= shortTermWindowFrames);

        latestMomentaryLufs = computeMaskedMomentaryLufs (detectMaskBits, numChannels);
        latestShortTermLufs = computeMaskedShortTermLufs (detectMaskBits, numChannels);
    }

    float LevelerModule::computeMaskedMomentaryLufs (uint16_t maskBits, int numChannels) const noexcept
    {
        numChannels = juce::jlimit (0, maxSupportedChannels, numChannels);

        const int count = juce::jmin (detectorFramesWritten, momentaryWindowFrames);
        if (count <= 0 || maskBits == 0)
            return -200.0f;

        float totalEnergy = 0.0f;

        for (int ch = 0; ch < numChannels; ++ch)
            if ((maskBits & bitForChannel (ch)) != 0)
                totalEnergy += (momentaryEnergySumsByChannel[(size_t) ch] / (float) count);

        return energyToLufsSafe (totalEnergy);
    }

    float LevelerModule::computeMaskedShortTermLufs (uint16_t maskBits, int numChannels) const noexcept
    {
        numChannels = juce::jlimit (0, maxSupportedChannels, numChannels);

        const int count = juce::jmin (detectorFramesWritten, shortTermWindowFrames);
        if (count <= 0 || maskBits == 0)
            return -200.0f;

        float totalEnergy = 0.0f;

        for (int ch = 0; ch < numChannels; ++ch)
            if ((maskBits & bitForChannel (ch)) != 0)
                totalEnergy += (shortTermEnergySumsByChannel[(size_t) ch] / (float) count);

        return energyToLufsSafe (totalEnergy);
    }

    float LevelerModule::computeMaskedCurrentMeasurementLufs (uint16_t maskBits,
                                                              int numChannels,
                                                              int measChoice) const noexcept
    {
        const float momentary = computeMaskedMomentaryLufs (maskBits, numChannels);
        const float shortTerm = computeMaskedShortTermLufs (maskBits, numChannels);

        switch (measChoice)
        {
            case measShortTerm:
                if (shortTermValid)
                    return shortTerm;
                if (momentaryValid)
                    return momentary;
                return shortTerm;

            case measMomentary:
                return momentary;

            case measAuto:
            default:
                return (shortTermValid ? shortTerm : momentary);
        }
    }

    float LevelerModule::computeDesiredGainDbForMeasurement (float measuredLufs,
                                                             float targetLufs,
                                                             float maxBoostDb,
                                                             float maxCutDb) const noexcept
    {
        if (! std::isfinite (measuredLufs) || measuredLufs <= -199.0f)
            return 0.0f;

        const float unclamped = targetLufs - measuredLufs;
        return juce::jlimit (-maxCutDb, maxBoostDb, unclamped);
    }

    float LevelerModule::applyRateLimitDb (float currentDb,
                                           float targetDb,
                                           float rateUpDbPerSec,
                                           float rateDownDbPerSec,
                                           int numSamples) const noexcept
    {
        const float sr = (float) juce::jmax (1.0, sampleRate);
        const float dt = (float) juce::jmax (0, numSamples) / sr;

        const float upMax   = juce::jmax (0.0f, rateUpDbPerSec)   * dt;
        const float downMax = juce::jmax (0.0f, rateDownDbPerSec) * dt;

        const float delta = targetDb - currentDb;

        if (delta > 0.0f)
            return currentDb + juce::jmin (delta, upMax);

        if (delta < 0.0f)
            return currentDb - juce::jmin (-delta, downMax);

        return currentDb;
    }

    void LevelerModule::commitLearnedGain() noexcept
    {
        for (int ch = 0; ch < maxSupportedChannels; ++ch)
            heldGainDbByChannel[(size_t) ch] = learnedCandidateGainDbByChannel[(size_t) ch];

        haveCommittedHoldGain = true;
    }

    void LevelerModule::updateMeteringFromAppliedGain (int numChannels,
                                                       int numSamples) noexcept
    {
        numChannels = juce::jlimit (0, maxSupportedChannels, numChannels);

        float rep = 0.0f;
        bool found = false;

        for (int ch = 0; ch < numChannels; ++ch)
        {
            if ((applyMaskBits & bitForChannel (ch)) == 0)
                continue;

            const float g = currentGainDbByChannel[(size_t) ch];

            if (! found || std::abs (g) > std::abs (rep))
            {
                rep = g;
                found = true;
            }
        }

        if (! found)
            rep = 0.0f;

        metering.gainDbCurrent.store   (rep, std::memory_order_relaxed);
        metering.gainDbBlockPeak.store (rep, std::memory_order_relaxed);

        float hold = metering.gainDbHold.load (std::memory_order_relaxed);

        const float dt = (float) juce::jmax (0, numSamples) / (float) juce::jmax (1.0, sampleRate);
        const float decay = 12.0f * dt;

        if (std::abs (rep) >= std::abs (hold))
        {
            hold = rep;
        }
        else
        {
            const float mag = juce::jmax (0.0f, std::abs (hold) - decay);
            hold = std::copysign (mag, hold);
        }

        metering.gainDbHold.store (hold, std::memory_order_relaxed);
    }

    // [BEGIN LVLR-PROCESS-INTERNAL-VS-HOSTGAIN]
    void LevelerModule::process (ProcessContext& ctx) noexcept
    {
        auto& audio = ctx.audio;

        const int numChannels = juce::jlimit (0, maxSupportedChannels, audio.getNumChannels());
        const int numSamples  = juce::jmax (0, ctx.numSamples);

        if (numChannels <= 0 || numSamples <= 0)
            return;

        updateRoutingMasksIfNeeded (ctx.channelSet, numChannels);

        const bool enabled = (readParamOrDefault (pEnabled01, lvlr::Defaults::enabled01) >= 0.5f);
        const bool bypass  = bypassed.load (std::memory_order_relaxed);

        const float targetLufs =
            clampFinite (readParamOrDefault (pTargetLufs, lvlr::Defaults::targetLufs),
                         lvlr::Ranges::targetMinLufs,
                         lvlr::Ranges::targetMaxLufs,
                         lvlr::Defaults::targetLufs);

        const float maxBoostDb =
            clampFinite (readParamOrDefault (pMaxBoostDb, lvlr::Defaults::maxBoostDb),
                         lvlr::Ranges::maxBoostMinDb,
                         lvlr::Ranges::maxBoostMaxDb,
                         lvlr::Defaults::maxBoostDb);

        const float maxCutDb =
            clampFinite (readParamOrDefault (pMaxCutDb, lvlr::Defaults::maxCutDb),
                         lvlr::Ranges::maxCutMinDb,
                         lvlr::Ranges::maxCutMaxDb,
                         lvlr::Defaults::maxCutDb);

        const int measChoice =
            juce::jlimit (0, 2, readChoiceOrDefault (pMeasChoice, lvlr::Defaults::measChoice));

        const int modeChoice =
            juce::jlimit (0, 1, readChoiceOrDefault (pModeChoice, lvlr::Defaults::modeChoice));

        const int controlModeChoice =
            juce::jlimit (0, 1, readChoiceOrDefault (pControlModeChoice, lvlr::Defaults::controlModeChoice));

        const bool useHostGain = (controlModeChoice == controlHostGain);

        const bool learnOn =
            (readParamOrDefault (pLearn01, lvlr::Defaults::learn01) >= 0.5f);

        const float rateUpDbPerSec =
            clampFinite (readParamOrDefault (pRateUpDbPerSec, lvlr::Defaults::rateUpDbPerSec),
                         lvlr::Ranges::rateUpMinDbPerSec,
                         lvlr::Ranges::rateUpMaxDbPerSec,
                         lvlr::Defaults::rateUpDbPerSec);

        const float rateDownDbPerSec =
            clampFinite (readParamOrDefault (pRateDownDbPerSec, lvlr::Defaults::rateDownDbPerSec),
                         lvlr::Ranges::rateDownMinDbPerSec,
                         lvlr::Ranges::rateDownMaxDbPerSec,
                         lvlr::Defaults::rateDownDbPerSec);

        const float hostGainDb =
            clampFinite (readParamOrDefault (pHostGainDb, lvlr::Defaults::hostGainDb),
                         lvlr::Ranges::hostGainMinDb,
                         lvlr::Ranges::hostGainMaxDb,
                         lvlr::Defaults::hostGainDb);

        // Learn-Hold commit edges are relevant only for the Internal algorithm path.
        if (! useHostGain && modeChoice == modeLearnHold)
        {
            if (learnWasOnLastBlock && ! learnOn)
                commitLearnedGain();

            if (wasPlayingLastBlock && ! ctx.isPlaying && learnOn)
                commitLearnedGain();
        }

        const int mcPolicyChoice = readChoiceOrDefault (pMcPolicyChoice, 0);

        for (int ch = 0; ch < numChannels; ++ch)
            desiredGainDbByChannel[(size_t) ch] = 0.0f;

        if (enabled && ! bypass)
        {
            if (useHostGain)
            {
                for (int ch = 0; ch < numChannels; ++ch)
                {
                    const uint16_t bit = bitForChannel (ch);
                    desiredGainDbByChannel[(size_t) ch] = ((applyMaskBits & bit) != 0 ? hostGainDb : 0.0f);
                }
            }
            else if (mcPolicyChoice == 2) // Unlinked
            {
                for (int ch = 0; ch < numChannels; ++ch)
                {
                    const uint16_t bit = bitForChannel (ch);

                    if ((applyMaskBits & bit) == 0)
                    {
                        desiredGainDbByChannel[(size_t) ch] = 0.0f;
                        continue;
                    }

                    const float measured =
                        ((detectMaskBits & bit) != 0
                            ? computeMaskedCurrentMeasurementLufs (bit, numChannels, measChoice)
                            : -200.0f);

                    const float candidate =
                        computeDesiredGainDbForMeasurement (measured, targetLufs, maxBoostDb, maxCutDb);

                    if (modeChoice == modeAdaptive)
                    {
                        desiredGainDbByChannel[(size_t) ch] = candidate;
                    }
                    else
                    {
                        if (learnOn)
                            learnedCandidateGainDbByChannel[(size_t) ch] = candidate;

                        desiredGainDbByChannel[(size_t) ch] =
                            (haveCommittedHoldGain ? heldGainDbByChannel[(size_t) ch] : 0.0f);
                    }
                }
            }
            else // Linked / Dialog-mask
            {
                const float sharedMeasured =
                    computeMaskedCurrentMeasurementLufs (detectMaskBits, numChannels, measChoice);

                const float sharedCandidate =
                    computeDesiredGainDbForMeasurement (sharedMeasured, targetLufs, maxBoostDb, maxCutDb);

                for (int ch = 0; ch < numChannels; ++ch)
                {
                    const uint16_t bit = bitForChannel (ch);

                    if ((applyMaskBits & bit) == 0)
                    {
                        desiredGainDbByChannel[(size_t) ch] = 0.0f;
                        continue;
                    }

                    if (modeChoice == modeAdaptive)
                    {
                        desiredGainDbByChannel[(size_t) ch] = sharedCandidate;
                    }
                    else
                    {
                        if (learnOn)
                            learnedCandidateGainDbByChannel[(size_t) ch] = sharedCandidate;

                        desiredGainDbByChannel[(size_t) ch] =
                            (haveCommittedHoldGain ? heldGainDbByChannel[(size_t) ch] : 0.0f);
                    }
                }
            }
        }

        std::array<float, maxSupportedChannels> blockStartGainDbByChannel {};
        std::array<float, maxSupportedChannels> blockEndGainDbByChannel   {};
        std::array<float*, maxSupportedChannels> writePtrs {};
        std::array<float,  maxSupportedChannels> gainLinByChannel {};

        for (int ch = 0; ch < numChannels; ++ch)
        {
            blockStartGainDbByChannel[(size_t) ch] = currentGainDbByChannel[(size_t) ch];

            blockEndGainDbByChannel[(size_t) ch] =
                (useHostGain
                    ? desiredGainDbByChannel[(size_t) ch]
                    : applyRateLimitDb (currentGainDbByChannel[(size_t) ch],
                                        desiredGainDbByChannel[(size_t) ch],
                                        rateUpDbPerSec,
                                        rateDownDbPerSec,
                                        numSamples));

            currentGainDbByChannel[(size_t) ch] = blockEndGainDbByChannel[(size_t) ch];
            writePtrs[(size_t) ch] = audio.getWritePointer (ch);

            if (! useHostGain)
                gainLinByChannel[(size_t) ch] = juce::Decibels::decibelsToGain (blockEndGainDbByChannel[(size_t) ch]);
        }

        if (ctx.isDiscontinuity)
            resetDetectorState();

        const bool runDetector = (enabled && ! bypass && ! ctx.freezeAnalysis);

        if (useHostGain)
        {
            for (int i = 0; i < numSamples; ++i)
            {
                if (runDetector)
                {
                    for (int ch = 0; ch < numChannels; ++ch)
                    {
                        const float x = writePtrs[(size_t) ch][i];
                        const float y = kWeight.processSample (ch, x);
                        frameEnergyAccumByChannel[(size_t) ch] += (double) y * (double) y;
                    }

                    if (--samplesUntilNextFrame <= 0)
                    {
                        samplesUntilNextFrame += frameSamples;
                        pushDetectorFrameFromAccumulatedEnergy (numChannels);
                    }
                }

                const float alpha = (numSamples > 1 ? (float) i / (float) (numSamples - 1) : 1.0f);

                for (int ch = 0; ch < numChannels; ++ch)
                {
                    const float startDb = blockStartGainDbByChannel[(size_t) ch];
                    const float endDb   = blockEndGainDbByChannel[(size_t) ch];
                    const float gDb     = startDb + (endDb - startDb) * alpha;
                    writePtrs[(size_t) ch][i] *= juce::Decibels::decibelsToGain (gDb);
                }
            }
        }
        else
        {
            for (int i = 0; i < numSamples; ++i)
            {
                if (runDetector)
                {
                    for (int ch = 0; ch < numChannels; ++ch)
                    {
                        const float x = writePtrs[(size_t) ch][i];
                        const float y = kWeight.processSample (ch, x);
                        frameEnergyAccumByChannel[(size_t) ch] += (double) y * (double) y;
                    }

                    if (--samplesUntilNextFrame <= 0)
                    {
                        samplesUntilNextFrame += frameSamples;
                        pushDetectorFrameFromAccumulatedEnergy (numChannels);
                    }
                }

                for (int ch = 0; ch < numChannels; ++ch)
                    writePtrs[(size_t) ch][i] *= gainLinByChannel[(size_t) ch];
            }
        }

        updateMeteringFromAppliedGain (numChannels, numSamples);

        wasPlayingLastBlock = ctx.isPlaying;
        learnWasOnLastBlock = learnOn;
    }
    // [END LVLR-PROCESS-INTERNAL-VS-HOSTGAIN]
    // [END LVLR-MODULE-IMPL]
} // namespace levelscope
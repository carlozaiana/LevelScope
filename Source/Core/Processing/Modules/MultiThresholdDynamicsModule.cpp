#include "MultiThresholdDynamicsModule.h"
// [BEGIN MTDM-UPWARD-STRATEGY-INCLUDES]
#include <cmath>
// [END MTDM-UPWARD-STRATEGY-INCLUDES]

namespace levelscope
{
    // [BEGIN MTDM-MC-POLICY-MASK-HELPER-IMPL]
    void MultiThresholdDynamicsModule::updateChannelMasksIfNeeded (const ProcessContext& ctx,
                                                                   int policyChoice,
                                                                   int dialogDetectorChoice,
                                                                   int dialogApplyChoice,
                                                                   bool lfeInDetector,
                                                                   bool lfeInApply) noexcept
    {
        // AUDIO THREAD ONLY. No allocations, no locks.

        const int numChBuf = ctx.audio.getNumChannels();
        const int numChSet = ctx.channelSet.size();
        const int numCh = juce::jlimit (0, kMaxPolicyChannels, juce::jmin (numChBuf, numChSet));

        const bool dirty =
            (policyChoice != lastPolicyChoice) ||
            (dialogDetectorChoice != lastDialogDetectorChoice) ||
            (dialogApplyChoice    != lastDialogApplyChoice) ||
            (lfeInDetector != lastLfeInDetector) ||
            (lfeInApply    != lastLfeInApply) ||
            (numCh != lastMaskNumChannels) ||
            (ctx.channelSet != lastMaskChannelSet);

        if (! dirty)
            return;

        lastPolicyChoice = policyChoice;
        lastDialogDetectorChoice = dialogDetectorChoice;
        lastDialogApplyChoice    = dialogApplyChoice;
        lastLfeInDetector = lfeInDetector;
        lastLfeInApply    = lfeInApply;
        lastMaskNumChannels = numCh;
        lastMaskChannelSet  = ctx.channelSet;

        auto bitFor = [] (int ch) noexcept -> uint16_t
        {
            if (ch < 0 || ch >= kMaxPolicyChannels) return 0;
            return (uint16_t) (1u << (unsigned) ch);
        };

        uint16_t nonLfeBits = 0;
        uint16_t lfeBits    = 0;
        uint16_t centreBits = 0;
        uint16_t lcrBits    = 0;
        uint16_t allBits    = 0;

        for (int ch = 0; ch < numCh; ++ch)
        {
            const uint16_t b = bitFor (ch);
            allBits |= b;

            const auto t = ctx.channelSet.getTypeOfChannel (ch);

            if (t == juce::AudioChannelSet::LFE)
                lfeBits |= b;
            else
                nonLfeBits |= b;

            if (t == juce::AudioChannelSet::centre)
                centreBits |= b;

            if (t == juce::AudioChannelSet::left
                || t == juce::AudioChannelSet::right
                || t == juce::AudioChannelSet::centre)
                lcrBits |= b;
        }

        // Base masks from policy/dialog choices (LFE excluded at base layer)
        uint16_t baseDetect = 0;
        uint16_t baseApply  = 0;

        if (policyChoice == 1) // Dialog-mask
        {
            baseDetect = (dialogDetectorChoice == 0 ? centreBits : lcrBits);
            baseApply  = (dialogApplyChoice    == 0 ? centreBits : lcrBits);

            // Fallbacks if layout lacks requested channels
            if (baseDetect == 0) baseDetect = nonLfeBits;
            if (baseApply  == 0) baseApply  = nonLfeBits;
        }
        else // Linked (0) or Unlinked (2): all non-LFE by default
        {
            baseDetect = nonLfeBits;
            baseApply  = nonLfeBits;

            if (baseDetect == 0) baseDetect = allBits;
            if (baseApply  == 0) baseApply  = allBits;
        }

        // Apply LFE override layer (add/remove LFE bit)
        uint16_t det = baseDetect;
        uint16_t app = baseApply;

        if (lfeInDetector) det |= lfeBits; else det = (uint16_t) (det & (uint16_t) ~lfeBits);
        if (lfeInApply)    app |= lfeBits; else app = (uint16_t) (app & (uint16_t) ~lfeBits);

        // Final safety fallback
        if (det == 0) det = allBits;
        if (app == 0) app = allBits;

        detectMaskBits = det;
        applyMaskBits  = app;
    }
    // [END MTDM-MC-POLICY-MASK-HELPER-IMPL]
    // [BEGIN MTDM-MODULE-IMPL]
    juce::String MultiThresholdDynamicsModule::getModuleID() const
    {
        return kModuleID;
    }

    juce::String MultiThresholdDynamicsModule::getDisplayName() const
    {
        return "Multi-threshold Dynamics";
    }

    void MultiThresholdDynamicsModule::prepare (const ModulePrepareSpec& spec)
    {
        preparedSampleRate = spec.sampleRate;
        preparedMaxBlockSize = spec.maxBlockSize;
        preparedChannelSet = spec.channelSet;

        // [BEGIN MTDM-ZONE-AUDITION-PREPARE-GATE-BUFFERS]
        // Allocate once (NON-RT) for zone audition gate alignment.
        // Worst-case latency: SUC up to 8192 + limiter lookahead up to ~50ms + headroom.
        const int maxUpwardLatency = 8192;
        const int maxLimiterLatency = (int) std::ceil (spec.sampleRate * 0.050) + 4096; // generous
        const int maxTotalLatency = maxUpwardLatency + maxLimiterLatency;

        gateLineSize = juce::jmax (1, maxTotalLatency + spec.maxBlockSize + 8);
        zoneGateDelayLine.assign ((size_t) gateLineSize, 1.0f);

        zoneGateScratch.assign ((size_t) juce::jmax (1, spec.maxBlockSize), 1.0f);

        gateWritePos = 0;
        gateDelaySamples = 0;
        gateReadPos = 0;
        // [END MTDM-ZONE-AUDITION-PREPARE-GATE-BUFFERS]

        // [BEGIN MTDM-PREPARE-UPWARD-STRATEGIES]
        spectralUpwardProcessor.prepare (spec);
        broadbandUpwardProcessor.prepare (spec);

        // Default active strategy is Spectral; actual mode is selected in process() from APVTS.
        activeUpward = &spectralUpwardProcessor;
        lastUpwardModeChoice = -1;
        // [BEGIN MTDM-PREPARE-DOWNWARD]
            downwardProcessor.prepare (spec);
        // [END MTDM-PREPARE-DOWNWARD]
        // [END MTDM-PREPARE-UPWARD-STRATEGIES]

        // [BEGIN MTDM-PREPARE-LIMITER]
            limiterStage.prepare (spec);
        // [END MTDM-PREPARE-LIMITER]

        // [BEGIN MTDM-PREPARE-CLEAR-EXTERNAL-MASKS]
        // Ensure DSP blocks start with no external mask override (they use legacy LFE policy until MTDM sets masks).
            spectralUpwardProcessor.suc.setChannelMasksAudioThread (0, 0, false);
            broadbandUpwardProcessor.buc.setChannelMasksAudioThread (0, 0, false);
            downwardProcessor.comp.setChannelMasksAudioThread (0, 0, false);
        // [END MTDM-PREPARE-CLEAR-EXTERNAL-MASKS]

        // [BEGIN MTDM-UNTOUCHED-AUDITION-PREPARE]
        // Fixed audition detector times (can be parameterized later if needed)
            const float sr = (float) juce::jmax (1.0, preparedSampleRate);

            const float detAttackS  = 0.010f; // 10 ms
            const float detReleaseS = 0.100f; // 100 ms

            auditionDetA = std::exp (-1.0f / (detAttackS  * sr));
            auditionDetR = std::exp (-1.0f / (detReleaseS * sr));

            const float gateS = 0.005f; // 5 ms gate smoothing to avoid clicks
            auditionGateA = std::exp (-1.0f / (gateS * sr));

            auditionEnvMS = 0.0f;
            auditionGateZ = 1.0f;
            auditionWasActive = false;
        // [END MTDM-UNTOUCHED-AUDITION-PREPARE]
    }

        // [BEGIN MTDM-RESET-UPWARD-STRATEGIES]
        void MultiThresholdDynamicsModule::reset()
        {
            spectralUpwardProcessor.reset();
            broadbandUpwardProcessor.reset();
            // [BEGIN MTDM-RESET-DOWNWARD]
            downwardProcessor.reset();
            // [END MTDM-RESET-DOWNWARD]
            // [BEGIN MTDM-RESET-LIMITER]
            limiterStage.reset();
            // [END MTDM-RESET-LIMITER]
            // [BEGIN MTDM-UNTOUCHED-AUDITION-RESET]
            auditionEnvMS = 0.0f;
            auditionGateZ = 1.0f;
            auditionWasActive = false;
            // [END MTDM-UNTOUCHED-AUDITION-RESET]
            // [BEGIN MTDM-UPWARD-METERING-RESET]
            upwardHoldDbInternal = 0.0f;
            upwardHoldSamplesLeft = 0;

            upwardMetering.boostDbCurrent.store   (0.0f, std::memory_order_relaxed);
            upwardMetering.boostDbBlockPeak.store (0.0f, std::memory_order_relaxed);
            upwardMetering.boostDbHold.store      (0.0f, std::memory_order_relaxed);
            // [END MTDM-UPWARD-METERING-RESET]

            // [BEGIN MTDM-ZONE-AUDITION-RESET-GATE-BUFFERS]
            if (! zoneGateDelayLine.empty())
                std::fill (zoneGateDelayLine.begin(), zoneGateDelayLine.end(), 1.0f);

            if (! zoneGateScratch.empty())
                std::fill (zoneGateScratch.begin(), zoneGateScratch.end(), 1.0f);

            gateWritePos = 0;
            gateDelaySamples = 0;
            gateReadPos = 0;
            // [END MTDM-ZONE-AUDITION-RESET-GATE-BUFFERS]
        }
        // [END MTDM-RESET-UPWARD-STRATEGIES]

    // [BEGIN MTDM-UPWARD-STRATEGY-IMPL]
    void MultiThresholdDynamicsModule::SpectralUpwardProcessor::prepare (const ModulePrepareSpec& spec)
    {
        const int numCh = juce::jmax (1, spec.channelSet.size());
        suc.prepare (spec.sampleRate, numCh, spec.channelSet, spec.maxBlockSize);
        prepared = true;
    }

    void MultiThresholdDynamicsModule::SpectralUpwardProcessor::reset() noexcept
    {
        if (prepared)
            suc.reset();
    }

    void MultiThresholdDynamicsModule::SpectralUpwardProcessor::process (juce::AudioBuffer<float>& audio,
                                                                         const UpwardRuntimeParams& rp) noexcept
    {
        if (! prepared)
            return;
        // [BEGIN MTDM-SUC-PASS-STAGE-E-MASKS]
        suc.setChannelMasksAudioThread (rp.detectMaskBits, rp.applyMaskBits, rp.unlinked);
        // [END MTDM-SUC-PASS-STAGE-E-MASKS]

        levelscope::dsp::SpectralUpwardCompressor::Parameters p;
        p.t0Lufs = rp.t0Lufs;
        p.t1Lufs = rp.t1Lufs;

        p.amount01   = rp.amount01;
        p.maxBoostDb = rp.maxBoostDb;
        p.curve      = rp.curve;
        p.lowKneeDb  = rp.lowKneeDb;
        p.highKneeDb = rp.highKneeDb;
        p.attackMs   = rp.attackMs;
        p.releaseMs  = rp.releaseMs;

        p.calibrationTrimDb = rp.calibrationTrimDb;

        p.fftSizeChoice     = rp.fftSizeChoice;
        p.bandsPerOctChoice = rp.bandsPerOctChoice;
        p.minFreqHz         = rp.minFreqHz;
        p.maxFreqHz         = rp.maxFreqHz;

        p.curveType = (rp.curveTypeChoice == 1
                         ? levelscope::dsp::SpectralUpwardCompressor::CurveType::bell
                         : levelscope::dsp::SpectralUpwardCompressor::CurveType::monotonic);
        
        // [BEGIN MTDM-SUC-SET-LFE-MASK]
        p.lfeInDetector = rp.lfeInDetector;
        p.lfeInApply    = rp.lfeInApply;
        // [END MTDM-SUC-SET-LFE-MASK]

        // [BEGIN MTDM-SUC-PASS-AUDITION-BYPASS]
        p.auditionBypass = rp.auditionBypass;
        // [END MTDM-SUC-PASS-AUDITION-BYPASS]

        suc.setParametersAudioThread (p);
        suc.process (audio);
    }

    int MultiThresholdDynamicsModule::SpectralUpwardProcessor::getLatencySamples() const noexcept
    {
        return suc.getLatencySamples();
    }

    void MultiThresholdDynamicsModule::BroadbandUpwardProcessor::prepare (const ModulePrepareSpec& spec)
    {
        const int numCh = juce::jmax (1, spec.channelSet.size());
        buc.prepare (spec.sampleRate, numCh, spec.channelSet, spec.maxBlockSize);
        prepared = true;
    }

    void MultiThresholdDynamicsModule::BroadbandUpwardProcessor::reset() noexcept
    {
        if (prepared)
            buc.reset();
    }

    void MultiThresholdDynamicsModule::BroadbandUpwardProcessor::process (juce::AudioBuffer<float>& audio,
                                                                          const UpwardRuntimeParams& rp) noexcept
    {
        if (! prepared)
            return;
        // [BEGIN MTDM-BUC-PASS-STAGE-E-MASKS]
        buc.setChannelMasksAudioThread (rp.detectMaskBits, rp.applyMaskBits, rp.unlinked);
        // [END MTDM-BUC-PASS-STAGE-E-MASKS]

        levelscope::dsp::BroadbandUpwardCompressor::Parameters p;
        p.t0Lufs = rp.t0Lufs;
        p.t1Lufs = rp.t1Lufs;

        p.amount01   = rp.amount01;
        p.maxBoostDb = rp.maxBoostDb;
        p.curve      = rp.curve;
        p.lowKneeDb  = rp.lowKneeDb;
        p.highKneeDb = rp.highKneeDb;

        // minimal set: reuse attack/release
        p.attackMs   = rp.attackMs;
        p.releaseMs  = rp.releaseMs;

        p.curveType = (rp.curveTypeChoice == 1
                         ? levelscope::dsp::BroadbandUpwardCompressor::CurveType::bell
                         : levelscope::dsp::BroadbandUpwardCompressor::CurveType::monotonic);

        // [BEGIN MTDM-BUC-SET-LFE-MASK]
        p.lfeInDetector = rp.lfeInDetector;
        p.lfeInApply    = rp.lfeInApply;
        // [END MTDM-BUC-SET-LFE-MASK]                 

        buc.setParametersAudioThread (p);
        buc.process (audio);
    }

    int MultiThresholdDynamicsModule::BroadbandUpwardProcessor::getLatencySamples() const noexcept
    {
        return 0;
    }
    // [END MTDM-UPWARD-STRATEGY-IMPL]

    // [BEGIN MTDM-DOWNWARD-STRATEGY-IMPL]
        void MultiThresholdDynamicsModule::BroadbandDownwardProcessor::prepare (const ModulePrepareSpec& spec)
        {
            const int numCh = juce::jmax (1, spec.channelSet.size());
            comp.prepare (spec.sampleRate, numCh, spec.channelSet, spec.maxBlockSize);
            prepared = true;
        }

        void MultiThresholdDynamicsModule::BroadbandDownwardProcessor::reset() noexcept
        {
            if (prepared)
                comp.reset();
        }

        void MultiThresholdDynamicsModule::BroadbandDownwardProcessor::process (juce::AudioBuffer<float>& audio,
                                                                                const DownwardRuntimeParams& rp) noexcept
        {
            if (! prepared)
                return;
            // [BEGIN MTDM-BDC-PASS-STAGE-E-MASKS]
            comp.setChannelMasksAudioThread (rp.detectMaskBits, rp.applyMaskBits, rp.unlinked);
            // [END MTDM-BDC-PASS-STAGE-E-MASKS]

            levelscope::dsp::BroadbandDownwardCompressor::Parameters p;
            p.enabled   = rp.enabled;
            p.t2Lufs    = rp.t2Lufs;
            p.t3Lufs    = rp.t3Lufs;
            p.ratio     = rp.ratio;
            p.kneeDb    = rp.kneeDb;
            p.attackMs  = rp.attackMs;
            p.releaseMs = rp.releaseMs;
            p.makeupDb  = rp.makeupDb;

            comp.setParametersAudioThread (p);
            // [BEGIN MTDM-BDC-REMOVE-OLD-LFE-CALL]
            // LFE policy is now included in detectMaskBits/applyMaskBits; no separate call needed.
            // comp.setLfePolicyAudioThread (rp.lfeInDetector, rp.lfeInApply);
            // [END MTDM-BDC-REMOVE-OLD-LFE-CALL]
            comp.process (audio);
        }
    // [END MTDM-DOWNWARD-STRATEGY-IMPL]

    // [BEGIN MTDM-LIMITER-STRATEGY-IMPL]
    void MultiThresholdDynamicsModule::LookaheadLimiterStage::prepare (const ModulePrepareSpec& spec)
    {
        const int numCh = juce::jmax (1, spec.channelSet.size());
        lim.prepare (spec.sampleRate, numCh, spec.channelSet, spec.maxBlockSize);
        prepared = true;
    }

    void MultiThresholdDynamicsModule::LookaheadLimiterStage::reset() noexcept
    {
        if (prepared)
            lim.reset();
    }

    void MultiThresholdDynamicsModule::LookaheadLimiterStage::process (juce::AudioBuffer<float>& audio,
                                                                       const LimiterRuntimeParams& rp) noexcept
    {
        if (! prepared)
            return;

        // [BEGIN MTDM-LIMITER-WRAPPER-TP]
        levelscope::dsp::LookaheadLimiter::Parameters p;
        p.enabled      = rp.enabled;
        p.ceilingDb    = rp.ceilingDb;
        p.lookaheadMs  = rp.lookaheadMs;
        p.attackMs     = rp.attackMs;
        p.releaseMs    = rp.releaseMs;
        p.driveDb      = rp.driveDb;
        p.oversamplingChoice = rp.oversamplingChoice;
        // [END MTDM-LIMITER-WRAPPER-TP]

        // [BEGIN MTDM-LIMITER-PASS-AUDITION-BYPASS]
        p.auditionBypass = rp.auditionBypass;
        // [END MTDM-LIMITER-PASS-AUDITION-BYPASS]

        lim.setParametersAudioThread (p);
        lim.process (audio);
    }

    int MultiThresholdDynamicsModule::LookaheadLimiterStage::getLatencySamples() const noexcept
    {
        return lim.getLatencySamples();
    }
    // [END MTDM-LIMITER-STRATEGY-IMPL]

    // [BEGIN MTDM-BINDPARAMS-FULL-IMPL]
    void MultiThresholdDynamicsModule::bindParameters (std::atomic<float>* enabled01,
                                                      std::atomic<float>* thresholdDb,
                                                      std::atomic<float>* ratio,

                                                      std::atomic<float>* t0Lufs,
                                                      std::atomic<float>* t1Lufs,

                                                      std::atomic<float>* sucAmount01,
                                                      std::atomic<float>* sucMaxBoostDb,
                                                      std::atomic<float>* sucCurve,
                                                      std::atomic<float>* sucLowKneeDb,
                                                      std::atomic<float>* sucHighKneeDb,
                                                      std::atomic<float>* sucAttackMs,
                                                      std::atomic<float>* sucReleaseMs,

                                                      std::atomic<float>* sucFftSizeChoice,
                                                      std::atomic<float>* sucBandsPerOctChoice,
                                                      std::atomic<float>* sucMinFreqHz,
                                                      std::atomic<float>* sucMaxFreqHz,

                                                      std::atomic<float>* sucCalTrimDb,
                                                      std::atomic<float>* sucCurveTypeChoice,

                                                      std::atomic<float>* upwardModeChoice,
                                                      std::atomic<float>* lfeInDetector01,
                                                      std::atomic<float>* lfeInApply01,

                                                      std::atomic<float>* t2Lufs,
                                                      std::atomic<float>* t3Lufs,
                                                      std::atomic<float>* downEnabled01,
                                                      std::atomic<float>* downRatio,
                                                      std::atomic<float>* downKneeDb,
                                                      std::atomic<float>* downAttackMs,
                                                      std::atomic<float>* downReleaseMs,
                                                      std::atomic<float>* downMakeupDb,

                                                      std::atomic<float>* limEnabled01,
                                                      std::atomic<float>* limCeilingDb,
                                                      std::atomic<float>* limLookaheadMs,
                                                      std::atomic<float>* limReleaseMs,
                                                      std::atomic<float>* limAttackMs,
                                                      std::atomic<float>* limDriveDb,
                                                      std::atomic<float>* limOversamplingChoice,

                                                      std::atomic<float>* zoneSoloChoice,
                                                      std::atomic<float>* zoneUpwardMute01,
                                                      std::atomic<float>* zoneDownwardMute01,
                                                      std::atomic<float>* zoneLimiterMute01,
                                                      std::atomic<float>* zoneUntouchedMute01,
                                                      // [BEGIN MTDM-MC-POLICY-BINDPARAMS-IMPL]
                                                      std::atomic<float>* mcPolicyChoice,
                                                      std::atomic<float>* dialogDetectorChoice,
                                                      std::atomic<float>* dialogApplyChoice,
                                                      // [BEGIN MTDM-ZONE-AUDITION-BINDPARAMS-IMPL]
                                                      std::atomic<float>* zoneAudBelowT0_01,
                                                      std::atomic<float>* zoneAudT0T1_01,
                                                      std::atomic<float>* zoneAudT1T2_01,
                                                      std::atomic<float>* zoneAudT2T3_01,
                                                      std::atomic<float>* zoneAudAboveT3_01,

                                                      // [BEGIN MTDM-STAGE-BYPASS-BINDPARAMS-IMPL]
                                                      std::atomic<float>* upEnabled01,
                                                      std::atomic<float>* upBypass01,
                                                      std::atomic<float>* downBypass01,
                                                      std::atomic<float>* limBypass01) noexcept
                                                      // [END MTDM-STAGE-BYPASS-BINDPARAMS-IMPL]
                                                      // [END MTDM-ZONE-AUDITION-BINDPARAMS-IMPL]
                                                      // [END MTDM-MC-POLICY-BINDPARAMS-IMPL]
    {
        pEnabled01   = enabled01;
        pThresholdDb = thresholdDb;
        pRatio       = ratio;

        pT0Lufs = t0Lufs;
        pT1Lufs = t1Lufs;

        pSucAmount01   = sucAmount01;
        pSucMaxBoostDb = sucMaxBoostDb;
        pSucCurve      = sucCurve;
        pSucLowKneeDb  = sucLowKneeDb;
        pSucHighKneeDb = sucHighKneeDb;
        pSucAttackMs   = sucAttackMs;
        pSucReleaseMs  = sucReleaseMs;

        pSucFftSizeChoice     = sucFftSizeChoice;
        pSucBandsPerOctChoice = sucBandsPerOctChoice;
        pSucMinFreqHz         = sucMinFreqHz;
        pSucMaxFreqHz         = sucMaxFreqHz;

        pSucCalTrimDb       = sucCalTrimDb;
        pSucCurveTypeChoice = sucCurveTypeChoice;

        pUpwardModeChoice = upwardModeChoice;
        pLfeInDetector01  = lfeInDetector01;
        pLfeInApply01     = lfeInApply01;

        pT2Lufs = t2Lufs;
        pT3Lufs = t3Lufs;

        pDownEnabled01 = downEnabled01;
        pDownRatio     = downRatio;
        pDownKneeDb    = downKneeDb;
        pDownAttackMs  = downAttackMs;
        pDownReleaseMs = downReleaseMs;
        pDownMakeupDb  = downMakeupDb;

        pLimEnabled01   = limEnabled01;
        pLimCeilingDb   = limCeilingDb;
        pLimLookaheadMs = limLookaheadMs;
        pLimReleaseMs   = limReleaseMs;
        pLimAttackMs    = limAttackMs;
        pLimDriveDb     = limDriveDb;
        pLimOversamplingChoice = limOversamplingChoice;

        pZoneSoloChoice        = zoneSoloChoice;
        pZoneUpwardMute01      = zoneUpwardMute01;
        pZoneDownwardMute01    = zoneDownwardMute01;
        pZoneLimiterMute01     = zoneLimiterMute01;
        pZoneUntouchedMute01 = zoneUntouchedMute01;
        // [BEGIN MTDM-MC-POLICY-BINDPARAMS-ASSIGN]
        pMcPolicyChoice       = mcPolicyChoice;
        pDialogDetectorChoice = dialogDetectorChoice;
        pDialogApplyChoice    = dialogApplyChoice;
        // [END MTDM-MC-POLICY-BINDPARAMS-ASSIGN]

        // [BEGIN MTDM-ZONE-AUDITION-BINDPARAMS-ASSIGN]
        pZoneAudBelowT0_01 = zoneAudBelowT0_01;
        pZoneAudT0T1_01    = zoneAudT0T1_01;
        pZoneAudT1T2_01    = zoneAudT1T2_01;
        pZoneAudT2T3_01    = zoneAudT2T3_01;
        pZoneAudAboveT3_01 = zoneAudAboveT3_01;
        // [END MTDM-ZONE-AUDITION-BINDPARAMS-ASSIGN]

        // [BEGIN MTDM-STAGE-BYPASS-BINDPARAMS-ASSIGN]
        pUpEnabled01  = upEnabled01;
        pUpBypass01   = upBypass01;

        pDownBypass01 = downBypass01;

        pLimBypass01  = limBypass01;
        // [END MTDM-STAGE-BYPASS-BINDPARAMS-ASSIGN]
    }
    // [END MTDM-BINDPARAMS-FULL-IMPL]

    // [BEGIN MTDM-PROCESS-STAGE-D1A]
        void MultiThresholdDynamicsModule::process (ProcessContext& ctx) noexcept
        {
            // Default is disabled => pass-through (no audible change, no latency insertion).
            const float enabled = (pEnabled01 != nullptr
                                     ? pEnabled01->load (std::memory_order_relaxed)
                                     : 0.0f);

            if (enabled < 0.5f)
                return;

            // [BEGIN MTDM-MC-POLICY-READ-AND-MASK-UPDATE]
            const int policyChoice = (int) std::lround (pMcPolicyChoice != nullptr
                                                          ? pMcPolicyChoice->load (std::memory_order_relaxed)
                                                          : (float) levelscope::mtdm::Defaults::mcPolicyChoice);

            const int dialogDetChoice = (int) std::lround (pDialogDetectorChoice != nullptr
                                                             ? pDialogDetectorChoice->load (std::memory_order_relaxed)
                                                             : (float) levelscope::mtdm::Defaults::dialogDetectorChoice);

            const int dialogApplyChoice = (int) std::lround (pDialogApplyChoice != nullptr
                                                               ? pDialogApplyChoice->load (std::memory_order_relaxed)
                                                               : (float) levelscope::mtdm::Defaults::dialogApplyChoice);

            const bool lfeDet = (pLfeInDetector01 != nullptr
                                   ? (pLfeInDetector01->load (std::memory_order_relaxed) >= 0.5f)
                                   : (levelscope::mtdm::Defaults::lfeInDetector01 >= 0.5f));

            const bool lfeApp = (pLfeInApply01 != nullptr
                                   ? (pLfeInApply01->load (std::memory_order_relaxed) >= 0.5f)
                                   : (levelscope::mtdm::Defaults::lfeInApply01 >= 0.5f));

            updateChannelMasksIfNeeded (ctx, policyChoice, dialogDetChoice, dialogApplyChoice, lfeDet, lfeApp);
            // [END MTDM-MC-POLICY-READ-AND-MASK-UPDATE]

            // [BEGIN MTDM-PROCESS-UPWARD-MODE-SWITCH]
            // Select upward mode (0=Spectral, 1=Broadband)
            const int upwardMode = (int) std::lround (pUpwardModeChoice != nullptr
                                                      ? pUpwardModeChoice->load (std::memory_order_relaxed)
                                                      : (float) levelscope::mtdm::Defaults::upwardModeChoice);

            if (activeUpward == nullptr)
                activeUpward = &spectralUpwardProcessor;

            if (upwardMode != lastUpwardModeChoice)
            {
                lastUpwardModeChoice = upwardMode;

                activeUpward = (upwardMode == 1 ? static_cast<IUpwardProcessor*> (&broadbandUpwardProcessor)
                                                : static_cast<IUpwardProcessor*> (&spectralUpwardProcessor));

                // Rare user action; reset to avoid stale delay lines / FIFOs during mode switch.
                spectralUpwardProcessor.reset();
                broadbandUpwardProcessor.reset();
            }

            UpwardRuntimeParams up;

            const float t0 = (pT0Lufs != nullptr ? pT0Lufs->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::t0Lufs);
            const float t1 = (pT1Lufs != nullptr ? pT1Lufs->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::t1Lufs);
            up.t0Lufs = juce::jmin (t0, t1);
            up.t1Lufs = juce::jmax (t0, t1);

            up.amount01   = (pSucAmount01   != nullptr ? pSucAmount01->load   (std::memory_order_relaxed) : levelscope::mtdm::Defaults::sucAmount01);
            up.maxBoostDb = (pSucMaxBoostDb != nullptr ? pSucMaxBoostDb->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::sucMaxBoostDb);
            up.curve      = (pSucCurve      != nullptr ? pSucCurve->load      (std::memory_order_relaxed) : levelscope::mtdm::Defaults::sucCurve);

            up.lowKneeDb  = (pSucLowKneeDb  != nullptr ? pSucLowKneeDb->load  (std::memory_order_relaxed) : levelscope::mtdm::Defaults::sucLowKneeDb);
            up.highKneeDb = (pSucHighKneeDb != nullptr ? pSucHighKneeDb->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::sucHighKneeDb);

            up.attackMs   = (pSucAttackMs   != nullptr ? pSucAttackMs->load   (std::memory_order_relaxed) : levelscope::mtdm::Defaults::sucAttackMs);
            up.releaseMs  = (pSucReleaseMs  != nullptr ? pSucReleaseMs->load  (std::memory_order_relaxed) : levelscope::mtdm::Defaults::sucReleaseMs);

            up.calibrationTrimDb = (pSucCalTrimDb != nullptr ? pSucCalTrimDb->load (std::memory_order_relaxed)
                                                             : levelscope::mtdm::Defaults::sucCalTrimDb);

            // curve type shared across strategies
            up.curveTypeChoice = (int) std::lround (pSucCurveTypeChoice != nullptr
                                                    ? pSucCurveTypeChoice->load (std::memory_order_relaxed)
                                                    : (float) levelscope::mtdm::Defaults::sucCurveTypeChoice);

            // Spectral-only params (ignored by Broadband)
            up.fftSizeChoice = juce::jlimit (0, 3, (int) std::lround (pSucFftSizeChoice != nullptr
                                                                      ? pSucFftSizeChoice->load (std::memory_order_relaxed)
                                                                      : (float) levelscope::mtdm::Defaults::sucFftSizeChoice));

            up.bandsPerOctChoice = juce::jlimit (0, 3, (int) std::lround (pSucBandsPerOctChoice != nullptr
                                                                          ? pSucBandsPerOctChoice->load (std::memory_order_relaxed)
                                                                          : (float) levelscope::mtdm::Defaults::sucBandsPerOctChoice));

            up.minFreqHz = (pSucMinFreqHz != nullptr ? pSucMinFreqHz->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::sucMinFreqHz);
            up.maxFreqHz = (pSucMaxFreqHz != nullptr ? pSucMaxFreqHz->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::sucMaxFreqHz);

            up.minFreqHz = juce::jlimit (10.0f, 2000.0f, up.minFreqHz);
            up.maxFreqHz = juce::jlimit (1000.0f, 24000.0f, up.maxFreqHz);
            if (up.maxFreqHz < up.minFreqHz + 10.0f)
                up.maxFreqHz = up.minFreqHz + 10.0f;

            // [BEGIN MTDM-UPWARD-RP-SET-LFE-MASK]
            up.lfeInDetector = (pLfeInDetector01 != nullptr
                                  ? (pLfeInDetector01->load (std::memory_order_relaxed) >= 0.5f)
                                  : (levelscope::mtdm::Defaults::lfeInDetector01 >= 0.5f));

            up.lfeInApply = (pLfeInApply01 != nullptr
                               ? (pLfeInApply01->load (std::memory_order_relaxed) >= 0.5f)
                               : (levelscope::mtdm::Defaults::lfeInApply01 >= 0.5f));
            // [END MTDM-UPWARD-RP-SET-LFE-MASK]

            // [BEGIN MTDM-UPWARD-RP-SET-MC-MASKBITS]
            up.detectMaskBits = detectMaskBits;
            up.applyMaskBits  = applyMaskBits;
            up.unlinked       = (policyChoice == 2);
            // [END MTDM-UPWARD-RP-SET-MC-MASKBITS]

            // [BEGIN MTDM-ZONE-AUDITION-COMPUTE-GATE-PRE]
            // Zone audition toggles (A/B behavior):
            // - if none enabled => audition OFF => full audio
            // - if any enabled  => audition ON  => only selected zones pass
            const bool z0 = (pZoneAudBelowT0_01 != nullptr ? (pZoneAudBelowT0_01->load (std::memory_order_relaxed) >= 0.5f)
                                                           : (levelscope::mtdm::Defaults::zoneAudBelowT0_01 >= 0.5f));
            const bool z1 = (pZoneAudT0T1_01    != nullptr ? (pZoneAudT0T1_01->load    (std::memory_order_relaxed) >= 0.5f)
                                                           : (levelscope::mtdm::Defaults::zoneAudT0T1_01 >= 0.5f));
            const bool z2 = (pZoneAudT1T2_01    != nullptr ? (pZoneAudT1T2_01->load    (std::memory_order_relaxed) >= 0.5f)
                                                           : (levelscope::mtdm::Defaults::zoneAudT1T2_01 >= 0.5f));
            const bool z3 = (pZoneAudT2T3_01    != nullptr ? (pZoneAudT2T3_01->load    (std::memory_order_relaxed) >= 0.5f)
                                                           : (levelscope::mtdm::Defaults::zoneAudT2T3_01 >= 0.5f));
            const bool z4 = (pZoneAudAboveT3_01 != nullptr ? (pZoneAudAboveT3_01->load (std::memory_order_relaxed) >= 0.5f)
                                                           : (levelscope::mtdm::Defaults::zoneAudAboveT3_01 >= 0.5f));

            const bool zoneAuditionActive = (z0 || z1 || z2 || z3 || z4);

            // Read T2/T3 early so the audition gate uses the same semantics as processing.
            float t2Aud = (pT2Lufs != nullptr ? pT2Lufs->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::t2Lufs);
            float t3Aud = (pT3Lufs != nullptr ? pT3Lufs->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::t3Lufs);
            if (t3Aud < t2Aud) std::swap (t2Aud, t3Aud);

            // Enforce untouched semantics: T2 must be >= T1
            t2Aud = juce::jmax (t2Aud, up.t1Lufs);
            t3Aud = juce::jmax (t3Aud, t2Aud);

            const float t0Aud = up.t0Lufs;
            const float t1Aud = up.t1Lufs;

            // Determine current output latency to align gate (Spectral Upward + Limiter).
            const int upwardLatency = (lastUpwardModeChoice == 0 ? spectralUpwardProcessor.getLatencySamples() : 0);

            const bool limEnabledNow = (pLimEnabled01 != nullptr
                                         ? (pLimEnabled01->load (std::memory_order_relaxed) >= 0.5f)
                                         : (levelscope::mtdm::Defaults::limEnabled01 >= 0.5f));

            const int limiterLatency = (limEnabledNow ? limiterStage.getLatencySamples() : 0);

            const int desiredGateDelay = juce::jmax (0, upwardLatency + limiterLatency);

            if (desiredGateDelay != gateDelaySamples && gateLineSize > 0)
            {
                gateDelaySamples = desiredGateDelay;

                // Re-anchor read pointer relative to current write pointer.
                gateReadPos = gateWritePos - gateDelaySamples;
                gateReadPos %= gateLineSize;
                if (gateReadPos < 0) gateReadPos += gateLineSize;
            }

            // Build per-sample gate values from PRE-processing audio, but read delayed gate for POST-processing output.
            // Use detectMaskBits for detector channel selection.
            float* gate = (! zoneGateScratch.empty() ? zoneGateScratch.data() : nullptr);

            const float* const* rd = ctx.audio.getArrayOfReadPointers();
            const int numChBuf = ctx.audio.getNumChannels();
            const int numChSet = ctx.channelSet.size();
            const int numChDet = juce::jlimit (0, 16, juce::jmin (numChBuf, numChSet));

            uint16_t detBits = (uint16_t) (detectMaskBits & (uint16_t) ((numChDet >= 16) ? 0xFFFFu : ((1u << (unsigned) numChDet) - 1u)));

            for (int i = 0; i < ctx.numSamples; ++i)
            {
                // Linked detector energy across detector channels (mask-aware)
                double sumSq = 0.0;
                int count = 0;

                for (int ch = 0; ch < numChDet; ++ch)
                {
                    const uint16_t b = (uint16_t) (1u << (unsigned) ch);
                    if ((detBits & b) == 0)
                        continue;

                    const float x = rd[ch][i];
                    sumSq += (double) x * (double) x;
                    ++count;
                }

                if (count <= 0)
                {
                    // Fallback: all channels
                    for (int ch = 0; ch < numChDet; ++ch)
                    {
                        const float x = rd[ch][i];
                        sumSq += (double) x * (double) x;
                    }
                    count = juce::jmax (1, numChDet);
                }

                const float e = (float) (sumSq / (double) juce::jmax (1, count));

                const float a = (e > auditionEnvMS ? auditionDetA : auditionDetR);
                auditionEnvMS = a * auditionEnvMS + (1.0f - a) * e;

                const float L = (float) (-0.691 + 10.0 * std::log10 ((double) auditionEnvMS + 1.0e-12));

                int zone = 0;
                if      (L <  t0Aud) zone = 0;
                else if (L <  t1Aud) zone = 1;
                else if (L <  t2Aud) zone = 2;
                else if (L <  t3Aud) zone = 3;
                else                 zone = 4;

                const bool zoneSelected =
                    (zone == 0 ? z0 :
                     zone == 1 ? z1 :
                     zone == 2 ? z2 :
                     zone == 3 ? z3 : z4);

                const float targetGate = (zoneAuditionActive ? (zoneSelected ? 1.0f : 0.0f) : 1.0f);

                auditionGateZ = auditionGateA * auditionGateZ + (1.0f - auditionGateA) * targetGate;

                // Write current gate into delay line
                if (gateLineSize > 0 && ! zoneGateDelayLine.empty())
                {
                    zoneGateDelayLine[(size_t) gateWritePos] = auditionGateZ;
                    if (++gateWritePos >= gateLineSize) gateWritePos = 0;

                    // Read delayed gate for application later
                    float delayedGate = zoneGateDelayLine[(size_t) gateReadPos];
                    if (++gateReadPos >= gateLineSize) gateReadPos = 0;

                    if (gate != nullptr)
                        gate[i] = delayedGate;
                }
                else
                {
                    if (gate != nullptr)
                        gate[i] = auditionGateZ;
                }
            }
            // [END MTDM-ZONE-AUDITION-COMPUTE-GATE-PRE]

            // [BEGIN MTDM-ZONE-SOLO-MUTE-LOGIC]
            // Stage audition controls are deprecated (Option B). For now:
            // - user controls stages by Enabled/Amount params
            // - we always run stages normally (subject to each stage's own enabled param)
            const bool runUp   = true;
            const bool runDown = true;
            const bool runLim  = true;

            // Keep old params unused but avoid warnings in some builds.
            juce::ignoreUnused (pZoneSoloChoice, pZoneUpwardMute01, pZoneDownwardMute01, pZoneLimiterMute01, pZoneUntouchedMute01);
            // [END MTDM-ZONE-SOLO-MUTE-LOGIC]

            // [BEGIN MTDM-RUN-UPWARD]
            const bool upEnabled = (pUpEnabled01 != nullptr
                                      ? (pUpEnabled01->load (std::memory_order_relaxed) >= 0.5f)
                                      : (levelscope::mtdm::Defaults::upEnabled01 >= 0.5f));

            const bool upBypass = (pUpBypass01 != nullptr
                                     ? (pUpBypass01->load (std::memory_order_relaxed) >= 0.5f)
                                     : (levelscope::mtdm::Defaults::upBypass01 >= 0.5f));

            if (! upEnabled)
            {
                // Structural disable: stage not in chain. (Latency reporting is handled in PluginProcessor.)
                // Do nothing here.
            }
            else
            {
                // Safe bypass:
                // - Spectral: must still run to preserve STFT latency, but force unity via auditionBypass.
                // - Broadband: no latency, so we can simply skip when bypassed.
                if (lastUpwardModeChoice == 0) // Spectral
                {
                    up.auditionBypass = upBypass;
                    activeUpward->process (ctx.audio, up);
                }
                else
                {
                    if (! upBypass)
                        activeUpward->process (ctx.audio, up);
                }
            }
            // [END MTDM-RUN-UPWARD]

            // [BEGIN MTDM-UPWARD-METERING-UPDATE]
            // Upward metering: report max boost across channels (as dB >= 0).
            auto gainToBoostDb = [] (float g) noexcept
            {
                g = juce::jmax (1.0f, g);
                const float db = (float) (20.0 * std::log10 ((double) g));
                return (db > 0.0f ? db : 0.0f);
            };

            bool upwardActiveForMeter = runUp;

            float gMax  = 1.0f;
            float gLast = 1.0f;

            if (lastUpwardModeChoice == 0) // Spectral (called always, but may be auditionBypass)
            {
                if (up.auditionBypass)
                {
                    upwardActiveForMeter = false;
                }
                else
                {
                    gMax  = spectralUpwardProcessor.suc.getLastBlockMaxLinearGain();
                    gLast = spectralUpwardProcessor.suc.getLastBlockLastLinearGain();
                }
            }
            else // Broadband
            {
                if (runUp)
                {
                    gMax  = broadbandUpwardProcessor.buc.getLastBlockMaxLinearGain();
                    gLast = broadbandUpwardProcessor.buc.getLastBlockLastLinearGain();
                }
                else
                {
                    upwardActiveForMeter = false;
                }
            }

            if (! upwardActiveForMeter)
            {
                upwardHoldDbInternal = 0.0f;
                upwardHoldSamplesLeft = 0;

                upwardMetering.boostDbCurrent.store   (0.0f, std::memory_order_relaxed);
                upwardMetering.boostDbBlockPeak.store (0.0f, std::memory_order_relaxed);
                upwardMetering.boostDbHold.store      (0.0f, std::memory_order_relaxed);
            }
            else
            {
                const float boostBlockPeak = gainToBoostDb (gMax);
                const float boostCurrent   = gainToBoostDb (gLast);

                upwardMetering.boostDbBlockPeak.store (boostBlockPeak, std::memory_order_relaxed);
                upwardMetering.boostDbCurrent.store   (boostCurrent,   std::memory_order_relaxed);

                const int holdSamples = (int) std::lround (preparedSampleRate * (double) upwardHoldTimeSeconds);
                const float decayDbThisBlock =
                    (preparedSampleRate > 1.0 ? upwardHoldDecayDbPerSecond * (float) ((double) ctx.numSamples / preparedSampleRate)
                                              : 0.0f);

                if (boostBlockPeak > upwardHoldDbInternal)
                {
                    upwardHoldDbInternal = boostBlockPeak;
                    upwardHoldSamplesLeft = holdSamples;
                }
                else
                {
                    if (upwardHoldSamplesLeft > 0)
                    {
                        upwardHoldSamplesLeft -= ctx.numSamples;
                        if (upwardHoldSamplesLeft < 0)
                            upwardHoldSamplesLeft = 0;
                    }
                    else
                    {
                        upwardHoldDbInternal = std::max (0.0f, upwardHoldDbInternal - decayDbThisBlock);
                    }
                }

                upwardMetering.boostDbHold.store (upwardHoldDbInternal, std::memory_order_relaxed);
            }
            // [END MTDM-UPWARD-METERING-UPDATE]

            // [BEGIN MTDM-PROCESS-DOWNWARD]
            DownwardRuntimeParams down;

            down.enabled = (pDownEnabled01 != nullptr
                              ? (pDownEnabled01->load (std::memory_order_relaxed) >= 0.5f)
                              : (levelscope::mtdm::Defaults::downEnabled01 >= 0.5f));

            down.t2Lufs = (pT2Lufs != nullptr ? pT2Lufs->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::t2Lufs);
            down.t3Lufs = (pT3Lufs != nullptr ? pT3Lufs->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::t3Lufs);

            // Safety ordering (user may set T2/T3 in any order)
            if (down.t3Lufs < down.t2Lufs)
                std::swap (down.t2Lufs, down.t3Lufs);

            // [BEGIN MTDM-T1T2-UNTOUCHED-CLAMP]
            // Enforce T1–T2 untouched semantics: do not allow downward engagement below T1.
            // (UI will constrain too; this is a safety clamp.)
            down.t2Lufs = juce::jmax (down.t2Lufs, up.t1Lufs);
            down.t3Lufs = juce::jmax (down.t3Lufs, down.t2Lufs);
            // [END MTDM-T1T2-UNTOUCHED-CLAMP]

            down.ratio = (pDownRatio != nullptr ? pDownRatio->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::downRatio);
            down.kneeDb = (pDownKneeDb != nullptr ? pDownKneeDb->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::downKneeDb);
            down.attackMs = (pDownAttackMs != nullptr ? pDownAttackMs->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::downAttackMs);
            down.releaseMs = (pDownReleaseMs != nullptr ? pDownReleaseMs->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::downReleaseMs);
            down.makeupDb = (pDownMakeupDb != nullptr ? pDownMakeupDb->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::downMakeupDb);

            // Reuse the existing LFE policy params (already in MTDM)
            down.lfeInDetector = (pLfeInDetector01 != nullptr
                                    ? (pLfeInDetector01->load (std::memory_order_relaxed) >= 0.5f)
                                    : (levelscope::mtdm::Defaults::lfeInDetector01 >= 0.5f));

            down.lfeInApply = (pLfeInApply01 != nullptr
                                 ? (pLfeInApply01->load (std::memory_order_relaxed) >= 0.5f)
                                 : (levelscope::mtdm::Defaults::lfeInApply01 >= 0.5f));

            // [BEGIN MTDM-DOWNWARD-RP-SET-MC-MASKBITS]
            down.detectMaskBits = detectMaskBits;
            down.applyMaskBits  = applyMaskBits;
            down.unlinked       = (policyChoice == 2);
            // [END MTDM-DOWNWARD-RP-SET-MC-MASKBITS]

            // [BEGIN MTDM-RUN-DOWNWARD]
            const bool downBypass = (pDownBypass01 != nullptr
                                       ? (pDownBypass01->load (std::memory_order_relaxed) >= 0.5f)
                                       : (levelscope::mtdm::Defaults::downBypass01 >= 0.5f));

            if (! downBypass)
                downwardProcessor.process (ctx.audio, down);
            // [END MTDM-RUN-DOWNWARD]
            // [END MTDM-PROCESS-DOWNWARD]

            // [BEGIN MTDM-DOWNWARD-METERING-UPDATE]
            // Downward compressor metering: use compressor gain excluding makeup.
            // [BEGIN MTDM-DOWNWARD-METERING-GAINTOGR-RENAME]
            auto gainToGrDbDown = [] (float g) noexcept
            // [END MTDM-DOWNWARD-METERING-GAINTOGR-RENAME]
            {
                g = juce::jlimit (1.0e-6f, 1.0f, g);
                return (float) (-20.0 * std::log10 ((double) g));
            };

            // If downward disabled, publish zeros and reset hold
            if (! down.enabled)
            {
                downwardHoldDbInternal = 0.0f;
                downwardHoldSamplesLeft = 0;

                downwardMetering.grDbCurrent.store   (0.0f, std::memory_order_relaxed);
                downwardMetering.grDbBlockPeak.store (0.0f, std::memory_order_relaxed);
                downwardMetering.grDbHold.store      (0.0f, std::memory_order_relaxed);
            }
            else
            {
                const float minGain  = downwardProcessor.comp.getLastBlockMinCompGain();
                const float lastGain = downwardProcessor.comp.getLastBlockLastCompGain();

                const float grBlockPeak = gainToGrDbDown (minGain);
                const float grCurrent   = gainToGrDbDown (lastGain);

                downwardMetering.grDbBlockPeak.store (grBlockPeak, std::memory_order_relaxed);
                downwardMetering.grDbCurrent.store   (grCurrent,   std::memory_order_relaxed);

                const int holdSamples = (int) std::lround (preparedSampleRate * (double) downwardHoldTimeSeconds);
                const float decayDbThisBlock =
                    (preparedSampleRate > 1.0 ? downwardHoldDecayDbPerSecond * (float) ((double) ctx.numSamples / preparedSampleRate)
                                              : 0.0f);

                if (grBlockPeak > downwardHoldDbInternal)
                {
                    downwardHoldDbInternal = grBlockPeak;
                    downwardHoldSamplesLeft = holdSamples;
                }
                else
                {
                    if (downwardHoldSamplesLeft > 0)
                    {
                        downwardHoldSamplesLeft -= ctx.numSamples;
                        if (downwardHoldSamplesLeft < 0)
                            downwardHoldSamplesLeft = 0;
                    }
                    else
                    {
                        downwardHoldDbInternal = std::max (0.0f, downwardHoldDbInternal - decayDbThisBlock);
                    }
                }

                downwardMetering.grDbHold.store (downwardHoldDbInternal, std::memory_order_relaxed);
            }
            // [END MTDM-DOWNWARD-METERING-UPDATE]

            // [BEGIN MTDM-PROCESS-LIMITER]
            LimiterRuntimeParams lim;

            lim.enabled = (pLimEnabled01 != nullptr
                             ? (pLimEnabled01->load (std::memory_order_relaxed) >= 0.5f)
                             : (levelscope::mtdm::Defaults::limEnabled01 >= 0.5f));

            lim.ceilingDb = (pLimCeilingDb != nullptr ? pLimCeilingDb->load (std::memory_order_relaxed)
                                                      : levelscope::mtdm::Defaults::limCeilingDb);

            lim.lookaheadMs = (pLimLookaheadMs != nullptr ? pLimLookaheadMs->load (std::memory_order_relaxed)
                                                          : levelscope::mtdm::Defaults::limLookaheadMs);

            lim.releaseMs = (pLimReleaseMs != nullptr ? pLimReleaseMs->load (std::memory_order_relaxed)
                                                      : levelscope::mtdm::Defaults::limReleaseMs);

            // [BEGIN MTDM-PROCESS-LIMITER-TP]
            lim.attackMs = (pLimAttackMs != nullptr ? pLimAttackMs->load (std::memory_order_relaxed)
                                                    : levelscope::mtdm::Defaults::limAttackMs);

            lim.driveDb = (pLimDriveDb != nullptr ? pLimDriveDb->load (std::memory_order_relaxed)
                                                  : levelscope::mtdm::Defaults::limDriveDb);

            lim.oversamplingChoice = (int) std::lround (pLimOversamplingChoice != nullptr
                                                        ? pLimOversamplingChoice->load (std::memory_order_relaxed)
                                                        : (float) levelscope::mtdm::Defaults::limOversamplingChoice);
            // [END MTDM-PROCESS-LIMITER-TP]

            // [BEGIN MTDM-RUN-LIMITER]
            const bool limBypass = (pLimBypass01 != nullptr
                                      ? (pLimBypass01->load (std::memory_order_relaxed) >= 0.5f)
                                      : (levelscope::mtdm::Defaults::limBypass01 >= 0.5f));

            if (lim.enabled)
            {
                // Safe bypass preserves limiter latency/delay pipeline.
                lim.auditionBypass = limBypass;
                limiterStage.process (ctx.audio, lim);
            }
            // [END MTDM-RUN-LIMITER]
            // [END MTDM-PROCESS-LIMITER]

            // [BEGIN MTDM-LIMITER-METERING-UPDATE]
            // Metering: compute GR dB from limiter-applied gain (block min and last).
            // GR dB is reported as positive numbers (0 = no reduction).
            auto gainToGrDb = [] (float g) noexcept
            {
                g = juce::jlimit (1.0e-6f, 1.0f, g);
                return (float) (-20.0 * std::log10 ((double) g));
            };

            const float minGain  = limiterStage.lim.getLastBlockMinGain();
            const float lastGain = limiterStage.lim.getLastBlockLastGain();

            const float grBlockPeak = gainToGrDb (minGain);
            const float grCurrent   = gainToGrDb (lastGain);

            limiterMetering.grDbBlockPeak.store (grBlockPeak, std::memory_order_relaxed);
            limiterMetering.grDbCurrent.store   (grCurrent,   std::memory_order_relaxed);

            // Peak-hold with hold time + decay (audio-thread-only state).
            const int holdSamples = (int) std::lround (preparedSampleRate * (double) limiterHoldTimeSeconds);
            const float decayDbThisBlock =
                (preparedSampleRate > 1.0 ? limiterHoldDecayDbPerSecond * (float) ((double) ctx.numSamples / preparedSampleRate)
                                          : 0.0f);

            if (grBlockPeak > limiterHoldDbInternal)
            {
                limiterHoldDbInternal = grBlockPeak;
                limiterHoldSamplesLeft = holdSamples;
            }
            else
            {
                if (limiterHoldSamplesLeft > 0)
                {
                    limiterHoldSamplesLeft -= ctx.numSamples;
                    if (limiterHoldSamplesLeft < 0)
                        limiterHoldSamplesLeft = 0;
                }
                else
                {
                    limiterHoldDbInternal = std::max (0.0f, limiterHoldDbInternal - decayDbThisBlock);
                }
            }

            limiterMetering.grDbHold.store (limiterHoldDbInternal, std::memory_order_relaxed);
            // [END MTDM-LIMITER-METERING-UPDATE]
            // [BEGIN MTDM-ZONE-AUDITION-APPLY]
            // Apply the precomputed (latency-aligned) zone audition gate AFTER all processing.
            if (! zoneGateScratch.empty())
            {
                float* const* wr = ctx.audio.getArrayOfWritePointers();
                const int nCh = ctx.audio.getNumChannels();
                const int nS  = ctx.audio.getNumSamples();

                const float* g = zoneGateScratch.data();

                for (int ch = 0; ch < nCh; ++ch)
                    for (int i = 0; i < nS; ++i)
                        wr[ch][i] *= g[i];
            }
            // [END MTDM-ZONE-AUDITION-APPLY]

            juce::ignoreUnused (pThresholdDb, pRatio);
            // [END MTDM-PROCESS-UPWARD-MODE-SWITCH]

            // NOTE: Other zones (T1–T2, etc.) are pass-through for Stage D1a by simply not existing yet.
            // thresholdDb/ratio are currently unused placeholders for future time-domain zones.
        }
    // [END MTDM-PROCESS-STAGE-D1A]

    juce::ValueTree MultiThresholdDynamicsModule::saveState() const
    {
        auto t = AudioModuleBase::saveState();
        t.setProperty ("schemaVersion", 1, nullptr);
        return t;
    }

    void MultiThresholdDynamicsModule::loadState (const juce::ValueTree& state)
    {
        AudioModuleBase::loadState (state);
        // schemaVersion currently unused; tolerate older/missing fields.
    }
    // [END MTDM-MODULE-IMPL]
}
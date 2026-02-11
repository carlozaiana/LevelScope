#include "MultiThresholdDynamicsModule.h"
// [BEGIN MTDM-UPWARD-STRATEGY-INCLUDES]
#include <cmath>
// [END MTDM-UPWARD-STRATEGY-INCLUDES]

namespace levelscope
{
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
            comp.setLfePolicyAudioThread (rp.lfeInDetector, rp.lfeInApply);
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

        levelscope::dsp::LookaheadLimiter::Parameters p;
        p.enabled      = rp.enabled;
        p.ceilingDb    = rp.ceilingDb;
        p.lookaheadMs  = rp.lookaheadMs;
        p.releaseMs    = rp.releaseMs;

        lim.setParametersAudioThread (p);
        lim.process (audio);
    }

    int MultiThresholdDynamicsModule::LookaheadLimiterStage::getLatencySamples() const noexcept
    {
        return lim.getLatencySamples();
    }
    // [END MTDM-LIMITER-STRATEGY-IMPL]

    // [BEGIN MTDM-BINDPARAMS-STAGE-D1A-IMPL]
            // [BEGIN MTDM-BINDPARAMS-STAGE-D2A-IMPL]
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
                                                              std::atomic<float>* limReleaseMs) noexcept
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

                pUpwardModeChoice  = upwardModeChoice;
                pLfeInDetector01   = lfeInDetector01;
                pLfeInApply01      = lfeInApply01;

                pT2Lufs = t2Lufs;
                pT3Lufs = t3Lufs;

                pDownEnabled01 = downEnabled01;
                pDownRatio     = downRatio;
                pDownKneeDb    = downKneeDb;
                pDownAttackMs  = downAttackMs;
                pDownReleaseMs = downReleaseMs;
                pDownMakeupDb  = downMakeupDb;

                // [BEGIN MTDM-BINDPARAMS-STORE-LIMITER]
                pLimEnabled01   = limEnabled01;
                pLimCeilingDb   = limCeilingDb;
                pLimLookaheadMs = limLookaheadMs;
                pLimReleaseMs   = limReleaseMs;
                // [END MTDM-BINDPARAMS-STORE-LIMITER]
            }
            // [END MTDM-BINDPARAMS-STAGE-D2A-IMPL]

    // [BEGIN MTDM-PROCESS-STAGE-D1A]
        void MultiThresholdDynamicsModule::process (ProcessContext& ctx) noexcept
        {
            // Default is disabled => pass-through (no audible change, no latency insertion).
            const float enabled = (pEnabled01 != nullptr
                                     ? pEnabled01->load (std::memory_order_relaxed)
                                     : 0.0f);

            if (enabled < 0.5f)
                return;

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

            activeUpward->process (ctx.audio, up);
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

            downwardProcessor.process (ctx.audio, down);
            // [END MTDM-PROCESS-DOWNWARD]
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

            limiterStage.process (ctx.audio, lim);
            // [END MTDM-PROCESS-LIMITER]

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
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
                                                      std::atomic<float>* zoneUntouchedMute01) noexcept
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

            // [BEGIN MTDM-ZONE-SOLO-MUTE-LOGIC]
            const int solo = (pZoneSoloChoice != nullptr
                                ? (int) std::lround (pZoneSoloChoice->load (std::memory_order_relaxed))
                                : levelscope::mtdm::Defaults::zoneSoloChoice);

            const bool muteUp = (pZoneUpwardMute01 != nullptr
                                   ? (pZoneUpwardMute01->load (std::memory_order_relaxed) >= 0.5f)
                                   : (levelscope::mtdm::Defaults::zoneUpwardMute01 >= 0.5f));

            const bool muteDown = (pZoneDownwardMute01 != nullptr
                                     ? (pZoneDownwardMute01->load (std::memory_order_relaxed) >= 0.5f)
                                     : (levelscope::mtdm::Defaults::zoneDownwardMute01 >= 0.5f));

            const bool muteLim = (pZoneLimiterMute01 != nullptr
                                    ? (pZoneLimiterMute01->load (std::memory_order_relaxed) >= 0.5f)
                                    : (levelscope::mtdm::Defaults::zoneLimiterMute01 >= 0.5f));

            const bool runUp   = (! muteUp)   && (solo == 0 || solo == 1);
            const bool runDown = (! muteDown) && (solo == 0 || solo == 2);
            const bool runLim  = (! muteLim)  && (solo == 0 || solo == 3);
            // [BEGIN MTDM-UNTOUCHED-AUDITION-FLAGS]
            const bool soloUntouched = (solo == 4);

            const bool muteUntouched = (pZoneUntouchedMute01 != nullptr
                                          ? (pZoneUntouchedMute01->load (std::memory_order_relaxed) >= 0.5f)
                                          : (levelscope::mtdm::Defaults::zoneUntouchedMute01 >= 0.5f));

            // Gate mode precedence:
            // - Solo Untouched overrides everything.
            // - Untouched Mute only applies when not soloing a stage (solo==0).
            const bool gateSoloUntouched = soloUntouched;
            const bool gateMuteUntouched = (! soloUntouched) && (solo == 0) && muteUntouched;

            const bool auditionGateActive = gateSoloUntouched || gateMuteUntouched;
            // [END MTDM-UNTOUCHED-AUDITION-FLAGS]
            // [END MTDM-ZONE-SOLO-MUTE-LOGIC]

            // [BEGIN MTDM-RUN-UPWARD]
            // If Upward mode is Spectral, always call it to preserve STFT delay.
            // If muted/unsoloed, we run it in auditionBypass mode (unity gains).
            if (lastUpwardModeChoice == 0) // Spectral
            {
                up.auditionBypass = ! runUp;
                activeUpward->process (ctx.audio, up);
            }
            else
            {
                // Broadband upward has no latency; safe to skip.
                if (runUp)
                    activeUpward->process (ctx.audio, up);
            }
            // [END MTDM-RUN-UPWARD]

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

            // [BEGIN MTDM-RUN-DOWNWARD]
            if (runDown)
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
            if (lim.enabled)
            {
                lim.auditionBypass = ! runLim;
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
            // [BEGIN MTDM-UNTOUCHED-AUDITION-APPLY]
            // Audition gate: either "solo untouched zone" or "mute untouched zone".
            // Untouched zone is defined as T1 <= L < T2 (broadband LUFS-ish proxy).
            if (auditionGateActive)
            {
                // Reset detector smoothly when gate mode toggles to avoid popping.
                if (! auditionWasActive)
                {
                    auditionEnvMS = 0.0f;
                    auditionGateZ = 0.0f;
                }
                auditionWasActive = true;

                const float t1 = up.t1Lufs;
                float t2 = (pT2Lufs != nullptr ? pT2Lufs->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::t2Lufs);
                if (t2 < t1) t2 = t1;

                auto** chans = ctx.audio.getArrayOfWritePointers();
                const int numCh = ctx.audio.getNumChannels();
                const int numS  = ctx.audio.getNumSamples();

                // Decide detector channel set (exclude LFE unless lfeInDetector is true)
                const bool lfeDet = (pLfeInDetector01 != nullptr
                                       ? (pLfeInDetector01->load (std::memory_order_relaxed) >= 0.5f)
                                       : (levelscope::mtdm::Defaults::lfeInDetector01 >= 0.5f));

                for (int i = 0; i < numS; ++i)
                {
                    double sumSq = 0.0;
                    int count = 0;

                    for (int ch = 0; ch < numCh; ++ch)
                    {
                        const bool isLfe = (preparedChannelSet.getTypeOfChannel (ch) == juce::AudioChannelSet::LFE);

                        if (! lfeDet && isLfe)
                            continue;

                        const float x = chans[ch][i];
                        sumSq += (double) x * (double) x;
                        ++count;
                    }

                    const float e = (float) (sumSq / (double) juce::jmax (1, count));

                    const float a = (e > auditionEnvMS ? auditionDetA : auditionDetR);
                    auditionEnvMS = a * auditionEnvMS + (1.0f - a) * e;

                    const float L = (float) (-0.691 + 10.0 * std::log10 ((double) auditionEnvMS + 1.0e-12));

                    const bool inUntouched = (L >= t1 && L < t2);

                    const float targetGate =
                        gateSoloUntouched ? (inUntouched ? 1.0f : 0.0f)
                                          : (inUntouched ? 0.0f : 1.0f);

                    // Smooth the gate to avoid clicks
                    auditionGateZ = auditionGateA * auditionGateZ + (1.0f - auditionGateA) * targetGate;

                    for (int ch = 0; ch < numCh; ++ch)
                        chans[ch][i] *= auditionGateZ;
                }
            }
            else
            {
                auditionWasActive = false;
            }
            // [END MTDM-UNTOUCHED-AUDITION-APPLY]

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
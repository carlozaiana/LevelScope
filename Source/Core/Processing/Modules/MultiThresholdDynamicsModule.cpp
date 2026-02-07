#include "MultiThresholdDynamicsModule.h"

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

        // [BEGIN MTDM-PREPARE-SUC]
                // Stage D1a: prepare spectral upward compressor (allocations happen here only).
                const int numCh = juce::jmax (1, spec.channelSet.size());
                spectralUpward.prepare (spec.sampleRate, numCh, spec.channelSet, spec.maxBlockSize);
                spectralPrepared = true;
        // [END MTDM-PREPARE-SUC]
    }

    // [BEGIN MTDM-RESET-STAGE-D1A]
        void MultiThresholdDynamicsModule::reset()
        {
            if (spectralPrepared)
                spectralUpward.reset();
        }
    // [END MTDM-RESET-STAGE-D1A]

    // [BEGIN MTDM-BINDPARAMS-STAGE-D1A-IMPL]
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
                                                              std::atomic<float>* upwardModeChoice) noexcept
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

                // [BEGIN MTDM-BINDPARAMS-STORE-UPWARD-MODE]
                pUpwardModeChoice = upwardModeChoice;
                // [END MTDM-BINDPARAMS-STORE-UPWARD-MODE]
            }
    // [END MTDM-BINDPARAMS-STAGE-D1A-IMPL]

    // [BEGIN MTDM-PROCESS-STAGE-D1A]
        void MultiThresholdDynamicsModule::process (ProcessContext& ctx) noexcept
        {
            // Default is disabled => pass-through (no audible change, no latency insertion).
            const float enabled = (pEnabled01 != nullptr
                                     ? pEnabled01->load (std::memory_order_relaxed)
                                     : 0.0f);

            if (enabled < 0.5f)
                return;

            // Must be prepared before processing.
            if (! spectralPrepared)
                return;

            // Read parameters (audio-thread-only contract).
            levelscope::dsp::SpectralUpwardCompressor::Parameters p;

            // UI-oriented thresholds (LUFS axis); internally translated by spectralUpward
            const float t0 = (pT0Lufs != nullptr ? pT0Lufs->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::t0Lufs);
            const float t1 = (pT1Lufs != nullptr ? pT1Lufs->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::t1Lufs);

            // Safety: enforce ordering (UI should constrain too).
            p.t0Lufs = juce::jmin (t0, t1);
            p.t1Lufs = juce::jmax (t0, t1);

            p.amount01   = (pSucAmount01   != nullptr ? pSucAmount01->load   (std::memory_order_relaxed) : levelscope::mtdm::Defaults::sucAmount01);
            p.maxBoostDb = (pSucMaxBoostDb != nullptr ? pSucMaxBoostDb->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::sucMaxBoostDb);
            p.curve      = (pSucCurve      != nullptr ? pSucCurve->load      (std::memory_order_relaxed) : levelscope::mtdm::Defaults::sucCurve);

            p.lowKneeDb  = (pSucLowKneeDb  != nullptr ? pSucLowKneeDb->load  (std::memory_order_relaxed) : levelscope::mtdm::Defaults::sucLowKneeDb);
            p.highKneeDb = (pSucHighKneeDb != nullptr ? pSucHighKneeDb->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::sucHighKneeDb);

            p.attackMs   = (pSucAttackMs   != nullptr ? pSucAttackMs->load   (std::memory_order_relaxed) : levelscope::mtdm::Defaults::sucAttackMs);
            p.releaseMs  = (pSucReleaseMs  != nullptr ? pSucReleaseMs->load  (std::memory_order_relaxed) : levelscope::mtdm::Defaults::sucReleaseMs);

            // Advanced choices (choice params come through as float indices)
            const int fftChoice = (int) std::lround (pSucFftSizeChoice != nullptr
                                                     ? pSucFftSizeChoice->load (std::memory_order_relaxed)
                                                     : (float) levelscope::mtdm::Defaults::sucFftSizeChoice);

            const int bpoChoice = (int) std::lround (pSucBandsPerOctChoice != nullptr
                                                     ? pSucBandsPerOctChoice->load (std::memory_order_relaxed)
                                                     : (float) levelscope::mtdm::Defaults::sucBandsPerOctChoice);

            p.fftSizeChoice     = fftChoice;
            p.bandsPerOctChoice = bpoChoice;

            p.minFreqHz = (pSucMinFreqHz != nullptr ? pSucMinFreqHz->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::sucMinFreqHz);
            p.maxFreqHz = (pSucMaxFreqHz != nullptr ? pSucMaxFreqHz->load (std::memory_order_relaxed) : levelscope::mtdm::Defaults::sucMaxFreqHz);

                p.minFreqHz = juce::jlimit (10.0f, 2000.0f, p.minFreqHz);
                p.maxFreqHz = juce::jlimit (1000.0f, 24000.0f, p.maxFreqHz);

            // Ensure min/max freq ordering
            if (p.maxFreqHz < p.minFreqHz + 10.0f)
                p.maxFreqHz = p.minFreqHz + 10.0f;

            // [BEGIN MTDM-SUC-TRIM-AND-CURVETYPE-APPLY]
                p.calibrationTrimDb =
                    (pSucCalTrimDb != nullptr ? pSucCalTrimDb->load (std::memory_order_relaxed)
                        : levelscope::mtdm::Defaults::sucCalTrimDb);

                const int curveTypeChoice = (int) std::lround (pSucCurveTypeChoice != nullptr
                        ? pSucCurveTypeChoice->load (std::memory_order_relaxed)
                        : (float) levelscope::mtdm::Defaults::sucCurveTypeChoice);

                p.curveType = (curveTypeChoice == 1
                        ? levelscope::dsp::SpectralUpwardCompressor::CurveType::bell
                        : levelscope::dsp::SpectralUpwardCompressor::CurveType::monotonic);

                p.fftSizeChoice = juce::jlimit (0, 3, fftChoice);
                p.bandsPerOctChoice = juce::jlimit (0, 3, bpoChoice);

            // [END MTDM-SUC-TRIM-AND-CURVETYPE-APPLY]

            // Apply
            spectralUpward.setParametersAudioThread (p);

            // Process in-place (safe)
            spectralUpward.process (ctx.audio);

            // NOTE: Other zones (T1–T2, etc.) are pass-through for Stage D1a by simply not existing yet.
            // thresholdDb/ratio are currently unused placeholders for future time-domain zones.
            juce::ignoreUnused (pThresholdDb, pRatio);
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
} // namespace levelscope
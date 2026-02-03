#pragma once

#include "../IAudioModule.h"
// [BEGIN MTDM-PARAM-HEADER-INCLUDE]
#include "MultiThresholdDynamicsParamIDs.h"
// [END MTDM-PARAM-HEADER-INCLUDE]
// [BEGIN MTDM-SUC-INCLUDE]
#include "../DSP/SpectralUpwardCompressor.h"
// [END MTDM-SUC-INCLUDE]

namespace levelscope
{
    // [BEGIN MTDM-MODULE-DECL]
    // MultiThresholdDynamicsModule (Stage B skeleton)
    // - No DSP yet: process() is a strict no-op (pass-through).
    // - RT-safety: process() does not allocate, lock, or touch ValueTree.
    class MultiThresholdDynamicsModule final : public AudioModuleBase
    {
    public:
        // Stable persistence ID (must not change once released).
        static constexpr const char* kModuleID = "MTDM"; // Multi-Threshold Dynamics Module

        juce::String getModuleID() const override;
        juce::String getDisplayName() const override;

        void prepare (const ModulePrepareSpec& spec) override;
        void reset() override;

        // [BEGIN MTDM-RT-PARAM-BINDING-DECL]
            void process (ProcessContext& ctx) noexcept override;

            // Parameter binding (non-owning pointers; RT-safe reads in process()).
            // Caller must ensure pointers remain valid for module lifetime (APVTS raw param pointers do).
            // [BEGIN MTDM-BINDPARAMS-STAGE-D1A-DECL]
                        void bindParameters (std::atomic<float>* enabled01,
                                             std::atomic<float>* thresholdDb,
                                             std::atomic<float>* ratio,

                                             // Stage D1a params (additive)
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
                                             std::atomic<float>* sucMaxFreqHz) noexcept;
            // [END MTDM-BINDPARAMS-STAGE-D1A-DECL]

            // Persistence (non-audio-thread only)
            juce::ValueTree saveState() const override;
            void loadState (const juce::ValueTree& state) override;
        // [END MTDM-RT-PARAM-BINDING-DECL]

        // [BEGIN MTDM-RT-PARAM-BINDING-MEMBERS]
    private:
        // Non-owning APVTS parameter pointers (RT-safe to read; never write from audio thread)
        std::atomic<float>* pEnabled01     = nullptr; // 0/1
        std::atomic<float>* pThresholdDb   = nullptr;
        std::atomic<float>* pRatio         = nullptr;

        // Stored from prepare() for debugging/future DSP wiring. Not used for DSP yet.
        double preparedSampleRate = 0.0;
        int preparedMaxBlockSize = 0;
        juce::AudioChannelSet preparedChannelSet;
        // [END MTDM-RT-PARAM-BINDING-MEMBERS]

        // [BEGIN MTDM-SUC-MEMBERS]
                // Stage D1a parameter pointers
                std::atomic<float>* pT0Lufs = nullptr;
                std::atomic<float>* pT1Lufs = nullptr;

                std::atomic<float>* pSucAmount01        = nullptr;
                std::atomic<float>* pSucMaxBoostDb      = nullptr;
                std::atomic<float>* pSucCurve           = nullptr;
                std::atomic<float>* pSucLowKneeDb       = nullptr;
                std::atomic<float>* pSucHighKneeDb      = nullptr;
                std::atomic<float>* pSucAttackMs        = nullptr;
                std::atomic<float>* pSucReleaseMs       = nullptr;

                std::atomic<float>* pSucFftSizeChoice      = nullptr; // choice index in float form
                std::atomic<float>* pSucBandsPerOctChoice  = nullptr; // choice index in float form
                std::atomic<float>* pSucMinFreqHz          = nullptr;
                std::atomic<float>* pSucMaxFreqHz          = nullptr;

                // Stage D1a DSP
                levelscope::dsp::SpectralUpwardCompressor spectralUpward;
                bool spectralPrepared = false;
        // [END MTDM-SUC-MEMBERS]
    };
    // [END MTDM-MODULE-DECL]
} // namespace levelscope
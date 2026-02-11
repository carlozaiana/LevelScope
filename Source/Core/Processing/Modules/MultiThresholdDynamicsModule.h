#pragma once

#include "../IAudioModule.h"
// [BEGIN MTDM-PARAM-HEADER-INCLUDE]
#include "MultiThresholdDynamicsParamIDs.h"
// [END MTDM-PARAM-HEADER-INCLUDE]
// [BEGIN MTDM-SUC-INCLUDE]
#include "../DSP/SpectralUpwardCompressor.h"
// [END MTDM-SUC-INCLUDE]
// [BEGIN MTDM-BUC-INCLUDE]
#include "../DSP/BroadbandUpwardCompressor.h"
// [END MTDM-BUC-INCLUDE]
// [BEGIN MTDM-BDC-INCLUDE]
#include "../DSP/BroadbandDownwardCompressor.h"
// [END MTDM-BDC-INCLUDE]

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
                                     std::atomic<float>* sucMaxFreqHz,

                                     // Stage D1a usability
                                     std::atomic<float>* sucCalTrimDb,
                                     std::atomic<float>* sucCurveTypeChoice,

                                     // Stage D1c ADD-UPWARD-MODE
                                     std::atomic<float>* upwardModeChoice,
                                     // [BEGIN MTDM-BINDPARAMS-ADD-LFE-MASK]
                                     std::atomic<float>* lfeInDetector01,
                                     std::atomic<float>* lfeInApply01,
                                     // [END MTDM-BINDPARAMS-ADD-LFE-MASK]
                                     // [BEGIN MTDM-BINDPARAMS-ADD-DOWNWARD]
                                     // Stage D2a: downward zone params
                                     std::atomic<float>* t2Lufs,
                                     std::atomic<float>* t3Lufs,
                                     std::atomic<float>* downEnabled01,
                                     std::atomic<float>* downRatio,
                                     std::atomic<float>* downKneeDb,
                                     std::atomic<float>* downAttackMs,
                                     std::atomic<float>* downReleaseMs,
                                     std::atomic<float>* downMakeupDb,

                                     // Stage D3a: limiter params
                                     std::atomic<float>* limEnabled01,
                                     std::atomic<float>* limCeilingDb,
                                     std::atomic<float>* limLookaheadMs,
                                     std::atomic<float>* limReleaseMs) noexcept;
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

        // [BEGIN MTDM-UPWARD-STRATEGY-TYPES]
        struct UpwardRuntimeParams
        {
            float t0Lufs = levelscope::mtdm::Defaults::t0Lufs;
            float t1Lufs = levelscope::mtdm::Defaults::t1Lufs;

            float amount01   = levelscope::mtdm::Defaults::sucAmount01;
            float maxBoostDb = levelscope::mtdm::Defaults::sucMaxBoostDb;
            float curve      = levelscope::mtdm::Defaults::sucCurve;

            float lowKneeDb  = levelscope::mtdm::Defaults::sucLowKneeDb;
            float highKneeDb = levelscope::mtdm::Defaults::sucHighKneeDb;

            float attackMs   = levelscope::mtdm::Defaults::sucAttackMs;
            float releaseMs  = levelscope::mtdm::Defaults::sucReleaseMs;

            // 0=Monotonic, 1=Bell (shared across strategies)
            int curveTypeChoice = levelscope::mtdm::Defaults::sucCurveTypeChoice;

            // Spectral-only:
            int fftSizeChoice     = levelscope::mtdm::Defaults::sucFftSizeChoice;
            int bandsPerOctChoice = levelscope::mtdm::Defaults::sucBandsPerOctChoice;
            float minFreqHz       = levelscope::mtdm::Defaults::sucMinFreqHz;
            float maxFreqHz       = levelscope::mtdm::Defaults::sucMaxFreqHz;
            float calibrationTrimDb = levelscope::mtdm::Defaults::sucCalTrimDb;

            // [BEGIN MTDM-UPWARD-RP-LFE-MASK]
            bool lfeInDetector = false;
            bool lfeInApply    = false;
            // [END MTDM-UPWARD-RP-LFE-MASK]
        };

        struct IUpwardProcessor
        {
            virtual ~IUpwardProcessor() = default;
            virtual void prepare (const ModulePrepareSpec& spec) = 0;
            virtual void reset() noexcept = 0;
            virtual void process (juce::AudioBuffer<float>& audio, const UpwardRuntimeParams& p) noexcept = 0;
            virtual int  getLatencySamples() const noexcept = 0;
        };

        struct SpectralUpwardProcessor final : IUpwardProcessor
        {
            void prepare (const ModulePrepareSpec& spec) override;
            void reset() noexcept override;
            void process (juce::AudioBuffer<float>& audio, const UpwardRuntimeParams& p) noexcept override;
            int  getLatencySamples() const noexcept override;

            levelscope::dsp::SpectralUpwardCompressor suc;
            bool prepared = false;
        };

        struct BroadbandUpwardProcessor final : IUpwardProcessor
        {
            void prepare (const ModulePrepareSpec& spec) override;
            void reset() noexcept override;
            void process (juce::AudioBuffer<float>& audio, const UpwardRuntimeParams& p) noexcept override;
            int  getLatencySamples() const noexcept override;

            levelscope::dsp::BroadbandUpwardCompressor buc;
            bool prepared = false;
        };

        SpectralUpwardProcessor  spectralUpwardProcessor;
        BroadbandUpwardProcessor broadbandUpwardProcessor;

        IUpwardProcessor* activeUpward = nullptr;
        int lastUpwardModeChoice = -1; // 0=Spectral, 1=Broadband
        
        // [BEGIN MTDM-DOWNWARD-STRATEGY-TYPES]
        struct DownwardRuntimeParams
        {
            bool enabled = false;

            float t2Lufs = levelscope::mtdm::Defaults::t2Lufs;
            float t3Lufs = levelscope::mtdm::Defaults::t3Lufs;

            float ratio    = levelscope::mtdm::Defaults::downRatio;
            float kneeDb   = levelscope::mtdm::Defaults::downKneeDb;
            float attackMs = levelscope::mtdm::Defaults::downAttackMs;
            float releaseMs = levelscope::mtdm::Defaults::downReleaseMs;
            float makeupDb = levelscope::mtdm::Defaults::downMakeupDb;

            bool lfeInDetector = false;
            bool lfeInApply    = false;
        };

        struct IDownwardProcessor
        {
            virtual ~IDownwardProcessor() = default;
            virtual void prepare (const ModulePrepareSpec& spec) = 0;
            virtual void reset() noexcept = 0;
            virtual void process (juce::AudioBuffer<float>& audio, const DownwardRuntimeParams& p) noexcept = 0;
        };

        struct BroadbandDownwardProcessor final : IDownwardProcessor
        {
            void prepare (const ModulePrepareSpec& spec) override;
            void reset() noexcept override;
            void process (juce::AudioBuffer<float>& audio, const DownwardRuntimeParams& p) noexcept override;

            levelscope::dsp::BroadbandDownwardCompressor comp;
            bool prepared = false;
        };

        BroadbandDownwardProcessor downwardProcessor;
        // [END MTDM-DOWNWARD-STRATEGY-TYPES]
        // [END MTDM-UPWARD-STRATEGY-TYPES]

        // [BEGIN MTDM-LIMITER-STRATEGY-TYPES]
        struct LimiterRuntimeParams
        {
            bool enabled = false;
            float ceilingDb   = levelscope::mtdm::Defaults::limCeilingDb;
            float lookaheadMs = levelscope::mtdm::Defaults::limLookaheadMs;
            float releaseMs   = levelscope::mtdm::Defaults::limReleaseMs;
        };

        struct ILimiter
        {
            virtual ~ILimiter() = default;
            virtual void prepare (const ModulePrepareSpec& spec) = 0;
            virtual void reset() noexcept = 0;
            virtual void process (juce::AudioBuffer<float>& audio, const LimiterRuntimeParams& p) noexcept = 0;
            virtual int  getLatencySamples() const noexcept = 0;
        };

        struct LookaheadLimiterStage final : ILimiter
        {
            void prepare (const ModulePrepareSpec& spec) override;
            void reset() noexcept override;
            void process (juce::AudioBuffer<float>& audio, const LimiterRuntimeParams& p) noexcept override;
            int  getLatencySamples() const noexcept override;

            levelscope::dsp::LookaheadLimiter lim;
            bool prepared = false;
        };

        LookaheadLimiterStage limiterStage;
        // [END MTDM-LIMITER-STRATEGY-TYPES]

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

        // [BEGIN MTDM-SUC-TRIM-AND-CURVETYPE-MEMBERS]
                std::atomic<float>* pSucCalTrimDb       = nullptr;
                std::atomic<float>* pSucCurveTypeChoice = nullptr; // 0=Monotonic, 1=Bell
        // [END MTDM-SUC-TRIM-AND-CURVETYPE-MEMBERS]

        // [BEGIN MTDM-UPWARD-MODE-MEMBER]
                std::atomic<float>* pUpwardModeChoice = nullptr; // 0=Spectral, 1=Broadband
        // [END MTDM-UPWARD-MODE-MEMBER]

        // [BEGIN MTDM-LFE-MASK-MEMBERS]
                std::atomic<float>* pLfeInDetector01 = nullptr; // 0/1
                std::atomic<float>* pLfeInApply01    = nullptr; // 0/1
        // [END MTDM-LFE-MASK-MEMBERS]

        // [BEGIN MTDM-DOWNWARD-PARAM-PTRS]
        std::atomic<float>* pT2Lufs = nullptr;
        std::atomic<float>* pT3Lufs = nullptr;

        std::atomic<float>* pDownEnabled01 = nullptr;
        std::atomic<float>* pDownRatio     = nullptr;
        std::atomic<float>* pDownKneeDb    = nullptr;
        std::atomic<float>* pDownAttackMs  = nullptr;
        std::atomic<float>* pDownReleaseMs = nullptr;
        std::atomic<float>* pDownMakeupDb  = nullptr;
        // [END MTDM-DOWNWARD-PARAM-PTRS]

        // [BEGIN MTDM-LIMITER-PARAM-PTRS]
        std::atomic<float>* pLimEnabled01   = nullptr;
        std::atomic<float>* pLimCeilingDb   = nullptr;
        std::atomic<float>* pLimLookaheadMs = nullptr;
        std::atomic<float>* pLimReleaseMs   = nullptr;
        // [END MTDM-LIMITER-PARAM-PTRS]
    };
    // [END MTDM-MODULE-DECL]
} // namespace levelscope
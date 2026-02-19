#pragma once

#include "../IAudioModule.h"
// [BEGIN MTDM-INCLUDE-CSTDINT]
#include <cstdint>
// [END MTDM-INCLUDE-CSTDINT]
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
// [BEGIN MTDM-LIMITER-INCLUDE]
#include "../DSP/LookaheadLimiter.h"
// [END MTDM-LIMITER-INCLUDE]

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

        // [BEGIN MTDM-LIMITER-METERING-PUBLIC]
        struct LimiterMetering
        {
            // All values are "gain reduction in dB" as positive numbers.
            // 0.0 = no reduction, 6.0 = 6 dB reduction, etc.
            std::atomic<float> grDbCurrent   { 0.0f }; // current/instant (updated once per block)
            std::atomic<float> grDbBlockPeak { 0.0f }; // peak GR seen in last processed block
            std::atomic<float> grDbHold      { 0.0f }; // peak-hold with decay
        };

        const LimiterMetering& getLimiterMetering() const noexcept { return limiterMetering; }
        // [END MTDM-LIMITER-METERING-PUBLIC]

        // [BEGIN MTDM-DOWNWARD-METERING-PUBLIC]
        struct DownwardMetering
        {
            // Gain reduction in dB as positive numbers (0 = no reduction).
            std::atomic<float> grDbCurrent   { 0.0f };
            std::atomic<float> grDbBlockPeak { 0.0f };
            std::atomic<float> grDbHold      { 0.0f };
        };

        const DownwardMetering& getDownwardMetering() const noexcept { return downwardMetering; }
        // [END MTDM-DOWNWARD-METERING-PUBLIC]

        void prepare (const ModulePrepareSpec& spec) override;
        void reset() override;

        // [BEGIN MTDM-RT-PARAM-BINDING-DECL]
            void process (ProcessContext& ctx) noexcept override;

            // Parameter binding (non-owning pointers; RT-safe reads in process()).
            // Caller must ensure pointers remain valid for module lifetime (APVTS raw param pointers do).
            // [BEGIN MTDM-BINDPARAMS-FULL-DECL]
            void bindParameters (std::atomic<float>* enabled01,
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
                                 // [BEGIN MTDM-MC-POLICY-BINDPARAMS-DECL]
                                 std::atomic<float>* mcPolicyChoice,
                                 std::atomic<float>* dialogDetectorChoice,
                                 std::atomic<float>* dialogApplyChoice) noexcept;
                                 // [END MTDM-MC-POLICY-BINDPARAMS-DECL]
            // [END MTDM-BINDPARAMS-FULL-DECL]

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

        // [BEGIN MTDM-MC-POLICY-TYPES-AND-MASK-STATE]
        enum class MCPolicy : int
        {
            linked     = 0,
            dialogMask = 1,
            unlinked   = 2
        };

        static constexpr int kMaxPolicyChannels = 16; // enough for 7.1.4 (12ch) + headroom

        // Current masks (recomputed on audio thread when dirty)
        uint16_t detectMaskBits = 0; // bit i => channel i included in detector set
        uint16_t applyMaskBits  = 0; // bit i => channel i included in apply set

        // Cached inputs to detect changes (RT-safe; updated only on audio thread)
        int lastPolicyChoice = -1;
        int lastDialogDetectorChoice = -1;
        int lastDialogApplyChoice    = -1;
        bool lastLfeInDetector = false;
        bool lastLfeInApply    = false;
        int lastMaskNumChannels = -1;
        juce::AudioChannelSet lastMaskChannelSet; // value type; safe to store/copy

        // AUDIO THREAD ONLY: recompute detect/apply masks if needed (no alloc)
        void updateChannelMasksIfNeeded (const ProcessContext& ctx,
                                         int policyChoice,
                                         int dialogDetectorChoice,
                                         int dialogApplyChoice,
                                         bool lfeInDetector,
                                         bool lfeInApply) noexcept;
        // [END MTDM-MC-POLICY-TYPES-AND-MASK-STATE]

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

            // [BEGIN MTDM-UPWARD-RP-MC-MASKBITS]
            // Stage E multichannel masks (interpreted by DSP blocks in Part 3/4)
            uint16_t detectMaskBits = 0;
            uint16_t applyMaskBits  = 0;
            bool     unlinked       = false;
            // [END MTDM-UPWARD-RP-MC-MASKBITS]

            // [BEGIN MTDM-UPWARD-AUDITION-BYPASS]
            bool auditionBypass = false; // delay-preserving unity mode (Spectral only)
            // [END MTDM-UPWARD-AUDITION-BYPASS]
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

            // [BEGIN MTDM-DOWNWARD-RP-MC-MASKBITS]
            uint16_t detectMaskBits = 0;
            uint16_t applyMaskBits  = 0;
            bool     unlinked       = false;
            // [END MTDM-DOWNWARD-RP-MC-MASKBITS]
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
            // [BEGIN MTDM-LIM-RUNTIME-TP]
            float attackMs = levelscope::mtdm::Defaults::limAttackMs;
            float driveDb  = levelscope::mtdm::Defaults::limDriveDb;
            int   oversamplingChoice = levelscope::mtdm::Defaults::limOversamplingChoice;
            // [END MTDM-LIM-RUNTIME-TP]
            // [BEGIN MTDM-LIMITER-AUDITION-BYPASS]
            bool auditionBypass = false; // delay-preserving unity mode
            // [END MTDM-LIMITER-AUDITION-BYPASS]
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

        // [BEGIN MTDM-LIMITER-METERING-MEMBERS]
        LimiterMetering limiterMetering;

        // [BEGIN MTDM-DOWNWARD-METERING-MEMBERS]
        DownwardMetering downwardMetering;

        float downwardHoldDbInternal = 0.0f;
        int   downwardHoldSamplesLeft = 0;

        static constexpr float downwardHoldTimeSeconds = 0.25f;
        static constexpr float downwardHoldDecayDbPerSecond = 12.0f;
        // [END MTDM-DOWNWARD-METERING-MEMBERS]

        // Audio-thread-only state for peak-hold behavior (no atomics needed here).
        float limiterHoldDbInternal = 0.0f;
        int   limiterHoldSamplesLeft = 0;

        static constexpr float limiterHoldTimeSeconds = 0.25f;   // hold before decay
        static constexpr float limiterHoldDecayDbPerSecond = 12.0f; // pro-style decay rate
        // [END MTDM-LIMITER-METERING-MEMBERS]

        // [BEGIN MTDM-UNTOUCHED-AUDITION-DETECTOR]
        // Audition gate detector (broadband LUFS-ish proxy) for "Untouched Solo/Mute".
        float auditionEnvMS = 0.0f;
        float auditionGateZ = 1.0f;

        float auditionDetA = 0.99f; // per-sample one-pole coeffs (set in prepare)
        float auditionDetR = 0.999f;
        float auditionGateA = 0.99f; // gate smoothing (set in prepare)

        bool auditionWasActive = false;
        // [END MTDM-UNTOUCHED-AUDITION-DETECTOR]

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

        // [BEGIN MTDM-MC-POLICY-PARAM-PTRS]
        std::atomic<float>* pMcPolicyChoice       = nullptr; // choice index in float form
        std::atomic<float>* pDialogDetectorChoice = nullptr; // 0=C, 1=LCR
        std::atomic<float>* pDialogApplyChoice    = nullptr; // 0=C, 1=LCR
        // [END MTDM-MC-POLICY-PARAM-PTRS]

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

        // [BEGIN MTDM-LIMITER-PARAM-PTRS-TP]
        std::atomic<float>* pLimAttackMs          = nullptr;
        std::atomic<float>* pLimDriveDb           = nullptr;
        std::atomic<float>* pLimOversamplingChoice = nullptr; // 0=Off,1=2x,2=4x
        // [END MTDM-LIMITER-PARAM-PTRS-TP]

        // [BEGIN MTDM-ZONE-SOLO-MUTE-PTRS]
        std::atomic<float>* pZoneSoloChoice = nullptr; // 0=None,1=Upward,2=Downward,3=Limiter
        std::atomic<float>* pZoneUpwardMute01 = nullptr;
        std::atomic<float>* pZoneDownwardMute01 = nullptr;
        std::atomic<float>* pZoneLimiterMute01 = nullptr;
        // [BEGIN MTDM-ZONE-UNTOUCHED-MUTE-PTR]
        std::atomic<float>* pZoneUntouchedMute01 = nullptr;
        // [END MTDM-ZONE-UNTOUCHED-MUTE-PTR]
        // [END MTDM-ZONE-SOLO-MUTE-PTRS]
    };
    // [END MTDM-MODULE-DECL]
} // namespace levelscope
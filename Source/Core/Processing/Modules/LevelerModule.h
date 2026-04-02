#pragma once

#include "../IAudioModule.h"
#include "../../BS1770KWeighting.h"
#include "LevelerParamIDs.h"

#include <array>
#include <atomic>
#include <cstdint>
#include <limits>

//==============================================================================
// LevelerModule
//
// v1 design goals:
// - zero added latency
// - disabled by default
// - K-weighted loudness-like detector
// - 60 Hz internal detector frames
// - supports Adaptive and Learn-Hold
// - supports Auto / Momentary / Short-term measurement selection
// - reuses shared multichannel routing policy params from MTDM
// - RT-safe processing: no allocations, no locks in process()
//==============================================================================

namespace levelscope
{
    // [BEGIN LVLR-MODULE-CLASS]
    class LevelerModule final : public IAudioModule
    {
    public:
        static constexpr int maxSupportedChannels = 16;
        static constexpr double detectorFrameRateHz = 60.0;
        static constexpr int momentaryWindowFrames = 24;   // 0.4 s @ 60 Hz
        static constexpr int shortTermWindowFrames = 180;  // 3.0 s @ 60 Hz

        enum MeasurementChoice
        {
            measAuto = 0,
            measMomentary,
            measShortTerm
        };

        enum ModeChoice
        {
            modeAdaptive = 0,
            modeLearnHold
        };

        enum ControlModeChoice
        {
            controlInternal = 0,
            controlHostGain
        };

        struct Metering
        {
            std::atomic<float> gainDbCurrent      { 0.0f };
            std::atomic<float> gainDbBlockPeak    { 0.0f };
            std::atomic<float> gainDbHold         { 0.0f };

            // Current measurement-selected detector loudness (LUFS-ish):
            // Auto -> ST if valid else M
            // Momentary -> M
            // Short-term -> ST with internal fallback until valid
            std::atomic<float> measuredLufsCurrent { -200.0f };
        };

        LevelerModule();
        ~LevelerModule() override = default;

        juce::String getModuleID() const override      { return "leveler"; }
        juce::String getDisplayName() const override   { return "Leveler"; }

        void prepare (const ModulePrepareSpec& spec) override;
        void reset() override;
        void process (ProcessContext& ctx) noexcept override;

        void setBypassed (bool shouldBeBypassed) noexcept override
        {
            bypassed.store (shouldBeBypassed, std::memory_order_relaxed);
        }

        bool isBypassed() const noexcept override
        {
            return bypassed.load (std::memory_order_relaxed);
        }

        // [BEGIN LVLR-VALUETREE-PERSISTENCE-SIGNATURE]
        // Leveler currently relies on APVTS-backed params only.
        // Module-local persistence is minimal for v1.
        juce::ValueTree saveState() const override;
        void loadState (const juce::ValueTree& state) override;
        // [END LVLR-VALUETREE-PERSISTENCE-SIGNATURE]

        void bindParameters (std::atomic<float>* enabledParam,
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
                             std::atomic<float>* lfeInApplyParam) noexcept;

        const Metering& getMetering() const noexcept { return metering; }

    private:
        //--------------------------------------------------------------------------
        // Raw APVTS parameter pointers (owned by plugin/APVTS, not by this module)
        std::atomic<float>* pEnabled01             = nullptr;
        std::atomic<float>* pTargetLufs            = nullptr;
        std::atomic<float>* pMaxBoostDb            = nullptr;
        std::atomic<float>* pMaxCutDb              = nullptr;
        std::atomic<float>* pMeasChoice            = nullptr;
        std::atomic<float>* pModeChoice            = nullptr;
        std::atomic<float>* pLearn01               = nullptr;
        std::atomic<float>* pRateUpDbPerSec        = nullptr;
        std::atomic<float>* pRateDownDbPerSec      = nullptr;
        std::atomic<float>* pControlModeChoice     = nullptr;
        std::atomic<float>* pHostGainDb            = nullptr;

        // Reused shared routing policy params (currently under mtdm.* namespace)
        std::atomic<float>* pMcPolicyChoice        = nullptr;
        std::atomic<float>* pDialogDetectorChoice  = nullptr;
        std::atomic<float>* pDialogApplyChoice     = nullptr;
        std::atomic<float>* pLfeInDetector         = nullptr;
        std::atomic<float>* pLfeInApply            = nullptr;

        //--------------------------------------------------------------------------
        // Runtime configuration
        double sampleRate = 48000.0;
        int preparedNumChannels = 0;
        juce::AudioChannelSet preparedChannelSet;

        std::atomic<bool> bypassed { false };

        //--------------------------------------------------------------------------
        // Detector path (self-contained; does not depend on plugin timeline history)
        BS1770KWeighting kWeight;

        int frameSamples = 800;              // sampleRate / 60, clamped in prepare()
        int samplesUntilNextFrame = 800;

        // Per-channel K-weighted energy accumulation for current 60 Hz frame
        std::array<double, maxSupportedChannels> frameEnergyAccumByChannel {};

        // Rolling detector-history at 60 Hz, per channel
        std::array<std::array<float, shortTermWindowFrames>, maxSupportedChannels> frameEnergyHistoryByChannel {};
        std::array<float, maxSupportedChannels> momentaryEnergySumsByChannel {};
        std::array<float, maxSupportedChannels> shortTermEnergySumsByChannel {};

        int detectorHistoryWriteIndex = 0;
        int detectorFramesWritten = 0;

        // Latest measured loudness estimates (for debug/metering/decision paths)
        float latestMomentaryLufs = -200.0f;
        float latestShortTermLufs = -200.0f;
        bool momentaryValid = false;
        bool shortTermValid = false;

        //--------------------------------------------------------------------------
        // Multichannel routing cache
        uint16_t detectMaskBits = 0;
        uint16_t applyMaskBits  = 0;

        int  cachedMcPolicyChoice       = std::numeric_limits<int>::min();
        int  cachedDialogDetectorChoice = std::numeric_limits<int>::min();
        int  cachedDialogApplyChoice    = std::numeric_limits<int>::min();
        int  cachedLfeInDetector01      = std::numeric_limits<int>::min();
        int  cachedLfeInApply01         = std::numeric_limits<int>::min();
        bool routingMasksDirty          = true;

        //--------------------------------------------------------------------------
        // Gain state
        std::array<float, maxSupportedChannels> currentGainDbByChannel {};
        std::array<float, maxSupportedChannels> desiredGainDbByChannel {};
        std::array<float, maxSupportedChannels> learnedCandidateGainDbByChannel {};
        std::array<float, maxSupportedChannels> heldGainDbByChannel {};

        bool wasPlayingLastBlock = false;
        bool learnWasOnLastBlock = false;
        bool haveCommittedHoldGain = false;

        Metering metering;

        //--------------------------------------------------------------------------
        // Helpers
        static float readParamOrDefault (const std::atomic<float>* p, float fallback) noexcept;
        static int   readChoiceOrDefault (const std::atomic<float>* p, int fallback) noexcept;
        static uint16_t bitForChannel (int ch) noexcept;

        void resetDetectorState() noexcept;
        void resetGainState() noexcept;

        void updateRoutingMasksIfNeeded (const juce::AudioChannelSet& channelSet,
                                         int numChannels) noexcept;

        void rebuildRoutingMasks (const juce::AudioChannelSet& channelSet,
                                  int numChannels,
                                  int mcPolicyChoice,
                                  int dialogDetectorChoice,
                                  int dialogApplyChoice,
                                  bool lfeInDetector,
                                  bool lfeInApply) noexcept;

        void pushDetectorFrameFromAccumulatedEnergy (int numChannels) noexcept;

        float computeMaskedMomentaryLufs (uint16_t maskBits, int numChannels) const noexcept;
        float computeMaskedShortTermLufs (uint16_t maskBits, int numChannels) const noexcept;
        float computeMaskedCurrentMeasurementLufs (uint16_t maskBits,
                                                   int numChannels,
                                                   int measChoice) const noexcept;

        float computeDesiredGainDbForMeasurement (float measuredLufs,
                                                  float targetLufs,
                                                  float maxBoostDb,
                                                  float maxCutDb) const noexcept;

        float applyRateLimitDb (float currentDb,
                                float targetDb,
                                float rateUpDbPerSec,
                                float rateDownDbPerSec,
                                int numSamples) const noexcept;

        void commitLearnedGain() noexcept;
        void updateMeteringFromAppliedGain (int numChannels,
                                            int numSamples) noexcept;
    };
    // [END LVLR-MODULE-CLASS]
} // namespace levelscope
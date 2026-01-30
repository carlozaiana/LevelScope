#pragma once

#include "../IAudioModule.h"
// [BEGIN MTDM-PARAM-HEADER-INCLUDE]
#include "MultiThresholdDynamicsParamIDs.h"
// [END MTDM-PARAM-HEADER-INCLUDE]

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
            void bindParameters (std::atomic<float>* enabled01,
                                 std::atomic<float>* thresholdDb,
                                 std::atomic<float>* ratio) noexcept;

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
    };
    // [END MTDM-MODULE-DECL]
} // namespace levelscope
#pragma once

#include "../IAudioModule.h"

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

        void process (ProcessContext& ctx) noexcept override;

        // Persistence (non-audio-thread only)
        juce::ValueTree saveState() const override;
        void loadState (const juce::ValueTree& state) override;

    private:
        // Stored from prepare() for debugging/future DSP wiring. Not used in process().
        double preparedSampleRate = 0.0;
        int preparedMaxBlockSize = 0;
        juce::AudioChannelSet preparedChannelSet;
    };
    // [END MTDM-MODULE-DECL]
} // namespace levelscope

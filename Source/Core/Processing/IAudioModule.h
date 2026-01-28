#pragma once

#include "ProcessContext.h"

#include <juce_core/juce_core.h>
#include <juce_data_structures/juce_data_structures.h>

#include <atomic>

namespace levelscope
{
    struct ModulePrepareSpec
    {
        double sampleRate = 48000.0;
        int maxBlockSize = 0;
        juce::AudioChannelSet channelSet;

        // Future: hooks for shared scratch buffers, latency reporting, oversampling context, etc.
    };

    // [BEGIN IAudioModule]
    // Contract:
    // - process() must be RT-safe: no locks, no heap allocations, no file I/O, no ValueTree access.
    // - saveState/loadState are non-audio-thread only.
    class IAudioModule
    {
    public:
        virtual ~IAudioModule() = default;

        // Stable identifier used for persistence and module graph reconstruction.
        // Must not change once released.
        virtual juce::String getModuleID() const = 0;

        // UI-facing name (may change without breaking persistence).
        virtual juce::String getDisplayName() const = 0;

        // Lifecycle
        virtual void prepare (const ModulePrepareSpec& spec) = 0;
        virtual void reset() = 0;

        // DSP
        virtual void process (ProcessContext& ctx) noexcept = 0;

        // Bypass (RT-safe)
        virtual void setBypassed (bool shouldBypass) noexcept = 0;
        virtual bool isBypassed() const noexcept = 0;

        // Persistence (non-audio-thread only)
        virtual juce::ValueTree saveState() const = 0;
        virtual void loadState (const juce::ValueTree& state) = 0;
    };
    // [END IAudioModule]

    // Optional convenience base: provides atomic bypass + minimal default state handling.
    class AudioModuleBase : public IAudioModule
    {
    public:
        void setBypassed (bool shouldBypass) noexcept override { bypassed.store (shouldBypass, std::memory_order_relaxed); }
        bool isBypassed() const noexcept override { return bypassed.load (std::memory_order_relaxed); }

        juce::ValueTree saveState() const override
        {
            // Type is module ID to keep tree self-describing.
            juce::ValueTree t (getModuleID());
            t.setProperty ("bypassed", isBypassed(), nullptr);
            return t;
        }

        void loadState (const juce::ValueTree& state) override
        {
            // Be tolerant of missing properties / older versions.
            const auto b = state.getProperty ("bypassed", false);
            setBypassed ((bool) b);
        }

    private:
        std::atomic<bool> bypassed { false };
    };
} // namespace levelscope
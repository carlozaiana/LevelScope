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
    }

    void MultiThresholdDynamicsModule::reset()
    {
        // No state yet.
    }

    void MultiThresholdDynamicsModule::process (ProcessContext& ctx) noexcept
    {
        // Stage B: intentional no-op (pass-through).
        juce::ignoreUnused (ctx);
    }

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
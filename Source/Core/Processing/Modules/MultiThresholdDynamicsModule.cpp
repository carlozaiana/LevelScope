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

    // [BEGIN MTDM-RT-PARAM-BINDING-IMPL]
        void MultiThresholdDynamicsModule::bindParameters (std::atomic<float>* enabled01,
                                                          std::atomic<float>* thresholdDb,
                                                          std::atomic<float>* ratio) noexcept
        {
            // Non-RT call expected (constructor / graph build thread).
            // Storing pointers is RT-safe; we only read them in process().
            pEnabled01   = enabled01;
            pThresholdDb = thresholdDb;
            pRatio       = ratio;
        }

        void MultiThresholdDynamicsModule::process (ProcessContext& ctx) noexcept
        {
            // Stage C1: still no DSP. We only demonstrate RT-safe parameter reads.
            // Default is disabled => pass-through.
            const float enabled = (pEnabled01 != nullptr
                                     ? pEnabled01->load (std::memory_order_relaxed)
                                     : 0.0f);

            if (enabled < 0.5f)
            {
                // Module disabled => guaranteed pass-through.
                juce::ignoreUnused (ctx);
                return;
            }

            // Enabled but still skeleton/no-op.
            // (Reading these avoids “unused member” warnings and validates binding path.)
            if (pThresholdDb != nullptr) (void) pThresholdDb->load (std::memory_order_relaxed);
            if (pRatio       != nullptr) (void) pRatio->load       (std::memory_order_relaxed);

            juce::ignoreUnused (ctx);
        }
    // [END MTDM-RT-PARAM-BINDING-IMPL]

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
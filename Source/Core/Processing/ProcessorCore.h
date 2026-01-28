#pragma once

#include "IAudioModule.h"

#include <atomic>
#include <memory>
#include <vector>

namespace levelscope
{
    struct ModuleGraph final
    {
        uint32_t revision = 0; // incremented by whoever constructs graphs; useful for UI/debug
        std::vector<std::shared_ptr<IAudioModule>> modules;
    };

    // ProcessorCore: lock-free(ish) module-chain host.
    // Audio thread:
    //   - atomic_load() a shared_ptr<const ModuleGraph>
    //   - iterates modules and calls process()
    //
    // Non-audio thread:
    //   - build a new ModuleGraph (shared_ptr)
    //   - ensure modules are prepared if needed
    //   - atomic_store() to activate
    class ProcessorCore
    {
    public:
        ProcessorCore() = default;

        void prepare (const ModulePrepareSpec& spec)
        {
            currentSpec = spec;

            // prepare currently active modules (caller guarantees audio thread is not concurrently processing)
            auto g = getActiveGraph();
            if (! g)
                return;

            for (auto& m : g->modules)
                if (m)
                    m->prepare (spec);
        }

        void reset()
        {
            auto g = getActiveGraph();
            if (! g)
                return;

            for (auto& m : g->modules)
                if (m)
                    m->reset();
        }

        void process (ProcessContext& ctx) noexcept
        {
            auto g = getActiveGraph();
            if (! g)
                return;

            for (auto& m : g->modules)
            {
                if (m && ! m->isBypassed())
                    m->process (ctx);
            }
        }

        // Graph management
        void setActiveGraph (std::shared_ptr<const ModuleGraph> newGraph) noexcept
        {
            // Contract: caller must ensure any required prepare() has already been done safely.
            std::atomic_store_explicit (&activeGraph, std::move (newGraph), std::memory_order_release);
        }

        std::shared_ptr<const ModuleGraph> getActiveGraph() const noexcept
        {
            return std::atomic_load_explicit (&activeGraph, std::memory_order_acquire);
        }

        const ModulePrepareSpec& getCurrentPrepareSpec() const noexcept { return currentSpec; }

    private:
        ModulePrepareSpec currentSpec {};
        mutable std::shared_ptr<const ModuleGraph> activeGraph;
    };
} // namespace levelscope
#pragma once
#include "StandaloneServices.h"
#include <cmath>

namespace levelscope::standalone
{
// The qualified static-gain-v1 operation. No channel-role weighting, clipping,
// monitor controls or plugin processing belongs here.
class StaticGainProcessor
{
public:
    static bool supports (const ProcessingPlan& p) noexcept
    {
        return p.schemaVersion == 1 && p.engineVersion == "static-gain-v1"
            && p.graphJson.isEmpty() && std::isfinite (p.staticGainDb)
            && p.staticGainDb >= -60 && p.staticGainDb <= 24;
    }
    explicit StaticGainProcessor (double db) : coefficient (std::pow (10.0, db / 20)) {}
    bool process (float* const* audio, int channels, int frames) const noexcept
    {
        for (int channel = 0; channel < channels; ++channel)
            for (int i = 0; i < frames; ++i)
            {
                const float x = static_cast<float> (audio[channel][i] * coefficient);
                if (! std::isfinite (x)) return false;
                audio[channel][i] = x;
            }
        return true;
    }
private:
    double coefficient;
};

// Injectable test sink; caller owns storage/lifetime. RT implementations must be
// bounded, allocation-free and noexcept. No observer is installed in the product.
struct ProcessedAudioObserver
{
    virtual ~ProcessedAudioObserver() = default;
    virtual void observe (const float* const*, int channels, int frames, int64_t sourceFrame) noexcept = 0;
};
}

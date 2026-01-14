#pragma once

#include <vector>
#include <cmath>
#include <juce_core/juce_core.h>

//==============================================================================
// [BS1770] K-weighting filter (BS.1770 / EBU R128)
// Header-only utility. Not wired into the plugin yet (Phase 2A part 2).
//
// Implementation uses two biquads per channel:
//   1) High-shelf (approx +4 dB above ~1.7 kHz)
//   2) High-pass  (low frequency rolloff)
//
// Coefficients are computed per sample rate using RBJ cookbook formulas.
//==============================================================================

class BS1770KWeighting
{
public:
    BS1770KWeighting() = default;

    void prepare (double sampleRate, int numChannels)
    {
        fs = (sampleRate > 0.0 ? sampleRate : 48000.0);
        channels = juce::jmax (1, numChannels);

        shelf.designHighShelf (fs, 1681.974450955533, 0.7071752369554196, 3.999843853973347);
        highPass.designHighPass (fs, 38.13547087602444, 0.5003270373238773);

        shelfState.assign ((size_t) channels, {});
        hpState.assign ((size_t) channels, {});
        reset();
    }

    void reset()
    {
        for (auto& s : shelfState) s.reset();
        for (auto& s : hpState)    s.reset();
    }

    float processSample (int ch, float x) noexcept
    {
        if (ch < 0 || ch >= channels)
            return x;

        const double y1 = shelf.process (shelfState[(size_t) ch], (double) x);
        const double y2 = highPass.process (hpState[(size_t) ch], y1);
        return (float) y2;
    }

private:
    struct BiquadState
    {
        double z1 = 0.0, z2 = 0.0;
        void reset() noexcept { z1 = 0.0; z2 = 0.0; }
    };

    struct Biquad
    {
        // Normalised (a0 = 1)
        double b0 = 1.0, b1 = 0.0, b2 = 0.0;
        double a1 = 0.0, a2 = 0.0;

        double process (BiquadState& st, double x) const noexcept
        {
            // DF2T
            const double y = b0 * x + st.z1;
            st.z1 = b1 * x - a1 * y + st.z2;
            st.z2 = b2 * x - a2 * y;
            return y;
        }

        void designHighPass (double sampleRate, double f0, double Q)
        {
            const double w0 = 2.0 * juce::MathConstants<double>::pi * f0 / sampleRate;
            const double cosw0 = std::cos (w0);
            const double sinw0 = std::sin (w0);
            const double alpha = sinw0 / (2.0 * Q);

            const double a0 = 1.0 + alpha;

            const double nb0 =  (1.0 + cosw0) * 0.5;
            const double nb1 = -(1.0 + cosw0);
            const double nb2 =  (1.0 + cosw0) * 0.5;

            const double na1 = -2.0 * cosw0;
            const double na2 =  1.0 - alpha;

            b0 = nb0 / a0;
            b1 = nb1 / a0;
            b2 = nb2 / a0;
            a1 = na1 / a0;
            a2 = na2 / a0;
        }

        void designHighShelf (double sampleRate, double f0, double Q, double gainDb)
        {
            const double A = std::pow (10.0, gainDb / 40.0);
            const double w0 = 2.0 * juce::MathConstants<double>::pi * f0 / sampleRate;
            const double cosw0 = std::cos (w0);
            const double sinw0 = std::sin (w0);
            const double alpha = sinw0 / (2.0 * Q);
            const double sqrtA = std::sqrt (A);

            const double a0 =        (A + 1.0) - (A - 1.0) * cosw0 + 2.0 * sqrtA * alpha;
            const double na1 =  2.0 * ((A - 1.0) - (A + 1.0) * cosw0);
            const double na2 =        (A + 1.0) - (A - 1.0) * cosw0 - 2.0 * sqrtA * alpha;

            const double nb0 =    A * ((A + 1.0) + (A - 1.0) * cosw0 + 2.0 * sqrtA * alpha);
            const double nb1 = -2.0 * A * ((A - 1.0) + (A + 1.0) * cosw0);
            const double nb2 =    A * ((A + 1.0) + (A - 1.0) * cosw0 - 2.0 * sqrtA * alpha);

            b0 = nb0 / a0;
            b1 = nb1 / a0;
            b2 = nb2 / a0;
            a1 = na1 / a0;
            a2 = na2 / a0;
        }
    };

    double fs = 48000.0;
    int channels = 2;

    Biquad shelf;
    Biquad highPass;

    std::vector<BiquadState> shelfState;
    std::vector<BiquadState> hpState;
};
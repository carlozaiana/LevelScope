#pragma once

// Shared upward gain-law functions (header-only, no JUCE dependency, no allocations).
// Intended to be used by DSP blocks and UI curve sampling.
//
// Semantics:
// - levelDb is a generic "level" axis in dB (often LUFS-ish proxy).
// - Returns an UPWARD boost mapping (gain >= 1).
// - Active zone: fades in approaching T0 (below T0), fades out approaching T1 (below T1).
// - Above/at T1 the *target* gain is unity (0 dB). Any smoothing is handled by the caller.

#include <algorithm> // std::clamp, std::min, std::max
#include <cmath>     // std::pow, std::log10

namespace levelscope::dsp::UpwardGainLaw
{
    enum class CurveType : int
    {
        monotonic = 0, // quiet end (near T0) => max boost, approaching T1 => 0 boost
        bell      = 1  // mid-zone boosted most
    };

    struct Params
    {
        float t0Db = -45.0f;
        float t1Db = -30.0f;

        float lowKneeDb  = 3.0f;   // fade-in width below T0
        float highKneeDb = 3.0f;   // fade-out width below T1

        float maxBoostDb = 8.0f;   // cap
        float curve01    = 0.5f;   // expo = 1 + curve01*3

        CurveType curveType = CurveType::monotonic;
    };

    inline float clamp01 (float x) noexcept
    {
        return std::clamp (x, 0.0f, 1.0f);
    }

    inline float dbToLin (float db) noexcept
    {
        return std::pow (10.0f, db / 20.0f);
    }

    inline float linToDb (float g) noexcept
    {
        g = std::max (1.0e-12f, g);
        return (float) (20.0 * std::log10 ((double) g));
    }

    inline float smoothstep01 (float t) noexcept
    {
        t = clamp01 (t);
        return t * t * (3.0f - 2.0f * t);
    }

    // One-sided knee: ramps from 0..1 over [threshold - kneeWidth, threshold]
    inline float kneeUpToThreshold01 (float levelDb, float thresholdDb, float kneeWidthDb) noexcept
    {
        kneeWidthDb = std::max (1.0e-4f, kneeWidthDb);
        const float start = thresholdDb - kneeWidthDb;

        if (levelDb <= start)      return 0.0f;
        if (levelDb >= thresholdDb) return 1.0f;

        const float t = (levelDb - start) / kneeWidthDb;
        return smoothstep01 (t);
    }

    // 0..1 "active zone" factor using one-sided knees
    inline float computeActiveZone01 (float levelDb, const Params& p) noexcept
    {
        const float t0 = std::min (p.t0Db, p.t1Db);
        const float t1 = std::max (p.t0Db, p.t1Db);

        const float inAroundT0  = kneeUpToThreshold01 (levelDb, t0, p.lowKneeDb);
        const float outAroundT1 = 1.0f - kneeUpToThreshold01 (levelDb, t1, p.highKneeDb);

        return clamp01 (inAroundT0 * outAroundT1);
    }

    // Returns boost in dB (>= 0), already including active-zone knees and maxBoost clamp.
    inline float computeBaseBoostDb (float levelDb, const Params& p) noexcept
    {
        const float t0 = std::min (p.t0Db, p.t1Db);
        const float t1 = std::max (p.t0Db, p.t1Db);

        const float maxBoostDb = std::max (0.0f, p.maxBoostDb);
        if (maxBoostDb <= 0.0f)
            return 0.0f;

        const float zone01 = computeActiveZone01 (levelDb, p);
        if (zone01 <= 1.0e-6f)
            return 0.0f;

        const float range = std::max (1.0f, t1 - t0);
        float pos = (levelDb - t0) / range; // 0..1 across zone
        pos = std::clamp (pos, 0.0f, 1.0f);

        const float curve01 = clamp01 (p.curve01);
        const float expo = 1.0f + curve01 * 3.0f;

        float shaped = 0.0f;

        if (p.curveType == CurveType::monotonic)
        {
            shaped = std::pow (std::max (0.0f, 1.0f - pos), expo);
        }
        else // bell
        {
            const float d = std::abs (pos - 0.5f) * 2.0f; // 0..1
            shaped = std::pow (std::max (0.0f, 1.0f - d), expo);
        }

        float boostDb = maxBoostDb * shaped * zone01;
        boostDb = std::clamp (boostDb, 0.0f, maxBoostDb);
        return boostDb;
    }

    // Base gain (>= 1.0), no amount fade.
    inline float computeBaseGainLin (float levelDb, const Params& p) noexcept
    {
        return dbToLin (computeBaseBoostDb (levelDb, p));
    }

    // Amount fade-in/out in linear domain (matches existing DSP style).
    inline float applyAmountLin (float baseGainLin, float amount01) noexcept
    {
        amount01 = clamp01 (amount01);
        baseGainLin = std::max (1.0f, baseGainLin);
        return 1.0f + (baseGainLin - 1.0f) * amount01;
    }

    inline float computeGainLin (float levelDb, const Params& p, float amount01) noexcept
    {
        return applyAmountLin (computeBaseGainLin (levelDb, p), amount01);
    }

    inline float computeGainDb (float levelDb, const Params& p, float amount01) noexcept
    {
        return linToDb (computeGainLin (levelDb, p, amount01));
    }

    // UI-friendly alias: "gain in dB" meaning upward boost dB (>= 0)
    inline float computeUpwardGainDb (float levelDb, const Params& p, float amount01) noexcept
    {
        const float gDb = computeGainDb (levelDb, p, amount01);
        return (gDb > 0.0f ? gDb : 0.0f);
    }
} // namespace levelscope::dsp::UpwardGainLaw
#include "AttributionEngine.h"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace levelscope::proposal
{
bool buildAttribution (ProgramAnalysis& a, const std::function<bool()>& cancelled)
{
    double total = 0;
    for (size_t n = 0; n < a.momentary.size(); ++n)
    {
        if ((n % 256) == 0 && cancelled()) return false;
        const auto& w = a.momentary[n];
        if (w.included) { total += w.energy; a.integratedRank.push_back (n); }
        else a.excludedIntegrated.push_back (n);
    }
    for (auto n : a.integratedRank) a.momentary[n].energyShare = a.momentary[n].energy / total;
    // Stable explicit tie break; all full-resolution entries are retained.
    auto order = [&] (const auto& windows, bool descending)
    {
        return [&windows, descending] (size_t x, size_t y)
        {
            if (windows[x].energy != windows[y].energy)
                return descending ? windows[x].energy > windows[y].energy : windows[x].energy < windows[y].energy;
            return windows[x].end < windows[y].end;
        };
    };
    std::sort (a.integratedRank.begin(), a.integratedRank.end(), order (a.momentary, true));
    if (cancelled()) return false;
    for (size_t n = 0; n < a.lraWindows.size(); ++n)
    {
        if ((n % 256) == 0 && cancelled()) return false;
        const auto& w = a.lraWindows[n];
        if (w.included) a.lraOrder.push_back (n);
        if (! w.lufs) ++a.missingLraWindows;
    }
    std::sort (a.lraOrder.begin(), a.lraOrder.end(), order (a.lraWindows, false));
    if (cancelled()) return false;
    for (auto n : a.lraOrder)
    {
        const auto v = *a.lraWindows[n].lufs;
        if (a.p10 && v <= *a.p10) a.lowerTail.push_back (n);
        if (a.p95 && v >= *a.p95) a.upperTail.push_back (n);
    }
    std::sort (a.upperTail.begin(), a.upperTail.end(), order (a.lraWindows, true));
    double low = 0, high = 0; bool any = false;
    for (const auto& w : a.lraWindows) if (w.lufs)
    {
        if (! any) low = high = *w.lufs;
        low = std::min (low, *w.lufs); high = std::max (high, *w.lufs); any = true;
    }
    if (any)
    {
        low = std::floor (low); high = std::floor (high);
        // Finite input energies bound this span (~6500 dB for binary64).
        for (double v = low; v <= high; ++v) a.histogram.push_back ({ v, 0, 0, 0 });
        for (const auto& w : a.lraWindows) if (w.lufs)
        {
            auto& bin = a.histogram[static_cast<size_t> (std::floor (*w.lufs) - low)];
            if (w.included) ++bin.included; else ++bin.excluded;
        }
        size_t count = 0;
        for (auto& b : a.histogram) { count += b.included; b.cumulativeIncluded = count; }
    }
    return ! cancelled();
}
}

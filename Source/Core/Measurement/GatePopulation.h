#pragma once
#include "OfflineMeasurement.h"
#include <cmath>

namespace levelscope::measurement
{
// Shared full-precision population selection. These are the existing M1
// expressions/order, extracted for descriptive analysis; never gate displayed LUFS.
inline double absoluteGateEnergy() { return std::pow (10.0, (-70 + 0.691) / 10); }
inline bool aboveAbsolute (double energy, bool lra)
{ return lra ? energy >= absoluteGateEnergy() : energy > absoluteGateEnergy(); }
inline bool includedInPopulation (double energy, double relative, bool lra)
{ return aboveAbsolute (energy, lra) && (lra ? energy >= relative : energy > relative); }
template <typename Cancel>
std::optional<double> relativeGateEnergy (const std::vector<Window>& windows, bool lra, Cancel cancelled)
{
    double sum = 0; size_t count = 0, index = 0;
    for (const auto& w : windows)
    {
        if ((index++ % 256) == 0 && cancelled()) return {};
        if (aboveAbsolute (w.energy, lra)) { sum += w.energy; ++count; }
    }
    if (count == 0) return {};
    return sum / static_cast<double> (count) / (lra ? 100 : 10);
}
inline std::optional<double> relativeGateEnergy (const std::vector<Window>& windows, bool lra)
{ return relativeGateEnergy (windows, lra, [] { return false; }); }
inline size_t percentileIndex (size_t count, double p)
{ return static_cast<size_t> (std::floor ((count - 1) * p + 0.5)); }
}

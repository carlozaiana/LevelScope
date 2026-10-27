#pragma once
#include <juce_core/juce_core.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <locale>
#include <sstream>
#include <vector>

namespace levelscope::processing
{
// Whole-programme processing, after composition and before the final global trim.
// Original components and annotation documents are deliberately not owned here.
struct GainInterval
{
    int64_t begin = 0, end = 0;
    double db = 0;
};
struct SectionGainPlan
{
    static constexpr const char* version = "section-gain-v1";
    int schema = 1, sampleRate = 0;
    int64_t frames = 0;
    juce::String inputIdentity;
    std::vector<GainInterval> intervals;

    // The canonical representation excludes annotation IDs/names and redundant
    // boundaries. The ramp rule and rounding order are part of version v1.
    SectionGainPlan canonical() const
    {
        auto p = *this; p.intervals.clear();
        for (auto row : intervals)
        {
            if (row.db == 0) row.db = 0; // signed zero is not a different plan
            if (!p.intervals.empty() && p.intervals.back().end == row.begin && p.intervals.back().db == row.db)
                p.intervals.back().end = row.end;
            else p.intervals.push_back(row);
        }
        return p;
    }
    int64_t radius(size_t right) const noexcept
    {
        return std::min({int64_t(std::llround(sampleRate * .050)),
            (intervals[right-1].end-intervals[right-1].begin)/2,
            (intervals[right].end-intervals[right].begin)/2});
    }
    juce::String validationError() const
    {
        if (schema != 1 || (sampleRate != 44100 && sampleRate != 48000 && sampleRate != 96000 && sampleRate != 192000)
            || frames <= 0 || frames > int64_t(sampleRate)*10800 || inputIdentity.length()!=64
            || !inputIdentity.containsOnly("0123456789abcdef") || intervals.empty() || intervals.size()>10000)
            return "Invalid section-gain version, input, rate, duration or interval count";
        int64_t end=0;
        for (const auto& row : intervals)
        {
            if (row.begin!=end || row.end<=row.begin || row.end>frames || !std::isfinite(row.db) || row.db < -12 || row.db > 12)
                return "Section gains require complete contiguous frame coverage and finite gains from -12 to +12 dB";
            end=row.end;
        }
        if (end!=frames) return "Section gains must cover the complete programme";
        const auto p=canonical();
        for (size_t n=1;n<p.intervals.size();++n)
            if (p.radius(n)==0) return "Differing gains need at least two frames on each side of a transition";
        return {};
    }
    bool isUnity() const noexcept
    { return std::all_of(intervals.begin(),intervals.end(),[](const auto& i){return i.db==0;}); }
    juce::String audioKey() const
    {
        const auto p=canonical();
        std::ostringstream out; out.imbue(std::locale::classic()); out << std::hexfloat;
        out << version << ";after-composition;before-static-gain-v1;float32-stage-rounding;" << schema << ';'
            << sampleRate << ';' << frames << ';' << inputIdentity << ';';
        for(const auto& row:p.intervals) out << row.begin << ':' << row.end << ':' << row.db << ';';
        return out.str();
    }
    // Call on a validated canonical plan. No mutable cursor: seek/chunk invariant.
    double gainAt(int64_t frame) const noexcept
    {
        const auto it=std::upper_bound(intervals.begin(),intervals.end(),frame,
            [](int64_t f,const GainInterval& i){return f<i.end;});
        const size_t n=static_cast<size_t>(it-intervals.begin());
        if(n>=intervals.size()) return 0; // callers reject out-of-range requests
        auto ramp=[&](size_t right)
        {
            const auto r=radius(right), begin=intervals[right].begin-r;
            if(frame==begin) return intervals[right-1].db;
            if(frame==begin+2*r-1) return intervals[right].db;
            const double t=double(frame-begin)/double(2*r-1);
            return intervals[right-1].db+(intervals[right].db-intervals[right-1].db)*t;
        };
        if(n>0 && frame<intervals[n].begin+radius(n)) return ramp(n);
        if(n+1<intervals.size() && frame>=intervals[n].end-radius(n+1)) return ramp(n+1);
        return intervals[n].db;
    }
};
}

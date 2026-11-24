#include "LraTimelineDiagnostics.h"
#include "Core/Measurement/GatePopulation.h"
#include <algorithm>
#include <stdexcept>

namespace levelscope::proposal::lra_timeline
{
namespace
{
void addLevel(Column& column, unsigned flags, std::optional<double> level)
{
    column.flags |= flags;
    if (!level) return;
    if (!column.minLufs || *level < *column.minLufs) column.minLufs = level;
    if (!column.maxLufs || *level > *column.maxLufs) column.maxLufs = level;
    if ((flags & eligible) && (!column.maxEligibleLufs || *level > *column.maxEligibleLufs))
        column.maxEligibleLufs = level;
}
void combine(Column& dst, const Column& src)
{
    dst.flags |= src.flags;
    if (src.minLufs && (!dst.minLufs || *src.minLufs < *dst.minLufs)) dst.minLufs = src.minLufs;
    if (src.maxLufs && (!dst.maxLufs || *src.maxLufs > *dst.maxLufs)) dst.maxLufs = src.maxLufs;
    if (src.maxEligibleLufs && (!dst.maxEligibleLufs || *src.maxEligibleLufs > *dst.maxEligibleLufs))
        dst.maxEligibleLufs = src.maxEligibleLufs;
}
bool exceedsThreshold(const Target& target, const juce::String& identity, std::optional<double> lufs)
{
    return target.state == TargetState::exceeds && target.snapshotIdentity == identity &&
           lufs && target.upperLevelLufs && *lufs > *target.upperLevelLufs;
}
}

std::shared_ptr<const Index> Index::build(std::shared_ptr<const ProgramAnalysis> measured,
                                          const std::function<bool()>& cancelled)
{
    if (!measured) throw std::invalid_argument("LRA timeline requires a validated measurement");
    if (cancelled()) return {};
    auto index = std::shared_ptr<Index>(new Index(std::move(measured)));
    const auto& a = index->measured();
    if (a.lraStable != (a.lra.value.has_value() && a.frames >= int64_t(a.sampleRate)*60))
        throw std::invalid_argument("LRA stability differs from measured duration and population");
    const auto shortProgramme = a.frames < int64_t(a.sampleRate) * 3;
    const auto n = a.lraWindows.size();
    index->flags.resize(n);
    // Mark only the tied windows whose rank actually follows the P95 observation.
    // M4 lraOrder sorts by exact energy, then earlier window end for equal energy.
    if (!a.lraOrder.empty())
    {
        const auto first = measurement::percentileIndex(a.lraOrder.size(), .95) + 1;
        for (size_t rank = first; rank < a.lraOrder.size(); ++rank)
        {
            if (rank % 256 == 0 && cancelled()) return {};
            const auto i = a.lraOrder[rank];
            if (a.lraWindows[i].lufs == a.p95) index->flags[i] |= upperRankTie;
        }
    }
    index->leafCount = 1;
    while (index->leafCount * blockSize < n) index->leafCount *= 2;
    index->tree.resize(index->leafCount * 2);
    for (size_t i = 0; i < n; ++i)
    {
        if (i % 256 == 0 && cancelled()) return {};
        const auto& w = a.lraWindows[i];
        auto& mask = index->flags[i];
        if (w.padding) mask |= padded;
        if (w.energy == 0) mask |= zeroEnergy;
        if (shortProgramme) mask |= insufficientDuration;
        else if (w.included)
        {
            mask |= eligible;
            if (*w.lufs > *a.p95) mask |= aboveP95;
            if (*w.lufs < *a.p10) mask |= belowP10;
            if (*w.lufs == *a.p10) mask |= atP10;
            if (*w.lufs == *a.p95) mask |= atP95;
        }
        else mask |= measurement::aboveAbsolute(w.energy,true) ? relativeRejected : absoluteRejected;
        addLevel(index->tree[index->leafCount + i/blockSize], mask, w.lufs);
    }
    for (size_t i = index->leafCount - 1; i > 0; --i)
    {
        if (i % 256 == 0 && cancelled()) return {};
        auto& node = index->tree[i];
        node = index->tree[2*i];
        combine(node,index->tree[2*i+1]);
    }
    return cancelled() ? nullptr : index;
}

size_t Index::upperRankCount() const
{
    const auto count = analysis->lraOrder.size();
    return count ? count - 1 - measurement::percentileIndex(count, .95) : 0;
}
Target Index::target(const TargetProfile& profile, Applicability applicability) const
{
    Target t;
    t.snapshotIdentity = analysis->snapshotIdentity;
    if (profile.jobRequirementsError.isNotEmpty()) { t.state = TargetState::inapplicable; return t; }
    if (!profile.limits.hasMaxLoudnessRange) return t;
    if (profile.validationError().isNotEmpty() || applicability != profile.applicability)
    { t.state = TargetState::inapplicable; return t; }
    t.ceilingLu = profile.limits.maxLoudnessRangeLu;
    if (!analysis->lra.value || !analysis->p10 || !analysis->p95) t.state = TargetState::unavailable;
    else if (!analysis->lraStable) t.state = TargetState::unstable;
    else
    {
        t.upperLevelLufs = *analysis->p10 + *t.ceilingLu;
        t.state = *analysis->lra.value > *t.ceilingLu ? TargetState::exceeds : TargetState::within;
    }
    return t;
}

std::pair<size_t,size_t> Index::windowRange(int64_t begin, int64_t end) const
{
    // The qualified windows begin at i*hop and end exactly 30 hops later.
    // Intersect two half-open intervals. Callers already clip the query to real audio.
    const int64_t hop = analysis->sampleRate / 10;
    const auto count = flags.size();
    const auto first = static_cast<size_t>(std::max<int64_t>(0, begin/hop - 29));
    const auto last = static_cast<size_t>((end-1)/hop + 1);
    return {std::min(first,count),std::min(last,count)};
}

Column Index::aggregate(size_t first, size_t last) const
{
    Column column;
    const auto add = [&](size_t i) { addLevel(column, flags[i], analysis->lraWindows[i].lufs); };
    while (first < last && first % blockSize) add(first++);
    while (first < last && last % blockSize) add(--last);
    size_t l = leafCount + first/blockSize, r = leafCount + last/blockSize;
    while (l < r)
    {
        if (l % 2) combine(column,tree[l++]);
        if (r % 2) combine(column,tree[--r]);
        l /= 2;
        r /= 2;
    }
    return column;
}

std::vector<Column> Index::columns(int64_t begin, int64_t end, int pixels, const Target& t) const
{
    begin = std::max<int64_t>(0,begin);
    end = std::min(analysis->frames,end);
    if (end <= begin || pixels <= 0) return {};
    pixels = std::min(4096,pixels);
    std::vector<Column> out(static_cast<size_t>(pixels));
    const auto width = end - begin;
    for (int x = 0; x < pixels; ++x)
    {
        const auto first = begin + width*x/pixels;
        const auto last = begin + (width*(x+1)+pixels-1)/pixels;
        const auto range = windowRange(first,last);
        auto& column = out[static_cast<size_t>(x)];
        column = aggregate(range.first,range.second);
        if (exceedsThreshold(t,analysis->snapshotIdentity,column.maxEligibleLufs)) column.flags |= targetPassage;
    }
    return out;
}

std::vector<Window> Index::inspect(int64_t frame, const Target& t) const
{
    if (frame < 0 || frame >= analysis->frames) return {};
    const auto [first,last] = windowRange(frame,frame+1);
    std::vector<Window> out;
    out.reserve(last-first);
    for (auto i = first; i < last; ++i)
    {
        const auto& w = analysis->lraWindows[i];
        unsigned mask = flags[i];
        if ((mask & eligible) && exceedsThreshold(t,analysis->snapshotIdentity,w.lufs)) mask |= targetPassage;
        out.push_back({i,w.begin,std::min(w.end,analysis->frames),w.end,w.lufs,mask});
    }
    return out;
}
size_t Index::indexBytes() const
{
    return sizeof(*this) + flags.capacity()*sizeof(flags.front()) + tree.capacity()*sizeof(Column);
}
}

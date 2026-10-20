#include "SectionAnalysis.h"
#include <juce_cryptography/juce_cryptography.h>
#include <set>
#include <limits>
#include <sstream>
#include <locale>
#include <cmath>

namespace levelscope::proposal
{
namespace
{
void require (bool yes, const char* message) { if (! yes) throw std::invalid_argument (message); }
juce::String hash (const std::string& bytes) { return juce::SHA256 (bytes.data(), bytes.size()).toHexString(); }
void field (std::ostringstream& out, const juce::String& s) { const auto b = s.toStdString(); out << b.size() << ':' << b << ';'; }
}
bool SectionMap::valid() const
{
    if (schemaVersion != 1 || inputIdentity.isEmpty() || (sampleRate != 44100 && sampleRate != 48000 && sampleRate != 96000 && sampleRate != 192000)
        || frames <= 0 || frames > static_cast<int64_t> (sampleRate) * 86400 || sections.empty() || sections.size() > 10000) return false;
    int64_t end = 0; std::set<uint64_t> ids;
    for (const auto& s : sections)
    {
        if (s.id == 0 || s.id >= nextId || ! ids.insert (s.id).second || s.begin != end || s.end <= s.begin || s.end > frames
            || s.name.trim().isEmpty() || s.name.length() > 128) return false;
        end = s.end;
    }
    return end == frames;
}
juce::String SectionMap::digest() const
{
    require (valid(), "Invalid section map");
    std::ostringstream s; s.imbue (std::locale::classic()); field (s, sectionsVersion); field (s, inputIdentity);
    s << schemaVersion << ';' << sampleRate << ';' << frames << ';' << revision << ';' << nextId << ';';
    for (const auto& row : sections) { s << row.id << ';' << row.begin << ';' << row.end << ';'; field (s, row.name); }
    return hash (s.str());
}
bool SectionDocument::bind (const juce::String& input, int rate, int64_t frames)
{
    if (input == value.inputIdentity && rate == value.sampleRate && frames == value.frames) return false;
    SectionMap next;
    if (input.isNotEmpty())
    {
        next.inputIdentity = input; next.sampleRate = rate; next.frames = frames; next.nextId = 2;
        next.sections.push_back ({ 1, 0, frames, "Whole programme" });
        require (next.valid(), "Unsupported programme section input");
    }
    value = std::move (next); return true;
}
bool SectionDocument::edit (uint64_t id, const std::function<void(SectionMap&, size_t)>& fn, juce::String& error)
{
    error.clear();
    try
    {
        require (value.valid(), "Confirm a valid source and channel interpretation first");
        auto it = std::find_if (value.sections.begin(), value.sections.end(), [id] (const auto& s) { return s.id == id; });
        require (it != value.sections.end(), "Select an existing section");
        auto next = value; fn (next, static_cast<size_t> (it - value.sections.begin()));
        require (next.valid(), "Sections must form non-empty contiguous source-frame intervals (maximum 10,000 sections)");
        if (next.digest() == value.digest()) return true;
        require (value.revision < std::numeric_limits<uint64_t>::max(), "Section revision exhausted");
        next.revision = value.revision + 1; value = std::move (next); return true;
    }
    catch (const std::exception& e) { error = e.what(); return false; }
}
bool SectionDocument::split (uint64_t id, int64_t frame, juce::String& error)
{
    return edit (id, [frame] (auto& map, size_t n)
    {
        auto& row = map.sections[n]; require (frame > row.begin && frame < row.end, "Split frame must lie strictly inside the selected section");
        require (map.nextId < std::numeric_limits<uint64_t>::max(), "Section IDs exhausted");
        const auto newId = map.nextId++;
        ProgrammeSection right { newId, frame, row.end, "Section " + juce::String (newId) };
        row.end = frame; map.sections.insert (map.sections.begin() + static_cast<ptrdiff_t> (n + 1), std::move (right));
    }, error);
}
bool SectionDocument::moveEnd (uint64_t id, int64_t frame, juce::String& error)
{
    return edit (id, [frame] (auto& map, size_t n)
    {
        require (n + 1 < map.sections.size(), "Programme end is fixed; select an internal boundary");
        require (frame > map.sections[n].begin && frame < map.sections[n + 1].end, "Boundary must leave both neighbouring sections non-empty");
        map.sections[n].end = frame; map.sections[n + 1].begin = frame;
    }, error);
}
bool SectionDocument::mergeNext (uint64_t id, juce::String& error)
{
    return edit (id, [] (auto& map, size_t n)
    {
        require (n + 1 < map.sections.size(), "No following section to merge");
        map.sections[n].end = map.sections[n + 1].end; map.sections.erase (map.sections.begin() + static_cast<ptrdiff_t> (n + 1));
    }, error);
}
bool SectionDocument::rename (uint64_t id, const juce::String& name, juce::String& error)
{ return edit (id, [&name] (auto& map, size_t n) { map.sections[n].name = name.trim(); }, error); }
std::optional<size_t> WindowSubset::index (size_t row) const
{
    if (row < first.size()) return first.begin + row;
    row -= first.size(); if (row < second.size()) return second.begin + row;
    return {};
}
size_t SectionAnalysis::storageBytes() const
{
    size_t bytes = sizeof (*this) + sections.capacity() * sizeof (SectionSummary);
    for (const auto& s : sections) bytes += s.section.name.getNumBytesAsUTF8() + s.programmeMaxima.capacity() * sizeof (SectionPeakLocation);
    return bytes + identity.getNumBytesAsUTF8() + mapIdentity.getNumBytesAsUTF8() + analysisIdentity.getNumBytesAsUTF8();
}
std::shared_ptr<const SectionAnalysis> analyzeSections (const ProgramAnalysis& a, const SectionMap& map,
    const juce::String& expected, const std::function<bool()>& cancelled)
{
    if (cancelled()) return {};
    require (map.valid() && map.inputIdentity == expected && a.sampleRate == map.sampleRate && a.frames == map.frames
             && a.identity.isNotEmpty() && a.snapshotIdentity.isNotEmpty(), "Section map and measured programme must describe matching inputs");
    auto result = std::make_shared<SectionAnalysis>(); result->mapIdentity = map.digest(); result->analysisIdentity = a.identity;
    result->identity = hash ((juce::String (sectionsVersion) + ";" + result->mapIdentity + ";" + a.identity).toStdString());
    result->sections.reserve (map.sections.size());
    auto summarize = [&] (const auto& windows, const ProgrammeSection& section, bool gated, bool lra, SectionWindowSummary& out)
    {
        const auto offset = [&] (auto it) { return static_cast<size_t> (it - windows.begin()); };
        const auto first = offset (std::upper_bound (windows.begin(), windows.end(), section.begin, [] (auto frame, const auto& w) { return frame < w.end; }));
        const auto last = offset (std::lower_bound (windows.begin(), windows.end(), section.end, [] (const auto& w, auto frame) { return w.begin < frame; }));
        const auto noPadding = std::clamp (offset (std::upper_bound (windows.begin(), windows.end(), a.frames, [] (auto frame, const auto& w) { return frame < w.end; })), first, last);
        auto begin = offset (std::lower_bound (windows.begin(), windows.end(), section.begin, [] (const auto& w, auto frame) { return w.begin < frame; }));
        auto end = offset (std::upper_bound (windows.begin(), windows.end(), section.end, [] (auto frame, const auto& w) { return frame < w.end; }));
        begin = std::clamp (begin, first, noPadding); end = std::clamp (end, first, noPadding);
        if (end < begin) begin = end = first;
        out.contained.first = { begin, end }; out.crossing = { { first, begin }, { end, noPadding } }; out.padding.first = { noPadding, last };
        if (gated && ! lra && ! a.integratedRank.empty()) out.programmeBlockEnergyShare = 0;
        std::vector<double> finite; finite.reserve (end - begin);
        for (size_t n = begin; n < end; ++n)
        {
            if ((n % 256) == 0 && cancelled()) return false;
            const auto& w = windows[n];
            if (w.lufs) { require (std::isfinite (*w.lufs), "Non-finite measured section window"); finite.push_back (*w.lufs); }
            else ++out.missingCount;
            if (gated)
            {
                if (w.included) { ++out.included; if (out.programmeBlockEnergyShare) *out.programmeBlockEnergyShare += w.energyShare; }
                else ++out.excluded;
            }
            if (lra && w.included && w.lufs)
            { if (a.p10 && *w.lufs <= *a.p10) ++out.lowerTail; if (a.p95 && *w.lufs >= *a.p95) ++out.upperTail; }
        }
        out.finiteCount = finite.size();
        if (! finite.empty())
        {
            std::sort (finite.begin(), finite.end()); out.minimum = finite.front(); out.maximum = finite.back();
            const auto middle = finite.size() / 2;
            out.median = finite.size() % 2 ? finite[middle] : finite[middle - 1] + (finite[middle] - finite[middle - 1]) / 2;
        }
        return ! cancelled();
    };
    for (const auto& section : map.sections)
    {
        if (cancelled()) return {};
        SectionSummary row; row.section = section;
        if (! summarize (a.momentary, section, true, false, row.momentary)
            || ! summarize (a.shortTerm, section, false, false, row.shortTerm)
            || ! summarize (a.lraWindows, section, true, true, row.lraWindows)) return {};
        for (size_t n = 0; n < a.peaks.size(); ++n)
        {
            const auto& p = a.peaks[n];
            const bool sp = p.sampleDbfs && p.sampleFrame >= section.begin && p.sampleFrame < section.end;
            const bool tp = p.trueDbtp && p.truePeakEighthFrame >= section.begin * 8 && p.truePeakEighthFrame < section.end * 8;
            if (sp || tp) row.programmeMaxima.push_back ({ n, sp, tp });
        }
        result->sections.push_back (std::move (row));
    }
    return cancelled() ? nullptr : result;
}
}

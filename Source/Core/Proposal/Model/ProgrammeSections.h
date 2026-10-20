#pragma once
#include "ProgramAnalysis.h"
#include <functional>

namespace levelscope::proposal
{
inline constexpr const char* sectionsVersion = "programme-sections-m5.1";
struct ProgrammeSection
{
    uint64_t id = 0;
    int64_t begin = 0, end = 0;
    juce::String name;
};
struct SectionMap
{
    int schemaVersion = 1, sampleRate = 0;
    int64_t frames = 0;
    juce::String inputIdentity;
    uint64_t revision = 0, nextId = 1;
    std::vector<ProgrammeSection> sections;
    bool valid() const;
    juce::String digest() const;
};

// Independent annotation document. No audio session, processing or job reference.
// Edits are transactional; rejected and no-op edits preserve identity/revision.
class SectionDocument
{
public:
    bool bind (const juce::String& inputIdentity, int sampleRate, int64_t frames);
    const SectionMap& map() const { return value; }
    bool split (uint64_t id, int64_t frame, juce::String& error);
    bool moveEnd (uint64_t id, int64_t frame, juce::String& error);
    bool mergeNext (uint64_t id, juce::String& error);
    bool rename (uint64_t id, const juce::String&, juce::String& error);
private:
    SectionMap value;
    bool edit (uint64_t, const std::function<void(SectionMap&, size_t)>&, juce::String&);
};

struct IndexRange { size_t begin = 0, end = 0; size_t size() const { return end - begin; } };
struct WindowSubset
{
    IndexRange first, second;
    size_t size() const { return first.size() + second.size(); }
    std::optional<size_t> index (size_t row) const;
};
struct SectionWindowSummary
{
    WindowSubset contained, crossing, padding;
    size_t finiteCount = 0, missingCount = 0, included = 0, excluded = 0, lowerTail = 0, upperTail = 0;
    std::optional<double> minimum, median, maximum;
    std::optional<double> programmeBlockEnergyShare; // unavailable without a programme denominator
};
struct SectionPeakLocation { size_t peakIndex = 0; bool sample = false, truePeak = false; };
struct SectionSummary
{
    ProgrammeSection section;
    SectionWindowSummary momentary, shortTerm, lraWindows;
    std::vector<SectionPeakLocation> programmeMaxima;
};
struct SectionAnalysis
{
    juce::String identity, mapIdentity, analysisIdentity;
    std::vector<SectionSummary> sections;
    size_t storageBytes() const;
};
}

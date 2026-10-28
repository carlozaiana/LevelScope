#pragma once
#include "Model/ProgrammeSections.h"
#include <array>

namespace levelscope::proposal
{
inline constexpr const char* sectionPeakVersion = "continuous-section-peaks-v1";
juce::String sectionGeometryKey(const SectionMap&);
struct LocatedPeak
{
    double magnitude = 0;
    int64_t eighthFrame = 0;
    bool padding = false;
};
struct SectionPeaks
{
    std::vector<LocatedPeak> sample, truePeak;
};
struct SectionPeakMap
{
    juce::String geometry, snapshotIdentity;
    int channels = 0;
    std::vector<SectionPeaks> sections;
    size_t storageBytes() const;
};

// True-peak attribution expands each section by nominal 50 ms plus the 64-frame
// interpolation support. This deliberately conservative neighbourhood is not an
// event counter or an isolated-section measurement. Ordinary sample maxima use
// exact half-open section intervals. Retain actual signed/eighth-frame locations.
class SectionPeakCollector final : public measurement::PeakObserver
{
public:
    SectionPeakCollector(const SectionMap&, int channels, const juce::String& snapshot);
    void sample(size_t, double, int64_t) override;
    void truePeak(size_t, double, int64_t, bool) override;
    std::shared_ptr<const SectionPeakMap> finish();
private:
    SectionMap map;
    SectionPeakMap result;
    std::vector<int64_t> edges;
    std::vector<std::pair<int64_t,int64_t>> neighbourhoods;
    std::vector<LocatedPeak> buckets;
    std::array<size_t,measurement::channelCapacity> sampleCursor {}, peakCursor {};
    bool finished = false;
};
}

#pragma once
#include "Model/ProgramAnalysis.h"
#include "Model/TargetProfile.h"
#include <functional>

namespace levelscope::proposal::lra_timeline
{
inline constexpr const char* version = "measured-lra-timeline-v1";

// Each bit describes at least one measured window with positive-length source
// support in the queried interval. Overlapping windows may set several bits.
enum Flag : unsigned
{
    eligible = 1u, aboveP95 = 2u, upperRankTie = 4u, belowP10 = 8u,
    absoluteRejected = 16u, relativeRejected = 32u, zeroEnergy = 64u,
    insufficientDuration = 128u, padded = 256u, targetPassage = 512u,
    atP10 = 1024u, atP95 = 2048u
};
enum class TargetState { noCeiling, inapplicable, unavailable, unstable, within, exceeds };
struct Target
{
    TargetState state = TargetState::noCeiling;
    std::optional<double> ceilingLu, upperLevelLufs;
    juce::String snapshotIdentity;
};
struct Column
{
    unsigned flags = 0;
    std::optional<double> minLufs, maxLufs; // all finite measured windows touching the column
    std::optional<double> maxEligibleLufs;
};
struct Window
{
    size_t index = 0; // chronological index in ProgramAnalysis::lraWindows
    int64_t begin = 0, end = 0, measuredEnd = 0; // source [begin,end); measuredEnd includes drain
    std::optional<double> lufs;
    unsigned flags = 0;
};
class Index
{
public:
    // Worker only. Supply the accepted, validated full-resolution M4 analysis.
    static std::shared_ptr<const Index> build(std::shared_ptr<const ProgramAnalysis>,
                                              const std::function<bool()>& cancelled = [] { return false; });
    const ProgramAnalysis& measured() const { return *analysis; }
    size_t upperRankCount() const; // counts strictly later ranks, including LUFS ties assigned above P95
    Target target(const TargetProfile&, Applicability) const;
    // Returns one entry per requested pixel (up to 4096); invalid/empty viewport returns none.
    std::vector<Column> columns(int64_t begin, int64_t end, int pixels, const Target&) const;
    // All overlapping windows at an actual source frame, in time order (at most 30).
    std::vector<Window> inspect(int64_t frame, const Target&) const;
    size_t indexBytes() const;
private:
    explicit Index(std::shared_ptr<const ProgramAnalysis> a) : analysis(std::move(a)) {}
    static constexpr size_t blockSize = 32;
    std::shared_ptr<const ProgramAnalysis> analysis;
    size_t leafCount = 0;
    std::vector<unsigned short> flags;
    std::vector<Column> tree;
    std::pair<size_t,size_t> windowRange(int64_t begin, int64_t end) const;
    Column aggregate(size_t begin, size_t end) const;
};
}

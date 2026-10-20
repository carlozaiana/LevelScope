#pragma once
#include "StandaloneDiagnostics.h"
#include "Core/Proposal/SectionAnalysis.h"

namespace levelscope::standalone
{
// Content, source time base, ordered roles and range only. No gain/profile/video.
juce::String sectionInputIdentity (const ServiceContext&);
class SectionsController
{
public:
    ~SectionsController();
    void select (std::shared_ptr<const MeasuredDiagnostics>, const proposal::SectionMap&);
    void poll();
    void cancel();
    void retry();
    bool busy() const { return future.valid(); }
    std::shared_ptr<const proposal::SectionAnalysis> result() const { return accepted; }
    juce::String status = "Measure Source or Current State to build section diagnostics.";
    uint64_t buildsStarted = 0;
private:
    std::shared_ptr<const MeasuredDiagnostics> selected;
    proposal::SectionMap selectedMap;
    juce::String selectedKey;
    std::shared_ptr<JobControl> control;
    std::future<std::shared_ptr<const proposal::SectionAnalysis>> future;
    std::shared_ptr<const proposal::SectionAnalysis> accepted;
    uint64_t generation = 0, runningGeneration = 0;
    bool attempted = false;
};
}

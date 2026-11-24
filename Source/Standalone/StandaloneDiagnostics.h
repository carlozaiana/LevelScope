#pragma once
#include "StandaloneCompliance.h"
#include "StandaloneMeasurementService.h"
#include "Core/Proposal/AnalysisEngine.h"
#include "Core/Proposal/LraTimelineDiagnostics.h"
#include <future>

namespace levelscope::standalone
{
class StandaloneSessionModel;
enum class TimelineMeasurementState { unavailable, unmeasured, unverified, measuring, measured };
struct TimelineMeasurement
{
    TimelineMeasurementState state;
    std::shared_ptr<const AnalysisSnapshot> snapshot;
};
TimelineMeasurement timelineMeasurement(const StandaloneSessionModel&, MeasurementDomain, bool measuringThisDomain = false);

struct MeasuredDiagnostics
{
    std::shared_ptr<const AnalysisSnapshot> snapshot;
    std::shared_ptr<const proposal::ProgramAnalysis> analysis;
    proposal::MeasuredProgramme measuredBinding;
    std::shared_ptr<const proposal::lra_timeline::Index> lraTimeline;
};
struct DiagnosticAnnotations
{
    juce::String identity, text;
    std::optional<double> ceiling;
    std::vector<std::optional<double>> peakHeadroom;
};
DiagnosticAnnotations annotateDiagnostics (const MeasuredDiagnostics&, const TargetProfile&, proposal::Applicability);

// Independent, read-only job: never changes SessionJob, revisions, measurements,
// proposals or transport. A changed immutable snapshot cancels stale completion.
class DiagnosticsController
{
public:
    explicit DiagnosticsController(bool includeLraTimeline = false) : buildLraTimeline(includeLraTimeline) {}
    ~DiagnosticsController();
    void select (std::shared_ptr<const AnalysisSnapshot>, const ServiceContext&, MeasurementDomain);
    void poll();
    void cancel();
    void retry();
    bool busy() const { return future.valid(); }
    std::shared_ptr<const MeasuredDiagnostics> result() const { return accepted; }
    juce::String status = "No valid measured snapshot selected.";
    uint64_t buildsStarted = 0;
private:
    const bool buildLraTimeline;
    std::shared_ptr<const AnalysisSnapshot> selected, running;
    ServiceContext expected;
    MeasurementDomain domain = MeasurementDomain::source;
    juce::String expectedIdentity;
    std::shared_ptr<JobControl> control;
    std::future<ServiceResult<MeasuredDiagnostics>> future;
    std::shared_ptr<const MeasuredDiagnostics> accepted;
    uint64_t generation = 0, runningGeneration = 0;
    bool attempted = false;
};
}

#pragma once

#include "StandaloneMedia.h"
#include "StandaloneTargetProfiles.h"
#include <optional>
#include <vector>

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-SERVICE-BOUNDARY]
// Worker-thread APIs only. Requests/results are value snapshots, never UI/APVTS
// references. Callers own cancellation and reject stale revision/job tickets.
struct FrameRange { juce::int64 begin = 0, end = 0; }; // half-open, source-rate frames
struct ProcessingPlan
{
    int schemaVersion = 1;
    juce::String engineVersion;
    // Versioned module IDs/order/parameters/automation/freezes; no plugin state blob.
    juce::String graphJson;
};

struct ServiceContext
{
    JobTicket ticket;
    juce::File sourceFile;
    MediaInfo media;
    ChannelLayout fileLayout;
    FrameRange range;
    ProcessingPlan processing;
    TargetProfile profile;
    juce::String profileVersion, buildVersion, platformAndFloatingPointMode;

    bool isFullProgram() const { return range.begin == 0 && range.end == media.frames; }
    bool isValid() const
    {
        return ! media.sha256.isEmpty() && media.channels == fileLayout.count
            && fileLayout.confirmed && fileLayout.isValid()
            && range.begin >= 0 && range.end > range.begin && range.end <= media.frames;
    }
};

struct MeasuredMetrics
{
    // Missing metrics (including silence/short-file cases) must stay missing.
    std::optional<double> integratedLufs, loudnessRangeLu, truePeakDbtp, samplePeakDbfs;
};
enum class MeasurementDomain { source, currentState, renderedOutput };
struct AnalysisSnapshot
{
    ServiceContext context;
    MeasurementDomain domain = MeasurementDomain::source;
    MeasuredMetrics measured;
    juce::String measurementVersion;
    // Analysis engine may extend with curves, gates, attribution and scene data.
};
struct ProposalRequest
{
    ServiceContext context;
    AnalysisSnapshot source, current;
    juce::String methodId, policyMode, freezeConstraintsJson;
};
struct ProposalDraft
{
    ProcessingPlan plan;
    MeasuredMetrics predicted; // Never accepted as measured Current State.
    juce::String explanation;
};
struct RenderRequest
{
    ServiceContext context;
    juce::File destination; // Must be distinct from source; never overwrite source.
    juce::String format;
    int bitDepth = 0;
};
struct RenderArtifact
{
    juce::File file;
    juce::String sha256;
    AnalysisSnapshot verifiedOutput;
};
enum class ComplianceVerdict { notEvaluated, pass, fail, indeterminate };
struct ReportRequest
{
    ServiceContext context;
    AnalysisSnapshot source, current;
    RenderArtifact output;
    ComplianceVerdict verdict = ComplianceVerdict::notEvaluated;
    juce::File destination;
};
struct ReportArtifact { juce::File file; };

template <typename Value> struct ServiceResult
{
    std::optional<Value> value;
    ErrorCode error = ErrorCode::unavailable;
    juce::String message = "Service is not implemented";
};

class AnalysisService
{
public:
    virtual ~AnalysisService() = default;
    virtual ServiceResult<AnalysisSnapshot> analyze (const ServiceContext&, MeasurementDomain, JobControl&) = 0;
};
class ProposalService
{
public:
    virtual ~ProposalService() = default;
    virtual ServiceResult<ProposalDraft> propose (const ProposalRequest&, JobControl&) = 0;
};
class RenderService
{
public:
    virtual ~RenderService() = default;
    virtual ServiceResult<RenderArtifact> render (const RenderRequest&, JobControl&) = 0;
};
class ReportService
{
public:
    virtual ~ReportService() = default;
    virtual ServiceResult<ReportArtifact> write (const ReportRequest&, JobControl&) = 0;
};
// No implementations registered in the scaffold. No pass/fail from placeholder profiles.
// [END LS-STANDALONE-SERVICE-BOUNDARY]
}

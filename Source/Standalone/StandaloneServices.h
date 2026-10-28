#pragma once

#include "StandaloneMedia.h"
#include "ProgrammeSource.h"
#include "ProgrammeComposition.h"
#include "StandaloneTargetProfiles.h"
#include "Core/Measurement/OfflineMeasurement.h"
#include "Core/Processing/SectionGainPlan.h"
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
    juce::String engineVersion = "static-gain-v1";
    double staticGainDb = 0.0; // M1 typed plan; all file channels, including LFE.
    // Versioned module IDs/order/parameters/automation/freezes; no plugin state blob.
    juce::String graphJson;
    std::optional<processing::SectionGainPlan> sectionGains;
};

struct ServiceContext
{
    JobTicket ticket;
    juce::String sessionId;
    juce::File sourceFile;
    MediaInfo media;
    std::optional<ProgrammeAssembly> assembly;
    std::optional<ProgrammeComposition> composition; // immutable original Source recipe
    ChannelLayout fileLayout;
    FrameRange range;
    ProcessingPlan processing;
    TargetProfile profile;
    juce::String profileVersion, buildVersion, platformAndFloatingPointMode;

    bool isFullProgram() const { return range.begin == 0 && range.end == media.frames; }
    bool sourceFilesMatch() const { return composition ? composition->matchesFileStamps() : assembly ? assembly->matchesFileStamps() : media.matchesFileStamp(sourceFile); }
    std::vector<ProgrammeMember> sourceMembers() const
    {
        if(composition)return composition->files();
        if(assembly)return assembly->members;
        return {{"audio",sourceFile,media,fileLayout.roles[0],0}};
    }
    bool isValid() const
    {
        if(composition)
        {
            if(assembly||composition->validate().isNotEmpty())return false;
            const auto m=composition->summary();const auto roles=composition->layout();
            if(media.sha256!=m.sha256||media.channels!=m.channels||media.bytes!=m.bytes||media.frames!=m.frames
                ||media.sampleRate!=m.sampleRate||media.bitsPerSample!=m.bitsPerSample||media.format!=m.format
                ||fileLayout.count!=roles.count)return false;
            for(int i=0;i<roles.count;++i)if(fileLayout.roles[i]!=roles.roles[i])return false;
        }
        if (assembly)
        {
            if (assembly->validate().isNotEmpty()) return false;
            const auto summary = assembly->summary();
            if (media.sha256 != summary.sha256 || media.channels != summary.channels || media.bytes != summary.bytes
                || media.frames != summary.frames || media.sampleRate != summary.sampleRate || fileLayout.count != summary.channels) return false;
            for (int i=0; i<fileLayout.count; ++i) if (fileLayout.roles[i] != assembly->members[static_cast<size_t>(i)].role) return false;
        }
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
enum class MeasurementDomain { source, currentState, renderedOutput, privateCandidate };
struct AnalysisSnapshot
{
    ServiceContext context;
    MeasurementDomain domain = MeasurementDomain::source;
    MeasuredMetrics measured;
    juce::String measurementVersion;
    juce::String identity, measuredContentHash;
    measurement::Result detail;
    // Original measurement provenance survives runtime rebinding on session reopen.
    juce::String originSessionId;
    uint64_t originJobId = 0, originSessionRevision = 0;
    bool restoredFromSession = false;
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

// Identity binding only, not a compliance verdict or permission to publish.
// A future exporter must additionally obtain a passing authoritative verdict,
// recheck the artifact bytes, and reject destination aliases of original media.
struct ExportApproval
{
    juce::String sessionId, analysisIdentity, outputHash, profileId, profileVersion;
    uint64_t revision = 0;
    bool matches (const ServiceContext&, const RenderArtifact&) const;
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
// OfflineAnalysisService implements analysis. Proposal/render/report public services
// remain unavailable. Placeholder profiles never produce a compliance verdict.
// [END LS-STANDALONE-SERVICE-BOUNDARY]
}

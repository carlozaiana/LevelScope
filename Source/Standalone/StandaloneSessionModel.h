#pragma once

#include <juce_core/juce_core.h>
#include <cmath>
#include <memory>

#include "StandaloneTargetProfiles.h"
#include "StandaloneWorkflowReadiness.h"
#include "StandaloneMedia.h"
#include "StandaloneMeasurementService.h"
#include "StandaloneVideo.h"
#include "StandaloneCompliance.h"

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-WORKFLOW-MODEL]

enum class WorkflowPage
{
    importSource = 0,
    analyze,
    currentState,
    exportResults
};

struct MeasurementSummary
{
    bool hasMeasurement = false;

    std::optional<double> integratedLufs, loudnessRangeLu, truePeakDbtp;
    std::shared_ptr<const AnalysisSnapshot> analysis;

    juce::String statusText = "Not measured";
};

struct SourceDocument
{
    bool hasSource = false;

    juce::File file;
    juce::String displayName = "No source imported";
    juce::String fullPath;
    juce::String extension;
    juce::int64 fileSizeBytes = 0;
    MediaInfo media;
    ChannelLayout layout;
    OutputRouting routing;

    juce::String importStatus = "No source imported";

    void setValidatedFile (const juce::File& newFile, const MediaInfo& info)
    {
        file = newFile;
        hasSource = file.existsAsFile();

        if (! hasSource)
        {
            clear();
            return;
        }

        displayName = file.getFileName();
        fullPath = file.getFullPathName();
        extension = file.getFileExtension().toLowerCase();

        if (extension.startsWithChar ('.'))
            extension = extension.substring (1);

        fileSizeBytes = info.bytes;
        media = info;
        layout = info.suggestedLayout;
        routing = OutputRouting();
        importStatus = "Source decoded; confirm channel roles.";
    }

    void clear()
    {
        hasSource = false;
        file = juce::File();
        displayName = "No source imported";
        fullPath = {};
        extension = {};
        fileSizeBytes = 0;
        media = MediaInfo();
        layout = ChannelLayout();
        routing = OutputRouting();
        importStatus = "No source imported";
    }

};

struct CurrentStateDocument
{
    bool isInitialized = false;

    juce::String statusText = "No Current State initialized";
    juce::String basedOnSourceName = "No source";
    juce::String sourcePathSnapshot;

    bool needsRemeasurement = false;
    int revision = 0;
    double staticGainDb = 0.0;

    bool initializeFromSource (const SourceDocument& sourceDocument)
    {
        if (! sourceDocument.hasSource || ! sourceDocument.layout.confirmed || ! sourceDocument.layout.isValid())
        {
            clear();
            return false;
        }

        isInitialized = true;
        basedOnSourceName = sourceDocument.displayName;
        sourcePathSnapshot = sourceDocument.fullPath;
        needsRemeasurement = true;
        ++revision;

        staticGainDb = 0;
        statusText = "Identity processing; ready to measure";

        return true;
    }

    void clear()
    {
        isInitialized = false;
        statusText = "No Current State initialized";
        basedOnSourceName = "No source";
        sourcePathSnapshot = juce::String();
        needsRemeasurement = false;
        revision = 0;
        staticGainDb = 0;
    }
};

class StandaloneSessionModel
{
public:
    static constexpr int numTargetProfiles = StandaloneTargetProfileCatalog::numProfiles;

    WorkflowPage selectedPage = WorkflowPage::importSource;

    SourceDocument sourceDocument;
    CurrentStateDocument currentStateDocument;
    video::Reference referenceVideo; // Independent of audio revision, processing and routing.

    MeasurementSummary source;
    MeasurementSummary currentState;

    int selectedTargetProfileIndex = 0;
    proposal::Applicability programmeApplicability = proposal::Applicability::unconfirmed;
    uint64_t revision = 0;
    const juce::String sessionId = juce::Uuid().toString();
    SessionJob job;

    ServiceContext measurementContext (MeasurementDomain domain) const
    {
        ServiceContext c;
        c.ticket = job.ticket; c.ticket.sessionRevision = revision; c.sessionId = sessionId;
        c.sourceFile = sourceDocument.file; c.media = sourceDocument.media;
        c.fileLayout = sourceDocument.layout; c.range = { 0, c.media.frames };
        c.processing.staticGainDb = domain == MeasurementDomain::source ? 0 : currentStateDocument.staticGainDb;
        c.profile = getSelectedTargetProfile();
        c.profileVersion = c.profile.definitionVersion;
#ifdef LEVELSCOPE_MEASUREMENT_BUILD
        c.buildVersion = LEVELSCOPE_MEASUREMENT_BUILD;
#else
        c.buildVersion = "unidentified-development-build";
#endif
        c.platformAndFloatingPointMode = juce::SystemStats::getOperatingSystemName()
#if defined(__aarch64__) || defined(__arm64__)
            + ";arm64";
#else
            + ";x86_64";
#endif
        c.platformAndFloatingPointMode += ";float32-samples;float64-analysis;no-fast-math;fp-contract-off";
        return c;
    }

    bool setCurrentGainDb (double gain)
    {
        if (! currentStateDocument.isInitialized || ! std::isfinite (gain) || gain < -60 || gain > 24) return false;
        if (gain == currentStateDocument.staticGainDb) return true;
        ++revision; ++currentStateDocument.revision;
        currentStateDocument.staticGainDb = gain;
        currentStateDocument.needsRemeasurement = true;
        currentStateDocument.statusText = "Gain changed; re-measure Current State";
        currentState = MeasurementSummary();
        return true;
    }

    bool acceptMeasurement (AnalysisSnapshot result)
    {
        if (result.domain == MeasurementDomain::renderedOutput || ! sourceDocument.hasSource
            || result.context.sessionId != sessionId || result.context.ticket.sessionRevision != revision
            || result.context.ticket.id != job.ticket.id || job.state != JobState::succeeded
            || (result.domain == MeasurementDomain::currentState && ! currentStateDocument.isInitialized)) return false;
        const auto expected = measurementContext (result.domain);
        if (result.measuredContentHash != sourceDocument.media.sha256
            || result.identity != OfflineAnalysisService::identityFor (expected, result.domain, sourceDocument.media.sha256)) return false;
        auto& summary = result.domain == MeasurementDomain::source ? source : currentState;
        summary.hasMeasurement = true;
        summary.integratedLufs = result.measured.integratedLufs;
        summary.loudnessRangeLu = result.measured.loudnessRangeLu;
        summary.truePeakDbtp = result.measured.truePeakDbtp;
        summary.statusText = "Measured snapshot; selected rules evaluated separately";
        if (result.domain == MeasurementDomain::currentState)
        {
            currentStateDocument.needsRemeasurement = false;
            currentStateDocument.statusText = "Measured current revision";
        }
        summary.analysis = std::make_shared<const AnalysisSnapshot> (std::move (result));
        return true;
    }

    bool setValidatedSource (const juce::File& file, const MediaInfo& media)
    {
        if (! media.matchesFileStamp (file)) return false;
        ++revision;
        sourceDocument.setValidatedFile (file, media);

        source = MeasurementSummary();
        currentState = MeasurementSummary();
        currentStateDocument.clear();
        programmeApplicability = proposal::Applicability::unconfirmed;

        if (sourceDocument.hasSource)
        {
            source.statusText = "Validated source; measurement not run";
            currentState.statusText = "Waiting for Current State initialization";
        }
        return true;
    }

    void clearSourceFile()
    {
        ++revision;
        sourceDocument.clear();
        currentStateDocument.clear();
        programmeApplicability = proposal::Applicability::unconfirmed;

        source = MeasurementSummary();
        currentState = MeasurementSummary();

        selectedPage = WorkflowPage::importSource;
    }

    bool initializeCurrentStateFromSource()
    {
        ++revision;
        const auto didInitialize = currentStateDocument.initializeFromSource (sourceDocument);

        currentState = MeasurementSummary();

        if (didInitialize)
        {
            currentState.statusText = "Current State initialized; measurement not run";
            selectedPage = WorkflowPage::currentState;
        }

        return didInitialize;
    }

    void clearCurrentState()
    {
        ++revision;
        currentStateDocument.clear();
        currentState = MeasurementSummary();
    }

    bool canInitializeCurrentStateFromSource() const
    {
        return sourceDocument.hasSource && sourceDocument.layout.confirmed && sourceDocument.layout.isValid();
    }

    bool hasCurrentStateInitialized() const
    {
        return currentStateDocument.isInitialized;
    }

    bool currentStateNeedsRemeasurement() const
    {
        return currentStateDocument.isInitialized
            && currentStateDocument.needsRemeasurement;
    }

    void setSelectedTargetProfileIndex (int newIndex)
    {
        // Profile selection does not change samples or cancel a compatible audio
        // measurement. Evaluations are derived from the current profile on demand.
        selectedTargetProfileIndex = juce::jlimit (0, numTargetProfiles - 1, newIndex);
    }

    void setProgrammeApplicability (proposal::Applicability value) { programmeApplicability = value; }

    proposal::Evaluation evaluation (MeasurementDomain domain) const
    {
        const auto& summary = domain == MeasurementDomain::source ? source : currentState;
        const bool available = summary.hasMeasurement
            && (domain == MeasurementDomain::source || (currentStateDocument.isInitialized && ! currentStateDocument.needsRemeasurement));
        return evaluateSnapshot (available ? summary.analysis.get() : nullptr, measurementContext (domain),
                                 domain, getSelectedTargetProfile(), programmeApplicability);
    }

    void setChannelRole (int channel, ChannelRole role)
    {
        if (! juce::isPositiveAndBelow (channel, sourceDocument.layout.count)) return;
        ++revision;
        sourceDocument.layout.roles[channel] = role;
        sourceDocument.layout.confirmed = false;
        sourceDocument.importStatus = "File-channel roles changed; confirmation required.";
        source = MeasurementSummary();
        currentState = MeasurementSummary();
        currentStateDocument.clear();
    }

    bool confirmLayout()
    {
        if (! sourceDocument.hasSource || ! sourceDocument.layout.isValid()) return false;
        ++revision;
        sourceDocument.layout.confirmed = true;
        sourceDocument.importStatus = "Source decoded; file-channel roles confirmed.";
        return true;
    }

    TargetProfile getSelectedTargetProfile() const
    {
        return getTargetProfile (selectedTargetProfileIndex);
    }

    bool hasTargetProfileFamilySelected() const
    {
        return selectedTargetProfileIndex > 0
            && getSelectedTargetProfile().isSelectableForWorkflow;
    }

    bool hasAuthoritativeTargetLimits() const
    {
        return getSelectedTargetProfile().validationError().isEmpty();
    }

    WorkflowReadiness getWorkflowReadiness() const
    {
        WorkflowReadiness readiness;

        readiness.sourceImported = sourceDocument.hasSource;
        readiness.channelRolesConfirmed = sourceDocument.layout.confirmed && sourceDocument.layout.isValid();
        readiness.sourceMeasured = source.hasMeasurement;

        readiness.targetFamilySelected = hasTargetProfileFamilySelected();
        readiness.targetLimitsAuthoritative = hasAuthoritativeTargetLimits();

        readiness.currentStateInitialized = currentStateDocument.isInitialized;
        readiness.currentStateNeedsRemeasurement = currentStateNeedsRemeasurement();
        readiness.currentStateMeasured = currentState.hasMeasurement;

        readiness.proposalEngineAvailable = false;
        readiness.renderExportAvailable = false;

        return readiness;
    }

    static TargetProfile getTargetProfile (int index)
    {
        return StandaloneTargetProfileCatalog::getProfile (index);
    }
};

// [END LS-STANDALONE-WORKFLOW-MODEL]
}

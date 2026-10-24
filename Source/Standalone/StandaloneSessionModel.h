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
#include "StandaloneProposal.h"
#include "StandaloneAuditionRequest.h"
#include <functional>

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
    // Live controller guard; programmatic actions obey the same stopped-only rule.
    std::function<bool()> auditionBusy;
    MeasurementDomain auditionDomain = MeasurementDomain::source;
    bool audioEditsBlocked() const { return auditionBusy && auditionBusy(); }
    bool selectAuditionDomain (MeasurementDomain domain)
    {
        if (domain != MeasurementDomain::source && domain != MeasurementDomain::currentState) return false;
        if (domain == auditionDomain) return true;
        if (audioEditsBlocked()) return false;
        auditionDomain = domain; return true;
    }
    bool setOutputRoute (int channel, int output)
    {
        if (audioEditsBlocked() || ! juce::isPositiveAndBelow (channel, sourceDocument.layout.count)
            || output < -1 || output >= maxChannels) return false;
        sourceDocument.routing.fileToDevice[channel] = output; return true;
    }
    AuditionRequest auditionRequest() const
    {
        AuditionRequest r; r.context = measurementContext (auditionDomain); r.domain = auditionDomain;
        r.currentInitialized = currentStateDocument.isInitialized && ! sourceContentInvalid;
        if (sourceContentInvalid) r.context.media.sha256.clear();
        return r;
    }
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
    bool sourceContentInvalid = false;
    std::optional<proposal::TrimProposal> trimDraft;
    juce::String trimNotice = "No proposal drafted.";

    proposal::TrimRequest currentTrimRequest() const
    {
        return trimRequest (source.hasMeasurement && ! sourceContentInvalid ? source.analysis.get() : nullptr,
                            currentState.hasMeasurement && ! currentStateNeedsRemeasurement() ? currentState.analysis.get() : nullptr,
                            measurementContext (MeasurementDomain::source), measurementContext (MeasurementDomain::currentState),
                            getSelectedTargetProfile(), programmeApplicability);
    }

    bool acceptTrimProposal (const proposal::TrimProposal& result)
    {
        if (job.isActive()) return false;
        const auto fresh = proposal::draftGlobalTrim (currentTrimRequest());
        if (result.identity != fresh.identity) return false;
        trimDraft = fresh; // Reconstructed plan; never trust editable result fields.
        trimNotice = proposal::trimProposalText (fresh);
        return true;
    }

    bool canApplyTrim() const
    {
        return ! audioEditsBlocked() && ! job.isActive() && trimDraft && trimDraft->status == proposal::TrimStatus::draft
            && trimDraft->identity == proposal::draftGlobalTrim (currentTrimRequest()).identity;
    }

    bool applyTrim()
    {
        if (audioEditsBlocked() || job.isActive() || ! trimDraft) return false;
        const auto fresh = proposal::draftGlobalTrim (currentTrimRequest());
        if (trimDraft->identity != fresh.identity) { discardTrim(); return false; }
        if (fresh.status == proposal::TrimStatus::noChange) return true; // No audio mutation or remeasurement.
        if (fresh.status != proposal::TrimStatus::draft || ! fresh.totalGainDb) return false;
        if (! setCurrentGainDb (*fresh.totalGainDb)) return false;
        trimNotice = "Applied explicitly. Current State is UNVERIFIED; re-measure Current State. No automatic retry or rollback.\n"
                   + proposal::trimProposalText (fresh);
        return true;
    }

    void discardTrim()
    {
        trimDraft.reset();
        trimNotice = "No current proposal. Draft again from valid measured Current State.";
    }

    void measurementUnsuccessful (MeasurementDomain domain, ErrorCode error)
    {
        if (error == ErrorCode::sourceChanged)
        {
            ++revision; sourceContentInvalid = true; source = MeasurementSummary();
            currentState = MeasurementSummary(); discardTrim();
            sourceDocument.importStatus = "Content changed; import the source again";
            if (currentStateDocument.isInitialized) currentStateDocument.needsRemeasurement = true;
        }
        if (currentStateDocument.isInitialized && currentStateDocument.needsRemeasurement
            && (domain == MeasurementDomain::currentState || error == ErrorCode::sourceChanged))
        {
            currentState = MeasurementSummary();
            currentState.statusText = "UNVERIFIED: measurement did not complete successfully";
            currentStateDocument.statusText = "UNVERIFIED; applied gain retained. Explicit user action required.";
            trimNotice = "Verification unsuccessful. Applied gain retained; no automatic rollback, retry or correction.";
        }
    }

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
        if (audioEditsBlocked()) return false;
        discardTrim();
        ++revision; ++currentStateDocument.revision;
        currentStateDocument.staticGainDb = gain;
        currentStateDocument.needsRemeasurement = true;
        currentStateDocument.statusText = "UNVERIFIED: gain changed; re-measure Current State";
        currentState = MeasurementSummary();
        currentState.statusText = "UNVERIFIED: audio changed; previous measurement discarded";
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
        discardTrim();
        if (summary.analysis->domain == MeasurementDomain::currentState)
        {
            // The selected-profile verdict is derived above the status text; do
            // not cache it here across profile or applicability changes.
            currentStateDocument.statusText = "Measured current revision; selected rules evaluated separately";
            trimNotice = "Measurement complete. Actual Current State evaluation is authoritative. Gain retained; no automatic correction.";
        }
        return true;
    }

    bool setValidatedSource (const juce::File& file, const MediaInfo& media)
    {
        if (audioEditsBlocked() || ! media.matchesFileStamp (file)) return false;
        ++revision;
        discardTrim(); sourceContentInvalid = false;
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
        if (audioEditsBlocked()) return;
        ++revision;
        discardTrim(); sourceContentInvalid = false;
        sourceDocument.clear();
        currentStateDocument.clear();
        programmeApplicability = proposal::Applicability::unconfirmed;

        source = MeasurementSummary();
        currentState = MeasurementSummary();

        selectedPage = WorkflowPage::importSource;
    }

    bool initializeCurrentStateFromSource()
    {
        if (audioEditsBlocked() || sourceContentInvalid) return false;
        discardTrim();
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
        if (audioEditsBlocked()) return;
        discardTrim();
        ++revision;
        currentStateDocument.clear();
        currentState = MeasurementSummary();
    }

    bool canInitializeCurrentStateFromSource() const
    {
        return ! sourceContentInvalid && sourceDocument.hasSource && sourceDocument.layout.confirmed && sourceDocument.layout.isValid();
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
        const auto index = juce::jlimit (0, numTargetProfiles - 1, newIndex);
        if (selectedTargetProfileIndex != index) discardTrim();
        selectedTargetProfileIndex = index;
    }

    void setProgrammeApplicability (proposal::Applicability value)
    {
        if (programmeApplicability != value) discardTrim();
        programmeApplicability = value;
    }

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
        if (audioEditsBlocked() || ! juce::isPositiveAndBelow (channel, sourceDocument.layout.count)) return;
        discardTrim();
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
        if (audioEditsBlocked() || ! sourceDocument.hasSource || ! sourceDocument.layout.isValid()) return false;
        discardTrim();
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

        readiness.proposalEngineAvailable = true;
        readiness.proposalInputsValidated = ! job.isActive()
            && proposal::draftGlobalTrim (currentTrimRequest()).status != proposal::TrimStatus::unavailable;
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

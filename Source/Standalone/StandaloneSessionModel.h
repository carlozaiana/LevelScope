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
#include "StandaloneSections.h"
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
    std::optional<ProgrammeAssembly> assembly;
    std::optional<ProgrammeComposition> composition;
    ChannelLayout layout;
    OutputRouting routing;

    juce::String importStatus = "No source imported";
    bool matchesFileStamps() const { return composition ? composition->matchesFileStamps() : assembly ? assembly->matchesFileStamps() : media.matchesFileStamp(file); }
    juce::String channelName (int channel) const
    {
        if(composition)return "Channel "+juce::String(channel+1)+" (M&E + dialogue composition)";
        return assembly && juce::isPositiveAndBelow(channel, static_cast<int>(assembly->members.size()))
            ? assembly->members[static_cast<size_t>(channel)].file.getFileName() : "Channel " + juce::String(channel+1);
    }
    void setValidatedProgramme (const ProgrammeAssembly& group)
    {
        clear(); assembly = group; media = group.summary(); layout = group.layout();
        hasSource = true; fileSizeBytes = media.bytes; extension = "programme";
        displayName = juce::String(media.channels) + "-file " + (media.channels == 2 ? "stereo" : "5.1") + " programme";
        fullPath = "Assembled programme";
        importStatus = "Complete aligned group confirmed; metadata offsets are not applied";
    }
    void setValidatedComposition (const ProgrammeComposition& recipe)
    {
        clear();composition=recipe;media=recipe.summary();layout=recipe.layout();hasSource=true;
        fileSizeBytes=media.bytes;extension="composition";displayName="Composed programme - M&E + dialogue";
        fullPath="Original confirmed composition";
        importStatus="Complementary, correctly balanced stems and sample-zero alignment confirmed; unity component gains";
    }

    void setValidatedFile (const juce::File& newFile, const MediaInfo& info)
    {
        assembly.reset();
        composition.reset();
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
        assembly.reset();
        composition.reset();
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
    std::optional<processing::SectionGainPlan> sectionGains;
    proposal::SectionMap sectionGainAnnotations; // provenance, never audible identity

    ProcessingPlan processingPlan() const
    {
        ProcessingPlan p; p.staticGainDb=staticGainDb; p.sectionGains=sectionGains;
        if(sectionGains) { p.schemaVersion=2; p.engineVersion=processing::SectionGainPlan::version; }
        return p;
    }

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

        staticGainDb = 0; sectionGains.reset(); sectionGainAnnotations={};
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
        staticGainDb = 0; sectionGains.reset(); sectionGainAnnotations={};
    }
};

class StandaloneSessionModel
{
public:
    // Live controller guard; programmatic actions obey the same stopped-only rule.
    std::function<bool()> auditionBusy;
    MeasurementDomain auditionDomain = MeasurementDomain::source;
    bool persistenceLoading = false;
    bool audioEditsBlocked() const { return persistenceLoading || (auditionBusy && auditionBusy()); }
    bool selectAuditionDomain (MeasurementDomain domain)
    {
        if (domain != MeasurementDomain::source && domain != MeasurementDomain::currentState) return false;
        if (domain == auditionDomain) return true;
        if (audioEditsBlocked()) return false;
        auditionDomain = domain; documentChanged(); return true;
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
    std::shared_ptr<const AnalysisSnapshot> unresolvedSourceCache, unresolvedCurrentCache;

    int selectedTargetProfileIndex = 0;
    proposal::Applicability programmeApplicability = proposal::Applicability::unconfirmed;
    uint64_t revision = 0;
    juce::String sessionId = juce::Uuid().toString(); // unique runtime opening, never loaded from disk
    juce::String documentId = juce::Uuid().toString(), programmeId = juce::Uuid().toString();
    uint64_t documentRevision = 0, savedDocumentRevision = 0;
    juce::File documentFile;
    int openedSessionSchema = 4;
    juce::String restoredProfileId, restoredProfileVersion, restoredProfileDigest;
    bool videoUnresolved = false;
    proposal::SectionDocument sectionDocument;
    struct SectionTrimDraft
    {
        processing::SectionGainPlan plan;
        proposal::SectionMap annotations;
        juce::String geometry, runtime;
        uint64_t audioRevision=0;
    };
    std::optional<SectionTrimDraft> sectionTrimDraft;
    juce::String sectionTrimNotice="No section-gain draft. Applied processing is independent of annotations.";
    static juce::String sectionGeometry(const proposal::SectionMap& map)
    {
        juce::String key=map.inputIdentity+":"+juce::String(map.sampleRate)+":"+juce::String(map.frames);
        for(const auto& row:map.sections) key += ";"+juce::String(row.id)+":"+juce::String(row.begin)+":"+juce::String(row.end);
        return key;
    }
    bool sectionDraftFresh() const
    {
        return sectionTrimDraft && sectionTrimDraft->runtime==sessionId && sectionTrimDraft->audioRevision==revision
            && sectionTrimDraft->geometry==sectionGeometry(sectionDocument.map()) && !sourceContentInvalid
            && sectionTrimDraft->plan.inputIdentity==sectionInputIdentity(measurementContext(MeasurementDomain::source));
    }
    bool beginSectionTrimDraft()
    {
        if(audioEditsBlocked() || job.isActive() || sectionTrimDraft || !currentStateDocument.isInitialized || sourceContentInvalid) return false;
        synchronizeSections(); const auto& map=sectionDocument.map();
        if(!map.valid()) return false;
        auto c=measurementContext(MeasurementDomain::currentState);
        // Validate the layout through the same adapter, using a non-unity probe.
        processing::SectionGainPlan probe; probe.inputIdentity=map.inputIdentity; probe.sampleRate=map.sampleRate;
        probe.frames=map.frames; probe.intervals={{0,map.frames,1}};
        c.processing.sectionGains=probe;c.processing.schemaVersion=2;c.processing.engineVersion=processing::SectionGainPlan::version;
        const auto error=ProgrammeGainProcessor::validationError(c);
        if(error.isNotEmpty()){sectionTrimNotice=error;return false;}
        SectionTrimDraft d;d.annotations=map;d.geometry=sectionGeometry(map);d.runtime=sessionId;d.audioRevision=revision;
        d.plan=probe;d.plan.intervals.clear();
        const auto& active=currentStateDocument.sectionGains;
        const bool sameGeometry=active && sectionGeometry(currentStateDocument.sectionGainAnnotations)==d.geometry;
        for(const auto& row:map.sections)
        {
            double gain=0;
            if(sameGeometry) for(const auto& applied:active->intervals)
                if(applied.begin<=row.begin && applied.end>=row.end) {gain=applied.db;break;}
            d.plan.intervals.push_back({row.begin,row.end,gain});
        }
        sectionTrimDraft=std::move(d);
        sectionTrimNotice=sameGeometry?"DRAFT copied from applied gains; no audio change.":"DRAFT initialized to 0 dB on current annotations; no audio change.";
        return true;
    }
    bool setSectionDraftGain(size_t row,double db)
    {
        if(audioEditsBlocked() || job.isActive() || !sectionDraftFresh() || row>=sectionTrimDraft->plan.intervals.size()
            || !std::isfinite(db) || db < -12 || db > 12) return false;
        sectionTrimDraft->plan.intervals[row].db=db;
        sectionTrimNotice="DRAFT only: Apply explicitly while stopped, then re-measure Current State.";return true;
    }
    void discardSectionTrimDraft()
    {sectionTrimDraft.reset();sectionTrimNotice="Draft discarded. Applied processing and measurements unchanged.";}
    bool applySectionTrimDraft()
    {
        if(audioEditsBlocked() || job.isActive() || !sectionDraftFresh())
        {sectionTrimNotice="Cannot apply: stop playback, finish jobs and discard/recreate any stale draft.";return false;}
        const auto error=sectionTrimDraft->plan.validationError();
        if(error.isNotEmpty()){sectionTrimNotice=error;return false;}
        auto next=sectionTrimDraft->plan.canonical();
        auto c=measurementContext(MeasurementDomain::currentState);
        c.processing.sectionGains=next.isUnity()?std::nullopt:std::optional<processing::SectionGainPlan>(next);
        c.processing.schemaVersion=next.isUnity()?1:2;c.processing.engineVersion=next.isUnity()?"static-gain-v1":processing::SectionGainPlan::version;
        const auto requestError=OfflineAnalysisService::validateRequest(c,MeasurementDomain::currentState);
        if(requestError.isNotEmpty()){sectionTrimNotice=requestError;return false;}
        const auto oldKey=currentStateDocument.sectionGains?currentStateDocument.sectionGains->audioKey():juce::String();
        const auto newKey=next.isUnity()?juce::String():next.audioKey();
        if(oldKey==newKey)
        {sectionTrimDraft.reset();sectionTrimNotice="NO CHANGE: applied samples, revision and valid measurements retained.";return true;}
        currentStateDocument.sectionGains=c.processing.sectionGains;
        currentStateDocument.sectionGainAnnotations=next.isUnity()?proposal::SectionMap{}:sectionTrimDraft->annotations;
        sectionTrimDraft.reset();invalidateCurrentProcessing();
        sectionTrimNotice="APPLIED; Current State UNVERIFIED. Re-measure explicitly; no automatic correction or rollback.";return true;
    }
    void invalidateCurrentProcessing()
    {
        discardTrim();++revision;++currentStateDocument.revision;documentChanged();unresolvedCurrentCache.reset();
        currentStateDocument.needsRemeasurement=true;
        currentStateDocument.statusText="UNVERIFIED: processing changed; re-measure Current State";
        currentState=MeasurementSummary();currentState.statusText="UNVERIFIED: audio changed; previous measurement discarded";
    }
    void documentChanged() { ++documentRevision; }
    bool documentDirty() const { return sectionTrimDraft.has_value() || documentRevision != savedDocumentRevision; }
    void synchronizeSections()
    {
        const auto c = measurementContext (MeasurementDomain::source);
        sectionDocument.bind (sectionInputIdentity (c), static_cast<int> (c.media.sampleRate), c.media.frames);
    }
    void newDocument()
    {
        const auto guard = auditionBusy;
        *this = StandaloneSessionModel(); auditionBusy = guard;
    }
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
            if (sourceDocument.assembly || sourceDocument.composition)
            {
                unresolvedSourceCache = source.analysis;
                unresolvedCurrentCache = currentState.analysis;
            }
            ++revision; documentChanged(); sourceContentInvalid = true; source = MeasurementSummary();
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
        c.sourceFile = sourceDocument.file; c.media = sourceDocument.media; c.assembly = sourceDocument.assembly;
        c.composition = sourceDocument.composition;
        c.fileLayout = sourceDocument.layout; c.range = { 0, c.media.frames };
        c.processing = domain == MeasurementDomain::source ? ProcessingPlan{} : currentStateDocument.processingPlan();
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
        currentStateDocument.staticGainDb = gain;
        invalidateCurrentProcessing();
        return true;
    }

    bool acceptMeasurement (AnalysisSnapshot result)
    {
        if (result.domain == MeasurementDomain::renderedOutput || ! sourceDocument.hasSource
            || result.context.sessionId != sessionId || result.context.ticket.sessionRevision != revision
            || result.context.ticket.id != job.ticket.id || job.state != JobState::succeeded
            || (result.domain == MeasurementDomain::currentState && ! currentStateDocument.isInitialized)) return false;
        if ((sourceDocument.assembly || sourceDocument.composition) && !sourceDocument.matchesFileStamps())
        { measurementUnsuccessful(result.domain, ErrorCode::sourceChanged); return false; }
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
        documentChanged();
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
        ++revision; documentChanged(); programmeId = juce::Uuid().toString();
        discardTrim(); sourceContentInvalid = false; unresolvedSourceCache.reset(); unresolvedCurrentCache.reset();
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

    bool setValidatedProgramme (const ProgrammeAssembly& group)
    {
        if (audioEditsBlocked() || group.validate().isNotEmpty() || !group.matchesFileStamps()) return false;
        ++revision; documentChanged(); programmeId = juce::Uuid().toString();
        discardTrim(); sourceContentInvalid = false; unresolvedSourceCache.reset(); unresolvedCurrentCache.reset(); sourceDocument.setValidatedProgramme(group);
        source = MeasurementSummary(); currentState = MeasurementSummary(); currentStateDocument.clear();
        programmeApplicability = proposal::Applicability::unconfirmed;
        source.statusText = "Complete aligned programme; measurement not run";
        return true;
    }

    bool setValidatedComposition (const ProgrammeComposition& recipe)
    {
        if(audioEditsBlocked()||job.isActive()||recipe.validate().isNotEmpty()||!recipe.matchesFileStamps())return false;
        const bool editing=sourceDocument.composition.has_value();
        if(editing&&!sourceContentInvalid&&sourceDocument.composition->definitionIdentity()==recipe.definitionIdentity())
        {
            // Relinking/labels do not change audible input or the Current plan.
            sourceDocument.composition=recipe;sourceDocument.media=recipe.summary();sourceDocument.layout=recipe.layout();
            documentChanged();return true;
        }
        const bool keepCurrent=editing&&currentStateDocument.isInitialized;
        const double retainedGain=currentStateDocument.staticGainDb;
        ++revision;documentChanged();if(!editing)programmeId=juce::Uuid().toString();
        discardTrim();sourceContentInvalid=false;unresolvedSourceCache.reset();unresolvedCurrentCache.reset();
        sourceDocument.setValidatedComposition(recipe);source=MeasurementSummary();currentState=MeasurementSummary();
        currentStateDocument.clear();programmeApplicability=proposal::Applicability::unconfirmed;
        if(keepCurrent){currentStateDocument.initializeFromSource(sourceDocument);currentStateDocument.staticGainDb=retainedGain;}
        source.statusText="Original confirmed composition; measurement not run";
        currentState.statusText=keepCurrent?"UNVERIFIED: composition changed; applied global gain retained":"Waiting for Current State initialization";
        sectionDocument=proposal::SectionDocument();
        return true;
    }

    void clearSourceFile()
    {
        if (audioEditsBlocked()) return;
        ++revision; documentChanged(); programmeId = juce::Uuid().toString();
        discardTrim(); sourceContentInvalid = false; unresolvedSourceCache.reset(); unresolvedCurrentCache.reset();
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
        ++revision; documentChanged();
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
        unresolvedCurrentCache.reset();
        discardTrim();
        ++revision; documentChanged();
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
        if (persistenceLoading) return;
        // Profile selection does not change samples or cancel a compatible audio
        // measurement. Evaluations are derived from the current profile on demand.
        const auto index = juce::jlimit (0, numTargetProfiles - 1, newIndex);
        if (selectedTargetProfileIndex != index || restoredProfileId.isNotEmpty())
        { discardTrim(); documentChanged(); restoredProfileId.clear(); restoredProfileVersion.clear(); restoredProfileDigest.clear(); }
        selectedTargetProfileIndex = index;
    }

    void setProgrammeApplicability (proposal::Applicability value)
    {
        if (persistenceLoading) return;
        if (programmeApplicability != value) { discardTrim(); documentChanged(); }
        programmeApplicability = value;
    }

    // Accepted analyses are immutable. Cache only their derived evaluation, keyed
    // by the exact live context and profile; never deserialize a saved verdict.
    struct EvaluationCache
    {
        std::shared_ptr<const AnalysisSnapshot> snapshot;
        juce::String key;
        proposal::Evaluation value;
    };
    mutable std::array<EvaluationCache, 2> evaluations;
    proposal::Evaluation evaluation (MeasurementDomain domain) const
    {
        const auto& summary = domain == MeasurementDomain::source ? source : currentState;
        const bool available = !sourceContentInvalid && summary.hasMeasurement
            && (domain == MeasurementDomain::source || (currentStateDocument.isInitialized && ! currentStateDocument.needsRemeasurement));
        const auto analysis = available ? summary.analysis : nullptr;
        const auto context = measurementContext(domain);
        const auto profile = getSelectedTargetProfile();
        const auto key = context.sessionId + ":" + OfflineAnalysisService::identityFor(context, domain, context.media.sha256)
            + ":" + profile.digest() + ":" + juce::String(static_cast<int>(programmeApplicability));
        auto& cached = evaluations[domain == MeasurementDomain::source ? 0 : 1];
        if (cached.key != key || cached.snapshot != analysis)
        {
            cached.value = evaluateSnapshot(analysis.get(), context, domain, profile, programmeApplicability);
            cached.key = key; cached.snapshot = analysis;
        }
        return cached.value;
    }

    void setChannelRole (int channel, ChannelRole role)
    {
        // Assembly role reassignment is staged and explicitly reconfirmed as a
        // whole group, never performed through the interleaved role editor.
        if (sourceDocument.assembly || sourceDocument.composition) return;
        if (audioEditsBlocked() || ! juce::isPositiveAndBelow (channel, sourceDocument.layout.count)) return;
        discardTrim();
        ++revision; documentChanged();
        sourceDocument.layout.roles[channel] = role;
        sourceDocument.layout.confirmed = false;
        sourceDocument.layout.manuallyInterpreted = true;
        sourceDocument.importStatus = "Manual file-channel interpretation; confirmation required.";
        source = MeasurementSummary();
        currentState = MeasurementSummary();
        currentStateDocument.clear();
    }

    bool confirmLayout()
    {
        if (sourceDocument.composition) return sourceDocument.composition->validate().isEmpty();
        if (sourceDocument.assembly) return sourceDocument.assembly->validate().isEmpty();
        if (audioEditsBlocked() || ! sourceDocument.hasSource || ! sourceDocument.layout.isValid()) return false;
        discardTrim();
        ++revision; documentChanged();
        sourceDocument.layout.confirmed = true;
        sourceDocument.importStatus = sourceDocument.layout.manuallyInterpreted
            ? "Source decoded; explicit manual file-channel interpretation confirmed."
            : "Source decoded; format-suggested file-channel roles confirmed.";
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

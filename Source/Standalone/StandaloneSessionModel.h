#pragma once

#include <juce_core/juce_core.h>

#include "StandaloneTargetProfiles.h"
#include "StandaloneWorkflowReadiness.h"
#include "StandaloneMedia.h"

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

    double integratedLufs = 0.0;
    double loudnessRangeLu = 0.0;
    double truePeakDbtp = 0.0;

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
        importStatus = "Source decoded; confirm channel roles. Measurement not run.";
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

        statusText = "Initialized from Source; re-measure not implemented";

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
    }
};

class StandaloneSessionModel
{
public:
    static constexpr int numTargetProfiles = StandaloneTargetProfileCatalog::numProfiles;

    WorkflowPage selectedPage = WorkflowPage::importSource;

    SourceDocument sourceDocument;
    CurrentStateDocument currentStateDocument;

    MeasurementSummary source;
    MeasurementSummary currentState;

    int selectedTargetProfileIndex = 0;
    uint64_t revision = 0;
    SessionJob job;

    bool setValidatedSource (const juce::File& file, const MediaInfo& media)
    {
        if (! media.matchesFileStamp (file)) return false;
        ++revision;
        sourceDocument.setValidatedFile (file, media);

        source = MeasurementSummary();
        currentState = MeasurementSummary();
        currentStateDocument.clear();

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
        ++revision;
        selectedTargetProfileIndex = juce::jlimit (0, numTargetProfiles - 1, newIndex);
    }

    void setChannelRole (int channel, ChannelRole role)
    {
        if (! juce::isPositiveAndBelow (channel, sourceDocument.layout.count)) return;
        ++revision;
        sourceDocument.layout.roles[channel] = role;
        sourceDocument.layout.confirmed = false;
        sourceDocument.importStatus = "File-channel roles changed; confirmation required. Not measured.";
        source = MeasurementSummary();
        currentState = MeasurementSummary();
        currentStateDocument.clear();
    }

    bool confirmLayout()
    {
        if (! sourceDocument.hasSource || ! sourceDocument.layout.isValid()) return false;
        ++revision;
        sourceDocument.layout.confirmed = true;
        sourceDocument.importStatus = "Source decoded; file-channel roles confirmed. Not measured.";
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
        return getSelectedTargetProfile().limits.valuesAreAuthoritative;
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

#pragma once

#include <juce_core/juce_core.h>

#include "StandaloneTargetProfiles.h"

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

    juce::String importStatus = "No source imported";

    void setFromFile (const juce::File& newFile)
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

        fileSizeBytes = file.getSize();

        importStatus = hasSupportedAudioExtension()
            ? "Source selected; measurement not run"
            : "Source selected; extension not recognized by scaffold list";
    }

    void clear()
    {
        hasSource = false;
        file = juce::File();
        displayName = "No source imported";
        fullPath = {};
        extension = {};
        fileSizeBytes = 0;
        importStatus = "No source imported";
    }

    bool hasSupportedAudioExtension() const
    {
        return extension == "wav"
            || extension == "wave"
            || extension == "aif"
            || extension == "aiff"
            || extension == "flac"
            || extension == "mp3"
            || extension == "m4a";
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
        if (! sourceDocument.hasSource)
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

    void setSourceFile (const juce::File& file)
    {
        sourceDocument.setFromFile (file);

        source = MeasurementSummary();
        currentState = MeasurementSummary();
        currentStateDocument.clear();

        if (sourceDocument.hasSource)
        {
            source.statusText = sourceDocument.importStatus;
            currentState.statusText = "Waiting for Current State initialization";
        }
    }

    void clearSourceFile()
    {
        sourceDocument.clear();
        currentStateDocument.clear();

        source = MeasurementSummary();
        currentState = MeasurementSummary();

        selectedPage = WorkflowPage::importSource;
    }

    bool initializeCurrentStateFromSource()
    {
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
        currentStateDocument.clear();
        currentState = MeasurementSummary();
    }

    bool canInitializeCurrentStateFromSource() const
    {
        return sourceDocument.hasSource;
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
        selectedTargetProfileIndex = juce::jlimit (0, numTargetProfiles - 1, newIndex);
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

    static TargetProfile getTargetProfile (int index)
    {
        return StandaloneTargetProfileCatalog::getProfile (index);
    }
};

// [END LS-STANDALONE-WORKFLOW-MODEL]
}
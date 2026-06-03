#pragma once

#include <juce_core/juce_core.h>

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

struct TargetProfile
{
    juce::String id;
    juce::String displayName;
    juce::String description;
};

class StandaloneSessionModel
{
public:
    static constexpr int numTargetProfiles = 4;

    WorkflowPage selectedPage = WorkflowPage::importSource;

    SourceDocument sourceDocument;

    MeasurementSummary source;
    MeasurementSummary currentState;

    int selectedTargetProfileIndex = 0;

    void setSourceFile (const juce::File& file)
    {
        sourceDocument.setFromFile (file);

        source = MeasurementSummary();
        currentState = MeasurementSummary();

        if (sourceDocument.hasSource)
        {
            source.statusText = sourceDocument.importStatus;
            currentState.statusText = "Waiting for source measurement";
        }
    }

    void clearSourceFile()
    {
        sourceDocument.clear();

        source = MeasurementSummary();
        currentState = MeasurementSummary();

        selectedPage = WorkflowPage::importSource;
    }

    void setSelectedTargetProfileIndex (int newIndex)
    {
        selectedTargetProfileIndex = juce::jlimit (0, numTargetProfiles - 1, newIndex);
    }

    TargetProfile getSelectedTargetProfile() const
    {
        return getTargetProfile (selectedTargetProfileIndex);
    }

    static TargetProfile getTargetProfile (int index)
    {
        switch (index)
        {
            case 1:
                return {
                    "ebu-r128-placeholder",
                    "EBU R128 family (placeholder)",
                    "Target profile scaffold only. No compliance values are enforced yet."
                };

            case 2:
                return {
                    "atsc-a85-placeholder",
                    "ATSC A/85 family (placeholder)",
                    "Target profile scaffold only. No compliance values are enforced yet."
                };

            case 3:
                return {
                    "custom-placeholder",
                    "Custom target (placeholder)",
                    "Reserved for future user-defined target settings."
                };

            case 0:
            default:
                return {
                    "none",
                    "No target selected",
                    "Choose a target profile before proposal or export features are enabled."
                };
        }
    }
};

// [END LS-STANDALONE-WORKFLOW-MODEL]
}
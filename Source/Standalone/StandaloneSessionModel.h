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

    juce::String sourceDisplayName = "No source imported";

    MeasurementSummary source;
    MeasurementSummary currentState;

    int selectedTargetProfileIndex = 0;

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
#pragma once

#include <juce_core/juce_core.h>

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-WORKFLOW-READINESS]

struct WorkflowReadiness
{
    bool sourceImported = false;
    bool sourceMeasured = false;

    bool targetFamilySelected = false;
    bool targetLimitsAuthoritative = false;

    bool currentStateInitialized = false;
    bool currentStateNeedsRemeasurement = false;
    bool currentStateMeasured = false;

    bool proposalEngineAvailable = false;
    bool renderExportAvailable = false;

    static juce::String yesNo (bool value)
    {
        return value ? "yes" : "no";
    }

    bool isSourceReadyForFutureMeasurement() const
    {
        return sourceImported;
    }

    bool areProposalInputsSatisfied() const
    {
        return sourceMeasured
            && targetFamilySelected
            && targetLimitsAuthoritative
            && currentStateInitialized
            && currentStateMeasured
            && ! currentStateNeedsRemeasurement;
    }

    bool canDraftProposal() const
    {
        return proposalEngineAvailable
            && areProposalInputsSatisfied();
    }

    bool canExport() const
    {
        return renderExportAvailable
            && sourceMeasured
            && targetFamilySelected
            && targetLimitsAuthoritative
            && currentStateInitialized
            && currentStateMeasured
            && ! currentStateNeedsRemeasurement;
    }

    juce::String buildChecklistText() const
    {
        juce::String text;

        text << "Readiness checklist:\n";
        text << "- Source imported: " << yesNo (sourceImported) << "\n";
        text << "- Source measured: " << yesNo (sourceMeasured) << "\n";
        text << "- Target family selected: " << yesNo (targetFamilySelected) << "\n";
        text << "- Authoritative target limits loaded: " << yesNo (targetLimitsAuthoritative) << "\n";
        text << "- Current State initialized: " << yesNo (currentStateInitialized) << "\n";
        text << "- Current State needs re-measurement: " << yesNo (currentStateNeedsRemeasurement) << "\n";
        text << "- Current State measured: " << yesNo (currentStateMeasured) << "\n";
        text << "- Proposal engine available: " << yesNo (proposalEngineAvailable) << "\n";
        text << "- Render/export available: " << yesNo (renderExportAvailable) << "\n\n";

        text << "Can draft proposal: " << yesNo (canDraftProposal()) << "\n";
        text << "Can export: " << yesNo (canExport()) << "\n";

        return text;
    }

    juce::String buildNextActionText() const
    {
        if (! sourceImported)
            return "Next action: import a Source file.";

        if (! sourceMeasured)
            return "Next action: implement and run Source measurement.";

        if (! targetFamilySelected)
            return "Next action: select a Target profile family.";

        if (! targetLimitsAuthoritative)
            return "Next action: load authoritative Target profile limits.";

        if (! currentStateInitialized)
            return "Next action: initialize Current State from Source.";

        if (currentStateNeedsRemeasurement)
            return "Next action: implement and run Current State re-measurement.";

        if (! currentStateMeasured)
            return "Next action: measure Current State.";

        if (! proposalEngineAvailable)
            return "Next action: proposal engine remains intentionally unavailable in this scaffold.";

        if (! renderExportAvailable)
            return "Next action: render/export remains intentionally unavailable in this scaffold.";

        return "Next action: workflow gates are satisfied.";
    }
};

// [END LS-STANDALONE-WORKFLOW-READINESS]
}
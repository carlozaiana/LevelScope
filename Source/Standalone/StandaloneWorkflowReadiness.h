#pragma once

#include <juce_core/juce_core.h>

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-WORKFLOW-READINESS]

struct WorkflowReadiness
{
    bool sourceImported = false;
    bool channelRolesConfirmed = false;
    bool sourceMeasured = false;

    bool targetFamilySelected = false;
    bool targetLimitsAuthoritative = false;

    bool currentStateInitialized = false;
    bool currentStateNeedsRemeasurement = false;
    bool currentStateMeasured = false;

    bool proposalEngineAvailable = false;
    bool proposalInputsValidated = false;
    bool renderExportAvailable = false;
    bool exactVerificationPassed = false;
    bool exportApproved = false;

    static juce::String yesNo (bool value)
    {
        return value ? "yes" : "no";
    }

    bool isSourceReadyForFutureMeasurement() const
    {
        return sourceImported && channelRolesConfirmed;
    }

    bool areProposalInputsSatisfied() const
    {
        return isSourceReadyForFutureMeasurement() && sourceMeasured
            && targetFamilySelected
            && targetLimitsAuthoritative
            && currentStateInitialized
            && currentStateMeasured
            && ! currentStateNeedsRemeasurement;
    }

    bool canDraftProposal() const
    {
        return proposalEngineAvailable
            && proposalInputsValidated
            && areProposalInputsSatisfied();
    }

    bool canExport() const
    {
        return renderExportAvailable
            && exactVerificationPassed && exportApproved
            && isSourceReadyForFutureMeasurement()
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
        text << "- File-channel roles confirmed: " << yesNo (channelRolesConfirmed) << "\n";
        text << "- Source measured: " << yesNo (sourceMeasured) << "\n";
        text << "- Target family selected: " << yesNo (targetFamilySelected) << "\n";
        text << "- Authoritative target limits loaded: " << yesNo (targetLimitsAuthoritative) << "\n";
        text << "- Current State initialized: " << yesNo (currentStateInitialized) << "\n";
        text << "- Current State needs re-measurement: " << yesNo (currentStateNeedsRemeasurement) << "\n";
        text << "- Current State measured: " << yesNo (currentStateMeasured) << "\n";
        text << "- Proposal engine available: " << yesNo (proposalEngineAvailable) << "\n";
        text << "- Matching measured proposal inputs and applicability: " << yesNo (proposalInputsValidated) << "\n";
        text << "- Render/export available: " << yesNo (renderExportAvailable) << "\n\n";
        text << "- Exact rendered-output verification passed: " << yesNo (exactVerificationPassed) << "\n";
        text << "- Export approved for this revision: " << yesNo (exportApproved) << "\n\n";

        text << "Can draft proposal: " << yesNo (canDraftProposal()) << "\n";
        text << "Can export: " << yesNo (canExport()) << "\n";

        return text;
    }

    juce::String buildNextActionText() const
    {
        if (! sourceImported)
            return "Next action: import a Source file.";

        if (! sourceMeasured)
        {
            if (! channelRolesConfirmed)
                return "Next action: assign unique speaker roles to file channels and confirm.";
            return "Next action: run Source measurement.";
        }

        if (! targetFamilySelected)
            return "Next action: select a Target profile family.";

        if (! targetLimitsAuthoritative)
            return "Next action: load authoritative Target profile limits.";

        if (! currentStateInitialized)
            return "Next action: initialize Current State from Source.";

        if (currentStateNeedsRemeasurement)
            return "Next action: re-measure Current State.";

        if (! currentStateMeasured)
            return "Next action: measure Current State.";

        if (! proposalEngineAvailable)
            return "Next action: proposal engine remains intentionally unavailable in this scaffold.";

        if (! proposalInputsValidated)
            return "Next action: confirm applicability and valid matching Source/Current State measurements.";

        if (! renderExportAvailable)
            return "Next action: draft a global trim, review it, then Apply explicitly. Export remains unavailable.";

        return "Next action: workflow gates are satisfied.";
    }
};

// [END LS-STANDALONE-WORKFLOW-READINESS]
}

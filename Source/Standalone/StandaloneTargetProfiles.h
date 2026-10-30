#pragma once

#include "Core/Proposal/Model/TargetProfile.h"

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-TARGET-PROFILES]

using TargetProfileLimits = proposal::TargetProfileLimits;
using TargetProfile = proposal::TargetProfile;

class StandaloneTargetProfileCatalog
{
public:
    static constexpr int numProfiles = 4;

    static TargetProfile getProfile (int index)
    {
        TargetProfile profile;

        switch (index)
        {
            case 1:
                return proposal::ebuR128ProgrammeQc();

            case 2:
                profile.id = "atsc-a85-family-placeholder";
                profile.displayName = "ATSC A/85 family";
                profile.familyName = "ATSC A/85";
                profile.description = "Broadcast loudness target family placeholder. Numeric compliance values are intentionally not authoritative in this scaffold.";
                profile.isSelectableForWorkflow = true;
                profile.limits.validationNote = "Profile family selected. Authoritative ATSC A/85 compliance values are not loaded yet.";
                return profile;

            case 3:
                profile.id = "custom-family-placeholder";
                profile.displayName = "Custom target";
                profile.familyName = "Custom";
                profile.description = "Reserved for a future user-defined target profile. Editing custom target values is not implemented yet.";
                profile.isSelectableForWorkflow = true;
                profile.limits.validationNote = "Custom target selected. User-editable compliance values are not implemented yet.";
                return profile;

            case 0:
            default:
                profile.id = "none";
                profile.displayName = "No target selected";
                profile.familyName = "None";
                profile.description = "Source and Current State measurement do not require a delivery target. Select a profile to evaluate its supported measured rules.";
                profile.isSelectableForWorkflow = false;
                profile.limits.validationNote = "No target profile family selected.";
                return profile;
        }
    }

    static juce::String buildLimitsText (const TargetProfile& profile)
    {
        juce::String text;

        text << "Authoritative base profile loaded: ";
        text << (profile.limits.valuesAreAuthoritative ? "yes" : "no");
        text << "\n";

        text << "Integrated loudness target: ";
        if (profile.limits.hasIntegratedLoudnessTarget)
        {
            text << profile.limits.integratedLufs << " LUFS";
            text << " ±" << profile.limits.integratedToleranceLu << " LU";
        }
        else
        {
            text << "--";
        }
        text << "\n";

        text << "Max loudness range: ";
        if (profile.limits.hasMaxLoudnessRange)
            text << profile.limits.maxLoudnessRangeLu << " LU";
        else
            text << "--";
        text << "\n";

        text << "True peak ceiling: ";
        if (profile.limits.hasTruePeakCeiling)
            text << profile.limits.truePeakCeilingDbtp << " dBTP";
        else
            text << "--";
        text << "\n\n";

        text << "Edition: " << profile.standardEdition << "\n"
             << "Definition: " << profile.definitionVersion << "\n"
             << "Applicability: pre-produced whole-programme linear-audio QC only.\n"
             << "The base EBU profile has no LRA ceiling; an explicit job rule is separate.\n"
             << profile.meterToleranceNote << "\n\n"
             << "Profile note: " << profile.limits.validationNote;

        return text;
    }
};

// [END LS-STANDALONE-TARGET-PROFILES]
}

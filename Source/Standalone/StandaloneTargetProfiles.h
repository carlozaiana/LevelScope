#pragma once

#include <juce_core/juce_core.h>

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-TARGET-PROFILES]

struct TargetProfileLimits
{
    bool hasIntegratedLoudnessTarget = false;
    double integratedLufs = 0.0;
    double integratedToleranceLu = 0.0;

    bool hasMaxLoudnessRange = false;
    double maxLoudnessRangeLu = 0.0;

    bool hasTruePeakCeiling = false;
    double truePeakCeilingDbtp = 0.0;

    bool valuesAreAuthoritative = false;

    juce::String validationNote = "No authoritative compliance values loaded.";
};

struct TargetProfile
{
    juce::String id;
    juce::String displayName;
    juce::String familyName;
    juce::String description;

    bool isSelectableForWorkflow = false;

    TargetProfileLimits limits;
};

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
                profile.id = "ebu-r128-family-placeholder";
                profile.displayName = "EBU R128 family";
                profile.familyName = "EBU R128";
                profile.description = "Broadcast loudness target family placeholder. Numeric compliance values are intentionally not authoritative in this scaffold.";
                profile.isSelectableForWorkflow = true;
                profile.limits.validationNote = "Profile family selected. Authoritative EBU R128 compliance values are not loaded yet.";
                return profile;

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
                profile.description = "Choose a target profile family before future measurement, proposal, or export checks are enabled.";
                profile.isSelectableForWorkflow = false;
                profile.limits.validationNote = "No target profile family selected.";
                return profile;
        }
    }

    static juce::String buildLimitsText (const TargetProfile& profile)
    {
        juce::String text;

        text << "Authoritative compliance values loaded: ";
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

        text << "Profile note: " << profile.limits.validationNote;

        return text;
    }
};

// [END LS-STANDALONE-TARGET-PROFILES]
}
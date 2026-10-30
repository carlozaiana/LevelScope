#pragma once
#include <juce_core/juce_core.h>

namespace levelscope::proposal
{
enum class Applicability
{
    unconfirmed, preProducedLinearProgramme, liveProgramme, deliberatelyLowerTarget,
    shortForm, distributionSpecific
};
enum class MeasurementMethod { unsupported, bs1770FullProgramme };

struct TargetProfileLimits
{
    bool hasIntegratedLoudnessTarget = false;
    double integratedLufs = 0, integratedToleranceLu = 0;
    bool hasMaxLoudnessRange = false;
    double maxLoudnessRangeLu = 0;
    bool hasTruePeakCeiling = false;
    double truePeakCeilingDbtp = 0;
    bool valuesAreAuthoritative = false;
    juce::String validationNote = "No authoritative profile loaded.";
};

// Value definition, copied into requests; the catalogue exposes no mutable shared
// object or general editor. Authority is supplied by trusted catalogue code only.
struct TargetProfile
{
    juce::String id, displayName, familyName, description;
    bool isSelectableForWorkflow = false;
    TargetProfileLimits limits;
    int schemaVersion = 1;
    juce::String definitionVersion, standardEdition, referenceUrl;
    juce::String integratedProvenance, truePeakProvenance, lraProvenance;
    MeasurementMethod method = MeasurementMethod::unsupported;
    Applicability applicability = Applicability::unconfirmed;
    // Documentary meter provision, NEVER added to an acceptance ceiling.
    double truePeakMeterToleranceDb = 0, truePeakMeterBandwidthHz = 0;
    juce::String meterToleranceNote;

    // Evaluated view only; the immutable catalogue definition leaves these empty.
    juce::String baseProfileDigest, jobRequirementsDigest, jobRequirementsError;
    juce::String validationError() const;
    juce::String digest() const;
};

TargetProfile ebuR128ProgrammeQc();
}

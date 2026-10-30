#include "TargetProfile.h"
#include <juce_cryptography/juce_cryptography.h>
#include <cmath>
#include <iomanip>
#include <locale>
#include <sstream>

namespace levelscope::proposal
{
juce::String TargetProfile::validationError() const
{
    if (jobRequirementsError.isNotEmpty()) return jobRequirementsError;
    if (! isSelectableForWorkflow || ! limits.valuesAreAuthoritative)
        return "No authoritative profile selected.";
    if (schemaVersion != 1 || id.isEmpty() || definitionVersion.isEmpty()
        || standardEdition.isEmpty() || referenceUrl.isEmpty()
        || integratedProvenance.isEmpty() || truePeakProvenance.isEmpty())
        return "Invalid profile identity or rule provenance.";
    if (method != MeasurementMethod::bs1770FullProgramme
        || applicability != Applicability::preProducedLinearProgramme)
        return "Profile method or applicability is unsupported.";
    const auto& l = limits;
    if (! l.hasIntegratedLoudnessTarget || ! l.hasTruePeakCeiling
        || ! std::isfinite (l.integratedLufs) || ! std::isfinite (l.integratedToleranceLu)
        || l.integratedToleranceLu < 0
        || ! std::isfinite (l.integratedLufs - l.integratedToleranceLu)
        || ! std::isfinite (l.integratedLufs + l.integratedToleranceLu)
        || ! std::isfinite (l.truePeakCeilingDbtp)
        || ! std::isfinite (l.maxLoudnessRangeLu)
        || (l.hasMaxLoudnessRange && (l.maxLoudnessRangeLu < 0 || lraProvenance.isEmpty()))
        || ! std::isfinite (truePeakMeterToleranceDb) || truePeakMeterToleranceDb < 0
        || ! std::isfinite (truePeakMeterBandwidthHz) || truePeakMeterBandwidthHz < 0
        || (truePeakMeterToleranceDb > 0 && (truePeakMeterBandwidthHz == 0 || meterToleranceNote.isEmpty())))
        return "Invalid profile limits or meter-tolerance documentation.";
    return {};
}

juce::String TargetProfile::digest() const
{
    std::ostringstream s; s.imbue (std::locale::classic()); s << std::hexfloat;
    const auto text = [&] (const juce::String& value)
    { const auto bytes = value.toStdString(); s << bytes.size() << ':' << bytes << ';'; };
    s << "target-profile-v1;" << schemaVersion << ';';
    for (const auto& v : { id, definitionVersion, standardEdition, referenceUrl, displayName,
                           familyName, description, integratedProvenance, truePeakProvenance,
                           lraProvenance, meterToleranceNote, limits.validationNote }) text (v);
    s << isSelectableForWorkflow << ';' << static_cast<int> (method) << ';' << static_cast<int> (applicability)
      << ';' << limits.valuesAreAuthoritative << ';' << limits.hasIntegratedLoudnessTarget
      << ';' << limits.integratedLufs << ';' << limits.integratedToleranceLu
      << ';' << limits.hasTruePeakCeiling << ';' << limits.truePeakCeilingDbtp
      << ';' << limits.hasMaxLoudnessRange << ';' << limits.maxLoudnessRangeLu
      << ';' << truePeakMeterToleranceDb << ';' << truePeakMeterBandwidthHz;
    if (jobRequirementsDigest.isNotEmpty() || jobRequirementsError.isNotEmpty())
    { text(baseProfileDigest); text(jobRequirementsDigest); text(jobRequirementsError); }
    const auto bytes = s.str(); return juce::SHA256 (bytes.data(), bytes.size()).toHexString();
}

TargetProfile ebuR128ProgrammeQc()
{
    TargetProfile p;
    p.id = "ebu-r128-2023-programme-linear-qc";
    p.definitionVersion = "1.0.0";
    p.standardEdition = "EBU R128-2023 (November 2023)";
    p.referenceUrl = "https://tech.ebu.ch/docs/r/r128.pdf";
    p.displayName = "EBU R128-2023 programme QC"; p.familyName = "EBU R128";
    p.description = "Selected programme loudness and true-peak rules for pre-produced whole-programme linear-audio QC. "
                    "Excludes live exceptions, deliberate lower targets, short-form supplements and distribution-specific requirements.";
    p.isSelectableForWorkflow = true;
    p.method = MeasurementMethod::bs1770FullProgramme;
    p.applicability = Applicability::preProducedLinearProgramme;
    p.integratedProvenance = "R128-2023 clauses h, i, k, l: -23 LUFS intended target; +/-0.2 LU QC allowance; whole-programme level gating.";
    p.truePeakProvenance = "R128-2023 clause m: maximum -1 dBTP for linear audio.";
    p.lraProvenance = "R128-2023 clause n: informative; no universal LRA ceiling specified.";
    p.truePeakMeterToleranceDb = 0.3; p.truePeakMeterBandwidthHz = 20000;
    p.meterToleranceNote = "Clause m separately permits +/-0.3 dB meter tolerance for signals band-limited to 20 kHz. "
                           "Bandwidth is not certified here; this provision never relaxes the -1 dBTP evaluation ceiling.";
    p.limits.hasIntegratedLoudnessTarget = true;
    p.limits.integratedLufs = -23; p.limits.integratedToleranceLu = 0.2;
    p.limits.hasTruePeakCeiling = true; p.limits.truePeakCeilingDbtp = -1;
    p.limits.valuesAreAuthoritative = true;
    p.limits.validationNote = "Only the selected scalar rules are evaluated. No metadata, complete broadcaster delivery or future encoded-export verification.";
    return p;
}
}

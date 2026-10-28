#include "VerificationEngine.h"
#include <juce_cryptography/juce_cryptography.h>
#include <cmath>
#include <iomanip>
#include <locale>
#include <sstream>

namespace levelscope::proposal
{
const char* verdictName (Verdict v)
{
    switch (v)
    {
        case Verdict::pass: return "PASS";
        case Verdict::fail: return "FAIL";
        case Verdict::indeterminate: return "INDETERMINATE";
        default: return "NOT EVALUATED";
    }
}

Evaluation evaluateMeasuredProgramme (const TargetProfile& profile, const MeasuredProgramme& input)
{
    Evaluation out;
    out.domain = input.domain; out.snapshotIdentity = input.snapshotIdentity;
    out.profileDigest = profile.digest(); out.profileId = profile.id; out.profileVersion = profile.definitionVersion;
    // Include actual evaluated values/status as well as the adapter's full snapshot
    // binding; two different measurements cannot share an evaluation identity.
    std::ostringstream s; s.imbue (std::locale::classic()); s << std::hexfloat;
    s << evaluationVersion << ';' << out.profileDigest.toStdString() << ';' << input.snapshotIdentity.toStdString()
      << ';' << input.measurementVersion.toStdString() << ';' << static_cast<int> (input.domain)
      << ';' << static_cast<int> (input.applicability) << ';' << input.currentIdentityMatches << ';' << input.lraStable;
    for (const auto& m : { input.integrated, input.truePeak, input.lra })
    { s << ';' << static_cast<int> (m.status) << ';' << m.value.has_value(); if (m.value) s << ';' << *m.value; }
    const auto bytes = s.str(); out.identity = juce::SHA256 (bytes.data(), bytes.size()).toHexString();
    out.message = profile.validationError();
    if (out.message.isNotEmpty()) return out;
    if (input.domain != SnapshotDomain::source && input.domain != SnapshotDomain::currentState
        && input.domain != SnapshotDomain::renderedOutput && input.domain != SnapshotDomain::privateCandidate)
    { out.message = "Only explicitly identified measured snapshots are supported; predictions are not measurements."; return out; }
    if (! input.currentIdentityMatches || input.snapshotIdentity.isEmpty())
    { out.message = "No valid current measured snapshot; measure this state first."; return out; }
    if (input.measurementVersion != measurement::version)
    { out.message = "Unsupported measurement method/version."; return out; }
    if (input.applicability != profile.applicability)
    { out.message = "Confirm applicability: pre-produced whole-programme linear-audio QC only."; return out; }

    const auto add = [&] (const char* id, const char* units, const juce::String& provenance,
                          const measurement::Metric& metric, std::optional<double> lower,
                          double upper, double reference, bool stable)
    {
        RuleEvaluation rule; rule.id = id; rule.units = units; rule.provenance = provenance;
        rule.lower = lower; rule.upper = upper;
        if (! metric.value || ! std::isfinite (*metric.value) || metric.status != measurement::Status::available || ! stable)
        {
            rule.verdict = Verdict::indeterminate;
            if (! stable) rule.reason = "Required LRA is unstable for a programme under 60 seconds.";
            else if (metric.value && ! std::isfinite (*metric.value)) rule.reason = "Required metric is non-finite.";
            else if (metric.status == measurement::Status::silence) rule.reason = "Digital silence; required metric unavailable.";
            else if (metric.status == measurement::Status::belowGate) rule.reason = "No eligible loudness blocks.";
            else rule.reason = "Required metric unavailable or inconsistent.";
        }
        else
        {
            rule.measured = metric.value; rule.delta = *metric.value - reference;
            rule.verdict = (! lower || *metric.value >= *lower) && *metric.value <= upper ? Verdict::pass : Verdict::fail;
            rule.reason = rule.verdict == Verdict::pass ? "Within selected bounds." : "Outside selected bounds.";
        }
        out.rules.push_back (std::move (rule));
    };
    const auto& l = profile.limits;
    add ("programme-loudness", "LUFS", profile.integratedProvenance, input.integrated,
         l.integratedLufs - l.integratedToleranceLu, l.integratedLufs + l.integratedToleranceLu, l.integratedLufs, true);
    add ("true-peak", "dBTP", profile.truePeakProvenance, input.truePeak, {}, l.truePeakCeilingDbtp, l.truePeakCeilingDbtp, true);
    if (l.hasMaxLoudnessRange)
        add ("loudness-range", "LU", profile.lraProvenance, input.lra, {}, l.maxLoudnessRangeLu, l.maxLoudnessRangeLu, input.lraStable);
    // A known violation is a fail even if another required metric is unavailable.
    // Otherwise any indeterminate rule prevents an aggregate pass.
    out.verdict = Verdict::pass;
    for (const auto& rule : out.rules)
        if (rule.verdict == Verdict::fail) out.verdict = Verdict::fail;
        else if (rule.verdict != Verdict::pass && out.verdict != Verdict::fail) out.verdict = Verdict::indeterminate;
    out.message = out.verdict == Verdict::pass
        ? (l.hasMaxLoudnessRange ? "Passes selected programme loudness, true-peak and LRA rules."
                               : "Passes selected programme loudness and true-peak rules.")
        : out.verdict == Verdict::fail ? "Fails one or more selected measured rules."
                                     : "Required measurements are unavailable; selected rules are indeterminate.";
    return out;
}
}

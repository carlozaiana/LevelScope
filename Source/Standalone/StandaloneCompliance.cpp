#include "StandaloneCompliance.h"
#include <juce_cryptography/juce_cryptography.h>
#include <iomanip>
#include <locale>
#include <sstream>

namespace levelscope::standalone
{
namespace
{
juce::String snapshotDigest (const AnalysisSnapshot& a)
{
    std::ostringstream s; s.imbue (std::locale::classic()); s << std::hexfloat;
    const auto text = [&] (const juce::String& v)
    { const auto b = v.toStdString(); s << b.size() << ':' << b << ';'; };
    s << "accepted-analysis-v1;";
    for (const auto& v : { a.identity, a.measuredContentHash, a.measurementVersion, a.context.sessionId }) text (v);
    s << a.context.ticket.id << ';' << a.context.ticket.sessionRevision << ';' << static_cast<int> (a.domain);
    const auto optional = [&] (std::optional<double> v) { s << ';' << v.has_value(); if (v) s << ';' << *v; };
    for (const auto& v : { a.measured.integratedLufs, a.measured.truePeakDbtp, a.measured.loudnessRangeLu, a.measured.samplePeakDbfs }) optional (v);
    const auto& d = a.detail;
    s << ';' << d.sampleRate << ';' << d.frames << ';' << d.lraStable << ';' << d.integratedIncluded << ';' << d.lraIncluded;
    for (const auto& m : { d.integrated, d.truePeak, d.lra, d.samplePeak, d.momentaryMax, d.shortTermMax })
    { s << ';' << static_cast<int> (m.status); optional (m.value); }
    for (const auto& v : { d.integratedRelativeGate, d.lraRelativeGate, d.p10, d.p95 }) optional (v);
    for (const auto* windows : { &d.momentary, &d.shortTerm, &d.lraWindows })
    { s << ';' << windows->size(); for (const auto& w : *windows) s << ';' << w.endFrame << ';' << w.energy; }
    for (const auto* peaks : { &d.channelSamplePeaks, &d.channelTruePeaks })
    { s << ';' << peaks->size(); for (auto v : *peaks) s << ';' << v; }
    for (const auto* frames : { &d.samplePeakFrames, &d.truePeakEighthFrames })
    { s << ';' << frames->size(); for (auto v : *frames) s << ';' << v; }
    const auto bytes = s.str(); return juce::SHA256 (bytes.data(), bytes.size()).toHexString();
}
bool sameValue (std::optional<double> a, std::optional<double> b)
{
    return a.has_value() == b.has_value() && (! a || *a == *b);
}
}

proposal::MeasuredProgramme measuredProgrammeForSnapshot (const AnalysisSnapshot* analysis, const ServiceContext& expected,
                                                          MeasurementDomain domain, proposal::Applicability applicability)
{
    proposal::MeasuredProgramme input;
    input.domain = domain == MeasurementDomain::source ? proposal::SnapshotDomain::source
                 : domain == MeasurementDomain::currentState ? proposal::SnapshotDomain::currentState
                                                           : proposal::SnapshotDomain::renderedOutput;
    input.applicability = applicability;
    if (analysis)
    {
        const auto& a = *analysis;
        input.snapshotIdentity = snapshotDigest (a);
        input.measurementVersion = a.measurementVersion;
        input.integrated = a.detail.integrated; input.truePeak = a.detail.truePeak;
        input.lra = a.detail.lra; input.lraStable = a.detail.lraStable;
        input.currentIdentityMatches = a.domain == domain && ! expected.sessionId.isEmpty()
            && a.context.sessionId == expected.sessionId && a.context.isFullProgram()
            && expected.isFullProgram() && a.detail.frames == expected.media.frames
            && a.detail.sampleRate == expected.media.sampleRate
            && OfflineAnalysisService::validateRequest (expected, domain).isEmpty()
            && OfflineAnalysisService::validateRequest (a.context, domain).isEmpty()
            && a.measuredContentHash == expected.media.sha256
            && a.identity == OfflineAnalysisService::identityFor (a.context, domain, a.measuredContentHash)
            && a.identity == OfflineAnalysisService::identityFor (expected, domain, expected.media.sha256)
            && sameValue (a.measured.integratedLufs, a.detail.integrated.value)
            && sameValue (a.measured.truePeakDbtp, a.detail.truePeak.value)
            && sameValue (a.measured.loudnessRangeLu, a.detail.lra.value)
            && sameValue (a.measured.samplePeakDbfs, a.detail.samplePeak.value);
    }
    return input;
}

proposal::Evaluation evaluateSnapshot (const AnalysisSnapshot* analysis, const ServiceContext& expected,
                                      MeasurementDomain domain, const TargetProfile& profile, proposal::Applicability applicability)
{
    return proposal::evaluateMeasuredProgramme (profile, measuredProgrammeForSnapshot (analysis, expected, domain, applicability));
}

juce::String evaluationText (const proposal::Evaluation& e)
{
    juce::String text;
    text << (e.domain == proposal::SnapshotDomain::source ? "Source" : "Current State")
         << " evaluation: " << proposal::verdictName (e.verdict) << "\n" << e.message << "\n";
    for (const auto& rule : e.rules)
    {
        text << rule.id << ": " << proposal::verdictName (rule.verdict);
        // Enough display precision for diagnosis; the decision above uses the
        // original binary64 values even when adjacent values display identically.
        if (rule.measured) text << " | " << juce::String (*rule.measured, 6) << " " << rule.units;
        if (rule.lower) text << " | min " << juce::String (*rule.lower, 6);
        if (rule.upper) text << " | max " << juce::String (*rule.upper, 6);
        if (rule.delta) text << " | delta " << juce::String (*rule.delta, 6);
        text << "\n" << rule.reason << "\n";
    }
    text << "LRA is informative when the profile has no LRA rule. Verdict uses full precision.\n"
         << "No complete delivery, metadata or future encoded-export verification.\n";
    return text;
}
}

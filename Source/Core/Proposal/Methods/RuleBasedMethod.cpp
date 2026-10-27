#include "RuleBasedMethod.h"
#include <juce_cryptography/juce_cryptography.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>

namespace levelscope::proposal
{
TrimProposal draftGlobalTrim (const TrimRequest& r)
{
    TrimProposal p; p.originalGainDb = r.currentGainDb;
    const auto source = evaluateMeasuredProgramme (r.profile, r.source);
    const auto current = evaluateMeasuredProgramme (r.profile, r.current);
    std::ostringstream binding; binding.imbue (std::locale::classic()); binding << std::hexfloat;
    for (const auto& s : { juce::String (trimMethodVersion), source.identity, current.identity,
                           r.sourceInputIdentity, r.currentInputIdentity, r.processingVersion, r.fixedUpstreamPlanIdentity })
    { const auto bytes = s.toStdString(); binding << bytes.size() << ':' << bytes << ';'; }
    binding << r.currentGainDb;
    const auto bytes = binding.str(); p.identity = juce::SHA256 (bytes.data(), bytes.size()).toHexString();
    const auto finiteMetric = [] (const measurement::Metric& m)
    { return m.status == measurement::Status::available && m.value && std::isfinite (*m.value); };
    if (r.source.domain != SnapshotDomain::source || r.current.domain != SnapshotDomain::currentState
        || source.verdict == Verdict::notEvaluated || current.verdict == Verdict::notEvaluated
        || ! finiteMetric (r.source.integrated) || ! finiteMetric (r.source.truePeak)
        || ! finiteMetric (r.current.integrated) || ! finiteMetric (r.current.truePeak))
    { p.explanation = "Valid measured Source and Current State, required metrics and confirmed profile applicability are required."; return p; }
    if (r.sourceInputIdentity.isEmpty() || r.sourceInputIdentity != r.currentInputIdentity)
    { p.explanation = "Source and Current State must describe identical content, channel interpretation and programme range."; return p; }
    const bool supportedPlan=(r.processingVersion=="static-gain-v1" && r.fixedUpstreamPlanIdentity.isEmpty())
        || (r.processingVersion=="section-gain-v1" && r.fixedUpstreamPlanIdentity.isNotEmpty());
    if (!supportedPlan || ! std::isfinite (r.currentGainDb)
        || r.currentGainDb < -60 || r.currentGainDb > 24 || r.profile.limits.hasMaxLoudnessRange)
    { p.explanation = "Global trim requires a supported measured plan, fixed upstream section gains and programme loudness/true-peak rules only."; return p; }

    const double i = *r.current.integrated.value, tp = *r.current.truePeak.value;
    if (current.verdict == Verdict::pass)
    {
        p.status = TrimStatus::noChange; p.totalGainDb = r.currentGainDb; p.deltaDb = 0;
        p.explanation = "Already passes selected measured rules. No change: retain gain, measurements and evaluation.";
        return p;
    }
    const auto& l = r.profile.limits;
    const double desired = l.integratedLufs - i;
    const double lower = std::max (-60 - r.currentGainDb, l.integratedLufs - l.integratedToleranceLu - i);
    const double upper = std::min ({ 24 - r.currentGainDb, l.integratedLufs + l.integratedToleranceLu - i,
                                    l.truePeakCeilingDbtp - tp });
    p.status = TrimStatus::noSolution;
    if (! std::isfinite (desired) || ! std::isfinite (lower) || ! std::isfinite (upper) || lower > upper)
    { p.explanation = "No predicted trim-only solution: loudness allowance, true-peak headroom and supported total gain do not intersect. No plan applied."; return p; }
    double gain = r.currentGainDb + std::clamp (desired, lower, upper);
    // Resolve representable-gain rounding using forward comparisons. This never
    // widens a rule and is not an audio verification/refinement loop.
    for (int step = 0; step < 16; ++step)
    {
        const double delta = gain - r.currentGainDb, predictedI = i + delta, predictedTP = tp + delta;
        const bool tooLow = predictedI < l.integratedLufs - l.integratedToleranceLu || gain < -60;
        const bool tooHigh = predictedI > l.integratedLufs + l.integratedToleranceLu
                             || predictedTP > l.truePeakCeilingDbtp || gain > 24;
        if (! tooLow && ! tooHigh && gain != r.currentGainDb)
        {
            p.status = TrimStatus::draft; p.totalGainDb = gain; p.deltaDb = delta;
            p.predictedIntegrated = predictedI; p.predictedTruePeak = predictedTP;
            p.explanation = delta == desired ? "Global trim toward the profile's intended loudness target."
                : "Global trim constrained by peak headroom or total gain; nearest predicted point to the intended target inside the QC allowance.";
            p.explanation += " Prediction only: absolute-gate membership and sample rounding can change the measured outcome. Apply explicitly, then re-measure; no automatic correction.";
            return p;
        }
        if (tooLow == tooHigh) break;
        gain = std::nextafter (gain, tooLow ? std::numeric_limits<double>::infinity()
                                          : -std::numeric_limits<double>::infinity());
    }
    p.explanation = "No representable trim candidate satisfies the scalar bounds. No plan applied.";
    return p;
}

juce::String trimProposalText (const TrimProposal& p)
{
    const char* label = p.status == TrimStatus::draft ? "DRAFT — prediction only"
                      : p.status == TrimStatus::noChange ? "NO CHANGE"
                      : p.status == TrimStatus::noSolution ? "NO PREDICTED SOLUTION" : "UNAVAILABLE";
    juce::String text; text << "Global-trim proposal: " << label << "\n" << p.explanation << "\n";
    if (p.totalGainDb) text << (p.status == TrimStatus::noChange ? "Retained total gain: " : "Proposed total gain: ")
                           << juce::String (*p.totalGainDb, 12) << " dB\n";
    if (p.deltaDb) text << "Correction from measured Current State: " << juce::String (*p.deltaDb, 12) << " dB\n";
    if (p.predictedIntegrated) text << "Predicted I: " << juce::String (*p.predictedIntegrated, 6) << " LUFS\n";
    if (p.predictedTruePeak) text << "Predicted TP: " << juce::String (*p.predictedTruePeak, 6) << " dBTP\n";
    if (p.status == TrimStatus::draft) text << "Post-change LRA requires measurement. ";
    if (p.status == TrimStatus::noChange) text << "No new measurement is required. ";
    text << "Source audition remains original audio.\n";
    return text;
}
}

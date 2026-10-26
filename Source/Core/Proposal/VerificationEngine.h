#pragma once
#include "Model/TargetProfile.h"
#include "Core/Measurement/OfflineMeasurement.h"
#include <vector>

namespace levelscope::proposal
{
inline constexpr const char* evaluationVersion = "measured-programme-rules-m10.1";
enum class Verdict { notEvaluated, pass, fail, indeterminate };
enum class SnapshotDomain { source, currentState, predicted, renderedOutput };
struct MeasuredProgramme
{
    SnapshotDomain domain = SnapshotDomain::predicted; // Fail closed.
    juce::String snapshotIdentity, measurementVersion;
    Applicability applicability = Applicability::unconfirmed;
    bool currentIdentityMatches = false;
    measurement::Metric integrated, truePeak, lra;
    bool lraStable = false;
};
struct RuleEvaluation
{
    juce::String id, units, provenance, reason;
    Verdict verdict = Verdict::notEvaluated;
    std::optional<double> measured, lower, upper, delta;
};
struct Evaluation
{
    Verdict verdict = Verdict::notEvaluated;
    SnapshotDomain domain = SnapshotDomain::predicted;
    juce::String identity, snapshotIdentity, profileDigest, profileId, profileVersion;
    juce::String message;
    std::vector<RuleEvaluation> rules;
};

// A measured-value evaluator, not an audio processor or exporter. No rounding,
// epsilon, test tolerance, or documentary meter tolerance enters comparisons.
Evaluation evaluateMeasuredProgramme (const TargetProfile&, const MeasuredProgramme&);
const char* verdictName (Verdict);
}

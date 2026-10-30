#pragma once
#include "Core/Proposal/VerificationEngine.h"

namespace levelscope::proposal
{
inline constexpr const char* trimMethodVersion = "global-trim-m16.1";
enum class TrimStatus { unavailable, noChange, draft, noSolution };

// Each snapshot is validated against its expected context by the adapter.
// The input identities additionally bind common content, roles and range.
struct TrimRequest
{
    TargetProfile profile;
    MeasuredProgramme source, current;
    juce::String sourceInputIdentity, currentInputIdentity;
    juce::String processingVersion = "static-gain-v1";
    juce::String fixedUpstreamPlanIdentity; // validated section plan, held fixed
    juce::String preferenceIdentity; // metadata binding; never part of measured audio
    bool preferencesAvailable = true;
    double currentGainDb = 0;
};
struct TrimProposal
{
    TrimStatus status = TrimStatus::unavailable;
    juce::String identity, explanation;
    double originalGainDb = 0;
    std::optional<double> totalGainDb, deltaDb, predictedIntegrated, predictedTruePeak;
};

// Scalar planning only. No audio, session mutation, measured verdict or retry.
TrimProposal draftGlobalTrim (const TrimRequest&);
// Candidate-domain numerical planning only; never a Current State Apply authority.
TrimProposal draftGlobalTrimForCandidate (const TrimRequest&);
juce::String trimProposalText (const TrimProposal&);
}

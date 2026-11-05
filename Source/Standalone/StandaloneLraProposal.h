#pragma once
#include "StandaloneSectionProposal.h"

namespace levelscope::standalone
{
// Planning estimates only; final complete PCM must be independently assessed.
std::optional<double> lraFinalGlobalGain(const AnalysisSnapshot&,const proposal::TargetProfile&);
// Missing/indeterminate required metrics are never competitive. Lexicographic
// ordering is the versioned method contract, not a new compliance verdict.
using LraAttemptRank=std::array<double,7>;
std::optional<LraAttemptRank> lraAttemptRank(const AnalysisSnapshot&,const proposal::Evaluation&,
    int64_t affectedFrames,double maximumAdditional,int assessmentIndex);
}

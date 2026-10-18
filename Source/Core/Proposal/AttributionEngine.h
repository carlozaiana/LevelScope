#pragma once
#include "Model/ProgramAnalysis.h"
#include <functional>

namespace levelscope::proposal
{
// Descriptive overlapping-block shares and percentile membership, not causal
// removal effects, independent scene energy or processing recommendations.
bool buildAttribution (ProgramAnalysis&, const std::function<bool()>& cancelled);
}

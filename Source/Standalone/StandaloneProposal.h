#pragma once
#include "StandaloneCompliance.h"
#include "Core/Proposal/Methods/RuleBasedMethod.h"

namespace levelscope::standalone
{
proposal::TrimRequest trimRequest (const AnalysisSnapshot* source, const AnalysisSnapshot* current,
                                  const ServiceContext& expectedSource, const ServiceContext& expectedCurrent,
                                  const TargetProfile&, proposal::Applicability);
}

#pragma once
#include "StandaloneMeasurementService.h"
#include "Core/Proposal/VerificationEngine.h"

namespace levelscope::standalone
{
// Call only with an accepted immutable snapshot. Worker completion/ticket checks
// belong to the session; profile changes do not invalidate audio measurements.
proposal::MeasuredProgramme measuredProgrammeForSnapshot (const AnalysisSnapshot*, const ServiceContext&,
                                                          MeasurementDomain, proposal::Applicability);
proposal::Evaluation evaluateSnapshot (const AnalysisSnapshot*, const ServiceContext& current,
                                      MeasurementDomain, const TargetProfile&, proposal::Applicability);
juce::String evaluationText (const proposal::Evaluation&);
}

#pragma once
#include "StandaloneMeasurementService.h"
#include "Core/Proposal/VerificationEngine.h"

namespace levelscope::standalone
{
// Call only with an accepted immutable snapshot. Worker completion/ticket checks
// belong to the session; profile changes do not invalidate audio measurements.
bool measurementMatchesContext (const AnalysisSnapshot&, const ServiceContext&, MeasurementDomain);
proposal::MeasuredProgramme measuredProgrammeForSnapshot (const AnalysisSnapshot*, const ServiceContext&,
                                                          MeasurementDomain, proposal::Applicability);
// Published snapshots are immutable. Reuse only the measurement payload digest;
// context validation and applicability still use the caller's current inputs.
proposal::MeasuredProgramme measuredProgrammeForSnapshot (const std::shared_ptr<const AnalysisSnapshot>&,
                                                          const ServiceContext&, MeasurementDomain, proposal::Applicability);
proposal::Evaluation evaluateSnapshot (const AnalysisSnapshot*, const ServiceContext& current,
                                      MeasurementDomain, const TargetProfile&, proposal::Applicability);
juce::String evaluationText (const proposal::Evaluation&);
}

#pragma once
#include "LraTimelineLanes.h"
#include "StandaloneDiagnostics.h"
namespace levelscope::standalone
{
inline void presentLraTimeline(LraTimelineLanes& lanes,const TimelineMeasurement& binding,
    const std::shared_ptr<const MeasuredDiagnostics>& result,const TargetProfile& profile,
    proposal::Applicability applicability,MeasurementDomain domain,const juce::String& diagnosticStatus)
{
    const juce::String name=domain==MeasurementDomain::currentState?"Current":"Original";
    if(binding.snapshot&&result&&result->snapshot==binding.snapshot&&result->lraTimeline)
    {lanes.measured(result->lraTimeline,profile,applicability,name);return;}
    switch(binding.state)
    {
        case TimelineMeasurementState::unavailable:lanes.unavailable(name,"Confirm roles and resolve the programme source before measured inspection.");break;
        case TimelineMeasurementState::unmeasured:lanes.unavailable(name,"No completed measurement. Measure this state first.");break;
        case TimelineMeasurementState::unverified:lanes.unavailable(name,"Needs remeasurement. Previous diagnostics are unavailable for edited audio.");break;
        case TimelineMeasurementState::measuring:lanes.unavailable(name,"Measurement in progress - no completed LRA diagnostics yet.");break;
        case TimelineMeasurementState::measured:lanes.unavailable(name,diagnosticStatus);break;
    }
}
}

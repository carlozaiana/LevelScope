#pragma once
#include "Core/Proposal/RegionPlanning.h"
namespace levelscope::standalone
{
inline juce::String proposalPlanningName(const proposal::RegionPlanningSettings& s)
{
    const auto method=juce::String(s.methodVersion());
    if(method==proposal::innerSearchVersion)return "Competing passage boundaries";
    if(method==proposal::eventSearchVersion)return "Complete-event preference";
    return "Measured-passage refinement";
}
inline juce::String proposalDefinitionText(const proposal::RegionPlanningSettings& s)
{
    return "Stored definition: settings "+juce::String(s.schema)+" / "+s.methodVersion()+" / "+s.formationVersion()
        +(s.schema>=2?" / "+s.eventFeatures.encode()["version"].toString()+" / "+s.coverageVersion()+" / "+s.scheduleVersion():juce::String())
        +". Definition digest "+s.digest()+".";
}
}

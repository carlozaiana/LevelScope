#pragma once
#include "Model/ProgrammeSections.h"

namespace levelscope::proposal
{
// Uses full-resolution M4 windows and original programme gates. Does not decode
// audio, recalculate local I/LRA, evaluate a profile, or infer local peak safety.
std::shared_ptr<const SectionAnalysis> analyzeSections (const ProgramAnalysis&, const SectionMap&,
    const juce::String& expectedInputIdentity, const std::function<bool()>& cancelled = [] { return false; });
}

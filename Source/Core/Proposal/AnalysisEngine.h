#pragma once
#include "Model/ProgramAnalysis.h"
#include <functional>

namespace levelscope::proposal
{
// No profile, audio decoding, processing, GUI, transport or device dependencies.
// Throws on malformed/inconsistent measurement details; cancellation returns null.
std::shared_ptr<const ProgramAnalysis> analyzeProgramme (const measurement::Result&, const juce::String& snapshotIdentity,
                                                       const std::function<bool()>& cancelled = [] { return false; });
TimelinePlot prepareTimeline (const ProgramAnalysis&, double beginFrame, double endFrame, int columns);
// Nearest sampled window end, stable earlier-window tie break; no display-bin data.
std::optional<size_t> inspectWindow (const std::vector<DiagnosticWindow>&, double sourceFrame);
}

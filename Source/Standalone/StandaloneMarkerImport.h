#pragma once
#include "StandaloneSessionModel.h"
namespace levelscope::standalone
{
bool preflightMarkerImport(const StandaloneSessionModel&, const proposal::markers::Preview&, juce::String&);
// The caller supplies the hash of the current file/paste byte snapshot. No audio job is started.
bool commitMarkerImport(StandaloneSessionModel&, const proposal::markers::Preview&,
                        const juce::String& currentBytesHash, juce::String&);
}

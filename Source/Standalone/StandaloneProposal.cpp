#include "StandaloneProposal.h"
#include <juce_cryptography/juce_cryptography.h>
#include <iomanip>
#include <locale>
#include <sstream>

namespace levelscope::standalone
{
namespace
{
juce::String inputIdentity (const AnalysisSnapshot* a)
{
    if (! a || ! a->context.isValid()) return {};
    const auto& c = a->context;
    std::ostringstream s; s.imbue (std::locale::classic()); s << std::hexfloat;
    s << "programme-input-v1;" << c.media.sha256.toStdString() << ';' << a->measuredContentHash.toStdString()
      << ';' << c.media.sampleRate << ';' << c.media.frames << ';' << c.media.channels
      << ';' << c.range.begin << ';' << c.range.end << ';' << c.fileLayout.count << ';' << c.fileLayout.confirmed;
    for (int n = 0; n < c.fileLayout.count; ++n) s << ';' << static_cast<int> (c.fileLayout.roles[n]);
    const auto bytes = s.str(); return juce::SHA256 (bytes.data(), bytes.size()).toHexString();
}
}
proposal::TrimRequest trimRequest (const AnalysisSnapshot* source, const AnalysisSnapshot* current,
                                  const ServiceContext& expectedSource, const ServiceContext& expectedCurrent,
                                  const TargetProfile& profile, proposal::Applicability applicability)
{
    proposal::TrimRequest r; r.profile = profile;
    r.source = measuredProgrammeForSnapshot (source, expectedSource, MeasurementDomain::source, applicability);
    r.current = measuredProgrammeForSnapshot (current, expectedCurrent, MeasurementDomain::currentState, applicability);
    r.sourceInputIdentity = inputIdentity (source); r.currentInputIdentity = inputIdentity (current);
    r.currentGainDb = expectedCurrent.processing.staticGainDb;
    r.processingVersion = expectedCurrent.processing.engineVersion;
    return r;
}
}

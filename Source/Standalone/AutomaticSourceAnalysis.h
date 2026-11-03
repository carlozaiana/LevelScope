#pragma once
#include "StandaloneSessionModel.h"

namespace levelscope::standalone
{
// Message-thread scheduling only. A failed/cancelled attempt is remembered for
// the runtime and exact input identity; unrelated target edits never re-arm it.
class AutomaticSourceAnalysis
{
public:
    bool eligible(const StandaloneSessionModel& session) const
    {
        if(session.source.hasMeasurement || !session.canInitializeCurrentStateFromSource()
            || session.job.isActive() || session.audioEditsBlocked())return false;
        const auto context=session.measurementContext(MeasurementDomain::source);
        return OfflineAnalysisService::validateRequest(context,MeasurementDomain::source).isEmpty()
            && attempted!=key(session,context);
    }
    void attemptedSource(const StandaloneSessionModel& session)
    {attempted=key(session,session.measurementContext(MeasurementDomain::source));}
private:
    juce::String attempted;
    static juce::String key(const StandaloneSessionModel& s,const ServiceContext& c)
    {return s.sessionId+OfflineAnalysisService::identityFor(c,MeasurementDomain::source,c.media.sha256);}
};
}

#pragma once
#include "StandaloneSessionModel.h"
#include "StandaloneAudition.h"

namespace levelscope::standalone
{
// Message-thread boundary for timeline seeking. Inspection never calls this.
struct TimelineTransport
{
    static juce::String blockedReason(const StandaloneSessionModel& session,const StandaloneAudition& audio)
    {
        if(audio.active()||audio.preparing())return "Stop playback and wait for preparation/cancellation before seeking.";
        if(session.persistenceLoading||session.job.isActive())return "Finish or cancel the active job before seeking.";
        if(!session.sourceDocument.hasSource||session.sourceContentInvalid||!session.sourceDocument.matchesFileStamps())return "A current, confirmed programme is required for seeking.";
        return session.auditionRequest().validationError();
    }
    static juce::String seek(int64_t frame,const juce::String& expected,StandaloneSessionModel& session,StandaloneAudition& audio)
    {
        const auto reason=blockedReason(session,audio);if(reason.isNotEmpty())return reason;
        const auto request=session.auditionRequest();
        if(expected.isEmpty()||request.identity()!=expected)return "Audition changed; select a new seek position.";
        if(frame<0||frame>request.context.media.frames)return "Seek is outside the programme.";
        audio.synchronize(request,session.sourceDocument.routing);
        return audio.seek(frame,request,session.sourceDocument.routing);
    }
};
}

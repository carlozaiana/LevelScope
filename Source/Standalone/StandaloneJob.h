#pragma once

#include <atomic>
#include <cstdint>
#include <juce_core/juce_core.h>

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-JOB-CONTRACT]
enum class JobStage { ingest, analysis, proposal, render, report };
enum class JobState { idle, running, cancelling, succeeded, cancelled, failed };
enum class ErrorCode
{
    none, cancelled, unavailable, missingFile, unsupportedMedia, invalidMetadata,
    unsupportedChannels, decodeFailed, nonFiniteAudio, sourceChanged, invalidLayout,
    staleRevision, invalidRequest, deviceUnavailable, sampleRateMismatch
};

struct JobControl
{
    std::atomic<bool> cancel { false };
    std::atomic<float> progress { 0.0f };
};

struct JobTicket
{
    uint64_t id = 0, sessionRevision = 0;
};

struct SessionJob
{
    JobTicket ticket;
    JobStage stage = JobStage::ingest;
    JobState state = JobState::idle;
    ErrorCode error = ErrorCode::none;
    juce::String message = "No job started";

    bool isActive() const { return state == JobState::running || state == JobState::cancelling; }
    bool begin (JobStage nextStage, uint64_t revision)
    {
        if (isActive()) return false;
        ticket = { ticket.id + 1, revision };
        stage = nextStage;
        state = JobState::running;
        error = ErrorCode::none;
        message = "Validating source";
        return true;
    }
    void requestCancel()
    {
        if (state == JobState::running) { state = JobState::cancelling; message = "Cancelling"; }
    }
    bool finish (JobTicket completed, uint64_t currentRevision, ErrorCode result, const juce::String& detail)
    {
        if (! isActive() || completed.id != ticket.id) return false;
        if (state == JobState::cancelling) result = ErrorCode::cancelled;
        else if (completed.sessionRevision != currentRevision) result = ErrorCode::staleRevision;
        error = result;
        state = result == ErrorCode::none ? JobState::succeeded
              : result == ErrorCode::cancelled || result == ErrorCode::staleRevision ? JobState::cancelled
              : JobState::failed;
        message = state == JobState::cancelled ? "Job cancelled; result discarded" : detail;
        return state == JobState::succeeded;
    }
};
// [END LS-STANDALONE-JOB-CONTRACT]
}

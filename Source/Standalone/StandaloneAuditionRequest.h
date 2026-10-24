#pragma once
#include "StandaloneMeasurementService.h"
#include "StandaloneProcessing.h"

namespace levelscope::standalone
{
struct AuditionRequest
{
    ServiceContext context;
    MeasurementDomain domain = MeasurementDomain::source;
    bool currentInitialized = false;

    juce::String validationError() const
    {
        if (! context.isValid() || ! context.isFullProgram() || context.sessionId.isEmpty()
            || ! std::isfinite (context.media.sampleRate) || context.media.sampleRate <= 0
            || context.media.sampleRate != std::floor (context.media.sampleRate))
            return "Validated source, confirmed ordered roles and full programme required";
        if (domain == MeasurementDomain::currentState)
        {
            if (! currentInitialized) return "Initialize Current State first";
            return OfflineAnalysisService::validateRequest (context, domain);
        }
        if (domain != MeasurementDomain::source || ! StaticGainProcessor::supports (context.processing)
            || context.processing.staticGainDb != 0)
            return "Source audition requires identity processing";
        return {};
    }
    juce::String identity() const
    {
        if (validationError().isNotEmpty()) return {};
        return context.sessionId + ":" + OfflineAnalysisService::identityFor (context, domain, context.media.sha256)
            + (domain == MeasurementDomain::currentState ? ":revision=" + juce::String (context.ticket.sessionRevision) : juce::String());
    }
    static AuditionRequest source (const juce::File& file, const MediaInfo& media, const ChannelLayout& layout)
    {
        AuditionRequest r;
        r.context.sessionId = "source-transport";
        r.context.sourceFile = file; r.context.media = media; r.context.fileLayout = layout;
        r.context.range = { 0, media.frames };
        return r;
    }
};
}

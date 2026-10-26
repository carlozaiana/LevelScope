#pragma once
#include "StandaloneServices.h"

namespace levelscope::standalone
{
struct ProcessedAudioObserver;
class OfflineAnalysisService final : public AnalysisService
{
public:
    explicit OfflineAnalysisService (ProcessedAudioObserver* sink = nullptr) : observer (sink) {}
    ServiceResult<AnalysisSnapshot> analyze (const ServiceContext&, MeasurementDomain, JobControl&) override;
    static juce::String validateRequest (const ServiceContext&, MeasurementDomain);
    static juce::String identityFor (const ServiceContext&, MeasurementDomain, const juce::String& measuredHash);
    static std::vector<measurement::Role> measurementRoles (const ChannelLayout&);
    // Test/verification primitive. Encodes to an owned temporary PCM24 WAV, reopens
    // it, and measures those bytes. Does not publish an output or enable export.
    ServiceResult<AnalysisSnapshot> verifyTemporaryRender (const ServiceContext&, JobControl&);
    // Encoded artifact has its own media/layout/hash context and zero gain. It is
    // read in owned export staging; no second full-length copy is manufactured.
    ServiceResult<AnalysisSnapshot> verifyEncodedArtifact (const ServiceContext&, JobControl&);
private:
    ProcessedAudioObserver* observer;
};
}

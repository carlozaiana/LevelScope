#pragma once
#include "StandaloneServices.h"
#include "Core/Waveform/WaveformOverview.h"

namespace levelscope::standalone
{
struct ProcessedAudioObserver;
class VerifiedSource;
class WaveformProgress;
class LoudnessProgress;
class OfflineAnalysisService final : public AnalysisService
{
public:
    explicit OfflineAnalysisService (ProcessedAudioObserver* sink = nullptr, WaveformProgress* preview = nullptr, LoudnessProgress* loudness = nullptr)
        : observer (sink), waveformProgress(preview), loudnessProgress(loudness) {}
    ServiceResult<AnalysisSnapshot> analyze (const ServiceContext&, MeasurementDomain, JobControl&) override;
    // Replaceable runtime overview only. Does not measure or modify a valid snapshot.
    ServiceResult<std::shared_ptr<const waveform::Overview>> buildWaveform(const ServiceContext&,MeasurementDomain,JobControl&);
    static juce::String validateRequest (const ServiceContext&, MeasurementDomain);
    static bool hasCompleteCoverage (const AnalysisSnapshot&);
    static juce::String identityFor (const ServiceContext&, MeasurementDomain, const juce::String& measuredHash);
    static std::vector<measurement::Role> measurementRoles (const ChannelLayout&);
    // Test/verification primitive. Encodes to an owned temporary PCM24 WAV, reopens
    // it, and measures those bytes. Does not publish an output or enable export.
    ServiceResult<AnalysisSnapshot> verifyTemporaryRender (const ServiceContext&, JobControl&);
    // Encoded artifact has its own media/layout/hash context and zero gain. It is
    // read in owned export staging; no second full-length copy is manufactured.
    ServiceResult<AnalysisSnapshot> verifyEncodedArtifact (const ServiceContext&, JobControl&);
    // Private assessment only, reusing one verified original across bounded passes.
    // Result domain cannot satisfy Current freshness or artifact approval.
    ServiceResult<AnalysisSnapshot> analyzeCandidate(const ServiceContext&, const VerifiedSource&,
        JobControl&, measurement::PeakObserver* = nullptr);
private:
    ProcessedAudioObserver* observer;
    WaveformProgress* waveformProgress;
    LoudnessProgress* loudnessProgress;
};
}

#include "StandaloneMeasurementService.h"
#include "StandaloneProcessing.h"
#include "StandaloneMeasurementLayout.h"
#include "VerifiedSource.h"
#include "Pcm24Wave.h"
#include <juce_cryptography/juce_cryptography.h>
#include <cmath>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace levelscope::standalone
{
namespace
{
struct Failure { ErrorCode code; juce::String message; };
void checkCancel (JobControl& control)
{
    if (control.isCancelled()) throw Failure { ErrorCode::cancelled, "Measurement cancelled" };
}
struct OwnedFile
{
    juce::File file = juce::File::getSpecialLocation (juce::File::tempDirectory)
        .getChildFile ("levelscope-" + juce::Uuid().toString() + ".wav");
    ~OwnedFile() { file.setReadOnly (false); file.deleteFile(); }
};

AnalysisSnapshot measure (const ServiceContext& context, MeasurementDomain domain,
                          std::unique_ptr<ProgrammeReader> reader, const juce::String& hash, double gainDb, JobControl& control, ProcessedAudioObserver* observer)
{
    if (!reader) throw Failure { ErrorCode::decodeFailed, "Cannot decode complete verified programme" };
    measurement::OfflineMeasurement engine (static_cast<int> (reader->sampleRate), OfflineAnalysisService::measurementRoles (context.fileLayout));
    constexpr int blockSize = 4096;
    juce::AudioBuffer<float> audio (context.media.channels, blockSize);
    const StaticGainProcessor processing (gainDb);
    for (juce::int64 frame = 0; frame < reader->lengthInSamples;)
    {
        checkCancel (control);
        const int count = static_cast<int> (std::min<juce::int64> (blockSize, reader->lengthInSamples - frame));
        if (! reader->read (audio.getArrayOfWritePointers(), audio.getNumChannels(), frame, count))
            throw Failure { ErrorCode::decodeFailed, "Measurement decode failed" };
        if (! processing.process (audio.getArrayOfWritePointers(), audio.getNumChannels(), count))
            throw Failure { ErrorCode::nonFiniteAudio, "Non-finite processed audio" };
        if (observer) observer->observe (audio.getArrayOfReadPointers(), audio.getNumChannels(), count, frame);
        engine.push (audio.getArrayOfReadPointers(), count);
        frame += count;
        control.progress.store (0.5f + 0.49f * static_cast<float> (static_cast<double> (frame) / reader->lengthInSamples));
    }
    AnalysisSnapshot result;
    result.context = context; result.domain = domain; result.measuredContentHash = hash;
    result.measurementVersion = measurement::version;
    result.identity = OfflineAnalysisService::identityFor (context, domain, hash);
    result.detail = engine.finish();
    if (context.assembly && !context.sourceFilesMatch())
        throw Failure { ErrorCode::sourceChanged, "A programme member changed during measurement; result discarded" };
    checkCancel (control);
    result.measured = { result.detail.integrated.value, result.detail.lra.value,
                        result.detail.truePeak.value, result.detail.samplePeak.value };
    control.progress.store (1);
    return result;
}
template <typename Fn> ServiceResult<AnalysisSnapshot> run (Fn&& function)
{
    try { return { function(), ErrorCode::none, "Measurement complete; delivery compliance not evaluated" }; }
    catch (const Failure& error) { return { {}, error.code, error.message }; }
    catch (const std::exception& error) { return { {}, ErrorCode::invalidRequest, error.what() }; }
}
}

std::vector<measurement::Role> OfflineAnalysisService::measurementRoles (const ChannelLayout& layout)
{
    using R = measurement::Role;
    using C = juce::AudioChannelSet;
    std::vector<R> roles;
    if (! layout.isValid()) return roles;
    if (layout.count == 8 && ! isSystemI (layout)) return roles;
    const auto has = [&] (ChannelRole role) { return std::find (layout.roles.begin(), layout.roles.begin() + layout.count, role) != layout.roles.begin() + layout.count; };
    if (has (C::leftSurround) != has (C::rightSurround)
        || has (C::leftSurroundSide) != has (C::rightSurroundSide)) return roles;
    for (int i = 0; i < layout.count; ++i)
    {
        switch (layout.roles[i])
        {
            case C::left: roles.push_back (R::left); break;
            case C::right: roles.push_back (R::right); break;
            case C::centre: roles.push_back (R::centre); break;
            case C::LFE: roles.push_back (R::lfe); break;
            case C::leftSurround: case C::leftSurroundSide: roles.push_back (R::surroundLeft); break;
            case C::rightSurround: case C::rightSurroundSide: roles.push_back (R::surroundRight); break;
            case C::leftSurroundRear: if (! isSystemI (layout)) return {}; roles.push_back (R::rearLeft); break;
            case C::rightSurroundRear: if (! isSystemI (layout)) return {}; roles.push_back (R::rearRight); break;
            default: return {};
        }
    }
    return roles;
}
juce::String OfflineAnalysisService::validateRequest (const ServiceContext& c, MeasurementDomain domain)
{
    if (! c.isValid() || ! c.isFullProgram() || c.sessionId.isEmpty()) return "Full program, session identity and confirmed roles required";
    if (c.media.measurementRestriction.isNotEmpty()) return c.media.measurementRestriction;
    if (! std::isfinite (c.media.sampleRate) || c.media.sampleRate < 1 || c.media.sampleRate > 192000
        || c.media.sampleRate != std::floor (c.media.sampleRate)
        || ! measurement::OfflineMeasurement::supports (static_cast<int> (c.media.sampleRate), measurementRoles (c.fileLayout)))
        return "Measurement supports mono, stereo, 5.0/5.1 and explicit BS.2051 System I 7.1 at 44.1/48/96/192 kHz";
    if (! StaticGainProcessor::supports (c.processing))
        return "Unsupported processing plan; M1 accepts static gain from -60 to +24 dB";
    if (domain == MeasurementDomain::source && c.processing.staticGainDb != 0)
        return "Source measurement requires identity processing";
    return {};
}
juce::String OfflineAnalysisService::identityFor (const ServiceContext& c, MeasurementDomain domain, const juce::String& hash)
{
    std::ostringstream key; key.imbue (std::locale::classic());
    if (c.assembly) key << "mono-channel-assembly-v2\n" << c.assembly->definitionIdentity() << '\n';
    key << measurement::version << '\n' << layoutMappingVersion << '\n'
        << (isSystemI (c.fileLayout) ? systemICapability : "m1-layouts-v1") << '\n'
        << c.media.sha256 << '\n' << hash << '\n'
        << static_cast<int> (domain) << '\n' << std::hexfloat << c.media.sampleRate << '\n'
        << c.range.begin << ':' << c.range.end << '\n' << c.processing.schemaVersion << '\n'
        << c.processing.engineVersion << '\n' << c.processing.staticGainDb << '\n'
        << c.buildVersion << '\n' << c.platformAndFloatingPointMode << "\njuce-7.0.12-float32-decode\n";
    for (int i = 0; i < c.fileLayout.count; ++i) key << static_cast<int> (c.fileLayout.roles[i]) << ',';
    const auto text = key.str();
    return juce::SHA256 (text.data(), text.size()).toHexString();
}
ServiceResult<AnalysisSnapshot> OfflineAnalysisService::analyze (const ServiceContext& context, MeasurementDomain domain, JobControl& control)
{
    return run ([&]
    {
        const auto error = validateRequest (context, domain);
        if (error.isNotEmpty() || domain == MeasurementDomain::renderedOutput)
            throw Failure { ErrorCode::invalidRequest, error.isNotEmpty() ? error : "Rendered output must be measured through artifact verification" };
        const auto prepared = prepareVerifiedSource (context, control);
        if (! prepared.value) throw Failure { prepared.error, prepared.message };
        const auto& source = *prepared.value;
        return measure (context, domain, source->createReader(), context.media.sha256, context.processing.staticGainDb, control, observer);
    });
}
bool ExportApproval::matches (const ServiceContext& context, const RenderArtifact& artifact) const
{
    const auto& measured = artifact.verifiedOutput;
    const auto& a = context.profile.limits;
    const auto& b = measured.context.profile.limits;
    return ! outputHash.isEmpty() && ! profileVersion.isEmpty()
        && context.profile.limits.valuesAreAuthoritative
        && a.hasIntegratedLoudnessTarget == b.hasIntegratedLoudnessTarget && a.integratedLufs == b.integratedLufs
        && a.integratedToleranceLu == b.integratedToleranceLu && a.hasMaxLoudnessRange == b.hasMaxLoudnessRange
        && a.maxLoudnessRangeLu == b.maxLoudnessRangeLu && a.hasTruePeakCeiling == b.hasTruePeakCeiling
        && a.truePeakCeilingDbtp == b.truePeakCeilingDbtp && a.valuesAreAuthoritative == b.valuesAreAuthoritative
        && OfflineAnalysisService::validateRequest (context, MeasurementDomain::currentState).isEmpty()
        && measured.context.isFullProgram() && measured.domain == MeasurementDomain::renderedOutput
        && measured.measurementVersion == measurement::version
        && measured.measured.integratedLufs && measured.measured.truePeakDbtp
        && analysisIdentity == OfflineAnalysisService::identityFor (context, MeasurementDomain::renderedOutput, outputHash)
        && analysisIdentity == measured.identity
        && outputHash == artifact.sha256 && outputHash == measured.measuredContentHash
        && ! sessionId.isEmpty() && sessionId == context.sessionId && sessionId == measured.context.sessionId
        && revision == context.ticket.sessionRevision && revision == measured.context.ticket.sessionRevision
        && profileId == context.profile.id && profileVersion == context.profileVersion
        && profileId == measured.context.profile.id && profileVersion == measured.context.profileVersion;
}
ServiceResult<AnalysisSnapshot> OfflineAnalysisService::verifyTemporaryRender (const ServiceContext& context, JobControl& control)
{
    return run ([&]
    {
        const auto error = validateRequest (context, MeasurementDomain::currentState);
        if (error.isNotEmpty()) throw Failure { ErrorCode::invalidRequest, error };
        const auto prepared = prepareVerifiedSource (context, control);
        if (! prepared.value) throw Failure { prepared.error, prepared.message };
        const auto& source = *prepared.value;
        auto reader = source->createReader();
        if (!reader) throw Failure { ErrorCode::decodeFailed, "Cannot decode verified programme" };
        OwnedFile output;
        auto stream = output.file.createOutputStream();
        if (! stream) throw Failure { ErrorCode::unavailable, "Cannot create private render" };
        juce::WavAudioFormat wav;
        std::unique_ptr<juce::AudioFormatWriter> writer (wav.createWriterFor (stream.get(), reader->sampleRate, reader->numChannels, 24, {}, 0));
        if (! writer) throw Failure { ErrorCode::unavailable, "Cannot create PCM24 encoder" };
        stream.release();
        juce::AudioBuffer<float> audio (context.media.channels, 4096);
        const StaticGainProcessor processing (context.processing.staticGainDb);
        for (juce::int64 frame = 0; frame < reader->lengthInSamples;)
        {
            checkCancel (control);
            const int count = static_cast<int> (std::min<juce::int64> (4096, reader->lengthInSamples-frame));
            if (! reader->read (audio.getArrayOfWritePointers(), audio.getNumChannels(), frame, count))
                throw Failure { ErrorCode::decodeFailed, "Render decode failed" };
            if (! processing.process (audio.getArrayOfWritePointers(), audio.getNumChannels(), count))
                throw Failure { ErrorCode::nonFiniteAudio, "Non-finite processed audio" };
            for (int channel = 0; channel < audio.getNumChannels(); ++channel)
                for (int i = 0; i < count; ++i)
                {
                    const auto x = audio.getSample (channel, i);
                    if (! std::isfinite (x) || std::abs (x) >= 1)
                        throw Failure { ErrorCode::invalidRequest, "PCM24 render would clip" };
                    audio.setSample (channel, i, x);
                }
            if (! writer->writeFromAudioSampleBuffer (audio, 0, count)) throw Failure { ErrorCode::decodeFailed, "Render write failed" };
            frame += count;
        }
        writer.reset();
        JobControl verificationControl; verificationControl.parentCancel = &control.cancel;
        const auto validation = validateMedia (output.file, verificationControl);
        checkCancel (control);
        if (! validation.succeeded() || validation.media.frames != context.media.frames
            || validation.media.channels != context.media.channels || validation.media.sampleRate != context.media.sampleRate
            || validation.media.bitsPerSample != 24)
            throw Failure { ErrorCode::decodeFailed, "Rendered artifact validation failed" };
        return measure (context, MeasurementDomain::renderedOutput, ProgrammeReader::open({output.file}, validation.media), validation.media.sha256, 0, control, observer);
    });
}
ServiceResult<AnalysisSnapshot> OfflineAnalysisService::verifyEncodedArtifact (const ServiceContext& context, JobControl& control)
{
    return run ([&]
    {
        const auto error = validateRequest (context, MeasurementDomain::renderedOutput);
        if (error.isNotEmpty() || context.assembly || context.processing.staticGainDb != 0)
            throw Failure { ErrorCode::invalidRequest, "Artifact measurement requires its own supported PCM context and zero additional gain. " + error };
        pcm24::Format format; juce::String formatError;
        if (!pcm24::inspect (context.sourceFile, format, formatError, control))
            throw Failure { ErrorCode::invalidMetadata, formatError };
        const auto validation = validateMedia (context.sourceFile, control);
        if (!validation.succeeded() || validation.media.sha256 != context.media.sha256
            || format.frames != context.media.frames || format.rate != context.media.sampleRate
            || format.channels != context.media.channels)
            throw Failure { ErrorCode::sourceChanged, "Encoded artifact bytes or format changed" };
        auto result = measure (context, MeasurementDomain::renderedOutput,
            ProgrammeReader::open ({context.sourceFile}, context.media), context.media.sha256, 0, control, observer);
        if (pcm24::hashFile (context.sourceFile, control) != context.media.sha256)
            throw Failure { ErrorCode::sourceChanged, "Artifact changed during measurement" };
        return result;
    });
}
}

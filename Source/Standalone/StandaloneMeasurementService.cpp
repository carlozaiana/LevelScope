#include "StandaloneMeasurementService.h"
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

// Hashing is cancellable; never call SHA256(File) on a potentially large source.
class HashInput final : public juce::InputStream
{
public:
    HashInput (juce::FileInputStream& s, JobControl& c) : stream (s), control (c) {}
    juce::int64 getTotalLength() override { return stream.getTotalLength(); }
    juce::int64 getPosition() override { return stream.getPosition(); }
    bool setPosition (juce::int64 p) override { return stream.setPosition (p); }
    bool isExhausted() override { return control.isCancelled() || stream.isExhausted(); }
    int read (void* dst, int n) override { checkCancel (control); return stream.read (dst, n); }
private:
    juce::FileInputStream& stream;
    JobControl& control;
};
juce::String hashFile (const juce::File& file, JobControl& control)
{
    auto stream = file.createInputStream();
    if (! stream) throw Failure { ErrorCode::missingFile, "Cannot open measurement input" };
    HashInput input (*stream, control);
    const auto hash = juce::SHA256 (input).toHexString();
    checkCancel (control);
    if (stream->getStatus().failed() || stream->getPosition() != stream->getTotalLength())
        throw Failure { ErrorCode::decodeFailed, "Incomplete content hash" };
    return hash;
}
void snapshot (const ServiceContext& context, OwnedFile& destination, JobControl& control)
{
    auto input = context.sourceFile.createInputStream();
    if (! input) throw Failure { ErrorCode::missingFile, "Source is unavailable" };
    if (input->getTotalLength() != context.media.bytes)
        throw Failure { ErrorCode::sourceChanged, "Source length changed; import again" };
    auto output = destination.file.createOutputStream();
    if (! output) throw Failure { ErrorCode::unavailable, "Cannot create private measurement snapshot" };
    std::array<char, 65536> buffer;
    juce::int64 copied = 0;
    while (copied < context.media.bytes)
    {
        checkCancel (control);
        const int wanted = static_cast<int> (std::min<juce::int64> (buffer.size(), context.media.bytes - copied));
        const int n = input->read (buffer.data(), wanted);
        if (n != wanted || ! output->write (buffer.data(), static_cast<size_t> (n)))
            throw Failure { ErrorCode::decodeFailed, "Snapshot read/write failed" };
        copied += n;
        control.progress.store (0.15f * static_cast<float> (static_cast<double> (copied) / context.media.bytes));
    }
    output->flush();
    if (output->getStatus().failed()) throw Failure { ErrorCode::decodeFailed, "Snapshot flush failed" };
    output.reset();
    if (! input->isExhausted() || input->getTotalLength() != context.media.bytes)
        throw Failure { ErrorCode::sourceChanged, "Source length changed during snapshot" };
    if (hashFile (destination.file, control) != context.media.sha256)
        throw Failure { ErrorCode::sourceChanged, "Source content changed; import again" };
    destination.file.setReadOnly (true);
    // Exact bytes match an imported file, but repeat payload validation so decoder
    // short reads/checksum errors cannot become measurements through zero padding.
    JobControl validation;
    validation.parentCancel = &control.cancel;
    const auto verified = validateMedia (destination.file, validation);
    checkCancel (control);
    if (! verified.succeeded()) throw Failure { verified.error, verified.message };
    if (verified.media.sha256 != context.media.sha256 || verified.media.frames != context.media.frames
        || verified.media.channels != context.media.channels || verified.media.sampleRate != context.media.sampleRate)
        throw Failure { ErrorCode::sourceChanged, "Snapshot metadata differs from imported source" };
    control.progress.store (0.5f);
}
std::unique_ptr<juce::AudioFormatReader> readerFor (const juce::File& file)
{
    juce::AudioFormatManager formats; registerSourceFormats (formats);
    // Snapshot names are private implementation details. Probe bytes, not suffix:
    // AIFF/FLAC snapshots deliberately share the same temporary naming convention.
    auto reader = std::unique_ptr<juce::AudioFormatReader> (formats.createReaderFor (file.createInputStream()));
    if (! reader) throw Failure { ErrorCode::decodeFailed, "Cannot decode verified snapshot" };
    return reader;
}
AnalysisSnapshot measure (const ServiceContext& context, MeasurementDomain domain,
                          const juce::File& file, const juce::String& hash, double gainDb, JobControl& control)
{
    auto reader = readerFor (file);
    measurement::OfflineMeasurement engine (static_cast<int> (reader->sampleRate), OfflineAnalysisService::measurementRoles (context.fileLayout));
    constexpr int blockSize = 4096;
    juce::AudioBuffer<float> audio (context.media.channels, blockSize);
    const double gain = std::pow (10.0, gainDb / 20);
    for (juce::int64 frame = 0; frame < reader->lengthInSamples;)
    {
        checkCancel (control);
        const int count = static_cast<int> (std::min<juce::int64> (blockSize, reader->lengthInSamples - frame));
        if (! reader->read (audio.getArrayOfWritePointers(), audio.getNumChannels(), frame, count))
            throw Failure { ErrorCode::decodeFailed, "Measurement decode failed" };
        for (int channel = 0; channel < audio.getNumChannels(); ++channel)
            for (int i = 0; i < count; ++i)
            {
                const float x = static_cast<float> (audio.getSample (channel, i) * gain);
                if (! std::isfinite (x)) throw Failure { ErrorCode::nonFiniteAudio, "Non-finite processed audio" };
                audio.setSample (channel, i, x);
            }
        engine.push (audio.getArrayOfReadPointers(), count);
        frame += count;
        control.progress.store (0.5f + 0.49f * static_cast<float> (static_cast<double> (frame) / reader->lengthInSamples));
    }
    AnalysisSnapshot result;
    result.context = context; result.domain = domain; result.measuredContentHash = hash;
    result.measurementVersion = measurement::version;
    result.identity = OfflineAnalysisService::identityFor (context, domain, hash);
    result.detail = engine.finish();
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
            default: return {};
        }
    }
    return roles;
}
juce::String OfflineAnalysisService::validateRequest (const ServiceContext& c, MeasurementDomain domain)
{
    if (! c.isValid() || ! c.isFullProgram() || c.sessionId.isEmpty()) return "Full program, session identity and confirmed roles required";
    if (! std::isfinite (c.media.sampleRate) || c.media.sampleRate < 1 || c.media.sampleRate > 192000
        || c.media.sampleRate != std::floor (c.media.sampleRate)
        || ! measurement::OfflineMeasurement::supports (static_cast<int> (c.media.sampleRate), measurementRoles (c.fileLayout)))
        return "Measurement supports mono, stereo, 5.0/5.1 at 44.1/48/96/192 kHz";
    if (c.processing.schemaVersion != 1 || c.processing.engineVersion != "static-gain-v1"
        || c.processing.graphJson.isNotEmpty() || ! std::isfinite (c.processing.staticGainDb)
        || c.processing.staticGainDb < -60 || c.processing.staticGainDb > 24)
        return "Unsupported processing plan; M1 accepts static gain from -60 to +24 dB";
    if (domain == MeasurementDomain::source && c.processing.staticGainDb != 0)
        return "Source measurement requires identity processing";
    return {};
}
juce::String OfflineAnalysisService::identityFor (const ServiceContext& c, MeasurementDomain domain, const juce::String& hash)
{
    std::ostringstream key; key.imbue (std::locale::classic());
    key << measurement::version << '\n' << c.media.sha256 << '\n' << hash << '\n'
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
        OwnedFile source; snapshot (context, source, control);
        return measure (context, domain, source.file, context.media.sha256, context.processing.staticGainDb, control);
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
        OwnedFile source; snapshot (context, source, control);
        auto reader = readerFor (source.file);
        OwnedFile output;
        auto stream = output.file.createOutputStream();
        if (! stream) throw Failure { ErrorCode::unavailable, "Cannot create private render" };
        juce::WavAudioFormat wav;
        std::unique_ptr<juce::AudioFormatWriter> writer (wav.createWriterFor (stream.get(), reader->sampleRate, reader->numChannels, 24, {}, 0));
        if (! writer) throw Failure { ErrorCode::unavailable, "Cannot create PCM24 encoder" };
        stream.release();
        juce::AudioBuffer<float> audio (context.media.channels, 4096);
        const double gain = std::pow (10.0, context.processing.staticGainDb / 20);
        for (juce::int64 frame = 0; frame < reader->lengthInSamples;)
        {
            checkCancel (control);
            const int count = static_cast<int> (std::min<juce::int64> (4096, reader->lengthInSamples-frame));
            if (! reader->read (audio.getArrayOfWritePointers(), audio.getNumChannels(), frame, count))
                throw Failure { ErrorCode::decodeFailed, "Render decode failed" };
            for (int channel = 0; channel < audio.getNumChannels(); ++channel)
                for (int i = 0; i < count; ++i)
                {
                    const auto x = static_cast<float> (audio.getSample (channel, i) * gain);
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
        return measure (context, MeasurementDomain::renderedOutput, output.file, validation.media.sha256, 0, control);
    });
}
}

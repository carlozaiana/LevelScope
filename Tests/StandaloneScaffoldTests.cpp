#include "Standalone/StandaloneAudition.h"
#include "Standalone/StandaloneSessionModel.h"
#include "Standalone/StandaloneServices.h"
#include <juce_cryptography/juce_cryptography.h>
#include <iostream>
#include <limits>
#include <thread>

using namespace levelscope::standalone;

namespace
{
int checks = 0;
void expect (bool value, const char* description)
{
    ++checks;
    if (! value) throw std::runtime_error (description);
}

struct TempDirectory
{
    juce::File path = juce::File::getCurrentWorkingDirectory()
        .getNonexistentChildFile ("levelscope-fixtures", {}, false);
    TempDirectory() { expect (path.createDirectory().wasOk(), "create fixture directory"); }
    ~TempDirectory() { path.deleteRecursively(); }
};

juce::File writeFixture (const juce::File& directory, const juce::String& name, int channels,
                         int frames = 9000, double rate = 48000, int format = 0)
{
    const auto file = directory.getChildFile (name);
    juce::WavAudioFormat wav;
    juce::AiffAudioFormat aiff;
    juce::FlacAudioFormat flac;
    auto* codec = format == 1 ? static_cast<juce::AudioFormat*> (&aiff)
                : format == 2 ? static_cast<juce::AudioFormat*> (&flac) : &wav;
    auto stream = file.createOutputStream();
    expect (stream != nullptr, "create fixture output");
    std::unique_ptr<juce::AudioFormatWriter> writer (
        codec->createWriterFor (stream.get(), rate, static_cast<unsigned int> (channels), 24, {}, 0));
    expect (writer != nullptr, "create fixture writer");
    stream.release();
    juce::AudioBuffer<float> samples (channels, frames);
    for (int channel = 0; channel < channels; ++channel)
        for (int frame = 0; frame < frames; ++frame)
            samples.setSample (channel, frame, static_cast<float> (channel + 1) * 0.01f);
    expect (writer->writeFromAudioSampleBuffer (samples, 0, frames), "write fixture audio");
    return file;
}

void mediaChecks (const juce::File& directory)
{
    for (int channels : { 1, 2, 6, 12 })
    {
        const auto file = writeFixture (directory, juce::String (channels) + ".wav", channels);
        const auto hashBefore = juce::SHA256 (file).toHexString();
        JobControl control;
        const auto result = validateMedia (file, control);
        expect (result.succeeded(), "valid WAV accepted");
        expect (result.media.channels == channels && result.media.frames == 9000, "decoded channel/frame counts");
        expect (result.media.sampleRate == 48000 && result.media.bitsPerSample == 24, "decoded rate/depth");
        expect (result.media.sha256 == hashBefore && juce::SHA256 (file).toHexString() == hashBefore, "source preserved and hash captured");
        expect (! result.media.suggestedLayout.confirmed, "decoder roles require confirmation");
        expect (control.progress.load() == 1.0f, "validation completes progress");
    }
    for (int format : { 1, 2 })
    {
        const auto file = writeFixture (directory, format == 1 ? "source.aiff" : "source.flac", 2, 2000, 44100, format);
        JobControl control;
        const auto result = validateMedia (file, control);
        expect (result.succeeded() && result.media.sampleRate == 44100, "AIFF/FLAC decoded");
    }
    JobControl control;
    const auto wide = writeFixture (directory, "too-wide.wav", 13);
    expect (validateMedia (wide, control).error == ErrorCode::unsupportedChannels, "reject 13 channels");
    const auto corrupt = directory.getChildFile ("corrupt.wav");
    corrupt.replaceWithText ("this is not audio");
    expect (validateMedia (corrupt, control).error == ErrorCode::unsupportedMedia, "reject fake WAV extension");
    expect (validateMedia (directory.getChildFile ("missing.wav"), control).error == ErrorCode::missingFile, "missing file failure");
    const auto empty = writeFixture (directory, "empty.wav", 2, 0);
    expect (! validateMedia (empty, control).succeeded(), "reject zero-frame audio");
    const auto renamed = directory.getChildFile ("audio.unknown");
    expect (directory.getChildFile ("2.wav").copyFileTo (renamed), "copy extension-independent fixture");
    auto result = validateMedia (renamed, control);
    expect (result.succeeded(), "probe content independent of extension");
    const auto truncated = directory.getChildFile ("truncated.wav");
    juce::MemoryBlock bytes;
    expect (renamed.loadFileAsData (bytes), "read fixture bytes");
    expect (truncated.replaceWithData (bytes.getData(), bytes.getSize() - 100), "truncate payload fixture");
    expect (validateMedia (truncated, control).error == ErrorCode::decodeFailed, "reject truncated PCM, no silent padding");
    const auto brokenFlac = directory.getChildFile ("truncated.flac");
    expect (directory.getChildFile ("source.flac").loadFileAsData (bytes), "read FLAC bytes");
    expect (brokenFlac.replaceWithData (bytes.getData(), bytes.getSize() / 2), "truncate FLAC fixture");
    expect (! validateMedia (brokenFlac, control).succeeded(), "reject truncated FLAC");
    control.cancel.store (true);
    expect (validateMedia (renamed, control).error == ErrorCode::cancelled, "cancel import");
    renamed.appendText ("changed");
    expect (! result.media.matchesFileStamp (renamed), "changed source rejected by stamp");

    const auto invalidFloat = directory.getChildFile ("nonfinite.wav");
    {
        juce::WavAudioFormat wav;
        auto stream = invalidFloat.createOutputStream();
        std::unique_ptr<juce::AudioFormatWriter> writer (wav.createWriterFor (stream.get(), 48000, 1, 32, {}, 0));
        expect (writer != nullptr, "create float fixture");
        stream.release();
        juce::AudioBuffer<float> samples (1, 64);
        samples.clear();
        samples.setSample (0, 32, std::numeric_limits<float>::quiet_NaN());
        expect (writer->writeFromAudioSampleBuffer (samples, 0, 64), "write nonfinite fixture");
    }
    JobControl finiteControl;
    expect (validateMedia (invalidFloat, finiteControl).error == ErrorCode::nonFiniteAudio, "reject nonfinite decoded audio");
    const auto longFile = writeFixture (directory, "cancel.wav", 2, 2000000);
    JobControl runningControl;
    MediaValidationResult cancelledResult;
    std::atomic<bool> finished { false };
    std::thread validator ([&] { cancelledResult = validateMedia (longFile, runningControl); finished.store (true); });
    while (runningControl.progress.load() == 0 && ! finished.load()) std::this_thread::yield();
    runningControl.cancel.store (true);
    validator.join();
    expect (cancelledResult.error == ErrorCode::cancelled, "cancel an in-flight decode job");
}

void routingChecks()
{
    ChannelLayout layout;
    layout.count = 2;
    expect (! layout.isValid(), "unknown roles invalid");
    layout.roles[0] = juce::AudioChannelSet::left;
    layout.roles[1] = juce::AudioChannelSet::left;
    expect (! layout.isValid(), "duplicate roles invalid");
    layout.roles[1] = juce::AudioChannelSet::right;
    expect (layout.isValid() && ! layout.confirmed, "roles distinct from confirmation");
    ChannelLayout immersive;
    const auto surround = juce::AudioChannelSet::create7point1point4();
    immersive.count = surround.size();
    for (int i = 0; i < immersive.count; ++i) immersive.roles[i] = surround.getTypeOfChannel (i);
    expect (immersive.count == 12 && immersive.isValid(), "all 7.1.4 roles supported");
    OutputRouting routing;
    expect (routing.validate (2, 12).isNotEmpty(), "routing muted by default");
    routing.fileToDevice[0] = 11;
    expect (routing.validate (2, 12).isEmpty(), "physical output 12 supported");
    expect (routing.validate (2, 2).isNotEmpty(), "device shrink invalidates route");
    routing.fileToDevice[1] = 11;
    expect (routing.validate (2, 12).isNotEmpty(), "duplicate outputs rejected");
    routing.fileToDevice[1] = 12;
    expect (routing.validate (2, 12).isNotEmpty(), "out-of-range output rejected");
}

void bufferChecks()
{
    auto buffer = std::make_unique<AuditionBuffer>();
    OutputRouting map;
    map.fileToDevice[0] = 11;
    map.fileToDevice[1] = 0;
    buffer->reset (2, map);
    buffer->gain.store (1.0f);
    juce::AudioBuffer<float> input (2, 4096), output (12, 4096);
    for (int frame = 0; frame < 4096; ++frame)
    {
        input.setSample (0, frame, 0.25f);
        input.setSample (1, frame, -0.5f);
    }
    for (int iteration = 0; iteration < 20; ++iteration)
    {
        buffer->push (input.getArrayOfReadPointers(), 4096);
        buffer->render (output.getArrayOfWritePointers(), 12, 4096);
        expect (std::abs (output.getSample (11, 4095) * 2 + output.getSample (0, 4095)) < 1e-6f,
                "routing keeps file channel order at wrap-around");
        expect (output.getMagnitude (1, 0, 4096) == 0.0f, "unused output silent");
    }
    expect (output.getSample (11, 4095) > 0.249f, "output 12 receives file channel 1");
    expect (buffer->playedFrames() == 81920 && buffer->drained(), "frame accounting after wrap");
    buffer->render (output.getArrayOfWritePointers(), 12, 4096);
    expect (output.getMagnitude (0, 4096) == 0.0f && buffer->underrunBlocks() == 1, "underrun is silence");
    buffer->endOfFile.store (true);
    buffer->render (output.getArrayOfWritePointers(), 12, 4096);
    expect (buffer->underrunBlocks() == 1, "EOF silence is not underrun");

    // Exercise actual concurrent producer/consumer memory ordering across wraps.
    map.fileToDevice[0] = 0;
    buffer->reset (1, map);
    constexpr int blocks = 1000, frames = 127;
    std::thread producer ([&]
    {
        juce::AudioBuffer<float> block (1, frames);
        for (int i = 0; i < blocks; ++i)
        {
            while (buffer->writable() < frames) std::this_thread::yield();
            for (int f = 0; f < frames; ++f) block.setSample (0, f, 0.125f);
            buffer->push (block.getArrayOfReadPointers(), frames);
        }
        buffer->endOfFile.store (true);
    });
    while (! buffer->endOfFile.load() || ! buffer->drained())
        buffer->render (output.getArrayOfWritePointers(), 12, frames);
    producer.join();
    expect (buffer->playedFrames() == blocks * frames, "concurrent ring preserves total frames");
}

void jobChecks()
{
    SessionJob job;
    expect (job.begin (JobStage::ingest, 2), "job starts");
    const auto ticket = job.ticket;
    expect (! job.begin (JobStage::analysis, 2), "overlapping jobs refused");
    expect (! job.finish ({ ticket.id + 1, 2 }, 2, ErrorCode::none, "stale"), "wrong job result ignored");
    expect (! job.finish (ticket, 3, ErrorCode::none, "stale"), "stale revision discarded");
    expect (job.state == JobState::cancelled, "stale job terminates");
    expect (job.begin (JobStage::ingest, 3), "new job after terminal state");
    job.requestCancel();
    expect (! job.finish (job.ticket, 3, ErrorCode::none, "late success"), "cancel wins completion race");
    expect (job.begin (JobStage::ingest, 3), "retry allowed");
    expect (job.finish (job.ticket, 3, ErrorCode::none, "complete"), "success accepted");
    expect (! job.finish (job.ticket, 3, ErrorCode::none, "duplicate"), "duplicate completion rejected");
}

void sessionChecks (const juce::File& directory)
{
    StandaloneSessionModel session;
    expect (! session.initializeCurrentStateFromSource(), "empty source cannot initialise current");
    JobControl control;
    const auto file = directory.getChildFile ("2.wav");
    const auto result = validateMedia (file, control);
    session.setValidatedSource (file, result.media);
    expect (! session.setValidatedSource (directory.getChildFile ("deleted.wav"), result.media)
            && session.sourceDocument.file == file, "changed source at commit keeps previous source");
    expect (! session.canInitializeCurrentStateFromSource(), "unconfirmed roles gate current state");
    session.setChannelRole (0, juce::AudioChannelSet::right);
    session.setChannelRole (1, juce::AudioChannelSet::left);
    expect (session.confirmLayout() && session.initializeCurrentStateFromSource(), "explicit reversed file order accepted");
    expect (session.currentStateNeedsRemeasurement(), "neutral current state still unmeasured");
    expect (session.sourceDocument.layout.roles[0] == juce::AudioChannelSet::right, "no sorted-set channel reorder");
    session.source.hasMeasurement = true;
    session.currentState.hasMeasurement = true;
    session.setChannelRole (0, juce::AudioChannelSet::centre);
    expect (! session.source.hasMeasurement && ! session.currentState.hasMeasurement
            && ! session.hasCurrentStateInitialized(), "role edit invalidates dependent truth");
    for (int profile = 0; profile < StandaloneSessionModel::numTargetProfiles; ++profile)
    {
        session.setSelectedTargetProfileIndex (profile);
        const auto gates = session.getWorkflowReadiness();
        expect (! gates.canDraftProposal() && ! gates.canExport()
                && gates.targetLimitsAuthoritative == (profile == 1), "only EBU scalar rules authoritative; proposal/export remain unavailable");
    }
    session.confirmLayout();
    session.initializeCurrentStateFromSource();
    session.clearSourceFile();
    expect (! session.sourceDocument.hasSource && ! session.hasCurrentStateInitialized(), "clear source clears dependent state");
    ServiceContext request;
    expect (! request.isValid(), "empty service request invalid");
    request.media = result.media;
    request.fileLayout = result.media.suggestedLayout;
    request.fileLayout.roles[0] = juce::AudioChannelSet::left;
    request.fileLayout.roles[1] = juce::AudioChannelSet::right;
    request.fileLayout.confirmed = true;
    request.range = { 0, result.media.frames };
    expect (request.isValid() && request.isFullProgram(), "full-program service snapshot valid");
    request.range = { 10, 100 };
    expect (request.isValid() && ! request.isFullProgram(), "segment cannot masquerade as full program");
    request.range.end = result.media.frames + 1;
    expect (! request.isValid(), "out-of-bounds analysis range rejected");
}

class FakeDevice final : public juce::AudioIODevice
{
public:
    FakeDevice() : juce::AudioIODevice ("Test outputs", "Fixture") {}
    juce::StringArray getOutputChannelNames() override
    {
        juce::StringArray names;
        for (int i = 0; i < 12; ++i) names.add (juce::String (i + 1));
        return names;
    }
    juce::StringArray getInputChannelNames() override { return {}; }
    juce::Array<double> getAvailableSampleRates() override { return { 48000 }; }
    juce::Array<int> getAvailableBufferSizes() override { return { 256 }; }
    int getDefaultBufferSize() override { return 256; }
    juce::String open (const juce::BigInteger& input, const juce::BigInteger& output, double rate, int) override
    {
        expect (input.isZero(), "audition never opens inputs");
        outputs = output; sampleRate = rate; opened = true;
        return {};
    }
    void close() override { stop(); opened = false; }
    bool isOpen() override { return opened; }
    void start (juce::AudioIODeviceCallback* next) override { callback = next; next->audioDeviceAboutToStart (this); }
    void stop() override { if (callback) callback->audioDeviceStopped(); callback = nullptr; }
    bool isPlaying() override { return callback != nullptr; }
    juce::String getLastError() override { return {}; }
    int getCurrentBufferSizeSamples() override { return 256; }
    double getCurrentSampleRate() override { return sampleRate; }
    int getCurrentBitDepth() override { return 24; }
    juce::BigInteger getActiveOutputChannels() const override { return outputs; }
    juce::BigInteger getActiveInputChannels() const override { return {}; }
    int getOutputLatencyInSamples() override { return 0; }
    int getInputLatencyInSamples() override { return 0; }
    void tick (juce::AudioBuffer<float>& block)
    {
        if (callback) callback->audioDeviceIOCallbackWithContext (nullptr, 0, block.getArrayOfWritePointers(),
                                                                block.getNumChannels(), block.getNumSamples(), {});
    }
    bool opened = false;
    double sampleRate = 48000;
    juce::BigInteger outputs;
    juce::AudioIODeviceCallback* callback = nullptr;
};

void auditionChecks (const juce::File& directory)
{
    auto device = std::make_unique<FakeDevice>();
    auto* fake = device.get();
    auto audio = std::make_unique<StandaloneAudition> (std::move (device));
    expect (! fake->opened && ! audio->connected(), "constructing engine does not open device");
    expect (audio->connect (44100).isNotEmpty() && ! audio->connected(), "sample rate mismatch refused");
    expect (audio->connect (48000).isEmpty() && audio->connected(), "explicit open accepted");
    expect (fake->outputs.countNumberOfSetBits() == 12, "opens all 12 output channels explicitly");
    const auto file = directory.getChildFile ("2.wav");
    JobControl control;
    const auto result = validateMedia (file, control);
    auto layout = result.media.suggestedLayout;
    layout.roles[0] = juce::AudioChannelSet::left;
    layout.roles[1] = juce::AudioChannelSet::right;
    OutputRouting map;
    map.fileToDevice[0] = 11; map.fileToDevice[1] = 0;
    expect (audio->play (file, result.media, layout, map).isNotEmpty(), "unconfirmed roles block playback");
    layout.confirmed = true;
    expect (audio->play (file, result.media, layout, map).isEmpty(), "source audition starts buffering");
    expect (audio->getState() == AuditionState::buffering && ! fake->isPlaying(), "prefill before audio callback");
    juce::AudioBuffer<float> block (12, 256);
    bool received = false;
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds (5);
    while (audio->active() && std::chrono::steady_clock::now() < deadline)
    {
        audio->poll();
        block.clear();
        fake->tick (block);
        received = received || block.getMagnitude (11, 0, block.getNumSamples()) > 0;
        expect (block.getMagnitude (1, 0, block.getNumSamples()) == 0, "unmapped device channel stays silent");
        std::this_thread::sleep_for (std::chrono::milliseconds (1));
    }
    expect (received && audio->getState() == AuditionState::finished, "audition reaches EOF");
    expect (audio->underruns() == 0, "prefilled fixture plays without underrun");
    expect (audio->play (file, result.media, layout, map).isEmpty(), "restart after EOF");
    audio->stop();
    expect (! fake->isPlaying() && audio->getState() == AuditionState::ready, "stop while buffering joins worker");
    expect (audio->play (file, result.media, layout, map).isEmpty(), "start before device loss");
    fake->close();
    audio->poll();
    expect (audio->getState() == AuditionState::failed && ! audio->active(), "device loss terminates playback");
    audio->disconnect();
    expect (! audio->connected(), "explicit close releases device");
}
}

int main()
{
    try
    {
        TempDirectory fixtures;
        mediaChecks (fixtures.path);
        routingChecks();
        bufferChecks();
        jobChecks();
        sessionChecks (fixtures.path);
        auditionChecks (fixtures.path);
        std::cout << "PASS: " << checks << " hardware-independent scaffold checks\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "FAIL after " << checks << " checks: " << error.what() << '\n';
        return 1;
    }
}

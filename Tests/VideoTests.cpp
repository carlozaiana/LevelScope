#include "Standalone/StandaloneVideo.h"
#include "Standalone/StandaloneSessionModel.h"
#include "Standalone/StandaloneAudition.h"
#include <iostream>
#include <chrono>

#if JUCE_MAC
void beginVideoReaderAudit();
int videoReaderAuditCount (int);
#else
void beginVideoReaderAudit() {}
int videoReaderAuditCount (int) { return 0; }
#endif
using namespace levelscope::standalone;
using namespace levelscope::standalone::video;
namespace
{
int checks = 0;
void expect (bool ok, const char* name) { ++checks; if (!ok) throw std::runtime_error (name); }
template <typename F> void rejects (F f, const char* name)
{ bool threw = false; try { f(); } catch (const std::exception&) { threw = true; } expect (threw, name); }
void timing()
{
    expect (Time (2002, 60000) == Time (1001, 30000), "rational reduction");
    expect (Time (-1, 3).floorTicks (2) == -1, "negative rounding is floor");
    expect (Time::parse ("-0.250") == Time (-1, 4), "negative decimal offset");
    expect (Time::parse ("1001/30000") == Time (1001, 30000), "exact fractional offset");
    rejects ([] { Time (1, 0); }, "zero denominator");
    rejects ([] { Time (INT64_MAX, 1) + Time (1, 1); }, "addition overflow");
    rejects ([] { Time (INT64_MIN, 1) - Time (1, 1); }, "subtraction overflow");
    rejects ([] { Time (INT64_MAX, 1).floorTicks (192000); }, "rescale overflow");
    rejects ([] { Time::parse ("0.1234567890"); }, "decimal precision bound");
    expect ((Time (INT64_MAX, 48000) - Time (INT64_MAX - 1, 48000)) == Time (1, 48000), "large adjacent sample positions");
    for (const auto rate : { Time (24, 1), Time (25, 1), Time (30, 1), Time (24000, 1001), Time (30000, 1001), Time (60000, 1001) })
    {
        std::vector<FrameTime> frames;
        const Time duration (rate.denominator, rate.numerator);
        for (int64_t f = 0; f < 100; ++f) frames.push_back ({ Time (f, 1) * duration, Time (f + 1, 1) * duration });
        expect (!frameAt (frames, Time (-1, 48000)), "before coverage blank");
        expect (!frameAt (frames, frames.back().end), "end-exclusive coverage");
        for (size_t f = 1; f < frames.size(); ++f)
        {
            expect (frameAt (frames, frames[f].pts) == f, "exact boundary chooses new frame");
            expect (frameAt (frames, frames[f].pts - Time (1, 192000)) == f - 1, "before boundary chooses prior frame");
            expect (frameAt (frames, frames[f].pts + Time (1, 192000)) == f, "after boundary chooses new frame");
        }
        // Independent integer formula across more than one hour. No accumulated double steps.
        for (int64_t sample = 0; sample <= 48000LL * 3601; sample += 4799)
        {
            const auto expected = static_cast<int64_t> ((__int128) sample * rate.numerator / (48000 * rate.denominator));
            expect ((Time (sample, 48000) * rate).floorTicks (1) == expected, "one-hour rational timeline has no drift");
        }
    }
    std::vector<FrameTime> vfr { { Time (0, 1), Time (1, 24) }, { Time (1, 24), Time (3, 24) }, { Time (4, 24), Time (5, 24) } };
    expect (frameAt (vfr, Time (2, 24)) == 1 && !frameAt (vfr, Time (3, 24)), "VFR durations and explicit gap");
    expect (frameAt (vfr, Time (1, 24) - Time (1, 24)) == 0, "positive offset delays picture");
    expect (frameAt (vfr, Time (0, 1) - Time (-1, 24)) == 1, "negative offset advances picture");
    expect (timecodeLabel (1799, Time (30000, 1001), true) == "00:00:59;29", "drop-frame preceding minute");
    expect (timecodeLabel (1800, Time (30000, 1001), true) == "00:01:00;02", "drop-frame minute");
    expect (timecodeLabel (17982, Time (30000, 1001), true) == "00:10:00;00", "drop-frame ten minutes");
    expect (timecodeLabel (107892, Time (30000, 1001), true) == "01:00:00;00", "drop-frame hour");
    expect (timecodeLabel (3600, Time (60000, 1001), true) == "00:01:00;04", "59.94 drop-frame minute");
    expect (timecodeLabel (35964, Time (60000, 1001), true) == "00:10:00;00", "59.94 ten minutes");
    expect (timecodeLabel (-25, Time (25, 1), false) == "-00:00:01:00", "negative timecode label");
    rejects ([] { timecodeLabel (0, Time (24, 1), true); }, "invalid drop-frame rate refused");
}
void clockChecks()
{
    auto clock = std::make_unique<AudioPresentationClock>();
    clock->submit ({ 7, 1000000000, 24000, 480 });
    expect (clock->position (999999999, 48000, 7, 24000) == 24000, "preroll does not advance");
    expect (clock->position (1005000000, 48000, 7, 24000) == 24240, "presented span interpolates");
    expect (clock->position (1020000000, 48000, 7, 24000) == 24480, "clock cannot run beyond submitted audio");
    clock->submit ({ 7, 1010000000, 24480, 0 });
    expect (clock->position (1015000000, 48000, 7, 24000) == 24480, "underrun silence holds source position");
    clock->submit ({ 8, 1020000000, 100, 480 });
    expect (clock->position (1025000000, 48000, 8, 100) == 340, "seek generation rejects old spans");
    expect (clock->position (1025000000, 48000, 9, 0) == 0, "no stale generation fallback");
    for (uint64_t i = 0; i < 5000; ++i) clock->submit ({ 10, 2000000000 + i * 10000000, static_cast<int64_t> (i * 480), 480 });
    expect (clock->position (2000000000 + 4999ULL * 10000000 + 5000000, 48000, 10, 0) == 4999 * 480 + 240, "bounded clock ring wraps");
    clock->reset(); expect (clock->position (90000000000, 48000, 10, 99) == 99, "flush resets submitted spans");
}
void independence()
{
    StandaloneSessionModel session;
    session.source.hasMeasurement = true; session.currentState.hasMeasurement = true;
    session.sourceDocument.routing.fileToDevice[0] = 3;
    const auto identity = OfflineAnalysisService::identityFor (session.measurementContext (MeasurementDomain::source), MeasurementDomain::source, "fixture");
    const auto revision = session.revision;
    session.referenceVideo.file = juce::File ("/tmp/reference.mov"); session.referenceVideo.sha256 = "picture";
    session.referenceVideo.startOnAudio = Time (-1001, 30000); session.referenceVideo.streamId = 2;
    expect (session.revision == revision && session.source.hasMeasurement && session.currentState.hasMeasurement, "video edits preserve measured state");
    expect (session.sourceDocument.routing.fileToDevice[0] == 3, "video edits preserve routing");
    expect (OfflineAnalysisService::identityFor (session.measurementContext (MeasurementDomain::source), MeasurementDomain::source, "fixture") == identity,
            "video attachment and offset do not enter audio measurement identity");
}
class FakeDevice final : public juce::AudioIODevice
{
public:
    FakeDevice() : AudioIODevice ("Video test", "Injected") {}
    juce::StringArray getOutputChannelNames() override { return { "L", "R" }; }
    juce::StringArray getInputChannelNames() override { return {}; }
    juce::Array<double> getAvailableSampleRates() override { return { 48000 }; }
    juce::Array<int> getAvailableBufferSizes() override { return { 256 }; }
    int getDefaultBufferSize() override { return 256; }
    juce::String open (const juce::BigInteger& inputs, const juce::BigInteger& outputs, double, int) override
    { expect (inputs.isZero(), "no input capture"); channels = outputs; opened = true; return {}; }
    void close() override { stop(); opened = false; }
    bool isOpen() override { return opened; }
    void start (juce::AudioIODeviceCallback* c) override { callback = c; c->audioDeviceAboutToStart (this); }
    void stop() override { if (callback) callback->audioDeviceStopped(); callback = nullptr; }
    bool isPlaying() override { return callback != nullptr; }
    juce::String getLastError() override { return {}; }
    int getCurrentBufferSizeSamples() override { return 256; }
    double getCurrentSampleRate() override { return 48000; }
    int getCurrentBitDepth() override { return 24; }
    juce::BigInteger getActiveOutputChannels() const override { return channels; }
    juce::BigInteger getActiveInputChannels() const override { return {}; }
    int getOutputLatencyInSamples() override { return 512; }
    int getInputLatencyInSamples() override { return 0; }
    void tick (juce::AudioBuffer<float>& block)
    { if (callback) callback->audioDeviceIOCallbackWithContext (nullptr, 0, block.getArrayOfWritePointers(), 2, block.getNumSamples(), {}); }
    juce::BigInteger channels;
    bool opened = false;
    juce::AudioIODeviceCallback* callback = nullptr;
};
void transportChecks(bool current = false)
{
    auto file = juce::File::getSpecialLocation (juce::File::tempDirectory).getNonexistentChildFile ("levelscope-video-transport", ".wav", false);
    struct Remove { juce::File file; ~Remove() { file.deleteFile(); } } cleanup { file };
    {
        juce::WavAudioFormat wav;
        auto stream = file.createOutputStream(); expect (stream != nullptr, "transport fixture file");
        std::unique_ptr<juce::AudioFormatWriter> writer (wav.createWriterFor (stream.get(), 48000, 2, 24, {}, 0));
        expect (writer != nullptr, "transport fixture writer"); stream.release();
        juce::AudioBuffer<float> data (2, 9000);
        for (int i = 0; i < 9000; ++i) { data.setSample (0, i, i * 0.0001f); data.setSample (1, i, -i * 0.0001f); }
        expect (writer->writeFromAudioSampleBuffer (data, 0, 9000), "transport fixture payload");
    }
    JobControl control; auto media = validateMedia (file, control);
    expect (media.succeeded(), "validated transport source");
    auto layout = media.media.suggestedLayout; layout.roles[0] = juce::AudioChannelSet::left; layout.roles[1] = juce::AudioChannelSet::right; layout.confirmed = true;
    OutputRouting routing; routing.fileToDevice[0] = 0; routing.fileToDevice[1] = 1;
    auto device = std::make_unique<FakeDevice>(); auto* fake = device.get();
    auto audio = std::make_unique<StandaloneAudition> (std::move (device));
    audio->sourceChanged (media.media);
    auto request = AuditionRequest::source(file, media.media, layout);
    request.currentInitialized = true;
    if (current) { request.domain = MeasurementDomain::currentState; request.context.processing.staticGainDb = -6.020599913279624; request.context.ticket.sessionRevision = 1; }
    auto play = [&](const OutputRouting& map, bool restart = false, bool muted = false)
    { audio->synchronize(request,map); return audio->play(request,map,restart,muted); };
    auto seek = [&](int64_t frame)
    { audio->synchronize(request,routing); return audio->seek(frame,request,routing); };

    expect (!audio->connected(), "source replacement does not open device");
    expect (seek (4000).isEmpty() && audio->positionFrames() == 4000 && !audio->active(), "seek while device closed retains stopped position");
    expect (audio->connect (48000).isEmpty(), "injected output open");
    expect (play (routing).isEmpty(), "resume at stopped seek");
    for (int i = 0; i < 5000 && audio->preparing(); ++i)
    { audio->poll (false); std::this_thread::sleep_for (std::chrono::milliseconds (1)); }
    expect (!fake->isPlaying() && audio->getState() == AuditionState::buffering, "picture preroll gates callback start");
    auto start = [&]
    {
        for (int i = 0; i < 200 && !fake->isPlaying(); ++i)
        { audio->poll(); std::this_thread::sleep_for (std::chrono::milliseconds (1)); }
        expect (fake->isPlaying(), "preroll starts callback");
    };
    start(); juce::AudioBuffer<float> block (2, 256); fake->tick (block);
    const auto scale = juce::Decibels::decibelsToGain (-30.0f) * 0.002f * (current ? 0.5f : 1.0f);
    expect (std::abs (block.getSample (0, 0) - 0.4f * scale) < 1e-8f, "seek decodes selected source frame");
    expect (block.getSample (0, 0) == -block.getSample (1, 0), "channel routing preserved");
    std::this_thread::sleep_for (std::chrono::milliseconds (25));
    expect (audio->positionFrames() == 4256, "clock stops at last submitted frame");
    const auto generation = audio->transportGeneration(); audio->stop();
    expect (audio->positionFrames() == 4256 && !fake->isPlaying() && audio->transportGeneration() > generation, "stop retains presented position and flushes generation");
    play (routing); start(); fake->tick (block);
    expect (std::abs (block.getSample (0, 0) - 0.4256f * scale) < 1e-8f, "resume starts at retained frame");
    expect (seek (1000).isEmpty() && audio->getState() == AuditionState::buffering, "playing seek enters preroll");
    start(); fake->tick (block);
    expect (std::abs (block.getSample (0, 0) - 0.1f * scale) < 1e-8f, "playing seek flushes old audio samples");
    for (int frame : { 2000, 7000, 200, 5000 }) seek (frame);
    start(); fake->tick (block);
    expect (std::abs (block.getSample (0, 0) - 0.5f * scale) < 1e-8f, "repeated seek publishes only final request");
    play (routing, true); start(); fake->tick (block);
    expect (block.getSample (0, 0) == 0, "restart seeks zero");
    seek (9000);
    expect (!audio->active() && audio->getState() == AuditionState::finished && audio->positionFrames() == 9000, "seek to EOF stops transport");
    expect (seek (-1).isNotEmpty(), "negative audio seek rejected");
    seek (8800);
    play (routing); start(); fake->tick (block);
    std::this_thread::sleep_for (std::chrono::milliseconds (25)); audio->poll();
    expect (audio->getState() == AuditionState::finished && audio->positionFrames() == 9000 && !fake->isPlaying(), "natural EOF drains presentation before stopping");
    OutputRouting muted;
    expect (play (muted).isNotEmpty(), "normal audition refuses all-muted routes");
    expect (play (muted, true, true).isEmpty(), "explicit silent timing run accepts all-muted routes");
    start(); fake->tick (block);
    expect (block.getMagnitude (0, 256) == 0, "silent timing run delivers only zero output");
    const auto beforeReplace = audio->transportGeneration(); audio->sourceChanged (media.media);
    expect (!fake->isPlaying() && audio->positionFrames() == 0 && audio->transportGeneration() > beforeReplace, "source replacement stops and resets the transport");
    play (routing); fake->close(); audio->poll();
    for (int i = 0; i < 5000 && audio->preparing(); ++i)
    { audio->poll(); std::this_thread::sleep_for (std::chrono::milliseconds (1)); }
    expect (audio->getState() == AuditionState::failed && !audio->active(), "device error stops and invalidates transport");
}
Player::View awaitPicture (Player& player, Time time, uint64_t generation, int timeoutSeconds = 15)
{
    const auto end = std::chrono::steady_clock::now() + std::chrono::seconds (timeoutSeconds);
    Player::View view;
    do
    {
        view = player.update (time, generation);
        if (view.state == Player::State::failed) throw std::runtime_error (view.message.toStdString());
        if (view.picture) return view;
        std::this_thread::sleep_for (std::chrono::milliseconds (2));
    } while (std::chrono::steady_clock::now() < end);
    throw std::runtime_error ("Timed out waiting for picture");
}
void mediaChecks (const juce::File& file)
{
    auto decoder = makeDecoder();
    auto asset = decoder->load (file, 0, [] { return false; });
    expect (asset->frames.size() == 120, "deterministic fixture has 120 frames");
    expect (!asset->sha256.isEmpty() && asset->selectedStream > 0, "video content and stream identity");
    const bool vfr = file.getFileName() == "vfr-exact.mov";
    for (int64_t i = 0; i < 120; ++i)
    {
        const auto ticks = i + (vfr ? i / 3 : 0);
        // AVAssetWriter gives the final VFR sample one base-frame duration;
        // endSession adds a trailing gap, which must not extend that picture.
        const auto endTicks = vfr && i == 119 ? ticks + 1 : i + 1 + (vfr ? (i + 1) / 3 : 0);
        expect (asset->frames[static_cast<size_t> (i)].pts == Time (ticks * 1001, 30000), "index matches generated presentation timeline");
        expect (asset->frames[static_cast<size_t> (i)].end == Time (endTicks * 1001, 30000), "index preserves generated frame coverage");
    }
    if (file.getFileName() == "timecode.mov")
    {
        std::cout << "TIMECODE " << asset->timecodeStatus << '\n';
        expect (asset->timecode.has_value(), "native timecode track recognized");
        expect (asset->timecode->startFrame == 107892 && asset->timecode->dropFrame
                && asset->timecode->rate == Time (30000, 1001), "native drop-frame origin and rational rate preserved");
    }
    for (size_t i : { 0u, 45u, 1u, 119u, 60u, 2u })
    {
        decoder->seek (i); const auto picture = decoder->next();
        expect (picture && picture->index == i && !picture->image.isNull(), "native random seek selects exact indexed PTS");
    }
    decoder->seek (0);
    for (size_t i = 0; i < 120; ++i)
    {
        auto frame = decoder->next(); expect (frame && frame->index == i, "native sequential PTS preserved");
        // Generator squares are near the bottom in the top-down decoded image.
        int burnIn = 0;
        for (int bit = 0; bit < 8; ++bit)
            if (frame->image.getPixelAt (20 + bit * 30, 160).getBrightness() > 0.5f) burnIn |= 1 << bit;
        expect (burnIn == static_cast<int> (i), "decoded burn-in agrees with PTS frame index");
    }
    expect (!decoder->next(), "native EOF is explicit");
    if (asset->streams.size() > 1)
    {
        auto second = makeDecoder(); auto secondAsset = second->load (file, asset->streams[1].id, [] { return false; });
        expect (secondAsset->selectedStream == asset->streams[1].id, "explicit second video stream selected");
        second->seek (45); expect (second->next()->index == 45, "second stream uses its own PTS index");
        expect (secondAsset->embeddedAudioTracks == 1, "embedded audio reported independently");
    }
    rejects ([&] { auto invalid = makeDecoder(); invalid->load (file, 9999, [] { return false; }); }, "missing stream rejected");
    Player player; player.attach (file);
    uint64_t generation = 100;
    for (size_t index : { 0u, 65u, 5u, 99u, 20u, 119u })
    {
        ++generation;
        auto view = awaitPicture (player, asset->frames[index].pts, generation);
        expect (view.picture->index == index && view.picture->generation == generation, "worker seek rejects stale generations");
        expect (view.queued <= 8, "bounded picture queue");
    }
    auto before = player.update (Time (-1, 48000), ++generation);
    expect (!before.picture, "before coverage has no stale picture");
    auto after = player.update (asset->frames.back().end, ++generation);
    expect (!after.picture, "after coverage has no stale picture");
    player.clear(); expect (!player.update ({}, ++generation).picture, "clear flushes picture");
    player.attach (file); player.attach (file); player.clear();
    expect (player.update ({}, ++generation).state == Player::State::empty, "replacement load cannot publish after clear");
    for (int i = 0; i < 5; ++i) { Player closing; closing.attach (file); }
    player.attach (file.getSiblingFile ("missing-reference.mov"));
    for (int i = 0; i < 200 && player.update ({}, generation).state != Player::State::failed; ++i)
        std::this_thread::sleep_for (std::chrono::milliseconds (2));
    const auto failed = player.update ({}, generation);
    expect (failed.state == Player::State::failed && !failed.picture && player.prerollReady(), "missing video reports failure without blocking audio");
    auto copy = juce::File::getSpecialLocation (juce::File::tempDirectory).getNonexistentChildFile ("levelscope-video-disappearing", ".mov", false);
    struct Remove { juce::File file; ~Remove() { file.deleteFile(); } } cleanup { copy };
    expect (file.copyFileTo (copy), "temporary video mutation fixture");
    Player disappearing; disappearing.attach (copy);
    (void) awaitPicture (disappearing, asset->frames[0].pts, 1);
    expect (copy.deleteFile(), "remove temporary video after loading");
    Player::View missing;
    for (int i = 0; i < 1000; ++i)
    {
        missing = disappearing.update (asset->frames[60].pts, 2);
        if (missing.state == Player::State::failed) break;
        std::this_thread::sleep_for (std::chrono::milliseconds (2));
    }
    expect (missing.state == Player::State::failed && !missing.picture && disappearing.prerollReady(), "file disappearance clears picture without blocking audio");
    expect (copy.replaceWithText ("This is deliberately not a media container."), "create corrupt-container fixture");
    rejects ([&] { auto corrupt = makeDecoder(); corrupt->load (copy, 0, [] { return false; }); }, "corrupt container rejected explicitly");
    std::cout << "MEDIA " << file.getFileName() << " sha256=" << asset->sha256 << " frames=" << asset->frames.size()
              << " vfr=" << asset->variableRate << " embedded_audio=" << asset->embeddedAudioTracks << '\n';
}
void longIndexChecks (const juce::File& file)
{
    auto decoder = makeDecoder();
    const auto asset = decoder->load (file, 0, [] { return false; });
    expect (asset->frames.size() == 54000, "long fixture frame count");
    for (int64_t i = 0; i < 54000; ++i)
    {
        expect (asset->frames[static_cast<size_t> (i)].pts == Time (i * 1001, 30000), "long fixture exact presentation timestamp");
        expect (asset->frames[static_cast<size_t> (i)].end == Time ((i + 1) * 1001, 30000), "long fixture exact coverage");
    }
    for (size_t i : { 0u, 27000u, 53999u, 45u })
    {
        decoder->seek (i); const auto frame = decoder->next();
        expect (frame && frame->index == i, "long fixture random seek PTS");
        int burnIn = 0;
        for (int bit = 0; bit < 8; ++bit)
            if (frame->image.getPixelAt (20 + bit * 30, 160).getBrightness() > 0.5f) burnIn |= 1 << bit;
        expect (burnIn == static_cast<int> (i & 255), "long fixture decoded burn-in");
    }
    std::cout << "LONG_INDEX " << file.getFileName() << " sha256=" << asset->sha256
              << " frames=" << asset->frames.size() << " exact_pts_and_coverage=PASS\n";
}
#include "ProResChecks.h"
}
int main (int argc, char** argv)
{
    try
    {
        if (argc >= 3 && juce::String (argv[1]).startsWith ("--"))
        {
            const juce::String mode (argv[1]); const juce::File file { juce::String (argv[2]) };
            if (mode == "--prores" && argc == 7)
                proResChecks (file, std::stoll (argv[3]), std::stoll (argv[4]), std::stoll (argv[5]), true, juce::String (argv[6]) == "vfr");
            else if (mode == "--actual") proResChecks (file, 2772, 25, 1, false, false);
            else if (mode == "--isolation") isolationChecks (file);
            else if (mode == "--loading") loadingChecks (file);
            else if (mode == "--reject")
            {
                auto decoder = makeDecoder(); bool rejected = false;
                try { decoder->load (file, 0, [] { return false; }); }
                catch (const std::exception& e) {
                    rejected = true; std::cout << "EXPECTED_REJECTION " << e.what() << '\n';
                    if (argc >= 4) expect (juce::String(e.what()).contains(argv[3]), "understandable expected rejection reason");
                }
                expect (rejected, "unsupported ProRes policy must reject");
            }
            else throw std::runtime_error ("Unknown ProRes test arguments");
            std::cout << "PASS: " << checks << " ProRes checks; controlled sink only, no physical audio device\n";
            return 0;
        }
        timing(); clockChecks(); independence(); transportChecks(); transportChecks(true);
        for (int i = 1; i < argc; ++i)
        {
            const juce::File file { juce::String (argv[i]) };
            if (file.getFileName() == "half-hour.mov") longIndexChecks (file);
            else mediaChecks (file);
        }
        std::cout << "PASS: " << checks << " video timing/transport checks; no audio device opened\n";
        return 0;
    }
    catch (const std::exception& error) { std::cerr << "FAIL after " << checks << ": " << error.what() << '\n'; return 1; }
}

#include "CompositionTestFixtures.h"
#include "Standalone/StandaloneWaveform.h"
#include <future>

namespace
{
std::shared_ptr<const VerifiedSource> acquire(const ServiceContext& c)
{
    JobControl control;
    auto result = prepareVerifiedSource(c, control);
    if (!result.value) throw std::runtime_error(result.message.toStdString());
    return *result.value;
}
struct PauseCopy : JobFeedback
{
    std::promise<void> entered, release;
    std::shared_future<void> released = release.get_future().share();
    bool paused = false;
    void publish(JobPhase phase, double fraction) override
    {
        JobFeedback::publish(phase, fraction);
        if (!paused && phase == JobPhase::copy)
        {
            paused = true;
            entered.set_value();
            released.wait();
        }
    }
};
enum class Cancellation { none, producer, waiter };
void concurrent(const AuditionRequest& request, Cancellation cancel)
{
    StandaloneSessionModel model;
    expect(model.setValidatedSource(request.context.sourceFile, request.context.media), "adopt concurrency fixture");
    const auto context = model.measurementContext(MeasurementDomain::source);
    JobControl first, second;
    auto pause = std::make_shared<PauseCopy>();
    first.feedback = pause;
    auto started = pause->entered.get_future();
    auto a = std::async(std::launch::async, [&] { return prepareVerifiedSource(context, first); });
    const bool began = started.wait_for(std::chrono::seconds(8)) == std::future_status::ready;
    auto b = std::async(std::launch::async, [&] { return prepareVerifiedSource(context, second); });
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(8);
    while (!model.verifiedAssets->counts().waits && std::chrono::steady_clock::now() < deadline)
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    const bool joined = model.verifiedAssets->counts().waits > 0;
    if (cancel == Cancellation::producer) first.cancel = true;
    if (cancel == Cancellation::waiter) second.cancel = true;
    const bool waiterLeft = cancel != Cancellation::waiter || b.wait_for(std::chrono::seconds(2)) == std::future_status::ready;
    pause->release.set_value();
    const auto producer = a.get(), waiter = b.get();
    expect(began && joined && waiterLeft, "concurrent acquisition waits without starting a second copy and waiter cancels promptly");
    if (cancel == Cancellation::none)
        expect(producer.value && waiter.value && *producer.value == *waiter.value, "two successful consumers coalesce one cold preparation");
    else if (cancel == Cancellation::producer)
        expect(producer.error == ErrorCode::cancelled && bool(waiter.value), "live waiter continues after producer cancellation");
    else
        expect(waiter.error == ErrorCode::cancelled && bool(producer.value), "producer continues after waiter cancellation");
    expect(model.verifiedAssets->counts().copies == 1, "one complete verification after competing requests");
    auto next = acquire(context);
    expect(next == *(producer.value ? producer.value : waiter.value), "subsequent consumer shares completed file");
}
void newDuringPreparation(const AuditionRequest& request)
{
    StandaloneSessionModel model;
    expect(model.setValidatedSource(request.context.sourceFile, request.context.media), "adopt source before background job");
    model.sourceDocument.layout = layoutFor(2);
    expect(model.job.begin(JobStage::analysis, model.revision), "start measured worker on original opening");
    const auto oldRequest = model.measurementContext(MeasurementDomain::source);
    JobControl worker;
    auto pause = std::make_shared<PauseCopy>();
    worker.feedback = pause;
    auto started = pause->entered.get_future();
    auto future = std::async(std::launch::async, [&] { return OfflineAnalysisService().analyze(oldRequest, MeasurementDomain::source, worker); });
    const bool began = started.wait_for(std::chrono::seconds(8)) == std::future_status::ready;
    model = StandaloneSessionModel();
    const bool adopted = model.setValidatedSource(request.context.sourceFile, request.context.media);
    model.sourceDocument.layout = layoutFor(2);
    pause->release.set_value();
    auto completed = future.get();
    expect(began && adopted && completed.value, "old worker safely finishes reading its leased source after New");
    expect(!model.acceptMeasurement(*completed.value) && model.verifiedAssets->counts().copies == 0 && !model.source.hasMeasurement,
           "old verified storage and measured result cannot attach to the new opening");
    auto fresh = acquire(model.measurementContext(MeasurementDomain::source));
    expect(fresh && model.verifiedAssets->counts().copies == 1, "same original in new session establishes its own trust");
}
#if JUCE_MAC
void writeWithMatchingStamps(const juce::File& file)
{
    auto writer = file.createOutputStream();
    expect(writer && writer->write("one", 3), "create same-size write fixture");
    writer->flush();
    auto captured = SourceFileVersion::capture(file);
    expect(captured.matches(file), "clean monitored file matches before writing");
    expect(writer->setPosition(0) && writer->write("two", 3), "overwrite three bytes through an already open file");
    writer->flush();

    // Isolate kernel change notification on both APFS and HFS+. Present all
    // the current stat fields as equal, just as for two same-tick HFS+ writes.
    const auto current = SourceFileVersion::read(file);
    expect(captured.sameFile(current) && captured.size == current.size, "same file and length after overwrite");
    captured.seconds = current.seconds;
    captured.nanos = current.nanos;
    auto secondConsumer = captured;
    expect(!captured.matches(file) && !secondConsumer.matches(file),
           "a same-size write with matching timestamps invalidates all consumers");
    expect(SourceFileVersion::capture(file).matches(file), "a fresh monitor starts clean after the write");
}
#endif
void lifecycle(const juce::File& root)
{
#if JUCE_MAC
    writeWithMatchingStamps(root.getChildFile("already-open.bin"));
#endif
    const auto request = fixture(root.getChildFile("source.wav"), 48000, 2, 192003);
    newDuringPreparation(request);
    concurrent(request, Cancellation::none);
    concurrent(request, Cancellation::waiter);
    concurrent(request, Cancellation::producer);
    StandaloneSessionModel model;
    expect(model.setValidatedSource(request.context.sourceFile, request.context.media), "adopt verified source");
    model.sourceDocument.layout = layoutFor(2);
    auto context = model.measurementContext(MeasurementDomain::source);
    auto first = acquire(context);
    auto copy = first->file();
    auto reader = first->createReader();
    expect(reader && copy != context.sourceFile, "reader leases a private file");
    auto changed = context;
    changed.processing.staticGainDb = -3;
    changed.ticket.sessionRevision++;
    expect(acquire(changed) == first && model.verifiedAssets->counts().copies == 1, "revision and gain share source storage");
    model.clearSourceFile();
    JobControl stale;
    expect(prepareVerifiedSource(context, stale).error == ErrorCode::staleRevision, "closed source generation cannot acquire from replacement");
    first.reset();
    juce::AudioBuffer<float> samples(2, 64);
    expect(copy.existsAsFile() && reader->read(samples.getArrayOfWritePointers(), 2, 948, 64), "active reader retains file after clearing source");
    reader.reset();
    expect(!copy.exists(), "last reader releases the old private file");

    expect(model.setValidatedSource(request.context.sourceFile, request.context.media), "new generation accepts original");
    model.sourceDocument.layout = layoutFor(2);
    context = model.measurementContext(MeasurementDomain::source);
    auto rebuilt = acquire(context);
    auto modified = rebuilt->file();
    expect(modified.setReadOnly(false), "owned fixture made writable for fault injection");
    {auto output = modified.createOutputStream(); expect(output && output->setPosition(80) && output->writeByte(77), "alter private copy in place");}
    expect(!rebuilt->matches(context.media) && !rebuilt->createReader(), "altered private bytes cannot open a reader");
    auto fresh = acquire(context);
    expect(fresh->file() != modified && fresh->intact() && model.verifiedAssets->counts().copies == 2, "altered private copy reverified from original");
    expect(fresh->file().deleteFile() && !fresh->intact(), "missing private copy invalidated");
    auto recovered = acquire(context);
    expect(recovered->intact() && model.verifiedAssets->counts().copies == 3, "missing private copy rebuilt");
    const auto path = recovered->file();
    const auto displaced = root.getChildFile("displaced-private.wav");
    expect(path.moveFileTo(displaced) && path.replaceWithText("unrelated replacement"), "replace the private pathname with an unrelated owned fixture");
    recovered.reset();
    auto replacement = acquire(context);
    expect(replacement->file() != path && path.loadFileAsString() == "unrelated replacement", "cleanup leaves an unrelated replacement untouched");
    const auto oldStamp = context.media.modified;
    {auto output = context.sourceFile.createOutputStream(); expect(output && output->setPosition(80) && output->writeByte(44), "mutate owned original fixture");}
    expect(context.sourceFile.setLastModificationTime(oldStamp) && context.media.matchesFileStamp(context.sourceFile), "mutation preserves legacy size/mtime");
    expect(!context.sourceFilesMatch() && prepareVerifiedSource(context, stale).error == ErrorCode::sourceChanged, "source mutation rejects reuse despite valid owned bytes");
}
void assembled(const juce::File& root)
{
    const auto f = fixtureComposition(root, 48000, CompositionPreset::stereoStereo, true, 144007);
    auto s = compositionModel(f.recipe);
    const auto context = s.measurementContext(MeasurementDomain::source);
    auto assets = acquire(context);
    expect(assets->files().size() == 4 && s.verifiedAssets->counts().copies == 4, "one private copy per physical composition asset");
    Capture source(2, f.expected.getNumSamples());
    measuredComposition(s, MeasurementDomain::source, &source);
    expect(sameSamples(source.data, f.expected), "shared composed Source samples equal independent oracle");
    expect(s.setCurrentGainDb(-3), "edit Current gain");
    auto current = measuredComposition(s, MeasurementDomain::currentState);
    expect(current.domain == MeasurementDomain::currentState && current.detail.frames == context.media.frames, "processed Current independently measured");
    const auto c = s.measurementContext(MeasurementDomain::currentState);
    WaveformCache waves;
    auto sink = waves.begin(WaveformBinding::from(c, MeasurementDomain::currentState));
    JobControl job;
    expect(bool(OfflineAnalysisService(nullptr, sink.get()).buildWaveform(c, MeasurementDomain::currentState, job).value), "independent waveform rebuild completes");
    auto artifact = OfflineAnalysisService().verifyTemporaryRender(c, job);
    expect(artifact.value && artifact.value->domain == MeasurementDomain::renderedOutput && artifact.value->measuredContentHash != c.media.sha256,
           "fresh rendered bytes independently verified");
    expect(s.verifiedAssets->counts().copies == 4 && s.verifiedAssets->counts().cacheHits >= 16, "waveform, measured domains and render reuse four physical files");
    // Exercise the same staged encoder, artifact measurement and final hash
    // checks used by the app, at a measured normalizing trim for this fixture.
    expect(s.setCurrentGainDb(-23.0 - *s.source.integratedLufs), "normalize known composition for delivery fixture");
    measuredComposition(s, MeasurementDomain::currentState);
    const auto exportRequest = ex::capture(s);
    auto encoded = ex::render(exportRequest, root, job);
    expect(bool(encoded.value), "export renderer uses session assets and measures the new encoded artifact");
    auto candidate = *encoded.value;
    const auto approval = ex::approvalIdentity(*candidate, exportRequest);
    expect(approval.isNotEmpty(), "independent artifact and Current pass delivery fixture");
    const auto destination = root.getChildFile("verified-export");
    auto delivered = ex::publish(*candidate, exportRequest, destination, approval, job);
    expect(delivered.value && juce::SHA256(destination.getChildFile("programme.wav")).toHexString() == candidate->hash(),
           "final exported PCM exactly matches the independently verified artifact");
    expect(s.verifiedAssets->counts().copies == 4, "normal export and finalization incur zero further source copies");
    auto saved = root.getChildFile("reopen.levelscope");
    expect(sf::save(sf::capture(s), saved, job).succeeded(), "save composition session without owned-file locators");
    auto document = sf::read(saved, job);
    expect(bool(document.document), "read saved session");
    auto verified = sf::verify(*document.document, job);
    expect(verified.prepared && verified.prepared->model.sessionId != s.sessionId, "opening verifies inputs in a distinct runtime");
    auto reopened = acquire(verified.prepared->model.measurementContext(MeasurementDomain::source));
    expect(verified.prepared->model.verifiedAssets->counts().copies == 4 && reopened->files()[0] != assets->files()[0],
           "reopened session establishes fresh trust and private storage");
    auto rerender = ex::render(ex::capture(verified.prepared->model), root, job, {}, nullptr, candidate->verifiedSource());
    expect(rerender.value && (*rerender.value)->verifiedSource()->files() == reopened->files()
           && (*rerender.value)->hash() == candidate->hash(), "explicit export hint from an old opening uses the receiving session files and produces identical PCM");

    auto group = *f.recipe.me.assembly;
    std::swap(group.members[0].file, group.members[1].file);
    std::swap(group.members[0].media, group.members[1].media);
    expect(group.confirmAlignment() && s.setValidatedProgramme(group), "replace composition with reordered separate-file assembly");
    auto reordered = acquire(s.measurementContext(MeasurementDomain::source));
    expect(s.verifiedAssets->counts().copies == 0 && reordered->files()[0] == assets->files()[1], "new assembly uses new order and retained unchanged assets");
}
}
int sessionSourceBenchmark(int argc, char** argv);
int main(int argc, char** argv)
{
    if (argc == 4 && juce::String(argv[1]) == "--benchmark") return sessionSourceBenchmark(argc - 1, argv + 1);
    if (argc != 2) return 2;
    const juce::File receipt(argv[1]);
    if (receipt.exists()) return 2;
    juce::String failure;
    try { Temp t; lifecycle(t.directory); assembled(t.directory.getChildFile("composed")); }
    catch (const std::exception& e) { failure = e.what(); }
    juce::var result(new juce::DynamicObject());
    auto* out = result.getDynamicObject();
    out->setProperty("passed", failure.isEmpty()); out->setProperty("failure", failure);
    out->setProperty("build", LEVELSCOPE_MEASUREMENT_BUILD); out->setProperty("checks", checks);
    std::cout << checks.size() << " checks " << failure << std::endl;
    return receipt.replaceWithText(juce::JSON::toString(result)) && failure.isEmpty() ? 0 : 1;
}

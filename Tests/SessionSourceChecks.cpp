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
void concurrent(const AuditionRequest& request, bool cancelProducer)
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
    if (cancelProducer) first.cancel = true;
    else second.cancel = true;
    const bool waiterLeft = cancelProducer || b.wait_for(std::chrono::seconds(2)) == std::future_status::ready;
    pause->release.set_value();
    const auto producer = a.get(), waiter = b.get();
    expect(began && joined && waiterLeft, "concurrent acquisition waits without starting a second copy and waiter cancels promptly");
    expect(cancelProducer ? producer.error == ErrorCode::cancelled && bool(waiter.value)
                          : waiter.error == ErrorCode::cancelled && bool(producer.value), "one cancellation preserves the other consumer");
    expect(model.verifiedAssets->counts().copies == 1, "one complete verification after competing requests");
    auto next = acquire(context);
    expect(next == *(cancelProducer ? waiter.value : producer.value), "subsequent consumer shares completed file");
}
void lifecycle(const juce::File& root)
{
    const auto request = fixture(root.getChildFile("source.wav"), 48000, 2, 192003);
    concurrent(request, false);
    concurrent(request, true);
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

#include "AuditionTestSupport.h"

// Device-free comparison driver. Compile the same driver against the qualified
// baseline and integrated source; the latter requests automatic preparation at
// the same point where the main window schedules it after Source analysis.
namespace
{
using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point start) { return std::chrono::duration<double>(Clock::now()-start).count(); }
struct CountedWriter final : SnapshotWriter
{
    CountedWriter(const juce::File& file,std::atomic<int64_t>& n): stream(file.createOutputStream()),bytes(n) {}
    bool write(const void* p,size_t n) override { if(!stream||!stream->write(p,n))return false;bytes+=int64_t(n);return true; }
    bool flush() override { if(!stream)return false;stream->flush();return stream->getStatus().wasOk(); }
    std::unique_ptr<juce::FileOutputStream> stream;std::atomic<int64_t>& bytes;
};
}
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI init;
    if(argc!=3)return 2;
    const juce::File file(argv[1]),output(argv[2]);
    auto result=juce::var(new juce::DynamicObject());auto* row=result.getDynamicObject();bool passed=false;
    try
    {
        expect(!output.exists(),"timing receipt is new");
        JobControl validation;auto begin=Clock::now();const auto media=validateMedia(file,validation);
        row->setProperty("validation_s",elapsed(begin));expect(media.succeeded(),"existing programme validated");
        expect(media.media.channels==2||media.media.channels==6,"timing programme is stereo or 5.1");
        row->setProperty("file",file.getFullPathName());row->setProperty("sha256",media.media.sha256);
        row->setProperty("channels",media.media.channels);row->setProperty("frames",media.media.frames);
        StandaloneSessionModel session;expect(session.setValidatedSource(file,media.media),"timing session imports confirmed input");
        session.sourceDocument.layout=layoutFor(media.media.channels);
        expect(session.job.begin(JobStage::analysis,session.revision),"initial analysis job begins");
        JobControl measurement;begin=Clock::now();auto measured=OfflineAnalysisService().analyze(session.measurementContext(MeasurementDomain::source),MeasurementDomain::source,measurement);
        row->setProperty("source_analysis_s",elapsed(begin));expect(measured.value.has_value(),"initial Source analysis completes");
        expect(session.job.finish(session.job.ticket,session.revision,measured.error,measured.message),"initial analysis job completes");
        expect(session.acceptMeasurement(std::move(*measured.value)),"initial measured Source accepted");
        auto device=std::make_unique<Device>();auto* d=device.get();auto audio=std::make_unique<StandaloneAudition>(std::move(device));
#if LS_BACKGROUND_TIMING
        const auto beforePlayback=session.verifiedAssets->counts();
        row->setProperty("session_copies_after_analysis",int64_t(beforePlayback.copies));
        row->setProperty("session_copied_bytes_after_analysis",int64_t(beforePlayback.copiedBytes));
#else
        std::atomic<int64_t> copiedBytes{0};std::atomic<int> copies{0};SnapshotStorage storage;
        storage.createWriter=[&](const juce::File& f){++copies;return std::make_unique<CountedWriter>(f,copiedBytes);};
        audio->setSnapshotStorage(storage);
#endif
        const auto map=routes(media.media.channels);auto request=session.auditionRequest();
        audio->synchronize(request,map);
#if LS_BACKGROUND_TIMING
        begin=Clock::now();expect(audio->prepare(request).isEmpty(),"automatic preparation accepts stopped request");
        waitPreparation(*audio);row->setProperty("automatic_preparation_s",elapsed(begin));
        expect(audio->playbackReady()&&!d->opened&&d->opens==0,"automatic preparation is ready without opening output");
        const auto afterWarmup=session.verifiedAssets->counts();
        row->setProperty("warmup_copy_count",int64_t(afterWarmup.copies-beforePlayback.copies));
        row->setProperty("warmup_copied_bytes",int64_t(afterWarmup.copiedBytes-beforePlayback.copiedBytes));
        row->setProperty("session_copies_after_warmup",int64_t(afterWarmup.copies));
        row->setProperty("session_copied_bytes_after_warmup",int64_t(afterWarmup.copiedBytes));
#else
        row->setProperty("automatic_preparation_s",0.0);
#endif
        expect(audio->connect(media.media.sampleRate).isEmpty(),"explicit injected output connection");
        audio->setMonitorGainDb(-40);
        auto play=[&](const char* name)
        {
            const auto start=Clock::now();expect(audio->play(request,map).isEmpty(),"explicit timed Play accepted");
            awaitPlaying(*audio,*d);juce::AudioBuffer<float> block(maxChannels,256);block.clear();d->tick(block);
            int first=0;for(;first<block.getNumSamples();++first){bool audible=false;for(int c=0;c<maxChannels;++c)audible|=block.getSample(c,first)!=0;if(audible)break;}
            expect(first<block.getNumSamples(),"non-zero sample delivered to injected output");
            row->setProperty(name,elapsed(start)+double(first)/media.media.sampleRate);audio->stop();
        };
        play("first_play_to_nonzero_output_s");play("subsequent_play_to_nonzero_output_s");
#if LS_BACKGROUND_TIMING
        const auto afterPlayback=session.verifiedAssets->counts();
        row->setProperty("audition_copy_count",int64_t(afterPlayback.copies-beforePlayback.copies));
        row->setProperty("audition_copied_bytes",int64_t(afterPlayback.copiedBytes-beforePlayback.copiedBytes));
        row->setProperty("session_total_copies",int64_t(afterPlayback.copies));row->setProperty("session_total_copied_bytes",int64_t(afterPlayback.copiedBytes));
        row->setProperty("playback_cache_hits",int64_t(afterPlayback.cacheHits-beforePlayback.cacheHits));row->setProperty("playback_coalesced_waits",int64_t(afterPlayback.waits-beforePlayback.waits));
#else
        row->setProperty("audition_copy_count",copies.load());row->setProperty("audition_copied_bytes",copiedBytes.load());
#endif
        row->setProperty("timing_scope","Device-free engine with 1 ms polling; first non-zero sample at injected app output, not physical or native 10 Hz UI latency");
        row->setProperty("physical_device_opened",false);row->setProperty("measurement_build",LEVELSCOPE_MEASUREMENT_BUILD);
        passed=true;
    }
    catch(const std::exception& e){row->setProperty("failure",e.what());std::cerr<<e.what()<<'\n';}
    row->setProperty("passed",passed);row->setProperty("checks",checks);
    if(!output.replaceWithText(juce::JSON::toString(result)))return 2;
    return passed?0:1;
}

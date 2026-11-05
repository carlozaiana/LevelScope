#include "AuditionTestSupport.h"
#include "Standalone/StandaloneMainComponent.h"
#include "Standalone/VerifiedSource.h"
#include <future>

namespace levelscope::standalone
{
struct PreparationFeedbackTestAccess
{
    static StandaloneSessionModel& model(StandaloneMainComponent& c) { return c.session; }
    static void tick(StandaloneMainComponent& c) { c.timerCallback(); }
    static void stopTimer(StandaloneMainComponent& c) { c.stopTimer(); }
    static void import(StandaloneMainComponent& c,const juce::File& f) { c.beginImport(f); }
    static void measure(StandaloneMainComponent& c) { c.beginMeasurement(MeasurementDomain::source); }
    static void cancel(StandaloneMainComponent& c) { c.cancelJobButton.onClick(); }
    static void refresh(StandaloneMainComponent& c) { c.refreshJobFeedback(); }
    static void present(StandaloneMainComponent& c,JobPhase phase,double fraction)
    {
        c.importControl=std::make_shared<JobControl>();c.importControl->feedback=std::make_shared<JobFeedback>();
        c.importControl->reportPhase(phase,fraction);c.refreshFromSession();c.refreshTimeline();
    }
    static bool stopped(StandaloneMainComponent& c) { return !c.audition.active()&&!c.audition.connected(); }
};
}
namespace
{
using Access=PreparationFeedbackTestAccess;
struct RecordingFeedback final : JobFeedback
{
    JobControl* cancelControl=nullptr;
    JobPhase cancelAt=JobPhase::starting;
    double cancelAfter=0;
    bool sawCancel=false,ordered=true;
    std::vector<JobPhase> phases;
    double previous=-1;
    std::function<void(JobPhase,double)> onPhase;
    void publish(JobPhase phase,double fraction) override
    {
        JobFeedback::publish(phase,fraction);
        if(phases.empty()||phases.back()!=phase){phases.push_back(phase);previous=-1;}
        ordered=ordered&&(fraction<0||fraction>=previous);previous=fraction;
        if(onPhase)onPhase(phase,fraction);
        if(cancelControl&&phase==cancelAt&&(fraction<0||fraction>=cancelAfter))
        {cancelControl->cancel.store(true);sawCancel=true;}
    }
};
void coherentPublication()
{
    JobFeedback feedback;std::atomic<bool> done{false};bool coherent=true;
    std::thread writer([&]{for(int n=0;n<100000;++n){feedback.publish(JobPhase::copy,.25);feedback.publish(JobPhase::measuring,.75);}done=true;});
    do{const auto s=feedback.snapshot();coherent&=(s.phase==JobPhase::starting&&s.fraction<0)||(s.phase==JobPhase::copy&&s.fraction==.25)||(s.phase==JobPhase::measuring&&s.fraction==.75);}while(!done);
    writer.join();expect(coherent,"concurrent stage and fraction publication is coherent");
    JobControl parent,child,grandchild;parent.feedback=std::make_shared<JobFeedback>();child.parent=&parent;grandchild.parent=&child;
    grandchild.reportPhase(JobPhase::samples,.375);expect(parent.feedback->snapshot().fraction==.375,"nested member workers report to the same bounded snapshot");
    expect(parent.progress.load()==0,"feedback never changes legacy numeric job progress");parent.cancel=true;expect(grandchild.isCancelled(),"parent cancellation still reaches nested verification");
    const auto pre=describeJobFeedback({JobPhase::copy,.5},true,false,false), measure=describeJobFeedback({JobPhase::measuring,.5},true,true,false);
    expect(pre.title.contains("preparation")&&pre.detail.contains("current file / pass")&&measure.title.contains("Current")&&measure.detail.contains("incomplete"),"preparation and measurement have explicit domains, units and incomplete authority");
    expect(describeJobFeedback({JobPhase::finalising,-1},true,false,false).fraction<0,"finalisation does not invent remaining time or a percentage");
}
void exact(const AnalysisSnapshot& a,const AnalysisSnapshot& b)
{
    expect(a.identity==b.identity&&a.measuredContentHash==b.measuredContentHash,"feedback preserves exact measurement identity");
    const auto& x=a.detail;const auto& y=b.detail;
    auto metric=[](const measurement::Metric& p,const measurement::Metric& q){return p.value==q.value&&p.status==q.status;};
    expect(metric(x.integrated,y.integrated)&&metric(x.lra,y.lra)&&metric(x.momentaryMax,y.momentaryMax)&&metric(x.shortTermMax,y.shortTermMax)&&metric(x.samplePeak,y.samplePeak)&&metric(x.truePeak,y.truePeak),"feedback preserves every metric and availability bit exactly");
    expect(x.frames==y.frames&&x.sampleRate==y.sampleRate&&x.lraStable==y.lraStable&&x.integratedRelativeGate==y.integratedRelativeGate&&x.lraRelativeGate==y.lraRelativeGate&&x.p10==y.p10&&x.p95==y.p95&&x.integratedIncluded==y.integratedIncluded&&x.lraIncluded==y.lraIncluded,"feedback preserves gates, populations, percentiles and coverage");
    expect(x.channelSamplePeaks==y.channelSamplePeaks&&x.channelTruePeaks==y.channelTruePeaks&&x.samplePeakFrames==y.samplePeakFrames&&x.truePeakEighthFrames==y.truePeakEighthFrames,"feedback preserves per-channel peaks and exact sample positions");
    auto windows=[](const auto& p,const auto& q){if(p.size()!=q.size())return false;for(size_t n=0;n<p.size();++n)if(p[n].endFrame!=q[n].endFrame||p[n].energy!=q[n].energy)return false;return true;};
    expect(windows(x.momentary,y.momentary)&&windows(x.shortTerm,y.shortTerm)&&windows(x.lraWindows,y.lraWindows),"feedback preserves every measured M, S and LRA window exactly");
}
void stagesAndCancellation(const juce::File& root,int channels)
{
    root.createDirectory();auto request=fixture(root.getChildFile("mix.wav"),48000,channels,240013);
    JobControl imported;auto importFeedback=std::make_shared<RecordingFeedback>();imported.feedback=importFeedback;
    auto valid=validateMedia(request.context.sourceFile,imported);
    expect(valid.succeeded()&&valid.media.sha256==request.context.media.sha256,"feedback-enabled import validates unchanged source identity");
    expect(importFeedback->phases==std::vector<JobPhase>{JobPhase::metadata,JobPhase::samples,JobPhase::sourceHash}&&importFeedback->ordered,"import reports real metadata, decoded-sample and identity stages with per-stage progress");
    auto storageDirectory=root.getChildFile("copies");storageDirectory.createDirectory();SnapshotStorage storage;storage.directory=storageDirectory;
    const auto original=juce::SHA256(request.context.sourceFile).toHexString();
    for(auto phase:{JobPhase::copy,JobPhase::copyHash,JobPhase::samples,JobPhase::sourceHash})
    {
        JobControl c;auto f=std::make_shared<RecordingFeedback>();c.feedback=f;f->cancelControl=&c;f->cancelAt=phase;f->cancelAfter=.25;
        const auto result=prepareVerifiedSource(request.context,c,storage);
        expect(f->sawCancel&&result.error==ErrorCode::cancelled&&!result.value,"cancellation inside each preparation pass rejects the private source");
        expect(storageDirectory.getNumberOfChildFiles(juce::File::findFiles)==0,"cancelled preparation deletes its owned partial source copy");
        expect(juce::SHA256(request.context.sourceFile).toHexString()==original,"preparation cancellation never changes original audio");
    }
    JobControl plain;const auto baseline=OfflineAnalysisService().analyze(request.context,MeasurementDomain::source,plain);expect(bool(baseline.value),"baseline measurement succeeds");
    const auto binding=WaveformBinding::from(request.context,MeasurementDomain::source);
    WaveformCache cache;auto wave=cache.begin(binding);LoudnessProgress loudness(binding);
    JobControl measured;auto f=std::make_shared<RecordingFeedback>();measured.feedback=f;bool partial=false,pending=false;
    f->onPhase=[&](JobPhase phase,double fraction)
    {
        if(phase==JobPhase::measuring&&fraction>.65&&fraction<1)
        {auto w=wave->overview();auto l=loudness.view(binding);partial|=w&&!w->complete()&&l&&l->coverage()>0;pending|=!cache.view(binding).accepted;}
    };
    const auto result=OfflineAnalysisService(nullptr,wave.get(),&loudness).analyze(request.context,MeasurementDomain::source,measured);
    expect(result.value&&f->ordered,"feedback-enabled measurement completes with ordered progress within each stage");
    expect(f->phases==std::vector<JobPhase>{JobPhase::copy,JobPhase::copyHash,JobPhase::metadata,JobPhase::samples,JobPhase::sourceHash,JobPhase::measuring,JobPhase::finalising},"preparation precedes measurement and finalisation without skipping verification");
    expect(partial&&pending,"progressive waveform and measured M/S remain visible but unaccepted during measurement");exact(*baseline.value,*result.value);
    expect(wave->overview()->complete()&&!cache.view(binding).accepted&&cache.accept(binding,binding),"EOF preview requires separate completed-job acceptance");
    for(auto phase:{JobPhase::measuring,JobPhase::finalising})
    {
        JobControl c;auto cancel=std::make_shared<RecordingFeedback>();c.feedback=cancel;cancel->cancelControl=&c;cancel->cancelAt=phase;cancel->cancelAfter=.65;
        WaveformCache abandoned;auto w=abandoned.begin(binding);LoudnessProgress l(binding);
        const auto outcome=OfflineAnalysisService(nullptr,w.get(),&l).analyze(request.context,MeasurementDomain::source,c);
        expect(cancel->sawCancel&&outcome.error==ErrorCode::cancelled&&!outcome.value&&!l.view(binding)&&w->state()==WaveformProgress::State::abandoned,"measurement/finalisation cancellation retires both previews without final authority");
    }
}
void flacChecksum(const juce::File& root)
{
    auto request=fixture(root.getChildFile("checksum-source.wav"),48000,2,48013);auto input=request.context.sourceFile;
    juce::AudioFormatManager formats;registerSourceFormats(formats);std::unique_ptr<juce::AudioFormatReader> reader(formats.createReaderFor(input));expect(bool(reader),"open checksum fixture");juce::AudioBuffer<float> pcm(2,48013);reader->read(&pcm,0,48013,0,true,true);
    auto file=root.getChildFile("checksum.flac");auto stream=file.createOutputStream();juce::FlacAudioFormat format;
    std::unique_ptr<juce::AudioFormatWriter> writer(format.createWriterFor(stream.get(),48000,2,24,{},5));expect(bool(writer),"create FLAC fixture");stream.release();expect(writer->writeFromAudioSampleBuffer(pcm,0,48013),"write FLAC fixture");writer.reset();
    JobControl c;auto f=std::make_shared<RecordingFeedback>();c.feedback=f;f->cancelControl=&c;f->cancelAt=JobPhase::checksum;f->cancelAfter=.25;
    const auto result=validateMedia(file,c);expect(f->sawCancel&&result.error==ErrorCode::cancelled,"FLAC checksum progress supports cancellation during decoded checksum pass");
}
void memberImport(const juce::File& root)
{
    const auto a=fixture(root.getChildFile("left.wav"),48000,1,48013),b=fixture(root.getChildFile("right.wav"),48000,1,48013);
    bool completed=false,adopted=false;ProgrammeImportPanel panel([&](auto result){completed=true;adopted=result.has_value();},root.getChildFile("role-memory"));
    panel.inspectFiles({a.context.sourceFile,b.context.sourceFile});const auto control=panel.preparationControl();
    expect(control&&control->feedback,"separate-file import exposes worker feedback to the main window");panel.cancel();
    const auto deadline=juce::Time::getMillisecondCounterHiRes()+5000;
    while(!completed&&juce::Time::getMillisecondCounterHiRes()<deadline)juce::MessageManager::getInstance()->runDispatchLoopUntil(10);
    expect(completed&&!adopted&&control->isCancelled(),"cancelling separate-file validation preserves previous programme and alignment confirmation boundary");
}
void imageOf(juce::Component& c,const juce::File& file)
{auto output=file.createOutputStream();expect(output&&juce::PNGImageFormat().writeImageToStream(c.createComponentSnapshot(c.getLocalBounds()),*output),"save main-window feedback image");}
void waitJob(StandaloneMainComponent& c)
{
    const auto limit=juce::Time::getMillisecondCounterHiRes()+15000;
    while(Access::model(c).job.isActive()&&juce::Time::getMillisecondCounterHiRes()<limit)
    {Access::tick(c);std::this_thread::sleep_for(std::chrono::milliseconds(2));}
    expect(!Access::model(c).job.isActive(),"main-window job finishes without blocking message-thread polling");
}
void shell(const juce::File& root)
{
    StandaloneMainComponent ui(emptyOutputCatalog(),noOutputSettings());Access::stopTimer(ui);auto& s=Access::model(ui);
    auto* stage=dynamic_cast<juce::Label*>(ui.findChildWithID("job-stage"));auto* detail=dynamic_cast<juce::Label*>(ui.findChildWithID("job-detail"));auto* progress=ui.findChildWithID("job-progress");
    expect(stage&&detail&&progress&&stage->getY()<110&&detail->getBottom()<130,"main stage, detail and progress are prominent above comparison rows");
    s.job.begin(JobStage::ingest,s.revision);Access::present(ui,JobPhase::samples,.37);
    expect(stage->getText().contains("source validation")&&detail->getText().contains("37%")&&progress->isVisible(),"main window immediately identifies validation work and pass progress");imageOf(ui,root.getChildFile("import.png"));
    Access::cancel(ui);expect(s.job.state==JobState::cancelling&&stage->getText().contains("Cancelling"),"main-window cancellation feedback is immediate");s.job.finish(s.job.ticket,s.revision,ErrorCode::none,"late result");
    auto request=fixture(root.getChildFile("ui.wav"),48000,2,240013);
    Access::import(ui,request.context.sourceFile);Access::cancel(ui);waitJob(ui);expect(!s.sourceDocument.hasSource&&s.job.state==JobState::cancelled,"cancelled import cannot adopt a late successful validation");
    expect(s.setValidatedSource(request.context.sourceFile,request.context.media),"adopt validated fixture for analysis lifecycle");
    s.job.begin(JobStage::ingest,s.revision);s.job.finish(s.job.ticket,s.revision,ErrorCode::invalidMetadata,"Invalid replacement source");Access::refresh(ui);
    expect(stage->getText()=="Invalid replacement source","import failure remains visible when previous source still needs role confirmation");
    s.sourceDocument.layout=request.context.fileLayout;
    Access::measure(ui);Access::cancel(ui);waitJob(ui);for(int n=0;n<10;++n)Access::tick(ui);
    expect(s.job.state==JobState::cancelled&&!s.source.hasMeasurement&&!s.currentState.hasMeasurement,"cancelled automatic analysis does not restart or grant Current measurement authority");
    s.job.begin(JobStage::analysis,s.revision);Access::present(ui,JobPhase::copy,.42);imageOf(ui,root.getChildFile("preparation.png"));
    Access::present(ui,JobPhase::measuring,.68);imageOf(ui,root.getChildFile("measurement.png"));
    s.job.requestCancel();s.job.finish(s.job.ticket,s.revision,ErrorCode::cancelled,"cancelled rendering fixture");
    Access::measure(ui);waitJob(ui);
    expect(s.source.hasMeasurement&&s.currentState.hasMeasurement&&s.currentState.analysis->identicalSourceMeasurementIdentity==s.source.analysis->identity,"explicit retry completes Source and identical initial Current with provenance");
    expect(Access::stopped(ui),"all import/analysis feedback interactions leave transport disconnected and stopped");
    auto ticket=s.job.ticket;expect(!s.job.finish(ticket,s.revision,ErrorCode::none,"late duplicate"),"late completion cannot replace completed job state");
    s.job.begin(JobStage::analysis,s.revision);Access::present(ui,JobPhase::finalising,-1);expect(detail->getText().contains("Finalising")&&!detail->getText().contains("100%"),"main banner distinguishes finishing work from completed results");
    auto old=s.job.ticket;++s.revision;expect(!s.job.finish(old,s.revision,ErrorCode::none,"stale result")&&s.job.error==ErrorCode::staleRevision,"stale revision completion remains rejected");
}
}
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI gui;if(argc!=2)return 2;const juce::File root(argv[1]);if(root.exists())return 2;root.createDirectory();bool passed=false;juce::String error;
    try{coherentPublication();stagesAndCancellation(root.getChildFile("stereo"),2);stagesAndCancellation(root.getChildFile("surround"),6);flacChecksum(root);memberImport(root);shell(root);passed=true;}
    catch(const std::exception& e){error=e.what();}
    auto data=juce::var(new juce::DynamicObject());auto* r=data.getDynamicObject();r->setProperty("passed",passed);r->setProperty("error",error);r->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);r->setProperty("checks",checks);r->setProperty("physical_devices_opened",false);root.getChildFile("result.json").replaceWithText(juce::JSON::toString(data,true));std::cout<<(passed?"PASS":"FAIL")<<" "<<checks.size()<<" "<<error<<std::endl;return passed?0:1;
}

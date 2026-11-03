// Dedicated combined qualification. Reuses oracle helpers, never the historical runner.
#define main unusedM14QualificationMain
#include "MonoTests.cpp"
#undef main
#include "Standalone/MeasurementTimeline.h"
#include "Standalone/StandaloneMarkerImport.h"
#include "Standalone/StandaloneWaveform.h"
namespace {
int64_t qFrames=2073600000LL;
double qBudgetMinutes=300;
int64_t qImpulse(){return qFrames/2+monoLongRate;}
struct QView:juce::Timer {
    MeasurementTimeline timeline;
    std::shared_ptr<waveform::Budget> budget=std::make_shared<waveform::Budget>();
    std::shared_ptr<WaveformProgress> wave;
    std::shared_ptr<LoudnessProgress> loudness;
    WaveformBinding binding;
    std::shared_ptr<const waveform::Overview> source;
    std::vector<double> preparation,painting;
    uint64_t ticks=0,partial=0;double first=0,start=monoNow();
    explicit QView(int64_t frames){timeline.setSize(1200,600);timeline.setProgramme(frames,monoLongRate,"long-qualification");startTimer(200);}
    void timerCallback() override {
        const auto t=monoNow();if(wave){auto data=wave->overview();timeline.setWaveform(data,data&&data->complete(),{"Centre"});if(data&&data->coveredFrames()>0&&!data->complete()){++partial;if(!first)first=monoNow()-start;}}
        if(loudness)timeline.setPreview(loudness->view(binding));
        // Alternate whole programme and bounded passages throughout the run.
        if(++ticks%25==0){const double at=double((ticks/25*7919)%1000)/1000*qFrames;timeline.setViewRange(at,std::min(double(qFrames),at+monoLongRate*30.));}
        else if(ticks%25==1)timeline.setViewRange(0,double(qFrames));
        preparation.push_back(monoNow()-t);const auto p=monoNow();auto image=timeline.createComponentSnapshot(timeline.getLocalBounds());painting.push_back(monoNow()-p);
        if(preparation.size()>100000){preparation.erase(preparation.begin(),preparation.begin()+50000);painting.erase(painting.begin(),painting.begin()+50000);}
    }
    juce::var report(){stopTimer();auto p95=[](std::vector<double> x){std::sort(x.begin(),x.end());return x.empty()?0.:x[size_t((x.size()-1)*.95)];};
        expect(partial>0&&ticks>10,"long live previews rendered and navigation exercised");expect(p95(preparation)<=50&&p95(painting)<=100,"timeline p95 preparation and paint budgets");expect(budget->used()<=waveform::combinedBudgetBytes,"combined Source Current overview allocation budget");
        auto r=juce::var(new juce::DynamicObject());auto* o=r.getDynamicObject();o->setProperty("updates",juce::String(ticks));o->setProperty("partial_updates",juce::String(partial));o->setProperty("first_waveform_ms",first);o->setProperty("prepare_p95_ms",p95(preparation));o->setProperty("paint_p95_ms",p95(painting));o->setProperty("rendering","actual timeline component offscreen at 5 Hz");return r;}
};
AnalysisSnapshot qMeasure(StandaloneSessionModel& s,MeasurementDomain d,QView& view,bool observe){
    expect(s.job.begin(JobStage::analysis,s.revision),"begin combined measurement");auto c=s.measurementContext(d);c.ticket=s.job.ticket;JobControl control;
    auto binding=WaveformBinding::from(c,d);auto wave=std::make_shared<WaveformProgress>(view.budget);auto loudness=std::make_shared<LoudnessProgress>(binding);
    if(observe){view.wave=wave;view.loudness=loudness;view.binding=binding;view.start=monoNow();}
    auto future=std::async(std::launch::async,[&]{return OfflineAnalysisService(nullptr,observe?wave.get():nullptr,observe?loudness.get():nullptr).analyze(c,d,control);});
    auto result=monoAwait(future);if(!result.value)std::cerr<<result.message<<std::endl;
    expect(result.value.has_value()&&s.job.finish(c.ticket,s.revision,result.error,result.message),"combined measurement finishes");
    if(observe){expect(s.acceptMeasurement(*result.value),"accept independent long measurement");auto data=wave->overview();expect(data&&data->complete(),"complete long PCM overview");
        if(d==MeasurementDomain::source)view.source=data;
        auto preview=loudness->view(binding);expect(preview&&preview->coverage()==qFrames,"long measured M S coverage reaches EOF");
        for(bool st:{false,true})for(const auto& w:st?result.value->detail.shortTerm:result.value->detail.momentary){auto x=preview->inspect(st,w.endFrame,qFrames);if(!x||x->endFrame!=w.endFrame||x->energy!=w.energy)throw std::runtime_error("long progressive loudness differs from completed data");}
    }
    return *result.value;
}
int32_t qCode(int64_t frame)
{
    if(frame==qImpulse())return 6710886;
    static const auto table=[]{std::array<int32_t,192> a{};for(int n=0;n<192;++n)a[size_t(n)]=int32_t(std::llround(.08*8388608*std::sin(n*2*3.14159265358979323846/192)));return a;}();
    return table[size_t(frame%192)];
}
StandaloneSessionModel qMake(const juce::File& root)
{
    const auto file=root.getChildFile("mono-duration-192k.wav");auto out=file.createOutputStream();expect(out!=nullptr,"owned sustained mono source stream");

    pcm::Format format;juce::String reason;
    expect(pcm::describe(monoLongRate,1,qFrames,format,reason,4),"bounded exact RF64 format");
    const auto header=pcm::header(format);expect(out->write(header.data(),header.size()),"RF64 ds64 header");
    constexpr int block=192*2048;std::vector<uint8_t> packed(size_t(block)*3);
    for(int n=0;n<block;++n){const uint32_t code=uint32_t(qCode(n));for(int b=0;b<3;++b)packed[size_t(n)*3+size_t(b)]=uint8_t(code>>(8*b));}
    for(int64_t first=0;first<qFrames;first+=block)
    {
        const int count=int(std::min<int64_t>(block,qFrames-first));
        const bool impulse=first<=qImpulse()&&qImpulse()<first+count;
        if(impulse)for(int b=0;b<3;++b)packed[size_t(qImpulse()-first)*3+size_t(b)]=uint8_t(uint32_t(qCode(qImpulse()))>>(8*b));
        if(!out->write(packed.data(),size_t(count)*3))throw std::runtime_error("Sustained source write failed");
        if(impulse)for(int b=0;b<3;++b)packed[size_t(qImpulse()-first)*3+size_t(b)]=uint8_t(uint32_t(qCode(qImpulse()%192))>>(8*b));
    }
    out->flush();expect(out->getStatus().wasOk(),"sustained source finalized");out.reset();JobControl job;auto v=validateMedia(file,job);
    StandaloneSessionModel s;expect(v.succeeded()&&v.media.frames==qFrames&&s.setValidatedSource(file,v.media),"complete sustained mono source validated");
    s.setChannelRole(0,juce::AudioChannelSet::centre);expect(s.confirmLayout()&&s.initializeCurrentStateFromSource(),"sustained explicit mono interpretation");
    s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);return s;
}
void qLong(const juce::File& root)
{
    const auto begin=monoNow(),initialFree=double(monoFree(root));expect(initialFree>=40.*1024*1024*1024,"40 GiB free before one sustained workflow");
    MonoPulse pulse;pulse.startTimer(10);auto stages=juce::var(new juce::DynamicObject());
    const auto stage=[&](const char* name,double start){stages.getDynamicObject()->setProperty(name,monoNow()-start);std::cout<<"M14 LONG "<<name<<" elapsed_s="<<(monoNow()-begin)/1000<<std::endl;};
    double t=monoNow();auto generating=std::async(std::launch::async,[&]{return qMake(root);});auto s=monoAwait(generating);stage("generate_import_ms",t);
    const auto baseline=monoRss();std::atomic<uint64_t> peak{baseline};std::atomic<int64_t> minimumFree{monoFree(root)};std::atomic<bool> stop{false};
    std::thread sampler([&]{while(!stop){peak=std::max(peak.load(),monoRss());minimumFree=std::min(minimumFree.load(),monoFree(root));std::this_thread::sleep_for(std::chrono::milliseconds(10));}});
    struct Join{std::atomic<bool>& stop;std::thread& worker;~Join(){stop=true;worker.join();}}join{stop,sampler};
    QView view(qFrames);t=monoNow();auto plain=qMeasure(s,MeasurementDomain::source,view,false);const auto plainMs=monoNow()-t;stage("plain_source_ms",t);
    const auto plainPeak=peak.load();peak=monoRss();t=monoNow();auto observed=qMeasure(s,MeasurementDomain::source,view,true);const auto observedMs=monoNow()-t;stage("observed_source_ms",t);
    // Input identity is not a numerical result digest. ExactLongSource.cpp performs
    // the complete metric/window comparison against the saved observed cache.
    expect(plain.identity==observed.identity&&plain.measuredContentHash==observed.measuredContentHash,"observed and plain three-hour input identities and content hashes match");
    expect(plain.detail.samplePeakFrames==observed.detail.samplePeakFrames&&plain.detail.truePeakEighthFrames==observed.detail.truePeakEighthFrames,"all exact long peak positions unchanged");
    expect(observedMs<=1.15*plainMs,"long observation overhead <=15 percent");
    const auto observedPeak=peak.load();expect(observedPeak<=plainPeak+256ULL*1024*1024,"incremental waveform peak RSS <=256 MiB above plain baseline");
    expect(s.currentState.hasMeasurement&&s.currentState.analysis->identicalSourceMeasurementIdentity==s.source.analysis->identity,"identical Current automatically initialized from completed Source");
    plain={};observed={};
    s.synchronizeSections();juce::String error;std::string markerBytes;
    for(int n=1;n<1000;++n){if(n>1)markerBytes+=",";markerBytes+=std::to_string(qFrames*n/1000);}
    proposal::markers::Options options;options.format=proposal::markers::Format::pastedPositions;options.unit=proposal::markers::Unit::audioFrames;options.confirmedUnits=options.confirmedPointMarkers=true;
    const auto imported=proposal::markers::prepare(markerBytes,"Qualification positions",options,s.sectionDocument.map(),s.sessionId);
    expect(commitMarkerImport(s,imported,proposal::markers::contentHash(markerBytes),error)&&s.sectionDocument.map().sections.size()==1000&&s.markerImport,"1000 long sections retain explicit marker provenance");

    JobControl job,proposalJob;const auto request=s.sectionProposalRequest();t=monoNow();auto proposing=std::async(std::launch::async,[&]{return generateSectionProposal(request,proposalJob);});auto proposed=monoAwait(proposing);stage("proposal_ms",t);
    if(proposed.value)std::cout<<(*proposed.value)->text()<<std::endl;else std::cerr<<proposed.message<<std::endl;
    expect(proposed.value&&(*proposed.value)->status()==SectionProposalStatus::ready&&(*proposed.value)->passes()<=3,"sustained mono exact candidate within three passes");
    const auto automatic=*proposed.value;expect(s.acceptSectionProposal(automatic)&&s.applySectionProposal()&&!s.currentState.hasMeasurement,"sustained explicit Apply leaves Current unverified");
    const auto planKey=s.currentStateDocument.sectionGains->audioKey();
    t=monoNow();qMeasure(s,MeasurementDomain::currentState,view,true);stage("independent_current_ms",t);
    expect(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"independently measured sustained Current passes selected rules");
    t=monoNow();auto analysis=std::async(std::launch::async,[&]{auto p=proposal::analyzeProgramme(s.currentState.analysis->detail,s.currentState.analysis->identity);return proposal::analyzeSections(*p,s.sectionDocument.map(),s.sectionDocument.map().inputIdentity);});
    expect(monoAwait(analysis)->sections.size()==1000,"sustained full-resolution diagnostics");stage("diagnostics_ms",t);
    expect(s.lastAppliedExplanation&&s.lastAppliedExplanation->candidate,"full measured section explanation retained");
    const auto history=s.lastAppliedExplanation->text();
    const auto doc=sf::capture(s);const auto session=root.getChildFile("mono-duration.levelscope");
    t=monoNow();auto saving=std::async(std::launch::async,[&]{return sf::save(doc,session,job);});expect(monoAwait(saving).succeeded(),"sustained coherent session save");const auto saveMs=monoNow()-t;stage("save_ms",t);
    t=monoNow();auto parsing=std::async(std::launch::async,[&]{return sf::read(session,job);});auto read=monoAwait(parsing);const auto parseMs=monoNow()-t;stage("parse_ms",t);
    expect(read.document&&read.document->current.sectionGains->audioKey()==planKey&&read.document->sections.sections.size()==1000,"sustained exact plan and annotation roundtrip");
    expect(read.document->lastAppliedExplanation&&read.document->lastAppliedExplanation->text()==history&&read.document->markerImport&&read.document->markerImport->identity()==s.markerImport->identity(),"complete long history and marker receipt survive save/parse");
    t=monoNow();auto opening=std::async(std::launch::async,[&]{return sf::verify(*read.document,job);});auto reopened=monoAwait(opening);expect(reopened.prepared&&reopened.prepared->model.source.hasMeasurement&&reopened.prepared->model.currentState.hasMeasurement,"sustained source and caches reverified on reopening");stage("reopen_ms",t);expect(reopened.prepared->model.lastAppliedExplanation&&reopened.prepared->model.lastAppliedExplanation->text()==history,"verified reopen preserves complete historical explanation");read={};reopened={};
    auto historical=s;auto requirements=historical.jobRequirements;requirements.maximumLraTenths=180;
    expect(historical.setJobRequirements(requirements,error)&&historical.appliedExplanationStatus().contains("Requirements changed"),"changed long requirements mark history separately");
    expect(historical.setCurrentGainDb(historical.currentStateDocument.staticGainDb-.1)&&!historical.currentState.hasMeasurement&&historical.lastAppliedExplanation->text()==history,"later processing leaves long Current unverified and history intact");
    expect(sf::save(sf::capture(historical),root.getChildFile("historical-unverified.levelscope"),job).succeeded(),"save independent stale-history case");
    JobControl cancel;t=monoNow();auto preparing=std::async(std::launch::async,[&]{return prepareVerifiedSource(s.measurementContext(MeasurementDomain::currentState),cancel);});
    while(cancel.progress.load()==0&&preparing.wait_for(std::chrono::milliseconds(1))!=std::future_status::ready)juce::MessageManager::getInstance()->runDispatchLoopUntil(10);
    t=monoNow();cancel.cancel=true;const auto cancelled=monoAwait(preparing);const auto cancelMs=monoNow()-t;expect(cancelled.error==ErrorCode::cancelled&&cancelMs<=1000,"sustained owned-copy cancellation within one second");
    auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));OutputRouting muted;
    s.selectAuditionDomain(MeasurementDomain::currentState);auto auditionRequest=s.auditionRequest();audition.synchronize(auditionRequest,muted);
    t=monoNow();expect(audition.prepare(auditionRequest).isEmpty(),"sustained device-closed preparation");const auto deadline=monoNow()+1800000;
    while(audition.preparing()&&monoNow()<deadline){audition.poll();juce::MessageManager::getInstance()->runDispatchLoopUntil(10);}
    expect(audition.preparedCopy()&&!injected->opened,"sustained copy prepared without physical device");const auto* copy=audition.preparedCopy();
    expect(audition.connect(monoLongRate).isEmpty(),"sustained injected device only");
    for(auto domain:{MeasurementDomain::source,MeasurementDomain::currentState})for(int64_t first:{int64_t(0),qFrames/1000-2048,qImpulse()-1500,qFrames-4096})
    {
        s.selectAuditionDomain(domain);auditionRequest=s.auditionRequest();audition.synchronize(auditionRequest,muted);expect(audition.seek(first,auditionRequest,muted).isEmpty(),"sustained exact source-frame seek");
        Capture captured(1,4096,first);audition.setCaptureObserver(&captured);expect(audition.play(auditionRequest,muted,false,true).isEmpty(),"sustained explicit injected playback");awaitPlaying(audition,*injected);juce::AudioBuffer<float> output(maxChannels,256);
        for(int n=0;n<16;++n){injected->tick(output);audition.poll();std::this_thread::sleep_for(std::chrono::milliseconds(1));}audition.stop();audition.setCaptureObserver(nullptr);
        juce::AudioBuffer<float> raw(1,4096);for(int n=0;n<4096;++n)raw.setSample(0,n,float(qCode(first+n))/8388608.f);
        expect(captured.valid&&captured.next==first+4096&&sameSamples(captured.data,monoOracle(raw,auditionRequest.context.processing,first))&&audition.preparedCopy()==copy&&audition.underruns()==0,"sustained Source Current captures exact at first transition impulse and EOF with copy reuse");
    }
    audition.disconnect();audition.sourceChanged({});stage("audition_capture_ms",t);
    const auto exportRequest=ex::capture(s);t=monoNow();auto rendering=std::async(std::launch::async,[&]{return ex::render(exportRequest,root,job);});auto rendered=monoAwait(rendering);stage("render_verify_ms",t);
    if(!rendered.value)std::cerr<<rendered.message<<'\n';expect(rendered.value.has_value(),"sustained mono render reopen and actual artifact measurement");
    const auto artifact=*rendered.value;expect(artifact->spec().format.rf64==(qFrames==2073600000LL)&&artifact->spec().format.frames==qFrames&&artifact->file().getSize()==int64_t(artifact->spec().format.fileBytes)&&ex::artifactEvaluation(*artifact,exportRequest).verdict==proposal::Verdict::pass,"sustained artifact has the required RIFF/RF64 container with exact frame count and measured pass");
    const auto destination=root.getChildFile("export");t=monoNow();auto publishing=std::async(std::launch::async,[&]{return ex::publish(*artifact,exportRequest,destination,ex::approvalIdentity(*artifact,exportRequest),job);});expect(monoAwait(publishing).value.has_value(),"explicit local verified sustained export publication");stage("publish_ms",t);
    pulse.stopTimer();std::sort(pulse.gaps.begin(),pulse.gaps.end());expect(pulse.gaps.size()>10,"sustained heartbeat sampled");const auto p95=pulse.gaps[size_t((pulse.gaps.size()-1)*.95)],growth=double(peak.load()-baseline);const auto reportBytes=destination.getChildFile("report.json").getSize()+destination.getChildFile("report.html").getSize();
    expect(saveMs<=5000&&parseMs<=5000&&session.getSize()<=128LL*1024*1024,"sustained session resource budgets");expect(p95<=100&&reportBytes<=16*1024*1024,"sustained RSS responsiveness and report budgets");expect(monoNow()-begin<=qBudgetMinutes*60000.,"one sustained workflow completes within host execution budget");
    longEvidence=juce::var(new juce::DynamicObject());auto* e=longEvidence.getDynamicObject();e->setProperty("frames",juce::String(qFrames));e->setProperty("rate",monoLongRate);e->setProperty("channels",1);e->setProperty("intervals",1000);e->setProperty("capture_cases",8);e->setProperty("private_passes",automatic->passes());e->setProperty("stages_ms",stages);e->setProperty("total_ms",monoNow()-begin);e->setProperty("cancel_ms",cancelMs);e->setProperty("ui_p95_ms",p95);e->setProperty("additional_peak_rss_bytes",growth);e->setProperty("disk_growth_bytes",initialFree-double(minimumFree.load()));e->setProperty("session_bytes",session.getSize());e->setProperty("report_bytes",reportBytes);e->setProperty("artifact_sha256",artifact->hash());e->setProperty("plain_source_ms",plainMs);e->setProperty("observed_source_ms",observedMs);e->setProperty("plain_peak_rss",juce::String(plainPeak));e->setProperty("observed_peak_rss",juce::String(observedPeak));e->setProperty("waveform_retained_bytes",juce::String(view.budget->used()));e->setProperty("timeline",view.report());
    e->setProperty("three_hour_qualified",qFrames==2073600000LL);e->setProperty("rf64_end_to_end_qualified",artifact->spec().format.rf64);
}
}
int main(int argc,char** argv){
    if(argc!=4)return 2;
    juce::ScopedJuceInitialiser_GUI gui;const juce::File root(argv[1]);qFrames=int64_t(std::stoll(argv[2]))*monoLongRate;qBudgetMinutes=std::stod(argv[3]);bool passed=false;juce::String error;
    try{expect(!root.exists()&&root.createDirectory().wasOk(),"fresh combined qualification root");expect(qFrames>=monoLongRate*10LL&&qFrames<=2073600000LL,"bounded requested qualification duration");qLong(root);passed=true;}
    catch(const std::exception& e){error=e.what();std::cerr<<error<<std::endl;}
    auto report=juce::var(new juce::DynamicObject());auto* r=report.getDynamicObject();r->setProperty("passed",passed);r->setProperty("failure",error);r->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);r->setProperty("checks",checks);r->setProperty("long_workflow",longEvidence);r->setProperty("physical_devices_opened",false);r->setProperty("os",juce::SystemStats::getOperatingSystemName());root.getChildFile("result.json").replaceWithText(juce::JSON::toString(report,true));return passed?0:1;
}

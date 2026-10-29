#pragma once
#include <sys/statvfs.h>
namespace
{
constexpr int monoLongRate=192000;
constexpr int64_t monoLongFrames=int64_t(monoLongRate)*90*60;
constexpr int64_t monoImpulse=monoLongFrames/2+monoLongRate;
double monoNow(){return juce::Time::getMillisecondCounterHiRes();}
uint64_t monoRss(){mach_task_basic_info_data_t i{};mach_msg_type_number_t n=MACH_TASK_BASIC_INFO_COUNT;return task_info(mach_task_self(),MACH_TASK_BASIC_INFO,reinterpret_cast<task_info_t>(&i),&n)==KERN_SUCCESS?i.resident_size:0;}
int64_t monoFree(const juce::File& f){struct statvfs v{};return statvfs(f.getFullPathName().toRawUTF8(),&v)==0?int64_t(v.f_bavail)*v.f_frsize:0;}
template<class T>T monoAwait(std::future<T>& pending)
{while(pending.wait_for(std::chrono::milliseconds(0))!=std::future_status::ready)juce::MessageManager::getInstance()->runDispatchLoopUntil(10);return pending.get();}
struct MonoPulse:juce::Timer
{double last=monoNow();std::vector<double> gaps;void timerCallback()override{const auto now=monoNow();gaps.push_back(now-last);last=now;}};
int32_t monoLongCode(int64_t frame)
{
    if(frame==monoImpulse)return 6710886;
    static const auto table=[]{std::array<int32_t,192> a{};for(int n=0;n<192;++n)a[size_t(n)]=int32_t(std::llround(.08*8388608*std::sin(n*2*3.14159265358979323846/192)));return a;}();
    return table[size_t(frame%192)];
}
StandaloneSessionModel makeLongMono(const juce::File& root)
{
    const auto file=root.getChildFile("mono-90min-192k.wav");auto out=file.createOutputStream();expect(out!=nullptr,"owned sustained mono source stream");
    const uint64_t bytes=uint64_t(monoLongFrames)*3;
    out->write("RIFF",4);out->writeInt(int(bytes+36));out->write("WAVEfmt ",8);out->writeInt(16);out->writeShort(1);out->writeShort(1);
    out->writeInt(monoLongRate);out->writeInt(monoLongRate*3);out->writeShort(3);out->writeShort(24);out->write("data",4);out->writeInt(int(bytes));
    constexpr int block=192*2048;std::vector<uint8_t> packed(size_t(block)*3);
    for(int n=0;n<block;++n){const uint32_t code=uint32_t(monoLongCode(n));for(int b=0;b<3;++b)packed[size_t(n)*3+size_t(b)]=uint8_t(code>>(8*b));}
    for(int64_t first=0;first<monoLongFrames;first+=block)
    {
        const int count=int(std::min<int64_t>(block,monoLongFrames-first));
        const bool impulse=first<=monoImpulse&&monoImpulse<first+count;
        if(impulse)for(int b=0;b<3;++b)packed[size_t(monoImpulse-first)*3+size_t(b)]=uint8_t(uint32_t(monoLongCode(monoImpulse))>>(8*b));
        if(!out->write(packed.data(),size_t(count)*3))throw std::runtime_error("Sustained source write failed");
        if(impulse)for(int b=0;b<3;++b)packed[size_t(monoImpulse-first)*3+size_t(b)]=uint8_t(uint32_t(monoLongCode(monoImpulse%192))>>(8*b));
    }
    out->flush();expect(out->getStatus().wasOk(),"sustained source finalized");out.reset();JobControl job;auto v=validateMedia(file,job);
    StandaloneSessionModel s;expect(v.succeeded()&&v.media.frames==monoLongFrames&&s.setValidatedSource(file,v.media),"complete 90-minute mono source validated");
    s.setChannelRole(0,juce::AudioChannelSet::centre);expect(s.confirmLayout()&&s.initializeCurrentStateFromSource(),"sustained explicit mono interpretation");
    s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);return s;
}
void monoLong(const juce::File& root)
{
    const auto begin=monoNow(),initialFree=double(monoFree(root));expect(initialFree>=20.*1024*1024*1024,"20 GiB free before one 90-minute workflow");
    MonoPulse pulse;pulse.startTimer(10);auto stages=juce::var(new juce::DynamicObject());
    const auto stage=[&](const char* name,double start){stages.getDynamicObject()->setProperty(name,monoNow()-start);std::cout<<"M14 LONG "<<name<<" elapsed_s="<<(monoNow()-begin)/1000<<std::endl;};
    double t=monoNow();auto generating=std::async(std::launch::async,[&]{return makeLongMono(root);});auto s=monoAwait(generating);stage("generate_import_ms",t);
    const auto baseline=monoRss();std::atomic<uint64_t> peak{baseline};std::atomic<int64_t> minimumFree{monoFree(root)};std::atomic<bool> stop{false};
    std::thread sampler([&]{while(!stop){peak=std::max(peak.load(),monoRss());minimumFree=std::min(minimumFree.load(),monoFree(root));std::this_thread::sleep_for(std::chrono::milliseconds(10));}});
    struct Join{std::atomic<bool>& stop;std::thread& worker;~Join(){stop=true;worker.join();}}join{stop,sampler};
    t=monoNow();auto source=std::async(std::launch::async,[&]{return measuredComposition(s,MeasurementDomain::source);});monoAwait(source);stage("source_measure_ms",t);
    s.synchronizeSections();juce::String error;for(int n=1;n<1000;++n)expect(s.sectionDocument.split(uint64_t(n),monoLongFrames*n/1000,error),"sustained manual annotations");s.documentChanged();
    t=monoNow();auto current=std::async(std::launch::async,[&]{return measuredComposition(s,MeasurementDomain::currentState);});monoAwait(current);stage("baseline_current_measure_ms",t);
    JobControl job,proposalJob;const auto request=s.sectionProposalRequest();t=monoNow();auto proposing=std::async(std::launch::async,[&]{return generateSectionProposal(request,proposalJob);});auto proposed=monoAwait(proposing);stage("proposal_ms",t);
    if(proposed.value)std::cout<<(*proposed.value)->text()<<std::endl;else std::cerr<<proposed.message<<std::endl;
    expect(proposed.value&&(*proposed.value)->status()==SectionProposalStatus::ready&&(*proposed.value)->passes()<=3,"sustained mono exact candidate within three passes");
    const auto automatic=*proposed.value;expect(s.acceptSectionProposal(automatic)&&s.applySectionProposal()&&!s.currentState.hasMeasurement,"sustained explicit Apply leaves Current unverified");
    const auto planKey=s.currentStateDocument.sectionGains->audioKey();
    t=monoNow();auto verification=std::async(std::launch::async,[&]{return measuredComposition(s,MeasurementDomain::currentState);});monoAwait(verification);stage("independent_current_ms",t);
    expect(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"independently measured 90-minute Current passes selected rules");
    t=monoNow();auto analysis=std::async(std::launch::async,[&]{auto p=proposal::analyzeProgramme(s.currentState.analysis->detail,s.currentState.analysis->identity);return proposal::analyzeSections(*p,s.sectionDocument.map(),s.sectionDocument.map().inputIdentity);});
    expect(monoAwait(analysis)->sections.size()==1000,"90-minute full-resolution diagnostics");stage("diagnostics_ms",t);
    const auto doc=sf::capture(s);const auto session=root.getChildFile("mono-90min.levelscope");
    t=monoNow();auto saving=std::async(std::launch::async,[&]{return sf::save(doc,session,job);});expect(monoAwait(saving).succeeded(),"sustained coherent session save");const auto saveMs=monoNow()-t;stage("save_ms",t);
    t=monoNow();auto parsing=std::async(std::launch::async,[&]{return sf::read(session,job);});auto read=monoAwait(parsing);const auto parseMs=monoNow()-t;stage("parse_ms",t);
    expect(read.document&&read.document->current.sectionGains->audioKey()==planKey&&read.document->sections.sections.size()==1000,"sustained exact plan and annotation roundtrip");
    t=monoNow();auto opening=std::async(std::launch::async,[&]{return sf::verify(*read.document,job);});auto reopened=monoAwait(opening);expect(reopened.prepared&&reopened.prepared->model.source.hasMeasurement&&reopened.prepared->model.currentState.hasMeasurement,"sustained source and caches reverified on reopening");stage("reopen_ms",t);read={};reopened={};
    JobControl cancel;t=monoNow();auto preparing=std::async(std::launch::async,[&]{return prepareVerifiedSource(s.measurementContext(MeasurementDomain::currentState),cancel);});
    while(cancel.progress.load()==0&&preparing.wait_for(std::chrono::milliseconds(1))!=std::future_status::ready)juce::MessageManager::getInstance()->runDispatchLoopUntil(10);
    t=monoNow();cancel.cancel=true;const auto cancelled=monoAwait(preparing);const auto cancelMs=monoNow()-t;expect(cancelled.error==ErrorCode::cancelled&&cancelMs<=1000,"sustained owned-copy cancellation within one second");
    auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));OutputRouting muted;
    s.selectAuditionDomain(MeasurementDomain::currentState);auto auditionRequest=s.auditionRequest();audition.synchronize(auditionRequest,muted);
    t=monoNow();expect(audition.prepare(auditionRequest).isEmpty(),"sustained device-closed preparation");const auto deadline=monoNow()+1800000;
    while(audition.preparing()&&monoNow()<deadline){audition.poll();juce::MessageManager::getInstance()->runDispatchLoopUntil(10);}
    expect(audition.preparedCopy()&&!injected->opened,"90-minute copy prepared without physical device");const auto* copy=audition.preparedCopy();
    expect(audition.connect(monoLongRate).isEmpty(),"sustained injected device only");
    for(auto domain:{MeasurementDomain::source,MeasurementDomain::currentState})for(int64_t first:{int64_t(0),monoLongFrames/1000-2048,monoImpulse-1500,monoLongFrames-4096})
    {
        s.selectAuditionDomain(domain);auditionRequest=s.auditionRequest();audition.synchronize(auditionRequest,muted);expect(audition.seek(first,auditionRequest,muted).isEmpty(),"sustained exact source-frame seek");
        Capture captured(1,4096,first);audition.setCaptureObserver(&captured);expect(audition.play(auditionRequest,muted,false,true).isEmpty(),"sustained explicit injected playback");awaitPlaying(audition,*injected);juce::AudioBuffer<float> output(maxChannels,256);
        for(int n=0;n<16;++n){injected->tick(output);audition.poll();std::this_thread::sleep_for(std::chrono::milliseconds(1));}audition.stop();audition.setCaptureObserver(nullptr);
        juce::AudioBuffer<float> raw(1,4096);for(int n=0;n<4096;++n)raw.setSample(0,n,float(monoLongCode(first+n))/8388608.f);
        expect(captured.valid&&captured.next==first+4096&&sameSamples(captured.data,monoOracle(raw,auditionRequest.context.processing,first))&&audition.preparedCopy()==copy&&audition.underruns()==0,"sustained Source Current captures exact at first transition impulse and EOF with copy reuse");
    }
    audition.disconnect();audition.sourceChanged({});stage("audition_capture_ms",t);
    const auto exportRequest=ex::capture(s);t=monoNow();auto rendering=std::async(std::launch::async,[&]{return ex::render(exportRequest,root,job);});auto rendered=monoAwait(rendering);stage("render_verify_ms",t);
    if(!rendered.value)std::cerr<<rendered.message<<'\n';expect(rendered.value.has_value(),"90-minute mono render reopen and actual artifact measurement");
    const auto artifact=*rendered.value;expect(!artifact->spec().format.rf64&&artifact->file().getSize()==3110400068LL&&ex::artifactEvaluation(*artifact,exportRequest).verdict==proposal::Verdict::pass,"90-minute artifact is ordinary RIFF with exact frame count and measured pass");
    const auto destination=root.getChildFile("export");t=monoNow();auto publishing=std::async(std::launch::async,[&]{return ex::publish(*artifact,exportRequest,destination,ex::approvalIdentity(*artifact,exportRequest),job);});expect(monoAwait(publishing).value.has_value(),"explicit local verified sustained export publication");stage("publish_ms",t);
    pulse.stopTimer();std::sort(pulse.gaps.begin(),pulse.gaps.end());expect(pulse.gaps.size()>10,"sustained heartbeat sampled");const auto p95=pulse.gaps[size_t((pulse.gaps.size()-1)*.95)],growth=double(peak.load()-baseline);const auto reportBytes=destination.getChildFile("report.json").getSize()+destination.getChildFile("report.html").getSize();
    expect(saveMs<=5000&&parseMs<=5000&&session.getSize()<=128LL*1024*1024,"sustained session resource budgets");expect(growth<=256.*1024*1024&&p95<=100&&reportBytes<=16*1024*1024,"sustained RSS responsiveness and report budgets");expect(monoNow()-begin<=120*60000.,"one 90-minute workflow completes within 120-minute execution budget");
    longEvidence=juce::var(new juce::DynamicObject());auto* e=longEvidence.getDynamicObject();e->setProperty("frames",juce::String(monoLongFrames));e->setProperty("rate",monoLongRate);e->setProperty("channels",1);e->setProperty("intervals",1000);e->setProperty("capture_cases",8);e->setProperty("private_passes",automatic->passes());e->setProperty("stages_ms",stages);e->setProperty("total_ms",monoNow()-begin);e->setProperty("cancel_ms",cancelMs);e->setProperty("ui_p95_ms",p95);e->setProperty("additional_peak_rss_bytes",growth);e->setProperty("disk_growth_bytes",initialFree-double(minimumFree.load()));e->setProperty("session_bytes",session.getSize());e->setProperty("report_bytes",reportBytes);e->setProperty("artifact_sha256",artifact->hash());e->setProperty("three_hour_qualified",false);e->setProperty("rf64_end_to_end_qualified",false);
}
}

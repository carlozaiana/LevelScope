#pragma once
#include "CompositionLongTests.h" // resource/fixture helpers only; never call the historical workflow
namespace
{
void longSectionProposal(const juce::File& root,const juce::File& originalSession)
{
    const auto begin=compositionNow();const auto initialFree=freeBytes(root);
    expect(initialFree>=int64_t(65)*1024*1024*1024,"65 GiB available before single M13 long workflow");
    JobControl job;auto historical=sf::read(originalSession,job);
    expect(historical.document&&historical.document->source.composition&&historical.document->reference.file==juce::File(),"reuse immutable M11 long input references with no video attached");
    auto recipe=*historical.document->source.composition;historical={};
    expect(recipe.summary().frames==longFrames&&recipe.summary().sampleRate==48000&&recipe.summary().channels==6,"three-hour 48 kHz 5.1 composition workload");
    // Reuse original non-LFE members. Two owned replacement LFE files add one
    // finite over-full-scale composed impulse; never patch historical media.
    constexpr int64_t impulseFrame=longFrames/2+48000;
    for(auto* component:{&recipe.me,&recipe.dialogue})
    {
        auto group=*component->assembly;
        auto& member=group.members[3];const auto file=root.getChildFile(component->id+"-m13-lfe.wav");
        expect(member.file.copyFileTo(file),"copy only owned long LFE fixture");
        auto stream=file.createOutputStream();expect(stream&&stream->setPosition(44+impulseFrame*3),"seek owned PCM24 impulse");
        const uint8_t value[3]={0,0,96};expect(stream->write(value,3),"write exact +0.75 component impulse");stream->flush();expect(stream->getStatus().wasOk(),"owned impulse finalized");stream.reset();
        auto inspected=validateMedia(file,job);expect(inspected.succeeded(),"modified long LFE member content reverified");member.file=file;member.media=inspected.media;
        expect(group.confirmAlignment(),"generator reconfirms changed long assembly identity");component->assembly=group;component->media=group.summary();component->layout=group.layout();
    }
    expect(recipe.confirmComplementaryBalancedAlignment(),"generator confirms changed complementary long composition");
    auto s=compositionModel(recipe);CompositionPulse pulse;pulse.startTimer(10);const auto baseline=compositionRss();
    std::atomic<uint64_t> peak{baseline};std::atomic<int64_t> minimumFree{initialFree};std::atomic<bool> stop{false};
    std::thread sampler([&]{while(!stop){peak=std::max(peak.load(),compositionRss());minimumFree=std::min(minimumFree.load(),freeBytes(root));std::this_thread::sleep_for(std::chrono::milliseconds(10));}});
    struct Join{std::atomic<bool>& stop;std::thread& thread;~Join(){stop=true;thread.join();}}join{stop,sampler};
    auto stages=juce::var(new juce::DynamicObject());auto stage=[&](const char* name,double start){stages.getDynamicObject()->setProperty(name,compositionNow()-start);std::cout<<"M13 LONG "<<name<<" elapsed_s="<<(compositionNow()-begin)/1000<<std::endl;};
    double t=compositionNow();auto source=std::async(std::launch::async,[&]{return measuredComposition(s,MeasurementDomain::source);});compositionAwait(source);stage("source_measure_ms",t);
    s.synchronizeSections();juce::String error;for(int n=1;n<1000;++n)expect(s.sectionDocument.split(n,longFrames*n/1000,error),"long manual section boundaries");s.documentChanged();
    t=compositionNow();auto current=std::async(std::launch::async,[&]{return measuredComposition(s,MeasurementDomain::currentState);});compositionAwait(current);stage("baseline_current_measure_ms",t);
    const auto requestProposal=s.sectionProposalRequest();JobControl proposalJob;
    t=compositionNow();auto proposing=std::async(std::launch::async,[&]{return generateSectionProposal(requestProposal,proposalJob);});auto proposalResult=compositionAwait(proposing);stage("private_proposal_ms",t);
    if(proposalResult.value)std::cout<<(*proposalResult.value)->text()<<std::endl;else std::cerr<<proposalResult.message<<std::endl;
    expect(proposalResult.value&&(*proposalResult.value)->status()==SectionProposalStatus::ready&&(*proposalResult.value)->passes()==3,"long bounded local proposal completes three private passes");
    const auto automatic=*proposalResult.value;expect(s.acceptSectionProposal(automatic)&&s.applySectionProposal(),"explicit long proposal Apply");
    const auto key=s.currentStateDocument.sectionGains->audioKey();
    t=compositionNow();auto verification=std::async(std::launch::async,[&]{return measuredComposition(s,MeasurementDomain::currentState);});compositionAwait(verification);stage("independent_current_verification_ms",t);
    expect(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"actual long processed Current passes selected rules");
    t=compositionNow();auto diagnostics=std::async(std::launch::async,[&]{auto a=proposal::analyzeProgramme(s.currentState.analysis->detail,s.currentState.analysis->identity);return proposal::analyzeSections(*a,s.sectionDocument.map(),s.sectionDocument.map().inputIdentity);});
    expect(compositionAwait(diagnostics)->sections.size()==1000,"long full-resolution Current section diagnostics");stage("diagnostics_ms",t);
    const auto doc=sf::capture(s);const auto session=root.getChildFile("section.levelscope");
    t=compositionNow();auto saving=std::async(std::launch::async,[&]{return sf::save(doc,session,job);});expect(compositionAwait(saving).succeeded(),"long coherent schema4 save");const auto saveMs=compositionNow()-t;stage("save_ms",t);
    t=compositionNow();auto reading=std::async(std::launch::async,[&]{return sf::read(session,job);});auto read=compositionAwait(reading);const auto parseMs=compositionNow()-t;stage("parse_ms",t);
    expect(read.document&&read.document->current.sectionGains->audioKey()==key&&read.document->sections.sections.size()==1000,"long exact sections and applied plan recovered");
    t=compositionNow();auto opening=std::async(std::launch::async,[&]{return sf::verify(*read.document,job);});auto reopened=compositionAwait(opening);
    expect(reopened.prepared&&reopened.prepared->model.currentState.hasMeasurement&&reopened.prepared->model.source.hasMeasurement,"long input reverified and same-build caches validated");stage("reopen_verify_ms",t);reopened={};read={};
    JobControl cancel;auto preparing=std::async(std::launch::async,[&]{return prepareVerifiedSource(s.measurementContext(MeasurementDomain::currentState),cancel);});
    while(cancel.progress.load()==0&&preparing.wait_for(std::chrono::milliseconds(1))!=std::future_status::ready)juce::MessageManager::getInstance()->runDispatchLoopUntil(10);
    t=compositionNow();cancel.cancel=true;const auto cancelled=compositionAwait(preparing);const auto cancelMs=compositionNow()-t;expect(cancelled.error==ErrorCode::cancelled&&cancelMs<=1000,"long owned-source cancellation <=1 second");
    auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));OutputRouting muted;
    s.selectAuditionDomain(MeasurementDomain::currentState);auto request=s.auditionRequest();audition.synchronize(request,muted);
    t=compositionNow();expect(audition.prepare(request).isEmpty(),"long device-closed preparation");const auto deadline=compositionNow()+1800000;
    while(audition.preparing()&&compositionNow()<deadline){audition.poll();juce::MessageManager::getInstance()->runDispatchLoopUntil(10);}
    expect(audition.preparedCopy()&&!injected->opened,"complete verified copies without physical output");const auto* copy=audition.preparedCopy();
    expect(audition.connect(48000).isEmpty(),"long audio-only injected device");
    for(auto domain:{MeasurementDomain::source,MeasurementDomain::currentState})for(int64_t first:{int64_t(0),longFrames/1000-2048,impulseFrame-1500,longFrames-4096})
    {
        s.selectAuditionDomain(domain);request=s.auditionRequest();audition.synchronize(request,muted);expect(audition.seek(first,request,muted).isEmpty(),"long exact first transition middle final seek");
        Capture capture(6,4096,first);audition.setCaptureObserver(&capture);expect(audition.play(request,muted,false,true).isEmpty(),"long injected section audition");awaitPlaying(audition,*injected);juce::AudioBuffer<float> out(maxChannels,256);
        for(int n=0;n<16;++n){injected->tick(out);audition.poll();std::this_thread::sleep_for(std::chrono::milliseconds(1));}audition.stop();audition.setCaptureObserver(nullptr);
        juce::AudioBuffer<float> raw(6,4096);for(int ch=0;ch<6;++ch)for(int n=0;n<4096;++n)raw.setSample(ch,n,(ch==3&&first+n==impulseFrame?1.5f:float(double(longOperand(0,ch,first+n))/32768+double(longOperand(1,ch,first+n))/32768)));
        expect(capture.valid&&capture.next==first+4096&&sameSamples(capture.data,oracleAudio(raw,request.context.processing,first))&&audition.preparedCopy()==copy&&audition.underruns()==0,"long all-channel oracle capture exact and continuous across transitions; copy reused");
    }
    audition.disconnect();audition.sourceChanged({});stage("prepare_capture_ms",t);
    const auto requestExport=ex::capture(s);t=compositionNow();auto rendering=std::async(std::launch::async,[&]{return ex::render(requestExport,root,job);});auto rendered=compositionAwait(rendering);
    if(!rendered.value)std::cerr<<rendered.message<<'\n';expect(rendered.value.has_value(),"one real three-hour section-processed RF64 render/reopen/measure");stage("render_verify_ms",t);
    const auto artifact=*rendered.value;expect(artifact->spec().format.rf64&&artifact->file().getSize()==9331200104LL&&ex::artifactEvaluation(*artifact,requestExport).verdict==proposal::Verdict::pass,"real RF64 artifact and measured rules verified");
    const auto destination=root.getChildFile("export");t=compositionNow();auto publishing=std::async(std::launch::async,[&]{return ex::publish(*artifact,requestExport,destination,ex::approvalIdentity(*artifact,requestExport),job);});expect(compositionAwait(publishing).value.has_value(),"explicit verified local long export");stage("publish_ms",t);
    pulse.stopTimer();std::sort(pulse.gaps.begin(),pulse.gaps.end());expect(pulse.gaps.size()>10,"long UI heartbeat sampled");const auto p95=pulse.gaps[size_t((pulse.gaps.size()-1)*.95)];const auto growth=peak.load()-baseline;const auto reports=destination.getChildFile("report.json").getSize()+destination.getChildFile("report.html").getSize();
    expect(saveMs<=5000&&parseMs<=5000&&session.getSize()<=int64_t(sf::maximumBytes),"long save/parse/session budgets");
    expect(growth<=256ULL*1024*1024&&p95<=100&&reports<=16*1024*1024,"long additional RSS UI responsiveness and report budgets");
#if defined(__arm64__) || defined(__aarch64__)
    constexpr double timeBudget=240*60000.;
#else
    constexpr double timeBudget=300*60000.;
#endif
    expect(compositionNow()-begin<=timeBudget,"single consolidated long workflow within actual-platform time budget");
    longEvidence=juce::var(new juce::DynamicObject());auto* e=longEvidence.getDynamicObject();e->setProperty("private_passes",automatic->passes());e->setProperty("impulse_frame",juce::String(impulseFrame));e->setProperty("frames",juce::String(longFrames));e->setProperty("rate",48000);e->setProperty("channels",6);e->setProperty("intervals",1000);e->setProperty("capture_cases",8);e->setProperty("stages_ms",stages);e->setProperty("total_ms",compositionNow()-begin);e->setProperty("cancel_ms",cancelMs);e->setProperty("ui_p95_ms",p95);e->setProperty("additional_peak_rss_bytes",juce::int64(growth));e->setProperty("disk_growth_bytes",initialFree-minimumFree.load());e->setProperty("session_bytes",session.getSize());e->setProperty("report_bytes",reports);e->setProperty("artifact_sha256",artifact->hash());
}
}

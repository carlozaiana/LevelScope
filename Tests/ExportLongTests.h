#pragma once
namespace
{
uint64_t exportRss(){mach_task_basic_info_data_t info{};mach_msg_type_number_t count=MACH_TASK_BASIC_INFO_COUNT;return task_info(mach_task_self(),MACH_TASK_BASIC_INFO,reinterpret_cast<task_info_t>(&info),&count)==KERN_SUCCESS?info.resident_size:0;}
double exportNow(){return juce::Time::getMillisecondCounterHiRes();}
struct ExportPulse:juce::Timer
{
    double previous=exportNow();std::vector<double> gaps;
    void timerCallback() override{const double now=exportNow();gaps.push_back(now-previous);previous=now;}
};
template<class T> T exportAwait(std::future<T>& work)
{while(work.wait_for(std::chrono::milliseconds(0))!=std::future_status::ready)juce::MessageManager::getInstance()->runDispatchLoopUntil(10);return work.get();}
double auditionLargeArtifact(const juce::File& file)
{
    const auto began=exportNow();JobControl control;const auto media=validateMedia(file,control);expect(media.succeeded(),"published RF64 reimports with validated metadata");
    auto request=AuditionRequest::source(file,media.media,layoutFor(6));request.context.sessionId="m10-rf64-audition";request.context.buildVersion=LEVELSCOPE_MEASUREMENT_BUILD;request.currentInitialized=true;
    auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));OutputRouting muted;
    audition.synchronize(request,muted);expect(audition.prepare(request).isEmpty(),"RF64 verified copy prepares while device closed");
    while(audition.preparing()&&exportNow()-began<900000){audition.poll();juce::MessageManager::getInstance()->runDispatchLoopUntil(10);}
    expect(!audition.preparing()&&audition.preparedCopy()&&!injected->opened,"RF64 source preparation completes without physical device");const auto* copy=audition.preparedCopy();
    auto reader=ProgrammeReader::open({file},media.media);expect(reader!=nullptr&&audition.connect(48000).isEmpty(),"RF64 reader and injected device ready");
    for(double gain:{0.,-3.305074118944})for(int64_t first:{int64_t(0),int64_t(259200019),int64_t(518400000-4096)})
    {
        request.domain=gain==0?MeasurementDomain::source:MeasurementDomain::currentState;request.context.processing.staticGainDb=gain;++request.context.ticket.sessionRevision;
        audition.synchronize(request,muted);expect(audition.seek(first,request,muted).isEmpty(),"RF64 exact sample seek");
        Capture captured(6,4096,first);audition.setCaptureObserver(&captured);expect(audition.play(request,muted,false,true).isEmpty(),"RF64 all-muted injected audition");awaitPlaying(audition,*injected);
        juce::AudioBuffer<float> block(maxChannels,256);for(int n=0;n<16;++n){injected->tick(block);audition.poll();std::this_thread::sleep_for(std::chrono::milliseconds(1));}
        audition.stop();audition.setCaptureObserver(nullptr);juce::AudioBuffer<float> expected(6,4096);expect(reader->read(expected.getArrayOfWritePointers(),6,first,4096),"RF64 expected seek samples read");
        StaticGainProcessor processing(gain);expect(processing.process(expected.getArrayOfWritePointers(),6,4096),"RF64 shared full-precision oracle");
        bool exact=captured.valid&&captured.next==first+4096;for(int ch=0;ch<6;++ch)exact&=std::memcmp(expected.getReadPointer(ch),captured.data.getReadPointer(ch),4096*sizeof(float))==0;
        expect(exact&&audition.preparedCopy()==copy&&audition.underruns()==0,"RF64 first/middle/final Source/Current all-channel capture exact and copy reused");
    }
    audition.disconnect();return exportNow()-began;
}
void longExport(const juce::File& input,const juce::File& destination)
{
    namespace sf=levelscope::standalone::sessionfile;
    JobControl control;const auto begin=exportNow();auto reading=sf::read(input,control);
    expect(reading.document&&reading.document->source.assembly,"long export grouped session parsed");auto document=*reading.document;
    const auto members=document.source.assembly->members;for(const auto& m:members)sf::relinkAudioMember(document,m.id,input.getParentDirectory().getChildFile(m.file.getFileName()));
    auto verified=sf::verify(document,control);expect(verified.document&&verified.document->audioResolved,"three-hour group files verified");
    StandaloneSessionModel s;juce::String error;expect(sf::adopt(*verified.document,*verified.prepared,s,error),"three-hour programme adopted");
    expect(s.sourceDocument.media.frames==518400000&&s.sourceDocument.media.channels==6&&s.sourceDocument.media.sampleRate==48000,"real three-hour 5.1 48 kHz scope");
    s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
    std::cout<<"LONG source measurement started"<<std::endl;measure(s,MeasurementDomain::source);
    expect(s.setCurrentGainDb(-23-*s.source.integratedLufs),"long exact normalization gain");
    std::cout<<"LONG Current measurement started; elapsed_ms="<<exportNow()-begin<<std::endl;measure(s,MeasurementDomain::currentState);
    std::cout<<"LONG measurements complete; elapsed_ms="<<exportNow()-begin<<std::endl;
    const auto uiSession=destination.getParentDirectory().getChildFile("long-ui-ready.levelscope");
    expect(sf::save(sf::capture(s),uiSession,control).succeeded(),"save coherent long fixture for separate silent UI qualification");
    s.synchronizeSections();auto request=ex::capture(s);request.sections.sections.clear();request.sections.nextId=1001;
    for(int i=0;i<1000;++i)request.sections.sections.push_back({uint64_t(i+1),int64_t(i)*518400,int64_t(i+1)*518400,"Annotation "+juce::String(i+1)});
    expect(request.sections.valid()&&ex::readiness(request).isEmpty(),"1000 annotations and measured export readiness");
    const double prerequisiteMs=exportNow()-begin;verified.prepared.reset();
    const auto baseline=exportRss();std::atomic<uint64_t> peak{baseline};std::atomic<bool> stop{false};
    std::thread sampling([&]{while(!stop){peak=std::max(peak.load(),exportRss());std::this_thread::sleep_for(std::chrono::milliseconds(5));}});
    struct Join{std::atomic<bool>& stop;std::thread& thread;~Join(){stop=true;thread.join();}}join{stop,sampling};
    ExportPulse pulse;pulse.startTimer(10);const auto parent=destination.getParentDirectory();
    expect(!destination.exists()&&parent.isDirectory(),"fresh long export destination");
    JobControl cancellation;auto cancelling=std::async(std::launch::async,[&]{return ex::render(request,parent,cancellation);});
    while(cancellation.progress.load()==0&&cancelling.wait_for(std::chrono::milliseconds(1))!=std::future_status::ready)juce::MessageManager::getInstance()->runDispatchLoopUntil(10);
    double start=exportNow();cancellation.cancel=true;auto cancelled=exportAwait(cancelling);const auto cancelMs=exportNow()-start;
    expect(cancelled.error==ErrorCode::cancelled&&cancelMs<=1000,"active three-hour preparation cancellation <=1 second");
    juce::var stageTimes(new juce::DynamicObject());const auto renderBegin=exportNow();pcm::Hooks hooks;
    hooks.reached=[&](const char* name){if(juce::String(name)=="audio-create"||juce::String(name)=="verify"||juce::String(name)=="verified"){stageTimes.getDynamicObject()->setProperty(name,exportNow()-renderBegin);std::cout<<"LONG "<<name<<" elapsed_ms="<<exportNow()-renderBegin<<std::endl;}};
    auto rendering=std::async(std::launch::async,[&]{return ex::render(request,parent,control,hooks);});auto result=exportAwait(rendering);
    if(!result.value)std::cerr<<result.message<<'\n';expect(result.value.has_value(),"actual 9.33 GB RF64 rendered and fully remeasured");auto candidate=*result.value;
    const auto renderMs=exportNow()-renderBegin;expect(candidate->spec().format.rf64&&candidate->file().getSize()==9331200104LL,"real RF64 complete physical size");
    expect(ex::artifactEvaluation(*candidate,request).verdict==proposal::Verdict::pass,"real long actual PCM passes selected rules");
    // Exercise the exact lightweight paths used on every export UI refresh.
    (void)s.evaluation(MeasurementDomain::currentState); // pre-existing session cache
    std::vector<double> inspection;
    for(int n=0;n<100;++n)
    {
        const auto tick=exportNow();const auto captured=ex::capture(s);
        expect(ex::preview(captured,false).isNotEmpty()&&candidate->matchesAudio(captured)
            &&ex::artifactEvaluation(*candidate,captured).verdict==proposal::Verdict::pass
            &&s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"long export UI uses accepted immutable caches");
        inspection.push_back(exportNow()-tick);
    }
    std::sort(inspection.begin(),inspection.end());const auto inspectionP95=inspection[94];
    expect(inspectionP95<=100,"long export UI inspection p95 <=100ms");
    const auto owned=candidate->folder().findChildFiles(juce::File::findFiles,false);expect(owned.size()==1&&owned[0].getSize()==9331200104LL,"single output payload before report; no second interleaved render");
    start=exportNow();auto publication=std::async(std::launch::async,[&]{return ex::publish(*candidate,request,destination,ex::approvalIdentity(*candidate,request),control);});const auto published=exportAwait(publication);const auto publishMs=exportNow()-start;
    expect(published.value.has_value(),"real large export atomically published with reports");pulse.stopTimer();
    std::sort(pulse.gaps.begin(),pulse.gaps.end());expect(pulse.gaps.size()>=10,"UI heartbeat during long offline work");const auto p95=pulse.gaps[static_cast<size_t>((pulse.gaps.size()-1)*.95)];
    const auto growth=peak.load()-baseline;const auto reportBytes=destination.getChildFile("report.json").getSize()+destination.getChildFile("report.html").getSize();
    expect(growth<=256ULL*1024*1024&&p95<=100&&reportBytes<=16*1024*1024,"three-hour memory/UI/report-size budgets");
    expect(prerequisiteMs+renderMs+publishMs<=10800000,"full three-hour preparation/render/verification within programme duration");
    const auto artifactHash=candidate->hash();const auto artifactBytes=candidate->file().getSize();const auto seed=candidate->spec().seedHex;const auto copies=candidate->verifiedSource()->files();
    start=exportNow();result.value.reset();candidate.reset();const auto cleanupMs=exportNow()-start;
    bool cleaned=true;for(const auto& file:copies)cleaned&=!file.exists();expect(cleaned,"successful long render releases only its owned private source copies");
    const auto reimportAuditionMs=auditionLargeArtifact(destination.getChildFile("programme.wav"));
    expect(prerequisiteMs+renderMs+publishMs+reimportAuditionMs<=10800000,"complete long export and additional RF64 audition checks within three hours");
    auto* evidence=new juce::DynamicObject();evidence->setProperty("rf64_audition_cases",6);evidence->setProperty("rf64_reimport_audition_ms",reimportAuditionMs);evidence->setProperty("prerequisite_ms",prerequisiteMs);evidence->setProperty("render_verify_ms",renderMs);evidence->setProperty("publish_ms",publishMs);evidence->setProperty("stage_ms",stageTimes);
    evidence->setProperty("cancel_ms",cancelMs);evidence->setProperty("ui_p95_ms",p95);evidence->setProperty("inspection_p95_ms",inspectionP95);evidence->setProperty("ui_callbacks",int(pulse.gaps.size()));evidence->setProperty("rss_growth_bytes",juce::int64(growth));evidence->setProperty("report_bytes",reportBytes);evidence->setProperty("artifact_bytes",artifactBytes);evidence->setProperty("artifact_sha256",artifactHash);evidence->setProperty("cleanup_ms",cleanupMs);evidence->setProperty("source_bytes",request.current.media.bytes);evidence->setProperty("exact_gain",request.current.processing.staticGainDb);evidence->setProperty("seed_hex",seed);evidence->setProperty("frames",juce::String(request.current.media.frames));evidence->setProperty("annotations",1000);longEvidence=juce::var(evidence);
}
}

#pragma once
// Included by ProgrammeAssemblyTests: shares its fixture, assertion and capture helpers.
namespace
{
int referenceCases=0;
juce::var longProgrammeEvidence;
uint64_t programmeRss(){mach_task_basic_info_data_t info{};mach_msg_type_number_t count=MACH_TASK_BASIC_INFO_COUNT;return task_info(mach_task_self(),MACH_TASK_BASIC_INFO,reinterpret_cast<task_info_t>(&info),&count)==KERN_SUCCESS?info.resident_size:0;}
double programmeNow(){return juce::Time::getMillisecondCounterHiRes();}
void nativeReferences(const juce::File& root,const juce::File& videoFolder,const juce::File& proresFolder)
{
    auto f=makeFixture(root,6,48000,192001,true);auto s=model(f.group);measure(s,MeasurementDomain::source);s.setCurrentGainDb(-3.305074118944);measure(s,MeasurementDomain::currentState);
    std::vector<juce::File> files{videoFolder.getChildFile("cfr-exact.mov"),videoFolder.getChildFile("vfr-exact.mov")};for(const char* name:{"proxy-with-tone.mov","lt-with-tone.mov","422-with-tone.mov"})files.push_back(proresFolder.getChildFile(name));
    beginVideoReaderAudit();
    for(auto domain:{MeasurementDomain::source,MeasurementDomain::currentState})for(const auto& file:files)
    {
        s.selectAuditionDomain(domain);const auto request=s.auditionRequest();const auto audioIdentity=request.identity();auto decoder=video::makeDecoder();const auto asset=decoder->load(file,0,[]{return false;});expect(asset!=nullptr,"M9 native reference opens");
        s.referenceVideo={file,asset->sha256,asset->selectedStream,video::Time(-1001,30000)};
        video::Player player;player.attach(file,asset->selectedStream);uint64_t generation=1;
        for(auto offset:{video::Time(-1001,30000),video::Time(7,3)})for(auto time:{video::Time(0,1),video::Time(3,2),video::Time(3,1),video::Time(0,1),video::Time(f.frames,48000)})
        {
            s.referenceVideo.startOnAudio=offset;const auto target=time-offset;const auto expected=video::frameAt(asset->frames,target);player.invalidate();++generation;bool matched=false;const auto deadline=programmeNow()+10000;
            while(programmeNow()<deadline){auto view=player.update(target,generation);if(view.asset&&((expected&&view.picture&&view.picture->index==*expected)||(!expected&&!view.picture))){matched=true;break;}std::this_thread::sleep_for(std::chrono::milliseconds(2));}
            expect(matched,"Source/Current group rational signed-offset seek/restart/EOF picture selection");
            expect(s.auditionRequest().identity()==audioIdentity,"reference attachment/offset never changes programme audition identity");
        }
        player.clear();++referenceCases;
    }
    expect(referenceCases==10&&videoReaderAuditCount(0)>0&&videoReaderAuditCount(1)==0,"ten native group reference cases; embedded audio remains undecoded");
}
void exchangeOut(const juce::File& root)
{
    expect(!root.exists()&&root.createDirectory().wasOk(),"new M9 exchange destination");
    for(int channels:{2,6})
    {
        auto dir=root.getChildFile(juce::String(channels));auto f=makeFixture(dir,channels,48000,192001,true);auto s=model(f.group);s.setCurrentGainDb(-3.305074118944);measure(s,MeasurementDomain::source);measure(s,MeasurementDomain::currentState);
        s.synchronizeSections();juce::String error;expect(s.sectionDocument.split(1,77777,error),"exchange sections");s.documentChanged();s.referenceVideo.startOnAudio=video::Time(-1001,30000);JobControl c;
        expect(sf::save(sf::capture(s),dir.getChildFile("programme.levelscope"),c).succeeded(),"exchange complete group saved");
    }
}
void exchangeIn(const juce::File& root)
{
    for(int channels:{2,6})
    {
        const auto dir=root.getChildFile(juce::String(channels));JobControl c;auto saved=sf::read(dir.getChildFile("programme.levelscope"),c);expect(saved.document&&saved.document->source.assembly,"foreign group parsed");const auto confirmation=saved.document->source.assembly->alignmentConfirmation;
        const auto original=saved.document->sourceCache;expect(original!=nullptr,"foreign original provenance cache available for compatibility check");
        const auto members=saved.document->source.assembly->members;for(const auto& m:members)expect(sf::relinkAudioMember(*saved.document,m.id,dir.getChildFile(m.file.getFileName())),"explicit exchange relink by file slot");
        auto verified=sf::verify(*saved.document,c);expect(verified.document&&verified.document->audioResolved,"foreign group every member verifies");StandaloneSessionModel s;juce::String notice;expect(sf::adopt(*verified.document,*verified.prepared,s,notice),"foreign grouped session adopted");
        const bool compatible=original->context.buildVersion==s.measurementContext(MeasurementDomain::source).buildVersion&&original->context.platformAndFloatingPointMode==s.measurementContext(MeasurementDomain::source).platformAndFloatingPointMode;
        expect(s.source.hasMeasurement==compatible&&s.currentState.hasMeasurement==compatible,"conservative exact build/platform cache policy");
        expect(s.sourceDocument.assembly->alignmentConfirmation==confirmation&&s.currentStateDocument.staticGainDb==-3.305074118944&&s.referenceVideo.startOnAudio==video::Time(-1001,30000)&&s.sectionDocument.map().sections.size()==2,"cross-platform settings and exact alignment confirmation retained");
        measure(s,MeasurementDomain::source);measure(s,MeasurementDomain::currentState);expect(s.source.hasMeasurement&&s.currentState.hasMeasurement,"reopened foreign group can be remeasured");
    }
}
juce::Array<juce::var> migrationEvidence;
void migrateM8(const juce::File& file,const juce::File& out)
{
    JobControl c;auto d=sf::read(file,c);expect(d.document&&d.document->loadedSchema==1&&!d.document->source.assembly,("M8 core migration remains interleaved: "+d.message).toRawUTF8());
    juce::FileInputStream before(file);const auto hash=juce::SHA256(before).toHexString();auto r=sf::verify(*d.document,c);expect(r.document.has_value(),"M8 settings recover even if media missing");
    StandaloneSessionModel s;juce::String notice;expect(sf::adopt(*r.document,*r.prepared,s,notice),"M8 migrated session adoption");expect(s.openedSessionSchema==1&&!s.source.hasMeasurement&&!s.currentState.hasMeasurement,"M8 caches not reused across changed build; upgrade remains explicit");
    expect(s.currentStateDocument.staticGainDb==d.document->current.staticGainDb&&s.referenceVideo.startOnAudio==d.document->reference.startOnAudio&&s.sectionDocument.map().digest()==d.document->sections.digest(),"M8 migration exact settings sections rational alignment");
    const auto profile=s.getSelectedTargetProfile();expect(profile.id==d.document->profile.id&&profile.definitionVersion==d.document->profile.version&&profile.digest()==d.document->profile.digest&&s.programmeApplicability==d.document->applicability,"M8 migration exact profile definition and applicability");
    expect(s.sourceDocument.layout.roles==d.document->source.layout.roles&&s.sourceDocument.layout.count==d.document->source.layout.count&&s.sourceDocument.layout.confirmed==d.document->source.layout.confirmed&&s.sourceDocument.media.sha256==d.document->source.media.sha256,"M8 migration preserves ordered channel interpretation and original content identity");
    expect(sf::save(sf::capture(s),out,c).succeeded(),"explicit M9 save as");auto newer=sf::read(out,c);expect(newer.document&&newer.document->loadedSchema==2,"new core explicitly requires M9 reader");juce::FileInputStream after(file);expect(juce::SHA256(after).toHexString()==hash,"original M8 file unchanged by migration and Save As");
    auto row=juce::var(new juce::DynamicObject());row.getDynamicObject()->setProperty("original_sha256",hash);row.getDynamicObject()->setProperty("channels",s.sourceDocument.media.channels);row.getDynamicObject()->setProperty("audio_resolved",r.document->audioResolved);row.getDynamicObject()->setProperty("original_unchanged",true);migrationEvidence.add(row);
}
juce::Array<juce::var> realProgrammeEvidence;
void confirmedRealProgrammes(const juce::File& input,const juce::File& outputFolder)
{
    const auto request=juce::JSON::parse(input);expect(static_cast<bool>(request["alignment_confirmed"]),"explicit real-media alignment authorization required");
    const auto* groups=request["groups"].getArray();expect(groups!=nullptr,"real-media group definitions");
    for(const auto& definition:*groups)
    {
        const auto* rows=definition["files"].getArray();expect(rows!=nullptr,"real-media explicit file mapping");std::vector<juce::File> files;for(const auto& row:*rows)files.emplace_back(row["path"].toString());JobControl c;
        auto inspected=inspectProgrammeFiles(files,c);expect(inspected.error==ErrorCode::none,"real-media all members fully validate");auto group=inspected.group;
        const juce::StringArray roleNames{ "L","R","C","LFE","Ls","Rs" };const auto roles=assemblyRoles(static_cast<int>(files.size()));
        for(size_t i=0;i<files.size();++i){const auto& row=(*rows)[static_cast<int>(i)];expect(group.members[i].media.sha256==row["sha256"].toString(),"real member still matches user-confirmed bytes");const int role=roleNames.indexOf(row["role"].toString());expect(role>=0&&role<static_cast<int>(roles.size()),"explicit supported real channel role");group.members[i].role=roles[static_cast<size_t>(role)];}
        expect(group.orderByRoles()&&group.confirmAlignment(),"user-confirmed exact real assembly");auto s=model(group);s.setCurrentGainDb(-3.305074118944);measure(s,MeasurementDomain::source);measure(s,MeasurementDomain::currentState);
        auto record=juce::var(new juce::DynamicObject());auto* e=record.getDynamicObject();e->setProperty("name",definition["name"]);e->setProperty("assembly_identity",group.definitionIdentity());e->setProperty("frames",s.sourceDocument.media.frames);e->setProperty("rate",s.sourceDocument.media.sampleRate);e->setProperty("channels",s.sourceDocument.media.channels);e->setProperty("processing_gain_db",s.currentStateDocument.staticGainDb);
        for(auto domain:{MeasurementDomain::source,MeasurementDomain::currentState}){const auto& a=domain==MeasurementDomain::source?s.source.analysis:s.currentState.analysis;auto metrics=juce::var(new juce::DynamicObject());metrics.getDynamicObject()->setProperty("integrated_lufs",a->measured.integratedLufs?juce::var(*a->measured.integratedLufs):juce::var());metrics.getDynamicObject()->setProperty("true_peak_dbtp",a->measured.truePeakDbtp?juce::var(*a->measured.truePeakDbtp):juce::var());metrics.getDynamicObject()->setProperty("identity",a->identity);e->setProperty(domain==MeasurementDomain::source?"source":"current",metrics);}
        auto document=outputFolder.getChildFile("real-"+definition["name"].toString()+".levelscope");expect(sf::save(sf::capture(s),document,c).succeeded(),"real programme saved to separate test document");auto read=sf::read(document,c);expect(read.document.has_value(),"real grouped document read");auto verified=sf::verify(*read.document,c);expect(verified.document&&verified.document->audioResolved,"real grouped session members verify unchanged");
        // Re-hash all originals after measurement/reopening; temporary snapshots never replace them.
        for(size_t i=0;i<files.size();++i){juce::FileInputStream stream(files[i]);expect(juce::SHA256(stream).toHexString()==(*rows)[static_cast<int>(i)]["sha256"].toString(),"original real media unchanged");}realProgrammeEvidence.add(record);
    }
}
void prepareProgramme(const juce::File& input,const juce::File& destination)
{
    JobControl c;auto read=sf::read(input,c);expect(read.document&&read.document->source.assembly,"prepare grouped qualification session");auto d=*read.document;const auto members=d.source.assembly->members;
    for(const auto& m:members)expect(sf::relinkAudioMember(d,m.id,input.getParentDirectory().getChildFile(m.file.getFileName())),"qualification explicit stable-slot relocation");
    auto verified=sf::verify(d,c);expect(verified.document&&verified.document->audioResolved,"qualification relocated group content verified");StandaloneSessionModel s;juce::String error;expect(sf::adopt(*verified.document,*verified.prepared,s,error),"qualification relocated session adoption");expect(sf::save(sf::capture(s),destination,c).succeeded(),"qualification relocated session saved separately");
}
void generateLong(const juce::File& root,int rate,int seconds)
{
    expect(!root.exists()&&root.createDirectory().wasOk(),"new long group destination");std::vector<juce::File> files;const int64_t frames=int64_t(rate)*seconds;
    for(int ch=0;ch<6;++ch)
    {
        auto file=root.getChildFile("mono-"+juce::String(ch)+".flac");auto out=file.createOutputStream();juce::FlacAudioFormat format;std::unique_ptr<juce::AudioFormatWriter> writer(format.createWriterFor(out.get(),rate,1,16,{},0));expect(writer!=nullptr,"long FLAC writer");out.release();juce::AudioBuffer<float> block(1,65536);
        for(int i=0;i<65536;++i)block.setSample(0,i,static_cast<float>(std::llround(500*(ch+1)*std::sin((i%128)*2*3.14159265358979323846/128)))/32768.f);
        for(int64_t first=0;first<frames;first+=65536){const int n=static_cast<int>(std::min<int64_t>(65536,frames-first));if(first+n==frames)block.setSample(0,n-1,(10000.f+ch*1000)/32768);if(!writer->writeFromAudioSampleBuffer(block,0,n))throw std::runtime_error("long fixture write");}
        writer.reset();files.push_back(file);
    }
    JobControl c;auto inspected=inspectProgrammeFiles(files,c);expect(inspected.error==ErrorCode::none,"real long group inspect");auto group=inspected.group;auto roles=assemblyRoles(6);for(int ch=0;ch<6;++ch)group.members[ch].role=roles[ch];expect(group.confirmAlignment(),"test-generator known alignment confirmed");auto s=model(group);s.setCurrentGainDb(-3.305074118944);expect(sf::save(sf::capture(s),root.getChildFile("programme.levelscope"),c).succeeded(),"long programme session saved");
}
void longProgramme(const juce::File& input,const juce::File& output)
{
    auto report=juce::var(new juce::DynamicObject());auto* evidence=report.getDynamicObject();JobControl c;auto read=sf::read(input,c);expect(read.document&&read.document->source.assembly,"real long grouped session read");
    auto d=*read.document;const auto members=d.source.assembly->members;for(const auto& m:members)sf::relinkAudioMember(d,m.id,input.getParentDirectory().getChildFile(m.file.getFileName()));
    const double begin=programmeNow();auto valid=sf::verify(d,c);expect(valid.document&&valid.document->audioResolved,"real three-hour member verification");evidence->setProperty("verify_ms",programmeNow()-begin);StandaloneSessionModel s;juce::String error;expect(sf::adopt(*valid.document,*valid.prepared,s,error),"long group install");
    expect(s.sourceDocument.media.frames==518400000&&s.sourceDocument.media.sampleRate==48000&&s.sourceDocument.media.channels==6,"real three-hour 48 kHz six-file scope");
    const auto baseline=programmeRss();std::atomic<uint64_t> peak{baseline};std::atomic<bool> stop{false};std::thread sampler([&]{while(!stop){peak=std::max(peak.load(),programmeRss());std::this_thread::sleep_for(std::chrono::milliseconds(5));}});
    struct Join{std::atomic<bool>& stop;std::thread& thread;~Join(){stop=true;thread.join();}} join{stop,sampler};
    for(auto domain:{MeasurementDomain::source,MeasurementDomain::currentState}){std::cout<<"Three-hour "<<(domain==MeasurementDomain::source?"Source":"Current")<<" measurement started"<<std::endl;auto start=programmeNow();measure(s,domain);const double elapsed=programmeNow()-start;evidence->setProperty(domain==MeasurementDomain::source?"source_measure_ms":"current_measure_ms",elapsed);std::cout<<"Measurement completed in "<<elapsed/1000<<" seconds"<<std::endl;}
    s.synchronizeSections();for(int n=1;n<1000;++n)expect(s.sectionDocument.split(n,s.sourceDocument.media.frames*n/1000,error),"long group section split");s.documentChanged();auto captured=sf::capture(s);
    struct Pulse final:juce::Timer{std::vector<double> gaps;double last=programmeNow();void timerCallback()override{auto now=programmeNow();gaps.push_back(now-last);last=now;}}pulse;
    pulse.startTimer(5);const auto await=[](auto& pending){while(pending.wait_for(std::chrono::milliseconds(0))!=std::future_status::ready)juce::MessageManager::getInstance()->runDispatchLoopUntil(5);return pending.get();};
    auto start=programmeNow();auto writing=std::async(std::launch::async,[&]{return sf::save(captured,output,c);});expect(await(writing).succeeded(),"real long group coherent save");const double saveMs=programmeNow()-start;start=programmeNow();auto reading=std::async(std::launch::async,[&]{return sf::read(output,c);});auto reopened=await(reading);const double parseMs=programmeNow()-start;
    pulse.stopTimer();std::sort(pulse.gaps.begin(),pulse.gaps.end());expect(pulse.gaps.size()>=5,"long group UI heartbeat observed during persistence");const double p95=pulse.gaps[static_cast<size_t>((pulse.gaps.size()-1)*.95)];expect(p95<=100,"long group persistence UI p95 under 100 ms");expect(reopened.document&&reopened.document->sourceCache&&reopened.document->currentCache,"real full-resolution long caches roundtrip");
    auto complete=sf::verify(*reopened.document,c);expect(complete.document&&complete.document->audioResolved,"real long saved group fully reopened and verified");StandaloneSessionModel restored;expect(sf::adopt(*complete.document,*complete.prepared,restored,error)&&restored.source.hasMeasurement&&restored.currentState.hasMeasurement,"real long measured caches restored after full member verification");
    expect(saveMs<=5000&&parseMs<=5000&&output.getSize()<=int64_t(sf::maximumBytes),"long group persistence time and size budgets");
    JobControl cancel;auto opening=std::async(std::launch::async,[&]{return sf::verify(d,cancel);});while(cancel.progress.load()==0&&opening.wait_for(std::chrono::milliseconds(1))!=std::future_status::ready){}start=programmeNow();cancel.cancel=true;const auto interrupted=opening.get();const double cancelMs=programmeNow()-start;expect(interrupted.error==ErrorCode::cancelled&&cancelMs<=1000,"active long group open cancellation within one second");
    auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));s.selectAuditionDomain(MeasurementDomain::currentState);auto request=s.auditionRequest();OutputRouting routing;audition.synchronize(request,routing);expect(audition.prepare(request).isEmpty(),"long group silent preparation");
    const auto preparationDeadline=programmeNow()+1800000;while(audition.preparing()&&programmeNow()<preparationDeadline){audition.poll();std::this_thread::sleep_for(std::chrono::milliseconds(2));}expect(audition.preparedCopy()!=nullptr&&!injected->opened,"long complete verified group prepared with device closed");
    JobControl prepareCancel;auto context=s.measurementContext(MeasurementDomain::currentState);auto copying=std::async(std::launch::async,[&]{return prepareVerifiedSource(context,prepareCancel);});
    while(prepareCancel.progress.load()==0&&copying.wait_for(std::chrono::milliseconds(1))!=std::future_status::ready){}start=programmeNow();prepareCancel.cancel=true;auto cancelledCopy=copying.get();const double preparationCancelMs=programmeNow()-start;expect(cancelledCopy.error==ErrorCode::cancelled&&preparationCancelMs<=1000,"active long group preparation cancellation within one second");
    const auto* copy=audition.preparedCopy();int64_t payload=0;for(auto file:copy->files())payload+=file.getSize();expect(payload==s.sourceDocument.media.bytes,"long group has no extra interleaved temporary payload");
    expect(audition.connect(48000).isEmpty(),"long injected device only");double maxSeek=0;
    for(int64_t frame:{int64_t(0),int64_t(123456789),s.sourceDocument.media.frames-513})
    {
        Capture played(6,static_cast<int>(std::min<int64_t>(1024,s.sourceDocument.media.frames-frame)),frame);audition.setCaptureObserver(&played);expect(audition.seek(frame,request,routing).isEmpty(),"long group seek request");start=programmeNow();expect(audition.play(request,routing,false,true).isEmpty(),"long muted injected play");awaitPlaying(audition,*injected);maxSeek=std::max(maxSeek,programmeNow()-start);
        juce::AudioBuffer<float> block(maxChannels,played.data.getNumSamples());injected->tick(block);audition.poll();audition.stop();
        juce::AudioBuffer<float> oracle(6,played.data.getNumSamples());for(int ch=0;ch<6;++ch)for(int i=0;i<oracle.getNumSamples();++i){const auto pos=frame+i;oracle.setSample(ch,i,pos==s.sourceDocument.media.frames-1?(10000.f+ch*1000)/32768:static_cast<float>(std::llround(500*(ch+1)*std::sin((pos%128)*2*3.14159265358979323846/128)))/32768.f);}StaticGainProcessor processing(s.currentStateDocument.staticGainDb);expect(processing.process(oracle.getArrayOfWritePointers(),6,oracle.getNumSamples()),"long exact processing oracle");bool equal=played.valid&&played.next==frame+oracle.getNumSamples();for(int ch=0;ch<6;++ch)equal&=std::memcmp(oracle.getReadPointer(ch),played.data.getReadPointer(ch),oracle.getNumSamples()*4)==0;expect(equal&&audition.preparedCopy()==copy,"long random/final-sample all-channel capture parity and copy reuse");
    }
    audition.setCaptureObserver(nullptr);audition.disconnect();expect(maxSeek<=2000,"prepared three-hour seek/prefill within two seconds");const uint64_t growth=peak.load()-baseline;expect(growth<=256ULL*1024*1024,"long group bounded reader/measurement memory including results under 256 MiB growth");
    evidence->setProperty("frames",s.sourceDocument.media.frames);evidence->setProperty("channels",6);evidence->setProperty("rate",48000);evidence->setProperty("programme_identity",s.sourceDocument.media.sha256);evidence->setProperty("save_ms",saveMs);evidence->setProperty("parse_ms",parseMs);evidence->setProperty("cancel_ms",cancelMs);evidence->setProperty("max_seek_prefill_ms",maxSeek);evidence->setProperty("rss_growth_bytes",static_cast<juce::int64>(growth));evidence->setProperty("session_bytes",output.getSize());evidence->setProperty("owned_payload_bytes",payload);evidence->setProperty("original_payload_bytes",s.sourceDocument.media.bytes);evidence->setProperty("preparation_cancel_ms",preparationCancelMs);evidence->setProperty("ui_callback_p95_ms",p95);evidence->setProperty("ui_callback_count",static_cast<int>(pulse.gaps.size()));longProgrammeEvidence=report;
}
}

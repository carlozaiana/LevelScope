#include "Standalone/SessionDocument.h"
#include "Standalone/StandaloneMeasurementLayout.h"
#include <iostream>
#include <juce_cryptography/juce_cryptography.h>
#include <fstream>
#include <mach/mach.h>
#include <future>
#include <sys/resource.h>
#include <thread>
using namespace levelscope; using namespace levelscope::standalone;
namespace sf = levelscope::standalone::sessionfile;
namespace
{
juce::var performance(new juce::DynamicObject());
juce::var mediaPerformance;
int checks=0;void expect(bool ok,const juce::String& text){++checks;if(!ok)throw std::runtime_error(text.toStdString());}

std::string bytes(const juce::File& f)
{juce::MemoryBlock b;expect(f.loadFileAsData(b),"read file bytes");return {static_cast<const char*>(b.getData()),b.getSize()};}
uint64_t little(const std::string& b,size_t pos){uint64_t n=0;for(int i=0;i<8;++i)n|=uint64_t(static_cast<uint8_t>(b[pos+i]))<<(8*i);return n;}
void coreEdit(const juce::File& from,const juce::File& to,const std::function<void(juce::var&)>& edit)
{
    auto b=bytes(from);auto length=little(b,8);auto v=juce::JSON::parse(juce::String::fromUTF8(b.data()+80,static_cast<int>(length)));edit(v);
    auto json=juce::JSON::toString(v,true).toStdString();auto hash=juce::SHA256(json.data(),json.size()).toHexString().toStdString();std::string out="LSSESS01";
    for(int i=0;i<8;++i)out.push_back(static_cast<char>(uint64_t(json.size())>>(8*i)));out+=hash;out+=json;out+=b.substr(80+length);expect(to.replaceWithData(out.data(),out.size()),"write edited document fixture");
}

void measure(StandaloneSessionModel& s,MeasurementDomain d)
{
    expect(s.job.begin(JobStage::analysis,s.revision),"begin measurement");auto c=s.measurementContext(d);JobControl control;
    auto r=OfflineAnalysisService().analyze(c,d,control);expect(r.value.has_value(),r.message);expect(s.job.finish(c.ticket,s.revision,r.error,r.message),"finish measurement");expect(s.acceptMeasurement(std::move(*r.value)),"accept measured snapshot");
}
juce::File fixture(const juce::File& folder,int channels,int rate,int frames=0,bool silent=false)
{
    auto f=folder.getChildFile("tone-"+juce::String(channels)+"-"+juce::String(rate)+".wav");auto out=f.createOutputStream();expect(out!=nullptr,"fixture output: "+f.getFullPathName());juce::WavAudioFormat format;
    std::unique_ptr<juce::AudioFormatWriter> writer(format.createWriterFor(out.get(),rate,channels,24,{},0));expect(writer!=nullptr,"fixture writer");out.release();juce::AudioBuffer<float> b(channels,frames>0?frames:rate*4);
    for(int c=0;c<channels;++c)for(int i=0;i<b.getNumSamples();++i)b.setSample(c,i,silent?0.f:static_cast<float>(.04*std::sin(2*3.141592653589793*(351+c*100)*i/rate)));
    expect(writer->writeFromAudioSampleBuffer(b,0,b.getNumSamples()),"fixture samples");return f;
}
StandaloneSessionModel model(const juce::File& file,int channels)
{
    StandaloneSessionModel s;JobControl c;auto m=validateMedia(file,c);expect(m.succeeded(),m.message);expect(s.setValidatedSource(file,m.media),"import source");
    if(channels==8){auto l=systemILayout();for(int n=0;n<8;++n)s.setChannelRole(n,l.roles[n]);}
    expect(s.confirmLayout(),"confirm roles");s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
    expect(s.initializeCurrentStateFromSource(),"initialize Current");expect(s.setCurrentGainDb(-6.020599913279624),"precise trim");measure(s,MeasurementDomain::source);measure(s,MeasurementDomain::currentState);s.synchronizeSections();juce::String error;expect(s.sectionDocument.split(1,s.sourceDocument.media.frames/2,error),error);s.documentChanged();s.referenceVideo.startOnAudio=video::Time(-1001,30000);s.documentChanged();return s;
}
void roundtrip(const juce::File& root,int channels,int rate)
{
    auto s=model(fixture(root,channels,rate),channels);expect(s.acceptTrimProposal(proposal::draftGlobalTrim(s.currentTrimRequest())),"live proposal prepared before saving");auto d=sf::capture(s);auto f=root.getChildFile("roundtrip.levelscope");JobControl control;
    auto written=sf::save(d,f,control);expect(written.succeeded(),written.message);expect(sf::saved(s,d,f)&&!s.documentDirty(),"save revision clean");
    auto read=sf::read(f,control);expect(read.document.has_value(),read.message);expect(read.document->sourceCache&&read.document->currentCache,read.document->notices.joinIntoString(";"));
    auto verified=sf::verify(*read.document,control);expect(verified.document&&verified.document->audioResolved,verified.message);
    StandaloneSessionModel opened;juce::String notice;expect(sf::install(*verified.document,opened,notice),notice);
    expect(opened.sessionId!=s.sessionId&&opened.documentId==s.documentId,"persistent/runtime identity separate");
    StandaloneSessionModel adopted;expect(verified.prepared&&sf::adopt(*verified.document,*verified.prepared,adopted,notice),"worker-prepared session adopted");
    expect(!sf::adopt(*verified.document,*verified.prepared,adopted,notice),"prepared opening cannot be reused");
    const auto stableEvaluation=adopted.evaluation(MeasurementDomain::currentState);auto stableAnalysis=adopted.currentState.analysis;
    adopted.setOutputRoute(0,0);adopted.referenceVideo.startOnAudio=video::Time(7,3);
    expect(adopted.currentState.analysis==stableAnalysis&&adopted.evaluation(MeasurementDomain::currentState).snapshotIdentity==stableEvaluation.snapshotIdentity,"monitor/video preserve reconstructed evaluation");
    adopted.setProgrammeApplicability(proposal::Applicability::unconfirmed);
    expect(adopted.currentState.analysis==stableAnalysis&&adopted.evaluation(MeasurementDomain::currentState).verdict!=proposal::Verdict::pass,"applicability reevaluation reuses measurements");
    adopted.setCurrentGainDb(-11.1);expect(adopted.evaluation(MeasurementDomain::currentState).verdict!=proposal::Verdict::pass&&!sf::capture(adopted).currentCache,"changed gain invalidates restored evaluation/cache");expect(opened.currentState.hasMeasurement&&opened.source.hasMeasurement,notice);expect(!opened.trimDraft,"saved proposals/predictions never restore as authority");
    expect(opened.currentStateDocument.staticGainDb==s.currentStateDocument.staticGainDb,"full precision gain");expect(opened.sectionDocument.map().digest()==s.sectionDocument.map().digest(),"exact sections");expect(opened.referenceVideo.startOnAudio==s.referenceVideo.startOnAudio,"rational offset");expect(opened.evaluation(MeasurementDomain::currentState).verdict==s.evaluation(MeasurementDomain::currentState).verdict,"evaluation reconstructed");
    expect(opened.currentState.analysis->originSessionId==s.sessionId,"original measurement provenance");for(auto route:opened.sourceDocument.routing.fileToDevice)expect(route==-1,"restored routes muted");
    auto stale=*s.currentState.analysis;expect(!opened.acceptMeasurement(stale),"old runtime completion rejected");
    auto snap=sf::capture(s);s.setCurrentGainDb(-2.345678901234567);expect(sf::saved(s,snap,f)&&s.documentDirty(),"edits during save stay dirty");
    auto unverified=sf::capture(s);expect(!unverified.currentCache&&unverified.current.staticGainDb==-2.345678901234567,"unverified gain does not save old cache");
    auto previous=bytes(f);for(const char* stage:{"space","create","write","flush","replace"})
    {JobControl c;sf::FileHooks h;h.fail=[&](const char* x){return juce::String(x)==stage;};auto r=sf::save(snap,f,c,h);expect(!r.succeeded(),"injected save failure");expect(bytes(f)==previous,"previous saved bytes retained");expect(root.findChildFiles(juce::File::findFiles,false,"*.partial").isEmpty(),"owned partial cleaned");}
    sf::FileHooks collision;juce::File occupied;collision.beforeCreate=[&](const juce::File& path){occupied=path;expect(path.replaceWithText("unowned sentinel"),"simulate exclusive-create collision");};
    JobControl collisionControl;expect(!sf::save(snap,f,collisionControl,collision).succeeded(),"exclusive-create collision rejected");expect(occupied.loadFileAsString()=="unowned sentinel"&&bytes(f)==previous,"failed create preserves unowned file and previous session");occupied.deleteFile();
    JobControl cancelled;sf::FileHooks hooks;hooks.reached=[&](const char* x){if(juce::String(x)=="write")cancelled.cancel=true;};auto cancelledSave=sf::save(snap,f,cancelled,hooks);expect(cancelledSave.error==ErrorCode::cancelled,"cancelled save");expect(bytes(f)==previous,"cancelled save preserves file");
    juce::MemoryBlock damagedBytes;expect(f.loadFileAsData(damagedBytes),"read saved bytes");auto damaged=root.getChildFile("damaged.levelscope");static_cast<uint8_t*>(damagedBytes.getData())[damagedBytes.getSize()-1]^=1;expect(damaged.replaceWithData(damagedBytes.getData(),damagedBytes.getSize()),"write damaged cache");JobControl c;auto recovered=sf::read(damaged,c);expect(recovered.document.has_value(),"damaged cache preserves valid core");expect(!recovered.document->currentCache&&recovered.document->sourceCache,"cache damage isolated by domain");
    auto incompatible=*verified.document;auto a=std::make_shared<AnalysisSnapshot>(*incompatible.currentCache);a->context.platformAndFloatingPointMode="another-platform";a->identity=OfflineAnalysisService::identityFor(a->context,a->domain,a->measuredContentHash);incompatible.currentCache=a;StandaloneSessionModel other;expect(sf::install(incompatible,other,notice),notice);expect(other.source.hasMeasurement&&!other.currentState.hasMeasurement&&other.currentStateNeedsRemeasurement(),"incompatible cache discarded, settings recovered");
    auto wrong=*verified.document;sf::relinkAudio(wrong,root.getChildFile("missing.wav"));auto missing=sf::verify(wrong,c);expect(missing.document&&!missing.document->audioResolved,"missing audio unresolved");expect(sf::install(*missing.document,other,notice),notice);expect(other.sourceContentInvalid&&!other.source.hasMeasurement&&other.sectionDocument.map().digest()==s.sectionDocument.map().digest(),"unresolved group preserves settings and sections without pass");
    const auto unresolvedSnapshot=sf::capture(other);expect(unresolvedSnapshot.sections.digest()==d.sections.digest()&&!unresolvedSnapshot.sourceCache&&!unresolvedSnapshot.currentCache,"saving unresolved audio preserves annotations and omits untrusted caches");
    auto unresolvedFile=root.getChildFile("unresolved-resave.levelscope");expect(sf::save(unresolvedSnapshot,unresolvedFile,c).succeeded(),"unresolved session can be saved");auto unresolvedRead=sf::read(unresolvedFile,c);expect(unresolvedRead.document&&unresolvedRead.document->sections.digest()==d.sections.digest(),"unresolved annotations survive reopening");
    auto moved=root.getChildFile("moved.wav");expect(d.source.file.copyFileTo(moved),"copy same-content relink fixture");

    auto movedDoc=*read.document;sf::relinkAudio(movedDoc,moved);auto movedVerified=sf::verify(movedDoc,c);expect(movedVerified.document&&movedVerified.document->audioResolved,"same content relink verifies");
    expect(sf::install(*movedVerified.document,other,notice)&&other.currentState.hasMeasurement,"same content relink reuses compatible measurements");
    auto changed=bytes(moved);changed.back()^=1;expect(moved.replaceWithData(changed.data(),changed.size()),"change same-name content");auto changedResult=sf::verify(movedDoc,c);expect(changedResult.document&&!changedResult.document->audioResolved,"changed same-name bytes cannot inherit pass");
    auto alias=root.getChildFile("source-alias.levelscope");expect(d.source.file.createSymbolicLink(alias,true),"source alias fixture");
    expect(!sf::save(d,alias,c).succeeded(),"session cannot overwrite aliased media");alias.deleteFile();
    auto afterVerify=sf::verify(*read.document,c);auto originalTime=d.source.file.getLastModificationTime();
    expect(d.source.file.setLastModificationTime(juce::Time(originalTime.toMilliseconds()+2000)),"change verification stamp");
    const auto previousId=other.sessionId;expect(!sf::adopt(*afterVerify.document,*afterVerify.prepared,other,notice)&&other.sessionId==previousId,"late media change cannot replace previous session");
    d.source.file.setLastModificationTime(originalTime);
    auto bad=root.getChildFile("invalid-core.levelscope");
    std::vector<std::function<void(juce::var&)>> invalid {
        [](auto& v){v.getDynamicObject()->setProperty("schema",3);},
        [](auto& v){v.getDynamicObject()->setProperty("documentId",juce::String::repeatedString("z",32));},
        [](auto& v){v["programme"]["channels"][1].getDynamicObject()->setProperty("fileChannel",0);},
        [](auto& v){v["programme"]["layout"]["roles"].getArray()->set(1,1);},
        [](auto& v){v["video"].getDynamicObject()->setProperty("offset","1/0");},
        [](auto& v){v["programme"].getDynamicObject()->setProperty("mode","stems");},
        [](auto& v){v["programme"].getDynamicObject()->setProperty("rangeEnd","1");},
        [](auto& v){v["programme"]["files"][0]["media"].getDynamicObject()->setProperty("frames","1152921504606846975");},
        [](auto& v){v["current"].getDynamicObject()->setProperty("gain","7ff0000000000000");},
        [](auto& v){v["current"].getDynamicObject()->setProperty("engine","unknown-processor");},
        [](auto& v){auto a=v["sections"]["rows"][1];a.getArray()->set(0,v["sections"]["rows"][0][0]);},
        [](auto& v){v["sections"].getDynamicObject()->setProperty("input",juce::String::repeatedString("a",64));}
    };
    for(auto edit:invalid){coreEdit(f,bad,edit);auto rejected=sf::read(bad,c);expect(!rejected.document,"invalid core rejected");}
    coreEdit(f,bad,[](auto& v){v["caches"].getArray()->set(1,juce::var("damaged cache entry"));});auto damagedDirectory=sf::read(bad,c);expect(damagedDirectory.document&&damagedDirectory.document->sourceCache&&!damagedDirectory.document->currentCache,"independently malformed cache entry recoverable");
    coreEdit(f,bad,[](auto& v){v["profile"].getDynamicObject()->setProperty("version","future-profile");});auto profileRead=sf::read(bad,c);expect(profileRead.document.has_value(),"unknown profile retains valid document");auto profileVerified=sf::verify(*profileRead.document,c);expect(profileVerified.document&&sf::install(*profileVerified.document,other,notice),"unknown profile opens");expect(other.restoredProfileVersion=="future-profile"&&other.evaluation(MeasurementDomain::currentState).verdict!=proposal::Verdict::pass,"unknown profile never passes");
    auto truncated=previous.substr(0,previous.size()-17);expect(bad.replaceWithData(truncated.data(),truncated.size()),"truncate last cache");auto truncRead=sf::read(bad,c);expect(truncRead.document&&truncRead.document->sourceCache&&!truncRead.document->currentCache,"truncated cache recovers core and other cache");
    auto coreDamage=previous;coreDamage[81]^=1;expect(bad.replaceWithData(coreDamage.data(),coreDamage.size()),"damage core");expect(!sf::read(bad,c).document,"damaged core checksum rejected");
    auto coherent=sf::capture(s);sf::FileHooks during;during.reached=[&](const char* x){if(juce::String(x)=="replace")s.setCurrentGainDb(-9.123456789012345);};
    expect(sf::save(coherent,f,c,during).succeeded(),"concurrent edit save succeeds");expect(sf::saved(s,coherent,f)&&s.documentDirty(),"concurrent edit remains dirty after atomic commit");auto coherentRead=sf::read(f,c);expect(coherentRead.document&&coherentRead.document->current.staticGainDb==coherent.current.staticGainDb,"file contains captured coherent gain");
    auto beforeMeasure=sf::capture(s);measure(s,MeasurementDomain::currentState);expect(sf::saved(s,beforeMeasure,f)&&s.documentDirty(),"background measurement completion remains unsaved");
    s.newDocument();expect(!sf::saved(s,beforeMeasure,f),"late save cannot mark newly opened document clean");
    auto oversized=root.getChildFile("oversized.levelscope");{std::ofstream huge(oversized.getFullPathName().toStdString(),std::ios::binary);huge.seekp(sf::maximumBytes);huge.put('x');}
    expect(!sf::read(oversized,c).document,"oversized document rejected before allocation");oversized.deleteFile();
    JobControl stop;stop.cancel=true;expect(sf::read(f,stop).error==ErrorCode::cancelled,"open cancellation before installation");

}
uint64_t rss(){mach_task_basic_info_data_t info{};mach_msg_type_number_t count=MACH_TASK_BASIC_INFO_COUNT;return task_info(mach_task_self(),MACH_TASK_BASIC_INFO,reinterpret_cast<task_info_t>(&info),&count)==KERN_SUCCESS?info.resident_size:0;}
double now(){return juce::Time::getMillisecondCounterHiRes();}
void longSnapshot(const juce::File& root)
{
    auto s=model(fixture(root,2,48000),2);s.sourceDocument.media.frames=int64_t(48000)*10800;
    s.synchronizeSections();
    juce::String error;for(int n=1;n<1000;++n){expect(s.sectionDocument.split(uint64_t(n),s.sourceDocument.media.frames*n/1000,error),error);}s.documentChanged();
    for(auto domain:{MeasurementDomain::source,MeasurementDomain::currentState})
    {
        auto& summary=domain==MeasurementDomain::source?s.source:s.currentState;auto a=std::make_shared<AnalysisSnapshot>(*summary.analysis);a->context=s.measurementContext(domain);a->identity=OfflineAnalysisService::identityFor(a->context,domain,a->measuredContentHash);
        a->detail.frames=s.sourceDocument.media.frames;
        for(auto* windows:{&a->detail.momentary,&a->detail.shortTerm,&a->detail.lraWindows})
        {
            const int64_t length=windows==&a->detail.momentary?19200:144000;
            const int64_t step=windows==&a->detail.lraWindows?48000:4800;
            windows->clear();for(int64_t end=length;end<=a->detail.frames;end+=step)windows->push_back({end,.0002});
        }
        summary.analysis=a;
    }
    const auto baseline=rss();std::atomic<bool> stop{false};std::atomic<uint64_t> peak{baseline};
    std::thread sampler([&]{while(!stop){auto n=rss();if(n>peak)peak=n;std::this_thread::sleep_for(std::chrono::milliseconds(1));}});
    struct Join{std::atomic<bool>& stop;std::thread& t;~Join(){stop=true;t.join();}} join{stop,sampler};
    auto t=now();auto d=sf::capture(s);double captureMs=now()-t;auto f=root.getChildFile("three-hour.levelscope");JobControl c;
    struct Pulse final : juce::Timer { std::vector<double> gaps;double last=now();void timerCallback() override{auto t=now();gaps.push_back(t-last);last=t;} } pulse;
    pulse.startTimer(5);
    const auto wait=[&](auto& future){while(future.wait_for(std::chrono::milliseconds(0))!=std::future_status::ready)juce::MessageManager::getInstance()->runDispatchLoopUntil(5);return future.get();};
    t=now();auto saving=std::async(std::launch::async,[&]{return sf::save(d,f,c);});auto written=wait(saving);double saveMs=now()-t;expect(written.succeeded(),written.message);
    t=now();auto reading=std::async(std::launch::async,[&]{return sf::read(f,c);});auto read=wait(reading);double readMs=now()-t;
    // Synthetic media metadata is intentional: this test exercises persistence
    // costs, not media validation. Real per-file verification has separate cases.
    auto prepared=std::make_shared<sf::Prepared>();auto staged=*read.document;staged.audioResolved=staged.videoResolved=true;
    t=now();auto preparing=std::async(std::launch::async,[&]{auto ok=sf::install(staged,prepared->model,prepared->notice);if(ok){prepared->model.evaluation(MeasurementDomain::source);prepared->model.evaluation(MeasurementDomain::currentState);}return ok;});
    expect(wait(preparing),"prepare three-hour session off UI thread");double preparationMs=now()-t;
    StandaloneSessionModel opened;t=now();expect(sf::adopt(staged,*prepared,opened,error),error);opened.evaluation(MeasurementDomain::source);opened.evaluation(MeasurementDomain::currentState);double adoptionMs=now()-t;
    pulse.stopTimer();std::sort(pulse.gaps.begin(),pulse.gaps.end());double pulse95=pulse.gaps.empty()?0:pulse.gaps[static_cast<size_t>(.95*(pulse.gaps.size()-1))];expect(read.document&&read.document->sourceCache&&read.document->currentCache,"three-hour caches survive exact parsing");
    expect(read.document->sections.sections.size()==1000,"1000 exact section annotations");expect(read.document->sourceCache->detail.momentary.size()==d.sourceCache->detail.momentary.size(),"full-resolution windows retained");
    JobControl cancel;sf::FileHooks hooks;double requested=0;hooks.reached=[&](const char* stage){if(juce::String(stage)=="write"&&requested==0){requested=now();cancel.cancel=true;}};
    auto cancelled=sf::save(d,root.getChildFile("cancel-long.levelscope"),cancel,hooks);double cancelMs=now()-requested;expect(cancelled.error==ErrorCode::cancelled,"active large save cancellation");
    const auto growth=peak.load()-baseline;
    auto* p=performance.getDynamicObject();p->setProperty("fixture","Synthetic three-hour full-resolution Source/Current caches and 1000 sections; media hashing separately qualified");p->setProperty("capture_ms",captureMs);p->setProperty("worker_preparation_ms",preparationMs);p->setProperty("ui_adoption_ms",adoptionMs);p->setProperty("save_ms",saveMs);p->setProperty("parse_ms",readMs);p->setProperty("active_cancel_ms",cancelMs);p->setProperty("rss_growth_bytes",static_cast<juce::int64>(growth));p->setProperty("session_bytes",f.getSize());p->setProperty("ui_callback_p95_ms",pulse95);p->setProperty("ui_callback_count",static_cast<int>(pulse.gaps.size()));
    expect(saveMs<=5000&&readMs<=5000,"five-second save/parse budget");expect(captureMs<=100&&adoptionMs<=100,"responsive snapshot capture and adoption");expect(cancelMs<=1000,"one-second active cancellation budget");expect(growth<=256ULL*1024*1024,"256 MiB persistence growth budget");expect(f.getSize()<=int64_t(sf::maximumBytes),"128 MiB session cap");expect(pulse.gaps.size()>=5&&pulse95<=100,"UI callback p95 under 100 ms during save/read");
}
void unavailableMetrics(const juce::File& root)
{
    for(bool silence:{false,true})
    {
        auto s=model(fixture(root,2,48000,silence?48000*4:100,silence),2);auto d=sf::capture(s);JobControl c;auto f=root.getChildFile("unavailable.levelscope");expect(sf::save(d,f,c).succeeded(),"unavailable metrics save");auto r=sf::read(f,c);expect(r.document&&r.document->sourceCache&&r.document->currentCache,"unavailable metrics cache preserved");
        const auto& a=s.source.analysis->detail;const auto& b=r.document->sourceCache->detail;
        expect(a.integrated.status==b.integrated.status&&a.integrated.value==b.integrated.value&&a.lra.status==b.lra.status&&a.lra.value==b.lra.value&&a.lraStable==b.lraStable,"silence/short status and optional LRA remain exact, never zero-filled");
    }
}
void longMedia(const juce::File& file,const juce::File& root)
{
    auto started=now();JobControl c;auto m=validateMedia(file,c);expect(m.succeeded(),m.message);double initialMs=now()-started;
    StandaloneSessionModel s;expect(s.setValidatedSource(file,m.media),"long media import");if(m.media.channels==8){auto l=systemILayout();for(int i=0;i<8;++i)s.setChannelRole(i,l.roles[i]);}expect(s.confirmLayout(),"long media roles");expect(s.initializeCurrentStateFromSource(),"long media Current");s.setCurrentGainDb(-3.305074118944);s.synchronizeSections();auto d=sf::capture(s);
    auto saved=root.getChildFile("long-media.levelscope");expect(sf::save(d,saved,c).succeeded(),"long media settings save");auto read=sf::read(saved,c);expect(read.document.has_value(),"long media core read");
    JobControl cancel;auto future=std::async(std::launch::async,[&]{return sf::verify(*read.document,cancel);});
    while(cancel.progress.load()<.01f&&future.wait_for(std::chrono::milliseconds(1))!=std::future_status::ready){}
    auto progress=cancel.progress.load();auto requested=now();cancel.cancel=true;auto interrupted=future.get();double cancelMs=now()-requested;
    expect(progress>0&&interrupted.error==ErrorCode::cancelled&&cancelMs<=1000,"active long-media verification cancellation/progress");
    started=now();auto complete=sf::verify(*read.document,c);double verifyMs=now()-started;
    expect(complete.document&&complete.document->audioResolved&&complete.prepared,"real three-hour media reopened and verified");
    expect(!complete.prepared->model.source.hasMeasurement&&complete.prepared->model.currentStateNeedsRemeasurement(),"settings-only long session stays unverified");
    mediaPerformance=juce::var(new juce::DynamicObject());auto* v=mediaPerformance.getDynamicObject();v->setProperty("source_sha256",m.media.sha256);v->setProperty("frames",m.media.frames);v->setProperty("rate",m.media.sampleRate);v->setProperty("channels",m.media.channels);v->setProperty("initial_verification_ms",initialMs);v->setProperty("reopen_verification_ms",verifyMs);v->setProperty("active_cancel_ms",cancelMs);v->setProperty("cancel_progress",progress);
}
void exchangeOut(const juce::File& directory)
{
    expect(!directory.exists(),"exchange destination must be new");expect(directory.createDirectory().wasOk(),"create exchange folder");auto s=model(fixture(directory,6,48000),6);auto d=sf::capture(s);JobControl c;auto target=directory.getChildFile("transfer.levelscope");auto r=sf::save(d,target,c);expect(r.succeeded(),r.message);
    auto raw=bytes(target);auto damaged=raw;damaged.back()^=1;
    expect(directory.getChildFile("damaged-cache.levelscope").replaceWithData(damaged.data(),damaged.size()),"UI damaged cache fixture");
    damaged=raw;damaged[81]^=1;expect(directory.getChildFile("invalid-core.levelscope").replaceWithData(damaged.data(),damaged.size()),"UI invalid core fixture");
    coreEdit(target,directory.getChildFile("missing-audio.levelscope"),[&](auto& v){v["programme"]["files"][0].getDynamicObject()->setProperty("locator",directory.getChildFile("missing.wav").getFullPathName());});
    auto changed=bytes(d.source.file);changed.back()^=1;expect(directory.getChildFile("changed-audio.wav").replaceWithData(changed.data(),changed.size()),"UI wrong-content relink fixture");
    s.setCurrentGainDb(-2.345678901234567);expect(sf::save(sf::capture(s),directory.getChildFile("unverified.levelscope"),c).succeeded(),"UI unverified fixture");
}
void exchangeIn(const juce::File& directory)
{
    JobControl c;auto r=sf::read(directory.getChildFile("transfer.levelscope"),c);expect(r.document&&r.document->sourceCache&&r.document->currentCache,"foreign-platform cache record remains readable");auto d=*r.document;sf::relinkAudio(d,directory.getChildFile(d.source.file.getFileName()));auto verified=sf::verify(d,c);expect(verified.document&&verified.document->audioResolved,"transferred audio verified per file");
    StandaloneSessionModel s;juce::String message;expect(sf::install(*verified.document,s,message),message);expect(!s.source.hasMeasurement&&!s.currentState.hasMeasurement,"different platform caches rejected conservatively");expect(s.currentStateDocument.staticGainDb==d.current.staticGainDb&&s.sectionDocument.map().digest()==d.sections.digest(),"cross-platform settings and annotations exact");expect(s.currentStateNeedsRemeasurement(),"cross-platform Current unverified");expect(s.currentStateDocument.isInitialized,"cross-platform processing retained");
}

}
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI init;juce::String output,exchangeWrite,exchangeRead,longFile;for(int n=1;n+1<argc;++n){auto arg=juce::String(argv[n]);if(arg=="--output")output=argv[++n];else if(arg=="--exchange-out")exchangeWrite=argv[++n];else if(arg=="--exchange-in")exchangeRead=argv[++n];else if(arg=="--long-media")longFile=argv[++n];}
    auto root=(output.isNotEmpty()?juce::File(output).getParentDirectory():juce::File::getCurrentWorkingDirectory()).getChildFile("levelscope-m8-"+juce::Uuid().toString());root.createDirectory();bool passed=false;juce::String failure;
    try{for(int channels:{2,6,8})roundtrip(root,channels,48000);for(int rate:{44100,96000,192000})roundtrip(root,2,rate);longSnapshot(root);unavailableMetrics(root);if(longFile.isNotEmpty())longMedia(juce::File(longFile),root);if(exchangeWrite.isNotEmpty()){exchangeOut(juce::File(exchangeWrite));if(longFile.isNotEmpty())expect(root.getChildFile("long-media.levelscope").copyFileTo(juce::File(exchangeWrite).getChildFile("long-media.levelscope")),"UI long opening/cancellation fixture");}if(exchangeRead.isNotEmpty())exchangeIn(juce::File(exchangeRead));passed=true;}catch(const std::exception& e){failure=e.what();std::cerr<<failure<<'\n';}
    auto v=juce::var(new juce::DynamicObject());auto* o=v.getDynamicObject();o->setProperty("suite","M8 sessions");o->setProperty("passed",passed);o->setProperty("checks",checks);o->setProperty("failure",failure);o->setProperty("performance",performance);o->setProperty("long_media",mediaPerformance);o->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);o->setProperty("os",juce::SystemStats::getOperatingSystemName());
#if defined(__arm64__)
    o->setProperty("architecture","arm64");
#else
    o->setProperty("architecture","x86_64");
#endif
    if(output.isNotEmpty())juce::File(output).replaceWithText(juce::JSON::toString(v));std::cout<<juce::JSON::toString(v)<<'\n';root.deleteRecursively();return passed?0:1;
}

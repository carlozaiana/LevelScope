#pragma once
#include <sys/statvfs.h>
namespace
{
constexpr int64_t longFrames=518400000;
double compositionNow(){return juce::Time::getMillisecondCounterHiRes();}
uint64_t compositionRss(){mach_task_basic_info_data_t i{};mach_msg_type_number_t n=MACH_TASK_BASIC_INFO_COUNT;return task_info(mach_task_self(),MACH_TASK_BASIC_INFO,reinterpret_cast<task_info_t>(&i),&n)==KERN_SUCCESS?i.resident_size:0;}
int64_t freeBytes(const juce::File& f){struct statvfs v{};return statvfs(f.getFullPathName().toRawUTF8(),&v)==0?int64_t(v.f_bavail)*v.f_frsize:0;}
template<class T> T compositionAwait(std::future<T>& pending)
{while(pending.wait_for(std::chrono::milliseconds(0))!=std::future_status::ready)juce::MessageManager::getInstance()->runDispatchLoopUntil(10);return pending.get();}
struct CompositionPulse:juce::Timer
{double last=compositionNow();std::vector<double> gaps;void timerCallback()override{const double now=compositionNow();gaps.push_back(now-last);last=now;}};
int longOperand(int component,int channel,int64_t frame)
{
    if(frame==longFrames-1)return component==0?1000+channel*100:500+channel*50;
    const int period=component==0?128:64;
    return static_cast<int>(std::llround((component==0?400:180)*(channel+1)*std::sin((frame%period)*2*3.14159265358979323846/period)));
}
ProgrammeComponent longComponent(const juce::File& root,int index)
{
    ProgrammeComponent c;c.id=index==0?"me":"dialogue";c.purpose=index==0?"music-effects":"dialogue";ProgrammeAssembly group;const auto roles=assemblyRoles(6);JobControl job;
    for(int ch=0;ch<6;++ch)
    {
        const auto file=root.getChildFile(c.id+"-"+juce::String(ch)+".wav");expect(!file.exists(),"long fixture never overwrites existing media");auto out=file.createOutputStream();expect(out!=nullptr,"long mono PCM24 stream");
        const uint32_t bytes=uint32_t(longFrames*3);out->write("RIFF",4);out->writeInt(static_cast<int>(bytes+36));out->write("WAVEfmt ",8);out->writeInt(16);out->writeShort(1);out->writeShort(1);out->writeInt(48000);out->writeInt(144000);out->writeShort(3);out->writeShort(24);out->write("data",4);out->writeInt(static_cast<int>(bytes));
        constexpr int block=65536;std::vector<uint8_t> packed(block*3);
        for(int n=0;n<block;++n){const uint32_t q=uint32_t(longOperand(index,ch,n)*256);for(int b=0;b<3;++b)packed[size_t(n)*3+b]=uint8_t(q>>(b*8));}
        for(int64_t first=0;first<longFrames;first+=block)
        {
            const int count=static_cast<int>(std::min<int64_t>(block,longFrames-first));
            if(first+count==longFrames){const uint32_t q=uint32_t(longOperand(index,ch,longFrames-1)*256);for(int b=0;b<3;++b)packed[size_t(count-1)*3+b]=uint8_t(q>>(b*8));}
            if(!out->write(packed.data(),size_t(count)*3))throw std::runtime_error("long fixture disk write failed");
        }
        out->flush();expect(out->getStatus().wasOk(),"long PCM24 fixture finalized");out.reset();auto v=validateMedia(file,job);expect(v.succeeded(),"long component member independently hashed and inspected");group.members.push_back({"channel-"+juce::String(ch),file,v.media,roles[ch],0});
    }
    expect(group.confirmAlignment(),"long generator confirms exact within-component alignment");c.assembly=group;c.media=group.summary();c.layout=group.layout();return c;
}
void longComposition(const juce::File& root)
{
    const auto begin=compositionNow();const auto initialFree=freeBytes(root);expect(initialFree>=int64_t(65)*1024*1024*1024,"65 GiB free before consolidated long workload");
    CompositionPulse pulse;pulse.startTimer(10);auto stages=juce::var(new juce::DynamicObject());
    auto stage=[&](const char* name,double start){stages.getDynamicObject()->setProperty(name,compositionNow()-start);std::cout<<"LONG "<<name<<" elapsed_s="<<(compositionNow()-begin)/1000<<std::endl;};
    double t=compositionNow();auto generating=std::async(std::launch::async,[&]{ProgrammeComposition c;c.preset=CompositionPreset::surroundSurround;c.me=longComponent(root,0);c.dialogue=longComponent(root,1);expect(c.confirmComplementaryBalancedAlignment(),"long complementary balanced generator recipe confirmed");return c;});auto recipe=compositionAwait(generating);stage("generate_verify_ms",t);
    auto s=compositionModel(recipe);const auto baseline=compositionRss();std::atomic<uint64_t> peak{baseline};std::atomic<int64_t> minimumFree{freeBytes(root)};std::atomic<bool> stop{false};
    std::thread sampler([&]{while(!stop){peak=std::max(peak.load(),compositionRss());minimumFree=std::min(minimumFree.load(),freeBytes(root));std::this_thread::sleep_for(std::chrono::milliseconds(10));}});
    struct Join{std::atomic<bool>& stop;std::thread& thread;~Join(){stop=true;thread.join();}}join{stop,sampler};
    t=compositionNow();auto source=std::async(std::launch::async,[&]{return measuredComposition(s,MeasurementDomain::source);});compositionAwait(source);stage("source_measure_ms",t);
    expect(s.setCurrentGainDb(-23-*s.source.integratedLufs),"explicit long Current normalization gain");
    t=compositionNow();auto current=std::async(std::launch::async,[&]{return measuredComposition(s,MeasurementDomain::currentState);});compositionAwait(current);stage("current_measure_ms",t);
    expect(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"long actual Current passes selected rules");
    s.synchronizeSections();juce::String error;for(int n=1;n<1000;++n)expect(s.sectionDocument.split(n,longFrames*n/1000,error),"long composition section split");s.documentChanged();auto captured=sf::capture(s);JobControl job;
    const auto session=root.getChildFile("composition.levelscope");t=compositionNow();auto saving=std::async(std::launch::async,[&]{return sf::save(captured,session,job);});expect(compositionAwait(saving).succeeded(),"long coherent composed save");const double saveMs=compositionNow()-t;stage("save_ms",t);
    t=compositionNow();auto reading=std::async(std::launch::async,[&]{return sf::read(session,job);});auto read=compositionAwait(reading);const double parseMs=compositionNow()-t;stage("parse_ms",t);
    expect(read.document&&read.document->sourceCache&&read.document->currentCache&&read.document->sections.sections.size()==1000,"long caches and annotations retained");
    t=compositionNow();auto opening=std::async(std::launch::async,[&]{return sf::verify(*read.document,job);});auto reopened=compositionAwait(opening);expect(reopened.prepared&&reopened.document->audioResolved&&reopened.prepared->model.source.hasMeasurement&&reopened.prepared->model.currentState.hasMeasurement,"all long component members reverified and measured caches restored");stage("reopen_verify_ms",t);reopened.prepared.reset();read.document.reset();
    JobControl cancel;auto preparing=std::async(std::launch::async,[&]{return prepareVerifiedSource(s.measurementContext(MeasurementDomain::source),cancel);});
    while(cancel.progress.load()==0&&preparing.wait_for(std::chrono::milliseconds(1))!=std::future_status::ready)juce::MessageManager::getInstance()->runDispatchLoopUntil(10);
    t=compositionNow();cancel.cancel=true;const auto cancelled=compositionAwait(preparing);const double cancelMs=compositionNow()-t;expect(cancelled.error==ErrorCode::cancelled&&cancelMs<=1000,"active long composition preparation cancellation <=1 second");
    auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));OutputRouting muted;s.selectAuditionDomain(MeasurementDomain::currentState);auto request=s.auditionRequest();audition.synchronize(request,muted);
    t=compositionNow();expect(audition.prepare(request).isEmpty(),"long composition device-closed preparation");
    const auto deadline=compositionNow()+1800000;while(audition.preparing()&&compositionNow()<deadline){audition.poll();juce::MessageManager::getInstance()->runDispatchLoopUntil(10);}
    expect(audition.preparedCopy()&&!injected->opened,"complete private component copies before any injected output");const auto* copy=audition.preparedCopy();int64_t payload=0;for(const auto& f:copy->files())payload+=f.getSize();expect(payload==s.sourceDocument.media.bytes,"no extra interleaved composed temporary payload");
    expect(audition.connect(48000).isEmpty(),"long injected output only");
    for(auto domain:{MeasurementDomain::source,MeasurementDomain::currentState})for(int64_t first:{int64_t(0),int64_t(259200019),longFrames-4096})
    {
        s.selectAuditionDomain(domain);request=s.auditionRequest();audition.synchronize(request,muted);expect(audition.seek(first,request,muted).isEmpty(),"long exact audio-source seek");Capture capturedAudio(6,4096,first);audition.setCaptureObserver(&capturedAudio);
        expect(audition.play(request,muted,false,true).isEmpty(),"long composed all-muted injected playback");awaitPlaying(audition,*injected);juce::AudioBuffer<float> block(maxChannels,256);
        for(int n=0;n<16;++n){injected->tick(block);audition.poll();std::this_thread::sleep_for(std::chrono::milliseconds(1));}audition.stop();audition.setCaptureObserver(nullptr);
        bool exact=capturedAudio.valid&&capturedAudio.next==first+4096;const double coefficient=std::pow(10.0,(domain==MeasurementDomain::source?0:s.currentStateDocument.staticGainDb)/20.0);
        for(int ch=0;ch<6;++ch)for(int n=0;n<4096;++n){const float sum=static_cast<float>(double(longOperand(0,ch,first+n))/32768+double(longOperand(1,ch,first+n))/32768);const float expected=static_cast<float>(double(sum)*coefficient);exact&=capturedAudio.data.getSample(ch,n)==expected;}
        expect(exact&&audition.preparedCopy()==copy&&audition.underruns()==0,"long Source/Current first/middle/final captures exact; private copies reused");
    }
    audition.disconnect();stage("prepare_seek_capture_ms",t);
    // Release audition copies before the single full RF64 rendering workflow.
    audition.sourceChanged(s.sourceDocument.media);
    const auto exportRequest=ex::capture(s);const auto destination=root.getChildFile("export");t=compositionNow();auto rendering=std::async(std::launch::async,[&]{return ex::render(exportRequest,root,job);});auto rendered=compositionAwait(rendering);
    if(!rendered.value)std::cerr<<rendered.message<<'\n';expect(rendered.value.has_value(),"one real long composition RF64 render and finalized-artifact remeasurement");stage("render_verify_ms",t);
    const auto candidate=*rendered.value;expect(candidate->file().getSize()==9331200104LL&&candidate->spec().format.rf64&&ex::artifactEvaluation(*candidate,exportRequest).verdict==proposal::Verdict::pass,"real RF64 size and independently measured artifact pass");
    t=compositionNow();auto publishing=std::async(std::launch::async,[&]{return ex::publish(*candidate,exportRequest,destination,ex::approvalIdentity(*candidate,exportRequest),job);});expect(compositionAwait(publishing).value.has_value(),"explicit long local export with verified reports");stage("publish_ms",t);
    pulse.stopTimer();std::sort(pulse.gaps.begin(),pulse.gaps.end());expect(pulse.gaps.size()>=10,"long UI event loop observed");const double p95=pulse.gaps[size_t((pulse.gaps.size()-1)*.95)];const uint64_t growth=peak.load()-baseline;const auto reports=destination.getChildFile("report.json").getSize()+destination.getChildFile("report.html").getSize();
    expect(saveMs<=5000&&parseMs<=5000&&session.getSize()<=int64_t(sf::maximumBytes),"long save/parse/session budgets");expect(growth<=256ULL*1024*1024&&p95<=100&&reports<=16*1024*1024,"long RSS UI and report budgets");expect(compositionNow()-begin<=10800000,"one complete consolidated three-hour workflow within 180 minutes");
    longEvidence=juce::var(new juce::DynamicObject());auto* e=longEvidence.getDynamicObject();e->setProperty("frames",juce::String(longFrames));e->setProperty("sample_rate",48000);e->setProperty("channels",6);e->setProperty("components",2);e->setProperty("annotations",1000);e->setProperty("capture_cases",6);e->setProperty("stages_ms",stages);e->setProperty("total_ms",compositionNow()-begin);e->setProperty("cancel_ms",cancelMs);e->setProperty("ui_p95_ms",p95);e->setProperty("ui_callbacks",int(pulse.gaps.size()));e->setProperty("rss_growth_bytes",juce::int64(growth));e->setProperty("disk_growth_bytes",initialFree-minimumFree.load());e->setProperty("session_bytes",session.getSize());e->setProperty("report_bytes",reports);e->setProperty("artifact_bytes",candidate->file().getSize());e->setProperty("artifact_sha256",candidate->hash());e->setProperty("composition_identity",recipe.definitionIdentity());
}
}

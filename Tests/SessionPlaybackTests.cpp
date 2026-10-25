#include "AuditionTestSupport.h"
#include "Standalone/SessionDocument.h"
#include "Standalone/StandaloneMeasurementLayout.h"
namespace sf = levelscope::standalone::sessionfile;
void beginVideoReaderAudit();
int videoReaderAuditCount(int);
namespace
{
int completedReferences=0;
void acceptMeasured(StandaloneSessionModel& s,MeasurementDomain domain,ProcessedAudioObserver* capture=nullptr)
{
    expect(s.job.begin(JobStage::analysis,s.revision),"session playback begin measurement");auto c=s.measurementContext(domain);JobControl control;
    auto result=OfflineAnalysisService(capture).analyze(c,domain,control);expect(result.value.has_value(),"offline measurement succeeds");expect(s.job.finish(c.ticket,s.revision,result.error,result.message),"measurement finish");expect(s.acceptMeasurement(std::move(*result.value)),"measurement accepted");
}
StandaloneSessionModel sessionFor(const AuditionRequest& input)
{
    StandaloneSessionModel s;expect(s.setValidatedSource(input.context.sourceFile,input.context.media),"session source");
    auto l=input.context.fileLayout;if(l.count==8)l=systemILayout();for(int c=0;c<l.count;++c)s.setChannelRole(c,l.roles[c]);expect(s.confirmLayout(),"session roles");expect(s.initializeCurrentStateFromSource(),"session current");return s;
}
StandaloneSessionModel reopened(StandaloneSessionModel& s,const juce::File& file)
{
    JobControl c;auto d=sf::capture(s);auto save=sf::save(d,file,c);expect(save.succeeded(),"playback session saved");auto read=sf::read(file,c);expect(read.document.has_value(),"playback session read");auto verified=sf::verify(*read.document,c);expect(verified.document&&verified.document->audioResolved&&verified.document->videoResolved,"playback session files verified");StandaloneSessionModel n;juce::String error;expect(sf::install(*verified.document,n,error),"playback session installed");return n;
}
void loadingGuards(const juce::File& root)
{
    auto input=fixture(root.getChildFile("loading.wav"),48000,2,8192);auto s=sessionFor(input);
    auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));StandaloneAudioPanel panel(s,audition,[]{});
    s.persistenceLoading=true;panel.refresh();
    for(const char* label:{"Open device","Prepare audition","Play selected","Restart","Seek"})
    {
        juce::TextButton* found=nullptr;for(auto* child:panel.getChildren())if(auto* b=dynamic_cast<juce::TextButton*>(child))if(b->getButtonText()==label)found=b;
        expect(found&&!found->isEnabled(),"opening disables audition lifecycle controls");found->onClick();
        expect(injected->opens==0&&!audition.active()&&!audition.preparing(),"queued opening action cannot open a device or start preparation/playback");
    }
    const auto gain=s.currentStateDocument.staticGainDb;expect(!s.setCurrentGainDb(1)&&s.currentStateDocument.staticGainDb==gain,"opening blocks processing changes");
    s.persistenceLoading=false;panel.refresh();expect(audition.connect(48000).isEmpty(),"explicit injected connection after load guard clears");audition.disconnect();
}
void sampleParity(const juce::File& root)
{
    for(int rate:{44100,48000,96000,192000})for(int channels:{2,6,8})for(double gain:{0.,-6.020599913279624,6.020599913279624,24.})
    {
        const int frames=8192;auto input=fixture(root.getChildFile("sample.wav"),rate,channels,frames);auto s=sessionFor(input);s.setCurrentGainDb(gain);Capture offline(channels,frames);acceptMeasured(s,MeasurementDomain::currentState,&offline);
        auto restored=reopened(s,root.getChildFile("sample.levelscope"));expect(restored.currentState.hasMeasurement,"current cache valid after reopen");restored.selectAuditionDomain(MeasurementDomain::currentState);
        auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));expect(injected->opens==0&&!audition.connected(),"reopen never opens a device");
        const int first=gain==0?0:1234;Capture captured(channels,frames-first,first);audition.setCaptureObserver(&captured);audition.sourceChanged(restored.sourceDocument.media);auto request=restored.auditionRequest();audition.synchronize(request,restored.sourceDocument.routing);
        expect(audition.connect(rate).isEmpty(),"explicit injected device only");expect(audition.seek(first,request,restored.sourceDocument.routing).isEmpty(),"restored audition exact seek");expect(audition.play(request,restored.sourceDocument.routing,false,true).isEmpty(),"restored audition all channels muted");awaitPlaying(audition,*injected);drain(audition,*injected,256);
        expect(captured.valid&&captured.next==frames,"all restored pre-monitor samples captured");for(int ch=0;ch<channels;++ch)expect(std::memcmp(captured.data.getReadPointer(ch),offline.data.getReadPointer(ch,first),static_cast<size_t>(frames-first)*4)==0,"restored processing bit-exact against qualified offline path");
        recordSamples(captured,rate,gain,"M8 restored Current",restored.currentState.analysis->detail);
        audition.disconnect();expect(!injected->opened,"injected lifecycle closed");
        input.context.sourceFile.deleteFile();
    }
}
void referenceRoundtrip(const juce::File& root,const juce::File& videos,const juce::File& prores)
{
    auto input=fixture(root.getChildFile("picture-audio.wav"),48000,2,48000*4);auto s=sessionFor(input);acceptMeasured(s,MeasurementDomain::source);acceptMeasured(s,MeasurementDomain::currentState);beginVideoReaderAudit();
    std::vector<juce::File> files{videos.getChildFile("cfr-exact.mov"),videos.getChildFile("vfr-exact.mov")};for(const char* name:{"proxy-with-tone.mov","lt-with-tone.mov","422-with-tone.mov"})files.push_back(prores.getChildFile(name));
    for(const auto& f:files)
    {
        auto decoder=video::makeDecoder();auto asset=decoder->load(f,0,[]{return false;});expect(asset!=nullptr,"reference asset available");
        s.referenceVideo={f,asset->sha256,asset->selectedStream,video::Time(-1001,30000)};auto restored=reopened(s,root.getChildFile("video.levelscope"));
        expect(restored.source.analysis->identity==s.source.analysis->identity&&restored.currentState.analysis->identity==s.currentState.analysis->identity,"restored reference preserves audio identities");expect(restored.referenceVideo.startOnAudio==s.referenceVideo.startOnAudio,"restored signed rational offset");
        video::Player player;player.attach(restored.referenceVideo.file,restored.referenceVideo.streamId);uint64_t generation=1;
        for(auto audioTime:{video::Time(0,1),video::Time(3,2),video::Time(3,1),video::Time(0,1)})
        {
            const auto pictureTime=audioTime-restored.referenceVideo.startOnAudio;auto expected=video::frameAt(asset->frames,pictureTime);auto deadline=juce::Time::getMillisecondCounterHiRes()+10000;bool found=false;++generation;player.invalidate();
            while(juce::Time::getMillisecondCounterHiRes()<deadline){auto view=player.update(pictureTime,generation);if(view.asset&&((expected&&view.picture&&view.picture->index==*expected)||(!expected&&!view.picture))){found=true;break;}std::this_thread::sleep_for(std::chrono::milliseconds(2));}
            expect(found,"restored stopped seek follows exact audio source coordinates");
        }
        player.clear();++completedReferences;
    }
    expect(videoReaderAuditCount(0)>0&&videoReaderAuditCount(1)==0,"embedded video audio never decoded by reopening");
}
}
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI gui;juce::String report,video,prores,samples;for(int n=1;n+1<argc;++n){auto a=juce::String(argv[n]);if(a=="--report")report=argv[++n];else if(a=="--video")video=argv[++n];else if(a=="--prores")prores=argv[++n];else if(a=="--samples")samples=argv[++n];}
    bool passed=false;juce::String error;try{if(samples.isNotEmpty()){sampleDirectory=juce::File(samples);expect(!sampleDirectory.exists()&&sampleDirectory.createDirectory().wasOk(),"fresh sample output");}Temp tmp;loadingGuards(tmp.directory);sampleParity(tmp.directory);if(video.isNotEmpty()&&prores.isNotEmpty())referenceRoundtrip(tmp.directory,juce::File(video),juce::File(prores));passed=true;}catch(const std::exception& e){error=e.what();std::cerr<<error<<'\n';}
    auto v=juce::var(new juce::DynamicObject());auto* o=v.getDynamicObject();o->setProperty("suite","M8 restored audition and reference picture");o->setProperty("passed",passed);o->setProperty("error",error);o->setProperty("checks",checks);o->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);o->setProperty("os",juce::SystemStats::getOperatingSystemName());o->setProperty("physical_devices_opened",false);o->setProperty("sample_cases",cases);o->setProperty("reference_cases",completedReferences);
#if defined(__arm64__)
    o->setProperty("architecture","arm64");
#else
    o->setProperty("architecture","x86_64");
#endif
    if(report.isNotEmpty())juce::File(report).replaceWithText(juce::JSON::toString(v));std::cout<<"M8 playback "<<(passed?"PASS":"FAIL")<<" ("<<checks.size()<<" checks)\n";return passed?0:1;
}

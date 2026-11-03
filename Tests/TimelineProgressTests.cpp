#include "AuditionTestSupport.h"
#include "Standalone/MeasurementTimeline.h"
#include "Standalone/StandaloneJobRequirements.h"
#include "Standalone/TimelineTransport.h"
#include "Standalone/TimelineLayout.h"
#include <future>

namespace
{
using Domain=MeasurementDomain;
void imageOf(juce::Component& component,const juce::File& path)
{
    auto stream=path.createOutputStream();expect(stream&&juce::PNGImageFormat().writeImageToStream(component.createComponentSnapshot(component.getLocalBounds()),*stream),"write rendered interaction evidence");
}
juce::MouseEvent mouse(juce::Component& c,juce::Point<float> p,bool left=true)
{
    const auto t=juce::Time::getCurrentTime();return {juce::Desktop::getInstance().getMainMouseSource(),p,juce::ModifierKeys(left?juce::ModifierKeys::leftButtonModifier:0),0,0,0,0,0,&c,&c,t,p,t,1,false};
}
void click(MeasurementTimeline& ui,juce::Point<float> p){ui.mouseDown(mouse(ui,p));ui.mouseUp(mouse(ui,p,false));}
juce::TextButton* button(juce::Component& c,const juce::String& name)
{for(auto* child:c.getChildren())if(auto* b=dynamic_cast<juce::TextButton*>(child);b&&b->getButtonText()==name)return b;return nullptr;}
void reducedWindows()
{
    constexpr int rate=48000,step=4800;LoudnessPreview preview(10800LL*rate,rate);
    expect(preview.storageBytes()<=LoudnessPreview::maximumBytes,"three-hour M/S preview fits its 4 MiB allocation bound");
    for(int k=4;k<=31;++k)
    {
        preview.window(false,{int64_t(k)*step,k==6?0.:.01*k});if(k>=30)preview.window(true,{int64_t(k)*step,.002*k});
        if(k==5)preview.publish(5*step,true);
    }
    const auto old=preview.coverage();expect(!preview.inspect(true,old,old),"no short-term window before three seconds");
    expect(!preview.inspect(false,old+1,old),"unread time is not populated by nearest earlier window");
    preview.publish(31*step,true);expect(preview.reduce(0,31*step,1,old).fullResolutionPoints==2,"old publication stays immutable after later measured windows arrive");
    auto small=preview.reduce(0,31*step,1,preview.coverage());const auto& bucket=small.momentary[0];
    expect(!bucket.included.count&&!bucket.excluded.count&&bucket.ungated.count==27&&small.shortTerm[0].ungated.count==2,"pending gate never classified; zero energy omitted; exact S cadence");
    expect(bucket.ungated.minimum==measurement::OfflineMeasurement::loudness(.04)&&bucket.ungated.maximum==measurement::OfflineMeasurement::loudness(.31),"reduced preview preserves both exact measured loudness extrema");
    auto wide=preview.reduce(0,31*step,100000,preview.coverage());expect(wide.momentary.size()==4096,"bounded preview columns");
    proposal::ValueExtent previous;bool gap=false,joined=false;
    for(const auto& col:wide.momentary)if(col.ungated.count){if(previous.count){if(col.ungated.firstEnd==7*step)gap=!col.ungated.follows(previous,step);else joined|=col.ungated.follows(previous,step);}previous=col.ungated;}
    expect(gap&&joined,"curve joins adjacent real windows without joining across silence");
    expect(preview.inspect(false,4.5*step,preview.coverage())->endFrame==4*step,"nearest preview window uses earlier tie");
    for(int k=32;k<=108000;++k){preview.window(false,{int64_t(k)*step,.001*(1+k%97)});preview.window(true,{int64_t(k)*step,.002*(1+k%43)});}
    preview.publish(preview.frames,true);std::vector<double> timing;
    for(int n=0;n<20;++n){const auto started=juce::Time::getMillisecondCounterHiRes();auto full=preview.reduce(0,double(preview.frames),4096,preview.frames);timing.push_back(juce::Time::getMillisecondCounterHiRes()-started);expect(full.fullResolutionPoints==215968,"synthetic three-hour display retains every complete M and S window");}
    std::sort(timing.begin(),timing.end());expect(timing[18]<50,"three-hour-sized synthetic viewport preparation p95 below 50 ms");
    auto* perf=new juce::DynamicObject();perf->setProperty("synthetic_three_hour_preview_bytes",juce::int64(preview.storageBytes()));perf->setProperty("synthetic_prepare_p95_ms",timing[18]);cases.add(juce::var(perf));
    WaveformBinding binding{"runtime","programme",{4,7},Domain::source};LoudnessProgress progress(binding);auto held=progress.begin(48000,rate);
    expect(held&&progress.view(binding),"job-bound preview available");
    for(int n=0;n<4;++n){auto stale=binding;if(n==0)++stale.ticket.id;if(n==1)++stale.ticket.sessionRevision;if(n==2)stale.runtime="reopened";if(n==3)stale.identity="other-processing";expect(!progress.view(stale),"stale runtime, job, revision or processing cannot view publication");}
    progress.abandon();held->window(false,{4*step,.1});held->publish(4*step,true);expect(!progress.view(binding),"late retired worker cannot republish");
    LoudnessProgress oversized(binding);expect(!oversized.begin(86400LL*rate,rate),"out-of-budget preview refused independently of measurement");
}
void concurrentPublication()
{
    constexpr int rate=48000,step=4800,steps=12000;LoudnessPreview data(int64_t(steps)*step,rate);
    std::atomic<bool> done{false};std::promise<void> first,start;auto ready=first.get_future();auto go=start.get_future();
    std::thread writer([&]{for(int k=4;k<=steps;++k){data.window(false,{int64_t(k)*step,.02});if(k>=30)data.window(true,{int64_t(k)*step,.01});
        if(k%8==0)data.publish(int64_t(k)*step,true);if(k==32){first.set_value();go.wait();}}
        data.publish(data.frames,true);done.store(true);});
    ready.wait();const auto initial=data.coverage();start.set_value();bool consistent=true;int reads=0;
    do {const auto token=data.coverage();const auto plot=data.reduce(0,double(data.frames),83,token);const auto k=token/step;
        consistent&=plot.fullResolutionPoints==size_t(k-3+k-29)&&!data.inspect(true,double(token+1),token);++reads;}while(!done.load());writer.join();
    expect(consistent&&reads>0&&data.coverage()==data.frames,"reader sees immutable, complete prefixes under concurrent window publication");
    expect(data.reduce(0,double(data.frames),1,initial).fullResolutionPoints==32,"old coverage token remains stable after final worker publication");
}
struct Pause final:ProcessedAudioObserver
{
    std::promise<void> ready,release;std::future<void> released=release.get_future();int stage=0;
    void observe(const float* const*,int,int,int64_t first) noexcept override
    {
        // Pause on the test worker only so the UI can examine an actual published prefix.
        if(stage==0&&first>=48000*3){std::this_thread::sleep_for(std::chrono::milliseconds(215));stage=1;}
        else if(stage==1){stage=2;ready.set_value();released.wait();}
    }
};
void partial(StandaloneSessionModel& s,const juce::File& root)
{
    expect(s.job.begin(JobStage::analysis,s.revision),"start analysis for cancellation evidence");auto c=s.measurementContext(Domain::source);const auto binding=WaveformBinding::from(c,Domain::source);
    LoudnessProgress progress(binding);WaveformCache cache;auto pcm=cache.begin(binding);JobControl control;Pause pause;auto ready=pause.ready.get_future();
    auto future=std::async(std::launch::async,[&]{return OfflineAnalysisService(&pause,pcm.get(),&progress).analyze(c,Domain::source,control);});
    const bool reached=ready.wait_for(std::chrono::seconds(10))==std::future_status::ready;
    struct Release {Pause& pause;bool done=false;~Release(){if(!done)pause.release.set_value();}} release{pause};
    expect(reached,"actual measured M/S available on worker before EOF");auto data=progress.view(binding);expect(data&&data->coverage()>=144000&&data->coverage()<c.media.frames,"sealed real prefix includes 3 s S but not EOF");
    MeasurementTimeline ui;ui.setSize(1150,470);ui.setWaveform(pcm->overview(),false,{"L","R"});ui.setPreview(data);ui.onRequestSeek=[](int64_t,const juce::String&){return juce::String("Unexpected seek");};ui.setTransport(0,false,"Finish analysis before seeking.",binding.identity);
    expect(ui.inspectSeconds(3)&&ui.inspectionText().contains("not a programme result")&&ui.inspectionText().contains("gate pending")&&ui.inspectionText().contains("ungated S"),"partial inspection distinguishes measured M/S from final gates/results");
    expect(ui.inspectSeconds(4.5)&&ui.inspectionText().contains("Pending coverage"),"unfinished passage explicitly pending");
    expect(!s.source.hasMeasurement&&!s.currentState.hasMeasurement,"partial preview does not initialize Current or create measured authority");
    imageOf(ui,root.getChildFile("partial.png"));auto started=juce::Time::getMillisecondCounterHiRes();control.cancel=true;pause.release.set_value();release.done=true;
    auto result=future.get();s.job.finish(c.ticket,s.revision,result.error,result.message);
    expect(!result.value&&result.error==ErrorCode::cancelled&&!progress.view(binding)&&!pcm->overview(),"cancellation discards both partial streams and no final snapshot");
    expect(juce::Time::getMillisecondCounterHiRes()-started<1000,"cancelled worker acknowledges within one second");
}
std::shared_ptr<const proposal::ProgramAnalysis> measure(StandaloneSessionModel& s,Domain d,std::shared_ptr<const waveform::Overview>* pcm=nullptr)
{
    expect(s.job.begin(JobStage::analysis,s.revision),"start fresh independent job");const auto c=s.measurementContext(d);LoudnessProgress progress(WaveformBinding::from(c,d));JobControl job;
    WaveformProgress waveform(std::make_shared<waveform::Budget>());auto result=OfflineAnalysisService(nullptr,pcm?&waveform:nullptr,&progress).analyze(c,d,job);if(pcm)*pcm=waveform.overview();expect(result.value.has_value(),"new measurement completes after cancellation");
    expect(s.job.finish(c.ticket,s.revision,result.error,result.message)&&s.acceptMeasurement(*result.value),"only accepted final snapshot provides programme result");
    return proposal::analyzeProgramme(result.value->detail,result.value->identity);
}
void navigation(StandaloneSessionModel& s,const juce::File& root)
{
    measure(s,Domain::source);expect(s.currentState.hasMeasurement,"initial identical Current ready only after final Source acceptance");
    expect(s.setCurrentGainDb(-4),"change independently inspected Current");std::shared_ptr<const waveform::Overview> currentPcm;const auto current=measure(s,Domain::currentState,&currentPcm);
    auto device=std::make_unique<Device>();auto* io=device.get();auto audio=std::make_unique<StandaloneAudition>(std::move(device));
    s.auditionBusy=[&]{return audio->active();};s.selectAuditionDomain(Domain::source);const auto req=s.auditionRequest();auto key=req.identity();const auto revision=s.revision,document=s.documentRevision;
    const auto sourceSnapshot=s.source.analysis,currentSnapshot=s.currentState.analysis;MeasurementTimeline ui;ui.setSize(1150,470);
    ui.setProgramme(req.context.media.frames,48000,"programme");ui.setWaveform({},false,{"L","R"});ui.setAnalysis(current);
    int seeks=0;ui.onRequestSeek=[&](int64_t frame,const juce::String& bound){++seeks;return TimelineTransport::seek(frame,bound,s,*audio);};ui.setTransport(0,true,{},key);
    auto mid=ui.rulerBounds().getCentre().toFloat();const auto graph=ui.graphBounds();const auto startGeneration=audio->transportGeneration();
    click(ui,graph.getCentre().toFloat());expect(seeks==0&&audio->transportGeneration()==startGeneration&&ui.inspectionText().contains("Measured loudness"),"graph click inspects Current without changing Source transport");
    const auto inspected=ui.inspectionText();click(ui,mid);expect(ui.inspectionText()==inspected,"seeking leaves the inspected measurement readout intact");const auto frame=(int64_t(mid.x)-graph.getX())*req.context.media.frames/(graph.getWidth()-1);
    expect(seeks==1&&audio->positionFrames()==frame&&io->opens==0&&!audio->active()&&s.auditionDomain==Domain::source,"stopped ruler seek uses source audition, exact floor and no device open");
    auto right=mid;right.x=float(graph.getRight()-1);click(ui,right);expect(audio->positionFrames()==req.context.media.frames&&audio->getState()==AuditionState::finished,"rightmost ruler pixel seeks exact EOF endpoint");
    auto* seconds=dynamic_cast<juce::TextEditor*>(ui.findChildWithID("seconds"));if(!seconds)for(auto* child:ui.getChildren())if(auto* e=dynamic_cast<juce::TextEditor*>(child))seconds=e;
    expect(seconds&&button(ui,"Seek"),"separate numeric Inspect and Seek controls");seconds->setText("1/7",false);seconds->onReturnKey();expect(seeks==2,"Enter inspects rational seconds without seeking");
    button(ui,"Seek")->onClick();expect(audio->positionFrames()==48000/7&&seeks==3,"explicit rational seek floors once through transport");
    const auto before=audio->transportGeneration();seconds->setText("-0.000001",false);button(ui,"Seek")->onClick();expect(audio->transportGeneration()==before,"negative sub-sample seek rejected");
    ui.setViewRange(48000,144000);auto original=ui.viewRange();auto drag=mid.translated(80,0);ui.mouseDown(mouse(ui,mid));ui.mouseDrag(mouse(ui,drag));ui.mouseUp(mouse(ui,drag,false));
    expect(ui.viewRange()!=original&&seeks==3&&audio->transportGeneration()==before,"ruler drag pans without seeking");
    button(ui,"Zoom +")->onClick();expect(ui.viewRange().second-ui.viewRange().first<original.second-original.first,"zoom changes viewport only");
    button(ui,"Whole programme")->onClick();expect(ui.viewRange().first==0&&ui.viewRange().second==req.context.media.frames,"whole-programme restores full extent");
    juce::ToggleButton* follow=nullptr;for(auto* child:ui.getChildren())if(auto* toggle=dynamic_cast<juce::ToggleButton*>(child);toggle&&toggle->getButtonText()=="Follow")follow=toggle;
    expect(follow!=nullptr,"explicit follow control");ui.setViewRange(0,48000);follow->setToggleState(true,juce::dontSendNotification);ui.setTransport(192000,true,{},key);
    expect(ui.viewRange().first<=192000&&ui.viewRange().second>=192000&&audio->transportGeneration()==before,"follow tracks observed playhead without issuing transport commands");
    juce::MouseWheelDetails wheel{};wheel.deltaX=.15f;ui.mouseWheelMove(mouse(ui,mid),wheel);expect(!follow->getToggleState()&&seeks==3,"manual wheel pan releases follow without seeking");button(ui,"Whole programme")->onClick();
    ui.mouseDown(mouse(ui,mid));ui.setTransport(audio->positionFrames(),false,"Source changed",{});ui.mouseUp(mouse(ui,mid,false));expect(seeks==3,"changed binding invalidates in-flight click");
    click(ui,mid);expect(seeks==3&&ui.navigationStatus().contains("Source changed"),"disabled seek gives clear reason");
    expect(TimelineTransport::seek(23,"stale",s,*audio).isNotEmpty()&&audio->transportGeneration()==before,"stale audition identity rejected before touching service");
    s.persistenceLoading=true;expect(TimelineTransport::seek(3,key,s,*audio).isNotEmpty(),"loading blocks seek");s.persistenceLoading=false;
    s.job.begin(JobStage::waveform,s.revision);expect(TimelineTransport::seek(3,key,s,*audio).isNotEmpty(),"active job blocks seek");s.job.finish(s.job.ticket,s.revision,ErrorCode::none,"done");
    expect(MeasurementTimeline::timecode(48000*(3600+62)+24000,48000)=="01:01:02.500","audio timecode uses programme clock without assumed video rate");
    expect(audio->connect(48000).isEmpty(),"only injected output opens");auto route=routes(2);audio->synchronize(req,route);
    expect(audio->play(req,route).isEmpty()&&audio->active(),"injected prepare/play active");const auto activeGen=audio->transportGeneration();
    expect(TimelineTransport::seek(1,key,s,*audio).isNotEmpty()&&audio->transportGeneration()==activeGen&&audio->active(),"timeline cannot stop or reposition preparing/buffering audition");
    awaitPlaying(*audio,*io);expect(TimelineTransport::seek(1,key,s,*audio).isNotEmpty()&&audio->getState()==AuditionState::playing,"playing seek rejected without stop/restart");
    JobRequirementsPanel targets(s,[]{},root.getChildFile("empty-presets"));juce::ComboBox* profiles=nullptr;for(auto* child:targets.getChildren())if(auto* b=dynamic_cast<juce::ComboBox*>(child))profiles=b;
    auto* save=button(targets,"Save as preset...");expect(profiles&&save&&!profiles->isEnabled()&&!save->isEnabled(),"Target controls intentionally locked during audition");
    audio->stop();targets.refreshAvailability();expect(profiles->isEnabled()&&save->isEnabled(),"Stop alone restores preset controls without touching applicability/draft");
    s.job.begin(JobStage::proposal,s.revision);targets.refreshAvailability();expect(!save->isEnabled(),"active proposal job locks editing");s.job.finish(s.job.ticket,s.revision,ErrorCode::none,"done");targets.refreshAvailability();expect(save->isEnabled()&&profiles->isEnabled(),"job completion alone restores Target editing");
    expect(s.revision==revision&&s.documentRevision==document&&s.source.analysis==sourceSnapshot&&s.currentState.analysis==currentSnapshot&&s.auditionDomain==Domain::source,"navigation preserves inspected/audition domains, revision and both accepted measurements");
    ui.setWaveform(currentPcm,true,{"L","R"});ui.setTransport(audio->positionFrames(),true,{},key);click(ui,graph.getCentre().toFloat());
    const auto preFit=audio->transportGeneration();button(ui,"Fit selection")->onClick();expect(ui.viewRange().second-ui.viewRange().first==currentPcm->binFrames&&audio->transportGeneration()==preFit,"Fit selection zooms to the inspected overview bin without seeking");
    button(ui,"Seek")->onClick();expect(audio->positionFrames()==frame,"explicit Seek uses the exact inspected frame regardless of rounded entry display");button(ui,"Whole programme")->onClick();ui.setTransport(audio->positionFrames(),true,{},key);imageOf(ui,root.getChildFile("navigation.png"));audio->disconnect();s.auditionBusy={};
}
void layout(StandaloneSessionModel& s,const juce::File& root)
{
    juce::Component shell;shell.setSize(1100,120);TimelineDivider divider("Comparison");shell.addAndMakeVisible(divider);divider.setBounds(0,90,1100,8);JobRequirementsPanel panel(s,[]{});shell.addAndMakeVisible(panel);shell.addAndMakeVisible(panel.informationLine());
    TimelineLayout layout;layout.moveComparison(-10000);expect(layout.comparison==90&&TimelineLayout::defaultComparison-layout.comparison==78,"upper divider allows 78 pixels compression from default");
    panel.setSize(1100,layout.comparison);panel.informationLine().setBounds(0,98,1100,TimelineLayout::informationHeight);panel.resized();
    for(auto* c:panel.getChildren())if(auto* field=dynamic_cast<juce::TextEditor*>(c))expect(field->getHeight()>=22&&field->getFont().getHeight()==16.f&&panel.getLocalBounds().contains(field->getBounds()),"compressed target fields retain font size and useful height");
    auto* presetButton=button(panel,"Save as preset...");expect(presetButton&&presetButton->getHeight()==22&&presetButton->getRight()<=230,"compact preset controls occupy the unused header cell without covering Target fields");
    expect(panel.informationLine().getY()==98&&panel.informationLine().getHeight()==22,"information line remains fixed outside compressed comparison above");
    imageOf(shell,root.getChildFile("compact.png"));
    layout.moveTimeline(-130,600,6);expect(layout.inspector==std::min(TimelineLayout::minimumInspector+130,600-233-TimelineLayout::dividerHeight),"lower divider releases waveform height to inspector");layout.moveTimeline(-1000,600,8);expect(600-layout.inspector-8>=137+8*16,"lower divider preserves controls, scrollbar and channel labels");
    layout.reset();expect(layout.inspector==TimelineLayout::minimumInspector&&layout.comparison==168,"Reset layout restores both dividers");
}
}
#include "TimelineGainLaneChecks.h"
#include "Standalone/AutomaticSourceAnalysis.h"
#include "TimelineSectionImportChecks.h"
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI gui;juce::String rootPath,output,error;bool passed=false;
    for(int n=1;n<argc;++n){juce::String arg(argv[n]);if(arg=="--root"&&n+1<argc)rootPath=argv[++n];else if(arg=="--output"&&n+1<argc)output=argv[++n];else return 2;}
    try
    {
        expect(rootPath.isNotEmpty(),"explicit owned fixture root");juce::File root(rootPath);expect(!root.exists()&&root.createDirectory().wasOk(),"create fresh test directory");
        reducedWindows();concurrentPublication();auto request=fixture(root.getChildFile("mix.wav"),48000,2,240013);StandaloneSessionModel s;expect(s.setValidatedSource(request.context.sourceFile,request.context.media),"verified interleaved source");s.sourceDocument.layout=request.context.fileLayout;
        partial(s,root);navigation(s,root);layout(s,root);gainLaneChecks(s,root);transportFocusChecks(s);gainDenseChecks(root);presetDismissal(root);drawSections(root);rememberedOrders(root);passed=true;
    }
    catch(const std::exception& e){error=e.what();std::cerr<<error<<'\n';}
    auto* report=new juce::DynamicObject();report->setProperty("passed",passed);report->setProperty("failure",error);report->setProperty("checks",juce::var(checks));report->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);report->setProperty("physical_devices_opened",false);report->setProperty("performance",juce::var(cases));
    if(output.isNotEmpty())juce::File(output).replaceWithText(juce::JSON::toString(juce::var(report),true));std::cout<<"M18 progress/navigation "<<(passed?"PASS":"FAIL")<<" "<<checks.size()<<" checks\n";return passed?0:1;
}

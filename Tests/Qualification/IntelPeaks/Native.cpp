// Reuse the frozen M18 native-window helpers; only this launcher relaxes access.
#define main m18_unused_native_main
#include "../../Profiling/ShortAcceptance/ShortAcceptance.cpp"
#undef main

namespace
{
struct Response : private juce::Timer
{
    std::vector<double> gaps;
    Clock::time_point previous=Clock::now();
    void begin(){previous=Clock::now();startTimer(100);}
    void end(){stopTimer();gaps.push_back(since(previous)*1000);}
    void timerCallback() override {gaps.push_back(since(previous)*1000);previous=Clock::now();}
    void record(juce::var& r)
    {
        stopTimer();std::sort(gaps.begin(),gaps.end());
        expect(!gaps.empty(),"native message loop heartbeat observed");
        const double p95=gaps[size_t(std::ceil(.95*gaps.size()))-1], maximum=gaps.back();
        put(r,"heartbeat_p95_ms",p95);put(r,"heartbeat_max_ms",maximum);
        expect(p95<=120 && maximum<=200,"native heartbeat within selected budget");
    }
};
void cancel(Main& ui,const char* key,juce::var& r)
{
    until("job did not reach measurement before cancellation",[&]{return ui.session.job.isActive()
        && ui.importControl && ui.importControl->progress.load()>.55f;});
    const auto start=Clock::now();click(ui.cancelJobButton);
    until("cancel acknowledgement exceeded one second",[&]{return !ui.session.job.isActive();},1);
    put(r,key,since(start)*1000);
    expect(ui.session.job.state==JobState::cancelled,"cancelled native job has explicit cancelled state");
    expect(!ui.audition.connected()&&!ui.audition.active(),"cancellation never opens a device or starts audition");
}
void checkNative(Main& ui,const juce::File& file,const juce::File& dir,juce::var& r)
{
    auto& s=ui.session;Response response;
    ui.audioPanel.roleMemory=RememberedChannelRoles(dir.getChildFile("owned-role-memory"));
    ui.chooseSourceFile();choose(ui.sourceFileChooser.get(),file);
    until("roles unavailable",[&]{return s.sourceDocument.hasSource&&ui.audioWindow&&ui.audioWindow->isVisible();});
    ui.jobRequirementsPanel.profiles.setSelectedId(2,juce::sendNotificationSync);
    ui.applicabilityConfirmation.setToggleState(true,juce::sendNotificationSync);
    response.begin();click(ui.audioPanel.confirmButton);
    cancel(ui,"source_cancel_ms",r);response.end();
    expect(!s.source.hasMeasurement&&!s.currentState.hasMeasurement,"cancelled first Source leaves both domains pending");
    response.begin();const auto source=Clock::now();click(ui.measureSourceButton);idle(ui);response.end();
    put(r,"source_seconds",since(source));
    expect(s.source.hasMeasurement&&s.currentState.hasMeasurement,"one successful Source initializes matching Current");
    const auto baseline=s.currentState.analysis;const auto original=s.source.analysis;
    expect(ui.timeline.waveformAccepted,"completed progressive timeline accepted");
    s.synchronizeSections();juce::String error;
    for(int n=1;n<20;++n)expect(s.sectionDocument.split(uint64_t(n),s.sourceDocument.media.frames*n/20,error),"split explicit 3-second annotation");
    ui.refreshFromSession();ui.refreshTimeline();
    response.begin();click(ui.proposalPanel.propose);cancel(ui,"proposal_cancel_ms",r);response.end();
    expect(s.currentState.analysis==baseline&&!s.canApplySectionProposal(),"cancelled proposal changes no audio and cannot apply");
    response.begin();const auto proposal=Clock::now();click(ui.proposalPanel.propose);idle(ui);response.end();
    put(r,"section_propose_seconds",since(proposal));
    expect(s.sectionProposal&&s.sectionProposal->status()==SectionProposalStatus::ready
        &&s.sectionProposal->passes()==3,"fresh section proposal completes all three private passes");
    click(ui.proposalPanel.apply);
    expect(!s.currentState.hasMeasurement&&s.currentStateNeedsRemeasurement(),"explicit Apply invalidates Current");
    const auto plan=s.currentStateDocument.processingPlan();
    response.begin();click(ui.remeasureCurrentStateButton);cancel(ui,"current_cancel_ms",r);response.end();
    expect(!s.currentState.hasMeasurement&&s.currentStateNeedsRemeasurement()
        &&proposalPlanIdentity(s.currentStateDocument.processingPlan())==proposalPlanIdentity(plan),
        "cancelled verification leaves applied audio unverified with no rollback");
    response.begin();const auto current=Clock::now();click(ui.remeasureCurrentStateButton);idle(ui);response.end();
    put(r,"current_seconds",since(current));
    expect(s.currentState.hasMeasurement&&s.currentState.analysis->identicalSourceMeasurementIdentity.isEmpty()
        &&s.evaluation(Domain::currentState).verdict==proposal::Verdict::pass,"independent Current verification passes");
    expect(s.source.analysis==original&&!ui.audition.connected()&&!ui.audition.active(),"Source preserved and device closed throughout");
    put(r,"channels",s.sourceDocument.media.channels);put(r,"input_sha256",s.sourceDocument.media.sha256);
    put(r,"heartbeat_scope","Active native analysis/proposal jobs, including completion and cancellation; excludes picker/setup/screenshot work");
    snapshot(ui,dir.getChildFile("verified-current.png"));response.record(r);
}
}
int main(int argc,char** argv)
{
    if(argc!=3)return 2;
    const juce::File dir(argv[2]);if(dir.exists()||dir.createDirectory().failed())return 2;
    juce::ScopedJuceInitialiser_GUI gui;juce::var r(new juce::DynamicObject());bool passed=false;
    try
    {
        juce::DocumentWindow window("Intel true-peak qualification",juce::Colours::darkgrey,juce::DocumentWindow::allButtons);
        auto content=std::make_unique<Main>();auto* ui=content.get();
        window.setUsingNativeTitleBar(true);window.setContentOwned(content.release(),true);window.centreWithSize(1240,1000);
        window.setVisible(true);window.toFront(true);pump();
        expect(window.getPeer()&&ui->isShowing(),"live native main component");
        checkNative(*ui,juce::File(argv[1]),dir,r);passed=true;
    }
    catch(const std::exception& e){put(r,"failure",e.what());std::cerr<<e.what()<<'\n';}
    put(r,"passed",passed);put(r,"checks",checks);put(r,"build",LEVELSCOPE_MEASUREMENT_BUILD);
    put(r,"physical_devices_opened",false);put(r,"video_loaded",false);
    if(!dir.getChildFile("result.json").replaceWithText(juce::JSON::toString(r)))return 2;
    return passed?0:1;
}

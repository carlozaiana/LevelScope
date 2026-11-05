// Test launcher only: same native component/product objects, controlled file selection.
#define main m18_unused_native_main
#include "../../Profiling/ShortAcceptance/ShortAcceptance.cpp"
#undef main

namespace
{
struct Heartbeat:private juce::Timer
{
    std::vector<double> gaps;Clock::time_point previous;
    void begin(){gaps.clear();previous=Clock::now();startTimer(20);}
    void timerCallback()override{gaps.push_back(since(previous)*1000);previous=Clock::now();}
    void finish(juce::var& out)
    {
        stopTimer();expect(!gaps.empty(),"native passage view heartbeat observed");std::sort(gaps.begin(),gaps.end());
        const auto p95=gaps[size_t(std::ceil(.95*gaps.size()))-1];put(out,"passages_p95_ms",p95);put(out,"passages_max_ms",gaps.back());
        expect(p95<=100,"passage inspection p95 under 100 ms on 20 ms timer");
    }
};
LraPassagesPanel& panel(Main& ui)
{expect(bool(ui.lraPassagesWindow),"passage window exists");auto* p=dynamic_cast<LraPassagesPanel*>(ui.lraPassagesWindow->getContentComponent());expect(p!=nullptr,"actual LRA passage panel");return *p;}
void ready(Main& ui)
{until("passage view did not complete",[&]{auto& p=panel(ui);return !p.future.valid()&&p.result.passages&&p.key==p.selectedKey();});}
void run(Main& ui,const juce::File& file,const juce::File& dir,juce::var& result)
{
    auto& s=ui.session;ui.audioPanel.roleMemory=RememberedChannelRoles(dir.getChildFile("owned-roles"));
    ui.chooseSourceFile();choose(ui.sourceFileChooser.get(),file);
    until("role confirmation unavailable",[&]{return s.sourceDocument.hasSource&&ui.audioWindow&&ui.audioWindow->isVisible();});
    ui.jobRequirementsPanel.profiles.setSelectedId(2,juce::sendNotificationSync);ui.applicabilityConfirmation.setToggleState(true,juce::sendNotificationSync);
    click(ui.audioPanel.confirmButton);idle(ui);expect(s.source.hasMeasurement&&s.currentState.hasMeasurement,"one-pass Source and untouched Current ready");
    const auto source=s.source.analysis,current=s.currentState.analysis;const auto revision=s.revision;const auto frame=ui.audition.positionFrames();
    ui.showLraPassages(false);ready(ui);expect(!panel(ui).result.passages->ceiling&&!panel(ui).assess.isEnabled(),"without a job LRA ceiling passage view stays descriptive");ui.lraPassagesWindow->setVisible(false);
    const auto target=int(std::floor((*current->detail.lra.value-1)*10));
    ui.jobRequirementsPanel.fields[2].setText(juce::String(target/10.,1),true);pump();
    expect(ui.jobRequirementsPanel.commitFields(),"inline LRA target committed");ui.refreshFromSession();ui.refreshTimeline();
    expect(s.currentState.analysis==current&&s.jobRequirements.maximumLraTenths==target,"requirements reuse valid measurement");
    click(ui.proposalPanel.propose);ready(ui);Heartbeat heartbeat;heartbeat.begin();
    expect(panel(ui).domain==Domain::currentState&&panel(ui).assess.isEnabled(),"Propose opens explicit Current estimates before any measurement pass");
    expect(!s.job.isActive()&&!s.sectionProposal&&panel(ui).result.selection&&panel(ui).result.selection->available,"estimate view does not launch private assessment or mutate audio");
    picture(panel(ui).createComponentSnapshot(panel(ui).getLocalBounds()),dir.getChildFile("passages.png"));
    panel(ui).list.selectRow(0);pump();
    expect(ui.timeline.selection&&ui.audition.positionFrames()==frame&&s.revision==revision,"passage selection highlights exact interval without seeking or audio mutation");
    click(panel(ui).source);ready(ui);expect(!panel(ui).assess.isEnabled()&&panel(ui).status.getText().startsWith("SOURCE"),"Source view descriptive; Current correction requires Current inspection");
    panel(ui).quiet.setToggleState(true,juce::sendNotificationSync);panel(ui).list.selectRow(0);pump();
    expect(!ui.audition.connected()&&!ui.audition.active()&&s.source.analysis==source&&s.currentState.analysis==current,"viewing quiet/high Source populations preserves closed stopped audio");
    click(panel(ui).current);ready(ui);panel(ui).quiet.setToggleState(false,juce::sendNotificationSync);
    const auto inspection=Clock::now();while(since(inspection)<2){panel(ui).list.selectRow(0);pump();}heartbeat.finish(result);
    click(panel(ui).assess);until("private assessment did not become active",[&]{return s.job.isActive()&&ui.importControl->progress>.15f;});
    const auto cancelStart=Clock::now();click(ui.cancelJobButton);until("native cancellation exceeds one second",[&]{return !s.job.isActive();},1);
    expect(s.job.state==JobState::cancelled&&!s.sectionProposal&&s.currentState.analysis==current,"cancelled assessment leaves valid original Current and no Apply authority");put(result,"cancel_ms",since(cancelStart)*1000);
    click(ui.proposalPanel.propose);ready(ui);click(panel(ui).assess);idle(ui);expect(s.sectionProposal&&s.canApplySectionProposal(),"native measured candidate enables explicit Apply");
    const auto candidate=s.sectionProposal;ui.refreshTimeline();expect(ui.gainLane.proposed&&ui.gainLane.rows().size()>1,"marker-free candidate has navigable processing-region preview");
    click(ui.proposalPanel.review);expect(ui.proposalInspectionWindow&&ui.proposalInspectionWindow->isVisible(),"complete measured rationale is accessible");
    ui.proposalInspectionWindow->setVisible(false);click(ui.proposalPanel.apply);
    expect(s.currentStateNeedsRemeasurement()&&!s.currentState.hasMeasurement&&s.lastAppliedExplanation,"explicit native Apply preserves historical rationale and marks Current unverified");
    click(ui.remeasureCurrentStateButton);idle(ui);expect(s.evaluation(Domain::currentState).verdict==proposal::Verdict::pass&&s.currentState.analysis->detail.lra.value==candidate->assessment()->detail.lra.value,"independent native Current remeasurement verifies exact candidate");
    click(ui.proposalPanel.history);expect(ui.proposalInspectionWindow->isVisible(),"historical explanation remains accessible after measurement");ui.proposalInspectionWindow->setVisible(false);
    snapshot(ui,dir.getChildFile("verified-current.png"));
    expect(s.source.analysis==source&&s.sectionDocument.map().sections.size()==1&&!ui.audition.connected()&&!ui.audition.active(),"original Source, whole-programme annotation and closed device preserved throughout");
    put(result,"channels",s.sourceDocument.media.channels);put(result,"input_sha256",s.sourceDocument.media.sha256);put(result,"private_passes",candidate->passes());put(result,"target_lra",target/10.);
}
}
int main(int argc,char** argv)
{
    if(argc!=3)return 2;juce::File dir(argv[2]);if(dir.exists()||dir.createDirectory().failed())return 2;
    juce::ScopedJuceInitialiser_GUI gui;juce::var result(new juce::DynamicObject());bool passed=false;
    try
    {
        juce::DocumentWindow window("M19 native silent qualification",juce::Colours::darkgrey,juce::DocumentWindow::allButtons);
        auto content=std::make_unique<Main>();auto* ui=content.get();window.setUsingNativeTitleBar(true);window.setContentOwned(content.release(),true);window.centreWithSize(1240,1000);window.setVisible(true);window.toFront(true);pump();
        expect(window.getPeer()&&ui->isShowing(),"live native M19 main component");run(*ui,juce::File(argv[1]),dir,result);passed=true;
    }
    catch(const std::exception& e){put(result,"failure",e.what());std::cerr<<e.what()<<'\n';}
    put(result,"passed",passed);put(result,"checks",checks);put(result,"build",LEVELSCOPE_MEASUREMENT_BUILD);put(result,"physical_devices_opened",false);put(result,"video_loaded",false);
    if(!dir.getChildFile("result.json").replaceWithText(juce::JSON::toString(result)))return 2;return passed?0:1;
}

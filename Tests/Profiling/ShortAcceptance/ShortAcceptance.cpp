#include "Standalone/StandaloneMainComponent.h"
#include "Standalone/StandaloneProcessing.h"
#include "AuditionTestSupport.h"
#include <algorithm>
#import <AppKit/AppKit.h>

namespace
{
using Clock = std::chrono::steady_clock;
using Domain = MeasurementDomain;
using Main = StandaloneMainComponent;

double since(Clock::time_point start) { return std::chrono::duration<double>(Clock::now() - start).count(); }
void put(juce::var& row, const char* key, juce::var value) { row.getDynamicObject()->setProperty(key, value); }
void pump()
{
    @autoreleasepool
    {
        NSEvent* event=[NSApp nextEventMatchingMask:NSEventMaskAny untilDate:[NSDate dateWithTimeIntervalSinceNow:0.002]
            inMode:NSDefaultRunLoopMode dequeue:YES];
        if(event) [NSApp sendEvent:event];
        [NSApp updateWindows];
    }
}
template<class Predicate> void until(const char* reason, Predicate done)
{
    const auto start = Clock::now();
    while (!done())
    {
        if (since(start) > 30) throw std::runtime_error(reason);
        pump();
    }
}
void click(juce::Button& button)
{
    if (!button.isEnabled()) throw std::runtime_error(("Disabled: " + button.getButtonText() + " " + button.getTooltip()).toStdString());
    button.onClick();
}
void choose(juce::FileChooser* picker, const juce::File& path)
{
    expect(picker && picker->asyncCallback, "file picker opened with an active completion");
    picker->finished({juce::URL(path)});
}
juce::MouseEvent mouse(juce::Component& c, juce::Point<float> p, bool down)
{
    const auto now=juce::Time::getCurrentTime();
    return {juce::Desktop::getInstance().getMainMouseSource(),p,juce::ModifierKeys(down ? juce::ModifierKeys::leftButtonModifier : 0),
            0,0,0,0,0,&c,&c,now,p,now,1,false};
}
float timeX(MeasurementTimeline& t,double seconds)
{
    const auto range=t.viewRange();
    return float(t.graphBounds().getX()+(seconds*48000-range.first)*(t.graphBounds().getWidth()-1)/(range.second-range.first));
}
void drag(MeasurementTimeline& t, double from, double to)
{
    const float y=float(t.graphBounds().getCentreY());
    t.mouseDown(mouse(t,{timeX(t,from),y},true));
    t.mouseDrag(mouse(t,{timeX(t,to),y},true));
    t.mouseUp(mouse(t,{timeX(t,to),y},false));
}
void picture(const juce::Image& image,const juce::File& file)
{
    auto stream=file.createOutputStream();
    expect(stream && image.isValid() && juce::PNGImageFormat().writeImageToStream(image,*stream), "write fresh UI capture");
}
void snapshot(Main& ui,const juce::File& file) {picture(ui.createComponentSnapshot(ui.getLocalBounds()),file);}
void idle(Main& ui)
{
    until("workflow job timeout",[&]{return !ui.session.job.isActive() && !ui.session.persistenceLoading;});
    if (ui.session.job.state==JobState::failed || ui.session.job.state==JobState::cancelled)
        throw std::runtime_error(ui.session.job.message.toStdString());
}
void selectOutput(Main& ui, Device& device, int channels)
{
    if(ui.audioPanel.deviceBox.getNumItems()==0) ui.audioPanel.deviceBox.addItem("Acceptance output",1);
    ui.audioPanel.deviceBox.setSelectedId(1,juce::dontSendNotification);
    ui.audioPanel.rebuildOutputs();
    ui.audioPanel.refresh();
    click(ui.audioPanel.openButton);
    expect(ui.audition.connected() && device.isOpen(),"injected output opened at programme rate");
    for(int i=0;i<channels;++i) ui.audioPanel.outputBoxes[size_t(i)].setSelectedId(i+2,juce::sendNotificationSync);
    for(int i=0;i<channels;++i) expect(ui.session.sourceDocument.routing.fileToDevice[size_t(i)]==i,"explicit one-to-one monitoring route");
}
void auditionFirst(Main& ui,Device& device,int channels,Clock::time_point import,juce::var& row)
{
    juce::AudioBuffer<float> block(channels,256);
    const auto play=Clock::now();
    click(ui.timelinePlay);
    expect(ui.audition.active(),"first Play started transport preparation");
    until("device did not start",[&]{return device.isPlaying();});
    block.clear();
    const auto callback=Clock::now();
    device.tick(block);
    int first=0;
    for(;first<block.getNumSamples();++first)
    {
        bool nonzero=false;
        for(int c=0;c<channels;++c) nonzero |= block.getSample(c,first)!=0.0f;
        if(nonzero) break;
    }
    expect(first==1,"first non-zero post-monitor PCM occurs at programme frame one");
    for(int c=0;c<channels;++c) expect(std::isfinite(block.getSample(c,1)) && block.getSample(c,1)!=0,"first audible frame present on each explicitly routed output");
    const double withinBlock=first/48000.;
    put(row,"play_to_first_output_s",std::chrono::duration<double>(callback-play).count()+withinBlock);
    put(row,"import_to_first_output_s",std::chrono::duration<double>(callback-import).count()+withinBlock);
    put(row,"first_nonzero_frame",first);
    ui.timeline.focusForTransport();
    expect(ui.getPeer()->handleKeyPress(juce::KeyPress(juce::KeyPress::spaceKey)) && !ui.audition.active(),"Space stops native timeline playback");
    click(ui.audioPanel.closeButton);
    expect(!ui.audition.connected(),"transport stopped and output disconnected");
}
void auditionCurrent(Main& ui,Device& device,int channels)
{
    auto& s=ui.session;
    selectOutput(ui,device,channels);
    ui.auditionDomainBox.setSelectedId(2,juce::sendNotificationSync);
    const auto boundary=s.sectionDocument.map().sections[1].begin;
    ui.timeline.setViewRange(double(boundary-4800),double(boundary+4800));
    auto point=ui.timeline.rulerBounds().getCentre().toFloat();
    ui.timeline.mouseDown(mouse(ui.timeline,point,true));
    ui.timeline.mouseUp(mouse(ui.timeline,point,false));
    const auto frame=ui.audition.positionFrames();
    expect(s.auditionDomain==Domain::currentState && std::abs(frame-boundary)<32 && !ui.audition.active(),
           "stopped main timecode seek lands at the edited boundary");
    Capture capture(channels,256,frame);
    expect(ui.audition.setCaptureObserver(&capture),"attach a bounded pre-monitor capture");
    const auto measured=s.currentState.analysis;
    click(ui.timelinePlay);
    until("Current device did not start",[&]{return device.isPlaying();});
    juce::AudioBuffer<float> block(channels,256);
    block.clear(); device.tick(block);click(ui.timelineStop);
    expect(ui.audition.setCaptureObserver(nullptr) && capture.valid && capture.next==frame+256,"Current audition delivered all selected frames");
    auto reader=ProgrammeReader::open({s.sourceDocument.file},s.sourceDocument.media);
    expect(reader && reader->read(block.getArrayOfWritePointers(),channels,frame,256),"read selected programme passage independently");
    ProgrammeGainProcessor processing(s.currentStateDocument.processingPlan());
    expect(processing.process(block.getArrayOfWritePointers(),channels,256,frame),"process the exact applied state for comparison");
    for(int c=0;c<channels;++c)
        expect(std::memcmp(block.getReadPointer(c),capture.data.getReadPointer(c),256*sizeof(float))==0,"Current audition matches the applied samples on the routed channel");
    expect(measured==s.currentState.analysis,"audition leaves independent measurement intact");
    click(ui.audioPanel.closeButton);
}
void run(Main& ui,const juce::File& input,const juce::File& dir,int channels,juce::var& result)
{
    auto& s=ui.session;
    ui.audioPanel.roleMemory=RememberedChannelRoles(dir.getChildFile("role-preferences"));
    auto device=std::make_unique<Device>(); device->outputCount=channels;
    auto* output=device.get(); ui.audition.device=std::move(device);
    ui.audioPanel.rebuildOutputs();
    ui.chooseSourceFile();
    const auto imported=Clock::now(); choose(ui.sourceFileChooser.get(),input);
    until("import did not present roles",[&]{return s.sourceDocument.hasSource && ui.audioWindow && ui.audioWindow->isVisible();});
    put(result,"import_to_role_review_s",since(imported));
    expect(!s.sourceDocument.layout.confirmed && !s.source.analysis,"review precedes measurement");
    expect(ui.audioPanel.confirmButton.getButtonText().containsIgnoreCase("analyse"),"one-click confirm and analyse action");
    // Programme applicability is explicitly confirmed for each newly imported job.
    ui.jobRequirementsPanel.profiles.setSelectedId(2,juce::sendNotificationSync);
    ui.applicabilityConfirmation.setToggleState(true,juce::sendNotificationSync);
    expect(s.selectedTargetProfileIndex==1 && s.programmeApplicability==proposal::Applicability::preProducedLinearProgramme,"EBU R128 programme context");
    selectOutput(ui,*output,channels);
    const auto confirmed=Clock::now(); click(ui.audioPanel.confirmButton);
    if(!s.job.isActive()) throw std::runtime_error(("Confirm roles did not start analysis: "+s.job.message+" / "+ui.audioPanel.actionError).toStdString());
    until("Play remained enabled during measurement",[&]{return !ui.timelinePlay.isEnabled();});
    expect(s.job.isActive(),"analysis starts and playback waits for the active job");
    juce::Image preview;
    until("no incomplete waveform painted",[&]{
        const auto& cols=ui.timeline.displayedWaveform();
        if(!ui.timeline.waveformAccepted && std::any_of(cols.begin(),cols.end(),[](const auto& c){return c.available();}))
        {
            preview=ui.timeline.createComponentSnapshot(ui.timeline.getLocalBounds());
            return preview.isValid();
        }
        if(s.source.hasMeasurement) throw std::runtime_error("Source completed before a partial waveform reached the visible timeline");
        return false;
    });
    put(result,"import_to_first_waveform_s",since(imported));
    put(result,"confirm_to_first_waveform_s",since(confirmed));
    expect(!s.currentState.hasMeasurement && !ui.timeline.waveformAccepted,"painted preview is incomplete and has no Current verdict");
    until("no discernible partial waveform reached the native timeline",[&]{
        const auto& cols=ui.timeline.displayedWaveform();
        const auto painted=std::count_if(cols.begin(),cols.end(),[](const auto& c){
            return c.available() && c.channels[0].maximum>c.channels[0].minimum;
        });
        return painted>=16;
    });
    const auto discernible=ui.timeline.createComponentSnapshot(ui.timeline.getLocalBounds());
    put(result,"import_to_16_signal_columns_s",since(imported));
    until("Source did not initialize Current",[&]{return s.currentState.hasMeasurement && !s.job.isActive();});
    put(result,"import_to_measurements_ready_s",since(imported));
    expect(s.source.hasMeasurement && s.currentState.analysis->identicalSourceMeasurementIdentity==s.source.analysis->identity
           && ui.timeline.waveformAccepted,"one accepted Source measurement populated untouched Current with exact origin");
    expect(s.sourceDocument.media.channels==channels && s.sourceDocument.media.frames==2880000 && s.sourceDocument.media.sampleRate==48000,"full short fixture metadata");
    put(result,"input_sha256",s.sourceDocument.media.sha256); put(result,"frames",s.sourceDocument.media.frames);
    put(result,"source_il_lufs",*s.source.integratedLufs); put(result,"source_true_peak_dbtp",*s.source.truePeakDbtp);
    auditionFirst(ui,*output,channels,imported,result);
    picture(preview,dir.getChildFile("first-waveform.png"));
    picture(discernible,dir.getChildFile("progressive-waveform.png"));
    const auto original=s.source.analysis,baseline=s.currentState.analysis;
    ui.inspectDomain.setSelectedId(2,juce::sendNotificationSync);
    expect(s.auditionDomain==Domain::source,"inspect Current leaves audition on Original");
    ui.inspectDomain.setSelectedId(1,juce::sendNotificationSync);
    ui.timeline.setViewRange(0,2880000);
    drag(ui.timeline,30,34);
    expect(s.sectionDocument.map().sections.size()==3,"one waveform gesture creates a local section");
    click(ui.sectionsButton); pump();
    auto& second=ui.sectionsPanel->timelineView();
    second.setViewRange(0,2880000); second.refreshSectionEditing();
    auto boundary=s.sectionDocument.map().sections[1].end;
    drag(second,double(boundary)/48000.,33);
    expect(s.sectionDocument.map().sections[1].end<boundary && s.source.analysis==original && s.currentState.analysis==baseline,"boundary adjusted in Sections without disturbing either measurement");
    auto ruler=second.rulerBounds().getCentre().toFloat();
    second.mouseDown(mouse(second,ruler,true)); second.mouseUp(mouse(second,ruler,false));
    expect(ui.audition.positionFrames()>0 && !ui.audition.active(),"stopped ruler seek in Sections");
    ui.sectionsWindow->setVisible(false);
    ui.timeline.setViewRange(28*48000,36*48000);
    const auto position=ui.audition.positionFrames();
    ui.gainLane.selectFrame(31*48000,true);
    expect(ui.audition.positionFrames()==position,"gain selection and zoom leave audition position unchanged");
    ui.gainLane.gain.setText("-0.75",true); click(ui.gainLane.stage);
    expect(s.sectionTrimDraft && s.currentState.analysis==baseline,"numeric gain remains an unapplied draft");
    click(ui.gainLane.discard);
    expect(!s.sectionTrimDraft && s.currentState.analysis==baseline,"Discard draft restores applied gains");
    const auto propose=Clock::now(); click(ui.proposalPanel.propose); idle(ui);
    put(result,"proposal_s",since(propose));
    put(result,"proposal_outcome_text",ui.proposalPanel.outcomeText());
    expect(s.sectionProposal && s.canApplySectionProposal() && s.currentState.analysis==baseline,"Propose offers an unapplied correction");
    const auto& proposed=*s.sectionProposal;
    const bool section=channels==6;
    expect(proposed.status()==(section?SectionProposalStatus::ready:SectionProposalStatus::globalOnly),"expected global or section-led outcome");
    expect(section ? (proposed.passes()==3 && proposed.assessment() && proposed.selection().corrections.size()==1) : proposed.passes()==0,"candidate evidence matches proposal method");
    put(result,"proposal_passes",proposed.passes());put(result,"proposal_outcome",proposed.outcomeTitle());
    const auto previewPlan=proposed.previewPlan();
    expect(previewPlan.has_value(),"proposal exposes the complete preview plan");
    put(result,"proposed_global_db",previewPlan->staticGainDb);
    if(section)
    {
        const auto& change=proposed.selection().corrections.front();
        put(result,"local_gain_db",change.resultingDb);put(result,"local_begin",juce::int64(change.begin));put(result,"local_end",juce::int64(change.end));
    }
    put(result,"proposed_local_count",int(proposed.selection().corrections.size()));
    click(ui.proposalPanel.review);
    expect(ui.proposalInspectionWindow && ui.proposalInspectionWindow->isShowing(),"complete proposal review opens from main window");
    ui.proposalInspectionWindow->setVisible(false);
    snapshot(ui,dir.getChildFile("proposal-ready.png"));
    click(ui.proposalPanel.apply);
    expect(s.currentStateNeedsRemeasurement() && !s.currentState.hasMeasurement && !ui.renderPcmButton.isEnabled()
           && s.lastAppliedExplanation,"explicit Apply makes Current unverified and blocks rendering");
    expect(proposalPlanIdentity(s.currentStateDocument.processingPlan())==proposalPlanIdentity(*previewPlan),"Apply installs exactly the reviewed global and local gains");
    auto explanation=s.lastAppliedExplanation->text();
    const auto remeasure=Clock::now();click(ui.remeasureCurrentStateButton);idle(ui);
    put(result,"current_remeasurement_s",since(remeasure));
    expect(s.currentState.analysis && s.currentState.analysis->domain==Domain::currentState && s.currentState.analysis->identicalSourceMeasurementIdentity.isEmpty()
           && s.evaluation(Domain::currentState).verdict==proposal::Verdict::pass,"independently measured applied Current passes");
    put(result,"current_il_lufs",*s.currentState.integratedLufs);put(result,"current_true_peak_dbtp",*s.currentState.truePeakDbtp);
    auditionCurrent(ui,*output,channels);
    const auto currentId=s.currentState.analysis->identity;
    const auto map=s.sectionDocument.map().digest();
    const auto plan=proposalPlanIdentity(s.currentStateDocument.processingPlan());
    const auto layout=s.sourceDocument.layout;
    auto file=dir.getChildFile("job.levelscope");
    click(ui.saveSessionButton); choose(ui.sessionChooser.get(),file);
    until("save did not finish",[&]{return !ui.documentResult.valid();});
    expect(file.existsAsFile() && s.documentFile==file && !s.documentDirty(),"session saved from main window");
    const auto oldSession=s.sessionId;
    click(ui.openSessionButton);choose(ui.sessionChooser.get(),file);
    until("session reopen did not finish",[&]{return !ui.documentResult.valid();});
    expect(s.sessionId!=oldSession && s.currentState.analysis && s.currentState.analysis->restoredFromSession
           && s.currentState.analysis->identity==currentId && s.evaluation(Domain::currentState).verdict==proposal::Verdict::pass,"reopened and revalidated independently measured Current");
    expect(s.lastAppliedExplanation && s.lastAppliedExplanation->text()==explanation && s.sectionDocument.map().digest()==map
           && proposalPlanIdentity(s.currentStateDocument.processingPlan())==plan,"reopened complete applied explanation, annotations and identical processing");
    expect(layout.roles==s.sourceDocument.layout.roles && s.sourceDocument.layout.confirmed && !ui.audition.connected()
           && !ui.audition.active() && std::all_of(s.sourceDocument.routing.fileToDevice.begin(),s.sourceDocument.routing.fileToDevice.end(),[](int r){return r==-1;}),"reopen restores project roles and leaves playback stopped and routes muted");
    click(ui.proposalPanel.history);expect(ui.proposalInspectionWindow->isShowing(),"saved explanation inspectable after reopen");ui.proposalInspectionWindow->setVisible(false);
    const auto dest=dir.getChildFile("export");
    const auto render=Clock::now();click(ui.renderPcmButton);choose(ui.exportChooser.get(),dest);idle(ui);
    put(result,"render_and_measure_artifact_s",since(render));
    expect(ui.exportCandidate && ui.publishPcmButton.isEnabled() && !dest.exists(),"independently measured artifact ready for explicit Publish");
    const auto request=exporting::capture(s);
    expect(exporting::artifactEvaluation(*ui.exportCandidate,request).verdict==proposal::Verdict::pass,"encoded artifact independently passes");
    put(result,"artifact_sha256",ui.exportCandidate->hash());
    put(result,"artifact_il_lufs",*ui.exportCandidate->analysis().detail.integrated.value);
    put(result,"artifact_true_peak_dbtp",*ui.exportCandidate->analysis().detail.truePeak.value);
    click(ui.publishPcmButton);idle(ui);
    expect(dest.getChildFile("programme.wav").existsAsFile() && dest.getChildFile("manifest.json").existsAsFile()
           && dest.getChildFile("report.html").existsAsFile() && dest.getChildFile("report.json").existsAsFile(),"published PCM24 audio and complete report bundle");
    snapshot(ui,dir.getChildFile("reopened-exported.png"));
}
}

int main(int argc,char** argv)
{
    if(argc!=4) return 2;
    const juce::File dir(argv[2]); if(dir.exists() || dir.createDirectory().failed())return 2;
    juce::ScopedJuceInitialiser_GUI gui;
    juce::var result(new juce::DynamicObject());
    bool passed=false;
    try
    {
        juce::DocumentWindow window("M18 short acceptance",juce::Colours::darkgrey,juce::DocumentWindow::allButtons);
        auto content=std::make_unique<Main>();auto* main=content.get();
        window.setUsingNativeTitleBar(true);window.setContentOwned(content.release(),true);window.centreWithSize(1240,1000);
        window.setVisible(true);window.toFront(true);pump();
        expect(window.getPeer() && main->isShowing(),"main window has a live native peer");
        run(*main,juce::File(argv[1]),dir,std::stoi(argv[3]),result);
        passed=true;
    }
    catch(const std::exception& e) {put(result,"failure",e.what());std::cerr<<e.what()<<'\n';}
    put(result,"passed",passed);put(result,"checks",checks);put(result,"measurement_build",LEVELSCOPE_MEASUREMENT_BUILD);
    if(!dir.getChildFile("result.json").replaceWithText(juce::JSON::toString(result)))return 2;
    return passed?0:1;
}

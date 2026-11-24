#include "CompositionTestFixtures.h"
#include "Standalone/StandaloneMainComponent.h"
#include "Standalone/TimelineTransport.h"
#include "Standalone/RegionReviewFormatting.h"
#include "Standalone/SessionSavePresentation.h"
#include "Standalone/ProposalPresentation.h"
#include "Core/Measurement/GatePopulation.h"

namespace levelscope::standalone
{
struct PreparationFeedbackTestAccess
{
    static void enableEventAlternatives(StandaloneMainComponent& c,bool enabled){c.eventAlternativeSelectionQualified=enabled;}
    static bool eventAlternativesEnabled(const StandaloneMainComponent& c){return c.eventAlternativeSelectionQualified;}
    static void review(StandaloneMainComponent& c,bool historical=false){c.reviewProposal(historical);}
    static RegionProposalReviewPanel* reviewPanel(StandaloneMainComponent& c){return c.proposalInspectionWindow?dynamic_cast<RegionProposalReviewPanel*>(c.proposalInspectionWindow->getContentComponent()):nullptr;}
    static juce::String historyStatus(StandaloneMainComponent& c){return c.session.appliedExplanationStatus();}
    static bool proposalPending(StandaloneMainComponent& c){return c.sectionProposalResult.valid();}
    static void closeReview(StandaloneMainComponent& c){c.proposalInspectionWindow.reset();}
    static void settings(StandaloneMainComponent& c){c.showRegionPlanning();}
    static RegionPlanningPanel* settingsPanel(StandaloneMainComponent& c){return c.regionPlanningWindow?dynamic_cast<RegionPlanningPanel*>(c.regionPlanningWindow->getContentComponent()):nullptr;}
    static void closeSettings(StandaloneMainComponent& c){c.regionPlanningWindow.reset();}
    static void save(StandaloneMainComponent& c,bool copy=false){c.beginSaveSession(copy);}
    static void open(StandaloneMainComponent& c,const juce::File& f){c.beginOpenSession(f);}
    static bool filePending(StandaloneMainComponent& c){return c.documentResult.valid();}
    static void deterministicChooser(StandaloneMainComponent& c){c.useNativeSessionFileChooser=false;}
    static juce::String fileNotice(StandaloneMainComponent& c){return c.documentNotice;}
    static StandaloneSessionModel& model(StandaloneMainComponent& c){return c.session;}
    static void stopTimer(StandaloneMainComponent& c){c.stopTimer();}
    static void tick(StandaloneMainComponent& c){c.timerCallback();}
    static void refresh(StandaloneMainComponent& c){c.refreshFromSession();c.refreshTimeline();c.resized();}
    // Historical explicit-Apply checks exercise the retained legacy service route.
    static void propose(StandaloneMainComponent& c){c.beginSectionProposal();}
    static void workflow(StandaloneMainComponent& c){c.proposeWithActiveMethod();}
    static ProposalOperation& operation(StandaloneMainComponent& c){return c.proposalOperation;}
    static bool workers(StandaloneMainComponent& c){return c.analysisResult.valid()||c.sectionProposalResult.valid();}
    static void cancel(StandaloneMainComponent& c){c.cancelJobButton.onClick();}
    static TimelineGainLane& lane(StandaloneMainComponent& c){return c.gainLane;}
    static ProposalWorkflowPanel& proposal(StandaloneMainComponent& c){return c.proposalPanel;}
    static bool stopped(StandaloneMainComponent& c){return !c.audition.active()&&!c.audition.connected();}
    static MeasurementTimeline& timeline(StandaloneMainComponent& c){return c.timeline;}
    static void moveDivider(StandaloneMainComponent& c,int delta){c.waveformDivider.move(delta);}
    static juce::Viewport& inspector(StandaloneMainComponent& c){return c.inspectorViewport;}
    static void inspectCurrent(StandaloneMainComponent& c){c.inspectDomain.setSelectedId(2,juce::sendNotificationSync);}
    static int inspection(StandaloneMainComponent& c){return c.inspectDomain.getSelectedId();}
    static StandaloneAudioPanel& setup(StandaloneMainComponent& c){return c.audioPanel;}
    static StandaloneAudition& audio(StandaloneMainComponent& c){return c.audition;}
    static void refreshTimeline(StandaloneMainComponent& c){c.refreshTimeline();}
    static juce::Label& warning(StandaloneMainComponent& c){return c.monitorWarning;}
    static juce::Label& scope(StandaloneMainComponent& c){return c.qcScope;}
    static juce::Slider& globalGain(StandaloneMainComponent& c){return c.currentGain;}
    static void measureSource(StandaloneMainComponent& c){c.beginMeasurement(MeasurementDomain::source);}
    static void play(StandaloneMainComponent& c){c.playTimeline();}
    static void clearSource(StandaloneMainComponent& c){c.clearSourceFile();}
    static juce::String stage(StandaloneMainComponent& c){return c.jobLabel.getText();}
    static juce::String detail(StandaloneMainComponent& c){return c.jobDetail.getText();}
    static bool cancellationAvailable(StandaloneMainComponent& c){return c.cancelJobButton.isEnabled();}
    static bool playAvailable(StandaloneMainComponent& c){return c.timelinePlay.isEnabled();}
    static bool overviewComplete(StandaloneMainComponent& c){return c.waveformCache.view(WaveformBinding::from(c.session.measurementContext(MeasurementDomain::source),MeasurementDomain::source)).accepted;}
};
}
namespace
{
using namespace levelscope::proposal;
using Access=PreparationFeedbackTestAccess;
struct Heartbeat:private juce::Timer
{
    std::vector<double> gaps;double previous=0;
    void begin(){gaps.clear();previous=juce::Time::getMillisecondCounterHiRes();startTimer(20);}
    void timerCallback()override{const auto now=juce::Time::getMillisecondCounterHiRes();gaps.push_back(now-previous);previous=now;}
    void finish()
    {
        stopTimer();std::sort(gaps.begin(),gaps.end());
        const auto p95=gaps.empty()?std::numeric_limits<double>::infinity():gaps[size_t(std::ceil(.95*gaps.size()))-1];auto v=juce::var(new juce::DynamicObject());
        v.getDynamicObject()->setProperty("native_heartbeat_p95_ms",p95);v.getDynamicObject()->setProperty("native_heartbeat_max_ms",gaps.empty()?0.:gaps.back());v.getDynamicObject()->setProperty("samples",int(gaps.size()));cases.add(v);
        expect(gaps.size()>20,"native M20 heartbeat collected meaningful samples");
        expect(p95<=100,"native M20 proposal/review heartbeat p95 <=100 ms");
    }
};
void inspectionPollingChecks()
{
    int standaloneCalls=0,ownedCalls=0;
    ProposalInspectionPanel standalone("Standalone report",[&]{++standaloneCalls;return "Live standalone status";});
    ProposalInspectionPanel owned("Embedded report",[&]{++ownedCalls;return "Parent-refreshed status";},1500,false);
    standalone.addToDesktop(juce::ComponentPeer::windowHasTitleBar);owned.addToDesktop(juce::ComponentPeer::windowHasTitleBar);
    standalone.setVisible(false);owned.setVisible(true);
    const auto hidden=standaloneCalls,initialOwned=ownedCalls;
    juce::MessageManager::getInstance()->runDispatchLoopUntil(650);
    expect(standaloneCalls==hidden,"hidden standalone inspection does not poll status");
    expect(ownedCalls==initialOwned,"embedded inspection has no duplicate status timer");
    owned.refreshStatus();expect(ownedCalls==initialOwned+1,"parent can synchronously refresh visible embedded status");
    standalone.setVisible(true);juce::MessageManager::getInstance()->runDispatchLoopUntil(650);
    expect(standaloneCalls>hidden,"visible standalone inspection retains live status polling");
    standalone.setVisible(false);const auto closed=standaloneCalls;
    juce::MessageManager::getInstance()->runDispatchLoopUntil(650);
    expect(standaloneCalls==closed,"closing inspection stops background status work");
}
template<class T> T* find(juce::Component& root,const juce::String& name)
{
    if(auto* c=dynamic_cast<T*>(&root))
    {if(c->getName()==name)return c;if constexpr(std::is_base_of_v<juce::Button,T>)if(c->getButtonText()==name)return c;}
    for(auto* child:root.getChildren())if(auto* c=find<T>(*child,name))return c;return nullptr;
}
template<class T> T& control(juce::Component& root,const juce::String& name)
{auto* c=find<T>(root,name);if(!c)throw std::runtime_error(("Missing UI control: "+name).toStdString());return *c;}
void click(juce::Component& root,const juce::String& name)
{auto& b=control<juce::TextButton>(root,name);expect(b.isEnabled(),"requested UI action is enabled");b.onClick();juce::MessageManager::getInstance()->runDispatchLoopUntil(2);}
void text(juce::Component& root,const juce::String& name,const juce::String& value)
{control<juce::TextEditor>(root,name).setText(value,true);}
void choose(juce::Component& root,const juce::String& name,int id)
{control<juce::ComboBox>(root,name).setSelectedId(id,juce::sendNotificationSync);}
void toggle(juce::Component& root,const juce::String& name,bool value)
{auto& b=control<juce::ToggleButton>(root,name);b.setToggleState(value,juce::dontSendNotification);if(b.onClick)b.onClick();}
void imageOf(juce::Component& c,const juce::File& file)
{auto image=c.createComponentSnapshot(c.getLocalBounds());auto stream=file.createOutputStream();juce::PNGImageFormat format;expect(stream&&format.writeImageToStream(image,*stream),"native UI image written");}
void waitJob(StandaloneMainComponent& c)
{
    const auto until=juce::Time::getMillisecondCounterHiRes()+90000;
    do{Access::tick(c);juce::MessageManager::getInstance()->runDispatchLoopUntil(5);}while(Access::model(c).job.isActive()&&juce::Time::getMillisecondCounterHiRes()<until);
    expect(!Access::model(c).job.isActive(),"main worker finishes within focused timeout");Access::refresh(c);
}
StandaloneSessionModel programme(const juce::File& root,int channels,bool heterogeneous=false)
{
    expect(root.createDirectory().wasOk(),"owned test root");const int rate=48000,seconds=65;
    const auto file=root.getChildFile("programme.wav");expect(!file.existsAsFile(),"never overwrite original or earlier fixture");
    auto stream=file.createOutputStream();juce::WavAudioFormat format;std::unique_ptr<juce::AudioFormatWriter> writer(format.createWriterFor(stream.get(),rate,unsigned(channels),24,{},0));expect(bool(writer),"PCM24 fixture writer");stream.release();
    juce::AudioBuffer<float> block(channels,4096);bool written=true;
    for(int first=0;first<rate*seconds;first+=4096)
    {
        const int count=std::min(4096,rate*seconds-first);
        for(int ch=0;ch<channels;++ch)for(int n=0;n<count;++n)
        {const int f=first+n;const bool high=heterogeneous?((f>=rate*15&&f<rate*23)||(f>=rate*35&&f<rate*43)):(f>=rate*30&&f<rate*36);const double amplitude=channels==6&&ch==3?.12:.03*(1+.07*ch);block.setSample(ch,n,float(amplitude*(high?2.:1.)*std::sin(2*juce::MathConstants<double>::pi*997*f/rate)));}
        written&=writer->writeFromAudioSampleBuffer(block,0,count);
    }
    expect(written,"complete fixture");writer.reset();JobControl job;const auto media=validateMedia(file,job);expect(media.succeeded(),"fixture validated");StandaloneSessionModel s;expect(s.setValidatedSource(file,media.media),"fresh programme");s.sourceDocument.layout=layoutFor(channels);s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(Applicability::preProducedLinearProgramme);measuredComposition(s,MeasurementDomain::source);if(!s.currentState.hasMeasurement){s.initializeCurrentStateFromSource();measuredComposition(s,MeasurementDomain::currentState);}s.synchronizeSections();return s;
}

void editorChecks(StandaloneSessionModel& s,const juce::File& root)
{
    RegionPlanningPanel panel(s);const auto revision=s.revision;const auto snapshot=s.currentState.analysis->identity;
    expect(control<juce::ComboBox>(panel,"Proposal method").getSelectedId()==1&&!s.regionPlanning,"M19 remains the default; opening editor is not experimental opt-in");
    expect(RegionPlanningPanel::exactFrame("9223372036854775807")==std::numeric_limits<int64_t>::max()&&!RegionPlanningPanel::exactFrame("9223372036854775808")&&!RegionPlanningPanel::exactFrame("4.5")&&!RegionPlanningPanel::exactFrame("-1"),"frame parser preserves full int64 precision and rejects overflow/fractions");
    choose(panel,"Proposal method",2);choose(panel,"Region mode",3);text(panel,"Temporary begin frame","576001");text(panel,"Temporary end frame","2207999");
    expect(!s.regionPlanning&&s.revision==revision,"method and numeric edits remain private before confirmation");
    std::optional<RegionPlanningSettings> draft;expect(panel.readDraft(draft).isEmpty()&&draft&&draft->selection.begin==576001&&!draft->allowGlobal,"selected mode defaults to global fixed and retains exact off-grid frame");
    toggle(panel,"Restrict processing start",true);toggle(panel,"Pin start",true);text(panel,"Processing start minimum frame","576001");
    text(panel,"Left margin seconds (0-3, step 0.1)","0.5");text(panel,"Right margin seconds (0-3, step 0.1)","1.5");
    expect(panel.readDraft(draft).isEmpty()&&draft->start->minimum==576001&&draft->start->maximum==576001&&draft->leftMarginMs==500&&draft->rightMarginMs==1500,"pins and decimal margins translate exactly without rounding");
    text(panel,"Need-change threshold LU (0.5-3)","1.01");expect(panel.readDraft(draft).isNotEmpty(),"UI rejects off-step discovery values rather than clamping");text(panel,"Need-change threshold LU (0.5-3)","1.0");
    auto& cuts=control<RegionFrameRows>(panel,"Internal nominal cuts - source frames; ranges or exact pins");
    text(cuts,cuts.getName()+" minimum / begin","1200001");toggle(cuts,"Exact pin",true);expect(panel.readDraft(draft).isNotEmpty(),"unfinished exact-frame row cannot be silently omitted");click(cuts,"Add row");
    auto& protectedRows=control<RegionFrameRows>(panel,"No additional local-envelope change - [begin, end) source frames");
    text(protectedRows,protectedRows.getName()+" minimum / begin","1400001");text(protectedRows,protectedRows.getName()+" maximum / end","1500007");click(protectedRows,"Add row");
    expect(panel.readDraft(draft).isEmpty()&&draft->internalCuts.size()==1&&draft->internalCuts[0].minimum==1200001&&draft->protectedLocal[0].end==1500007,"internal pin and half-open protected rows preserved exactly");
    juce::MessageManager::getInstance()->runDispatchLoopUntil(30);expect(panel.readDraft(draft).isEmpty(),"deferred editor notifications cannot dirty saved frame rows");
    click(panel,"Confirm settings");expect(s.regionPlanning&&s.regionPlanning->digest()==draft->digest()&&s.revision==revision&&s.currentState.analysis->identity==snapshot,"explicit confirmation changes metadata only, preserving baseline measurement/audio");
    imageOf(panel,root.getChildFile("region-settings.png"));
    for(auto* c:panel.getChildren())if(auto* view=dynamic_cast<juce::Viewport*>(c)){view->setViewPosition(0,view->getViewedComponent()->getHeight());break;}
    imageOf(panel,root.getChildFile("region-settings-frames.png"));
    s.auditionBusy=[] {return true;};text(panel,"Maximum splits per original root (0-4)","4");control<juce::TextButton>(panel,"Confirm settings").onClick();expect(s.regionPlanning->maximumSplits==2,"confirmation checks stopped transport even before next UI refresh");s.auditionBusy={};
    click(panel,"Discard edits / reload");expect(control<juce::TextEditor>(panel,"Maximum splits per original root (0-4)").getText()=="2","discard reloads active values");
    auto changed=*s.regionPlanning;changed.maximumSplits=3;juce::String error;expect(s.setRegionPlanning(changed,error)&&panel.readDraft(draft).isNotEmpty(),"external settings change makes private editor stale");
    click(panel,"Discard edits / reload");choose(panel,"Region mode",2);expect(panel.readDraft(draft).isEmpty()&&draft->sectionId==s.sectionDocument.map().sections.front().id,"selected section binds original ID, not a renamed processing region");
    expect(s.sectionDocument.split(s.sectionDocument.map().sections.front().id,48000,error)&&panel.readDraft(draft).isNotEmpty(),"edited selected annotation bounds require explicit editor reload");
    click(panel,"Discard edits / reload");choose(panel,"Proposal method",1);click(panel,"Confirm settings");expect(!s.regionPlanning&&s.currentState.analysis->identity==snapshot,"explicit default-method selection removes M20 without changing audio/measurement");
    const auto runtime=s.sessionId;s.sessionId="new-runtime";expect(panel.readDraft(draft).isNotEmpty(),"old editor cannot apply into a different session");s.sessionId=runtime;
}
void workflowChecks(StandaloneMainComponent& ui,const juce::File& root)
{
    auto& s=Access::model(ui);auto requirements=s.jobRequirements;requirements.maximumLraTenths=49;juce::String error;expect(s.setJobRequirements(requirements,error),"UI workflow has explicit LRA requirement");
    RegionPlanningPanel editor(s);choose(editor,"Proposal method",2);click(editor,"Confirm settings");Access::refresh(ui);
    const auto source=s.source.analysis->identity,baseline=s.currentState.analysis->identity;const auto revision=s.revision;
    Heartbeat heartbeat;heartbeat.begin();Access::propose(ui);expect(s.job.isActive(),"experimental Propose goes directly through asynchronous existing worker");waitJob(ui);
    expect(s.sectionProposal&&s.sectionProposal->regionReview()&&s.canApplySectionProposal(),"main completion presents measured region recommendation through existing Apply boundary");
    expect(!s.beginAppliedBoundary(1),"pending measured proposal protects Current from competing boundary edits");
    auto result=s.sectionProposal;const auto review=result->regionReview();expect(review->retainedAttempt<review->attempts.size(),"retained measured attempt index is checked before UI inspection");
    expect(s.revision==revision&&s.currentState.analysis->identity==baseline&&s.source.analysis->identity==source,"proposal presentation never installs private candidate as Current");
    auto& lane=Access::lane(ui);expect(lane.preferredHeight()>TimelineGainLane::height,"M20 shows separate geometry bands above applied/proposed envelopes");
    expect(lane.getHeight()==lane.preferredHeight()&&lane.plotBounds().getBottom()<=lane.getHeight()-28,"geometry bands and gain controls occupy separate layout space");
    lane.inspectRegionAttempt(review->retainedAttempt);lane.selectFrame(review->attempts[review->retainedAttempt].plan.activity.support.front().begin,false);
    expect(lane.inspectionText().contains("Changed-support")&&lane.inspectionText().contains("nominal plateau")&&lane.inspectionText().contains("Not Current verification"),"gain inspection distinguishes support, nominal cuts and measurement domains");
    const auto oldDomain=s.auditionDomain;RegionProposalReviewPanel panel(result,[&]{return s.sectionProposalStale()?"STALE":"Current binding";});
    panel.onInspect=[&](size_t n,int64_t begin,int64_t){lane.inspectRegionAttempt(n);lane.selectFrame(begin,false);};
    choose(panel,"Measured region attempt",int(review->retainedAttempt+2));click(panel,"Inspect in timeline");
    expect(s.auditionDomain==oldDomain&&Access::stopped(ui)&&s.revision==revision,"attempt/region inspection does not switch audition, connect, seek or apply");
    imageOf(panel,root.getChildFile("region-review.png"));imageOf(ui,root.getChildFile("main-regions.png"));
    panel.addToDesktop(juce::ComponentPeer::windowHasTitleBar);panel.setVisible(true);
    const auto inspectUntil=juce::Time::getMillisecondCounterHiRes()+1000;std::vector<double> completeSwitch,actualSwitch;
    while(juce::Time::getMillisecondCounterHiRes()<inspectUntil)
    {
        auto begin=juce::Time::getMillisecondCounterHiRes();choose(panel,"Measured region attempt",1);completeSwitch.push_back(juce::Time::getMillisecondCounterHiRes()-begin);
        begin=juce::Time::getMillisecondCounterHiRes();choose(panel,"Measured region attempt",int(review->retainedAttempt+2));actualSwitch.push_back(juce::Time::getMillisecondCounterHiRes()-begin);juce::MessageManager::getInstance()->runDispatchLoopUntil(20);
    }
    std::sort(completeSwitch.begin(),completeSwitch.end());std::sort(actualSwitch.begin(),actualSwitch.end());
    auto switchTiming=juce::var(new juce::DynamicObject());switchTiming.getDynamicObject()->setProperty("complete_review_switch_p95_ms",completeSwitch[size_t(std::ceil(.95*completeSwitch.size()))-1]);switchTiming.getDynamicObject()->setProperty("actual_review_switch_p95_ms",actualSwitch[size_t(std::ceil(.95*actualSwitch.size()))-1]);cases.add(switchTiming);
    const auto visibleReport=[&]() -> ProposalInspectionPanel&
    {for(auto* child:panel.getChildren())if(auto* report=dynamic_cast<ProposalInspectionPanel*>(child);report&&report->isVisible())return *report;throw std::runtime_error("No visible review page");};
    const auto pageLocation=[](ProposalInspectionPanel& report)
    {for(auto* child:report.getChildren())if(auto* label=dynamic_cast<juce::Label*>(child);label&&label->getText().startsWith("Complete inspection - page"))return label->getText();return juce::String();};
    choose(panel,"Measured region attempt",1);expect(visibleReport().pageCount()>1,"complete explanation has multiple fully accessible pages");visibleReport().nextPage();const auto location=pageLocation(visibleReport());
    choose(panel,"Measured region attempt",int(review->retainedAttempt+2));choose(panel,"Measured region attempt",1);
    expect(location.isNotEmpty()&&pageLocation(visibleReport())==location,"switching between immutable reports preserves the user's inspection page");
    panel.removeFromDesktop();heartbeat.finish();
    const auto originalMap=s.sectionDocument.map();auto dense=originalMap;dense.sections.clear();dense.nextId=10001;++dense.revision;
    for(int64_t n=0;n<10000;++n)dense.sections.push_back({uint64_t(n+1),n*dense.frames/10000,(n+1)*dense.frames/10000,"Annotation "+juce::String(n+1)});
    expect(s.sectionDocument.restore(dense,dense.inputIdentity)&&!s.sectionProposalStale(),"automatic measured region proposal survives unrelated dense annotation subdivisions");Access::refresh(ui);
    std::vector<double> paints;
    for(int n=0;n<12;++n){const auto start=juce::Time::getMillisecondCounterHiRes();auto screenshot=lane.createComponentSnapshot(lane.getLocalBounds());paints.push_back(juce::Time::getMillisecondCounterHiRes()-start);}
    std::sort(paints.begin(),paints.end());expect(paints[10]<=100,"all M20 geometry bands with 10000 annotations paint within 100 ms p95");
    auto denseTiming=juce::var(new juce::DynamicObject());denseTiming.getDynamicObject()->setProperty("m20_dense_geometry_p95_ms",paints[10]);cases.add(denseTiming);
    expect(s.sectionDocument.restore(originalMap,originalMap.inputIdentity),"restore original annotations after isolated paint stress");Access::refresh(ui);
    const auto appliedBeforeDelete=proposalPlanIdentity(s.currentStateDocument.processingPlan());const auto proposalBeforeDelete=s.sectionProposal;
    choose(lane,"Gain lane selection",3);lane.selectAt(48000*32,lane.plotBounds().getCentreY());
    expect(lane.deleteSelected()&&s.sectionProposal==proposalBeforeDelete&&!s.sectionTrimDraft&&proposalPlanIdentity(s.currentStateDocument.processingPlan())==appliedBeforeDelete,"Backspace on proposed geometry cannot mutate the measured candidate or applied processing");
    const auto report=result->text();expect(report.contains("Total splits")&&report.contains("maximum cumulative attenuation")&&report.contains("actual transition count")&&report.contains("separate whole-programme change")&&report.contains("omitted cut"),"complete explanation retains independent counts, cumulative gains, global scope and omitted search work");
    auto failedRequirements=requirements;failedRequirements.maximumLraTenths=1;expect(s.setJobRequirements(failedRequirements,error),"diagnostic request has a strict deliberately unreachable LRA rule");
    Access::propose(ui);waitJob(ui);auto failed=s.sectionProposal;
    expect(failed&&failed->status()==SectionProposalStatus::noSolution&&failed->regionReview()&&!failed->regionReview()->attempts.empty()&&!s.canApplySectionProposal(),"bounded measured failure remains inspectable with Apply unavailable");
    RegionProposalReviewPanel failedPanel(failed,[&]{return s.sectionProposalDisplayText().substring(0,100);});
    failedPanel.onInspect=[&](size_t n,int64_t begin,int64_t){lane.inspectRegionAttempt(n);lane.selectFrame(begin,false);Access::refresh(ui);};
    for(size_t n=0;n<failed->regionReview()->attempts.size();++n)
    {choose(failedPanel,"Measured region attempt",int(n+2));expect(lane.diagnosticRegionInspection()&&!control<juce::TextButton>(Access::proposal(ui),"Apply proposal").isEnabled(),"every displayed measured failure is diagnostic and cannot Apply");}
    imageOf(failedPanel,root.getChildFile("region-failure.png"));
    auto activeSettings=*s.regionPlanning;auto quiet=RegionPlanningSettings::defaults(RegionMode::selection,activeSettings.inputIdentity);quiet.selection={0,480000};
    expect(s.setRegionPlanning(quiet,error),"explicit quiet selection has global fixed");
    JobControl quietJob;auto quietResult=generateSectionProposal(s.sectionProposalRequest(),quietJob);
    expect(quietResult.value&&(*quietResult.value)->regionReview()&&(*quietResult.value)->regionReview()->attempts.empty(),"selection with no high correction need has no complete candidate or retained attempt");
    RegionProposalReviewPanel emptyPanel(*quietResult.value,[] {return "No complete measured attempt";});
    expect(!control<juce::TextButton>(emptyPanel,"Inspect in timeline").isEnabled(),"missing retained index cannot enable a no-op or invalid timeline inspection");
    expect(s.setRegionPlanning(activeSettings,error)&&s.setJobRequirements(requirements,error)&&s.acceptSectionProposal(result),"restore the exact original measured recommendation binding after diagnostic review");Access::refresh(ui);
    auto settings=*s.regionPlanning;settings.maximumSplits=1;expect(s.setRegionPlanning(settings,error)&&!s.canApplySectionProposal(),"settings edit stales measured proposal without changing baseline");Access::refresh(ui);expect(Access::proposal(ui).outcomeText().contains("Stale"),"main stale outcome is prominent");
    settings.maximumSplits=2;expect(s.setRegionPlanning(settings,error)&&s.canApplySectionProposal(),"exact original binding restores recommendation eligibility");Access::refresh(ui);click(Access::proposal(ui),"Apply proposal");
    expect(s.currentStateNeedsRemeasurement()&&!s.currentState.hasMeasurement&&s.lastAppliedExplanation,"explicit UI Apply leaves Current unverified and stores historical explanation");
    const auto appliedPlan=s.manualState();const auto historyBytes=encodeAppliedExplanation(*s.lastAppliedExplanation);
    const auto token=s.beginAppliedBoundary(1);expect(token&&s.previewAppliedBoundary(token,appliedPlan.gains->intervals[1].begin+20)&&s.commitAppliedBoundary(token),"boundary editing works immediately after explicit proposal Apply");
    expect(!s.canApplySectionProposal()&&s.appliedExplanationStatus().contains("processing edited since application")&&encodeAppliedExplanation(*s.lastAppliedExplanation)==historyBytes,"editing an applied proposal retires its applicability and preserves its measured explanation");
    expect(s.undoManual()&&s.manualState()==appliedPlan&&!s.canApplySectionProposal(),"Undo restores audio without resurrecting old proposal authority");
    expect(s.lastAppliedExplanation->narrative.contains("actual transition count"),"saved historical narrative includes added presentation accounting");
    auto current=measuredComposition(s,MeasurementDomain::currentState);expect(current.detail.lra.value==result->assessment()->detail.lra.value,"independent Current reproduces the exact privately measured plan");
    JobControl job;const auto path=root.getChildFile("ui-session.levelscope");expect(sf::save(sf::capture(s),path,job).succeeded(),"UI settings and applied history saved in schema9");auto document=sf::read(path,job);expect(document.succeeded(),"UI session read");auto reopened=sf::verify(*document.document,job);
    expect(reopened.prepared&&reopened.prepared->model.regionPlanning&&reopened.prepared->model.lastAppliedExplanation&&!reopened.prepared->model.sectionProposal,"reopened session restores settings/history, never a private candidate");
    expect(reopened.prepared->model.lastAppliedExplanation->narrative==s.lastAppliedExplanation->narrative,"full historical explanation survives exact save/reopen");
    TimelineGainLane restoredLane(reopened.prepared->model);restoredLane.setSize(1100,restoredLane.preferredHeight());restoredLane.synchronize("reopened-region");
    const auto& restoredPlan=*reopened.prepared->model.currentStateDocument.sectionGains;
    control<juce::ComboBox>(restoredLane,"Gain lane selection").setSelectedId(2,juce::sendNotificationSync);
    restoredLane.selectAt(restoredPlan.intervals[1].begin,restoredLane.plotBounds().getBottom()+5);
    expect(restoredLane.inspectionText().contains("Applied processing region")&&!control<juce::TextEditor>(restoredLane,"Timeline section gain dB").isEnabled(),"saved processing regions remain directly selectable after reopening without a live proposal");
    auto& removal=reopened.prepared->model;
    const auto oldHistory=encodeAppliedExplanation(*removal.lastAppliedExplanation);const auto oldGlobal=removal.currentStateDocument.staticGainDb;
    const auto nonzero=std::find_if(removal.currentStateDocument.sectionGains->intervals.begin(),removal.currentStateDocument.sectionGains->intervals.end(),[](const auto& r){return r.db!=0;});
    expect(nonzero!=removal.currentStateDocument.sectionGains->intervals.end()&&removal.draftRemoveAppliedRegion(nonzero->begin,nonzero->end),"reopened applied region stages removal independently of annotation geometry");
    expect(removal.applySectionTrimDraft()&&removal.currentStateNeedsRemeasurement()&&!removal.currentState.hasMeasurement&&removal.currentStateDocument.staticGainDb==oldGlobal,"reopened region removal explicitly applies and requires independent Current measurement");
    expect(encodeAppliedExplanation(*removal.lastAppliedExplanation)==oldHistory&&removal.appliedExplanationStatus().containsIgnoreCase("historical"),"region removal retains complete explanation as historical");
    const auto removedFile=root.getChildFile("removed-region.levelscope");expect(sf::save(sf::capture(removal),removedFile,job).succeeded(),"save applied region removal");auto removedRead=sf::read(removedFile,job);auto removedOpen=sf::verify(*removedRead.document,job);
    expect(removedOpen.prepared&&!removedOpen.prepared->model.currentState.hasMeasurement&&removedOpen.prepared->model.currentStateNeedsRemeasurement()&&encodeAppliedExplanation(*removedOpen.prepared->model.lastAppliedExplanation)==oldHistory,"reopened removal retains unverified Current and intact historical explanation");
    // Stress the new documentary text at the supported interval capacity. These
    // synthetic plans are presentation fixtures, never candidates for Apply.
    auto large=*review;large.constraints.frames=int64_t(10000)*48000;large.attempts.resize(3,review->attempts.front());
    for(auto& attempt:large.attempts)
    {
        auto& plan=attempt.plan;plan.assignments.clear();plan.envelope.frames=large.constraints.frames;plan.envelope.intervals.clear();
        for(int64_t n=0;n<10000;++n)plan.envelope.intervals.push_back({n*48000,(n+1)*48000,n%2?-1.:0.});
        plan.activity.intervalCount=10000;plan.activity.transitionCount=9999;
    }
    const auto largeText=large.text();auto historical=*s.lastAppliedExplanation;historical.narrative=largeText;
    const auto attachment=encodeAppliedExplanation(historical);
    expect(attachment.size()<AppliedProposalExplanation::maximumBytes&&decodeAppliedExplanation(attachment,job)->narrative==largeText,"three complete 10000-interval narratives fit the unchanged history budget and round trip exactly");
    Heartbeat pageHeartbeat;pageHeartbeat.begin();ProposalInspectionPanel largePanel(largeText,{},RegionProposalReviewPanel::pageBytes);largePanel.addToDesktop(juce::ComponentPeer::windowHasTitleBar);largePanel.setVisible(true);
    expect(largePanel.pageCount()>100&&largePanel.allText()==largeText,"large M20 explanation remains completely accessible through bounded pages");
    const auto pageUntil=juce::Time::getMillisecondCounterHiRes()+3000;
    while(juce::Time::getMillisecondCounterHiRes()<pageUntil){largePanel.nextPage();juce::MessageManager::getInstance()->runDispatchLoopUntil(20);}
    largePanel.removeFromDesktop();pageHeartbeat.finish();
    auto storage=juce::var(new juce::DynamicObject());storage.getDynamicObject()->setProperty("three_attempt_10000_interval_attachment_bytes",int64_t(attachment.size()));storage.getDynamicObject()->setProperty("complete_pages",largePanel.pageCount());cases.add(storage);
    // A fresh request cancellation retains already-applied processing.
    s.setCurrentGainDb(s.currentStateDocument.staticGainDb+.1);measuredComposition(s,MeasurementDomain::currentState);requirements.maximumLraTenths=30;expect(s.setJobRequirements(requirements,error),"force non-no-op request for cancellation check");
    const auto applied=proposalPlanIdentity(s.currentStateDocument.processingPlan());Access::propose(ui);const auto cancelledAt=juce::Time::getMillisecondCounterHiRes();Access::cancel(ui);waitJob(ui);
    expect(!s.sectionProposal&&s.job.state==JobState::cancelled&&proposalPlanIdentity(s.currentStateDocument.processingPlan())==applied,"main cancellation cannot apply a partial region result or roll back applied audio");
    const auto cancelMs=juce::Time::getMillisecondCounterHiRes()-cancelledAt;expect(cancelMs<=1000,"native region cancellation completes within one second");
    auto cancellation=juce::var(new juce::DynamicObject());cancellation.getDynamicObject()->setProperty("cancel_ms",cancelMs);cases.add(cancellation);
}

void overloadBannerChecks(const juce::File& root)
{
    const auto file=root.getChildFile("centre-overload.wav");juce::AudioBuffer<float> pcm(6,48000);pcm.clear();
    for(int n=0;n<48000;++n)pcm.setSample(2,n,.9f);floatFile(file,48000,pcm);JobControl job;const auto media=validateMedia(file,job);expect(media.succeeded(),"known Centre-only overload fixture");
    auto output=std::make_unique<Device>();auto* device=output.get();
    StandaloneMainComponent ui(emptyOutputCatalog(),noOutputSettings(),std::move(output));Access::stopTimer(ui);ui.setSize(1320,1100);ui.addToDesktop(juce::ComponentPeer::windowHasTitleBar);ui.setVisible(true);
    expect(!find<juce::ToggleButton>(ui,"Confirm pre-produced whole-programme linear-audio QC")&&Access::scope(ui).isVisible()&&Access::scope(ui).getText().contains("this file"),"main UI replaces scope confirmation with automatic file-limit explanation");
    auto& s=Access::model(ui);expect(s.setValidatedSource(file,media.media),"overload UI owns source fixture");s.sourceDocument.layout=layoutFor(6);expect(s.initializeCurrentStateFromSource(),"overload UI Current initialized");s.setCurrentGainDb(12);s.selectAuditionDomain(MeasurementDomain::currentState);
    Access::refresh(ui);for(int ch=0;ch<6;++ch)s.sourceDocument.routing.fileToDevice[ch]=ch;
    auto& a=Access::audio(ui);a.sourceChanged(media.media);const auto request=s.auditionRequest();a.synchronize(request,s.sourceDocument.routing);expect(a.connect(48000).isEmpty(),"connect injected output only");a.setMonitorGainDb(0);
    expect(a.play(request,s.sourceDocument.routing).isEmpty(),"explicit injected overload playback");awaitPlaying(a,*device);juce::AudioBuffer<float> out(12,256);device->tick(out);a.poll();Access::refreshTimeline(ui);
    auto& banner=Access::warning(ui);expect(banner.isVisible()&&banner.getText().contains("Centre (channel 3, output 3)")&&banner.getText().contains("Lower monitor volume")&&ui.getLocalBounds().contains(banner.getBounds()),"main window prominently exposes actual overload role/output/source position and recovery advice");
    expect(out.getMagnitude(0,256)==0&&!a.active()&&s.currentStateDocument.staticGainDb==12,"banner accompanies stopped silent output without changing applied processing");imageOf(ui,root.getChildFile("monitor-overload.png"));
    a.setMonitorGainDb(-30);a.stop();Access::refreshTimeline(ui);expect(banner.isVisible()&&!a.active(),"main warning survives volume adjustment and Stop without automatic playback");
    expect(a.play(request,s.sourceDocument.routing).isEmpty(),"explicit reduced-monitor recovery");awaitPlaying(a,*device);device->tick(out);a.poll();Access::refreshTimeline(ui);expect(!banner.isVisible()&&a.active(),"successful explicit recovery clears banner");a.stop();ui.removeFromDesktop();
}

void convenienceChecks(StandaloneMainComponent& ui,const juce::File& root)
{
    auto& s=Access::model(ui);const auto revision=s.revision;const auto source=s.source.analysis,current=s.currentState.analysis;
    Access::timeline(ui).focusForTransport();
    const juce::KeyPress shortcut('A',juce::ModifierKeys{},'a');
    expect(ui.keyPressed(shortcut)&&s.auditionDomain==MeasurementDomain::currentState&&Access::inspection(ui)==2&&Access::stopped(ui),"A pairs inspection and audition Current without connecting or starting playback");
    expect(ui.keyPressed(shortcut)&&s.auditionDomain==MeasurementDomain::source&&Access::inspection(ui)==1,"shortcut toggles both domains back to Original");
    expect(!ui.keyPressed(juce::KeyPress('A',juce::ModifierKeys::commandModifier|juce::ModifierKeys::shiftModifier,0))&&s.auditionDomain==MeasurementDomain::source,"old modified shortcut no longer toggles inspection/audition");
    juce::TextEditor typing;ui.addAndMakeVisible(typing);typing.setBounds(20,20,150,24);typing.grabKeyboardFocus();
    expect(!ui.keyPressed(shortcut)&&s.auditionDomain==MeasurementDomain::source,"typing A in a text field never switches inspection/audition");
    ui.removeChildComponent(&typing);Access::timeline(ui).focusForTransport();
    s.job.begin(JobStage::waveform,s.revision);ui.keyPressed(shortcut);expect(s.auditionDomain==MeasurementDomain::source,"active job blocks domain toggle");s.job.finish(s.job.ticket,s.revision,ErrorCode::none,"done");
    const auto busy=s.auditionBusy;s.auditionBusy=[]{return true;};ui.keyPressed(shortcut);expect(s.auditionDomain==MeasurementDomain::source,"active transport/preparation blocks shortcut without changing audition");s.auditionBusy=busy;
    expect(s.revision==revision&&s.source.analysis==source&&s.currentState.analysis==current,"inspection controls cannot alter measurements or processing");
    auto& viewport=Access::inspector(ui);const int before=viewport.getHeight();Access::moveDivider(ui,100);
    expect(viewport.getHeight()<before&&viewport.getViewedComponent()->getHeight()>=TimelineLayout::defaultInspector,"downward divider gives graph room while preserving all lower controls in viewport");
    viewport.setViewPosition(0,viewport.getViewedComponent()->getHeight());expect(viewport.getViewPositionY()>0,"compressed inspector controls remain scroll-accessible");viewport.setViewPosition(0,0);Access::moveDivider(ui,-100);
    expect(control<juce::Label>(ui,"Programme filename").getText()==s.sourceDocument.file.getFileName(),"programme filename appears above session status");
    imageOf(ui,root.getChildFile("inspection-convenience.png"));
}

void directCommandChecks(StandaloneMainComponent& ui,const juce::File& root)
{
    auto& s=Access::model(ui);s.discardSectionProposal();s.discardSectionTrimDraft();s.synchronizeSections();Access::refresh(ui);auto& timeline=Access::timeline(ui);
    expect(ui.getMenuBarNames()==juce::StringArray({"File","Edit","View","Transport","Analysis"}),"native menu model exposes the five command groups");
    expect(ui.getMenuForIndex(0,{}).getNumItems()>3&&ui.getMenuForIndex(4,{}).getNumItems()>1,"File and Analysis menus expose working existing actions");
    const auto editMenu=ui.getMenuForIndex(1,{});bool undoText=false,redoText=false;
    for(juce::PopupMenu::MenuItemIterator i(editMenu,false);i.next();){const auto& item=i.getItem();if(item.itemID==StandaloneMainComponent::undoEdit)undoText=item.text==juce::String::fromUTF8("Undo (⌘Z)");if(item.itemID==StandaloneMainComponent::redoEdit)redoText=item.text==juce::String::fromUTF8("Redo (⇧⌘Z)");}
    expect(undoText&&redoText,"native menu labels retain Command and Shift Unicode symbols");
    timeline.focusForTransport();const auto source=s.source.analysis;auto before=s.manualState();const auto id=before.map.sections.front().id;
    expect(s.setManualSectionGains({id},-1),"direct UI test applies a local gain");Access::refresh(ui);
    expect(!s.currentState.hasMeasurement&&!s.currentState.analysis&&s.currentStateNeedsRemeasurement()&&s.source.analysis==source,"direct change clears all Current measurement values, preserving Original");
    int unavailable=0;std::function<void(juce::Component&)> countUnavailable=[&](juce::Component& c){if(auto* label=dynamic_cast<juce::Label*>(&c);label&&label->getText()==juce::String::charToString(0x2014))++unavailable;for(auto* child:c.getChildren())countUnavailable(*child);};countUnavailable(ui);
    expect(unavailable==5,"all five unverified Current cells display the intended Unicode dash");
    expect(ui.keyPressed(juce::KeyPress('Z',juce::ModifierKeys::commandModifier,0))&&s.manualState()==before,"main Command-Z dispatches the shared edit history");
    expect(ui.keyPressed(juce::KeyPress('Z',juce::ModifierKeys::commandModifier|juce::ModifierKeys::shiftModifier,0))&&s.currentStateDocument.sectionGains,"main Shift-Command-Z redoes the committed edit");
    const auto warnings=timeline.peakVisibility().getToggleState();expect(ui.performCommand(StandaloneMainComponent::peakWarnings)&&timeline.peakVisibility().getToggleState()!=warnings,"View menu toggles the same peak-warning control");ui.performCommand(StandaloneMainComponent::peakWarnings);
    expect(ui.performCommand(StandaloneMainComponent::measureCurrent),"Analysis command starts explicit Current measurement");waitJob(ui);
    expect(s.currentState.hasMeasurement&&!s.currentStateNeedsRemeasurement()&&Access::stopped(ui),"explicit measurement restores verification without connecting playback");
    expect(ui.performCommand(StandaloneMainComponent::undoEdit)&&s.currentStateNeedsRemeasurement(),"Undo after measurement invalidates Current independently");
    timeline.focusForTransport();expect(ui.keyPressed(juce::KeyPress('A',juce::ModifierKeys{},'a'))&&Access::inspection(ui)==2&&timeline.isGainPreview()&&!timeline.displayedWaveform().empty(),"unverified Current keeps a gain-adjusted waveform from its last completed cache");
    imageOf(ui,root.getChildFile("direct-current-status.png"));
    Access::globalGain(ui).setValue(4,juce::sendNotificationSync);Access::refresh(ui);
    Access::globalGain(ui).setValue(5,juce::sendNotificationSync);Access::refresh(ui);
    expect(s.currentStateDocument.staticGainDb==5&&timeline.isGainPreview(),"global gain 5 dB leaves an unverified Current waveform preview");
    auto& remeasure=control<juce::TextButton>(ui,"Re-measure Current State");expect(remeasure.isEnabled(),"Current remeasurement is available after global gain entry");
    remeasure.setWantsKeyboardFocus(true);remeasure.grabKeyboardFocus();
    expect(remeasure.hasKeyboardFocus(false),"regression starts with keyboard focus on Re-measure");remeasure.onClick();
    // Paint synchronously before a timer refresh can replace the cleared overview.
    // This is the transition from the reported native paint-thread crash.
    imageOf(ui,root.getChildFile("remeasure-after-global-gain.png"));
    expect(s.job.isActive()&&Access::stopped(ui)&&!timeline.isGainPreview()&&timeline.displayedWaveform().empty(),"immediate remeasurement repaint has retired the preview and preserves stopped transport");
    waitJob(ui);expect(s.currentState.hasMeasurement&&!s.currentStateNeedsRemeasurement()&&s.currentState.analysis->context.processing.staticGainDb==5,"global-gain remeasurement completes and verifies the exact applied 5 dB setting");
    expect(timeline.hasKeyboardFocus(false),"Re-measure transfers focus to the stable timeline before disabling its button");
    const auto undo=juce::KeyPress('Z',juce::ModifierKeys::commandModifier,0),redo=juce::KeyPress('Z',juce::ModifierKeys::commandModifier|juce::ModifierKeys::shiftModifier,0);
    const auto endpoint=s.manualState();const auto history=s.manualHistory->undoCount();
    expect(history>=2&&ui.getPeer()->handleKeyPress(undo)&&ui.getPeer()->handleKeyPress(undo)&&s.manualHistory->undoCount()==history-2,"two native Command-Z events after measurement work without another mouse click");
    expect(s.currentStateNeedsRemeasurement()&&!s.currentState.hasMeasurement&&s.source.analysis==source,"Undo preserves Original and invalidates Current authority");
    expect(ui.getPeer()->handleKeyPress(redo)&&ui.getPeer()->handleKeyPress(redo)&&s.manualState()==endpoint,"native Shift-Command-Z restores both edits after measurement");
    const bool shading=timeline.transitionShading();expect(ui.performCommand(StandaloneMainComponent::transitionBands)&&timeline.transitionShading()!=shading,"View menu toggles waveform gain transition bands");ui.performCommand(StandaloneMainComponent::transitionBands);
    expect(s.stageBoundary(0,s.sectionDocument.map().sections.front().end+10),"main window can stage a shared boundary batch");Access::refresh(ui);
    expect(!remeasure.isEnabled()&&!ui.performCommand(StandaloneMainComponent::measureCurrent)&&!ui.performCommand(StandaloneMainComponent::newDocument)&&!ui.performCommand(StandaloneMainComponent::saveDocument)&&Access::stopped(ui),"pending boundaries block measurement, document replacement and saving while playback stays stopped");
    bool boundaryReason=false;for(auto* c:Access::proposal(ui).getChildren())if(auto* label=dynamic_cast<juce::Label*>(c);label&&label->getText().startsWith("Propose unavailable."))boundaryReason=label->getText().contains("Finish boundary edits with Enter, or cancel with Esc.")&&!label->getText().contains("Stop playback");
    expect(boundaryReason,"disabled proposal actions name the pending boundary edit instead of stopped playback");
    imageOf(ui,root.getChildFile("main-boundary-pending.png"));s.cancelBoundaryBatch();Access::refresh(ui);
    Access::globalGain(ui).setValue(6,juce::sendNotificationSync);Access::refresh(ui);remeasure.onClick();Access::cancel(ui);
    imageOf(ui,root.getChildFile("remeasure-cancel-after-gain.png"));waitJob(ui);
    expect(s.job.state==JobState::cancelled&&s.currentStateNeedsRemeasurement()&&!s.currentState.hasMeasurement&&Access::stopped(ui),"cancelling preparation after a gain edit paints safely and cannot verify Current");
}

juce::var gateInspectionReport(const measurement::Result& original,const measurement::Result& current)
{
    auto report=juce::var(new juce::DynamicObject());
    for(bool isCurrent:{false,true})
    {
        const auto& r=isCurrent?current:original;const auto gate=measurement::relativeGateEnergy(r.lraWindows,true);
        auto row=juce::var(new juce::DynamicObject());auto* d=row.getDynamicObject();int visible=0,excluded=0;std::optional<double> minimum;juce::Array<juce::var> excludedEnds;
        for(const auto& w:r.lraWindows)if(w.endFrame<=r.frames&&w.energy>0)
        {
            ++visible;const auto loudness=measurement::OfflineMeasurement::loudness(w.energy);minimum=minimum?std::min(*minimum,loudness):loudness;
            if(!gate||!measurement::includedInPopulation(w.energy,*gate,true)){++excluded;excludedEnds.add(double(w.endFrame)/r.sampleRate);}
        }
        d->setProperty("relative_lra_gate_lufs",r.lraRelativeGate?juce::var(*r.lraRelativeGate):juce::var());d->setProperty("minimum_visible_s_lufs",minimum?juce::var(*minimum):juce::var());
        d->setProperty("finite_visible_s_windows",visible);d->setProperty("excluded_visible_s_windows",excluded);d->setProperty("excluded_s_window_end_seconds",excludedEnds);
        MeasurementTimeline timeline;timeline.setSize(1150,540);timeline.setProgramme(r.frames,r.sampleRate,"gate-report");timeline.setAnalysis(proposal::analyzeProgramme(r,"gate-report"));
        size_t plotted=0,plottedExcluded=0;for(const auto& bin:timeline.shortTermGatePlot()){plotted+=bin.included.count+bin.excluded.count;plottedExcluded+=bin.excluded.count;}
        expect(plotted==size_t(visible)&&plottedExcluded==size_t(excluded),"actual comparison short-term plotted membership exactly matches measured population");
        report.getDynamicObject()->setProperty(isCurrent?"current":"original",row);
    }
    const auto og=measurement::relativeGateEnergy(original.momentary,false),cg=measurement::relativeGateEnergy(current.momentary,false);
    expect(og&&cg&&original.momentary.size()==current.momentary.size(),"paired comparison has matching measured momentary windows and gates");
    juce::Array<juce::var> changed;
    MeasurementTimeline timeline;timeline.setSize(1150,540);timeline.setProgramme(current.frames,current.sampleRate,"gate-changes");timeline.setAnalysis(proposal::analyzeProgramme(current,"gate-changes"));timeline.setReferenceMeasurement(std::make_shared<measurement::Result>(original),true);
    for(size_t n=0;n<current.momentary.size();++n)
    {
        const auto& o=original.momentary[n];const auto& c=current.momentary[n];
        const bool oi=measurement::includedInPopulation(o.energy,*og,false),ci=measurement::includedInPopulation(c.energy,*cg,false);if(oi==ci)continue;
        auto row=juce::var(new juce::DynamicObject());auto* d=row.getDynamicObject();const auto begin=c.endFrame-current.sampleRate/10;
        d->setProperty("begin_seconds",double(begin)/current.sampleRate);d->setProperty("end_seconds",double(c.endFrame)/current.sampleRate);d->setProperty("original_included",oi);d->setProperty("current_included",ci);
        d->setProperty("original_m_lufs",measurement::OfflineMeasurement::loudness(o.energy));d->setProperty("current_m_lufs",measurement::OfflineMeasurement::loudness(c.energy));changed.add(row);
        timeline.setViewRange(begin+current.sampleRate*.02,c.endFrame-current.sampleRate*.02);
        const auto& bins=timeline.changedMomentaryColumns();expect(!bins.empty()&&bins[bins.size()/2],"actual changed-membership cell remains red when zoomed inside its sampled cadence");
    }
    report.getDynamicObject()->setProperty("original_il_gate_lufs",original.integratedRelativeGate?juce::var(*original.integratedRelativeGate):juce::var());
    report.getDynamicObject()->setProperty("current_il_gate_lufs",current.integratedRelativeGate?juce::var(*current.integratedRelativeGate):juce::var());
    report.getDynamicObject()->setProperty("membership_change_cells",changed);return report;
}

void refreshListening(const juce::File& from,const juce::File& root)
{
    for(auto method:{"m19","m20"})
    {
        JobControl preflight;const auto read=sf::read(from.getChildFile(method).getChildFile("inspection.levelscope"),preflight);
        expect(read.document&&read.document->lastAppliedExplanation&&read.document->source.file.getFileName()=="MusicInBackAndForegroundSomeDialog.wav","preflight rejects changed comparison documents before remeasurement or writing");
    }
    for(auto method:{"m19","m20"})
    {
        const auto input=from.getChildFile(method).getChildFile("inspection.levelscope");JobControl job;auto read=sf::read(input,job);expect(read.succeeded(),"preserved comparison session reads");
        expect(read.document&&read.document->lastAppliedExplanation&&read.document->source.file.getFileName()=="MusicInBackAndForegroundSomeDialog.wav","comparison still contains the original music programme and applied explanation; changed user session must be preserved");
        const auto before=read.document;const auto plan=proposalPlanIdentity(before->current.processingPlan());
        const auto history=encodeAppliedExplanation(*before->lastAppliedExplanation);const auto required=before->jobRequirements.digest();
        const auto expectedSource=before->sourceCache,expectedCurrent=before->currentCache;
        expect(bool(expectedSource)&&bool(expectedCurrent),"original comparison carries independent reference results");
        auto prepared=sf::verify(*before,job);StandaloneSessionModel s;juce::String notice;
        expect(prepared.prepared&&sf::adopt(*prepared.document,*prepared.prepared,s,notice),"normal session adoption validates original media");
        measuredComposition(s,MeasurementDomain::source);measuredComposition(s,MeasurementDomain::currentState);
        const auto same=[](const auto& a,const auto& b)
        {return a.integrated.value==b.integrated.value&&a.lra.value==b.lra.value&&a.truePeak.value==b.truePeak.value&&a.samplePeak.value==b.samplePeak.value&&a.momentaryMax.value==b.momentaryMax.value&&a.shortTermMax.value==b.shortTermMax.value&&a.samplePeakFrames==b.samplePeakFrames&&a.truePeakEighthFrames==b.truePeakEighthFrames;};
        expect(same(s.source.analysis->detail,expectedSource->detail)&&same(s.currentState.analysis->detail,expectedCurrent->detail),"UI-only build independently reproduces all six scalars and peak positions for Original and Current");
        expect(proposalPlanIdentity(s.currentStateDocument.processingPlan())==plan&&s.jobRequirements.digest()==required&&encodeAppliedExplanation(*s.lastAppliedExplanation)==history,"comparison gains, requirements and complete historical explanation preserved exactly");
        expect(s.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass,"same applied comparison passes independently under new build identity");
        const auto folder=root.getChildFile(method);folder.createDirectory();const auto output=folder.getChildFile("inspection.levelscope");
        folder.getChildFile("gate-inspection.json").replaceWithText(juce::JSON::toString(gateInspectionReport(s.source.analysis->detail,s.currentState.analysis->detail),true));
        expect(sf::save(sf::capture(s),output,job).succeeded(),"save new verified comparison copy without replacing preserved original");
        const auto check=sf::read(output,job);const auto verified=sf::verify(*check.document,job);
        expect(verified.prepared&&verified.prepared->model.currentState.hasMeasurement,"new comparison reopens verified in matching build");
    }
}


#include "Qualification/M20/HoldoutChecks.h"

void recoverMusic(const juce::File& inputs,const juce::File& evidence,const juce::File& root)
{
    controlledHoldouts(inputs,root,"MusicInBackAndForegroundSomeDialog.wav");
    for(auto method:{"m19","m20"})
    {
        const auto fresh=juce::JSON::parse(root.getChildFile("MusicInBackAndForegroundSomeDialog").getChildFile(method).getChildFile("result.json"));
        const auto original=juce::JSON::parse(evidence.getChildFile(method).getChildFile("result.json"));
        for(auto key:{"content_sha256","requirements","preferences","settings","baseline","outcome","private_passes","apply_available","retained_actual_candidate","independent_current","inspection_only_plan"})
            expect(juce::JSON::toString(fresh[key],true)==juce::JSON::toString(original[key],true),(juce::String("Recovered fixed comparison exactly matches archived ")+method+" "+key).toRawUTF8());
        const auto file=root.getChildFile("MusicInBackAndForegroundSomeDialog").getChildFile(method).getChildFile("inspection.levelscope");
        JobControl job;const auto read=sf::read(file,job);expect(read.document&&read.document->lastAppliedExplanation,"recovered comparison has complete newly generated explanation");
        const auto verified=sf::verify(*read.document,job);expect(verified.prepared&&verified.prepared->model.currentState.hasMeasurement,"recovered comparison reopens independently verified in new UI build");
    }
}

void exchangeChecks(const juce::File& from,const juce::File& root)
{
    for(auto layout:{"stereo","surround"})
    {
        const auto folder=from.getChildFile(layout);JobControl job;auto document=sf::read(folder.getChildFile("session.levelscope"),job);
        expect(document.document&&document.document->loadedSchema==9&&document.document->regionPlanning&&document.document->lastAppliedExplanation&&document.document->sourceCache,"foreign schema9 settings, history and source provenance are present");
        const auto& before=*document.document;const auto settings=before.regionPlanning->digest(),processing=proposalPlanIdentity(before.current.processingPlan());
        const auto requirements=before.jobRequirements.digest(),preferences=before.proposalPolicy.digest();const auto history=encodeAppliedExplanation(*before.lastAppliedExplanation);const auto roles=before.source.layout;
        const auto origin=before.sourceCache->context.platformAndFloatingPointMode;
        expect(sf::relinkAudio(*document.document,folder.getChildFile("programme.wav")),"foreign owned PCM is explicitly relinked for byte verification");
        auto prepared=sf::verify(*document.document,job);StandaloneSessionModel s;juce::String notice;
        expect(prepared.prepared&&sf::adopt(*prepared.document,*prepared.prepared,s,notice),"foreign programme opens through normal verification/adoption boundary");
        const auto receiver=s.measurementContext(MeasurementDomain::source).platformAndFloatingPointMode;
        expect(origin!=receiver,"session exchange crosses actual platform identity");
        expect(!s.source.hasMeasurement&&!s.currentState.hasMeasurement&&!s.sectionProposal,"foreign caches and private candidates confer no current verification or Apply authority");
        expect(s.regionPlanning->digest()==settings&&proposalPlanIdentity(s.currentStateDocument.processingPlan())==processing,"exact M20 settings and applied plan survive platform change");
        expect(s.jobRequirements.digest()==requirements&&s.proposalPolicy.digest()==preferences&&encodeAppliedExplanation(*s.lastAppliedExplanation)==history,"requirements, preferences and complete historical explanation are byte exact across Macs");
        expect(s.sourceDocument.layout.confirmed&&s.sourceDocument.layout.roles==roles.roles&&s.sourceDocument.layout.count==roles.count,"programme speaker roles survive independently of physical routing");
        expect(std::all_of(s.sourceDocument.routing.fileToDevice.begin(),s.sourceDocument.routing.fileToDevice.end(),[](int route){return route==-1;}),"portable session cannot restore another Mac's physical outputs");
        measuredComposition(s,MeasurementDomain::source);measuredComposition(s,MeasurementDomain::currentState);
        expect(s.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass,"foreign applied M20 plan independently passes on receiving Mac");
        const auto saved=root.getChildFile(juce::String(layout)+".levelscope");expect(sf::save(sf::capture(s),saved,job).succeeded(),"receiving platform saves independently verified state");auto read=sf::read(saved,job);auto reopened=sf::verify(*read.document,job);
        expect(reopened.prepared&&reopened.prepared->model.currentState.hasMeasurement&&reopened.prepared->model.lastAppliedExplanation->narrative==s.lastAppliedExplanation->narrative,"receiving platform save/reopen preserves verification and complete explanation separately");
        auto record=juce::var(new juce::DynamicObject());auto* r=record.getDynamicObject();r->setProperty("layout",layout);r->setProperty("origin",origin);r->setProperty("receiver",receiver);r->setProperty("history_sha256",juce::SHA256(history.data(),history.size()).toHexString());r->setProperty("settings",settings);r->setProperty("processing",processing);cases.add(record);
    }
}

void controlVariants(const juce::File& root)
{
    auto s=programme(root.getChildFile("variant-fixture"),2,true);
    auto requirements=s.jobRequirements;requirements.maximumLraTenths=49;
    juce::String error;expect(s.setJobRequirements(requirements,error),"variant requirement set without changing baseline");
    auto policy=s.proposalPolicy;policy.durationThousandths=600;expect(s.setProposalPolicy(policy,error),"explicit 60 percent allowance for two high portions");
    const auto baseline=s.currentState.analysis;const auto processing=proposalPlanIdentity(s.currentStateDocument.processingPlan());
    auto base=RegionPlanningSettings::defaults(RegionMode::selection,sectionInputIdentity(s.measurementContext(MeasurementDomain::source)));
    base.selection={int64_t(48000)*12,int64_t(48000)*46};base.allowGlobal=true;
    std::vector<std::pair<juce::String,RegionPlanningSettings>> variants{{"approved-defaults",base}};
    const auto add=[&](juce::String name,auto change){auto settings=base;change(settings);variants.emplace_back(name,std::move(settings));};
    for(int ms:{1000,6000})add("minimum-"+juce::String(ms),[&](auto& v){v.minimumIntervalMs=ms;});
    for(int splits:{0,1,4})add("splits-"+juce::String(splits),[&](auto& v){v.maximumSplits=splits;});
    for(int ms:{500,1500,3000})add("margin-"+juce::String(ms),[&](auto& v){v.leftMarginMs=v.rightMarginMs=ms;});
    for(int tenths:{5,30})add("need-tenths-"+juce::String(tenths),[&](auto& v){v.needThresholdTenths=tenths;});
    for(int ms:{300,3000})add("persistence-"+juce::String(ms),[&](auto& v){v.persistenceMs=ms;});
    add("boundary-refinement-off",[](auto& v){v.allowRefinement=false;});
    add("whole-programme-global-fixed",[](auto& v){v.allowGlobal=false;});
    for(const auto& variant:variants)
    {
        expect(s.setRegionPlanning(variant.second,error),"variant accepted through immutable settings boundary");
        RegionPlanningPanel editor(s);std::optional<RegionPlanningSettings> draft;
        expect(editor.readDraft(draft).isEmpty()&&draft&&draft->digest()==variant.second.digest(),"variant survives native editor reload exactly");
        JobControl job;const auto outcome=generateSectionProposal(s.sectionProposalRequest(),job);
        expect(outcome.value&&(*outcome.value)->regionReview(),"bounded variant returns complete inspectable outcome");
        const auto& result=**outcome.value;const auto review=result.regionReview();
        expect(result.passes()<=3&&review->proxyEvaluations<=64&&review->planningBytes<=64*1024*1024,"variant retains measurement, proxy and storage bounds");
        expect(s.currentState.analysis==baseline&&proposalPlanIdentity(s.currentStateDocument.processingPlan())==processing,"variant cannot mutate Current or applied processing");
        expect(s.acceptSectionProposal(*outcome.value),"variant completion has exact current request binding");
        expect(s.canApplySectionProposal()==(result.status()==SectionProposalStatus::ready),"only measured ready outcome authorizes Apply for every variant");
        auto record=juce::var(new juce::DynamicObject());auto* row=record.getDynamicObject();row->setProperty("variant",variant.first);row->setProperty("settings",variant.second.encode());row->setProperty("outcome",result.outcomeTitle());row->setProperty("passes",result.passes());row->setProperty("proxy_count",int(review->proxyEvaluations));row->setProperty("planning_seconds",review->planningSeconds);row->setProperty("planning_bytes",int64_t(review->planningBytes));
        std::vector<RegionReviewScore> scores;juce::Array<juce::var> attempts;
        for(const auto& attempt:review->attempts)
        {
            scores.push_back(regionMeasuredScore(attempt.evaluation,attempt.plan.activity,attempt.plan.identity));
            auto actual=juce::var(new juce::DynamicObject());auto* a=actual.getDynamicObject();a->setProperty("identity",attempt.plan.identity);a->setProperty("pass",attempt.evaluation.verdict==Verdict::pass);a->setProperty("lra",attempt.measured[2].value?juce::var(*attempt.measured[2].value):juce::var());a->setProperty("affected_frames",attempt.plan.activity.affectedFrames);a->setProperty("connected_regions",int(attempt.plan.activity.support.size()));a->setProperty("transitions",int(attempt.plan.activity.transitionCount));attempts.add(actual);
            if(attempt.evaluation.verdict==Verdict::pass)expect(attempt.measured[2].value&&*attempt.measured[2].value<=4.9,"variant pass is strict actual LRA, not a prediction or tolerance");
        }
        row->setProperty("actual_attempts",attempts);juce::Array<juce::var> bands;
        for(double band:{0.,.001,.01,.05})
        {auto record=juce::var(new juce::DynamicObject());auto* b=record.getDynamicObject();b->setProperty("diagnostic_band_lu",band);auto selected=selectRegionReview(scores,band);b->setProperty("attempt",selected?int(*selected):-1);bands.add(record);}
        row->setProperty("diagnostic_only_band_comparison",bands);row->setProperty("retained_attempt",review->retainedAttempt<review->attempts.size()?int(review->retainedAttempt):-1);cases.add(record);
        root.getChildFile(variant.first+".txt").replaceWithText(result.text());
    }
}

#include "EventUIIntegrationChecks.h"
#include "InnerUIIntegrationChecks.h"
#include "AppliedBoundaryUIChecks.h"
#include "ProposalOperationChecks.h"
#include "MonoLraWorkflowChecks.h"
#include "PlaybackPreparationChecks.h"
#include "GraphBoundaryUIChecks.h"
}
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI gui;const bool monoLra=argc==3&&juce::String(argv[2])=="--mono-lra",graph=argc==4&&juce::String(argv[2])=="--graph-ui",playback=argc==3&&juce::String(argv[2])=="--playback-ui",operation=argc==4&&juce::String(argv[2])=="--operation-ui",applied=argc==4&&juce::String(argv[2])=="--applied-ui",inner=argc==4&&juce::String(argv[2])=="--inner-ui",events=argc==3&&juce::String(argv[2])=="--event-ui",variants=argc==3&&juce::String(argv[2])=="--variants",exchange=argc==4&&juce::String(argv[2])=="--exchange",holdouts=argc==4&&juce::String(argv[2])=="--holdouts",listening=argc==4&&juce::String(argv[2])=="--refresh-listening",recover=argc==5&&juce::String(argv[2])=="--recover-music";if(argc!=2&&!variants&&!exchange&&!holdouts&&!listening&&!recover&&!events&&!inner&&!applied&&!operation&&!playback&&!graph&&!monoLra)return 2;const juce::File root(argv[1]);if(root.exists())return 2;root.createDirectory();bool passed=false;juce::String error;
    try
    {
        if(monoLra)monoLraWorkflowChecks(root);else if(graph){appliedBoundaryUIChecks(root,juce::File(argv[3]),true);for(int channels:{2,6})graphBoundaryPlayback(root,channels);}else if(playback)playbackPreparationChecks(root);else if(operation)proposalOperationChecks(root,juce::File(argv[3]));else if(applied)appliedBoundaryUIChecks(root,juce::File(argv[3]));else if(inner)innerUIChecks(root,juce::File(argv[3]));else if(events)eventUIChecks(root);else if(recover)recoverMusic(juce::File(argv[3]),juce::File(argv[4]),root);else if(listening)refreshListening(juce::File(argv[3]),root);else if(holdouts)controlledHoldouts(juce::File(argv[3]),root);else if(exchange)exchangeChecks(juce::File(argv[3]),root);else if(variants)controlVariants(root);else
        {
        overloadBannerChecks(root);auto s=programme(root.getChildFile("fixture"),2);editorChecks(s,root);
        auto ui=std::make_unique<StandaloneMainComponent>(emptyOutputCatalog(),noOutputSettings());Access::stopTimer(*ui);ui->setSize(1320,1100);ui->setName("LevelScope - isolated M20 integration test");ui->addToDesktop(juce::ComponentPeer::windowHasTitleBar);ui->setVisible(true);
        auto busy=Access::model(*ui).auditionBusy;Access::model(*ui)=s;Access::model(*ui).auditionBusy=busy;Access::refresh(*ui);convenienceChecks(*ui,root);workflowChecks(*ui,root);directCommandChecks(*ui,root);
        }
        passed=true;
    }
    catch(const std::exception& e){error=e.what();}
    auto data=juce::var(new juce::DynamicObject());auto* d=data.getDynamicObject();d->setProperty("passed",passed);d->setProperty("error",error);d->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);d->setProperty("checks",checks);d->setProperty("cases",cases);d->setProperty("physical_devices_opened",false);
    root.getChildFile("result.json").replaceWithText(juce::JSON::toString(data,true));std::cout<<(passed?"PASS":"FAIL")<<" "<<checks.size()<<" "<<error<<std::endl;return passed?0:1;
}

#include "CompositionTestFixtures.h"
#include "Standalone/StandaloneMainComponent.h"
#include "Standalone/RegionReviewFormatting.h"

namespace levelscope::standalone
{
struct PreparationFeedbackTestAccess
{
    static StandaloneSessionModel& model(StandaloneMainComponent& c){return c.session;}
    static void stopTimer(StandaloneMainComponent& c){c.stopTimer();}
    static void tick(StandaloneMainComponent& c){c.timerCallback();}
    static void refresh(StandaloneMainComponent& c){c.refreshFromSession();c.refreshTimeline();c.resized();}
    static void propose(StandaloneMainComponent& c){c.proposeWithActiveMethod();}
    static void cancel(StandaloneMainComponent& c){c.cancelJobButton.onClick();}
    static TimelineGainLane& lane(StandaloneMainComponent& c){return c.gainLane;}
    static ProposalWorkflowPanel& proposal(StandaloneMainComponent& c){return c.proposalPanel;}
    static bool stopped(StandaloneMainComponent& c){return !c.audition.active()&&!c.audition.connected();}
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
    expect(s.lastAppliedExplanation->narrative.contains("actual transition count"),"saved historical narrative includes added presentation accounting");
    auto current=measuredComposition(s,MeasurementDomain::currentState);expect(current.detail.lra.value==result->assessment()->detail.lra.value,"independent Current reproduces the exact privately measured plan");
    JobControl job;const auto path=root.getChildFile("ui-session.levelscope");expect(sf::save(sf::capture(s),path,job).succeeded(),"UI settings and applied history saved in schema9");auto document=sf::read(path,job);expect(document.succeeded(),"UI session read");auto reopened=sf::verify(*document.document,job);
    expect(reopened.prepared&&reopened.prepared->model.regionPlanning&&reopened.prepared->model.lastAppliedExplanation&&!reopened.prepared->model.sectionProposal,"reopened session restores settings/history, never a private candidate");
    expect(reopened.prepared->model.lastAppliedExplanation->narrative==s.lastAppliedExplanation->narrative,"full historical explanation survives exact save/reopen");
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
    Heartbeat pageHeartbeat;pageHeartbeat.begin();ProposalInspectionPanel largePanel(largeText);largePanel.addToDesktop(juce::ComponentPeer::windowHasTitleBar);largePanel.setVisible(true);
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

}
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI gui;const bool variants=argc==3&&juce::String(argv[2])=="--variants",exchange=argc==4&&juce::String(argv[2])=="--exchange";if(argc!=2&&!variants&&!exchange)return 2;const juce::File root(argv[1]);if(root.exists())return 2;root.createDirectory();bool passed=false;juce::String error;
    try
    {
        if(exchange)exchangeChecks(juce::File(argv[3]),root);else if(variants)controlVariants(root);else
        {
        auto s=programme(root.getChildFile("fixture"),2);editorChecks(s,root);
        auto ui=std::make_unique<StandaloneMainComponent>(emptyOutputCatalog(),noOutputSettings());Access::stopTimer(*ui);ui->setSize(1320,1100);ui->setName("LevelScope - isolated M20 integration test");ui->addToDesktop(juce::ComponentPeer::windowHasTitleBar);ui->setVisible(true);
        auto busy=Access::model(*ui).auditionBusy;Access::model(*ui)=s;Access::model(*ui).auditionBusy=busy;Access::refresh(*ui);workflowChecks(*ui,root);
        }
        passed=true;
    }
    catch(const std::exception& e){error=e.what();}
    auto data=juce::var(new juce::DynamicObject());auto* d=data.getDynamicObject();d->setProperty("passed",passed);d->setProperty("error",error);d->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);d->setProperty("checks",checks);d->setProperty("cases",cases);d->setProperty("physical_devices_opened",false);
    root.getChildFile("result.json").replaceWithText(juce::JSON::toString(data,true));std::cout<<(passed?"PASS":"FAIL")<<" "<<checks.size()<<" "<<error<<std::endl;return passed?0:1;
}

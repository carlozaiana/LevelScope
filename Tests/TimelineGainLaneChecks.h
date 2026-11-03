// Included by TimelineProgressTests.cpp to exercise the same real UI/controller.
#include "Standalone/TimelineGainLane.h"
#include "Standalone/TimelineGainOutline.h"
#include "Standalone/StandaloneSectionsPanel.h"

namespace
{
juce::TextEditor* laneEditor(TimelineGainLane& lane)
{
    for(auto* c:lane.getChildren())if(auto* e=dynamic_cast<juce::TextEditor*>(c))return e;
    return nullptr;
}
void stageGain(TimelineGainLane& lane,const juce::String& value)
{
    auto* field=laneEditor(lane);expect(field&&field->isEnabled(),"main timeline gain field available");
    field->setText(value,false);field->onTextChange();field->onReturnKey();
}
juce::String laneNotice(TimelineGainLane& lane)
{
    juce::String text;
    for(auto* c:lane.getChildren())if(auto* l=dynamic_cast<juce::Label*>(c))text+=l->getText()+"\n";
    return text;
}
void gainLaneChecks(StandaloneSessionModel& s,const juce::File& root)
{
    s.synchronizeSections();juce::String error;auto& document=s.sectionDocument;
    expect(document.split(document.map().sections.front().id,96000,error)&&document.split(document.map().sections.back().id,192000,error),"prepare independent annotation sections");
    const auto original=s.source.analysis,verified=s.currentState.analysis;const auto revision=s.revision;
    TimelineGainLane lane(s);MeasurementTimeline ui;ui.setSize(1150,540);
    ui.attachGainLane(lane);ui.setProgramme(s.sourceDocument.media.frames,48000,"gain-source");
    WaveformProgress originalPcm(std::make_shared<waveform::Budget>());JobControl waveformJob;
    expect(OfflineAnalysisService(nullptr,&originalPcm).buildWaveform(s.measurementContext(Domain::source),Domain::source,waveformJob).value.has_value(),"gain lane shares a real Original waveform and accepted loudness view");
    ui.setWaveform(originalPcm.overview(),true,{"L","R"});ui.setAnalysis(proposal::analyzeProgramme(original->detail,original->identity));
    lane.synchronize("gain-source");
    auto point=lane.plotBounds().getCentre().toFloat()+lane.getPosition().toFloat();
    int seeks=0;ui.onRequestSeek=[&](int64_t,const juce::String&){++seeks;return juce::String{};};ui.setTransport(0,true,{},"unused");
    click(ui,point);expect(seeks==0&&s.revision==revision&&s.currentState.analysis==verified,"gain lane click only inspects");
    expect(lane.inspectionText().contains("[96000, 192000)")&&lane.inspectionText().contains("Applied: 0"),"lane click selects exact section and applied frame value");
    button(ui,"Fit selection")->onClick();expect(ui.viewRange()==std::make_pair(96000.,192000.),"gain selection fits the section on the shared waveform viewport");
    expect(lane.plotBounds().getX()==ui.graphBounds().getX()&&lane.plotBounds().getWidth()==ui.graphBounds().getWidth(),"gain lane shares exact horizontal scale with waveform");
    const auto view=ui.viewRange();ui.mouseDown(mouse(ui,point));ui.mouseDrag(mouse(ui,point.translated(60,0)));ui.mouseUp(mouse(ui,point.translated(60,0),false));
    expect(ui.viewRange()!=view&&seeks==0&&s.revision==revision,"gain drag pans without editing or seeking");button(ui,"Whole programme")->onClick();
    lane.selectFrame(96000,true);const double precise=-1.234567891234567;
    stageGain(lane,"-1.234567891234567");
    expect(s.sectionTrimDraft&&s.sectionTrimDraft->plan.intervals[1].db==precise&&s.currentState.analysis==verified&&s.revision==revision,"Return stages exact precision without applying or invalidating");
    expect(button(lane,"Apply gains")->isEnabled()&&laneNotice(lane).contains("DRAFT"),"draft visibly identified and explicit Apply enabled");
    laneEditor(lane)->setText("13",false);laneEditor(lane)->onTextChange();button(lane,"Set draft")->onClick();
    expect(s.sectionTrimDraft->plan.intervals[1].db==precise&&!button(lane,"Apply gains")->isEnabled(),"invalid unstaged value cannot accidentally apply previous draft");
    button(lane,"Discard gains")->onClick();expect(!s.sectionTrimDraft&&s.currentState.analysis==verified&&s.revision==revision,"Discard leaves exact applied audio and measurement intact");
    for(const auto* invalid:{"nan","inf","1 trailing","12.00001","-12.00001",""}){stageGain(lane,invalid);expect(!s.sectionTrimDraft,"invalid numerical text creates no draft");}
    stageGain(lane,"-1.234567891234567");lane.selectFrame(192000,true);stageGain(lane,"3.0000000000000004");
    expect(s.sectionTrimDraft->plan.intervals[2].db==3.0000000000000004,"next section precision retained");
    imageOf(ui,root.getChildFile("gain-draft.png"));
    button(lane,"Apply gains")->onClick();
    expect(!s.sectionTrimDraft&&s.currentStateNeedsRemeasurement()&&!s.currentState.hasMeasurement&&s.source.analysis==original&&s.revision==revision+1,"Apply marks only Current unverified and preserves Source");
    const auto appliedKey=s.currentStateDocument.sectionGains->audioKey();const auto& plan=*s.currentStateDocument.sectionGains;
    expect(plan.gainAt(96000-2400)==0&&plan.gainAt(96000+2399)==precise,"applied lane uses unchanged 100 ms transition endpoints");
    lane.selectFrame(96000,false);expect(lane.inspectionText().contains("at frame 96000")&&laneEditor(lane)->getText().getDoubleValue()==precise,"transition inspection and stored full-precision section setting remain distinct");
    s.measurementUnsuccessful(Domain::currentState,ErrorCode::cancelled);lane.synchronize("gain-source");
    expect(s.currentStateNeedsRemeasurement()&&s.currentStateDocument.sectionGains->audioKey()==appliedKey,"cancelled verification retains applied gain unverified");
    measure(s,Domain::currentState);lane.synchronize("gain-source");const auto remeasured=s.currentState.analysis;
    expect(s.currentState.hasMeasurement&&!s.currentStateNeedsRemeasurement(),"independent measurement verifies applied section processing");
    stageGain(lane,"-1.234567891234567");button(lane,"Apply gains")->onClick();
    expect(s.revision==revision+1&&s.currentState.analysis==remeasured,"identical Apply is a genuine no-op preserving verification");
    stageGain(lane,"-2");s.auditionBusy=[]{return true;};lane.synchronize("gain-source");
    expect(!button(lane,"Apply gains")->isEnabled()&&!laneEditor(lane)->isEnabled(),"active transport disables numerical mutation");
    button(lane,"Apply gains")->onClick();expect(s.currentState.analysis==remeasured,"active Apply callback cannot mutate");
    s.auditionBusy={};s.job.begin(JobStage::waveform,s.revision);lane.synchronize("gain-source");expect(!button(lane,"Apply gains")->isEnabled(),"jobs block Apply");
    s.job.finish(s.job.ticket,s.revision,ErrorCode::none,"done");lane.synchronize("gain-source");
    expect(document.moveEnd(document.map().sections[0].id,95000,error),"annotation-only boundary edit accepted");
    lane.synchronize("gain-source");
    expect(!button(lane,"Apply gains")->isEnabled()&&laneNotice(lane).contains("Stale")&&s.currentStateDocument.sectionGains->audioKey()==appliedKey&&s.currentState.analysis==remeasured,"geometry edits stale draft without moving applied ramps or losing measurement");
    button(lane,"Discard gains")->onClick();lane.selectFrame(100000,true);stageGain(lane,"-2");
    expect(s.sectionTrimDraft->plan.intervals[2].db==0,"changed annotations start a new zero baseline rather than moving old gains");
    button(lane,"Discard gains")->onClick();
    lane.selectFrame(document.map().frames,true);expect(lane.inspectionText().contains("at frame "+juce::String(document.map().frames-1)),"EOF lane click inspects last actual sample");
    imageOf(ui,root.getChildFile("gain-applied.png"));
    // Updates arriving from the existing Sections window also refresh this lane.
    expect(s.beginSectionTrimDraft()&&s.setSectionDraftGain(2,4.123456789123456),"external window stages exact draft");
    lane.synchronize("gain-source");lane.selectFrame(200000,false);
    expect(laneEditor(lane)->getText().getDoubleValue()==4.123456789123456,"main field observes existing window's draft revision");
    StandaloneSectionsPanel other(s);other.suspend();juce::TextEditor* otherField=nullptr;
    for(auto* child:other.getChildren())if(child->getName()=="Section draft gain dB")otherField=dynamic_cast<juce::TextEditor*>(child);
    lane.selectFrame(0,false);stageGain(lane,"1.234567891234567");other.synchronize();
    expect(otherField&&otherField->getText().getDoubleValue()==1.234567891234567,"existing Sections editor observes the main timeline's exact staged value");
    s.discardSectionTrimDraft();lane.synchronize("gain-source");other.synchronize();
    expect(otherField->getText().isEmpty(),"discard clears the old draft field in the Sections window");
    auto previousRuntime=s.sessionId;s.sessionId=juce::Uuid().toString();
    stageGain(lane,"2");expect(!s.sectionTrimDraft,"stale runtime cannot stage through an old visible binding");s.sessionId=previousRuntime;
    stageGain(lane,"1");
    s.sourceContentInvalid=true;lane.synchronize({});expect(!laneEditor(lane)->isEnabled()&&!button(lane,"Apply gains")->isEnabled(),"invalid source clears editable lane");
    button(lane,"Discard gains")->onClick();expect(!s.sectionTrimDraft,"explicit Discard can retire a displayed stale draft after source loss");
    s.sourceContentInvalid=false;lane.synchronize("gain-source");
}
void gainDenseChecks(const juce::File& root)
{
    auto req=fixture(root.getChildFile("dense.wav"),48000,2,240013);StandaloneSessionModel s;expect(s.setValidatedSource(req.context.sourceFile,req.context.media),"dense gain fixture");
    s.sourceDocument.layout=req.context.fileLayout;expect(s.initializeCurrentStateFromSource(),"manual current for dense lane");
    s.synchronizeSections();auto map=s.sectionDocument.map();map.sections.clear();map.nextId=10001;++map.revision;
    for(int n=0;n<10000;++n)map.sections.push_back({uint64_t(n+1),int64_t(n)*24,n==9999?map.frames:int64_t(n+1)*24,"Dense "+juce::String(n+1)});
    expect(s.sectionDocument.restore(map,map.inputIdentity)&&s.beginSectionTrimDraft(),"bounded 10000-section annotation and draft");
    for(size_t n=0;n<10000;++n)s.sectionTrimDraft->plan.intervals[n].db=n%2?12:-12;
    ++s.sectionTrimDraft->editSerial;
    TimelineGainLane lane(s);MeasurementTimeline ui;ui.setSize(1150,540);ui.attachGainLane(lane);ui.setProgramme(map.frames,48000,"dense");ui.setWaveform({},false,{"L","R"});lane.synchronize("dense");
    // Both sides are only 24 frames: the existing shortened-ramp rule is retained.
    const auto canonical=s.sectionTrimDraft->plan.canonical();expect(canonical.radius(1)==12&&canonical.gainAt(12)==-12&&canonical.gainAt(35)==12,"short-section transition keeps exact endpoints");
    auto reduced=gainOutline(&canonical,map.frames,0,map.frames,1);
    expect(reduced.size()<=4&&std::any_of(reduced.begin(),reduced.end(),[](auto k){return k.db==12;})&&std::any_of(reduced.begin(),reduced.end(),[](auto k){return k.db==-12;}),"one display column preserves both actual gain extremes");
    auto clip=gainOutline(&canonical,map.frames,24,30,9000);
    expect(clip.size()==2&&clip.front().db==canonical.gainAt(24)&&clip.back().db==canonical.gainAt(30),"viewport wholly inside a transition follows the processor at its actual sample positions");
    auto longView=gainOutline(&canonical,map.frames,0,map.frames,50000);
    expect(longView.size()<=4096*4&&std::is_sorted(longView.begin(),longView.end(),[](auto a,auto b){return a.frame<b.frame;}),"large display bounded to 4096 columns in frame order");
    auto start=juce::Time::getMillisecondCounterHiRes();auto image=ui.createComponentSnapshot(ui.getLocalBounds());const auto elapsed=juce::Time::getMillisecondCounterHiRes()-start;
    expect(elapsed<100,"dense 10000-section software preparation/paint under 100 ms");
    auto* perf=new juce::DynamicObject();perf->setProperty("dense_10000_sections_prepare_paint_ms",elapsed);cases.add(juce::var(perf));
    ui.setViewRange(10,38);lane.selectFrame(24,true);imageOf(ui,root.getChildFile("gain-short-transition.png"));
    auto seven=fixture(root.getChildFile("system-i.wav"),48000,8,48013);StandaloneSessionModel unsupported;
    expect(unsupported.setValidatedSource(seven.context.sourceFile,seven.context.media),"System I fixture");unsupported.sourceDocument.layout=seven.context.fileLayout;expect(unsupported.initializeCurrentStateFromSource(),"7.1 current global gain available");
    TimelineGainLane refused(unsupported);refused.synchronize("unsupported");expect(!laneEditor(refused)->isEnabled()&&laneNotice(refused).contains("5.1 only"),"7.1 section editing stays unsupported with a visible reason");
}
void transportFocusChecks(StandaloneSessionModel& s)
{
    auto device=std::make_unique<Device>();auto* io=device.get();StandaloneAudition audio(std::move(device));auto req=s.auditionRequest();auto routing=routes(2);
    expect(audio.connect(48000).isEmpty(),"focus test uses injected output only");audio.synchronize(req,routing);
    MeasurementTimeline ui;ui.setSize(1150,440);ui.setProgramme(req.context.media.frames,48000,"focus");ui.setWaveform({},false,{"L","R"});
    ui.onRequestSeek=[&](int64_t f,const juce::String& identity){return TimelineTransport::seek(f,identity,s,audio);};
    auto refresh=[&]{auto reason=TimelineTransport::blockedReason(s,audio);ui.setTransport(audio.positionFrames(),reason.isEmpty(),reason,audio.boundIdentity());};
    int toggles=0;ui.onTogglePlayback=[&]{++toggles;if(audio.active())audio.stop();else expect(audio.play(req,routing).isEmpty(),"Space starts existing transport");refresh();};
    s.auditionBusy=[&]{return audio.active();};refresh();
    ui.addToDesktop(0);ui.setVisible(true);ui.toFront(true);juce::MessageManager::getInstance()->runDispatchLoopUntil(20);
    auto* peer=ui.getPeer();expect(peer!=nullptr,"native JUCE focus surface available");
    auto* seek=button(ui,"Seek");seek->grabKeyboardFocus();seek->onClick();
    expect(ui.hasKeyboardFocus(false),"Seek transfers focus to a stable timeline surface");
    expect(peer->handleKeyPress(juce::KeyPress(juce::KeyPress::spaceKey))&&toggles==1,"Space after Seek starts playback through actual JUCE key dispatch");awaitPlaying(audio,*io);refresh();
    expect(!seek->isEnabled()&&ui.hasKeyboardFocus(false),"disabling Seek does not move focus into numeric editor");
    expect(peer->handleKeyPress(juce::KeyPress(juce::KeyPress::spaceKey))&&toggles==2&&!audio.active(),"second Space stops without mouse interaction");
    juce::TextEditor* text=nullptr;for(auto* c:ui.getChildren())if(auto* e=dynamic_cast<juce::TextEditor*>(c))text=e;
    expect(text!=nullptr,"numeric inspection editor");text->grabKeyboardFocus();peer->handleKeyPress(juce::KeyPress(juce::KeyPress::spaceKey));expect(toggles==2,"Space in focused text field remains text editing");
    // Starting from a keyboard-tabbed Seek button also moves focus before it disables.
    seek->grabKeyboardFocus();peer->handleKeyPress(juce::KeyPress(juce::KeyPress::spaceKey));awaitPlaying(audio,*io);refresh();peer->handleKeyPress(juce::KeyPress(juce::KeyPress::spaceKey));
    expect(toggles==4&&!audio.active(),"Space on a focused Seek button can start and stop twice");
    audio.disconnect();s.auditionBusy={};ui.removeFromDesktop();
}
}

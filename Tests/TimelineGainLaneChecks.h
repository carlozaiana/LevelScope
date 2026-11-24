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
    const auto original=s.source.analysis,verified=s.currentState.analysis;
    TimelineGainLane lane(s);MeasurementTimeline ui;ui.setSize(1150,540);ui.attachGainLane(lane);ui.setProgramme(s.sourceDocument.media.frames,48000,"gain-source");lane.synchronize("gain-source");
    lane.selectAt(120000,lane.plotBounds().getCentreY());
    expect(!button(lane,"Apply gains")&&!button(lane,"Discard gains"),"manual lane has no Apply or Discard buttons");
    const double precise=-1.234567891234567;stageGain(lane,"-1.234567891234567");
    expect(!s.sectionTrimDraft&&s.currentStateDocument.sectionGains&&s.currentStateDocument.sectionGains->gainAt(120000)==precise&&s.currentStateNeedsRemeasurement()&&!s.currentState.analysis&&s.source.analysis==original,"Return applies exact local gain and clears only Current measurement authority");
    auto applied=s.currentStateDocument.sectionGains->audioKey();
    for(const auto* invalid:{"nan","inf","1 trailing","12.00001","-12.00001",""}){stageGain(lane,invalid);expect(s.currentStateDocument.sectionGains->audioKey()==applied,"invalid entry never replaces applied processing");laneEditor(lane)->onEscapeKey();}
    lane.selectFrame(192000,true);stageGain(lane,"3.0000000000000004");
    expect(s.currentStateDocument.sectionGains->gainAt(220000)==3.0000000000000004&&s.currentStateDocument.sectionGains->gainAt(120000)==precise,"next edit preserves both exact settings without an intervening measurement");
    expect(s.undoManual()&&s.currentStateDocument.sectionGains->gainAt(220000)==0&&s.currentStateDocument.sectionGains->gainAt(120000)==precise,"one undo restores only the last edit");
    expect(s.undoManual(true)&&s.currentStateDocument.sectionGains->gainAt(220000)==3.0000000000000004,"redo restores exact processing");
    measure(s,Domain::currentState);lane.synchronize("gain-source");const auto measured=s.currentState.analysis;const auto revision=s.revision;
    lane.selectFrame(120000,false);stageGain(lane,"-1.234567891234567");expect(s.currentState.analysis==measured&&s.revision==revision,"same numeric setting is a no-op preserving verification");
    s.auditionBusy=[]{return true;};lane.synchronize("gain-source");expect(!laneEditor(lane)->isEnabled()&&!s.undoManual(),"playback blocks editing and undo");s.auditionBusy={};
    s.job.begin(JobStage::analysis,s.revision);lane.synchronize("gain-source");expect(!laneEditor(lane)->isEnabled()&&!s.setManualSectionGains({document.map().sections[1].id},-2),"active jobs block direct edits");s.job.finish(s.job.ticket,s.revision,ErrorCode::none,"done");lane.synchronize("gain-source");
    lane.selectRange(1,192000);expect(lane.selectedCount()==1,"range selection excludes partially enclosed left section");stageGain(lane,"-3");expect(s.currentStateDocument.sectionGains->gainAt(100)==0&&s.currentStateDocument.sectionGains->gainAt(120000)==-3,"range edit affects only fully selected sections");
    lane.selectRange(0,192000);expect(lane.selectedCount()==2&&laneEditor(lane)->getTextToShowWhenEmpty()=="Mixed","heterogeneous selection displays Mixed");stageGain(lane,"-2");
    expect(s.currentStateDocument.sectionGains->gainAt(100)==-2&&s.currentStateDocument.sectionGains->gainAt(120000)==-2&&s.currentStateDocument.sectionGains->gainAt(220000)==3.0000000000000004,"batch value replaces each selected local gain and preserves outside passage");
    expect(s.undoManual()&&s.currentStateDocument.sectionGains->gainAt(100)==0&&s.currentStateDocument.sectionGains->gainAt(120000)==-3,"batch gain change is one Undo step");lane.synchronize("gain-source");
    lane.selectRange(10,20);lane.synchronize("gain-source");expect(lane.selectedCount()==0&&!laneEditor(lane)->isEnabled(),"empty marquee cannot edit the prior selection after a refresh");
    lane.selectAt(120000,lane.plotBounds().getBottom()+5);expect(laneEditor(lane)->isEnabled()&&lane.inspectionText().contains("Section"),"readout click retains section mode");
    auto* field=laneEditor(lane);field->setText("-2.5",false);field->onTextChange();field->onFocusLost();expect(s.currentStateDocument.sectionGains->gainAt(120000)==-3,"focus loss cancels unconfirmed gain");
    StandaloneSectionsPanel other(s);other.suspend();juce::TextEditor* otherField=nullptr;
    for(auto* c:other.getChildren())if(c->getName()=="Section gain dB")otherField=dynamic_cast<juce::TextEditor*>(c);
    expect(otherField&&otherField->isEnabled()&&!button(other,"Apply section trims")&&!button(other,"Discard draft"),"Sections window also uses direct entry without manual Apply/Discard");
    otherField->setText("1.25",false);otherField->onTextChange();otherField->onFocusLost();expect(s.currentStateDocument.sectionGains->gainAt(100)==0,"separate Sections editor cancels on focus loss");
    otherField->setText("1.25",false);otherField->onTextChange();otherField->onReturnKey();expect(s.currentStateDocument.sectionGains->gainAt(100)==1.25,"separate Sections editor confirms on Enter");
    lane.synchronize("gain-source");lane.selectFrame(0,false);expect(laneEditor(lane)->getText()=="1.25","main lane observes the separate window's committed change");
    imageOf(ui,root.getChildFile("direct-gain-lane.png"));
    const auto prior=s.currentStateDocument.sectionGains->audioKey();s.sessionId=juce::Uuid().toString();stageGain(lane,"2");expect(s.currentStateDocument.sectionGains->audioKey()==prior&&!s.canUndoManual(),"stale runtime neither commits an old field nor exposes another session's undo");
    s.sourceContentInvalid=true;lane.synchronize({});expect(!laneEditor(lane)->isEnabled(),"invalid source clears editable lane");s.sourceContentInvalid=false;
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
    s.discardSectionTrimDraft();const auto editStart=juce::Time::getMillisecondCounterHiRes();
    expect(s.setManualSectionGains({map.sections[5000].id},-2)&&s.canUndoManual(),"dense 10000-section direct edit is journalled without audio copies");
    auto* editPerf=new juce::DynamicObject();editPerf->setProperty("dense_10000_sections_commit_ms",juce::Time::getMillisecondCounterHiRes()-editStart);cases.add(juce::var(editPerf));
    expect(s.undoManual()&&!s.currentStateDocument.sectionGains,"dense edit undo restores unity exactly");
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

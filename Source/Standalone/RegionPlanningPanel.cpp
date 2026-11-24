#include "RegionPlanningPanel.h"
#include "ProposalPresentation.h"
#include <limits>

namespace levelscope::standalone
{
namespace
{
juce::String settingsKey(const StandaloneSessionModel& s){return s.regionPlanning?s.regionPlanning->digest():"M19";}
juce::String spanText(proposal::RegionSpan s){return "["+juce::String(s.begin)+", "+juce::String(s.end)+")";}
proposal::RegionConstraints geometry(const proposal::RegionPlanningSettings& s,const StandaloneSessionModel& session)
{
    const auto c=session.measurementContext(MeasurementDomain::currentState);
    processing::SectionGainPlan before;before.inputIdentity=s.inputIdentity;before.sampleRate=c.media.sampleRate;before.frames=c.media.frames;
    before.intervals={{0,before.frames,0}};if(c.processing.sectionGains)before=*c.processing.sectionGains;
    // Only coordinate validation, never a measured result or an automatic root.
    auto metadata=std::make_shared<proposal::ProgramAnalysis>();metadata->frames=before.frames;metadata->sampleRate=before.sampleRate;
    proposal::LraPassages p;p.analysis=metadata;
    return proposal::resolveRegionConstraints(s,session.sectionDocument.map(),p,before);
}
}
std::optional<int64_t> RegionPlanningPanel::exactFrame(juce::String text)
{
    const auto bytes=text.trim().toStdString();if(bytes.empty()||bytes.size()>19)return {};
    int64_t n=0;for(char c:bytes){if(c<'0'||c>'9'||n>(std::numeric_limits<int64_t>::max()-(c-'0'))/10)return {};n=n*10+(c-'0');}return n;
}
RegionFrameRows::RegionFrameRows(juce::String name,bool pins,size_t maximum):allowPins(pins),limit(maximum)
{
    setName(name);title.setText(name,juce::dontSendNotification);rows.setRowHeight(24);
    for(auto* c:std::initializer_list<juce::Component*>{&title,&rows,&first,&last,&pin,&add,&update,&remove,&status})addAndMakeVisible(c);
    first.setName(name+" minimum / begin");last.setName(name+" maximum / end");
    first.setTextToShowWhenEmpty("Minimum / begin frame",juce::Colours::grey);last.setTextToShowWhenEmpty("Maximum / end frame",juce::Colours::grey);
    first.setInputRestrictions(19);last.setInputRestrictions(19);pin.setVisible(pins);
    pin.onClick=[this]{last.setEnabled(!pin.getToggleState());};
    add.onClick=[this]{store(true);};update.onClick=[this]{store(false);};
    remove.onClick=[this]{const int n=rows.getSelectedRow();if(n<0||size_t(n)>=values.size())return;values.erase(values.begin()+n);load(values);if(onChanged)onChanged();};
    load({});
}
void RegionFrameRows::load(std::vector<proposal::RegionEdgeRange> v)
{values=std::move(v);rows.updateContent();rows.deselectAllRows();first.setText({},false);last.setText({},false);pin.setToggleState(false,juce::dontSendNotification);last.setEnabled(true);markClean();status.setText("Select a row to edit; Add / Update commits it to this private draft.",juce::dontSendNotification);}
void RegionFrameRows::markClean(){cleanFirst=first.getText().trim();cleanLast=last.getText().trim();cleanPin=pin.getToggleState();}
juce::String RegionFrameRows::pendingError()const
{
    // Inspect actual editor values synchronously. Deferred TextEditor callbacks
    // cannot turn an already-stored row or an internal reload into a new edit.
    const bool pending=first.getText().trim()!=cleanFirst||last.getText().trim()!=cleanLast||pin.getToggleState()!=cleanPin;
    return pending?"Finish the pending frame row with Add / Update, or Discard edits.":juce::String();
}
void RegionFrameRows::store(bool append)
{
    auto a=RegionPlanningPanel::exactFrame(first.getText()),b=allowPins&&pin.getToggleState()?a:RegionPlanningPanel::exactFrame(last.getText());
    const int n=rows.getSelectedRow();
    if(!a||!b||*b<*a||(!allowPins&&*b==*a)){status.setText("Use exact non-negative frames in order; protections must have positive length.",juce::dontSendNotification);return;}
    if(append&&values.size()>=limit){status.setText("Row limit: "+juce::String(limit)+". Edit or remove an existing row.",juce::dontSendNotification);return;}
    if(!append&&(n<0||size_t(n)>=values.size())){status.setText("Select a row to update.",juce::dontSendNotification);return;}
    if(append)values.push_back({*a,*b});else values[size_t(n)]={*a,*b};
    rows.updateContent();rows.selectRow(append?int(values.size())-1:n);markClean();status.setText("Private draft row saved; Confirm settings is still required.",juce::dontSendNotification);if(onChanged)onChanged();
}
void RegionFrameRows::selectedRowsChanged(int n)
{
    if(n<0||size_t(n)>=values.size())return;
    first.setText(juce::String(values[size_t(n)].minimum),false);last.setText(juce::String(values[size_t(n)].maximum),false);
    pin.setToggleState(allowPins&&values[size_t(n)].minimum==values[size_t(n)].maximum,juce::dontSendNotification);last.setEnabled(!pin.getToggleState());markClean();
}
void RegionFrameRows::paintListBoxItem(int n,juce::Graphics& g,int w,int h,bool selected)
{
    if(n<0||size_t(n)>=values.size())return;if(selected)g.fillAll(juce::Colour(0xff325965));g.setColour(juce::Colours::white);g.setFont(13.f);
    const auto r=values[size_t(n)];g.drawText(juce::String(n+1)+": "+juce::String(r.minimum)+" - "+juce::String(r.maximum)+(allowPins?(r.minimum==r.maximum?"  exact pin":"  closed range"):"  end exclusive"),4,0,w-8,h,juce::Justification::centredLeft);
}
void RegionFrameRows::resized()
{
    auto r=getLocalBounds();title.setBounds(r.removeFromTop(25));status.setBounds(r.removeFromBottom(24));auto edit=r.removeFromBottom(30);
    for(auto* b:{&remove,&update,&add})b->setBounds(edit.removeFromRight(100).reduced(2));
    if(allowPins)pin.setBounds(edit.removeFromRight(95));first.setBounds(edit.removeFromLeft(edit.getWidth()/2).reduced(2));last.setBounds(edit.reduced(2));rows.setBounds(r.reduced(0,2));
}
RegionPlanningPanel::RegionPlanningPanel(StandaloneSessionModel& s):session(s)
{
    setName("Experimental region planning - private settings draft");viewport.setViewedComponent(&body,false);addAndMakeVisible(viewport);
    heading.setText("Proposal settings  |  Draft settings do not change audio",juce::dontSendNotification);heading.setFont(18.f);
    methodLabel.setText("Proposal method",juce::dontSendNotification);modeLabel.setText("Planning area",juce::dontSendNotification);
    sectionLabel.setText("Original section (annotation remains unchanged)",juce::dontSendNotification);selectionLabel.setText("Temporary selection [begin, end) - exact audio-source frames",juce::dontSendNotification);
    edgeLabel.setText("Processing start/end includes ramps. Unchecked = full legal range. Equal bounds = exact pin.",juce::dontSendNotification);
    footnote.setText("Local protections do not prevent separately permitted global gain. Gain / reserve / affected-duration limits remain in Proposal limits.",juce::dontSendNotification);
    method.setName("Proposal method");method.addItem("Standard proposals (default)",1);method.addItem("Measured-passage refinement",2);method.addItem("Prefer complete-event treatment",3);method.addItem("Competing passage boundaries",4);
    mode.setName("Region mode");mode.addItem("Automatic - refine measured seeds inward",1);mode.addItem("Inside one selected section (+ optional margins)",2);mode.addItem("Temporary time selection (+ optional margins)",3);
    section.setName("Original planning section");
    const char* names[]={"Maximum splits per original root (0-4)","Left margin seconds (0-3, step 0.1)","Right margin seconds (0-3, step 0.1)","Minimum nominal interval seconds (0.5-10)","Minimum unchanged separation seconds (0-3)","Need-change threshold LU (0.5-3)","Need persistence seconds (0.3-3)"};
    for(size_t n=0;n<fields.size();++n){labels[n].setText(names[n],juce::dontSendNotification);fields[n].setName(names[n]);body.addAndMakeVisible(labels[n]);body.addAndMakeVisible(fields[n]);fields[n].setInputRestrictions(8);}
    for(auto* c:std::initializer_list<juce::Component*>{&heading,&methodLabel,&modeLabel,&sectionLabel,&selectionLabel,&edgeLabel,&footnote,&method,&mode,&section,&selectionBegin,&selectionEnd,&startRange,&startMin,&startMax,&endRange,&endMin,&endMax,&startPin,&endPin,&splitting,&refinement,&fillGap,&global,&compare,&automaticMargins,&eventNote,&cuts,&protection})body.addAndMakeVisible(c);
    for(auto* c:std::initializer_list<juce::Component*>{&status,&confirm,&assess,&discard,&cancel})addAndMakeVisible(c);
    status.setMultiLine(true);status.setReadOnly(true);status.setScrollbarsShown(true);status.setFont(13.f);status.setName("Region draft validation and effective permissions");
    const char* frames[]={"Temporary begin frame","Temporary end frame","Processing start minimum frame","Processing start maximum frame","Processing end minimum frame","Processing end maximum frame"};
    int n=0;for(auto* f:{&selectionBegin,&selectionEnd,&startMin,&startMax,&endMin,&endMax}){f->setName(frames[n++]);f->setInputRestrictions(19);f->setTextToShowWhenEmpty("Exact source frame",juce::Colours::grey);}
    automaticMargins.setName("Allow automatic outward margins");
    eventNote.setFont(13.f);eventNote.setMinimumHorizontalScale(1.f);
    automaticMargins.onClick=[this]{if(!automaticMargins.getToggleState()){fields[1].setText("0.0",false);fields[2].setText("0.0",false);}feedback="Automatic margin permission changed in this private draft. Context reading alone grants no processing permission.";updateEnabled();};
    method.onChange=[this]
    {
        if(method.getSelectedId()>=3)compare.setToggleState(true,juce::dontSendNotification);
        else if(method.getSelectedId()==2)
        {automaticMargins.setToggleState(false,juce::dontSendNotification);if(mode.getSelectedId()==1){fields[1].setText("0.0",false);fields[2].setText("0.0",false);}}
        feedback="Planning approach changed in this private draft. Review permissions, then Confirm settings.";updateEnabled();
    };
    mode.onChange=[this]
    {
        // A mode change is a visible private draft reset of mode-specific permissions.
        global.setToggleState(mode.getSelectedId()==1,juce::dontSendNotification);automaticMargins.setToggleState(false,juce::dontSendNotification);
        startRange.setToggleState(false,juce::dontSendNotification);endRange.setToggleState(false,juce::dontSendNotification);cuts.load({});
        fields[1].setText("0.0",false);fields[2].setText("0.0",false);feedback="Mode changed in the private draft; margins and edge/cut restrictions reset. Review before confirming.";updateEnabled();
    };
    for(auto* t:{&startRange,&endRange,&startPin,&endPin,&splitting,&refinement})t->onClick=[this]{feedback.clear();updateEnabled();};
    confirm.onClick=[this]{commit(false);};assess.onClick=[this]{commit(true);};discard.onClick=[this]{load();feedback="Private edits discarded; active method and audio unchanged.";timerCallback();};
    cancel.onClick=[this]{load();if(auto* w=findParentComponentOfClass<ProposalWindow>())w->setVisible(false);};
    load();setSize(1060,860);startTimerHz(4);timerCallback();
}
void RegionPlanningPanel::showWindow(std::unique_ptr<ProposalWindow>& w,StandaloneSessionModel& s,juce::Component* owner,std::function<void()> changed,std::function<void()> assessCallback)
{
    auto* p=w?dynamic_cast<RegionPlanningPanel*>(w->getContentComponent()):nullptr;
    if(!p||p->runtime!=s.sessionId||!(w->isVisible()||w->isMinimised()))
    {auto content=std::make_unique<RegionPlanningPanel>(s);content->onChanged=std::move(changed);content->onAssess=std::move(assessCallback);w=std::make_unique<ProposalWindow>("Proposal settings - method and permissions",std::move(content),owner);}
    w->present();
}
void RegionPlanningPanel::load()
{
    runtime=session.sessionId;input=sectionInputIdentity(session.measurementContext(MeasurementDomain::source));baseDigest=settingsKey(session);annotations=session.sectionDocument.map();
    auto s=session.regionPlanning.value_or(proposal::RegionPlanningSettings::defaults(proposal::RegionMode::automatic,input));
    loadedSettings=s;automaticMargins.setToggleState(s.allowAutomaticMargins,juce::dontSendNotification);
    method.setSelectedId(session.regionPlanning?s.schema+1:1,juce::dontSendNotification);mode.setSelectedId(int(s.mode)+1,juce::dontSendNotification);
    section.clear(juce::dontSendNotification);for(size_t n=0;n<annotations.sections.size();++n){const auto& r=annotations.sections[n];section.addItem(r.name+" "+spanText({r.begin,r.end}),int(n+1));if(r.id==s.sectionId)section.setSelectedId(int(n+1),juce::dontSendNotification);}
    if(section.getSelectedId()==0&&!annotations.sections.empty())section.setSelectedId(1,juce::dontSendNotification);
    selectionBegin.setText(juce::String(s.selection.begin),false);selectionEnd.setText(juce::String(s.selection.end),false);
    const int values[]={s.maximumSplits,s.leftMarginMs,s.rightMarginMs,s.minimumIntervalMs,s.separationMs,s.needThresholdTenths,s.persistenceMs};
    for(size_t n=0;n<fields.size();++n)fields[n].setText(n==0?juce::String(values[n]):juce::String(values[n]/(n==5?10.:1000.),1),false);
    splitting.setToggleState(s.allowSplitting,juce::dontSendNotification);refinement.setToggleState(s.allowRefinement,juce::dontSendNotification);fillGap.setToggleState(s.fillShortGap,juce::dontSendNotification);global.setToggleState(s.allowGlobal,juce::dontSendNotification);compare.setToggleState(s.compareAlternatives,juce::dontSendNotification);
    const auto edge=[](const auto& r,auto& enabled,auto& pin,auto& a,auto& b){enabled.setToggleState(bool(r),juce::dontSendNotification);pin.setToggleState(r&&r->minimum==r->maximum,juce::dontSendNotification);a.setText(r?juce::String(r->minimum):juce::String(),false);b.setText(r?juce::String(r->maximum):juce::String(),false);};
    edge(s.start,startRange,startPin,startMin,startMax);edge(s.end,endRange,endPin,endMin,endMax);cuts.load(s.internalCuts);
    std::vector<proposal::RegionEdgeRange> protectedRows;for(auto p:s.protectedLocal)protectedRows.push_back({p.begin,p.end});protection.load(std::move(protectedRows));feedback.clear();updateEnabled();
}
juce::String RegionPlanningPanel::readDraft(std::optional<proposal::RegionPlanningSettings>& out)const
{
    out.reset();if(runtime!=session.sessionId||input!=sectionInputIdentity(session.measurementContext(MeasurementDomain::source)))return "Programme/session changed. Discard edits / reload before confirming.";
    if(baseDigest!=settingsKey(session))return "Active method/settings changed. Discard edits / reload before confirming.";
    if(method.getSelectedId()==1)return {};
    if(!session.sourceDocument.hasSource||!session.sourceDocument.layout.confirmed||session.sourceContentInvalid)return "Import and confirm programme roles first.";
    if(auto error=cuts.pendingError();error.isNotEmpty())return error;if(auto error=protection.pendingError();error.isNotEmpty())return error;
    const auto schema=method.getSelectedId()-1;const auto area=proposal::RegionMode(mode.getSelectedId()-1);
    if(schema<1||schema>3)return "Select a supported planning approach.";
    auto s=schema==loadedSettings.schema?loadedSettings:schema==3?proposal::RegionPlanningSettings::innerDefaults(area,input):schema==2?proposal::RegionPlanningSettings::eventDefaults(area,input):proposal::RegionPlanningSettings::defaults(area,input);
    s.mode=area;s.inputIdentity=input;
    if(s.schema>=2)s.allowAutomaticMargins=s.mode==proposal::RegionMode::automatic&&automaticMargins.getToggleState();
    if(s.mode==proposal::RegionMode::selectedSection)
    {
        const int n=section.getSelectedId()-1;if(n<0||size_t(n)>=annotations.sections.size())return "Select an original section.";
        const auto& selected=annotations.sections[size_t(n)];s.sectionId=selected.id;
        const auto& now=session.sectionDocument.map().sections;
        if(std::none_of(now.begin(),now.end(),[&](const auto& r){return r.id==selected.id&&r.begin==selected.begin&&r.end==selected.end;}))return "Selected section bounds changed. Discard edits / reload before confirming.";
    }
    if(s.mode==proposal::RegionMode::selection)
    {auto a=exactFrame(selectionBegin.getText()),b=exactFrame(selectionEnd.getText());if(!a||!b)return "Use exact integer source frames for the temporary selection.";s.selection={*a,*b};}
    int* values[]={&s.maximumSplits,&s.leftMarginMs,&s.rightMarginMs,&s.minimumIntervalMs,&s.separationMs,&s.needThresholdTenths,&s.persistenceMs};
    const int minimum[]={0,0,0,5,0,5,3},maximum[]={4,30,30,100,30,30,30};
    for(size_t n=0;n<fields.size();++n)
    {
        auto v=ProposalPreferencesPanel::exactUnits(fields[n].getText(),n!=0);
        if(!v)return "Use integer splits and decimal controls with at most one fractional digit. Values are never rounded.";
        if(*v<minimum[n]||*v>maximum[n])return labels[n].getText()+": value outside the approved range.";
        *values[n]=*v*(n==0||n==5?1:100);
    }
    s.allowSplitting=splitting.getToggleState();s.allowRefinement=refinement.getToggleState();s.fillShortGap=fillGap.getToggleState();s.allowGlobal=global.getToggleState();s.compareAlternatives=s.schema>=2||compare.getToggleState();
    const auto edge=[&](bool enabled,bool pin,const auto& a,const auto& b,std::optional<proposal::RegionEdgeRange>& r)->bool
    {if(!enabled)return true;auto x=exactFrame(a.getText()),y=pin?x:exactFrame(b.getText());if(!x||!y)return false;r=proposal::RegionEdgeRange{*x,*y};return true;};
    s.start.reset();s.end.reset();
    if(!edge(startRange.getToggleState(),startPin.getToggleState(),startMin,startMax,s.start)||!edge(endRange.getToggleState(),endPin.getToggleState(),endMin,endMax,s.end))return "Use exact integer source frames for edge ranges/pins.";
    s.internalCuts=cuts.values;s.protectedLocal.clear();for(auto p:protection.values)s.protectedLocal.push_back({p.minimum,p.maximum});
    auto error=s.validationError();if(error.isNotEmpty())return error;
    const auto resolved=geometry(s,session);if(resolved.error.isNotEmpty())return resolved.error;out=std::move(s);return {};
}
juce::String RegionPlanningPanel::geometryText(const proposal::RegionPlanningSettings& s)const
{
    const auto c=geometry(s,session);juce::String text;
    if(s.schema>=2)text=proposalPlanningName(s)+" (experimental). Comparative assessment is required; partial/uncertain passes need explicit review and selection. "+juce::String(s.eventFeatures.hopMs)+" ms feature cadence is evaluation spacing, not proven onset accuracy.\n";
    text+=proposalDefinitionText(s)+"\n";
    if(s.mode==proposal::RegionMode::automatic)text+=s.allowAutomaticMargins?
        "Automatic margins explicitly permitted: left "+juce::String(s.leftMarginMs/1000.,1)+" s, right "+juce::String(s.rightMarginMs/1000.,1)+" s. Assessment freezes each root's ownership and shows requested versus effective clipped areas/ranges. Context grants no permission.\n":
        "Automatic: measured roots will be resolved by assessment; no outward movement beyond each seed. Protections exclude additional local change.\n";
    for(const auto& r:c.roots)
    {
        const int64_t left=int64_t(s.leftMarginMs)*c.sampleRate/1000,right=int64_t(s.rightMarginMs)*c.sampleRate/1000;
        text+="Hint "+spanText(r.hint)+"; requested area "+spanText({r.hint.begin-left,r.hint.end+right})+"; effective clipped area "+spanText(r.allowed)+".\n";
        text+="Effective processing-start range ["+juce::String(r.start.minimum)+", "+juce::String(r.start.maximum)+"]; end ["+juce::String(r.end.minimum)+", "+juce::String(r.end.maximum)+"]. Exact pins are not rounded.\n";
    }
    return text+(s.allowGlobal?"Global trim permitted: the whole programme, including protected/local-excluded audio, may change globally.":"Global gain fixed. Local protection preserves the existing local envelope.");
}
void RegionPlanningPanel::commit(bool after)
{
    if(session.audioEditsBlocked()||session.job.isActive()||session.persistenceLoading){feedback="Stop playback and finish/cancel the current job before confirming.";return;}
    std::optional<proposal::RegionPlanningSettings> value;auto error=readDraft(value);
    if(error.isNotEmpty()){feedback=error;return;}
    if(!session.setRegionPlanning(value,error)){feedback=error;return;}baseDigest=settingsKey(session);
    feedback="Settings confirmed. Compatible measurements/audio unchanged. Previous proposals may be stale; Apply remains separate.";
    if(onChanged)onChanged();if(after&&onAssess)onAssess();timerCallback();
}
void RegionPlanningPanel::updateEnabled()
{
    const bool experimental=method.getSelectedId()>=2,selected=mode.getSelectedId()!=1,event=method.getSelectedId()>=3;
    automaticMargins.setEnabled(experimental&&event&&!selected);
    const auto selectedMode=mode.getSelectedId();
    mode.changeItemText(1,event?"Automatic - measured seeds with explicit permissions":"Automatic - refine measured seeds inward");
    // JUCE does not refresh the displayed text in changeItemText. Keep both
    // text and ID aligned without issuing a mode-change permission reset.
    mode.setSelectedId(selectedMode,juce::dontSendNotification);
    eventNote.setText(method.getSelectedId()==4?"Compares a continuing passage with a possible inner beginning and earlier/later returns. Unresolved evidence stays visible; a measured pass may still need explicit selection.":event?"Prefers supported-complete measured passes. Partial or uncertain coverage needs explicit selection; this is not phrase recognition.":"Standard proposals retain all supported jobs. Boundary approaches require measured IL/true peak/LRA and confirmed mono (Left, Right or Centre), stereo or conventional 5.1.",juce::dontSendNotification);
    mode.setEnabled(experimental);section.setEnabled(experimental&&mode.getSelectedId()==2);
    selectionBegin.setEnabled(experimental&&mode.getSelectedId()==3);selectionEnd.setEnabled(selectionBegin.isEnabled());
    for(size_t n=0;n<fields.size();++n)fields[n].setEnabled(experimental&&((n!=1&&n!=2)||selected||(event&&automaticMargins.getToggleState())));
    for(auto* t:{&splitting,&refinement,&fillGap,&global})t->setEnabled(experimental);compare.setEnabled(experimental&&!event);
    startRange.setEnabled(experimental&&selected);endRange.setEnabled(experimental&&selected);
    startPin.setEnabled(startRange.isEnabled()&&startRange.getToggleState());endPin.setEnabled(endRange.isEnabled()&&endRange.getToggleState());
    startMin.setEnabled(startPin.isEnabled());startMax.setEnabled(startPin.isEnabled()&&!startPin.getToggleState());endMin.setEnabled(endPin.isEnabled());endMax.setEnabled(endPin.isEnabled()&&!endPin.getToggleState());
    cuts.setEnabled(experimental&&selected);protection.setEnabled(experimental);
    timerCallback();
}
void RegionPlanningPanel::timerCallback()
{
    const bool idle=!session.audioEditsBlocked()&&!session.job.isActive()&&!session.persistenceLoading;
    confirm.setEnabled(idle);assess.setEnabled(idle&&session.source.hasMeasurement&&session.currentState.hasMeasurement&&!session.currentStateNeedsRemeasurement()&&!session.sectionTrimDraft);
    std::optional<proposal::RegionPlanningSettings> draft;const auto error=readDraft(draft);
    juce::String text=idle?"Private draft. Confirming settings changes no audio; assessment and Apply are separate.\n":"Inspection/draft only. Stop playback and finish/cancel the job to confirm.\n";
    if(feedback.isNotEmpty())text+=feedback+"\n";
    text+=error.isNotEmpty()?"Draft unavailable: "+error:draft?geometryText(*draft):"Standard proposals: existing global-first and section/passage paths for all supported jobs. No boundary-method opt-in.";
    if(!assess.isEnabled())text+="\nAssessment needs idle transport, matching measured Source/Current and no manual gain draft.";
    if(status.getText()!=text)status.setText(text,false);
}
void RegionPlanningPanel::resized()
{
    auto r=getLocalBounds().reduced(12);auto actions=r.removeFromBottom(34);
    confirm.setBounds(actions.removeFromLeft(175).reduced(2));assess.setBounds(actions.removeFromLeft(190).reduced(2));discard.setBounds(actions.removeFromLeft(200).reduced(2));cancel.setBounds(actions.removeFromLeft(100).reduced(2));
    status.setBounds(r.removeFromBottom(140).reduced(0,4));viewport.setBounds(r);body.setSize(std::max(760,viewport.getWidth()-18),1220);
    auto b=body.getLocalBounds().reduced(4);heading.setBounds(b.removeFromTop(32));
    auto row=b.removeFromTop(34);methodLabel.setBounds(row.removeFromLeft(180));method.setBounds(row.reduced(2));
    eventNote.setBounds(b.removeFromTop(50));automaticMargins.setBounds(b.removeFromTop(28));
    row=b.removeFromTop(34);modeLabel.setBounds(row.removeFromLeft(180));mode.setBounds(row.reduced(2));
    sectionLabel.setBounds(b.removeFromTop(24));section.setBounds(b.removeFromTop(30).reduced(2));
    selectionLabel.setBounds(b.removeFromTop(24));row=b.removeFromTop(30);selectionBegin.setBounds(row.removeFromLeft(row.getWidth()/2).reduced(2));selectionEnd.setBounds(row.reduced(2));
    row=b.removeFromTop(30);splitting.setBounds(row.removeFromLeft(row.getWidth()/2));refinement.setBounds(row);
    for(size_t n=0;n<fields.size();++n){row=b.removeFromTop(30);fields[n].setBounds(row.removeFromRight(140).reduced(2));labels[n].setBounds(row);}
    fillGap.setBounds(b.removeFromTop(28));global.setBounds(b.removeFromTop(28));compare.setBounds(b.removeFromTop(28));
    edgeLabel.setBounds(b.removeFromTop(28));
    const auto edge=[&](auto& label,auto& pin,auto& a,auto& z){auto line=b.removeFromTop(32);label.setBounds(line.removeFromLeft(215));pin.setBounds(line.removeFromLeft(100));a.setBounds(line.removeFromLeft(line.getWidth()/2).reduced(2));z.setBounds(line.reduced(2));};
    edge(startRange,startPin,startMin,startMax);edge(endRange,endPin,endMin,endMax);
    cuts.setBounds(b.removeFromTop(150));protection.setBounds(b.removeFromTop(185));footnote.setBounds(b.removeFromTop(42));
}
}

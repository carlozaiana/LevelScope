#include "TimelineGainLane.h"
#include "TimelineGainOutline.h"
#include "TimelineAppearance.h"
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>

namespace levelscope::standalone
{
namespace
{
juce::String exact(double value)
{
    std::ostringstream s;s.imbue(std::locale::classic());s<<std::setprecision(std::numeric_limits<double>::max_digits10)<<value;return s.str();
}
template<class Rows> size_t intervalAt(const Rows& rows,int64_t frame)
{
    return size_t(std::upper_bound(rows.begin(),rows.end(),frame,[](int64_t f,const auto& r){return f<r.end;})-rows.begin());
}
}
TimelineGainLane::TimelineGainLane(StandaloneSessionModel& s):session(s)
{
    setName("Current section gain lane");setInterceptsMouseClicks(false,true);setBufferedToImage(true);
    for(auto* c:std::initializer_list<juce::Component*>{&legend,&notice,&detail,&previous,&next,&gain,&dbLabel,&stage,&apply,&discard,&undoDelete})addAndMakeVisible(c);
    for(auto* label:{&legend,&notice,&detail,&dbLabel}){label->setFont(13.f);label->setInterceptsMouseClicks(false,false);}
    legend.setText("CURRENT section gain (dB)   |   Solid: applied   /   Dashed: draft",juce::dontSendNotification);
    legend.setTooltip("All programme channels, before global trim. This lane always describes Current processing.");
    dbLabel.setText("Draft dB",juce::dontSendNotification);gain.setName("Timeline section draft gain dB");
    gain.setTooltip("Enter -12 to +12 dB, then Return or Set draft. Apply gains changes Current.");
    gain.setInputRestrictions(36);gain.setSelectAllWhenFocused(true);
    gain.onTextChange=[this]{dirty=true;feedback.clear();updateControls();};
    gain.onReturnKey=[this]{stageValue();};
    gain.onEscapeKey=[this]{dirty=false;feedback.clear();showSelection();updateControls();};
    stage.onClick=[this]{stageValue();};
    undoDelete.onClick=[this]
    {
        if(!deletedAnnotation||!bindingMatches()||session.sectionDocument.map().revision!=deletedRevision
            ||session.audioEditsBlocked()||session.job.isActive()||session.sectionTrimDraft)return;
        auto restored=*deletedAnnotation;restored.revision=deletedRevision+1;
        if(session.sectionDocument.restore(restored,input)){session.documentChanged();deletedAnnotation.reset();afterAction();}
    };
    apply.onClick=[this]
    {
        if(!bindingMatches()||dirty||!apply.isEnabled())return;
        session.applySectionTrimDraft();afterAction();
    };
    discard.onClick=[this]
    {
        // Discard is also valid when the file disappeared or a draft became
        // stale. It must still belong to this opening and displayed draft.
        if(runtime!=session.sessionId||audioRevision!=session.revision||hasDraft!=bool(session.sectionTrimDraft)
            ||(hasDraft&&editSerial!=session.sectionTrimDraft->editSerial)||session.audioEditsBlocked()||session.job.isActive()||!discard.isEnabled())return;
        if(session.sectionTrimDraft){session.discardSectionTrimDraft();afterAction();}
        else {dirty=false;showSelection();feedback="Field edit discarded. Applied gains unchanged.";updateControls();}
    };
    previous.setTooltip("Inspect previous section");next.setTooltip("Inspect next section");
    previous.onClick=[this]{adjacent(-1);};
    next.onClick=[this]{adjacent(1);};
    addAndMakeVisible(selectionKind);selectionKind.setName("Gain lane selection");
    selectionKind.addItem("Select sections",1);selectionKind.addItem("Select applied regions",2);selectionKind.addItem("Select proposed regions",3);
    selectionKind.setSelectedId(1,juce::dontSendNotification);
    selectionKind.setTooltip("Sections are editable annotations. Processing regions inspect fixed gain intervals, including after reopening. Click a gain readout below the lane to select its region.");
    selectionKind.onChange=[this]{deleteArmed=false;rebuildProposalRows();selectFrame(cursor,true);};
}
const std::vector<proposal::ProgrammeSection>& TimelineGainLane::rows() const
{return !lraRows.empty()?lraRows:session.sectionDocument.map().sections;}
const RegionAttempt* TimelineGainLane::inspectedAttempt() const
{
    const auto review=automatic?automatic->regionReview():nullptr;
    return review&&regionAttempt<review->attempts.size()?&review->attempts[regionAttempt]:nullptr;
}
bool TimelineGainLane::diagnosticRegionInspection() const
{
    const auto review=automatic?automatic->regionReview():nullptr;
    return review&&regionAttempt<review->attempts.size()&&(regionAttempt!=review->retainedAttempt||review->attempts[regionAttempt].evaluation.verdict!=proposal::Verdict::pass);
}
void TimelineGainLane::rebuildProposalRows()
{
    lraRows.clear();
    const bool proposal=selectionKind.getSelectedId()==3;
    const auto* plan=proposal?(proposed&&proposed->sectionGains?&*proposed->sectionGains:nullptr):selectionKind.getSelectedId()==2&&applied?&*applied:nullptr;
    if(plan){uint64_t id=1;for(const auto& r:plan->intervals)lraRows.push_back({id++,r.begin,r.end,proposal?"Proposed processing region":"Applied processing region"});}
}
void TimelineGainLane::inspectRegionAttempt(size_t n)
{
    const auto review=automatic?automatic->regionReview():nullptr;if(!review||n>=review->attempts.size())return;
    regionAttempt=n;const auto& plan=review->attempts[n].plan;proposed=ProcessingPlan{};proposed->staticGainDb=plan.globalDb;proposed->sectionGains=plan.envelope;
    selectionKind.setSelectedId(3,juce::dontSendNotification);
    rebuildProposalRows();row=intervalAt(rows(),cursor);showSelection();updateControls();repaint();
}
void TimelineGainLane::adjacent(int direction)
{
    if(!bindingMatches()||!selected())return;
    if(selectionKind.getSelectedId()==1&&automatic&&!automatic->selection().corrections.empty())
    {
        const auto& changes=automatic->selection().corrections;
        if(direction>0){for(const auto& c:changes)if(c.begin>cursor){selectFrame(c.begin,true);return;}}
        else {for(auto i=changes.rbegin();i!=changes.rend();++i)if(i->begin<cursor){selectFrame(i->begin,true);return;}}
    }
    else if(direction<0&&row>0)selectFrame(rows()[row-1].begin,true);
    else if(direction>0&&row+1<rows().size())selectFrame(rows()[row+1].begin,true);
}
const proposal::ProgrammeSection* TimelineGainLane::selected() const
{
    const auto& rows=this->rows();
    return usableMap&&row<rows.size()?&rows[row]:nullptr;
}
bool TimelineGainLane::bindingMatches() const
{
    const auto& d=session.sectionTrimDraft;
    return programme.isNotEmpty()&&session.sourceDocument.hasSource&&!session.sourceContentInvalid&&session.sourceDocument.matchesFileStamps()
        &&runtime==session.sessionId&&audioRevision==session.revision&&mapRevision==session.sectionDocument.map().revision
        &&initialized==session.hasCurrentStateInitialized()&&hasDraft==bool(d)&&(!d||editSerial==d->editSerial);
}
void TimelineGainLane::synchronize(const juce::String& key)
{
    const bool replaced=key!=programme||runtime!=session.sessionId;
    if(replaced)
    {
        deleteArmed=false;deletedAnnotation.reset();
        programme=key;runtime=session.sessionId;cursor=0;
        input.clear();
        if(key.isNotEmpty()&&session.sourceDocument.layout.confirmed)
        {
            const auto c=session.measurementContext(MeasurementDomain::source);
            if(c.media.frames<=int64_t(c.media.sampleRate)*86400)input=sectionInputIdentity(c);
            if(input.isNotEmpty())session.synchronizeSections();
        }
    }
    const auto& map=session.sectionDocument.map();const auto& d=session.sectionTrimDraft;
    if(deletedAnnotation&&map.revision!=deletedRevision)deletedAnnotation.reset();
    const bool stale=session.sectionProposalStale();
    if(replaced||automatic!=session.sectionProposal||automaticStale!=stale||audioRevision!=session.revision||mapRevision!=map.revision||initialized!=session.hasCurrentStateInitialized()
        ||hasDraft!=bool(d)||(d&&editSerial!=d->editSerial))
    {
        const bool newProposal=automatic!=session.sectionProposal;
        deleteArmed=false;
        automatic=session.sectionProposal;automaticStale=stale;proposed=automatic?automatic->previewPlan():std::nullopt;
        if(newProposal)regionAttempt=automatic&&automatic->regionReview()?automatic->regionReview()->retainedAttempt:std::numeric_limits<size_t>::max();
        if(const auto* attempt=inspectedAttempt()){proposed=ProcessingPlan{};proposed->staticGainDb=attempt->plan.globalDb;proposed->sectionGains=attempt->plan.envelope;}
        audioRevision=session.revision;mapRevision=map.revision;initialized=session.hasCurrentStateInitialized();
        usableMap=input.isNotEmpty()&&map.inputIdentity==input&&map.valid();
        hasDraft=bool(d);editSerial=d?d->editSerial:0;fresh=usableMap&&session.sectionDraftFresh();
        applied=usableMap&&initialized?session.currentStateDocument.sectionGains:std::nullopt;
        if(newProposal&&proposed&&proposed->sectionGains)selectionKind.setSelectedId(3,juce::dontSendNotification);
        if(selectionKind.getSelectedId()==3&&(!proposed||!proposed->sectionGains))selectionKind.setSelectedId(applied?2:1,juce::dontSendNotification);
        if(selectionKind.getSelectedId()==2&&!applied)selectionKind.setSelectedId(1,juce::dontSendNotification);
        selectionKind.setItemEnabled(2,bool(applied));selectionKind.setItemEnabled(3,proposed&&proposed->sectionGains.has_value());
        rebuildProposalRows();
        geometryMatches=!applied||session.sectionGeometry(map)==session.sectionGeometry(session.currentStateDocument.sectionGainAnnotations);
        supportError=usableMap&&initialized?session.sectionDraftSupportError():juce::String{};
        draftError=fresh?d->plan.validationError():juce::String{};
        draft=fresh&&draftError.isEmpty()?std::optional<processing::SectionGainPlan>(d->plan.canonical()):std::nullopt;
        cursor=std::clamp<int64_t>(cursor,0,std::max<int64_t>(0,map.frames-1));row=intervalAt(rows(),cursor);
        dirty=false;feedback.clear();showSelection();resized();repaint();
    }
    updateControls();
}
void TimelineGainLane::showSelection()
{
    const auto* section=selected();
    if(!section){detail.setText("No confirmed sections",juce::dontSendNotification);detail.setTooltip({});gain.setText({},false);return;}
    double stored=0;
    if(fresh&&session.sectionTrimDraft->processingRegions)stored=session.sectionTrimDraft->plan.intervals[intervalAt(session.sectionTrimDraft->plan.intervals,cursor)].db;
    else if(fresh&&selectionKind.getSelectedId()==1)stored=session.sectionTrimDraft->plan.intervals[row].db;
    else if(applied&&(geometryMatches||selectionKind.getSelectedId()==2))stored=applied->intervals[intervalAt(applied->intervals,cursor)].db;
    gain.setText(exact(stored),false);
    auto description=section->name+"   ["+juce::String(section->begin)+", "+juce::String(section->end)+")";
    if(initialized)description+="   Applied: "+exact(applied?applied->gainAt(cursor):0.)+" dB";
    if(selectionKind.getSelectedId()!=1)description+=" (at cursor; nominal interval shown, transitions may cross its edges).";
    if(fresh&&session.sectionTrimDraft->processingRegions)description+=" | removal draft local "+exact(draft?draft->gainAt(cursor):0.)+" dB";
    if(initialized)description+=" | applied global "+exact(session.currentStateDocument.staticGainDb)+" dB";
    if(proposed)
    {
        const double local=proposed->sectionGains?proposed->sectionGains->gainAt(cursor):0.;
        description+="   Proposed local "+exact(local)+" dB / global "+exact(proposed->staticGainDb)+" dB / combined "+exact(local+proposed->staticGainDb)+" dB";
        if(!automatic->previewHasFinalGlobal())description+=" (rejected local selection with baseline global gain; combined candidate not measured)";
        if(const auto* attempt=inspectedAttempt())
        {
            const auto& review=*automatic->regionReview();
            description="Measured attempt "+juce::String(regionAttempt+1)+(diagnosticRegionInspection()?" — DIAGNOSTIC ONLY. ":" — retained recommendation. ")+description;
            for(const auto& a:attempt->plan.assignments)if(a.geometry.support.begin<=cursor&&cursor<a.geometry.support.end)
            {
                description+="\nRoot "+(a.geometry.root<review.constraints.roots.size()?review.constraints.roots[a.geometry.root].id:"unavailable")+". "+a.geometry.reason;
                description+="\nChanged-support outer edges ["+juce::String(a.geometry.support.begin)+", "+juce::String(a.geometry.support.end)+"); nominal plateau cuts "+juce::String(a.nominal.begin);
                for(auto cut:a.geometry.cuts)description+=", "+juce::String(cut);description+=", "+juce::String(a.nominal.end)+". All values are exact source frames; ramps are included only in changed support.";
                const size_t part=size_t(std::upper_bound(a.geometry.cuts.begin(),a.geometry.cuts.end(),cursor)-a.geometry.cuts.begin());
                if(part<a.reductions.size())description+="\nAdditional plateau attenuation "+exact(a.reductions[part])+" dB";
                if(part<a.requested.size()&&part<a.reserve.size()&&part<a.uncovered.size())description+="; requested / reserve achieved / deficit "+exact(a.requested[part])+" / "+exact(a.reserve[part])+" / "+exact(a.uncovered[part])+" dB.";
            }
            description+="\nActual private result: "+evaluationText(attempt->evaluation)+". Not Current verification. Complete attempts, support and transition coordinates are in Review.";
        }
        if(automatic->lraSelection())for(const auto& c:automatic->lraSelection()->selected)if(c.passage.begin<=cursor&&cursor<c.passage.end)
        {
            description+="\nLRA region support ["+juce::String(c.passage.begin)+", "+juce::String(c.passage.end)+") frames; score "+exact(c.passage.score)+" LU*windows. Estimated requested/reduction "+exact(c.requested)+" / "+exact(c.reduction)+" dB; reserve "+exact(c.reserve)+", uncovered overshoot "+exact(c.uncovered)+". See Review for actual complete measured results.";break;
        }
        description=automaticStale?"STALE preview. "+description:automatic->status()==SectionProposalStatus::noSolution?"REJECTED preview. "+description:description;
        for(const auto& c:automatic->selection().corrections)if(c.begin==section->begin&&c.id==section->id)
        {
            description+="\nReason — channel "+juce::String(c.channel+1)+" TP "+juce::String(c.peakDbtp,9)+" dBTP at "+juce::String(c.peakEighthFrame)+"/8 frames. Original local "+exact(c.originalDb)+" dB -> proposed plateau "+exact(c.resultingDb)+" dB. Requested reduction "+exact(c.requestedDb)+", planned "+exact(c.reductionDb)+" dB. Reserve achieved "+exact(c.reserveDb)+" / "+juce::String(automatic->selection().policy.reserveTenths/10.,1)+" dB. Uncovered deficit "+exact(c.uncoveredDb)+" dB.";break;
        }
    }
    detail.setText(description.upToFirstOccurrenceOf("\n",false,false),juce::dontSendNotification);
    detail.setTooltip(description+" at frame "+juce::String(cursor)+(fresh&&session.sectionTrimDraft->processingRegions
        ?". The disabled field shows the draft region plateau; the applied readout includes transitions."
        :". The field shows the section setting; the applied readout includes transitions."));
}
void TimelineGainLane::updateControls()
{
    if(measured!=session.currentState.hasMeasurement){measured=session.currentState.hasMeasurement;feedback.clear();}
    const bool idle=!session.audioEditsBlocked()&&!session.job.isActive();
    const bool editable=selectionKind.getSelectedId()==1&&selected()&&initialized&&supportError.isEmpty()&&idle&&!session.sectionProposal&&(!hasDraft||(fresh&&!session.sectionTrimDraft->processingRegions));
    previous.setEnabled(selected()&&row>0);next.setEnabled(selected()&&row+1<rows().size());
    if(selectionKind.getSelectedId()==1&&automatic&&!automatic->selection().corrections.empty())
    {
        const auto& changes=automatic->selection().corrections;
        previous.setEnabled(selected()&&cursor>changes.front().begin);next.setEnabled(selected()&&cursor<changes.back().begin);
    }
    previous.setTooltip("Inspect previous section or processing region");next.setTooltip("Inspect next section or processing region");
    for(auto* c:std::initializer_list<juce::Component*>{&gain,&dbLabel,&stage,&apply,&discard})c->setVisible(!automatic);
    gain.setEnabled(editable);stage.setEnabled(editable);
    apply.setEnabled(initialized&&idle&&!session.sectionProposal&&supportError.isEmpty()&&fresh&&!dirty&&draftError.isEmpty());discard.setEnabled((hasDraft||dirty)&&idle);
    undoDelete.setVisible(deletedAnnotation.has_value());undoDelete.setEnabled(deletedAnnotation&&idle&&!hasDraft&&bindingMatches());
    juce::String message;
    if(!usableMap)message="Confirm programme audio and channel roles to inspect sections.";
    else if(!initialized)message="Current is not initialised.";
    else if(hasDraft&&!fresh)message="Stale draft. Discard it and draft on the current sections.";
    else if(supportError.isNotEmpty())message=supportError;
    else if(!idle)message="Inspection available. Stop playback or finish the job to edit.";
    else if(proposed)message=(automaticStale?"STALE. ":diagnosticRegionInspection()?"DIAGNOSTIC ATTEMPT — Apply unavailable. ":automatic->status()==SectionProposalStatus::noSolution?"REJECTED. ":"Preview only. ")
        +juce::String("Proposed global ")+juce::String(proposed->staticGainDb,9)+" dB. Click/arrow through changes for local + global gain and reasons. "
        +(automaticStale?"Apply unavailable. Generate a fresh proposal.":automatic->status()==SectionProposalStatus::noSolution?"Apply unavailable. Review the rejection.":"Apply proposal explicitly in the main controls.")
        +(automatic->previewHasFinalGlobal()?juce::String():" Global shown is baseline only; no final combined solution.");
    else if(session.sectionProposal)message="This proposal has no gain preview. Review the outcome, or Discard to edit manually.";
    else if(feedback.isNotEmpty())message=feedback;
    else if(dirty)message="Unstaged field. Return or Set draft stages it; Apply is separate.";
    else if(draftError.isNotEmpty())message=draftError;
    else if(hasDraft)message=session.sectionTrimDraft->processingRegions?session.sectionTrimNotice:"Draft only. Apply gains changes Current; Discard keeps applied audio.";
    else if(selectionKind.getSelectedId()!=1)message="Click an applied region, then Backspace to draft removal of its local gain. Apply/Discard remain explicit. Proposed regions are read-only.";
    else if(!geometryMatches)message="Annotations changed. Applied ramps stay fixed; a new draft starts at 0 dB.";
    else if(session.currentStateNeedsRemeasurement())message="Current unverified. Re-measure Current State to verify the applied gain.";
    else message="Click a section to inspect. Backspace removes its annotation into a neighbour; applied audio stays unchanged. Enter a gain and choose Set draft.";
    if(!geometryMatches&&fresh&&!session.sectionTrimDraft->processingRegions)message+=" Other sections start at 0 dB on changed annotations.";
    legend.setText(automatic?"CURRENT local gain (dB)   |   Green solid = applied   /   Cyan dotted = proposal":
        "CURRENT local gain (dB)   |   Green solid = applied   /   Amber dashed = manual draft",juce::dontSendNotification);
    notice.setText(message,juce::dontSendNotification);notice.setTooltip(message);
}
void TimelineGainLane::stageValue()
{
    if(!bindingMatches()||!stage.isEnabled()||session.audioEditsBlocked()||session.job.isActive()||session.sectionProposal)return;
    std::istringstream in(gain.getText().trim().toStdString());in.imbue(std::locale::classic());double value=0;
    if(!(in>>value)||!(in>>std::ws).eof()||!std::isfinite(value)||value < -12||value > 12)
    {feedback="Enter one finite value from -12 to +12 dB.";updateControls();return;}
    if(!hasDraft&&!session.beginSectionTrimDraft()){feedback=session.sectionTrimNotice;updateControls();return;}
    session.setSectionDraftGain(row,value);afterAction();
}
void TimelineGainLane::afterAction()
{
    dirty=false;synchronize(programme);showSelection();feedback=session.sectionTrimNotice;updateControls();
    if(onChanged)onChanged();
}
void TimelineGainLane::selectFrame(int64_t frame,bool notify)
{
    const auto& map=session.sectionDocument.map();if(!usableMap||map.sections.empty()||!bindingMatches())return;
    cursor=std::clamp<int64_t>(frame,0,map.frames-1);row=intervalAt(rows(),cursor);
    dirty=false;feedback.clear();showSelection();updateControls();repaint();
    if(notify&&onSelection){const auto& s=rows()[row];onSelection(s.begin,s.end,cursor);}
}
void TimelineGainLane::selectAt(int64_t frame,int y)
{
    const auto base=plotBounds().getBottom();
    if(y>=base&&y<base+20&&applied)selectionKind.setSelectedId(2,juce::dontSendNotification);
    else if(y>=base+20&&y<base+42)selectionKind.setSelectedId(proposed&&proposed->sectionGains?3:hasDraft&&fresh&&session.sectionTrimDraft->processingRegions?2:1,juce::dontSendNotification);
    rebuildProposalRows();selectFrame(frame,true);deleteArmed=selected()!=nullptr;
}
bool TimelineGainLane::deleteSelected()
{
    if(!deleteArmed||!selected())return false;
    deleteArmed=false; // One explicit selection per deletion; held Backspace cannot cascade.
    if(!bindingMatches()||session.audioEditsBlocked()||session.job.isActive()||session.persistenceLoading)
    {feedback="Stop playback and finish the job before deleting.";updateControls();return true;}
    if(dirty){feedback="Stage or discard the gain field edit before deleting.";updateControls();return true;}
    if(selectionKind.getSelectedId()==3)
    {feedback="Proposed regions belong to the complete measured proposal. Apply or Discard the proposal first.";notice.setText(feedback,juce::dontSendNotification);notice.setTooltip(feedback);return true;}
    if(selectionKind.getSelectedId()==2)
    {
        const auto region=*selected();
        if(session.draftRemoveAppliedRegion(region.begin,region.end)){afterAction();return true;}
        feedback="Cannot remove this region: it has no local gain, or an unapplied proposal/stale or section draft needs Discard first.";
        notice.setText(feedback,juce::dontSendNotification);notice.setTooltip(feedback);return true;
    }
    if(hasDraft){feedback="Apply or Discard the gain draft before deleting annotations.";updateControls();return true;}
    auto map=session.sectionDocument.map();if(map.sections.size()<2)
    {feedback="The final whole-programme annotation must remain; there are no section boundaries left to delete.";updateControls();return true;}
    deletedAnnotation=map;
    if(row>0)map.sections[row-1].end=map.sections[row].end;
    else map.sections[1].begin=0;
    map.sections.erase(map.sections.begin()+ptrdiff_t(row));++map.revision;deletedRevision=map.revision;
    if(!session.sectionDocument.restore(map,input)){deletedAnnotation.reset();return true;}
    session.documentChanged();afterAction();
    feedback="Section annotation deleted into its neighbour. Audio/gains/measurements unchanged. Undo delete restores the annotation.";
    notice.setText(feedback,juce::dontSendNotification);notice.setTooltip(feedback);resized();return true;
}
void TimelineGainLane::setViewport(double begin,double end,int left,int width)
{
    if(begin==viewBegin&&end==viewEnd&&left==plotLeft&&width==plotWidth)return;
    viewBegin=begin;viewEnd=end;plotLeft=left;plotWidth=width;repaint();
}
void TimelineGainLane::resized()
{
    auto b=getLocalBounds();auto top=b.removeFromTop(22);selectionKind.setBounds(top.removeFromRight(195));legend.setBounds(top);notice.setBounds(b.removeFromTop(24));
    auto controls=b.removeFromBottom(28).reduced(2,1);
    previous.setBounds(controls.removeFromLeft(27));next.setBounds(controls.removeFromLeft(27));
    if(deletedAnnotation)undoDelete.setBounds(controls.removeFromLeft(100));
    if(!automatic)
    {
        discard.setBounds(controls.removeFromRight(105));apply.setBounds(controls.removeFromRight(100));stage.setBounds(controls.removeFromRight(86));
        gain.setBounds(controls.removeFromRight(156));dbLabel.setBounds(controls.removeFromRight(60));
    }
    detail.setBounds(controls);
}
void TimelineGainLane::paint(juce::Graphics& g)
{
    g.fillAll(palette::panel());
    if(hasRegionLayers())paintRegionLayers(g);
    paintReadouts(g);
    const auto box=plotBounds();if(box.getWidth()<2||viewEnd<=viewBegin)return;
    const auto x=[&](double f){return float(box.getX()+(f-viewBegin)*(box.getWidth()-1)/(viewEnd-viewBegin));};
    const auto y=[&](double db){return float(box.getCentreY()-db*(box.getHeight()-6)/24.);};
    g.setFont(11.f);g.setColour(juce::Colours::lightgrey);
    for(int db:{12,0,-12})g.drawText(juce::String(db)+" dB",0,int(y(db))-7,box.getX()-5,14,juce::Justification::centredRight);
    g.reduceClipRegion(box);g.setColour(palette::card());g.fillRect(box);
    const auto* s=selected();if(!s)return;
    g.setColour(juce::Colour(0xff263440));
    if(s->end>viewBegin&&s->begin<viewEnd)g.fillRect(x(std::max<double>(s->begin,viewBegin)),float(box.getY()),x(std::min<double>(s->end,viewEnd))-x(std::max<double>(s->begin,viewBegin)),float(box.getHeight()));
    const auto& rows=this->rows();
    int previousX=-1;
    for(size_t i=intervalAt(rows,int64_t(viewBegin));i<rows.size()&&rows[i].begin<viewEnd;++i)
    {
        const auto& r=rows[i];if(r.begin<viewBegin)continue;const int px=int(x(r.begin));
        if(px!=previousX){g.setColour(juce::Colour(0xff566570));g.drawVerticalLine(px,float(box.getY()),float(box.getBottom()));previousX=px;}
    }
    if(!initialized)return;
    // Clip in frame space first, including ramps whose centres are outside the view.
    // Canonical transition endpoints match SectionGainPlan::gainAt, including shortened ramps.
    auto draw=[&](const processing::SectionGainPlan* p,int style)
    {
        const auto knots=gainOutline(p,rows.back().end,viewBegin,viewEnd,box.getWidth());
        if(knots.empty())return;
        g.setColour(style==2?juce::Colour(0xff83d9ff):style==1?juce::Colour(0xffffc171):juce::Colour(0xff81d5b4));
        for(size_t i=1;i<knots.size();++i)
        {
            const float x0=x(knots[i-1].frame),x1=x(knots[i].frame),y0=y(knots[i-1].db),y1=y(knots[i].db);
            // Solid stroke for applied values. Dash the draft in screen space:
            // tessellating a long dashed JUCE Path through thousands of steep
            // sub-pixel transitions can take far longer than painting the plot.
            if(style==0||x1-x0<1.f){g.drawLine(x0,y0,x1,y1,1.6f);continue;}
            const auto at=[&](float px){return y0+(px-x0)/(x1-x0)*(y1-y0);};
            for(int phase=int(std::floor((x0-box.getX())/11))*11+box.getX();phase<x1;phase+=11)
            {
                const auto left=std::max(x0,float(phase)),right=std::min(x1,float(phase+(style==2?2:6)));
                if(right>left)g.drawLine(left,at(left),right,at(right),1.8f);
            }
        }
    };
    draw(applied?&*applied:nullptr,0);if(draft)draw(&*draft,1);if(proposed)draw(proposed->sectionGains?&*proposed->sectionGains:nullptr,2);
}
void TimelineGainLane::paintReadouts(juce::Graphics& g)
{
    if(!usableMap||plotWidth<2||viewEnd<=viewBegin)return;
    const auto x=[&](int64_t frame){return plotLeft+int((std::clamp(double(frame),viewBegin,viewEnd)-viewBegin)*(plotWidth-1)/(viewEnd-viewBegin));};
    const auto draw=[&](const processing::SectionGainPlan* plan,int line,const juce::String& title,juce::Colour colour)
    {
        const int y=plotBounds().getBottom()+line*20;g.setColour(colour);g.setFont(11.f);
        g.drawText(title,0,y,plotLeft-4,20,juce::Justification::centredRight);
        const auto label=[&](int64_t begin,int64_t end,double db)
        {
            if(end<=viewBegin||begin>=viewEnd)return;
            const auto left=x(begin),right=x(end);g.drawVerticalLine(left,float(y),float(y+17));
            if(right-left>=24)g.drawText(juce::String(db,2)+" dB",left+2,y,right-left-4,18,juce::Justification::centred);
        };
        if(plan)for(size_t i=intervalAt(plan->intervals,int64_t(viewBegin));i<plan->intervals.size()&&plan->intervals[i].begin<viewEnd;++i)
        {const auto& r=plan->intervals[i];label(r.begin,r.end,r.db);}
        else for(const auto& r:session.sectionDocument.map().sections)label(r.begin,r.end,0);
    };
    draw(applied?&*applied:nullptr,0,"Applied",juce::Colour(0xff81d5b4));
    if(proposed)draw(proposed->sectionGains?&*proposed->sectionGains:nullptr,1,automaticStale?"Stale prop.":diagnosticRegionInspection()||automatic->status()==SectionProposalStatus::noSolution?"Rejected prop.":"Proposed",juce::Colour(0xff83d9ff));
    else if(draft)draw(&*draft,1,"Draft",juce::Colour(0xffffc171));
}
void TimelineGainLane::paintRegionLayers(juce::Graphics& g)
{
    if(plotWidth<2||viewEnd<=viewBegin)return;
    const auto x=[&](int64_t f){return float(plotLeft+(std::clamp(double(f),viewBegin,viewEnd)-viewBegin)*(plotWidth-1)/(viewEnd-viewBegin));};
    const auto band=[&](int row,const char* label,juce::Colour colour,const std::vector<proposal::RegionSpan>& spans)
    {
        const float y=float(47+row*15);g.setFont(10.f);g.setColour(colour);g.drawText(label,0,int(y),plotLeft-4,14,juce::Justification::centredRight);
        juce::Graphics::ScopedSaveState save(g);g.reduceClipRegion({plotLeft,int(y),plotWidth,14});g.setColour(juce::Colour(0xff111a23));g.fillRect(plotLeft,int(y),plotWidth,14);
        for(auto s:spans)if(s.end>viewBegin&&s.begin<viewEnd){g.setColour(colour.withAlpha(.3f));g.fillRect(x(s.begin),y+1,std::max(1.f,x(s.end)-x(s.begin)),12.f);g.setColour(colour);g.drawVerticalLine(int(x(s.begin)),y+1,y+13);g.drawVerticalLine(int(x(s.end)),y+1,y+13);}
    };
    std::vector<proposal::RegionSpan> annotation,roots,allowed,protectedAreas,support;
    for(const auto& row:session.sectionDocument.map().sections)annotation.push_back({row.begin,row.end});
    const auto review=automatic?automatic->regionReview():nullptr;
    const auto settings=review?std::optional<proposal::RegionPlanningSettings>(review->constraints.settings):session.regionPlanning;
    if(review){for(const auto& r:review->constraints.roots){roots.push_back(r.hint);allowed.push_back(r.allowed);}if(const auto* a=inspectedAttempt())support=a->plan.activity.support;}
    else if(settings)
    {
        const auto& map=session.sectionDocument.map();
        if(settings->mode==proposal::RegionMode::automatic)allowed.push_back({0,map.frames});
        else
        {
            auto hint=settings->selection;
            if(settings->mode==proposal::RegionMode::selectedSection)for(const auto& s:map.sections)if(s.id==settings->sectionId){hint={s.begin,s.end};break;}
            if(hint.end>hint.begin){roots.push_back(hint);allowed.push_back({std::max<int64_t>(0,hint.begin-int64_t(settings->leftMarginMs)*map.sampleRate/1000),std::min(map.frames,hint.end+int64_t(settings->rightMarginMs)*map.sampleRate/1000)});}
        }
    }
    if(settings)protectedAreas=settings->protectedLocal;
    band(0,"Annotations",juce::Colour(0xff9caabd),annotation);
    band(1,settings&&settings->mode==proposal::RegionMode::automatic?"Measured roots":"Root hint",juce::Colour(0xffd7b5ff),roots);
    band(2,"Allowed area",juce::Colour(0xff87c9b0),allowed);
    band(3,"Protected",juce::Colour(0xffff9f85),protectedAreas);
    band(4,"Changed support",juce::Colour(0xff83d9ff),support);
}
}

#include "TimelineGainLane.h"
#include "TimelineGainOutline.h"
#include "TimelineAppearance.h"
#include "RegionReviewFormatting.h"
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <unordered_set>

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
    for(auto* c:std::initializer_list<juce::Component*>{&legend,&notice,&detail,&previous,&next,&gain,&dbLabel})addAndMakeVisible(c);
    for(auto* label:{&legend,&notice,&detail,&dbLabel}){label->setFont(13.f);label->setInterceptsMouseClicks(false,false);}
    legend.setText("CURRENT section gain (dB)   |   Solid: applied   /   Dashed: draft",juce::dontSendNotification);
    legend.setTooltip("All programme channels, before global trim. This lane always describes Current processing.");
    dbLabel.setText("Local dB",juce::dontSendNotification);gain.setName("Timeline section gain dB");
    gain.setTooltip("Absolute local gain before global trim: -12 to +12 dB. This replaces the section setting; it is not added to an existing region gain. Return or leaving the field applies it. Command-Z undoes the edit.");
    gain.setInputRestrictions(36);gain.setSelectAllWhenFocused(true);gain.setJustification(juce::Justification::centred);
    gain.onTextChange=[this]{dirty=true;feedback.clear();updateControls();};
    gain.onReturnKey=[this]{if(stageValue())if(auto* parent=getParentComponent();parent&&parent->isShowing()&&parent->getPeer()){parent->setWantsKeyboardFocus(true);parent->grabKeyboardFocus();}};
    gain.onEscapeKey=[this]{dirty=false;feedback.clear();showSelection();updateControls();};
    gain.onFocusLost=[this]{if(dirty)stageValue();};
    undoDelete.onClick=[this]
    {
        if(!deletedAnnotation||!bindingMatches()||session.sectionDocument.map().revision!=deletedRevision
            ||session.audioEditsBlocked()||session.job.isActive()||session.sectionTrimDraft)return;
        auto restored=*deletedAnnotation;restored.revision=deletedRevision+1;
        if(session.sectionDocument.restore(restored,input)){session.documentChanged();deletedAnnotation.reset();afterAction();}
    };
    apply.onClick=[this]
    {
        if(dirty&&!stageValue())return;
        if(!bindingMatches()||!apply.isEnabled())return;
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
    addAndMakeVisible(createSections);createSections.setName("Create sections from applied regions");
    createSections.setMouseClickGrabsKeyboardFocus(false);
    createSections.onClick=[this]
    {
        const auto blocked=session.appliedRegionSectionsBlocked();
        if(dirty||!bindingMatches()||blocked.isNotEmpty())
        {feedback=dirty?"Finish or cancel the local gain entry first.":blocked.isNotEmpty()?blocked:"Timeline changed; review the current sections and retry.";updateControls();return;}
        const auto previousRevision=session.sectionDocument.map().revision;
        std::unordered_set<uint64_t> previousIds;for(const auto& s:session.sectionDocument.map().sections)previousIds.insert(s.id);
        juce::String message;
        if(!session.createSectionsFromAppliedRegions(message)){feedback=message;updateControls();return;}
        if(session.sectionDocument.map().revision==previousRevision)
        {feedback=message;updateControls();return;}
        int64_t selectedFrame=cursor;
        for(const auto& s:session.sectionDocument.map().sections)if(!previousIds.count(s.id)){selectedFrame=s.begin;break;}
        synchronize(programme);selectionKind.setSelectedId(1,juce::dontSendNotification);rebuildProposalRows();selectFrame(selectedFrame,true);
        feedback=message;updateControls();repaint();
        if(onChanged)onChanged();
    };
    addAndMakeVisible(selectionKind);selectionKind.setName("Gain lane selection");
    selectionKind.addItem("Select sections",1);selectionKind.addItem("Select applied regions",2);selectionKind.addItem("Select proposed regions",3);
    selectionKind.setSelectedId(1,juce::dontSendNotification);
    selectionKind.setTooltip("Sections are editable annotations. Applied regions have draggable shared-edge handles in the Current graph and at the top of the gain plot, including after reopening. Choose the item type here; clicks keep this selection mode.");
    selectionKind.onChange=[this]{if(dirty&&!stageValue()){selectionKind.setSelectedId(1,juce::dontSendNotification);return;}deleteArmed=false;rebuildProposalRows();selectFrame(cursor,true);};
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
    if(!review||regionAttempt>=review->attempts.size())return false;
    if(review->eventMode)return !review->selectedAttempt||*review->selectedAttempt!=regionAttempt||automatic->status()!=SectionProposalStatus::ready||review->attempts[regionAttempt].evaluation.verdict!=proposal::Verdict::pass;
    return regionAttempt!=review->retainedAttempt||review->attempts[regionAttempt].evaluation.verdict!=proposal::Verdict::pass;
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
        deleteArmed=false;deletedAnnotation.reset();selectedIds.clear();rangeSelection=false;
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
    if(replaced||boundarySerial!=session.boundaryEditSerial||appliedBoundarySerial!=session.appliedBoundarySerial||automatic!=session.sectionProposal||automaticStale!=stale||audioRevision!=session.revision||mapRevision!=map.revision||initialized!=session.hasCurrentStateInitialized()
        ||hasDraft!=bool(d)||(d&&editSerial!=d->editSerial))
    {
        const bool newProposal=automatic!=session.sectionProposal;
        const auto oldReview=automatic?automatic->regionReview():nullptr,newReview=session.sectionProposal?session.sectionProposal->regionReview():nullptr;
        const bool sameInvocation=oldReview&&newReview&&oldReview->eventMode&&newReview->eventMode&&oldReview->resultIdentity==newReview->resultIdentity;
        const auto inspectedId=sameInvocation&&regionAttempt<oldReview->attempts.size()?oldReview->attempts[regionAttempt].attemptIdentity:juce::String();
        deleteArmed=false;
        automatic=session.sectionProposal;automaticStale=stale;proposed=automatic?automatic->previewPlan():std::nullopt;
        if(newProposal&&!sameInvocation)regionAttempt=automatic&&automatic->regionReview()?automatic->regionReview()->retainedAttempt:std::numeric_limits<size_t>::max();
        if(sameInvocation&&inspectedId.isNotEmpty())
        {const auto it=std::find_if(newReview->attempts.begin(),newReview->attempts.end(),[&](const auto& a){return a.attemptIdentity==inspectedId;});regionAttempt=size_t(it-newReview->attempts.begin());}
        if(const auto* attempt=inspectedAttempt()){proposed=ProcessingPlan{};proposed->staticGainDb=attempt->plan.globalDb;proposed->sectionGains=attempt->plan.envelope;}
        audioRevision=session.revision;mapRevision=map.revision;initialized=session.hasCurrentStateInitialized();
        boundarySerial=session.boundaryEditSerial;appliedBoundarySerial=session.appliedBoundarySerial;
        usableMap=input.isNotEmpty()&&map.inputIdentity==input&&map.valid();
        hasDraft=bool(d);editSerial=d?d->editSerial:0;fresh=usableMap&&session.sectionDraftFresh();
        applied=usableMap&&initialized?session.currentStateDocument.sectionGains:std::nullopt;
        // Keep the user-selected inspection mode when proposals arrive.
        if(selectionKind.getSelectedId()==3&&(!proposed||!proposed->sectionGains))selectionKind.setSelectedId(applied?2:1,juce::dontSendNotification);
        if(selectionKind.getSelectedId()==2&&!applied)selectionKind.setSelectedId(1,juce::dontSendNotification);
        selectionKind.setItemEnabled(2,bool(applied));selectionKind.setItemEnabled(3,proposed&&proposed->sectionGains.has_value());
        rebuildProposalRows();
        geometryMatches=!applied||session.sectionGeometry(map)==session.sectionGeometry(session.currentStateDocument.sectionGainAnnotations);
        supportError=usableMap&&initialized?session.sectionDraftSupportError():juce::String{};
        draftError=fresh?d->plan.validationError():juce::String{};
        draft=fresh&&draftError.isEmpty()?std::optional<processing::SectionGainPlan>(d->plan.canonical()):std::nullopt;
        if(session.pendingBoundaryEdits())draft=session.displayedSectionGains();
        cursor=std::clamp<int64_t>(cursor,0,std::max<int64_t>(0,map.frames-1));row=intervalAt(rows(),cursor);
        selectedIds.erase(std::remove_if(selectedIds.begin(),selectedIds.end(),[&](auto id){return std::none_of(map.sections.begin(),map.sections.end(),[&](const auto& r){return r.id==id;});}),selectedIds.end());
        dirty=false;feedback.clear();showSelection();resized();repaint();
    }
    updateControls();
}
void TimelineGainLane::showSelection()
{
    const auto* section=selected();
    if(!section){detail.setText("No confirmed sections",juce::dontSendNotification);detail.setTooltip({});gain.setText({},false);return;}
    const auto* displayed=fresh?&session.sectionTrimDraft->plan:applied?&*applied:nullptr;
    auto stored=section_edit::uniform(displayed,section->begin,section->end);
    if(selectedIds.size()>1)for(const auto& r:session.sectionDocument.map().sections)if(std::find(selectedIds.begin(),selectedIds.end(),r.id)!=selectedIds.end())
        if(section_edit::uniform(displayed,r.begin,r.end)!=stored)stored.reset();
    gain.setText(stored?exact(*stored):juce::String(),false);
    gain.setTextToShowWhenEmpty(stored?juce::String():"Mixed",juce::Colours::lightgrey);
    int64_t selectedEnd=section->end;
    if(selectedIds.size()>1)for(const auto& r:session.sectionDocument.map().sections)if(std::find(selectedIds.begin(),selectedIds.end(),r.id)!=selectedIds.end())selectedEnd=std::max(selectedEnd,r.end);
    auto description=selectedIds.size()>1?juce::String(selectedIds.size())+" selected sections in "+juce::String(selectedSpans().size())+" passage(s)":section->name+"   ["+juce::String(section->begin)+", "+juce::String(selectedEnd)+")";
    if(initialized)description+="   Applied: "+exact(applied?applied->gainAt(cursor):0.)+" dB";
    if(selectionKind.getSelectedId()!=1)description+=" (at cursor; nominal interval shown, transitions may cross its edges).";
    if(fresh&&session.sectionTrimDraft->processingRegions)description+=" | draft local "+exact(draft?draft->gainAt(cursor):0.)+" dB";
    if(initialized)description+=" | applied global "+exact(session.currentStateDocument.staticGainDb)+" dB";
    if(proposed)
    {
        const double local=proposed->sectionGains?proposed->sectionGains->gainAt(cursor):0.;
        description+="   Proposed local "+exact(local)+" dB / global "+exact(proposed->staticGainDb)+" dB / combined "+exact(local+proposed->staticGainDb)+" dB";
        if(!automatic->previewHasFinalGlobal())description+=" (rejected local selection with baseline global gain; combined candidate not measured)";
        if(const auto* attempt=inspectedAttempt())
        {
            const auto& review=*automatic->regionReview();
            description="Measured attempt "+juce::String(regionAttempt+1)+" — "+regionAttemptAuthority(*automatic,regionAttempt)+". "+description;
            if(review.eventMode)description+="\n"+eventCoverageSummary(attempt->coverage);
            for(const auto& a:attempt->plan.assignments)if(a.geometry.support.begin<=cursor&&cursor<a.geometry.support.end)
            {
                description+="\nRoot "+(a.geometry.root<review.constraints.roots.size()?review.constraints.roots[a.geometry.root].id:"unavailable")+". "+a.geometry.reason;
                if(review.constraints.settings.schema==3)description+="\n"+regionBoundaryFamily(a.geometry.family)+". Review shared boundary evidence for competing beginnings/returns and unresolved coverage.";
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
    const auto creationBlock=session.appliedRegionSectionsBlocked();
    createSections.setEnabled(!dirty&&bindingMatches()&&creationBlock.isEmpty());
    createSections.setTooltip(dirty?"Finish or cancel the local gain entry first.":creationBlock.isNotEmpty()?creationBlock:
        "Copy actual applied nominal cuts into sections once. Existing names and cuts, gains, sound and measured Current stay unchanged. Sections are not linked to later region moves. Undo: Command-Z.");
    const bool editable=(!rangeSelection||!selectedIds.empty())&&selectionKind.getSelectedId()==1&&selected()&&initialized&&supportError.isEmpty()&&idle&&session.manualEditAllowed();
    previous.setEnabled(selected()&&row>0);next.setEnabled(selected()&&row+1<rows().size());
    if(selectionKind.getSelectedId()==1&&automatic&&!automatic->selection().corrections.empty())
    {
        const auto& changes=automatic->selection().corrections;
        previous.setEnabled(selected()&&cursor>changes.front().begin);next.setEnabled(selected()&&cursor<changes.back().begin);
    }
    previous.setTooltip("Inspect previous section or processing region");next.setTooltip("Inspect next section or processing region");
    for(auto* c:std::initializer_list<juce::Component*>{&gain,&dbLabel})c->setVisible(!automatic||session.lastAppliedProposalBinding==automatic->identity());
    apply.setVisible(false);discard.setVisible(false);undoDelete.setVisible(false);
    gain.setEnabled(editable);
    std::istringstream field(gain.getText().trim().toStdString());field.imbue(std::locale::classic());double entered=0;
    const bool validField=(field>>entered)&&(field>>std::ws).eof()&&std::isfinite(entered)&&entered>=-12&&entered<=12;
    apply.setEnabled((!dirty||validField)&&initialized&&idle&&!session.sectionProposal&&supportError.isEmpty()&&(fresh||dirty)&&draftError.isEmpty());discard.setEnabled((hasDraft||dirty)&&idle);
    undoDelete.setVisible(false);undoDelete.setEnabled(deletedAnnotation&&idle&&!hasDraft&&bindingMatches());
    juce::String message;
    if(!usableMap)message="Confirm programme audio and channel roles to inspect sections.";
    else if(!initialized)message="Current is not initialised.";
    else if(hasDraft&&!fresh)message="Stale draft. Discard it and draft on the current sections.";
    else if(supportError.isNotEmpty())message=supportError;
    else if(!idle)message="Inspection available. Stop playback or finish the job to edit.";
    else if(proposed&&automatic->regionReview()&&automatic->regionReview()->eventMode)message=automaticStale?"STALE event proposal; Apply unavailable.":regionAttemptAuthority(*automatic,regionAttempt)+(diagnosticRegionInspection()?". Apply unavailable for this inspection; use Review to select an owned measured pass.":". Apply explicitly while stopped, then re-measure Current.");
    else if(proposed)message=(automaticStale?"STALE. ":diagnosticRegionInspection()?"DIAGNOSTIC ATTEMPT — Apply unavailable. ":automatic->status()==SectionProposalStatus::noSolution?"REJECTED. ":"Preview only. ")
        +juce::String("Proposed global ")+juce::String(proposed->staticGainDb,9)+" dB. Click/arrow through changes for local + global gain and reasons. "
        +(automaticStale?"Apply unavailable. Generate a fresh proposal.":automatic->status()==SectionProposalStatus::noSolution?"Apply unavailable. Review the rejection.":"Apply proposal explicitly in the main controls.")
        +(automatic->previewHasFinalGlobal()?juce::String():" Global shown is baseline only; no final combined solution.");
    else if(session.sectionProposal)message="This proposal has no gain preview. Review the outcome, or Discard to edit manually.";
    else if(feedback.isNotEmpty())message=feedback;
    else if(dirty)message="Editing local gain. Return or leaving the field applies it. Command-Z undoes it.";
    else if(draftError.isNotEmpty())message=draftError;
    else if(hasDraft)message=session.sectionTrimDraft->processingRegions?session.sectionTrimNotice:"Draft only. Apply gains changes Current; Discard keeps applied audio.";
    else if(selectionKind.getSelectedId()!=1)message=selectionKind.getSelectedId()==2?appliedBoundaryHelp():"Proposed regions are read-only. Choose Select applied regions to edit Current.";
    else if(!geometryMatches)message="Sections and gain regions differ. New edits retain existing gains outside the selected passage.";
    else if(session.currentStateNeedsRemeasurement())message="Current needs re-measurement. Continue editing, or re-measure when ready.";
    else message="Click a gain readout to edit. Drag this lane to select whole sections. Backspace resets/merges; Command-Z undoes.";

    legend.setText(automatic?"CURRENT local gain (dB)   |   Green solid = applied   /   Cyan dotted = proposal":
        "CURRENT local gain (dB)   |   Shaded: transition   /   Command-Z: Undo",juce::dontSendNotification);
    if(session.pendingBoundaryEdits()){message=session.pendingBoundaryMessage();legend.setText("CURRENT local gain | Solid: applied / Amber dashed: pending boundaries",juce::dontSendNotification);}
    notice.setText(message,juce::dontSendNotification);notice.setTooltip(message);
}
bool TimelineGainLane::stageValue()
{
    if(!dirty)return true;
    if(!bindingMatches()||!gain.isEnabled()||session.audioEditsBlocked()||session.job.isActive()||!session.manualEditAllowed())return false;
    std::istringstream in(gain.getText().trim().toStdString());in.imbue(std::locale::classic());double value=0;
    if(!(in>>value)||!(in>>std::ws).eof()||!std::isfinite(value)||value < -12||value > 12)
    {feedback="Enter one finite value from -12 to +12 dB. Nothing applied.";updateControls();return false;}
    auto ids=selectedIds;if(ids.empty()&&selected())ids.push_back(selected()->id);
    if(!session.setManualSectionGains(ids,value)){feedback=session.sectionTrimNotice;updateControls();return false;}
    afterAction();return true;
}
void TimelineGainLane::afterAction()
{
    dirty=false;synchronize(programme);showSelection();feedback=session.sectionTrimNotice;updateControls();
    if(selectionKind.getSelectedId()==1)deleteArmed=!selectedIds.empty();
    if(onChanged)onChanged();
}
void TimelineGainLane::selectFrame(int64_t frame,bool notify)
{
    if(dirty&&!stageValue())return;
    const auto& map=session.sectionDocument.map();if(!usableMap||map.sections.empty()||!bindingMatches())return;
    cursor=std::clamp<int64_t>(frame,0,map.frames-1);row=intervalAt(rows(),cursor);
    rangeSelection=false;selectedIds.clear();if(selectionKind.getSelectedId()==1&&selected())selectedIds.push_back(selected()->id);
    dirty=false;feedback.clear();showSelection();updateControls();resized();repaint();
    if(notify&&onSelection){const auto& s=rows()[row];onSelection(s.begin,s.end,cursor);}
}
void TimelineGainLane::selectAt(int64_t frame,int /*y*/)
{
    rebuildProposalRows();selectFrame(frame,true);deleteArmed=selected()!=nullptr;
}
void TimelineGainLane::selectRange(int64_t begin,int64_t end)
{
    if(dirty&&!stageValue())return;
    if(!bindingMatches()||selectionKind.getSelectedId()!=1)return;
    if(end<begin)std::swap(begin,end);rangeSelection=true;selectedIds.clear();
    for(const auto& r:rows())if(r.begin>=begin&&r.end<=end)selectedIds.push_back(r.id);
    deleteArmed=!selectedIds.empty();
    if(deleteArmed)
    {
        const auto first=std::find_if(rows().begin(),rows().end(),[&](const auto& r){return r.id==selectedIds.front();});
        const auto last=std::find_if(rows().begin(),rows().end(),[&](const auto& r){return r.id==selectedIds.back();});
        row=size_t(first-rows().begin());cursor=first->begin;
        if(onSelection)onSelection(first->begin,last->end,cursor);
    }
    showSelection();updateControls();gain.setEnabled(deleteArmed&&session.manualEditAllowed());resized();repaint();
    if(!deleteArmed)notice.setText("No fully enclosed sections selected.",juce::dontSendNotification);
}
std::vector<std::pair<int64_t,int64_t>> TimelineGainLane::selectedSpans() const
{
    std::vector<std::pair<int64_t,int64_t>> spans;if(selectionKind.getSelectedId()!=1)return spans;
    const std::unordered_set<uint64_t> ids(selectedIds.begin(),selectedIds.end());
    for(const auto& r:session.sectionDocument.map().sections)if(ids.count(r.id))
    {if(!spans.empty()&&spans.back().second==r.begin)spans.back().second=r.end;else spans.emplace_back(r.begin,r.end);}
    return spans;
}
void TimelineGainLane::toggleSection(int64_t frame)
{
    if(!commitPending()||!bindingMatches()||selectionKind.getSelectedId()!=1||session.boundaryBatch)return;
    const auto& map=session.sectionDocument.map();auto i=intervalAt(map.sections,std::clamp<int64_t>(frame,0,map.frames-1));if(i>=map.sections.size())return;
    const auto id=map.sections[i].id;auto found=std::find(selectedIds.begin(),selectedIds.end(),id);
    if(found==selectedIds.end())selectedIds.push_back(id);else selectedIds.erase(found);
    rangeSelection=true;deleteArmed=!selectedIds.empty();
    if(deleteArmed){row=0;while(row<map.sections.size()&&std::find(selectedIds.begin(),selectedIds.end(),map.sections[row].id)==selectedIds.end())++row;cursor=map.sections[row].begin;}
    dirty=false;feedback.clear();showSelection();updateControls();resized();repaint();
    const auto spans=selectedSpans();if(onSelection)onSelection(spans.empty()?0:spans.front().first,spans.empty()?0:spans.back().second,cursor);
    if(spans.empty())notice.setText("No sections selected. Command-click a section to add it.",juce::dontSendNotification);
}
void TimelineGainLane::editReadout(int64_t frame)
{
    const auto spans=selectedSpans();
    if(selectedIds.size()<2||std::none_of(spans.begin(),spans.end(),[&](const auto& s){return frame>=s.first&&frame<s.second;}))selectAt(frame,0);
    resized();if(gain.isEnabled()&&gain.isShowing()){gain.grabKeyboardFocus();gain.selectAll();}
}
bool TimelineGainLane::deleteSelected()
{
    if(!deleteArmed||!selected())return false;
    if(dirty&&!stageValue())return true;
    if(selectionKind.getSelectedId()==3){feedback="Apply or Discard the proposal before manual editing.";updateControls();return true;}
    if(!bindingMatches()||!session.manualEditAllowed())return true;
    const auto spans=selectedSpans();
    if(spans.size()>1)
    {
        deleteArmed=false;sectionEdit.requestResetSeparate(*this,spans,[safe=juce::Component::SafePointer<TimelineGainLane>(this)](bool,bool){if(safe){safe->selectedIds.clear();safe->afterAction();}});return true;
    }
    auto begin=selected()->begin,end=selected()->end;
    if(selectionKind.getSelectedId()==1&&!selectedIds.empty())for(const auto& r:rows())if(std::find(selectedIds.begin(),selectedIds.end(),r.id)!=selectedIds.end())
    {begin=std::min(begin,r.begin);end=std::max(end,r.end);}
    deleteArmed=false;
    sectionEdit.requestDelete(*this,begin,end,[safe=juce::Component::SafePointer<TimelineGainLane>(this)](bool,bool){if(safe){safe->selectedIds.clear();safe->afterAction();}});
    return true;
}
bool TimelineGainLane::removeSelectedBoundaries()
{
    if(selectionKind.getSelectedId()!=1||selectedIds.empty()||!bindingMatches()||!session.manualEditAllowed()||!stageValue())return false;
    if(selectedSpans().size()>1){feedback="Remove boundaries only needs one contiguous passage. Separate selections can reset their gains with Delete.";updateControls();return false;}
    auto next=session.manualState();auto& rows=next.map.sections;size_t first=rows.size(),last=0;
    for(size_t n=0;n<rows.size();++n)if(std::find(selectedIds.begin(),selectedIds.end(),rows[n].id)!=selectedIds.end()){first=std::min(first,n);last=n;}
    if(first==rows.size()||rows.size()==1)return false;
    if(first>0)--first;if(last+1<rows.size())++last;
    rows[first].end=rows[last].end;rows.erase(rows.begin()+ptrdiff_t(first+1),rows.begin()+ptrdiff_t(last+1));
    const auto accepted=session.commitManualEdit(std::move(next),"Remove section boundaries only");
    if(accepted){selectedIds.clear();afterAction();}return accepted;
}
void TimelineGainLane::setViewport(double begin,double end,int left,int width)
{
    if(begin==viewBegin&&end==viewEnd&&left==plotLeft&&width==plotWidth)return;
    viewBegin=begin;viewEnd=end;plotLeft=left;plotWidth=width;resized();repaint();
}
void TimelineGainLane::resized()
{
    auto b=getLocalBounds();auto top=b.removeFromTop(22);selectionKind.setBounds(top.removeFromRight(195));top.removeFromRight(6);createSections.setBounds(top.removeFromRight(264));top.removeFromRight(6);legend.setBounds(top);notice.setBounds(b.removeFromTop(24));
    auto controls=b.removeFromBottom(28).reduced(2,1);
    previous.setBounds(controls.removeFromLeft(27));next.setBounds(controls.removeFromLeft(27));
    if(!automatic||session.lastAppliedProposalBinding==automatic->identity())
    {
        gain.setBounds(controls.removeFromRight(110));dbLabel.setVisible(true);dbLabel.setBounds(controls.removeFromRight(60));
        if(const auto* s=selected();s&&selectedIds.size()<=1&&selectionKind.getSelectedId()==1&&viewEnd>viewBegin)
        {
            const auto left=plotLeft+int((std::max(viewBegin,double(s->begin))-viewBegin)*plotWidth/(viewEnd-viewBegin));
            const auto right=plotLeft+int((std::min(viewEnd,double(s->end))-viewBegin)*plotWidth/(viewEnd-viewBegin));
            if(right-left>=55){const auto width=std::min(140,right-left-4);gain.setBounds((left+right-width)/2,plotBounds().getBottom(),width,20);dbLabel.setVisible(false);}
        }
    }
    detail.setBounds(controls);
}
std::optional<size_t> TimelineGainLane::appliedBoundaryHandle(juce::Point<int> point) const
{
    const auto box=plotBounds();
    return point.y>=box.getY()&&point.y<=box.getY()+12?appliedBoundaryAt(float(point.x)):std::nullopt;
}
std::optional<size_t> TimelineGainLane::appliedBoundaryAt(float pointerX) const
{
    const auto box=plotBounds();const auto& plan=session.displayedSectionGains();
    if(selectionKind.getSelectedId()!=2||!plan||viewEnd<=viewBegin||box.getWidth()<2
        ||pointerX<box.getX()||pointerX>box.getRight())return {};
    std::optional<size_t> result;double distance=7;
    for(size_t n=1;n<plan->intervals.size();++n)
    {
        const auto frame=plan->intervals[n].begin;if(frame<viewBegin||frame>viewEnd)continue;
        const auto px=box.getX()+(frame-viewBegin)*(box.getWidth()-1)/(viewEnd-viewBegin);
        const auto d=std::abs(pointerX-px);if(d<distance){distance=d;result=n;}
    }
    return result;
}
juce::String TimelineGainLane::appliedBoundaryHelp() const
{
    if(session.pendingBoundaryEdits())return session.pendingBoundaryMessage();
    if(!session.manualEditAllowed())return "Applied-region handles unavailable: stop playback, finish the job and discard any pending proposal or gain draft first.";
    return "Edit Current: drag a boundary in the Current graph or a top-edge gain tab to move the shared cut and its ramp. Both neighbouring gains stay fixed. Release applies; Esc cancels; Command-Z undoes. Programme ends are fixed.";
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
    if((!rangeSelection||!selectedIds.empty())&&s->end>viewBegin&&s->begin<viewEnd)g.fillRect(x(std::max<double>(s->begin,viewBegin)),float(box.getY()),x(std::min<double>(s->end,viewEnd))-x(std::max<double>(s->begin,viewBegin)),float(box.getHeight()));
    const auto& rows=this->rows();
    int previousX=-1;
    for(size_t i=intervalAt(rows,int64_t(viewBegin));i<rows.size()&&rows[i].begin<viewEnd;++i)
    {
        const auto& r=rows[i];if(r.begin<viewBegin)continue;const int px=int(x(r.begin));
        if(px!=previousX){g.setColour(juce::Colour(0xff566570));g.drawVerticalLine(px,float(box.getY()),float(box.getBottom()));previousX=px;}
    }
    if(selectedIds.size()>1)for(const auto& r:session.sectionDocument.map().sections)
        if(std::find(selectedIds.begin(),selectedIds.end(),r.id)!=selectedIds.end()&&r.end>viewBegin&&r.begin<viewEnd)
        {g.setColour(juce::Colour(0xff81d5b4).withAlpha(.2f));g.fillRect(x(std::max<double>(r.begin,viewBegin)),float(box.getY()),x(std::min<double>(r.end,viewEnd))-x(std::max<double>(r.begin,viewBegin)),float(box.getHeight()));}
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
    const auto& visible=session.displayedSectionGains();
    if(visible)for(size_t n=1;n<visible->intervals.size();++n)
    {
        const auto radius=visible->radius(n),frame=visible->intervals[n].begin;
        const auto a=x(frame-radius),b=x(frame+radius);g.setColour(juce::Colour(0xff81d5b4).withAlpha(.12f));g.fillRect(a,float(box.getY()),std::max(1.f,b-a),float(box.getHeight()));
    }
    draw(applied?&*applied:nullptr,0);if(draft||session.pendingBoundaryEdits())draw(draft?&*draft:nullptr,1);if(proposed)draw(proposed->sectionGains?&*proposed->sectionGains:nullptr,2);
    if(selectionKind.getSelectedId()==2&&visible)
    {
        // Tabs sit on nominal cuts; the translucent band is the processor's
        // actual ramp support on either side. Programme endpoints stay fixed.
        g.setColour(session.manualEditAllowed()||session.appliedBoundaryDrag?juce::Colour(0xff81d5b4):juce::Colours::grey);
        for(size_t n=1;n<visible->intervals.size();++n)
        {
            const auto f=visible->intervals[n].begin;if(f<viewBegin||f>viewEnd)continue;
            const auto px=x(f);g.drawLine(px,float(box.getY()+8),px,float(box.getBottom()),1.f);
            g.fillRoundedRectangle(px-5,float(box.getY()),10.f,8.f,2.f);
        }
    }
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
    if(selectionKind.getSelectedId()==1)
    {
        const auto y=plotBounds().getBottom();g.setFont(11.f);
        for(const auto& r:session.sectionDocument.map().sections)if(r.end>viewBegin&&r.begin<viewEnd)
        {
            const auto l=x(r.begin),w=x(r.end)-l;const bool chosen=std::find(selectedIds.begin(),selectedIds.end(),r.id)!=selectedIds.end();
            if(chosen){g.setColour(juce::Colour(0xff81d5b4).withAlpha(.2f));g.fillRect(l,plotBounds().getY(),w,plotBounds().getHeight());}
            g.setColour(juce::Colour(0xff81d5b4));g.drawVerticalLine(l,float(y),float(y+18));
            const auto value=section_edit::uniform(applied?&*applied:nullptr,r.begin,r.end);
            if(w>=24)g.drawText(value?juce::String(*value,2)+" dB":"Mixed",l+2,y,w-4,18,juce::Justification::centred);
        }
    }
    else draw(applied?&*applied:nullptr,0,"Applied",juce::Colour(0xff81d5b4));
    if(proposed)draw(proposed->sectionGains?&*proposed->sectionGains:nullptr,1,automaticStale?"Stale prop.":automatic->regionReview()&&automatic->regionReview()->eventMode?(diagnosticRegionInspection()?"Review only":"Selected prop."):diagnosticRegionInspection()||automatic->status()==SectionProposalStatus::noSolution?"Rejected prop.":"Proposed",juce::Colour(0xff83d9ff));
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

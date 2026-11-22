#include "TimelineAppearance.h"
#include "ProposalWorkflowPanel.h"
#include "ProposalPresentation.h"
namespace levelscope::standalone
{
ProposalWorkflowPanel::ProposalWorkflowPanel(StandaloneSessionModel& s):session(s)
{
    for(auto* b:{&propose,&review,&apply,&discard,&history,&preferences,&regions,&recover,&boundaries})addAndMakeVisible(b);
    addAndMakeVisible(alternatives);alternatives.setName("Measured passing alternatives");
    alternatives.onChange=[this]{if(alternatives.isEnabled()&&onAlternative&&alternatives.getSelectedId()>0)onAlternative(size_t(alternatives.getSelectedId()-1));};
    recover.onClick=[this]{if(recover.isEnabled()&&onRecover)onRecover();};
    boundaries.onClick=[this]{if(boundaries.isEnabled()&&onBoundaries)onBoundaries();};
    addAndMakeVisible(method);method.setFont(13.f);method.setMinimumHorizontalScale(1.f);
    propose.setColour(juce::TextButton::buttonColourId,juce::Colour(0xff346967));
    headline.setFont(juce::Font(17.f,juce::Font::bold));headline.setColour(juce::Label::textColourId,juce::Colour(0xffb7e6db));
    availability.setFont(12.f);availability.setMinimumHorizontalScale(1.f);
    reason.setMultiLine(true);reason.setReadOnly(true);reason.setScrollbarsShown(true);reason.setFont(14.f);reason.setName("Proposal outcome and next action");
    reason.setColour(juce::TextEditor::backgroundColourId,palette::card());
    for(auto* c:std::initializer_list<juce::Component*>{&headline,&reason,&availability})addAndMakeVisible(c);
    propose.onClick=[this]{if(propose.isEnabled()&&onPropose)onPropose();};
    review.onClick=[this]{if(review.isEnabled()&&onReview)onReview();};
    history.onClick=[this]{if(history.isEnabled()&&onHistory)onHistory();};
    preferences.onClick=[this]{if(preferences.isEnabled()&&onPreferences)onPreferences();};
    regions.onClick=[this]{if(onRegions)onRegions();};
    apply.onClick=[this]{if(!apply.isEnabled())return;session.applySectionProposal();refresh();if(onChanged)onChanged();};
    discard.onClick=[this]{if(!discard.isEnabled())return;session.discardSectionProposal();refresh();if(onChanged)onChanged();};
    refresh();
}
void ProposalWorkflowPanel::refresh(float progress)
{
    const bool busy=session.audioEditsBlocked()||session.job.isActive()||session.persistenceLoading;
    const bool stale=session.sectionProposalStale();
    juce::String blocked;
    if(session.persistenceLoading)blocked="Finish opening the session.";
    else if(session.job.isActive())blocked="Finish or cancel the current job.";
    else if(session.pendingBoundaryEdits())blocked=session.pendingBoundaryMessage();
    else if(session.audioEditsBlocked())blocked="Stop playback to change processing.";
    else if(session.sourceContentInvalid||!session.sourceDocument.matchesFileStamps())blocked="Resolve the programme audio before proposing.";
    else if(!session.source.hasMeasurement)blocked="Measure Source first.";
    else if(!session.currentState.hasMeasurement||session.currentStateNeedsRemeasurement())blocked="Re-measure Current State first.";
    else if(session.sectionTrimDraft)blocked="Apply or discard the manual gain draft first.";
    if(operation)blocked=ProposalOperation::blocked(session);
    propose.setEnabled(blocked.isEmpty());review.setEnabled(bool(session.sectionProposal)||(operation&&bool(operation->result)));
    apply.setEnabled(session.canApplySectionProposal()&&!diagnosticInspection);discard.setEnabled(bool(session.sectionProposal)&&!busy);
    method.setText(session.regionPlanning?"Planning: "+proposalPlanningName(*session.regionPlanning)+" | "+juce::String(session.regionPlanning->mode==proposal::RegionMode::automatic?"automatic":session.regionPlanning->mode==proposal::RegionMode::selectedSection?"selected section":"temporary selection")
        +(session.regionPlanning->allowGlobal?" | whole-programme global trim permitted":" | global gain fixed"):"Planning: Standard proposals (default)",juce::dontSendNotification);
    method.setTooltip(method.getText()+". Settings remain private until confirmed."+(session.regionPlanning?"\n"+proposalDefinitionText(*session.regionPlanning):juce::String()));
    history.setEnabled(bool(session.lastAppliedExplanation));preferences.setEnabled(!busy);regions.setEnabled(!busy);
    juce::String title="Proposal workflow",detail=session.sectionProposalNotice;
    if(session.job.stage==JobStage::proposal&&session.job.isActive())
    {
        title=session.job.state==JobState::cancelling?"Cancelling proposal":"Proposing — Current remains unchanged";
        if(progress>=0)title+="  "+juce::String(juce::roundToInt(100*progress))+"%";
        detail=session.regionPlanning?
            "Experimental region search: one auxiliary energy pass, at most 64 proxy plans and three complete private measurements. Forecasts do not establish compliance. Global permission and complete old/new transition support are enforced; Current remains unchanged.":session.getSelectedTargetProfile().limits.hasMaxLoudnessRange?
            "Measured LRA passages guide bounded attenuation; at most three private assessments. Estimates are not measurements. Final global trim can be positive; Apply requires an exact complete candidate pass. Cancel keeps applied settings.":
            "Global gain is considered first. If it cannot meet IL and true peak together, one section envelope is assessed with up to three private passes. Cancel job keeps the applied settings.";
    }
    else if(stale)
    {title="Stale inputs — Apply unavailable";detail="Requirements, preferences, annotations, measured inputs or source files changed. Review shows the previous candidate. Choose Propose to replace it, or Discard proposal.";}
    else if(session.sectionProposalError!=ErrorCode::none)
    {
        title=session.sectionProposalError==ErrorCode::cancelled?"Proposal cancelled":
            session.sectionProposalError==ErrorCode::staleRevision||session.sectionProposalError==ErrorCode::sourceChanged?"Stale inputs — proposal discarded":"Proposal failed";
        detail=session.sectionProposalNotice+"\nApplied processing is preserved. Resolve the reported cause, then Propose again.";
    }
    else if(session.sectionProposal)
    {title=session.sectionProposal->outcomeTitle();detail=session.sectionProposal->outcomeReason();}
    else if(session.sectionProposalNotice.startsWith("APPLIED"))
    {title="Proposal applied — Current unverified";detail="Choose Re-measure Current State for independent verification. Last applied explanation retains the baseline, gains, reasons and available private candidate results.";}
    juce::String actions;
    if(!propose.isEnabled())actions="Propose unavailable. "+blocked;
    if(!apply.isEnabled())
    {
        if(actions.isNotEmpty())actions+="  ";
        actions+="Apply unavailable. ";
        actions+=session.pendingBoundaryEdits()?session.pendingBoundaryMessage():busy?"Stop playback or finish the job.":stale?"Generate a fresh proposal.":session.sectionProposal&&session.sectionProposal->status()==SectionProposalStatus::reviewRequired?"Review a passing attempt, then choose Use this measured alternative. Inspection alone does not select it.":diagnosticInspection?"A different or unselected attempt is displayed. Inspect the selected measured attempt in Review; selecting an alternative is a separate action.":session.sectionProposal?"This outcome has no applicable correction.":"No ready proposal.";
    }
    if(!review.isEnabled())actions+="  Review needs a completed proposal.";
    if(!discard.isEnabled())actions+=session.sectionProposal?"  Discard needs an idle transport/job.":"  No proposal to discard.";
    const bool normal=operation&&operation->phase!=ProposalOperation::Phase::idle;
    const bool legacy=bool(session.sectionProposal);
    apply.setVisible(legacy);discard.setVisible(legacy);
    recover.setVisible(!legacy);recover.setEnabled(operation&&operation->canRecover(session));
    boundaries.setVisible(!legacy&&operation&&operation->advisory);boundaries.setEnabled(operation&&!operation->active()&&operation->ownedBy(session));
    alternatives.setVisible(!legacy&&operation&&operation->choices.size()>1);
    alternatives.setEnabled(operation&&operation->canChoose(session));
    if(operation)
    {
        bool changed=alternatives.getNumItems()!=int(operation->choices.size());
        for(size_t n=0;!changed&&n<operation->choices.size();++n)changed=alternatives.getItemText(int(n))!=operation->choices[n].label;
        if(changed){alternatives.clear(juce::dontSendNotification);for(size_t n=0;n<operation->choices.size();++n)alternatives.addItem(operation->choices[n].label,int(n+1));}
        alternatives.setSelectedId(operation->choices.empty()?0:int(operation->selectedChoice+1),juce::dontSendNotification);
    }
    if(normal){title=operation->title;detail=operation->detail;}
    if(operation&&!legacy)actions=blocked.isEmpty()?"Propose assesses, applies an eligible measured pass, then independently verifies Current. Playback starts only with Play.":"Propose unavailable. "+blocked;
    if(normal&&operation->choices.size()>1)actions+="  Alternatives use the same pre-Propose baseline. Order does not rank sound quality.";
    headline.setText(title,juce::dontSendNotification);headline.setTooltip(title);
    if(reason.getText()!=detail){reason.setText(detail,false);reason.moveCaretToTop(false);}
    availability.setText(actions,juce::dontSendNotification);availability.setTooltip(actions);resized();
}
void ProposalWorkflowPanel::resized()
{
    auto r=getLocalBounds();auto row=r.removeFromTop(32);
    propose.setBounds(row.removeFromLeft(112).reduced(2));review.setBounds(row.removeFromLeft(142).reduced(2));
    apply.setBounds(row.removeFromLeft(132).reduced(2));discard.setBounds(row.removeFromLeft(142).reduced(2));
    recover.setBounds(apply.getX(),apply.getY(),206,apply.getHeight());
    preferences.setBounds(row.removeFromRight(146).reduced(2));history.setBounds(row.removeFromRight(204).reduced(2));
    row=r.removeFromTop(34);regions.setBounds(row.removeFromLeft(228).reduced(2));
    if(alternatives.isVisible())alternatives.setBounds(row.removeFromRight(280).reduced(2));
    if(boundaries.isVisible())boundaries.setBounds(row.removeFromRight(200).reduced(2));
    method.setBounds(row);
    headline.setBounds(r.removeFromTop(26));availability.setBounds(r.removeFromBottom(34));reason.setBounds(r.reduced(0,2));
}
}

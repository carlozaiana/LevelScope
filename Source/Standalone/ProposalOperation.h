#pragma once
#include "StandaloneSessionModel.h"

namespace levelscope::standalone
{
// One message-thread operation. Worker futures and device ownership remain in
// the shell. Nothing here is serialized as measurement or selection authority.
class ProposalOperation
{
public:
    enum class Phase { idle,source,current,search,apply,verify,verified,noChange,noPass,cancelled,failed,stale };
    struct Choice {size_t attempt=0;juce::String resultId,attemptId,label;};
    Phase phase=Phase::idle;
    uint64_t generation=0,ownedRevision=0;
    juce::String runtime,configuration,requirements,title,detail;
    bool applied=false,cancelRequested=false,advisory=false;
    ManualEditState original,owned;
    ProposalEditMetadata originalMetadata;
    SectionProposalRequest baseline;
    std::shared_ptr<const SectionProposalResult> result,selected,pending;
    size_t pendingChoice=0;bool pendingAdvisory=false;
    std::vector<Choice> choices;
    size_t selectedChoice=0;

    bool active() const {return phase==Phase::source||phase==Phase::current||phase==Phase::search||phase==Phase::apply||phase==Phase::verify;}
    static juce::String configurationKey(const StandaloneSessionModel& s)
    {
        const auto source=s.measurementContext(MeasurementDomain::source),current=s.measurementContext(MeasurementDomain::currentState);
        const auto text=s.programmeId+";"+sectionInputIdentity(source)+";"+sectionInputIdentity(current)+";"+source.sourceFile.getFullPathName()
            +";"+source.profile.digest()+";"+s.jobRequirements.digest()+";"+s.proposalPolicy.digest()+";"+juce::String(int(s.programmeApplicability))
            +";"+(s.regionPlanning?s.regionPlanning->digest():juce::String())+";"+source.buildVersion+";"+source.platformAndFloatingPointMode;
        return juce::SHA256(text.toRawUTF8(),text.getNumBytesAsUTF8()).toHexString();
    }
    static juce::String blocked(const StandaloneSessionModel& s)
    {
        if(s.proposalOperationBusy)return "Finish or cancel the current Propose operation.";
        if(s.persistenceLoading||s.job.isActive())return "Finish or cancel the current job.";
        if(s.pendingBoundaryEdits())return s.pendingBoundaryMessage();
        if(s.auditionBusy&&s.auditionBusy())return "Stop playback before proposing.";
        if(s.sectionTrimDraft)return "Apply or discard the pending manual gain draft.";
        if(s.sourceContentInvalid||!s.sourceDocument.matchesFileStamps())return "Resolve the programme audio before proposing.";
        const auto input=OfflineAnalysisService::validateRequest(s.measurementContext(MeasurementDomain::source),MeasurementDomain::source);
        if(input.isNotEmpty())return input;
        if(!s.currentStateDocument.isInitialized&&!s.awaitsInitialCurrent())return "Initialise Current before proposing.";
        return {};
    }
    bool ownedBy(const StandaloneSessionModel& s) const
    {return runtime==s.sessionId&&ownedRevision==s.revision&&owned==s.manualState()&&configuration==configurationKey(s)&&s.sourceDocument.matchesFileStamps()&&!s.sourceContentInvalid&&!s.sectionTrimDraft&&!s.pendingBoundaryEdits();}
    bool canRecover(const StandaloneSessionModel& s) const
    {return phase!=Phase::stale&&applied&&!active()&&ownedBy(s)&&!s.audioEditsBlocked()&&!s.job.isActive();}
    bool canChoose(const StandaloneSessionModel& s) const
    {return canRecover(s)&&choices.size()>1&&s.proposalOperationMatches(baseline,owned,ownedRevision,requirements);}
    void finish(StandaloneSessionModel& s,Phase next,const juce::String& heading,const juce::String& message)
    {
        if(runtime==s.sessionId)s.proposalOperationBusy=false;
        pending.reset();phase=next;title=heading;detail=message;
        if(next==Phase::stale)choices.clear();
    }
    bool begin(StandaloneSessionModel& s)
    {
        const auto error=blocked(s);if(error.isNotEmpty())return false;
        const auto nextGeneration=generation+1;*this=ProposalOperation();generation=nextGeneration;
        s.synchronizeSections();runtime=s.sessionId;configuration=configurationKey(s);requirements=s.jobRequirements.digest();
        original=owned=s.manualState();ownedRevision=s.revision;originalMetadata=s.proposalEditMetadata();
        s.discardSectionProposal();s.proposalOperationBusy=true;prepare(s);return true;
    }
    void prepare(StandaloneSessionModel& s)
    {
        const auto scalar=s.sectionProposalRequest().scalar();
        if(!s.source.hasMeasurement||!scalar.source.currentIdentityMatches)
        {phase=Phase::source;title="Preparing Propose - measuring Original";detail="Measuring the confirmed programme. Your applied Current processing is preserved.";}
        else if(!s.currentState.hasMeasurement||s.currentStateNeedsRemeasurement()||!scalar.current.currentIdentityMatches)
        {phase=Phase::current;title="Preparing Propose - remeasuring Current";detail="Measuring your actual applied gains and transitions before proposing. Manual edits are retained.";}
        else
        {baseline=s.sectionProposalRequest();phase=Phase::search;title="Finding a proposal";detail="Using the current processing and selected method/preferences. Current remains unchanged during the bounded search.";}
    }
    void observe(StandaloneSessionModel& s)
    {
        if(phase==Phase::idle||phase==Phase::stale)return;
        if(runtime!=s.sessionId){const auto next=generation+1;*this=ProposalOperation();generation=next;return;}
        if(s.pendingBoundaryEdits())return; // Escape can still restore the unchanged plan.
        if(!ownedBy(s))finish(s,Phase::stale,"Current or proposal inputs changed", "Previous alternatives are historical and cannot replace your edits. Re-measure or choose Propose again when ready.");
    }
    void measurementComplete(StandaloneSessionModel& s,bool accepted,const juce::String& message)
    {
        if(!accepted||cancelRequested)
        {finish(s,cancelRequested?Phase::cancelled:Phase::failed,applied?"Proposal applied - Current unverified":cancelRequested?"Propose cancelled - processing unchanged":"Propose preparation did not complete",
            message+(applied?"\nApplied processing is retained. Re-measure Current or Restore before Propose; no automatic rollback.":"\nProcessing was not changed."));return;}
        // Only a current-domain measurement accepted by the existing service can
        // supply the final verdict. Candidate receipts never enter live Current.
        if(phase==Phase::verify)
        {
            const auto evaluation=s.evaluation(MeasurementDomain::currentState);
            const bool pass=evaluation.verdict==proposal::Verdict::pass;
            finish(s,Phase::verified,pass?(advisory?"Measured requirements passed - check region boundaries":"Current independently verified - measured requirements passed"):
                "Current independently measured - requirements not passed",
                evaluationText(evaluation)+(advisory?"\n"+boundaryAdvice():juce::String())+"\nListening review and finalized-export verification remain separate.");
        }
        else prepare(s);
    }
    bool candidateComplete(StandaloneSessionModel& s,const SectionProposalOutcome& outcome)
    {
        if(cancelRequested){finish(s,Phase::cancelled,"Propose cancelled - processing unchanged","No proposal was applied.");return false;}
        if(!outcome.value)
        {finish(s,Phase::failed,"Proposal did not complete",outcome.message+"\nProcessing remains unchanged. This does not prove that no solution exists.");return false;}
        result=*outcome.value;
        if(result->identity()!=baseline.binding()){finish(s,Phase::stale,"Proposal inputs changed","The completed result cannot apply to this operation.");return false;}
        if(result->status()==SectionProposalStatus::noChange)
        {
            const bool pass=s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass;
            finish(s,pass?Phase::noChange:Phase::failed,pass?"Current already meets the measured requirements":"Current result could not be verified",result->outcomeReason());return false;
        }
        choices.clear();bool preferredFound=false;
        if(const auto& review=result->regionReview())
        {
            const auto preferred=review->recommendedAttempt.value_or(review->retainedAttempt);
            for(size_t n=0;n<review->attempts.size();++n)
            {
                const auto& attempt=review->attempts[n];if(!attempt.receipt)continue;
                const auto checked=selectWorkflowProposal(baseline,result,review->resultIdentity,attempt.attemptIdentity,false);if(!checked.value)continue;
                juce::String label="Attempt "+juce::String(n+1);
                if(attempt.measured[2].value)label+=" | LRA "+juce::String(*attempt.measured[2].value,2)+" LU";
                if(attempt.measured[1].value)label+=" | TP "+juce::String(*attempt.measured[1].value,2)+" dBTP";
                choices.push_back({n,review->resultIdentity,attempt.attemptIdentity,label});
                if(n==preferred){selectedChoice=choices.size()-1;preferredFound=true;}
            }
        }
        else if(measuredProposalApplyError(baseline,*result).isEmpty()){choices.push_back({0,{},{},"Measured proposed plan"});preferredFound=true;}
        if(choices.empty())
        {finish(s,Phase::noPass,"No eligible measured pass found in this search",result->outcomeTitle()+"\n"+result->outcomeReason()+"\nCurrent processing is unchanged. Search/preference limits do not prove that no solution exists.");return false;}
        if(!preferredFound){finish(s,Phase::failed,"Recommended plan is not eligible","The existing planner's selected plan has no valid passing receipt. Current is unchanged; inspect the complete report.");return false;}
        return select(s,selectedChoice,true);
    }
    bool select(StandaloneSessionModel& s,size_t index,bool automatic)
    {
        if(index>=choices.size()||!ownedBy(s)||(!automatic&&!canChoose(s)))return false;
        const auto& c=choices[index];const auto chosen=selectWorkflowProposal(baseline,result,c.resultId,c.attemptId,automatic);
        if(!chosen.value){finish(s,Phase::failed,"Measured alternative could not be selected",chosen.message);return false;}
        pending=*chosen.value;pendingChoice=index;cancelRequested=false;s.proposalOperationBusy=true;phase=Phase::apply;
        pendingAdvisory=false;if(const auto& review=pending->regionReview();review&&review->eventMode&&review->selectedAttempt)
        {const auto& a=review->attempts[*review->selectedAttempt];pendingAdvisory=!a.plan.assignments.empty()&&(!a.coverage.complete()||review->searchTruncated);}
        title="Measured candidate passes - applying proposal";detail=pending->outcomeReason()+"\nCurrent still needs independent verification after application.";return true;
    }
    bool apply(StandaloneSessionModel& s)
    {
        if(cancelRequested){finish(s,Phase::cancelled,applied?"Alternative switch cancelled":"Propose cancelled - processing unchanged","The pending choice was not applied.");return false;}
        if(!ownedBy(s)||!s.applyWorkflowProposal(baseline,result,pending,owned,ownedRevision,requirements))
        {finish(s,Phase::failed,"Proposal not applied",s.sectionTrimNotice+"\nProcessing was retained. Resolve the reported cause before proposing again.");return false;}
        selected=std::move(pending);selectedChoice=pendingChoice;advisory=pendingAdvisory;
        owned=s.manualState();ownedRevision=s.revision;applied=true;phase=Phase::verify;
        title="Proposal applied - independently verifying Current";detail="The candidate passed privately. Current remains unverified until this new measurement completes. Cancel retains the applied plan and recovery.";return true;
    }
    bool recover(StandaloneSessionModel& s)
    {
        if(!canRecover(s)||!s.commitProposalProcessing(original,originalMetadata,"Restore before Propose"))return false;
        owned=s.manualState();ownedRevision=s.revision;choices.clear();applied=false;
        finish(s,Phase::stale,"Processing restored to before Propose","Prior gains and explanation restored. Current needs independent remeasurement; no old PASS or proposal authority was restored.");return true;
    }
    juce::String boundaryAdvice() const
    {
        if(!selected||!selected->regionReview()||!selected->regionReview()->selectedAttempt)return {};
        const auto& review=*selected->regionReview();const auto& attempt=review.attempts[*review.selectedAttempt];
        juce::String text=attempt.coverage.notice;
        size_t shown=0;for(const auto& group:attempt.coverage.groups)if(group.knownPartial||group.incomplete)
        {if(shown++<6)text+="\n"+juce::String(double(group.extent.begin)/baseline.currentContext.media.sampleRate,3)+"-"+juce::String(double(group.extent.end)/baseline.currentContext.media.sampleRate,3)+" s: "+group.reason;}
        return text+"\nInspect the applied region handles and transition bands. Previous/next in the gain lane reaches every region; the full explanation retains all boundary evidence. Audition the entrances and returns; edit if needed, then remeasure.";
    }
};
}

#include "StandaloneSectionProposal.h"
#include "StandaloneProcessing.h"
#include "StandaloneSections.h"
#include "AppliedProposalExplanation.h"
#include <juce_cryptography/juce_cryptography.h>

namespace levelscope::standalone
{
proposal::TrimRequest SectionProposalRequest::scalar() const
{ auto r=trimRequest(source,current,sourceContext,currentContext,currentContext.profile,applicability);r.preferenceIdentity=policy.digest();r.preferencesAvailable=policy.validationError().isEmpty();return r; }
juce::String SectionProposalRequest::binding() const
{
    const auto s=proposal::draftGlobalTrim(scalar()).identity+";"+proposal::sectionTrimMethodVersion
        +";"+proposal::sectionGeometryKey(annotations)+";"+proposal::sectionGeometryKey(appliedAnnotations)
        +";"+currentContext.sessionId+";"+juce::String(currentContext.ticket.sessionRevision)
        +";"+juce::String(int(pendingManualDraft))+";"+policy.digest();
    return juce::SHA256(s.toRawUTF8(),s.getNumBytesAsUTF8()).toHexString();
}
std::optional<ProcessingPlan> SectionProposalResult::previewPlan() const
{
    if(state==SectionProposalStatus::ready)return proposed;
    if(state==SectionProposalStatus::globalOnly&&scalarResult.totalGainDb){auto p=proposed;p.staticGainDb=*scalarResult.totalGainDb;return p;}
    if(candidate)return candidate->context.processing;
    if(!seed.plan.intervals.empty()&&seed.plan.validationError().isEmpty())
    {auto p=proposed;p.sectionGains=seed.plan.canonical();p.schemaVersion=2;p.engineVersion=processing::SectionGainPlan::version;return p;}
    return {};
}
SectionProposalOutcome generateSectionProposal(const SectionProposalRequest& r,JobControl& control,const SnapshotStorage& storage)
{
    auto out=std::shared_ptr<SectionProposalResult>(new SectionProposalResult);
    out->bound=r.binding();out->proposed=r.currentContext.processing;out->annotations=r.annotations;
    const auto complete=[&]() -> SectionProposalOutcome
    {
        if(control.isCancelled())return {{},ErrorCode::cancelled,"Proposal cancelled; applied settings and measurements preserved."};
        out->narrative=out->explanation;
        out->title=out->state==SectionProposalStatus::ready||out->state==SectionProposalStatus::globalOnly?"Ready to apply":
            out->state==SectionProposalStatus::noChange?"Already meets requirements":
            out->problem==SectionProposalIssue::unsupportedRequirements?"Unsupported correction":
            out->problem==SectionProposalIssue::preferences?"Blocked by proposal preferences":
            out->problem==SectionProposalIssue::candidateRejected?"Candidate rejected":"Proposal unavailable";
        out->reason=out->state==SectionProposalStatus::globalOnly?"Global gain alone has a predicted solution. Apply explicitly, then Re-measure Current State.":
            out->state==SectionProposalStatus::ready?"The measured private candidate passes the selected requirements. Apply explicitly, then Re-measure Current State.":
            out->state==SectionProposalStatus::noChange?"Measured Current already passes every selected requirement. No Apply or remeasurement is needed.":out->explanation;
        if(out->problem==SectionProposalIssue::unsupportedRequirements)out->reason=out->scalarResult.explanation;
        if(out->seed.reason.isNotEmpty()&&!out->seed.available&&out->state==SectionProposalStatus::noSolution)out->reason=out->seed.reason;
        if(out->state==SectionProposalStatus::noSolution&&!out->seed.available&&out->seed.affectedFrames>=0&&r.annotations.frames>0)
        {
            out->reason="Local changes affect "+juce::String(100.*double(out->seed.affectedFrames)/r.annotations.frames,2)
                +"% of the programme against a "+juce::String(r.policy.durationThousandths/10.,1)+"% duration allowance. "+out->reason;
            auto longest=out->seed.corrections;
            std::sort(longest.begin(),longest.end(),[](const auto& a,const auto& b){return a.end-a.begin>b.end-b.begin;});
            out->reason+=" Largest selected sections (section lengths, excluding transition overlap): ";
            for(size_t n=0;n<std::min<size_t>(3,longest.size());++n)
            {
                const auto& c=longest[n];const auto it=std::lower_bound(r.annotations.sections.begin(),r.annotations.sections.end(),c.begin,[](const auto& a,int64_t f){return a.begin<f;});
                out->reason+=(n?", ":"")+(it!=r.annotations.sections.end()?it->name:juce::String(c.id))+" "+juce::String(double(c.end-c.begin)/r.annotations.sampleRate,3)+" s";
            }
            out->reason+=". Review all changes, edit sections or proposal limits, then Propose again.";
        }
        if(out->state==SectionProposalStatus::globalOnly)
            out->reason="Proposed global gain "+juce::String(*out->scalarResult.totalGainDb,9)+" dB. Predicted IL "+juce::String(*out->scalarResult.predictedIntegrated,6)+" LUFS / TP "+juce::String(*out->scalarResult.predictedTruePeak,6)+" dBTP. "+out->reason;
        if(out->problem==SectionProposalIssue::candidateRejected&&out->candidate)
        {
            const auto& d=out->candidate->detail;const auto& l=r.currentContext.profile.limits;
            out->reason="Measured private local/combined candidate IL "+(d.integrated.value?juce::String(*d.integrated.value,6):juce::String("unavailable"))+" LUFS against "+juce::String(l.integratedLufs,1)+" +/- "+juce::String(l.integratedToleranceLu,1)+" LU. TP "+(d.truePeak.value?juce::String(*d.truePeak.value,6):juce::String("unavailable"))+" dBTP against "+juce::String(l.truePeakCeilingDbtp,1)+" dBTP. "
                +"This bounded attempt found no passing combined plan. Review the measured results and uncovered deficits, then revise sections, limits or manual processing and remeasure.";
        }
        const char* outcome=out->state==SectionProposalStatus::ready?"READY FOR EXPLICIT APPLY - exact private candidate passes selected rules; Current unchanged":
            out->state==SectionProposalStatus::noChange?"NO CHANGE - measured Current already passes":
            out->state==SectionProposalStatus::globalOnly?"GLOBAL-ONLY DRAFT - prediction, explicit Apply and remeasurement required":
            out->state==SectionProposalStatus::noSolution?"NO SOLUTION FOUND WITHIN THESE PREFERENCES":"UNAVAILABLE";
        out->explanation=out->title+"\n"+out->reason+"\n\n"+juce::String(outcome)+"\nSection-led / "+juce::String(proposal::sectionTrimMethodVersion)+"\nBaseline binding: "+out->bound
            +"\nProfile: "+r.currentContext.profile.id+" / "+r.currentContext.profile.definitionVersion+" / "+r.currentContext.profile.digest()
            +"\n"+r.policy.text()+"\n"+out->explanation+"\nPrivate full-programme passes: "+juce::String(out->passCount)+" / 3.\n";
        if(!out->seed.reason.isEmpty())out->explanation+=proposal::sectionSeedText(out->seed)+"\n";
        out->explanation+=proposalTargetText(r.currentContext.profile)+"\n";
        if(r.source)out->explanation+=proposalMeasurementText(proposalMeasurementRecord(*r.source,r.currentContext.profile,r.applicability),"Source baseline")+"\n";
        if(r.current)out->explanation+=proposalMeasurementText(proposalMeasurementRecord(*r.current,r.currentContext.profile,r.applicability),"Current baseline before proposal")+"\n";
        if(out->candidate)out->explanation+=proposalMeasurementText(proposalMeasurementRecord(*out->candidate,r.currentContext.profile,r.applicability),"Measured private candidate");
        else out->explanation+="No private candidate measurement was performed.\n";
        control.setProgress(1);return {std::shared_ptr<const SectionProposalResult>(out),ErrorCode::none,"Proposal complete; no processing applied."};
    };
    try
    {
        if(r.policy.validationError().isNotEmpty()){out->problem=SectionProposalIssue::preferences;out->explanation=r.policy.validationError();return complete();}
        const auto scalar=r.scalar();out->scalarResult=proposal::draftGlobalTrim(scalar);
        out->explanation=proposal::trimProposalText(out->scalarResult);
        if(control.isCancelled())return complete();
        if(out->scalarResult.status==proposal::TrimStatus::unavailable)
        {
            if(out->scalarResult.explanation.startsWith("Automatic correction"))out->problem=SectionProposalIssue::unsupportedRequirements;
            return complete();
        }
        if(!r.sourceContext.sourceFilesMatch()||!r.currentContext.sourceFilesMatch())
            return {{},ErrorCode::sourceChanged,"Source members changed; proposal discarded."};
        if(out->scalarResult.status==proposal::TrimStatus::noChange)
        {out->state=SectionProposalStatus::noChange;out->problem=SectionProposalIssue::none;return complete();}
        if(r.pendingManualDraft){out->explanation="Resolve the unapplied manual draft explicitly before generating a changing proposal.";return complete();}
        if(out->scalarResult.status==proposal::TrimStatus::draft)
        {out->state=SectionProposalStatus::globalOnly;out->problem=SectionProposalIssue::none;return complete();}
        const auto& map=r.annotations;
        if(!map.valid()||map.inputIdentity!=sectionInputIdentity(r.currentContext))
        {out->explanation="Valid current programme annotations are required. No sections generated automatically.";return complete();}
        if(r.currentContext.processing.sectionGains
            &&proposal::sectionGeometryKey(r.appliedAnnotations)!=proposal::sectionGeometryKey(map))
        {out->explanation="Annotation geometry differs from the applied envelope provenance. Manual gains cannot be remapped.";return complete();}
        auto probe=r.currentContext;
        processing::SectionGainPlan plan;plan.inputIdentity=map.inputIdentity;plan.sampleRate=map.sampleRate;plan.frames=map.frames;plan.intervals={{0,map.frames,-1}};
        probe.processing.sectionGains=plan;probe.processing.schemaVersion=2;probe.processing.engineVersion=processing::SectionGainPlan::version;
        const auto error=OfflineAnalysisService::validateRequest(probe,MeasurementDomain::privateCandidate);
        if(error.isNotEmpty()){out->explanation=error;return complete();}
        JobControl prepare;prepare.parent=&control;prepare.parentProgress=&control.progress;prepare.progressScale=.1f;
        const auto copy=prepareVerifiedSource(r.currentContext,prepare,storage);
        if(!copy.value)return {{},copy.error,copy.message};
        OfflineAnalysisService worker;
        auto context=r.currentContext;
        const auto measure=[&](measurement::PeakObserver* observer=nullptr)
        {
            JobControl pass;pass.parent=&control;pass.parentProgress=&control.progress;
            pass.progressBase=.1f+float(out->passCount)*.3f;pass.progressScale=.3f;
            if(out->passCount>=3)return ServiceResult<AnalysisSnapshot>{{},ErrorCode::invalidRequest,"Private pass budget exhausted"};
            ++out->passCount;
            return worker.analyzeCandidate(context,**copy.value,pass,observer);
        };
        if(r.attribution&&r.attribution->snapshotIdentity==scalar.current.snapshotIdentity
            &&r.attribution->geometry==proposal::sectionGeometryKey(map)
            &&r.attribution->channels==r.currentContext.fileLayout.count)out->attribution=r.attribution;
        else
        {
            proposal::SectionPeakCollector collector(map,r.currentContext.fileLayout.count,scalar.current.snapshotIdentity);
            const auto baseline=measure(&collector);
            if(!baseline.value)return {{},baseline.error,baseline.message};
            // This is the identical processing stream; mismatches are an identity/numerical failure, not a correction opportunity.
            const auto& a=baseline.value->detail;const auto& b=r.current->detail;
            if(a.integrated.value!=b.integrated.value||a.truePeak.value!=b.truePeak.value||a.samplePeak.value!=b.samplePeak.value
                ||a.channelTruePeaks!=b.channelTruePeaks||a.channelSamplePeaks!=b.channelSamplePeaks)
                return {{},ErrorCode::staleRevision,"Private baseline disagrees with measured Current. Remeasure; no proposal offered."};
            out->attribution=collector.finish();
        }
        out->seed=proposal::selectSectionTrims(scalar,map,r.currentContext.processing.sectionGains,*out->attribution,r.policy);
        out->state=SectionProposalStatus::noSolution;out->problem=SectionProposalIssue::preferences;
        if(!out->seed.available)return complete();
        out->problem=SectionProposalIssue::candidateRejected;
        context.processing.sectionGains=out->seed.plan.canonical();context.processing.schemaVersion=2;
        context.processing.engineVersion=processing::SectionGainPlan::version;
        auto local=measure();if(!local.value)return {{},local.error,local.message};
        out->candidate=std::make_shared<const AnalysisSnapshot>(std::move(*local.value));
        auto adjusted=scalar;adjusted.current=measuredProgrammeForSnapshot(out->candidate.get(),context,MeasurementDomain::privateCandidate,r.applicability);
        adjusted.processingVersion=context.processing.engineVersion;
        const auto key=context.processing.sectionGains->audioKey();
        adjusted.fixedUpstreamPlanIdentity=juce::SHA256(key.toRawUTF8(),key.getNumBytesAsUTF8()).toHexString();
        const auto finalGlobal=proposal::draftGlobalTrimForCandidate(adjusted);
        out->explanation+="Private local-plan global calculation: "+finalGlobal.explanation+"\n";
        if(finalGlobal.status!=proposal::TrimStatus::draft&&finalGlobal.status!=proposal::TrimStatus::noChange)
        {out->explanation+="NO SOLUTION FOUND WITHIN THIS POLICY. One envelope assessed; no retry.";return complete();}
        context.processing.staticGainDb=*finalGlobal.totalGainDb;
        if(context.processing.staticGainDb!=r.currentContext.processing.staticGainDb)
        {
            auto final=measure();if(!final.value)return {{},final.error,final.message};
            out->candidate=std::make_shared<const AnalysisSnapshot>(std::move(*final.value));
        }
        const auto verdict=evaluateSnapshot(out->candidate.get(),context,MeasurementDomain::privateCandidate,r.currentContext.profile,r.applicability);
        if(verdict.verdict!=proposal::Verdict::pass)
        {out->explanation+="Exact combined candidate did not pass. Actual measured result shown; no second envelope or automatic retry.";return complete();}
        out->state=SectionProposalStatus::ready;out->problem=SectionProposalIssue::none;out->proposed=context.processing;
        out->explanation+="READY FOR EXPLICIT APPLY: private exact candidate passes selected programme rules. Current State remains unchanged.\n"
            "After Apply Current is UNVERIFIED until independent remeasurement. Global gain "+juce::String(context.processing.staticGainDb,12)+" dB affects the entire programme.\n";
        double minimumLocal=0,maximumLocal=0;
        if(context.processing.sectionGains){minimumLocal=12;maximumLocal=-12;for(const auto& row:context.processing.sectionGains->intervals){minimumLocal=std::min(minimumLocal,row.db);maximumLocal=std::max(maximumLocal,row.db);}}
        out->explanation+="Proposed local gain range "+juce::String(minimumLocal,12)+" to "+juce::String(maximumLocal,12)+" dB; combined local-plus-global range "+juce::String(minimumLocal+context.processing.staticGainDb,12)+" to "+juce::String(maximumLocal+context.processing.staticGainDb,12)+" dB. Local attenuation preferences do not cap whole-chain attenuation.\n";
        return complete();
    }
    catch(const std::exception& e){return {{},ErrorCode::invalidRequest,"Proposal failed safely: "+juce::String(e.what())};}
}
}

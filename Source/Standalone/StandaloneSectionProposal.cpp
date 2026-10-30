#include "StandaloneSectionProposal.h"
#include "StandaloneProcessing.h"
#include "StandaloneSections.h"
#include <juce_cryptography/juce_cryptography.h>

namespace levelscope::standalone
{
proposal::TrimRequest SectionProposalRequest::scalar() const
{ auto r=trimRequest(source.get(),current.get(),sourceContext,currentContext,currentContext.profile,applicability);r.preferenceIdentity=policy.digest();r.preferencesAvailable=policy.validationError().isEmpty();return r; }
juce::String SectionProposalRequest::binding() const
{
    const auto s=proposal::draftGlobalTrim(scalar()).identity+";"+proposal::sectionTrimMethodVersion
        +";"+proposal::sectionGeometryKey(annotations)+";"+proposal::sectionGeometryKey(appliedAnnotations)
        +";"+currentContext.sessionId+";"+juce::String(currentContext.ticket.sessionRevision)
        +";"+juce::String(int(pendingManualDraft))+";"+policy.digest();
    return juce::SHA256(s.toRawUTF8(),s.getNumBytesAsUTF8()).toHexString();
}
SectionProposalOutcome generateSectionProposal(const SectionProposalRequest& r,JobControl& control,const SnapshotStorage& storage)
{
    auto out=std::shared_ptr<SectionProposalResult>(new SectionProposalResult);
    out->bound=r.binding();out->proposed=r.currentContext.processing;
    const auto complete=[&]() -> SectionProposalOutcome
    {
        if(control.isCancelled())return {{},ErrorCode::cancelled,"Proposal cancelled; applied settings and measurements preserved."};
        const char* outcome=out->state==SectionProposalStatus::ready?"READY FOR EXPLICIT APPLY - exact private candidate passes selected rules; Current unchanged":
            out->state==SectionProposalStatus::noChange?"NO CHANGE - measured Current already passes":
            out->state==SectionProposalStatus::globalOnly?"GLOBAL-ONLY DRAFT - prediction, explicit Apply and remeasurement required":
            out->state==SectionProposalStatus::noSolution?"NO SOLUTION FOUND WITHIN THESE PREFERENCES":"UNAVAILABLE";
        out->explanation=juce::String(outcome)+"\nSection-led / "+juce::String(proposal::sectionTrimMethodVersion)+"\nBaseline binding: "+out->bound
            +"\nProfile: "+r.currentContext.profile.id+" / "+r.currentContext.profile.definitionVersion+" / "+r.currentContext.profile.digest()
            +"\n"+r.policy.text()+"\n"+out->explanation+"\nPrivate full-programme passes: "+juce::String(out->passCount)+" / 3.\n";
        if(r.current) out->explanation+="Baseline measured I / TP: "+(r.current->measured.integratedLufs?juce::String(*r.current->measured.integratedLufs,9):juce::String("unavailable"))
            +" LUFS / "+(r.current->measured.truePeakDbtp?juce::String(*r.current->measured.truePeakDbtp,9):juce::String("unavailable"))+" dBTP.\n";
        if(out->candidate)out->explanation+=evaluationText(evaluateSnapshot(out->candidate.get(),out->candidate->context,
            MeasurementDomain::privateCandidate,r.currentContext.profile,r.applicability));
        else out->explanation+="No private candidate measurement was performed.\n";
        control.setProgress(1);return {std::shared_ptr<const SectionProposalResult>(out),ErrorCode::none,"Proposal complete; no processing applied."};
    };
    try
    {
        if(r.policy.validationError().isNotEmpty()){out->explanation=r.policy.validationError();return complete();}
        const auto scalar=r.scalar();out->scalarResult=proposal::draftGlobalTrim(scalar);
        out->explanation=proposal::trimProposalText(out->scalarResult);
        if(control.isCancelled())return complete();
        if(out->scalarResult.status==proposal::TrimStatus::unavailable)return complete();
        if(!r.sourceContext.sourceFilesMatch()||!r.currentContext.sourceFilesMatch())
            return {{},ErrorCode::sourceChanged,"Source members changed; proposal discarded."};
        if(out->scalarResult.status==proposal::TrimStatus::noChange)
        {out->state=SectionProposalStatus::noChange;return complete();}
        if(r.pendingManualDraft){out->explanation="Resolve the unapplied manual draft explicitly before generating a changing proposal.";return complete();}
        if(out->scalarResult.status==proposal::TrimStatus::draft)
        {out->state=SectionProposalStatus::globalOnly;return complete();}
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
        out->state=SectionProposalStatus::noSolution;out->explanation+=proposal::sectionSeedText(out->seed);
        if(!out->seed.available)return complete();
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
        out->state=SectionProposalStatus::ready;out->proposed=context.processing;
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

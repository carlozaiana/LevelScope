#include "StandaloneLraProposal.h"
#include "AppliedProposalExplanation.h"
#include "StandaloneProcessing.h"
#include <limits>

namespace levelscope::standalone
{
namespace
{
bool finite(const measurement::Metric& m){return m.status==measurement::Status::available&&m.value&&std::isfinite(*m.value);}
}
std::optional<double> lraFinalGlobalGain(const AnalysisSnapshot& a,const proposal::TargetProfile& p)
{
    const auto& d=a.detail;const auto& l=p.limits;const double before=a.context.processing.staticGainDb;
    if(p.validationError().isNotEmpty()||!l.hasIntegratedLoudnessTarget||!l.hasTruePeakCeiling
        ||!finite(d.integrated)||!finite(d.truePeak)||(l.hasMaxMomentary&&!finite(d.momentaryMax))
        ||(l.hasMaxShortTerm&&!finite(d.shortTermMax))||!std::isfinite(before))return {};
    const double lower=std::max(-60-before,l.integratedLufs-l.integratedToleranceLu-*d.integrated.value);
    double upper=std::min({24-before,l.integratedLufs+l.integratedToleranceLu-*d.integrated.value,l.truePeakCeilingDbtp-*d.truePeak.value});
    if(l.hasMaxMomentary)upper=std::min(upper,l.maxMomentaryLufs-*d.momentaryMax.value);
    if(l.hasMaxShortTerm)upper=std::min(upper,l.maxShortTermLufs-*d.shortTermMax.value);
    if(lower>upper)return {};
    double gain=before+std::clamp(l.integratedLufs-*d.integrated.value,lower,upper);
    for(int n=0;n<16;++n)
    {
        const double delta=gain-before,i=*d.integrated.value+delta;
        const bool low=gain< -60||i<l.integratedLufs-l.integratedToleranceLu;
        const bool high=gain>24||i>l.integratedLufs+l.integratedToleranceLu||*d.truePeak.value+delta>l.truePeakCeilingDbtp
            ||(l.hasMaxMomentary&&*d.momentaryMax.value+delta>l.maxMomentaryLufs)
            ||(l.hasMaxShortTerm&&*d.shortTermMax.value+delta>l.maxShortTermLufs);
        if(!low&&!high)return gain;
        if(low==high)break;
        gain=std::nextafter(gain,low?std::numeric_limits<double>::infinity():-std::numeric_limits<double>::infinity());
    }
    return {};
}
std::optional<LraAttemptRank> lraAttemptRank(const AnalysisSnapshot& a,const proposal::Evaluation& e,int64_t frames,double maximum,int index)
{
    if(e.verdict==proposal::Verdict::notEvaluated||e.verdict==proposal::Verdict::indeterminate||!finite(a.detail.lra))return {};
    double excess=0,count=0,total=0;bool found=false;
    for(const auto& rule:e.rules)
    {
        if(!rule.measured||!rule.upper||!std::isfinite(*rule.measured))return {};
        const double above=std::max(0.,*rule.measured-*rule.upper),below=rule.lower?std::max(0.,*rule.lower-*rule.measured):0.;
        if(rule.id=="loudness-range"){found=true;excess=above;}
        else {count+=rule.verdict==proposal::Verdict::fail?1:0;total+=above+below;}
    }
    if(!found)return {};
    return LraAttemptRank{e.verdict==proposal::Verdict::pass?0.:1.,excess,count,total,double(frames),maximum,double(index)};
}
SectionProposalOutcome generateLraProposal(const SectionProposalRequest& r,JobControl& control,const SnapshotStorage& storage)
{
    auto out=std::shared_ptr<SectionProposalResult>(new SectionProposalResult);
    out->bound=r.binding();out->method=proposal::lraTrimMethodVersion;out->proposed=r.currentContext.processing;
    out->annotations=r.annotations;out->seed.policy=r.policy;
    const auto complete=[&]() -> SectionProposalOutcome
    {
        if(control.isCancelled())return {{},ErrorCode::cancelled,"LRA assessment cancelled; applied audio and measurements preserved."};
        if(r.current&&!r.currentContext.sourceFilesMatch())return {{},ErrorCode::sourceChanged,"Source changed during LRA assessment; candidate discarded."};
        out->title=out->state==SectionProposalStatus::ready?"Measured LRA candidate ready":out->state==SectionProposalStatus::noChange?"Already meets requirements":
            out->state==SectionProposalStatus::noSolution?"No acceptable candidate in this bounded attempt":"LRA correction unavailable";
        out->reason=out->explanation;
        out->narrative=out->title+"\n"+out->explanation+"\n"+out->narrative;
        out->explanation=out->narrative+"\nMethod "+out->method+" / binding "+out->bound+"\n"+r.policy.text()+"\n"+proposalTargetText(r.currentContext.profile);
        if(r.source)out->explanation+=proposalMeasurementText(proposalMeasurementRecord(*r.source,r.currentContext.profile,r.applicability),"Source baseline");
        if(r.current)out->explanation+=proposalMeasurementText(proposalMeasurementRecord(*r.current,r.currentContext.profile,r.applicability),"Current baseline");
        out->explanation+="\nPrivate full-programme assessments "+juce::String(out->passCount)+" / 3. Estimates never satisfy measured rules.\n";
        if(out->candidate)out->explanation+=proposalMeasurementText(proposalMeasurementRecord(*out->candidate,r.currentContext.profile,r.applicability),"Retained measured private attempt");
        control.setProgress(1);return {std::shared_ptr<const SectionProposalResult>(out),ErrorCode::none,"LRA assessment complete; Current unchanged."};
    };
    try
    {
        const auto scalar=r.scalar();out->scalarResult=proposal::draftGlobalTrim(scalar);
        if(r.policy.validationError().isNotEmpty()){out->problem=SectionProposalIssue::preferences;out->explanation=r.policy.validationError();return complete();}
        const auto source=proposal::evaluateMeasuredProgramme(r.currentContext.profile,scalar.source);
        const auto baseline=proposal::evaluateMeasuredProgramme(r.currentContext.profile,scalar.current);
        if(!r.source||!r.current||scalar.sourceInputIdentity.isEmpty()||scalar.sourceInputIdentity!=scalar.currentInputIdentity
            ||!scalar.source.currentIdentityMatches||!scalar.current.currentIdentityMatches
            ||source.verdict==proposal::Verdict::notEvaluated||baseline.verdict==proposal::Verdict::notEvaluated
            ||!r.sourceContext.sourceFilesMatch()||!r.currentContext.sourceFilesMatch())
        {out->explanation="Valid matching measured Source/Current, confirmed roles/range and applicability are required. Stale inputs cannot be assessed.";return complete();}
        if(out->scalarResult.status==proposal::TrimStatus::noChange)
        {out->state=SectionProposalStatus::noChange;out->problem=SectionProposalIssue::none;out->explanation="Measured Current passes every selected requirement. Genuine no-op: retain audio revision, settings and measurements.";return complete();}
        const auto& limits=r.currentContext.profile.limits;
        if(!limits.hasMaxLoudnessRange||!limits.hasIntegratedLoudnessTarget||!limits.hasTruePeakCeiling
            ||baseline.verdict==proposal::Verdict::indeterminate||!r.current->detail.lraStable||!finite(r.current->detail.lra)
            ||!finite(r.source->detail.integrated)||!finite(r.source->detail.truePeak))
        {out->problem=SectionProposalIssue::unsupportedRequirements;out->explanation="Method unavailable: select IL, true peak and LRA requirements with stable, available measured metrics. Additional M/S limits must also be measurable. This does not assert that no solution exists.";return complete();}
        if(r.pendingManualDraft){out->explanation="Apply or discard the pending manual draft before a changing assessment.";return complete();}
        if(r.currentContext.fileLayout.count!=2&&r.currentContext.fileLayout.count!=6)
        {out->problem=SectionProposalIssue::unsupportedRequirements;out->explanation="This LRA method supports stereo and conventional 5.1 only; other existing workflows are preserved.";return complete();}
        auto context=r.currentContext;
        processing::SectionGainPlan before;before.inputIdentity=sectionInputIdentity(context);before.sampleRate=int(context.media.sampleRate);before.frames=context.media.frames;
        before.intervals={{0,before.frames,0}};if(context.processing.sectionGains)before=context.processing.sectionGains->canonical();
        auto probe=context;probe.processing.sectionGains=before;probe.processing.sectionGains->intervals={{0,before.frames,-1}};
        probe.processing.schemaVersion=2;probe.processing.engineVersion=processing::SectionGainPlan::version;
        const auto error=OfflineAnalysisService::validateRequest(probe,MeasurementDomain::privateCandidate);
        if(error.isNotEmpty()){out->explanation=error;return complete();}
        const auto cancel=[&]{return control.isCancelled();};
        auto analysis=proposal::analyzeProgramme(r.current->detail,r.current->identity,cancel);if(!analysis)return complete();
        const auto passages=proposal::findLraPassages(analysis,limits.maxLoudnessRangeLu,cancel);if(!passages)return complete();
        const bool globalFirst=*r.current->detail.lra.value<=limits.maxLoudnessRangeLu&&bool(lraFinalGlobalGain(*r.current,r.currentContext.profile));
        proposal::LraTrimSelection selection;
        if(!globalFirst)
        {
            selection=proposal::selectLraTrims(*passages,before,r.policy,cancel);if(control.isCancelled())return complete();out->lraSeed=selection;
            out->narrative=selection.text();
            if(!selection.available)
            {out->state=SectionProposalStatus::noSolution;out->problem=SectionProposalIssue::preferences;out->explanation="No envelope selected. "+selection.text()+"No candidate was measured. Method/policy restrictions do not establish that the programme has no solution.";return complete();}
        }
        JobControl preparing;preparing.parent=&control;preparing.parentProgress=&control.progress;preparing.progressScale=.1f;
        const auto copy=prepareVerifiedSource(r.currentContext,preparing,storage);if(!copy.value)return {{},copy.error,copy.message};
        OfflineAnalysisService worker;std::optional<LraAttemptRank> bestRank;int bestIndex=0;
        std::optional<proposal::LraTrimSelection> bestSelection;
        const auto measure=[&](const std::optional<proposal::LraTrimSelection>& selected)->ServiceResult<AnalysisSnapshot>
        {
            if(out->passCount>=3)return {{},ErrorCode::invalidRequest,"Private LRA assessment budget exhausted"};
            JobControl pass;pass.parent=&control;pass.parentProgress=&control.progress;pass.progressBase=.1f+float(out->passCount)*.3f;pass.progressScale=.3f;++out->passCount;
            auto measured=worker.analyzeCandidate(context,**copy.value,pass);if(!measured.value)return measured;
            const auto evaluation=evaluateSnapshot(&*measured.value,context,MeasurementDomain::privateCandidate,r.currentContext.profile,r.applicability);
            auto rank=lraAttemptRank(*measured.value,evaluation,selected?selected->affectedFrames:0,selected?selected->maximumAdditional:0,out->passCount);
            out->narrative+="\nAssessment "+juce::String(out->passCount)+" exact plan "+proposalPlanIdentity(context.processing)+"; global "+juce::String(context.processing.staticGainDb,12)+" dB.\n";
            out->narrative+=proposalMeasurementText(proposalMeasurementRecord(*measured.value,r.currentContext.profile,r.applicability),"Actual complete private assessment");
            if(rank&&(!bestRank||*rank<*bestRank))
            {
                bestRank=rank;bestIndex=out->passCount;bestSelection=selected;
                out->candidate=std::make_shared<const AnalysisSnapshot>(*measured.value);
            }
            else out->narrative+="Earlier measured attempt retained: this attempt is worse or has unavailable required metrics under the versioned selection rule.\n";
            return measured;
        };
        if(globalFirst)
        {
            context.processing.staticGainDb=*lraFinalGlobalGain(*r.current,r.currentContext.profile);
            out->narrative+="LRA already meets the ceiling; assess explicit global trim first. Prediction alone cannot enable Apply.\n";
            auto result=measure({});if(!result.value)return {{},result.error,result.message};
        }
        else
        {
            context.processing.sectionGains=selection.plan.canonical();context.processing.schemaVersion=2;context.processing.engineVersion=processing::SectionGainPlan::version;
            auto first=measure(selection);if(!first.value)return {{},first.error,first.message};
            if(bestRank&&(*bestRank)[0]!=0&&finite(first.value->detail.lra)&&*first.value->detail.lra.value>limits.maxLoudnessRangeLu)
            {
                const auto refined=proposal::refineLraTrims(selection,*first.value->detail.lra.value-limits.maxLoudnessRangeLu,r.policy,cancel);
                if(control.isCancelled())return complete();
                out->narrative+="\nRefinement estimate: "+refined.text();
                if(refined.available){context.processing.sectionGains=refined.plan;auto second=measure(refined);if(!second.value)return {{},second.error,second.message};}
            }
            if(out->candidate&&bestRank&&(*bestRank)[0]!=0)
            {
                const auto global=lraFinalGlobalGain(*out->candidate,r.currentContext.profile);
                if(global&&*global!=out->candidate->context.processing.staticGainDb)
                {
                    context.processing=out->candidate->context.processing;context.processing.staticGainDb=*global;
                    out->narrative+="Final global trim estimate "+juce::String(*global,12)+" dB. A positive value raises the complete programme after local attenuation; it is not a selective quiet-passage boost.\n";
                    auto final=measure(bestSelection);if(!final.value)return {{},final.error,final.message};
                }
                else if(!global)out->narrative+="Method: no scalar global gain intersects the selected IL/TP/M/S constraints for the retained local plan. No other strategy attempted.\n";
            }
        }
        out->lraSeed=bestSelection?bestSelection:out->lraSeed;
        if(bestRank&&(*bestRank)[0]==0&&out->candidate)
        {
            out->state=SectionProposalStatus::ready;out->problem=SectionProposalIssue::none;out->proposed=out->candidate->context.processing;
            out->explanation="Exact complete measured assessment "+juce::String(bestIndex)+" passes all selected requirements. Global gain "+juce::String(out->proposed.staticGainDb,12)+" dB. Apply explicitly while stopped, then independently remeasure Current. Private measurements do not verify Current or an exported file.";
        }
        else
        {
            out->state=SectionProposalStatus::noSolution;out->problem=SectionProposalIssue::candidateRejected;
            out->explanation="Measured failure: this bounded attempt found no acceptable complete plan. ";
            if(out->candidate)out->explanation+="Retained measured attempt "+juce::String(bestIndex)+" has LRA "+juce::String(*out->candidate->detail.lra.value,6)+" LU against "+juce::String(limits.maxLoudnessRangeLu,6)+" LU. ";
            out->explanation+="Method/search limit: attenuation of selected high regions only, at most three complete assessments; no claim that the programme has no solution. Preferences were not relaxed. See Review for measured rule failures, gain/duration restrictions and estimates.";
        }
        if(bestSelection)out->narrative+="\nRetained selection estimates/provenance\n"+bestSelection->text();
        return complete();
    }
    catch(const std::exception& e){return {{},ErrorCode::invalidRequest,"LRA proposal failed safely: "+juce::String(e.what())};}
}
}

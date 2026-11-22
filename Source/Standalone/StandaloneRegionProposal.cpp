#include "StandaloneRegionProposal.h"
#include "RegionReviewFormatting.h"
#include "StandaloneLraProposal.h"
#include "StandaloneProcessing.h"
#include "AppliedProposalExplanation.h"
#include <set>
#include <numeric>
#include <thread>

namespace levelscope::standalone
{
namespace
{
using namespace proposal;
using Clock=std::chrono::steady_clock;
double seconds(Clock::time_point start){return std::chrono::duration<double>(Clock::now()-start).count();}
bool finite(const measurement::Metric& m){return m.status==measurement::Status::available&&m.value&&std::isfinite(*m.value);}
processing::SectionGainPlan baselinePlan(const ServiceContext& c)
{
    if(c.processing.sectionGains)return c.processing.sectionGains->canonical();
    processing::SectionGainPlan p;p.inputIdentity=sectionInputIdentity(c);p.frames=c.media.frames;p.sampleRate=int(c.media.sampleRate);p.intervals={{0,p.frames,0}};return p;
}
ProcessingPlan processingPlan(const RegionPlan& p)
{
    ProcessingPlan out;out.staticGainDb=p.globalDb;
    if(!p.envelope.isUnity()){out.sectionGains=p.envelope;out.schemaVersion=2;out.engineVersion=processing::SectionGainPlan::version;}return out;
}
struct ProxyChoice {RegionPlan plan;RegionPrediction forecast;RegionReviewScore score;bool measured=false;};
struct Alternative {size_t root=0;RegionAssignment assignment;};
constexpr size_t planningBudget=64*1024*1024;
size_t assignmentBytes(const RegionAssignment& a)
{
    return sizeof(a)+a.geometry.cuts.capacity()*sizeof(int64_t)
        +(a.reductions.capacity()+a.requested.capacity()+a.reserve.capacity()+a.uncovered.capacity())*sizeof(double)
        +a.geometry.reason.getNumBytesAsUTF8()*4;
}
size_t planBytes(const RegionPlan& p)
{
    size_t total=sizeof(p)+p.envelope.intervals.capacity()*sizeof(processing::GainInterval)
        +p.activity.support.capacity()*sizeof(RegionSpan)+p.activity.splits.capacity()*sizeof(size_t);
    for(const auto& a:p.assignments)total+=assignmentBytes(a);return total;
}
std::vector<Alternative> alternatives(const RegionConstraints& c,const RegionFormation& formed,const ProgramAnalysis& a,
    const processing::SectionGainPlan& before,const ProposalPolicy& policy,double ceiling,const std::function<bool()>& cancel)
{
    std::vector<std::vector<RegionAssignment>> perRoot(c.roots.size());
    std::vector<std::vector<std::vector<RegionAssignment>>> geometries(c.roots.size());
    for(size_t root=0;root<c.roots.size();++root)for(const auto& g:formed.alternatives[root])
    {
        if(cancel())return {};
        geometries[root].emplace_back();auto& variants=geometries[root].back();
        std::vector<int64_t> cuts{g.support.begin};cuts.insert(cuts.end(),g.cuts.begin(),g.cuts.end());cuts.push_back(g.support.end);
        std::vector<std::vector<double>> gains;RegionAssignment seed;seed.geometry=g;
        for(size_t part=1;part<cuts.size();++part)
        {
            auto need=regionGainNeed(before,a,{cuts[part-1],cuts[part]},policy,ceiling);
            seed.reductions.push_back(need.choices.size()==1?0:std::min(need.headroom,need.requested));gains.push_back(need.choices);
            seed.requested.push_back(need.choices.size()==1?0:need.requested);
        }
        if(seed.reductions.size()==3&&g.reason.startsWith("Reserved two-cut unchanged-middle"))seed.reductions[1]=0;
        variants.push_back(seed);
        for(int variant=0;variant<2;++variant)
        {
            auto trial=seed;for(size_t part=0;part<gains.size();++part)trial.reductions[part]=variant==0?gains[part][std::min<size_t>(1,gains[part].size()-1)]:gains[part].back();variants.push_back(std::move(trial));
        }
        for(size_t part=0;part<gains.size();++part)for(double gain:gains[part])if(gain!=seed.reductions[part])
        {auto trial=seed;trial.reductions[part]=gain;variants.push_back(std::move(trial));}
        if(c.settings.fillShortGap&&seed.reductions.size()>1)
        {
            auto filled=seed;const auto common=regionGainNeed(before,a,g.support,policy,ceiling);const double gain=std::min(common.headroom,common.requested);bool changed=false;
            const int64_t shortGap=int64_t(c.settings.separationMs)*c.sampleRate/1000+2*std::llround(c.sampleRate*.05);
            for(size_t n=1;n+1<seed.reductions.size();++n)if(seed.reductions[n]==0&&seed.reductions[n-1]>0&&seed.reductions[n+1]>0&&cuts[n+1]-cuts[n]<shortGap)
            {filled.reductions[n]=gain;filled.requested[n]=common.requested;changed=true;}
            if(changed){filled.geometry.reason+="; explicitly permitted separately assessed short-gap fill";variants.push_back(std::move(filled));}
        }
        for(auto& trial:variants)
        {
            const auto radius=int64_t(std::llround(c.sampleRate*.05));
            // Zero exterior parts refine the support, not the original root or its
            // split budget. Pins are still checked against the resulting support.
            if(c.settings.allowRefinement)
            {
                while(trial.reductions.size()>1&&trial.reductions.front()==0)
                {trial.geometry.support.begin=trial.geometry.cuts.front()-radius;trial.geometry.cuts.erase(trial.geometry.cuts.begin());trial.reductions.erase(trial.reductions.begin());trial.requested.erase(trial.requested.begin());}
                while(trial.reductions.size()>1&&trial.reductions.back()==0)
                {trial.geometry.support.end=trial.geometry.cuts.back()+radius;trial.geometry.cuts.pop_back();trial.reductions.pop_back();trial.requested.pop_back();}
            }
            for(size_t n=0;n<trial.reductions.size();++n)
            {
                const double deficit=std::max(0.,trial.requested[n]-policy.reserveTenths/10.);
                trial.reserve.push_back(std::min(policy.reserveTenths/10.,std::max(0.,trial.reductions[n]-deficit)));
                trial.uncovered.push_back(std::max(0.,deficit-trial.reductions[n]));
            }
        }
    }
    for(size_t root=0;root<geometries.size();++root)for(size_t variant=0;;++variant)
    {bool any=false;for(auto& geometry:geometries[root])if(variant<geometry.size()){perRoot[root].push_back(std::move(geometry[variant]));any=true;}if(!any)break;}
    std::vector<Alternative> out;for(size_t n=0;;++n){bool any=false;for(size_t r=0;r<perRoot.size();++r)if(n<perRoot[r].size()){out.push_back({r,std::move(perRoot[r][n])});any=true;}if(!any)break;}return out;
}
}

juce::String RegionProposalReview::text()const
{
    juce::String s="Experimental region planning / "+juce::String(constraints.settings.methodVersion())+" / "+proposal::regionProxyVersion+" / "+(eventMode?proposal::eventRankVersion:proposal::regionRankVersion)+"\n";
    s+="Effective immutable settings: "+juce::JSON::toString(constraints.settings.encode(),true)+"\nConstraint binding "+constraints.identity+"\n";
    s+="Roots "+juce::String(constraints.roots.size())+" (omitted "+juce::String(constraints.omittedRoots)+"); proxy plans "+juce::String(proxyEvaluations)+" / 64; duplicate audio plans "+juce::String(duplicatePlans)+"; rejected "+juce::String(rejectedPlans)+".\n";
    s+="Geometry constructions "+juce::String(formation.constructions)+"; omitted cut nominations "+juce::String(formation.omittedCuts)+"; omitted geometries "+juce::String(formation.omittedGeometries)+". Search truncation never proves impossibility.\n";
    s+=constraints.settings.allowGlobal?"Whole-programme global trim permitted separately: it may change audio outside the local area.\n":"Whole-programme global gain is fixed. Local protections do not claim a general stateful output freeze.\n";
    for(const auto& r:constraints.roots)s+="Original root "+r.id+"; hint ["+juce::String(r.hint.begin)+", "+juce::String(r.hint.end)+"); permitted ["+juce::String(r.allowed.begin)+", "+juce::String(r.allowed.end)+"); support-edge ranges ["+juce::String(r.start.minimum)+", "+juce::String(r.start.maximum)+"] / ["+juce::String(r.end.minimum)+", "+juce::String(r.end.maximum)+"].\n";
    if(eventMode)
    {
        s+="Event result "+resultIdentity+"; feature reserved/high-water bytes "+juce::String(featureBytes)+" / "+juce::String(featureHighWaterBytes)+".\n";
        s+="Recommended attempt "+(recommendedAttempt?juce::String(*recommendedAttempt+1):"none")+"; selected attempt "+(selectedAttempt?juce::String(*selectedAttempt+1):"none")+". Inspection alone never selects Apply.\n";
        for(const auto& root:constraints.roots)s+="Permission: "+root.permissionNotice+"\n";
        if(eventEvidence)s+="Current-only rise cells (processing-associated context, not Source event evidence): "+juce::String(std::count_if(eventEvidence->cells.begin(),eventEvidence->cells.end(),[](const auto& cell){return (cell.flags&2)!=0;}))+"\n";
        if(eventEvidence&&constraints.settings.schema==3)s+=eventEvidence->innerDescription();
        if(eventEvidence)for(const auto& g:eventEvidence->groups)s+="Source event "+g.id+" onset ["+juce::String(g.onset.begin)+", "+juce::String(g.onset.end)+") decay ["+juce::String(g.decay.begin)+", "+juce::String(g.decay.end)+"); onset raw-peak witness channel "+juce::String(g.onsetPeakChannel+1)+" frame "+juce::String(g.onsetPeakFrame)+"; baseline included need-window indices "+juce::String(g.firstNeedWindow)+".."+juce::String(g.lastNeedWindow)+" ("+juce::String(g.needWindowCount)+" included, exact inclusion/ceiling predicate); "+g.reason+"\n";
    }
    for(size_t n=0;n<attempts.size();++n)
    {
        const auto& a=attempts[n];s+="\nComplete measured attempt "+juce::String(n+1)+" / "+a.plan.identity+(n==retainedAttempt?" (retained)":"")+"\n";
        s+="Actual "+evaluationText(a.evaluation)+"\n";
        if(eventMode)s+="Family "+a.family+" / attempt "+a.attemptIdentity+"\n"+a.coverage.text();
        const char* labels[]={"IL","TP","LRA","M max","S max","sample peak"};for(size_t k=0;k<a.measured.size();++k)s+=juce::String(labels[k])+": "+(a.measured[k].value?juce::String(*a.measured[k].value,12):"unavailable/status "+juce::String(int(a.measured[k].status)))+"\n";
        s+="Forecast IL "+(a.prediction.integrated?juce::String(*a.prediction.integrated,12):"unknown")+"; LRA "+(a.prediction.lra?juce::String(*a.prediction.lra,12):"unknown")+"; TP unknown (not a measured failure). "+a.prediction.limitations+"\n";
        s+="Forecast gate entry/exit: I "+juce::String(a.prediction.integratedEntered)+"/"+juce::String(a.prediction.integratedLeft)+"; LRA "+juce::String(a.prediction.lraEntered)+"/"+juce::String(a.prediction.lraLeft)+".\n";
        s+="Actual minus forecast (IL/LRA/P10/P95/sampled M):";for(const auto& error:a.signedPredictionError)s+=" "+(error?juce::String(*error,12):juce::String("unknown"));s+=". Empirical errors, not bounds.\n";
        s+="Changed support "+juce::String(a.plan.activity.affectedFrames)+" programme frames; "+juce::String(a.plan.activity.support.size())+" connected regions; max added "+juce::String(a.plan.activity.maximumAdditional,12)+" dB; absolute local dB*frames "+juce::String(a.plan.activity.absoluteDbFrames,12)+"; global "+juce::String(a.plan.globalDb,12)+" dB.\n";
        s+=regionActivityDescription(a.plan,constraints.sampleRate,constraints.frames);
        for(const auto& part:a.plan.assignments)
        {
            s+="Root "+constraints.roots[part.geometry.root].id+"; "+part.geometry.reason+"; support ["+juce::String(part.geometry.support.begin)+", "+juce::String(part.geometry.support.end)+"); nominal cuts "+juce::String(part.nominal.begin);
            for(auto cut:part.geometry.cuts)s+=", "+juce::String(cut);s+=", "+juce::String(part.nominal.end)+"; additional reductions";for(auto gain:part.reductions)s+=" "+juce::String(gain,12);s+=" dB.\n";
            s+="Requested / achieved reserve / uncovered deficit by part:";
            for(size_t k=0;k<part.requested.size();++k)s+=" ["+juce::String(part.requested[k],12)+" / "+juce::String(part.reserve[k],12)+" / "+juce::String(part.uncovered[k],12)+"]";s+=" dB.\n";
        }
    }
    for(const auto& why:restrictions)s+=why+"\n";
    for(const auto& why:formation.reasons)s+="Formation limitation: "+why+"\n";
    return s+"Private assessment budget 3. Failed shortlist 0.01 LU is diagnostic ranking only, never a rule allowance.\nPreparation "+juce::String(preparationSeconds,3)+" s; auxiliary pass "+juce::String(sidecarSeconds,3)+" s; cached planning "+juce::String(planningSeconds,3)+" s; sidecar "+juce::String(sidecarBytes)+" bytes; bounded additional planning allocation estimate "+juce::String(planningBytes)+" bytes.\n";
}

juce::String validateRegionProposalForApply(const SectionProposalRequest& r,const SectionProposalResult& result)
{
    if(proposal::isEventMethod(result.methodVersion()))return validateEventRegionProposalForApply(r,result);
    if(result.methodVersion()!=proposal::regionSearchVersion)return {};
    if(!r.regionSettings||result.identity()!=r.binding()||!result.regionReview())return "Stale or missing region constraints";
    const auto& review=*result.regionReview();
    if(review.retainedAttempt>=review.attempts.size())return "Missing exact measured region attempt";
    auto plan=review.attempts[review.retainedAttempt].plan;
    if(proposalPlanIdentity(result.plan())!=proposalPlanIdentity(processingPlan(plan)))return "Measured plan differs from recorded region geometry";
    if(review.constraints.settings.digest()!=r.regionSettings->digest())return "Region permissions changed";
    return proposal::validateRegionPlan(baselinePlan(r.currentContext),r.currentContext.processing.staticGainDb,review.constraints,r.policy,plan);
}

SectionProposalOutcome generateRegionProposal(const SectionProposalRequest& r,JobControl& control,const SnapshotStorage& storage)
{
    using namespace proposal;
    auto out=std::shared_ptr<SectionProposalResult>(new SectionProposalResult);out->bound=r.binding();out->method=regionSearchVersion;out->proposed=r.currentContext.processing;out->annotations=r.annotations;out->seed.policy=r.policy;
    auto review=std::make_shared<RegionProposalReview>();review->resultIdentity=juce::Uuid().toString();out->regions=review;
    const auto finish=[&]()->SectionProposalOutcome
    {
        if(control.isCancelled())return {{},ErrorCode::cancelled,"Region assessment cancelled; Current unchanged"};
        if(r.current&&!r.currentContext.sourceFilesMatch())return {{},ErrorCode::sourceChanged,"Programme member changed during region assessment"};
        out->title=out->state==SectionProposalStatus::ready?"Measured region candidate ready":out->state==SectionProposalStatus::noChange?"Already meets requirements":out->state==SectionProposalStatus::noSolution?"No acceptable candidate in bounded region search":"Region planning unavailable";
        out->reason=out->explanation;out->narrative=out->title+"\n"+out->explanation+"\n";
        if(review->constraints.settings.validationError().isEmpty())out->narrative+=review->text();
        out->explanation=out->narrative+"\n"+r.policy.text()+"\n"+proposalTargetText(r.currentContext.profile);
        if(r.current)out->explanation+=proposalMeasurementText(proposalMeasurementRecord(*r.current,r.currentContext.profile,r.applicability),"Frozen measured Current baseline");
        control.setProgress(1);return {std::shared_ptr<const SectionProposalResult>(out),ErrorCode::none,"Region assessment complete; Current unchanged"};
    };
    try
    {
        const auto scalar=r.scalar();out->scalarResult=proposal::draftGlobalTrim(scalar);
        if(!r.source||!r.current||scalar.sourceInputIdentity.isEmpty()||scalar.sourceInputIdentity!=scalar.currentInputIdentity||!scalar.source.currentIdentityMatches||!scalar.current.currentIdentityMatches||!r.sourceContext.sourceFilesMatch()||!r.currentContext.sourceFilesMatch())
        {out->explanation="Valid matching measured Source/Current and confirmed ordered roles/range are required";return finish();}
        const auto baseline=evaluateMeasuredProgramme(r.currentContext.profile,scalar.current);
        if(baseline.verdict==Verdict::pass){out->state=SectionProposalStatus::noChange;out->problem=SectionProposalIssue::none;out->explanation="Every selected measured requirement passes. No plan/revision/measurement change or private analysis.";return finish();}
        if(!r.regionSettings||r.regionSettings->validationError().isNotEmpty()||r.policy.validationError().isNotEmpty())
        {out->problem=SectionProposalIssue::preferences;out->explanation="Invalid explicit region settings or proposal preferences";return finish();}
        const auto& limits=r.currentContext.profile.limits;
        if(baseline.verdict!=Verdict::fail||!limits.hasMaxLoudnessRange||!limits.hasIntegratedLoudnessTarget||!limits.hasTruePeakCeiling||!r.current->detail.lraStable||!finite(r.current->detail.lra)||!finite(r.source->detail.integrated)||!finite(r.source->detail.truePeak))
        {out->problem=SectionProposalIssue::unsupportedRequirements;out->explanation="Region method needs available stable LRA, selected IL/TP/LRA and every additional selected metric. Unavailable automatic correction does not prove no solution.";return finish();}
        if(r.pendingManualDraft){out->explanation="Resolve the pending manual draft before region assessment";return finish();}
        if(r.currentContext.fileLayout.count!=2&&r.currentContext.fileLayout.count!=6){out->explanation="Experimental region planning supports confirmed stereo/conventional 5.1 only";return finish();}
        auto context=r.currentContext;const auto before=baselinePlan(context);
        auto probe=context;probe.processing.sectionGains=before;probe.processing.sectionGains->intervals={{0,before.frames,-1}};probe.processing.schemaVersion=2;probe.processing.engineVersion=processing::SectionGainPlan::version;
        const auto supported=OfflineAnalysisService::validateRequest(probe,MeasurementDomain::privateCandidate);if(supported.isNotEmpty()){out->explanation=supported;return finish();}
        const auto cancelled=[&]{return control.isCancelled();};auto analysis=analyzeProgramme(r.current->detail,r.current->identity,cancelled);if(!analysis)return finish();
        const auto passages=findLraPassages(analysis,limits.maxLoudnessRangeLu,cancelled);if(!passages)return finish();
        review->constraints=resolveRegionConstraints(*r.regionSettings,r.annotations,*passages,before);
        if(review->constraints.error.isNotEmpty()){out->explanation=review->constraints.error;return finish();}
        const auto globalFirst=*r.current->detail.lra.value<=limits.maxLoudnessRangeLu&&r.regionSettings->allowGlobal?lraFinalGlobalGain(*r.current,r.currentContext.profile):std::nullopt;
        if(!globalFirst&&r.regionSettings->mode!=RegionMode::automatic)
        {
            const auto h=review->constraints.roots.front().hint;
            const bool need=analysis->p10&&std::any_of(analysis->lraWindows.begin(),analysis->lraWindows.end(),[&](const auto& w)
            {const auto centre=w.begin+(w.end-w.begin)/2;return centre>=h.begin&&centre<h.end&&w.included&&w.lufs&&*w.lufs>*analysis->p10+limits.maxLoudnessRangeLu;});
            if(!need){out->state=SectionProposalStatus::noSolution;out->explanation="No included high correction need is identified inside the original selection. Margins do not authorize correction of unrelated neighbouring material.";return finish();}
        }
        if(review->constraints.roots.empty()&&!globalFirst){out->state=SectionProposalStatus::noSolution;out->explanation="No high analysis seeds in this bounded automatic method; no candidate measured. This is not an impossibility claim.";return finish();}
        auto start=Clock::now();JobControl preparing;preparing.parent=&control;preparing.parentProgress=&control.progress;preparing.progressScale=.1f;
        const size_t binCount=size_t((before.frames+before.sampleRate/10-1)/(before.sampleRate/10)+17);
        // Reserve concurrent model scratch, rank vectors, bounded strings/metadata
        // and allocation capacity before allocating bins or evaluating a proxy.
        const size_t modelScratch=binCount*2*sizeof(double)
            +(analysis->momentary.size()+analysis->shortTerm.size()+analysis->lraWindows.size())*sizeof(measurement::Window)
            +analysis->lraWindows.size()*2*sizeof(double)+1024*1024;
        const size_t fixedPlanning=passages->storageBytes()+binCount*sizeof(PredictionBin)+modelScratch;
        if(fixedPlanning>planningBudget)return {{},ErrorCode::unavailable,"Region model cannot fit the 64 MiB planning budget; no candidate authorized"};
        const auto copy=prepareVerifiedSource(context,preparing,storage);review->preparationSeconds=seconds(start);if(!copy.value)return {{},copy.error,copy.message};
        std::shared_ptr<const PredictionSidecar> sidecar;
        if(!globalFirst)
        {
            start=Clock::now();auto reader=(**copy.value).createReader();if(!reader)return {{},ErrorCode::decodeFailed,"Cannot open verified programme for auxiliary analysis"};
            PredictionCollector collector(before.sampleRate,OfflineAnalysisService::measurementRoles(context.fileLayout),before.frames,out->bound,r.current->identity,regionPlanIdentity(before,r.currentContext.processing.staticGainDb));
            ProgrammeGainProcessor processor(context.processing);juce::AudioBuffer<float> audio(context.media.channels,4096);
            for(int64_t frame=0;frame<before.frames;)
            {
                if(control.isCancelled())return finish();const int count=int(std::min<int64_t>(4096,before.frames-frame));
                if(!reader->read(audio.getArrayOfWritePointers(),audio.getNumChannels(),frame,count))return {{},ErrorCode::decodeFailed,"Auxiliary analysis decode failed"};
                if(!processor.process(audio.getArrayOfWritePointers(),audio.getNumChannels(),count,frame))return {{},ErrorCode::nonFiniteAudio,"Auxiliary analysis non-finite audio"};
                collector.push(audio.getArrayOfReadPointers(),count);frame+=count;control.setProgress(.1f+.15f*float(double(frame)/before.frames));
            }
            sidecar=collector.finish();review->sidecarSeconds=seconds(start);review->sidecarBytes=sidecar->storageBytes();
        }
        std::vector<ProxyChoice> pool;std::set<juce::String> seen;std::vector<Alternative> choices;size_t cursor=0;
        double planningSeconds=0;bool exhausted=false;const double watchdog=
#if defined(__x86_64__)
            60;
#else
            30;
#endif
        Clock::time_point planningStart;bool planning=false;
        const auto abortPlanning=[&]{if(planning&&planningSeconds+seconds(planningStart)>watchdog)exhausted=true;return control.isCancelled()||exhausted;};
        const auto beginPlanning=[&]{planningStart=Clock::now();planning=true;};
        const auto endPlanning=[&]{planningSeconds+=seconds(planningStart);planning=false;review->planningSeconds=planningSeconds;};
        beginPlanning();
        if(!globalFirst){review->formation=formRegions(review->constraints,*analysis,limits.maxLoudnessRangeLu,abortPlanning);choices=alternatives(review->constraints,review->formation,*analysis,before,r.policy,limits.maxLoudnessRangeLu,abortPlanning);}endPlanning();
        size_t geometryBytes=0;for(const auto& choice:choices)geometryBytes+=assignmentBytes(choice.assignment);
        for(const auto& root:review->formation.alternatives)for(const auto& geometry:root)geometryBytes+=sizeof(geometry)+geometry.cuts.capacity()*sizeof(int64_t)+geometry.reason.getNumBytesAsUTF8()*4;
        review->planningBytes=fixedPlanning+geometryBytes;
        if(review->planningBytes>planningBudget)return {{},ErrorCode::unavailable,"Region geometry exceeds the 64 MiB planning budget; no candidate authorized"};
        auto add=[&](RegionPlan plan)->std::optional<size_t>
        {
            if(abortPlanning()||review->proxyEvaluations>=64)return {};
            if(plan.error.isNotEmpty())
            {++review->rejectedPlans;if(review->restrictions.size()<48&&std::find(review->restrictions.begin(),review->restrictions.end(),plan.error)==review->restrictions.end())review->restrictions.push_back(plan.error);return {};}
            if(!seen.insert(plan.identity).second){++review->duplicatePlans;return {};}
            if(plan.envelope.audioKey()==before.audioKey()&&plan.globalDb==r.currentContext.processing.staticGainDb)return {};
            size_t used=fixedPlanning+geometryBytes+2*planBytes(plan);
            for(const auto& p:pool)used+=sizeof(ProxyChoice)+planBytes(p.plan)+p.forecast.limitations.getNumBytesAsUTF8()*4;
            for(const auto& attempt:review->attempts)used+=sizeof(RegionAttempt)+planBytes(attempt.plan)+(attempt.receipt?attempt.receipt->storageBytes():0);
            if(used>planningBudget)throw std::runtime_error("Region model/search 64 MiB allocation budget exceeded; no partial Apply authority");
            review->planningBytes=std::max(review->planningBytes,used);
            ProxyChoice choice;choice.plan=std::move(plan);
            if(sidecar)choice.forecast=predictRegionPlan(*analysis,*sidecar,before,r.currentContext.processing.staticGainDb,choice.plan,abortPlanning);
            else
            {const double delta=choice.plan.globalDb-r.currentContext.processing.staticGainDb;choice.forecast.integrated=*r.current->detail.integrated.value+delta;choice.forecast.lra=r.current->detail.lra.value;choice.forecast.limitations="Global-only estimate; full candidate measurement remains required";}
            if(abortPlanning()||choice.forecast.cancelled)return {};
            choice.score=regionPredictedScore(choice.forecast,r.currentContext.profile,choice.plan.activity,choice.plan.identity);++review->proxyEvaluations;pool.push_back(std::move(choice));return pool.size()-1;
        };
        const auto bestProxy=[&]()->std::optional<size_t>
        {
            std::vector<RegionReviewScore> ranks;for(const auto& p:pool){auto s=p.score;if(p.measured)s.complete=false;ranks.push_back(s);}
            const auto best=selectRegionReview(ranks);if(best)return best;
            // An unavailable forecast is not a measured rejection. Keep one
            // deterministic eligible challenger if all forecasts are unavailable.
            for(size_t n=0;n<pool.size();++n)if(!pool[n].measured)return n;return {};
        };
        const auto beam=[&]
        {
            std::vector<size_t> result;std::vector<RegionReviewScore> ranks;for(const auto& p:pool)ranks.push_back(p.score);
            for(size_t n=0;n<4;++n){const auto selected=selectRegionReview(ranks);if(!selected)break;result.push_back(*selected);ranks[*selected].complete=false;}
            return result;
        };
        const auto globalVariants=[&](size_t index,size_t stageEnd)
        {
            if(!r.regionSettings->allowGlobal||!pool[index].forecast.integrated)return;
            const double current=pool[index].plan.globalDb,il=*pool[index].forecast.integrated;
            for(double desired:{limits.integratedLufs,limits.integratedLufs-limits.integratedToleranceLu,limits.integratedLufs+limits.integratedToleranceLu})
            {
                if(review->proxyEvaluations>=stageEnd||abortPlanning())break;auto plan=pool[index].plan;plan.globalDb=std::clamp(current+desired-il,-60.,24.);plan.error=validateRegionPlan(before,r.currentContext.processing.staticGainDb,review->constraints,r.policy,plan);add(std::move(plan));
            }
        };
        if(globalFirst)
        {
            auto plan=constructRegionPlan(before,r.currentContext.processing.staticGainDb,review->constraints,r.policy,{},*globalFirst);beginPlanning();add(std::move(plan));endPlanning();
        }
        else
        {
            // Preserve the earlier method's full-resolution gain recipe as a
            // complete unsplit challenger, with the NEW permissions enforced.
            // Its nominal M19 outer cuts may need to move inward for ramps to fit.
            beginPlanning();std::vector<RegionAssignment> seeds;
            for(size_t root=0;root<review->constraints.roots.size();++root)
            {
                auto found=std::find_if(choices.begin(),choices.end(),[&](const auto& a){return a.root==root&&a.assignment.geometry.cuts.empty()&&a.assignment.geometry.support==review->constraints.roots[root].hint;});
                if(found==choices.end())continue;
                auto candidate=seeds;candidate.push_back(found->assignment);
                auto plan=constructRegionPlan(before,r.currentContext.processing.staticGainDb,review->constraints,r.policy,candidate,r.currentContext.processing.staticGainDb);
                if(plan.error.isEmpty())seeds=std::move(candidate);
            }
            if(!seeds.empty()){auto index=add(constructRegionPlan(before,r.currentContext.processing.staticGainDb,review->constraints,r.policy,seeds,r.currentContext.processing.staticGainDb));if(index)globalVariants(*index,32);}
            endPlanning();
        }
        OfflineAnalysisService worker;std::vector<RegionReviewScore> measuredRanks;

        for(int pass=0;pass<3;++pass)
        {
            beginPlanning();const size_t stageEnd=pass==0?32:pass==1?48:64;
            if(!globalFirst)
            {
                // Feedback can nominate an exact complete global alternative, but it
                // consumes another proxy and another private assessment if selected.
                if(pass&&out->eventReceipt&&r.regionSettings->allowGlobal)
                {
                    const auto gain=lraFinalGlobalGain(out->eventReceipt->measured(),out->eventReceipt->context().processing.staticGainDb,r.currentContext.profile);
                    if(gain){auto plan=review->attempts[review->retainedAttempt].plan;plan.globalDb=*gain;plan.error=validateRegionPlan(before,r.currentContext.processing.staticGainDb,review->constraints,r.policy,plan);add(std::move(plan));}
                }
                size_t rejectedBudget=0;
                while(review->proxyEvaluations<stageEnd&&cursor<choices.size()&&!abortPlanning()&&rejectedBudget++<1024)
                {
                    const auto& alt=choices[cursor++];std::vector<RegionAssignment> selected{alt.assignment};
                    // At most four ranked complete plans form the beam. Root visits
                    // are round-robin; replacing one root never resets its split cap.
                    const auto leaders=beam();if(!leaders.empty()&&choices.size()>1)
                    {
                        const auto parent=leaders[(cursor-1)%leaders.size()];
                        for(const auto& a:pool[parent].plan.assignments)if(a.geometry.root!=alt.root)selected.push_back(a);
                    }
                    auto p=constructRegionPlan(before,r.currentContext.processing.staticGainDb,review->constraints,r.policy,std::move(selected),r.currentContext.processing.staticGainDb);
                    const auto index=add(std::move(p));if(index)globalVariants(*index,stageEnd);
                }
            }
            auto next=bestProxy();
            // If two refined attempts failed, reserve the final assessment for a
            // legal unmeasured unsplit challenger rather than deepening forever.
            if(pass==2&&std::none_of(measuredRanks.begin(),measuredRanks.end(),[](auto s){return s.allPass;}))
            {
                std::vector<RegionReviewScore> ranks;for(const auto& p:pool){auto rank=p.score;if(p.measured||std::any_of(p.plan.activity.splits.begin(),p.plan.activity.splits.end(),[](auto n){return n!=0;}))rank.complete=false;ranks.push_back(rank);}
                const auto challenger=selectRegionReview(ranks);if(challenger)next=challenger;
            }
            endPlanning();
            if(exhausted)return {{},ErrorCode::unavailable,"Region planning watchdog exhausted; no partial Apply authority"};
            if(control.isCancelled())return finish();if(!next)break;
            auto chosen=pool[*next];pool[*next].measured=true;
            context.processing=processingPlan(chosen.plan);JobControl assessment;assessment.parent=&control;assessment.parentProgress=&control.progress;assessment.progressBase=.25f+float(pass)*.24f;assessment.progressScale=.24f;
            start=Clock::now();++out->passCount;const auto actual=worker.analyzeCandidate(context,**copy.value,assessment);if(!actual.value)return {{},actual.error,actual.message};
            RegionAttempt attempt;attempt.plan=chosen.plan;attempt.prediction=chosen.forecast;attempt.seconds=seconds(start);attempt.snapshotIdentity=actual.value->identity;
            attempt.evaluation=evaluateSnapshot(&*actual.value,context,MeasurementDomain::privateCandidate,r.currentContext.profile,r.applicability);
            const auto& d=actual.value->detail;attempt.measured={d.integrated,d.truePeak,d.lra,d.momentaryMax,d.shortTermMax,d.samplePeak};
            std::array<std::optional<double>,5> predicted{chosen.forecast.integrated,chosen.forecast.lra,chosen.forecast.p10,chosen.forecast.p95,chosen.forecast.momentarySampled};
            std::array<std::optional<double>,5> values{d.integrated.value,d.lra.value,d.p10,d.p95,d.momentaryMax.value};
            for(size_t n=0;n<values.size();++n)if(values[n]&&predicted[n])attempt.signedPredictionError[n]=*values[n]-*predicted[n];
            attempt.receipt=std::shared_ptr<const RegionAssessmentReceipt>(new RegionAssessmentReceipt(*actual.value,r));
            const auto receiptKey=review->resultIdentity+";"+attempt.plan.identity+";"+attempt.receipt->identity();
            attempt.attemptIdentity=juce::SHA256(receiptKey.toRawUTF8(),receiptKey.getNumBytesAsUTF8()).toHexString();
            review->planningBytes+=sizeof(RegionAttempt)+planBytes(attempt.plan)+attempt.receipt->storageBytes();
            if(review->planningBytes>planningBudget)return {{},ErrorCode::unavailable,"Region receipts exceed the 64 MiB planning budget; no candidate authorized"};
            measuredRanks.push_back(regionMeasuredScore(attempt.evaluation,attempt.plan.activity,attempt.plan.identity));review->attempts.push_back(std::move(attempt));
            const auto retained=selectRegionReview(measuredRanks);
            if(retained)
            {
                if(*retained==review->attempts.size()-1)out->candidate=std::make_shared<const AnalysisSnapshot>(*actual.value);
                else if(*retained!=review->retainedAttempt)out->candidate.reset();
                review->retainedAttempt=*retained;out->eventReceipt=review->attempts[*retained].receipt;
            }
            if(retained&&measuredRanks[*retained].allPass&&!r.regionSettings->compareAlternatives)break;
            if(globalFirst)break;
        }
        review->searchTruncated=cursor<choices.size()||review->constraints.omittedRoots||review->formation.omittedCuts||review->formation.omittedGeometries;
        const auto best=selectRegionReview(measuredRanks);
        if(best&&measuredRanks[*best].allPass&&out->eventReceipt)
        {out->state=SectionProposalStatus::ready;out->problem=SectionProposalIssue::none;out->proposed=out->eventReceipt->context().processing;out->explanation="Exact complete candidate passes every selected measured requirement. Apply explicitly while stopped, then independently remeasure Current. Predictions do not verify Current or an export.";}
        else
        {out->state=SectionProposalStatus::noSolution;out->problem=SectionProposalIssue::candidateRejected;out->explanation="No acceptable complete candidate in this finite region search. All measured failures remain inspectable; Apply unavailable. Limits/uncertain forecasts/omitted alternatives do not prove no solution exists.";}
        return finish();
    }
    catch(const std::exception& e){return {{},ErrorCode::invalidRequest,"Region assessment failed safely: "+juce::String(e.what())};}
}
}

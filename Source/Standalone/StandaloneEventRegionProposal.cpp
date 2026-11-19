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
struct ProxyChoice {RegionPlan plan;RegionPrediction forecast;RegionReviewScore score;EventCoverage coverage;bool measured=false,original=false;};
struct Alternative {size_t root=0;RegionAssignment assignment;};
juce::String digest(const juce::String& s){return juce::SHA256(s.toRawUTF8(),s.getNumBytesAsUTF8()).toHexString();}
size_t coverageBytes(const EventCoverage& c)
{
    size_t n=sizeof(c)+c.groups.capacity()*sizeof(GroupTreatment)+c.signature.capacity()*sizeof(int64_t);
    for(const auto& g:c.groups)n+=(g.groupId.getNumBytesAsUTF8()+g.reason.getNumBytesAsUTF8())*4;
    return n+(c.identity.getNumBytesAsUTF8()+c.evidenceIdentity.getNumBytesAsUTF8()+c.notice.getNumBytesAsUTF8())*4;
}
juce::String authority(const SectionProposalRequest& r,const RegionProposalReview& review,const RegionAttempt& a,bool explicitSelection)
{return digest("region-selection-v1;"+r.binding()+";"+review.resultIdentity+";"+a.attemptIdentity+";"+a.receipt->identity()+";"+a.coverage.identity+";"+juce::String(int(explicitSelection)));}
constexpr size_t planningBudget=64*1024*1024;
#if defined(__x86_64__)
constexpr double planningWatchdogSeconds=60;
#else
constexpr double planningWatchdogSeconds=30;
#endif
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

SectionProposalOutcome generateEventRegionProposal(const SectionProposalRequest& r,JobControl& control,const SnapshotStorage& storage)
{
    using namespace proposal;
    auto out=std::shared_ptr<SectionProposalResult>(new SectionProposalResult);out->bound=r.binding();out->method=eventSearchVersion;out->proposed=r.currentContext.processing;out->annotations=r.annotations;out->seed.policy=r.policy;
    auto review=std::make_shared<RegionProposalReview>();out->regions=review;review->eventMode=true;review->resultIdentity=juce::Uuid().toString();
    const auto finish=[&]()->SectionProposalOutcome
    {
        if(control.isCancelled())return {{},ErrorCode::cancelled,"Region assessment cancelled; Current unchanged"};
        if(r.current&&!r.currentContext.sourceFilesMatch())return {{},ErrorCode::sourceChanged,"Programme member changed during region assessment"};
        out->title=out->state==SectionProposalStatus::ready?"Measured region candidate ready":out->state==SectionProposalStatus::reviewRequired?"Review required - passing alternatives, complete treatment not established":out->state==SectionProposalStatus::noChange?"Already meets requirements":out->state==SectionProposalStatus::noSolution?"No acceptable candidate in bounded region search":"Region planning unavailable";
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
        size_t fixedPlanning=passages->storageBytes()+binCount*sizeof(PredictionBin)+modelScratch;
        if(fixedPlanning>planningBudget)return {{},ErrorCode::unavailable,"Region model cannot fit the 64 MiB planning budget; no candidate authorized"};
        const auto copy=prepareVerifiedSource(context,preparing,storage);review->preparationSeconds=seconds(start);if(!copy.value)return {{},copy.error,copy.message};
        std::shared_ptr<const PredictionSidecar> sidecar;
        std::shared_ptr<const EventEvidence> events;double coarsePlanningSeconds=0;
        if(!globalFirst)
        {
            start=Clock::now();auto reader=(**copy.value).createReader();if(!reader)return {{},ErrorCode::decodeFailed,"Cannot open verified programme for auxiliary analysis"};
            PredictionCollector collector(before.sampleRate,OfflineAnalysisService::measurementRoles(context.fileLayout),before.frames,out->bound,r.current->identity,regionPlanIdentity(before,r.currentContext.processing.staticGainDb));
            const auto coarseStart=Clock::now();bool coarseExpired=false;
            const auto coarse=formRegions(review->constraints,*analysis,limits.maxLoudnessRangeLu,[&]{coarseExpired=seconds(coarseStart)>planningWatchdogSeconds;return cancelled()||coarseExpired;});
            coarsePlanningSeconds=seconds(coarseStart);
            if(cancelled())return finish();if(coarseExpired)return {{},ErrorCode::unavailable,"Coarse event nominations exceeded planning watchdog; no Apply authority"};
            start=Clock::now();
            EventFeatureCollector features(review->constraints,coarse,*analysis,OfflineAnalysisService::measurementRoles(context.fileLayout),limits.maxLoudnessRangeLu,out->bound,r.source->identity,regionPlanIdentity(before,context.processing.staticGainDb));
            if(fixedPlanning+features.reservedBytes()+size_t(context.media.channels)*4096*2*sizeof(float)>planningBudget)return {{},ErrorCode::unavailable,"Event features exceed combined planning budget"};
            ProgrammeGainProcessor processor(context.processing);juce::AudioBuffer<float> audio(context.media.channels,4096),dry(context.media.channels,4096);
            for(int64_t frame=0;frame<before.frames;)
            {
                if(control.isCancelled())return finish();const int count=int(std::min<int64_t>(4096,before.frames-frame));
                if(!reader->read(audio.getArrayOfWritePointers(),audio.getNumChannels(),frame,count))return {{},ErrorCode::decodeFailed,"Auxiliary analysis decode failed"};
                for(int channel=0;channel<audio.getNumChannels();++channel)dry.copyFrom(channel,0,audio,channel,0,count);
                if(!processor.process(audio.getArrayOfWritePointers(),audio.getNumChannels(),count,frame))return {{},ErrorCode::nonFiniteAudio,"Auxiliary analysis non-finite audio"};
                features.push(dry.getArrayOfReadPointers(),audio.getArrayOfReadPointers(),count,cancelled);
                collector.push(audio.getArrayOfReadPointers(),count);frame+=count;control.setProgress(.1f+.15f*float(double(frame)/before.frames));
            }
            events=features.finish();review->eventEvidence=events;review->featureBytes=events->storageBytes();review->featureHighWaterBytes=events->highWaterBytes;
            fixedPlanning+=events->highWaterBytes;
            sidecar=collector.finish();review->sidecarSeconds=seconds(start);review->sidecarBytes=sidecar->storageBytes();
        }
        std::vector<ProxyChoice> pool;pool.reserve(64);review->attempts.reserve(3);std::set<juce::String> seen;std::vector<Alternative> choices;size_t cursor=0;
        double planningSeconds=coarsePlanningSeconds;bool exhausted=false;const double watchdog=planningWatchdogSeconds;
        Clock::time_point planningStart;bool planning=false;
        const auto abortPlanning=[&]{if(planning&&planningSeconds+seconds(planningStart)>watchdog)exhausted=true;return control.isCancelled()||exhausted;};
        const auto beginPlanning=[&]{planningStart=Clock::now();planning=true;};
        const auto endPlanning=[&]{planningSeconds+=seconds(planningStart);planning=false;review->planningSeconds=planningSeconds;};
        beginPlanning();
        if(!globalFirst){review->formation=formEventRegions(review->constraints,*analysis,limits.maxLoudnessRangeLu,*events,abortPlanning);choices=alternatives(review->constraints,review->formation,*analysis,before,r.policy,limits.maxLoudnessRangeLu,abortPlanning);}endPlanning();
        size_t geometryBytes=choices.capacity()*sizeof(Alternative)+review->formation.alternatives.capacity()*sizeof(std::vector<RegionGeometry>);for(const auto& choice:choices)geometryBytes+=assignmentBytes(choice.assignment);
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
            // Reserve the next bounded coverage/receipt and concurrent schedule copies
            // before construction. Strings are bounded by the fixed reason vocabulary.
            size_t used=fixedPlanning+geometryBytes+2*planBytes(plan)+pool.capacity()*sizeof(ProxyChoice)
                +(events?events->groups.size()*(sizeof(GroupTreatment)+4096):0)+65536;
            for(const auto& p:pool)used+=sizeof(ProxyChoice)+planBytes(p.plan)+coverageBytes(p.coverage)+p.forecast.limitations.getNumBytesAsUTF8()*4;
            for(const auto& attempt:review->attempts)used+=sizeof(RegionAttempt)+planBytes(attempt.plan)+coverageBytes(attempt.coverage)+(attempt.receipt?attempt.receipt->storageBytes():0);
            if(used>planningBudget)throw std::runtime_error("Region model/search 64 MiB allocation budget exceeded; no partial Apply authority");
            review->planningBytes=std::max(review->planningBytes,used);
            ProxyChoice choice;choice.plan=std::move(plan);
            if(events)choice.coverage=assessEventCoverage(before,choice.plan,review->constraints,*events,*analysis,limits.maxLoudnessRangeLu);
            choice.original=std::all_of(choice.plan.assignments.begin(),choice.plan.assignments.end(),[&](const auto& a){return a.geometry.cuts.empty()&&a.geometry.support==review->constraints.roots[a.geometry.root].hint;});
            if(sidecar)choice.forecast=predictRegionPlan(*analysis,*sidecar,before,r.currentContext.processing.staticGainDb,choice.plan,abortPlanning);
            else
            {const double delta=choice.plan.globalDb-r.currentContext.processing.staticGainDb;choice.forecast.integrated=*r.current->detail.integrated.value+delta;choice.forecast.lra=r.current->detail.lra.value;choice.forecast.limitations="Global-only estimate; full candidate measurement remains required";}
            if(abortPlanning()||choice.forecast.cancelled)return {};
            choice.score=regionPredictedScore(choice.forecast,r.currentContext.profile,choice.plan.activity,choice.plan.identity);++review->proxyEvaluations;pool.push_back(std::move(choice));return pool.size()-1;
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
            // Complete family shares the same finite proxy budget. Preserve its
            // per-root nominees before the generic beam can favour short supports.
            seeds.clear();
            for(size_t root=0;root<review->constraints.roots.size();++root)
            {
                auto found=std::find_if(choices.begin(),choices.end(),[&](const auto& a){return a.root==root&&a.assignment.geometry.reason.startsWith("Supported complete-event");});
                if(found==choices.end())continue;
                auto candidate=seeds;candidate.push_back(found->assignment);
                auto plan=constructRegionPlan(before,context.processing.staticGainDb,review->constraints,r.policy,candidate,context.processing.staticGainDb);
                if(plan.error.isEmpty())seeds=std::move(candidate);
                else if(review->restrictions.size()<48)review->restrictions.push_back("Complete-event alternative: "+plan.error);
            }
            if(!seeds.empty()){auto index=add(constructRegionPlan(before,context.processing.staticGainDb,review->constraints,r.policy,seeds,context.processing.staticGainDb));if(index)globalVariants(*index,32);}
            endPlanning();
        }
        OfflineAnalysisService worker;std::vector<RegionReviewScore> measuredRanks;
        std::vector<EventCoverage> measuredCoverage;measuredCoverage.reserve(3);bool priorityRepairUsed=false;
        for(int pass=0;pass<3;++pass)
        {
            beginPlanning();const size_t stageEnd=pass==0?32:pass==1?48:64;
            std::optional<size_t> priority;
            if(pass&&!priorityRepairUsed&&r.regionSettings->allowGlobal)
            {
                const auto& last=review->attempts.back();
                if(last.coverage.complete()&&last.evaluation.verdict==Verdict::fail)
                {
                    const auto gain=lraFinalGlobalGain(last.receipt->measured(),last.plan.globalDb,r.currentContext.profile);
                    if(gain){auto plan=last.plan;plan.globalDb=*gain;plan.error=validateRegionPlan(before,context.processing.staticGainDb,review->constraints,r.policy,plan);priority=add(std::move(plan));if(priority){priorityRepairUsed=true;review->restrictions.push_back("One measured-feedback scalar repair has priority; it consumes a comparison slot and preserves local geometry.");}}
                }
            }
            if(!globalFirst)
            {
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
            std::vector<EventAssessmentChoice> schedule;
            for(const auto& p:pool)
            {
                const bool screening=p.score.complete&&p.score.failedRules==0
                    &&(!limits.hasMaxMomentary||p.forecast.momentarySampled.has_value())
                    &&(!limits.hasMaxShortTerm||p.forecast.shortTermSampled.has_value());
                schedule.push_back({p.score,p.coverage,p.measured,p.original,screening});
            }
            auto next=priority?priority:chooseEventAssessment(schedule,pass==0,selectCompleteEventPass(measuredRanks,measuredCoverage).has_value());
            endPlanning();
            if(exhausted)return {{},ErrorCode::unavailable,"Region planning watchdog exhausted; no partial Apply authority"};
            if(control.isCancelled())return finish();if(!next)break;
            auto chosen=pool[*next];pool[*next].measured=true;
            context.processing=processingPlan(chosen.plan);JobControl assessment;assessment.parent=&control;assessment.parentProgress=&control.progress;assessment.progressBase=.25f+float(pass)*.24f;assessment.progressScale=.24f;
            start=Clock::now();++out->passCount;const auto actual=worker.analyzeCandidate(context,**copy.value,assessment);if(!actual.value)return {{},actual.error,actual.message};
            RegionAttempt attempt;attempt.plan=chosen.plan;attempt.prediction=chosen.forecast;attempt.seconds=seconds(start);attempt.snapshotIdentity=actual.value->identity;
            attempt.receipt=std::shared_ptr<const RegionAssessmentReceipt>(new RegionAssessmentReceipt(*actual.value,r));
            attempt.coverage=chosen.coverage;attempt.family=globalFirst?"measured global-only":chosen.coverage.complete()?"supported-complete":chosen.original?"original boundaries":"tighter/spared/other";
            attempt.attemptIdentity=digest(review->resultIdentity+";"+chosen.plan.identity+";"+attempt.receipt->identity());
            attempt.evaluation=attempt.receipt->evaluation(r,context.processing);
            const auto& d=actual.value->detail;attempt.measured={d.integrated,d.truePeak,d.lra,d.momentaryMax,d.shortTermMax,d.samplePeak};
            std::array<std::optional<double>,5> predicted{chosen.forecast.integrated,chosen.forecast.lra,chosen.forecast.p10,chosen.forecast.p95,chosen.forecast.momentarySampled};
            std::array<std::optional<double>,5> values{d.integrated.value,d.lra.value,d.p10,d.p95,d.momentaryMax.value};
            for(size_t n=0;n<values.size();++n)if(values[n]&&predicted[n])attempt.signedPredictionError[n]=*values[n]-*predicted[n];
            measuredRanks.push_back(regionMeasuredScore(attempt.evaluation,attempt.plan.activity,attempt.plan.identity));measuredCoverage.push_back(attempt.coverage);review->attempts.push_back(std::move(attempt));
            // actual's full populations die at this iteration boundary; only the
            // privately minted compact receipt remains. No private cache is Current.
            if(globalFirst)break;
        }
        if(review->attempts.size()==3)review->restrictions.push_back("Three exact assessments used; remaining gain/boundary alternatives were not measured. No optimality claim.");
        review->searchTruncated=cursor<choices.size()||review->constraints.omittedRoots||review->formation.omittedCuts||review->formation.omittedGeometries;
        auto best=selectCompleteEventPass(measuredRanks,measuredCoverage);
        if(globalFirst&&!measuredRanks.empty()&&measuredRanks.front().allPass)best=0;
        if(best)
        {
            review->recommendedAttempt=review->selectedAttempt=best;review->retainedAttempt=*best;
            const auto& a=review->attempts[*best];out->eventReceipt=a.receipt;out->proposed=processingPlan(a.plan);
            out->selectionAuthority=authority(r,*review,a,false);out->state=SectionProposalStatus::ready;out->problem=SectionProposalIssue::none;
            out->explanation=globalFirst?"Measured global-only candidate passes; no local event treatment required.":"Supported complete-event candidate passes every selected measured requirement.";
            out->explanation+=" Apply explicitly while stopped, then independently remeasure Current. No artistic-quality or export verification claim.";
        }
        else
        {
            const auto diagnostic=selectRegionReview(measuredRanks);if(diagnostic)review->retainedAttempt=*diagnostic;
            const bool passing=std::any_of(measuredRanks.begin(),measuredRanks.end(),[](const auto& s){return s.allPass;});
            out->state=passing?SectionProposalStatus::reviewRequired:SectionProposalStatus::noSolution;
            out->problem=passing?SectionProposalIssue::none:SectionProposalIssue::candidateRejected;
            out->explanation=passing?"Measured passing alternatives exist, but complete event treatment is not established. Inspection does not select Apply. Explicitly select a measured alternative with its partial/uncertain limitations.":"No acceptable candidate in this bounded search. Measured failures and limits remain visible; this does not establish that no solution exists.";
        }
        return finish();
    }
    catch(const std::exception& e){return {{},control.isCancelled()?ErrorCode::cancelled:ErrorCode::invalidRequest,"Event assessment failed safely: "+juce::String(e.what())};}
}

juce::String validateEventRegionProposalForApply(const SectionProposalRequest& r,const SectionProposalResult& result)
{
    if(result.methodVersion()!=proposal::eventSearchVersion||!r.regionSettings||r.regionSettings->schema!=2||r.pendingManualDraft||result.identity()!=r.binding()||!result.regionReview())return "Stale/missing event request";
    const auto& review=*result.regionReview();
    if(!review.eventMode||!review.selectedAttempt||*review.selectedAttempt>=review.attempts.size()||review.constraints.settings.digest()!=r.regionSettings->digest())return "No explicitly authorized measured event attempt";
    const auto& a=review.attempts[*review.selectedAttempt];
    if(!a.receipt||a.receipt!=result.eventAssessment()||result.selectedAuthority()!=authority(r,review,a,result.explicitlySelectedAlternative())||proposalPlanIdentity(result.plan())!=proposalPlanIdentity(processingPlan(a.plan)))return "Selected result/plan/receipt identity mismatch";
    const bool globalOnly=a.plan.envelope.audioKey()==baselinePlan(r.currentContext).audioKey();
    if(!result.explicitlySelectedAlternative()&&!globalOnly&&(!review.recommendedAttempt||*review.recommendedAttempt!=*review.selectedAttempt||!a.coverage.complete()))return "Partial/uncertain treatment requires explicit alternative selection";
    auto plan=a.plan;const auto error=proposal::validateRegionPlan(baselinePlan(r.currentContext),r.currentContext.processing.staticGainDb,review.constraints,r.policy,plan);if(error.isNotEmpty())return error;
    const auto receiptError=a.receipt->validationError(r,result.plan());if(receiptError.isNotEmpty())return receiptError;
    if(a.receipt->evaluation(r,result.plan()).verdict!=proposal::Verdict::pass)return "Selected measured candidate does not pass current bound requirements";
    return {};
}
SectionProposalOutcome selectMeasuredEventAlternative(const SectionProposalRequest& r,const std::shared_ptr<const SectionProposalResult>& source,const juce::String& resultId,const juce::String& attemptId)
{
    if(!source||source->methodVersion()!=proposal::eventSearchVersion||source->identity()!=r.binding()||r.pendingManualDraft||!source->regionReview()||source->regionReview()->resultIdentity!=resultId)return {{},ErrorCode::staleRevision,"Stale/missing event result; inspection cannot authorize Apply"};
    const auto& old=*source->regionReview();const auto i=std::find_if(old.attempts.begin(),old.attempts.end(),[&](const auto& a){return a.attemptIdentity==attemptId;});
    if(i==old.attempts.end()||!i->receipt||i->receipt->evaluation(r,processingPlan(i->plan)).verdict!=proposal::Verdict::pass)return {{},ErrorCode::invalidRequest,"Only an owned, actual passing measured alternative may be selected"};
    auto result=std::shared_ptr<SectionProposalResult>(new SectionProposalResult(*source));auto review=std::make_shared<RegionProposalReview>(old);
    review->selectedAttempt=size_t(i-old.attempts.begin());review->retainedAttempt=*review->selectedAttempt;result->regions=review;
    result->proposed=processingPlan(i->plan);result->eventReceipt=i->receipt;result->userSelectedAlternative=true;
    result->selectionAuthority=authority(r,*review,*i,true);result->state=SectionProposalStatus::ready;result->problem=SectionProposalIssue::none;
    result->title="Selected passing alternative - "+juce::String(proposal::eventCoverageName(i->coverage.category));
    result->reason="Explicit selection, not an artistic recommendation. Current is unchanged until stopped-only Apply; independent measurement remains required.";
    result->narrative=result->title+"\n"+result->reason+"\n"+review->text();result->explanation=result->narrative+"\n"+r.policy.text();
    const auto error=validateEventRegionProposalForApply(r,*result);if(error.isNotEmpty())return {{},ErrorCode::staleRevision,error};
    return {std::shared_ptr<const SectionProposalResult>(result),ErrorCode::none,"Measured alternative selected; audio and measurements unchanged"};
}
}

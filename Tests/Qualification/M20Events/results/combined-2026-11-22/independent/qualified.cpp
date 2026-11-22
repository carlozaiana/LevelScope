#include "CompositionTestFixtures.h"
#include "Standalone/StandaloneRegionProposal.h"
#include "EventAllocationProbe.h"
#include <sys/resource.h>

using namespace levelscope::proposal;
using J = juce::var;
using F = juce::File;
using S = juce::String;
using Clock = std::chrono::steady_clock;

namespace {
J object() { return J(new juce::DynamicObject()); }
void put(J& target, const juce::Identifier& name, J value) { target.getDynamicObject()->setProperty(name, value); }
double seconds(Clock::time_point from) { return std::chrono::duration<double>(Clock::now()-from).count(); }
J optional(std::optional<double> value) { return value ? J(*value) : J(); }
template<typename T> J array(const std::vector<T>& values) {
    juce::Array<J> result;
    for (auto value : values) result.add(J(value));
    return result;
}
S json(const J& value) { return juce::JSON::toString(value, false, 17); }
J metrics(const measurement::Result& m) {
    J j=object();
    put(j,"il_lufs",optional(m.integrated.value)); put(j,"tp_dbtp",optional(m.truePeak.value));
    put(j,"lra_lu",optional(m.lra.value)); put(j,"maximum_m_lufs",optional(m.momentaryMax.value));
    put(j,"maximum_s_lufs",optional(m.shortTermMax.value)); put(j,"sample_peak_dbfs",optional(m.samplePeak.value));
    put(j,"p10",optional(m.p10)); put(j,"p95",optional(m.p95));
    put(j,"il_relative_gate",optional(m.integratedRelativeGate)); put(j,"lra_relative_gate",optional(m.lraRelativeGate));
    put(j,"il_included",(int)m.integratedIncluded); put(j,"lra_included",(int)m.lraIncluded);
    put(j,"sample_peaks_linear",array(m.channelSamplePeaks)); put(j,"sample_peak_frames",array(m.samplePeakFrames));
    put(j,"true_peaks_linear",array(m.channelTruePeaks)); put(j,"true_peak_eighth_frames",array(m.truePeakEighthFrames));
    return j;
}
J roles(const ChannelLayout& l) {
    juce::Array<J> result;
    for (int i=0;i<l.count;++i) result.add(juce::AudioChannelSet::getAbbreviatedChannelTypeName(l.roles[i]));
    return result;
}
J plan(const ProcessingPlan& p) {
    J j=object(); put(j,"identity",proposalPlanIdentity(p)); put(j,"global_db",p.staticGainDb);
    juce::Array<J> spans;
    if (p.sectionGains) {
        put(j,"gain_audio_key",p.sectionGains->audioKey());
        for (const auto& part:p.sectionGains->canonical().intervals) {
            J row=object(); put(row,"begin_frame",part.begin); put(row,"end_frame_exclusive",part.end); put(row,"plateau_db",part.db); spans.add(row);
        }
    }
    put(j,"intervals",spans); return j;
}
void selectAndVerify(const StandaloneSessionModel& baseline, const std::shared_ptr<const SectionProposalResult>& result,
                     size_t index, const F& dir, J& row) {
    auto s=baseline; JobControl job;
    auto review=result->regionReview();
    const auto& a=review->attempts.at(index);
    expect(s.acceptSectionProposal(result),"accept immutable applicable result in disposable copy");
    if (review->selectedAttempt!=index || !s.canApplySectionProposal())
        expect(s.selectMeasuredRegionAlternative(review->resultIdentity,a.attemptIdentity),"explicit selection of exact measured pass");
    expect(s.canApplySectionProposal()&&s.applySectionProposal()&&!s.currentState.hasMeasurement,"stopped Apply invalidates Current");
    const auto fresh=measuredComposition(s,MeasurementDomain::currentState);
    expect(s.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass,"independent receiving Current passes selected requirements");
    expect(json(metrics(fresh.detail))==json(metrics(a.receipt->measured())),"same-platform Current reproduces all measured scalars, gates and per-channel peaks and locations");
    put(row,"independent_current",metrics(fresh.detail));
    expect(sf::save(sf::capture(s),dir.getChildFile("applied.levelscope"),job).succeeded(),"save independent passing Current");
    const auto request=ex::capture(s); auto artifact=ex::render(request,dir,job);
    if(!artifact.value) std::cerr<<artifact.message<<std::endl;
    expect(bool(artifact.value),"finalized PCM24 rendered and independently read back");
    auto& c=**artifact.value; const bool passed=ex::artifactEvaluation(c,request).verdict==Verdict::pass;
    put(row,"final_pcm_pass",passed); put(row,"final_pcm_sha256",c.hash()); put(row,"final_pcm",metrics(c.analysis().detail));
    expect(c.analysis().context.processing.staticGainDb==0,"readback uses unity processing");
    if (passed) {
        auto saved=ex::publish(c,request,dir.getChildFile("verified-pcm24"),ex::approvalIdentity(c,request),job);
        expect(bool(saved.value),"retain verified finalized export in local evidence folder");
    } else {
        expect(c.file().copyFileTo(dir.getChildFile("failed-readback.wav")),"retain finalized nonpassing PCM as diagnostic evidence");
    }
    expect(passed,"finalized file independently meets selected rules");
}
void cancellation(const SectionProposalRequest& request,J& out) {
    JobControl control; std::atomic<bool> done{false}; Clock::time_point requested; bool issued=false;
    std::thread stop([&] {
        while(!done.load()) {
            if (control.progress.load()>=.15f) { requested=Clock::now(); issued=true; control.cancel.store(true); return; }
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
    });
    auto result=generateSectionProposal(request,control);
    auto ended=Clock::now(); done=true; stop.join();
    expect(issued&&!result.value&&result.error==ErrorCode::cancelled,"cancellation during auxiliary pass discards result");
    const double elapsed=std::chrono::duration<double>(ended-requested).count();
    put(out,"cancel_seconds_after_request",elapsed);
    expect(elapsed>=0&&elapsed<=1,"auxiliary cancellation completes within one second");
}
void holdout(const F& audio,const F& root,J& out) {
    const bool stereo=audio.getFileName()=="DialogAmbiance.wav";
    expect(stereo||audio.getFileName()=="5.1_Busy.wav","one of the two predeclared files");
    const S expected=stereo?"602485eef95774826a575b0069d358e6723a31f6928dbd6e04dd0df44c98eb90":"4fa30bfd9a8fbb64e68fc00eaf11f8adf81ead8502d05ce639b22cc07645d7b6";
    JobControl job; auto media=validateMedia(audio,job);
    expect(media.succeeded()&&media.media.sha256==expected&&media.media.frames==(stereo?7248001:8923042)&&media.media.sampleRate==48000&&media.media.channels==(stereo?2:6),"exact precommitted PCM, channels, length and rate");
    StandaloneSessionModel s; expect(s.setValidatedSource(audio,media.media),"adopt full source");
    s.sourceDocument.layout=layoutFor(media.media.channels); s.setSelectedTargetProfileIndex(1);
    s.setProgrammeApplicability(Applicability::preProducedLinearProgramme);
    auto baseline=measuredComposition(s,MeasurementDomain::source); s.synchronizeSections();
    expect(s.currentState.hasMeasurement&&s.currentState.analysis->identicalSourceMeasurementIdentity==baseline.identity,"unity Current uses identical Source");
    put(out,"file",audio.getFileName()); put(out,"sha256",expected); put(out,"roles",roles(s.sourceDocument.layout));
    put(out,"frames",media.media.frames); put(out,"baseline",metrics(baseline.detail)); put(out,"platform",baseline.context.platformAndFloatingPointMode);
    if(!baseline.detail.lraStable||!baseline.detail.lra.value||std::floor((*baseline.detail.lra.value-1)*10)<=0) {
        put(out,"outcome","unassessed: no stable positive reduced LRA ceiling"); return;
    }
    auto requirements=s.jobRequirements; requirements.integratedTenths=-230; requirements.integratedToleranceTenths=2;
    requirements.truePeakTenths=-10; requirements.maximumLraTenths=int(std::floor((*baseline.detail.lra.value-1)*10));
    requirements.maximumMomentaryTenths.reset(); requirements.maximumShortTermTenths.reset();
    requirements.provenance="Independent engineering comparison; no delivery-compliance claim";
    S error; expect(s.setJobRequirements(requirements,error),"install predeclared holdout target");
    const auto policy=s.proposalPolicy;
    expect(policy.durationThousandths==200&&policy.additionalTenths==30&&policy.cumulativeTenths==120&&policy.reserveTenths==5&&!policy.regionLimit,"default preservation policy exact");
    auto settings=RegionPlanningSettings::innerDefaults(RegionMode::automatic,sectionInputIdentity(s.measurementContext(MeasurementDomain::source)));
    settings.allowAutomaticMargins=true; settings.leftMarginMs=settings.rightMarginMs=1000;
    expect(settings.schema==3&&settings.allowGlobal&&s.setRegionPlanning(settings,error),"schema3 automatic and permitted 1s margins");
    put(out,"requirements",requirements.encode()); put(out,"policy",policy.encode()); put(out,"settings",settings.encode());
    auto request=s.sectionProposalRequest(); const auto unchanged=s.currentState.analysis;
    const auto originalPlan=proposalPlanIdentity(s.currentStateDocument.processingPlan());
    const auto start=Clock::now(); event_allocation::Scope alloc; auto r=generateSectionProposal(request,job);
    const auto heap=alloc.stop(); put(out,"assessment_seconds",seconds(start)); put(out,"peak_heap_bytes",juce::int64(heap));
    if(!r.value) std::cerr<<r.message<<std::endl;
    expect(bool(r.value),"completed honest bounded search"); auto result=*r.value;
    expect(s.currentState.analysis==unchanged&&proposalPlanIdentity(s.currentStateDocument.processingPlan())==originalPlan,"search preserves Current and baseline processing");
    auto review=result->regionReview(); expect(bool(review),"retained region review");
    put(out,"outcome",result->outcomeTitle()); put(out,"reason",result->outcomeReason()); put(out,"private_passes",result->passes());
    put(out,"proxy_count",int(review->proxyEvaluations)); put(out,"planning_bytes",juce::int64(review->planningBytes));
    put(out,"planning_seconds",review->planningSeconds); put(out,"sidecar_seconds",review->sidecarSeconds);
    put(out,"feature_high_water_bytes",juce::int64(review->featureHighWaterBytes));
    expect(root.getChildFile("proposal.txt").replaceWithText(result->text()),"retain complete measured proposal and boundary explanations");
    expect(alloc.valid()&&heap<=64*1024*1024&&review->planningBytes<=64*1024*1024&&review->proxyEvaluations<=64&&result->passes()<=3,"concurrent heap, planning storage and search counts within frozen caps");
    juce::Array<J> attempts; int passed=0;
    for(size_t i=0;i<review->attempts.size();++i) {
        const auto& a=review->attempts[i]; J row=object(); put(row,"attempt",int(i+1)); put(row,"family",a.family);
        put(row,"coverage",eventCoverageName(a.coverage.category)); put(row,"plan_identity",a.plan.identity);
        put(row,"affected_frames",a.plan.activity.affectedFrames); put(row,"pass",a.evaluation.verdict==Verdict::pass);
        expect(bool(a.receipt),"retained actual attempt owns measurement provenance");
        put(row,"metrics",metrics(a.receipt->measured())); put(row,"plan",plan(a.receipt->context().processing));
        if(a.evaluation.verdict==Verdict::pass) {
            ++passed; const auto dir=root.getChildFile("attempt-"+S(int(i+1)));
            expect(dir.createDirectory().wasOk(),"fresh attempt output"); attempts.add(row); put(out,"attempts",attempts);
            selectAndVerify(s,result,i,dir,row);
        } else attempts.add(row);
        put(out,"attempts",attempts);
    }
    put(out,"measured_passing_attempts",passed);
    if(!passed) {
        expect(s.acceptSectionProposal(result)&&!s.canApplySectionProposal()&&!s.applySectionProposal(),"no measured pass grants Apply");
        for(const auto& a:review->attempts) expect(!s.selectMeasuredRegionAlternative(review->resultIdentity,a.attemptIdentity),"failed candidate cannot be selected");
        expect(sf::save(sf::capture(s),root.getChildFile("baseline.levelscope"),job).succeeded(),"save verified neutral baseline for review");
    }
    if(stereo) {
        cancellation(request,out);
        expect(s.currentState.analysis==unchanged&&proposalPlanIdentity(s.currentStateDocument.processingPlan())==originalPlan,"cancel leaves processing and Current unchanged");
    }
    expect(juce::SHA256(audio).toHexString()==expected,"input unchanged after all offline checks");
}
void exchange(const F& session,const F& audio,const F& root,bool requireForeign,J& out) {
    JobControl job; auto input=sf::read(session,job);
    expect(input.document&&input.document->loadedSchema==11&&input.document->lastAppliedExplanation&&input.document->regionPlanning&&input.document->sourceCache&&input.document->currentCache,"read document11/programme9 with applied history and both origin caches");
    const auto before=*input.document;
    const auto history=encodeAppliedExplanation(*before.lastAppliedExplanation);
    expect(history.size()>8&&history.front()==4&&std::all_of(history.begin()+1,history.begin()+8,[](auto b){return b==0;}),"fourth applied explanation attachment format");
    const auto settings=before.regionPlanning->digest(), geometry=before.sections.digest();
    const auto processing=proposalPlanIdentity(before.current.processingPlan());
    const auto hash=before.source.media.sha256; const auto original=juce::SHA256(session).toHexString();
    expect(before.regionPlanning->schema==3&&before.regionPlanning->mode==RegionMode::automatic&&before.regionPlanning->allowAutomaticMargins&&before.regionPlanning->leftMarginMs==1000&&before.regionPlanning->rightMarginMs==1000,"expected competing method and permissions");
    expect(audio.existsAsFile()&&juce::SHA256(audio).toHexString()==hash,"receiving audio matches source content");
    expect(sf::relinkAudio(*input.document,audio),"explicitly substitute verified locator only in memory");
    auto verified=sf::verify(*input.document,job); StandaloneSessionModel s; S notice;
    expect(verified.prepared&&sf::adopt(*verified.document,*verified.prepared,s,notice),"normal session preparation and adoption");
    const auto context=s.measurementContext(MeasurementDomain::source);
    const auto origin=before.sourceCache->context.platformAndFloatingPointMode;
    expect((origin!=context.platformAndFloatingPointMode)==requireForeign,"actual sender/receiver platform relationship as requested");
    expect(context.buildVersion==LEVELSCOPE_MEASUREMENT_BUILD&&before.currentCache->context.buildVersion==context.buildVersion&&before.sourceCache->context.buildVersion==context.buildVersion,"same frozen measurement build on both ends");
    expect(before.currentCache->context.platformAndFloatingPointMode==origin,"same origin for Source and Current");
    if(requireForeign)
        expect(!s.source.hasMeasurement&&!s.currentState.hasMeasurement&&s.currentStateNeedsRemeasurement(),"foreign analysis discarded until independent measurements");
    auto invariants=[&](const StandaloneSessionModel& m) {
        expect(!m.sectionProposal&&!m.canApplySectionProposal(),"reopened document grants no live proposal authority");
        expect(m.lastAppliedExplanation&&encodeAppliedExplanation(*m.lastAppliedExplanation)==history,"exact history bytes preserved");
        expect(m.regionPlanning&&m.regionPlanning->digest()==settings&&proposalPlanIdentity(m.currentStateDocument.processingPlan())==processing&&m.sectionDocument.map().digest()==geometry,"exact settings, processing and annotation geometry");
        expect(m.jobRequirements.digest()==before.jobRequirements.digest()&&m.proposalPolicy.digest()==before.proposalPolicy.digest(),"target and preservation policy identical");
        expect(m.sourceDocument.layout.confirmed&&m.sourceDocument.layout.count==before.source.layout.count&&m.sourceDocument.layout.roles==before.source.layout.roles&&m.sourceDocument.layout.manuallyInterpreted==before.source.layout.manuallyInterpreted,"same confirmed ordered channel assignment");
        expect(std::all_of(m.sourceDocument.routing.fileToDevice.begin(),m.sourceDocument.routing.fileToDevice.end(),[](int n){return n==-1;}),"all physical routes muted");
    };
    invariants(s);
    put(out,"origin",origin); put(out,"receiver",context.platformAndFloatingPointMode);
    put(out,"source_cache_adopted",s.source.hasMeasurement); put(out,"current_cache_adopted",s.currentState.hasMeasurement);
    put(out,"source_sha256",hash); put(out,"input_session_sha256",original);
    put(out,"settings",s.regionPlanning->encode()); put(out,"settings_digest",settings);
    put(out,"requirements",s.jobRequirements.encode()); put(out,"policy",s.proposalPolicy.encode()); put(out,"roles",roles(s.sourceDocument.layout));
    put(out,"history_sha256",juce::SHA256(history.data(),history.size()).toHexString()); put(out,"sections_digest",geometry);
    put(out,"canonical_processing",plan(s.currentStateDocument.processingPlan()));
    put(out,"origin_source",metrics(before.sourceCache->detail)); put(out,"origin_current",metrics(before.currentCache->detail));
    const auto source=measuredComposition(s,MeasurementDomain::source); const auto current=measuredComposition(s,MeasurementDomain::currentState);
    put(out,"receiver_source",metrics(source.detail)); put(out,"receiver_current",metrics(current.detail));
    expect(s.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass,"independent receiving Current passes declared rules");
    const auto dest=root.getChildFile("receiver.levelscope");expect(sf::save(sf::capture(s),dest,job).succeeded(),"save independently measured receiving document");
    auto reopened=sf::read(dest,job);expect(reopened.document&&reopened.document->loadedSchema==11,"reread receiving schema11");
    auto prepared=sf::verify(*reopened.document,job); StandaloneSessionModel same;
    expect(prepared.prepared&&sf::adopt(*prepared.document,*prepared.prepared,same,notice),"normal same-platform reopen");
    invariants(same);
    expect(same.source.hasMeasurement&&same.currentState.hasMeasurement&&same.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass,"receiving caches valid after independent remeasure and reopen");
    expect(same.currentState.analysis->identity==current.identity&&same.source.analysis->identity==source.identity,"receiving platform cache identities retained");
    expect(juce::SHA256(session).toHexString()==original&&juce::SHA256(audio).toHexString()==hash,"exchange source session and audio byte-exact after test");
    put(out,"saved_session_sha256",juce::SHA256(dest).toHexString());
}
}
int main(int argc,char** argv) {
    if(argc<4)return 2;
    const S mode(argv[1]); bool independent=mode=="holdout";
    if(independent?argc!=4:(argc!=5||(mode!="receive"&&mode!="local-control")))return 2;
    const F root(argv[argc-1]);
    if(root.exists()||!root.createDirectory().wasOk()){std::cerr<<"Output must be new"<<std::endl;return 2;}
    juce::ScopedJuceInitialiser_GUI init; J out=object(); bool success=false;
    try {
        put(out,"mode",mode); put(out,"measurement_build",LEVELSCOPE_MEASUREMENT_BUILD);
        if(independent)holdout(F(argv[2]),root,out);
        else exchange(F(argv[2]),F(argv[3]),root,mode=="receive",out);
        success=true;
    } catch(const std::exception& e) {std::cerr<<e.what()<<std::endl;put(out,"error",S(e.what()));}
    rusage usage{}; getrusage(RUSAGE_SELF,&usage); put(out,"process_max_rss_bytes",juce::int64(usage.ru_maxrss));
    put(out,"checks",checks); put(out,"passed",success);
    if(!root.getChildFile("result.json").replaceWithText(json(out)))return 3;
    return success?0:1;
}

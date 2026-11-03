#include "Trace.h"
#include "Standalone/StandaloneSessionModel.h"
#include "Standalone/ProposalWorkflowPanel.h"
#include "Standalone/TimelineGainLane.h"
#include "Standalone/MeasurementTimeline.h"
#include <ctime>
#include <future>
#include <iostream>
#include <pthread.h>
#include <juce_cryptography/juce_cryptography.h>
#include <mach/mach.h>

using namespace levelscope;
using namespace levelscope::standalone;
using namespace lsprofile;
namespace
{
juce::var object() { return juce::var(new juce::DynamicObject()); }
void put(juce::var& o,const juce::Identifier& key,const juce::var& v) { o.getDynamicObject()->setProperty(key,v); }
void require(bool ok,const juce::String& detail) { if(!ok) throw std::runtime_error(detail.toStdString()); }
juce::String peaksDigest(const proposal::SectionPeakMap& map)
{
    juce::MemoryOutputStream bytes;
    bytes.writeInt(map.channels); bytes.writeInt64(map.sections.size());
    for(const auto& row:map.sections)
        for(const auto* column:{&row.sample,&row.truePeak})
            for(const auto& peak:*column)
            {
                bytes.writeDouble(peak.magnitude); bytes.writeInt64(peak.eighthFrame); bytes.writeBool(peak.padding);
            }
    return juce::SHA256(bytes.getData(),bytes.getDataSize()).toHexString();
}
double threadCpuMs()
{
    thread_basic_info_data_t info{};
    mach_msg_type_number_t n=THREAD_BASIC_INFO_COUNT;
    require(thread_info(pthread_mach_thread_np(pthread_self()),THREAD_BASIC_INFO,
                        reinterpret_cast<thread_info_t>(&info),&n)==KERN_SUCCESS,"thread CPU counter");
    return 1000.*(info.user_time.seconds+info.system_time.seconds)+(info.user_time.microseconds+info.system_time.microseconds)/1000.;
}

class Heartbeat : private juce::Timer
{
public:
    explicit Heartbeat(std::function<void()> render):paint(std::move(render)),last(Clock::now()) { startTimerHz(10); }
    void stop(juce::var& row)
    {
        stopTimer();
        put(row,"heartbeat_gaps_ms",gaps);
        put(row,"heartbeat_tail_ms",ms(last,Clock::now()));
        put(row,"timeline_paint_ms",paints);
    }
private:
    std::function<void()> paint;
    Clock::time_point last;
    juce::Array<juce::var> gaps,paints;
    void timerCallback() override
    {
        const auto now=Clock::now();
        gaps.add(ms(last,now)); last=now;
        paint(); paints.add(ms(now,Clock::now()));
    }
};

struct Trial
{
    StandaloneSessionModel s;
    WaveformCache waves;
    MeasurementTimeline timeline;
    TimelineGainLane lane{s};
    ProposalWorkflowPanel panel{s};
    std::function<void()> updatePreview;
    juce::var operations=object(), identities=object(), outcomes=object();
    juce::StringArray roles;
    juce::String sha;
    void render()
    {
        if(updatePreview) updatePreview();
        lane.synchronize(sha); panel.refresh();
        auto image=timeline.createComponentSnapshot(timeline.getLocalBounds());
        require(image.isValid(),"timeline paint failed");
    }
    template<class F> auto timeJob(const juce::String& name,bool background,F&& f)
    {
        juce::var row=object();
        auto run=[&]{
            if(background) pthread_setname_np(("ls-"+name).toRawUTF8());
            values={}; phase=0;
            const auto cpu=std::clock();
            const auto threadCpu=threadCpuMs();
            const auto begin=Clock::now();
            auto result=f();
            put(row,"wall_ms",ms(begin,Clock::now()));
            put(row,"worker_cpu_ms",threadCpuMs()-threadCpu);
            put(row,"process_cpu_ms",1000.*(std::clock()-cpu)/CLOCKS_PER_SEC);
            auto timings=object();
            for(size_t i=0;i<values.size();++i)
            {
                auto group=object(); bool present=false;
                for(size_t k=0;k<count;++k) if(values[i][k]>0)
                {put(group,names[k],values[i][k]);present=true;}
                if(present) put(timings,phases[i],group);
            }
            put(row,"timings_ms",timings);
            return result;
        };
        Heartbeat heartbeat([&]{render();});
        auto result=[&]{
            if(!background) return run();
            auto worker=std::async(std::launch::async,run);
            while(worker.wait_for(std::chrono::seconds(0))!=std::future_status::ready)
                juce::MessageManager::getInstance()->runDispatchLoopUntil(5);
            return worker.get();
        }();
        heartbeat.stop(row);
        put(operations,name,row);
        std::cout << name << " " << row["wall_ms"].toString() << " ms\n" << std::flush;
        return result;
    }

    void measure(MeasurementDomain domain,const juce::String& name)
    {
        require(s.job.begin(JobStage::analysis,s.revision),"begin analysis");
        const auto c=s.measurementContext(domain);
        const auto binding=WaveformBinding::from(c,domain);
        auto wave=waves.begin(binding);
        LoudnessProgress loudness(binding);
        updatePreview=[&]{
            const auto view=waves.view(binding);
            timeline.setWaveform(view.data,view.accepted,roles);
            timeline.setPreview(loudness.view(binding));
        };
        JobControl job;
        auto result=timeJob(name,true,[&]{return OfflineAnalysisService(nullptr,wave.get(),&loudness).analyze(c,domain,job);});
        updatePreview={}; timeline.setPreview({});
        require(result.value.has_value(),result.message);
        require(s.job.finish(c.ticket,s.revision,result.error,result.message) && s.acceptMeasurement(std::move(*result.value)),"measurement acceptance");
        require(waves.accept(binding,binding),"overview missing");
        const auto view=waves.view(binding); timeline.setWaveform(view.data,true,roles);
        const auto a=domain==MeasurementDomain::source?s.source.analysis:s.currentState.analysis;
        put(identities,name,measuredProgrammeForSnapshot(a,c,domain,s.programmeApplicability).snapshotIdentity);
        if(domain==MeasurementDomain::source)
        {
            require(s.currentState.analysis && s.currentState.analysis->identicalSourceMeasurementIdentity.isNotEmpty(),"initial Current origin");
            require(waves.shareInitialCurrent(*a,*s.currentState.analysis),"shared initial overview");
        }
    }

    void propose(const juce::String& name,int passes,SectionProposalStatus expected)
    {
        s.synchronizeSections();
        const auto request=s.sectionProposalRequest();
        const auto pre=Clock::now();
        const auto first=proposal::draftGlobalTrim(request.scalar());
        const double dispatchMs=ms(pre,Clock::now());
        const bool background=first.status==proposal::TrimStatus::noSolution;
        s.discardSectionProposal();
        if(background) require(s.job.begin(JobStage::proposal,s.revision),"begin proposal");
        JobControl job;
        auto res=timeJob(name,background,[&]{return generateSectionProposal(request,job);});
        if(background) require(s.job.finish(s.job.ticket,s.revision,res.error,res.message),res.message);
        require(res.value.has_value(),res.message);
        const auto acceptStart=Clock::now();
        require(s.acceptSectionProposal(*res.value),"proposal acceptance");
        auto out=*res.value;
        require(out->passes()==passes && out->status()==expected,out->text());
        auto fields=object();
        put(fields,"main_preflight_ms",dispatchMs);put(fields,"main_accept_ms",ms(acceptStart,Clock::now()));
        put(fields,"passes",out->passes());put(fields,"outcome",out->outcomeTitle());put(fields,"changes",int(out->selection().corrections.size()));
        put(fields,"status",int(out->status()));put(fields,"issue",int(out->issue()));put(fields,"can_apply",s.canApplySectionProposal());
        put(fields,"target_il",request.currentContext.profile.limits.integratedLufs);
        put(fields,"target_tp",request.currentContext.profile.limits.truePeakCeilingDbtp);
        put(outcomes,name,fields);
        put(identities,name,out->identity());
        put(identities,name+"_decision",juce::SHA256(out->text().toRawUTF8(),out->text().getNumBytesAsUTF8()).toHexString());
        if(out->peaks()) put(identities,name+"_section_peaks",peaksDigest(*out->peaks()));
        if(out->assessment())
        {
            auto a=out->assessment();
            put(identities,name+"_candidate",measuredProgrammeForSnapshot(a,a->context,MeasurementDomain::privateCandidate,s.programmeApplicability).snapshotIdentity);
            put(identities,name+"_plan",proposalPlanIdentity(out->plan()));
        }
    }

    juce::var run(const juce::File& file,bool decisionChecks)
    {
        s.sessionId="11111111111111111111111111111111";
        timeline.setSize(1160,460); timeline.attachGainLane(lane); panel.setSize(1160,ProposalWorkflowPanel::height);
        JobControl job;
        auto v=timeJob("import",true,[&]{return validateMedia(file,job);});
        require(v.succeeded() && s.setValidatedSource(file,v.media),v.message);
        require((v.media.channels==2||v.media.channels==6) && v.media.frames==60*48000 && v.media.bitsPerSample==24,"unexpected fixture");
        sha=v.media.sha256;
        auto order=assemblyRoles(v.media.channels);
        for(int i=0;i<v.media.channels;++i){s.setChannelRole(i,order.at(size_t(i)));roles.add(juce::AudioChannelSet::getAbbreviatedChannelTypeName(order.at(size_t(i))));}
        require(s.confirmLayout(),"confirm speaker roles");
        s.setSelectedTargetProfileIndex(1); s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
        s.synchronizeSections();
        juce::String error;
        for(int n=1;n<20;++n) require(s.sectionDocument.split(uint64_t(n),v.media.frames*n/20,error),error);
        timeline.setProgramme(v.media.frames,48000,sha);
        measure(MeasurementDomain::source,"source");
        auto normal=s.jobRequirements;
        auto fast=proposal::JobRequirements::editable(s.getSelectedTargetProfile());
        fast.id="22222222222222222222222222222222";
        fast.integratedTenths=juce::roundToInt(10*(*s.source.integratedLufs+.7));
        require(s.setJobRequirements(fast,error),error);
        propose("global",0,SectionProposalStatus::globalOnly);
        require(s.setJobRequirements(normal,error),error);
        propose("fresh",3,SectionProposalStatus::ready);
        if(decisionChecks)
        {
            const auto baseline=s.currentState.analysis;
            const auto policy=s.proposalPolicy;
            auto restrictive=policy;
            restrictive.durationThousandths=0;
            require(s.setProposalPolicy(restrictive,error),error);
            require(s.sectionProposalStale() && !s.applySectionProposal(),"stale proposal was actionable");
            propose("preference_blocked",0,SectionProposalStatus::noSolution);
            require(!s.canApplySectionProposal() && s.sectionProposal->issue()==SectionProposalIssue::preferences,"duration gate");
            require(s.setProposalPolicy(policy,error),error);
            auto strict=proposal::JobRequirements::editable(s.getSelectedTargetProfile());
            strict.id="33333333333333333333333333333333"; strict.integratedTenths=-180;
            require(s.setJobRequirements(strict,error),error);
            propose("measured_rejection",1,SectionProposalStatus::noSolution);
            require(!s.canApplySectionProposal() && s.sectionProposal->assessment()
                    && s.sectionProposal->issue()==SectionProposalIssue::candidateRejected,"candidate rejection");
            strict.maximumLraTenths=180;
            require(s.setJobRequirements(strict,error),error);
            propose("unsupported",0,SectionProposalStatus::unavailable);
            require(!s.canApplySectionProposal() && s.currentState.analysis==baseline,"rejection changed Current");
            require(s.setJobRequirements(normal,error),error);
        }
        propose("reused",2,SectionProposalStatus::ready);
        require(identities["fresh_plan"]==identities["reused_plan"],"reused attribution changed the applied plan");
        require(s.applySectionProposal() && s.currentStateNeedsRemeasurement(),"explicit Apply");
        measure(MeasurementDomain::currentState,"current");
        require(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"independent Current verdict");
        if(decisionChecks) propose("already_meets",0,SectionProposalStatus::noChange);
        auto result=object();
        put(result,"operations",operations);put(result,"identities",identities);put(result,"outcomes",outcomes);
        put(result,"roles",roles.joinIntoString(","));put(result,"file_sha256",sha);put(result,"frames",v.media.frames);
        put(result,"source_lufs",*s.source.integratedLufs);put(result,"current_lufs",*s.currentState.integratedLufs);
        put(result,"channels",v.media.channels);put(result,"base_build",LEVELSCOPE_MEASUREMENT_BUILD);
        return result;
    }
};
}
int main(int argc,char** argv)
{
    try
    {
        require(argc==3 || (argc==4 && std::string(argv[3])=="--decisions"),"usage: profile INPUT.wav OUTPUT.json [--decisions]");
        juce::ScopedJuceInitialiser_GUI gui;
        Trial trial;
        auto result=trial.run(juce::File(argv[1]),argc==4);
        const juce::File dest(argv[2]); require(!dest.exists() && dest.replaceWithText(juce::JSON::toString(result)),"write fresh report");
        return 0;
    }
    catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}

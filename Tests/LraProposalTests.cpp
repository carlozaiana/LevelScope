#include "CompositionTestFixtures.h"
#include "Standalone/StandaloneLraProposal.h"
#include "Standalone/LraPassagesPanel.h"
#include "Standalone/AppliedProposalExplanation.h"
#include "Core/Measurement/GatePopulation.h"
#include <future>
#include "LraListeningChecks.h"

namespace
{
juce::var resources(new juce::DynamicObject());
StandaloneSessionModel lraFixture(const juce::File& root,int channels,int rate=48000,int seconds=65)
{
    expect(root.createDirectory().wasOk(),"owned LRA fixture directory");auto file=root.getChildFile("programme.wav");
    auto stream=file.createOutputStream();juce::WavAudioFormat format;
    std::unique_ptr<juce::AudioFormatWriter> writer(format.createWriterFor(stream.get(),rate,unsigned(channels),24,{},0));expect(bool(writer),"LRA PCM24 fixture writer");stream.release();
    juce::AudioBuffer<float> block(channels,4096);const int frames=rate*seconds;
    bool written=true;for(int first=0;first<frames;first+=4096)
    {
        const int count=std::min(4096,frames-first);
        for(int ch=0;ch<channels;++ch)for(int n=0;n<count;++n)
        {
            const int f=first+n;const double high=f>=rate*30&&f<rate*36?2.:1.;
            const double amplitude=channels==6&&ch==3?.12:.03*(1+.07*ch);
            block.setSample(ch,n,float(amplitude*high*std::sin(2*juce::MathConstants<double>::pi*997*f/rate)));
        }
        written&=writer->writeFromAudioSampleBuffer(block,0,count);
    }
    expect(written,"complete LRA fixture payload written");
    writer.reset();JobControl control;const auto v=validateMedia(file,control);expect(v.succeeded(),"LRA fixture validated");
    StandaloneSessionModel s;expect(s.setValidatedSource(file,v.media),"LRA source adopted");s.sourceDocument.layout=layoutFor(channels);
    s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
    measuredComposition(s,MeasurementDomain::source);expect(s.initializeCurrentStateFromSource(),"LRA Current initialized");
    if(!s.currentState.hasMeasurement)measuredComposition(s,MeasurementDomain::currentState);
    s.synchronizeSections();return s;
}
void setLra(StandaloneSessionModel& s,int tenths)
{auto r=s.jobRequirements;r.maximumLraTenths=tenths;r.provenance="M19 synthetic test requirement, not a delivery recommendation";juce::String error;expect(s.setJobRequirements(r,error),"explicit LRA requirement adopted without remeasurement");}
void requireLraReduction(StandaloneSessionModel& s)
{expect(s.currentState.analysis->detail.lra.value.has_value(),"fixture LRA available");setLra(s,int(std::floor((*s.currentState.analysis->detail.lra.value-1.)*10)));}
void pureCases()
{
    auto a=std::make_shared<proposal::ProgramAnalysis>();a->identity="test-population";a->snapshotIdentity="snapshot";a->sampleRate=48000;a->frames=48000*65;a->p10=-30;a->p95=-20;a->lra={10,measurement::Status::available};a->lraStable=true;
    for(int n=0;n<620;++n){proposal::DiagnosticWindow w;w.begin=n*4800;w.end=w.begin+144000;w.lufs=-30;w.included=true;if(n==80||n==200)w.lufs=-20;a->lraWindows.push_back(w);}
    auto p=proposal::findLraPassages(a,5);expect(p&&p->high.size()==2&&p->high[0].begin<p->high[1].begin,"equal score/count regions break ties by earliest exact source frame");
    expect(p->high[0].score==5&&p->high[0].windows==1&&p->high[0].end-p->high[0].begin==144000,"full support and LU-times-window score explicit");
    auto padded=*a;padded.lraWindows.back().lufs=-20;padded.lraWindows.back().end=a->frames+4800;padded.lraWindows.back().padding=true;
    auto pad=proposal::findLraPassages(std::make_shared<proposal::ProgramAnalysis>(padded),5);
    expect(std::any_of(pad->high.begin(),pad->high.end(),[&](const auto& r){return r.end==a->frames&&r.paddedWindows==1;}),"padding labelled without extending programme duration");
    auto excluded=*a;excluded.lraWindows[80].included=false;excluded.lraWindows[200].lufs.reset();
    auto missing=proposal::findLraPassages(std::make_shared<proposal::ProgramAnalysis>(excluded),5);
    expect(missing->high.empty()&&missing->excluded==1&&missing->unavailable==1,"excluded and unavailable windows never become zero-valued contributions");
    expect(!proposal::findLraPassages(a,5,[]{return true;}),"population cancellation discards view");
    processing::SectionGainPlan before;before.inputIdentity=juce::String::repeatedString("a",64);before.sampleRate=48000;before.frames=a->frames;before.intervals={{0,before.frames,0}};
    proposal::ProposalPolicy policy;policy.regionLimit=1;auto selected=proposal::selectLraTrims(*p,before,policy);
    expect(selected.available&&selected.selected.size()==1&&selected.selected[0].passage.begin==p->high[0].begin,"greedy deterministic selection obeys actual region-count cap");
    expect(selected.selected[0].reduction==3&&selected.selected[0].uncovered==2&&selected.selected[0].reserve==0,"additional cap and uncovered overshoot explicit, reserve remains inside cap");
    auto refinement=proposal::refineLraTrims(selected,1,policy);expect(!refinement.available,"saturated caps do not launch duplicate refinement");
    auto cumulative=before;cumulative.intervals[0].db=-11.5;auto nearCap=proposal::selectLraTrims(*p,cumulative,policy);
    expect(nearCap.available&&nearCap.maximumAdditional==.5&&nearCap.selected.front().reserve==0,"existing attenuation leaves only actual cumulative headroom, with no reserve outside caps");
    cumulative.intervals[0].db=-12;expect(!proposal::selectLraTrims(*p,cumulative,policy).available,"saturated cumulative correction prevents further changes");
    policy.durationThousandths=0;expect(!proposal::selectLraTrims(*p,before,policy).available,"zero duration disables changes without relaxing preference");
    policy={};auto edge=before;edge.frames=10005;edge.intervals={{0,10005,0}};auto after=edge;after.intervals={{0,1334,-1},{1334,10005,0}};
    int64_t affected=0;size_t regions=0;double maximum=0;
    expect(proposal::validateLraEnvelope(edge,after,policy,affected,regions,maximum).isEmpty()&&affected==2001,"exact duration boundary includes transition support");
    --edge.frames;--after.frames;--edge.intervals.back().end;--after.intervals.back().end;
    expect(proposal::validateLraEnvelope(edge,after,policy,affected,regions,maximum).startsWith("Policy:"),"one frame beyond duration rejected without tolerance");
    after=edge;after.intervals={{0,edge.frames,std::nextafter(-3.,-4.)}};policy.durationThousandths=1000;
    expect(proposal::validateLraEnvelope(edge,after,policy,affected,regions,maximum).contains("additional"),"immediately excessive attenuation rejected at full precision");
    after.intervals[0].db=-3;expect(proposal::validateLraEnvelope(edge,after,policy,affected,regions,maximum).isEmpty(),"exact additional gain cap accepted");
    after.intervals[0].db=.1;expect(proposal::validateLraEnvelope(edge,after,policy,affected,regions,maximum).startsWith("Method:"),"attenuation-only method rejects local boost");
    AnalysisSnapshot snapshot;snapshot.detail.lra={6,measurement::Status::available};proposal::Evaluation e;e.verdict=proposal::Verdict::fail;
    proposal::RuleEvaluation lr;lr.id="loudness-range";lr.measured=6;lr.upper=5;lr.verdict=proposal::Verdict::fail;e.rules.push_back(lr);
    auto rank1=lraAttemptRank(snapshot,e,100,1,1);e.rules[0].measured=7;snapshot.detail.lra.value=7;auto rank2=lraAttemptRank(snapshot,e,100,2,2);
    expect(rank1&&rank2&&*rank1<*rank2,"gating-worsened refinement ranks worse and retains earlier actually measured attempt");
    e.rules[0].measured=4;snapshot.detail.lra.value=4;e.verdict=proposal::Verdict::pass;e.rules[0].verdict=proposal::Verdict::pass;
    expect(*lraAttemptRank(snapshot,e,100,2,2)<*rank1,"all-selected measured pass outranks failing attempt");
    const auto first=lraAttemptRank(snapshot,e,100,2,1),later=lraAttemptRank(snapshot,e,100,2,2);expect(*first<*later,"otherwise equal measured attempts retain earlier assessment");
    e.verdict=proposal::Verdict::indeterminate;expect(!lraAttemptRank(snapshot,e,0,0,1),"unavailable metrics cannot win best-candidate selection");
}
void workflow(const juce::File& root,int channels,int rate=48000)
{
    auto s=lraFixture(root,channels,rate);requireLraReduction(s);const auto source=s.source.analysis,current=s.currentState.analysis;const auto revision=s.revision;
    auto request=s.sectionProposalRequest();JobControl control;int preparations=0;SnapshotStorage store;store.directory=root.getChildFile("private");store.directory.createDirectory();store.availableBytes=[&](const juce::File&){++preparations;return int64_t(1)<<40;};
    const auto start=juce::Time::getMillisecondCounterHiRes();auto result=generateSectionProposal(request,control,store);
    if(!result.value)std::cerr<<result.message<<'\n';else std::cout<<(*result.value)->text()<<'\n';
    expect(result.value&&(*result.value)->status()==SectionProposalStatus::ready,"feasible LRA candidate actually passes complete selected requirements");
    auto candidate=*result.value;expect(candidate->passes()>=1&&candidate->passes()<=3&&preparations==1,"at most three full assessments reuse exactly one verified source copy");
    if(channels==2&&rate==48000)boundaryCaptures(root.getChildFile("listening"),request,*candidate);
    expect(store.directory.findChildFiles(juce::File::findFiles,true).isEmpty(),"private verified source cleaned after candidate generation");
    expect(s.source.analysis==source&&s.currentState.analysis==current&&s.revision==revision,"assessment changes no applied state or measurement");
    expect(candidate->plan().staticGainDb>0,"local attenuation followed by explicit positive global gain");
    auto dense=s;juce::String error;for(int n=1;n<65;++n)expect(dense.sectionDocument.split(n,n*rate,error),"additional annotations only");
    JobControl repeat;auto other=generateSectionProposal(dense.sectionProposalRequest(),repeat);
    expect(other.value&&(*other.value)->status()==candidate->status()&&proposalPlanIdentity((*other.value)->plan())==proposalPlanIdentity(candidate->plan())
        &&(*other.value)->assessment()->detail.lra.value==candidate->assessment()->detail.lra.value,"no-markers and dense markers yield identical audio plan and measured outcome");
    expect(s.acceptSectionProposal(candidate)&&s.canApplySectionProposal(),"exact ready candidate accepted with matching request binding");
    auto stale=s;expect(stale.setCurrentGainDb(.1)&&!stale.canApplySectionProposal(),"gain mutation invalidates LRA proposal");
    stale=s;auto requirements=stale.jobRequirements;requirements.maximumLraTenths=*requirements.maximumLraTenths+1;expect(stale.setJobRequirements(requirements,error)&&!stale.canApplySectionProposal(),"requirement change reuses measurement but invalidates candidate");
    stale=s;stale.proposalPolicy.durationThousandths=300;expect(!stale.canApplySectionProposal(),"preference identity invalidates candidate without changing audio");
    stale=s;stale.sessionId="different-runtime";expect(!stale.acceptSectionProposal(candidate),"late result cannot enter another reopened runtime");
    s.setOutputRoute(channels==6?3:0,-1);
    expect(s.canApplySectionProposal(),"monitor muting preserves all-channel candidate identity");
    stale=s;stale.auditionBusy=[]{return true;};expect(!stale.applySectionProposal(),"Apply requires stopped transport");
    expect(s.applySectionProposal()&&!s.currentState.hasMeasurement&&s.currentStateNeedsRemeasurement()&&s.revision==revision+1,"explicit Apply invalidates Current and advances audio revision once");
    const auto plan=s.currentStateDocument.processingPlan();s.measurementUnsuccessful(MeasurementDomain::currentState,ErrorCode::cancelled);
    expect(s.currentStateNeedsRemeasurement()&&proposalPlanIdentity(s.currentStateDocument.processingPlan())==proposalPlanIdentity(plan),"cancelled verification retains gain and unverified state without rollback");
    const auto verified=measuredComposition(s,MeasurementDomain::currentState);
    expect(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass&&verified.detail.integrated.value==candidate->assessment()->detail.integrated.value
        &&verified.detail.lra.value==candidate->assessment()->detail.lra.value&&verified.detail.channelTruePeaks==candidate->assessment()->detail.channelTruePeaks
        &&verified.detail.truePeakEighthFrames==candidate->assessment()->detail.truePeakEighthFrames,"independent Current exactly agrees with candidate metrics and all-channel peak locations");
    expect(s.lastAppliedExplanation&&s.lastAppliedExplanation->method==proposal::lraTrimMethodVersion&&s.lastAppliedExplanation->text().contains("Assessment 1"),"M18 history retains new method, reasons and actual assessment provenance");
    const auto saved=sf::capture(s);JobControl io;auto path=root.getChildFile("verified.levelscope");expect(sf::save(saved,path,io).succeeded(),"save verified LRA processing and historical explanation");
    auto read=sf::read(path,io);expect(bool(read.document),"read saved LRA session");auto reopened=sf::verify(*read.document,io);
    expect(reopened.prepared&&reopened.prepared->model.lastAppliedExplanation&&proposalPlanIdentity(reopened.prepared->model.currentStateDocument.processingPlan())==proposalPlanIdentity(plan)
        &&!reopened.prepared->model.sectionProposal,"reopen exact plan and history without private proposal authority");
    const auto m=s.currentState.analysis;const auto rev=s.revision;JobControl noop;auto nochange=generateSectionProposal(s.sectionProposalRequest(),noop);
    expect(nochange.value&&(*nochange.value)->status()==SectionProposalStatus::noChange&&(*nochange.value)->passes()==0&&s.acceptSectionProposal(*nochange.value)&&s.applySectionProposal()
        &&s.revision==rev&&s.currentState.analysis==m,"all selected requirements already pass: genuine no-op preserves state and measurement");s.discardSectionProposal();
    auto exportRoot=root.getChildFile("export");exportRoot.createDirectory();auto rendered=ex::render(ex::capture(s),exportRoot,io);
    if(!rendered.value)std::cerr<<rendered.message<<'\n';expect(bool(rendered.value),"strict final PCM24 artifact rendered and independently measured");
    expect(ex::artifactEvaluation(**rendered.value,ex::capture(s)).verdict==proposal::Verdict::pass,"final reopened artifact passes all selected requirements independently of Current");
    if(channels==6)
    {const auto& peaks=verified.detail.channelTruePeaks;expect(std::max_element(peaks.begin(),peaks.end())==peaks.begin()+3,"muted LFE still supplies highest true peak and governs complete measurement");}
    auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));OutputRouting muted;
    s.selectAuditionDomain(MeasurementDomain::currentState);auto auditionRequest=s.auditionRequest();audition.synchronize(auditionRequest,muted);
    expect(audition.connect(rate).isEmpty(),"M19 injected device only, no physical output");
    JobControl copying;auto copy=prepareVerifiedSource(s.measurementContext(MeasurementDomain::source),copying);expect(bool(copy.value),"source copy for independent pre-monitor oracle");auto reader=(*copy.value)->createReader();
    auto spans=proposal::changedEnvelopeSupport(candidate->lraSelection()->before,candidate->lraSelection()->plan);
    for(int64_t first:{int64_t(0),std::max<int64_t>(0,spans.front().first-1024),std::max<int64_t>(0,spans.front().second-1024),s.sourceDocument.media.frames-4096})
    {
        first=std::clamp<int64_t>(first,0,s.sourceDocument.media.frames-4096);expect(audition.seek(first,auditionRequest,muted).isEmpty(),"exact seek into LRA transition and EOF");Capture capture(channels,4096,first);audition.setCaptureObserver(&capture);
        expect(audition.play(auditionRequest,muted,false,true).isEmpty(),"muted injected pre-monitor capture");awaitPlaying(audition,*injected);juce::AudioBuffer<float> output(maxChannels,256);
        for(int n=0;n<16;++n){injected->tick(output);audition.poll();std::this_thread::sleep_for(std::chrono::milliseconds(1));}audition.stop();audition.setCaptureObserver(nullptr);
        juce::AudioBuffer<float> raw(channels,4096);expect(reader->read(raw.getArrayOfWritePointers(),channels,first,4096),"bounded unprocessed oracle read");
        // Independently evaluate documented plateaus/linear-dB ramp endpoints,
        // then the two float32 processing stages; never call the product gainAt.
        for(int n=0;n<4096;++n)
        {
            const auto frame=first+n;double db=0;const auto& rows=plan.sectionGains->intervals;
            for(const auto& row:rows)if(frame>=row.begin&&frame<row.end){db=row.db;break;}
            for(size_t k=1;k<rows.size();++k)
            {
                const auto half=std::min({int64_t(std::llround(rate*.05)),(rows[k-1].end-rows[k-1].begin)/2,(rows[k].end-rows[k].begin)/2});const auto begin=rows[k].begin-half,end=rows[k].begin+half;
                if(frame==begin)db=rows[k-1].db;else if(frame==end-1)db=rows[k].db;
                else if(frame>begin&&frame<end-1)db=rows[k-1].db+(rows[k].db-rows[k-1].db)*(double(frame-begin)/double(end-begin-1));
            }
            for(int ch=0;ch<channels;++ch){const float local=db==0?raw.getSample(ch,n):float(double(raw.getSample(ch,n))*std::pow(10.,db/20.));raw.setSample(ch,n,float(double(local)*std::pow(10.,plan.staticGainDb/20.)));}
        }
        expect(capture.valid&&capture.next==first+4096&&sameSamples(capture.data,raw)&&audition.underruns()==0,"all-channel pre-monitor samples match independent processing oracle at seek/transitions/EOF");
    }
    audition.disconnect();
    auto global=s;expect(global.setCurrentGainDb(plan.staticGainDb-1),"explicit changed baseline for global-first LRA case");measuredComposition(global,MeasurementDomain::currentState);JobControl only;
    auto normal=generateSectionProposal(global.sectionProposalRequest(),only);
    expect(normal.value&&(*normal.value)->status()==SectionProposalStatus::ready&&(*normal.value)->passes()==1
        &&(*normal.value)->plan().sectionGains->audioKey()==plan.sectionGains->audioKey(),"LRA already meets: global-only exact assessment precedes additional local changes");
    auto row=juce::var(new juce::DynamicObject());row.getDynamicObject()->setProperty("channels",channels);row.getDynamicObject()->setProperty("seconds",(juce::Time::getMillisecondCounterHiRes()-start)/1000);
    row.getDynamicObject()->setProperty("passes",candidate->passes());row.getDynamicObject()->setProperty("before_lra",*current->detail.lra.value);row.getDynamicObject()->setProperty("after_lra",*verified.detail.lra.value);row.getDynamicObject()->setProperty("gain_db",plan.staticGainDb);cases.add(row);
    expect(root.getChildFile("proposal.txt").replaceWithText(candidate->text()),"retain full test explanation without audio in Git");
}
void failureCases(const juce::File& root)
{
    auto s=lraFixture(root,2);requireLraReduction(s);const auto request=s.sectionProposalRequest();JobControl token;
    auto limited=request;limited.policy.durationThousandths=0;auto blocked=generateSectionProposal(limited,token);
    expect(blocked.value&&(*blocked.value)->status()==SectionProposalStatus::noSolution&&(*blocked.value)->passes()==0&&(*blocked.value)->text().contains("Policy:"),"policy-limited attempt is explicit and does not measure or claim impossibility");
    limited=request;limited.policy.additionalTenths=0;limited.policy.reserveTenths=0;blocked=generateSectionProposal(limited,token);
    expect(blocked.value&&(*blocked.value)->passes()==0&&(*blocked.value)->text().contains("gain cap"),"zero additional cap blocks changing method with clear reason");
    auto missing=request;missing.sourceContext.range.begin=1;auto invalid=generateSectionProposal(missing,token);expect(invalid.value&&(*invalid.value)->status()==SectionProposalStatus::unavailable&&(*invalid.value)->passes()==0,"mismatched programme range rejected before audio work");
    missing=request;missing.currentContext.fileLayout.roles[0]=juce::AudioChannelSet::centre;invalid=generateSectionProposal(missing,token);expect(invalid.value&&(*invalid.value)->passes()==0&&(*invalid.value)->status()==SectionProposalStatus::unavailable,"changed ordered roles reject stale snapshots");
    missing=request;missing.current.reset();invalid=generateSectionProposal(missing,token);expect(invalid.value&&(*invalid.value)->passes()==0,"missing Current never becomes a candidate");
    JobControl pre;pre.cancel=true;expect(generateSectionProposal(request,pre).error==ErrorCode::cancelled,"pre-cancelled assessment publishes no result");
    SnapshotStorage disk;disk.availableBytes=[](const juce::File&){return int64_t(0);};expect(generateSectionProposal(request,token,disk).error==ErrorCode::insufficientDiskSpace,"preparation disk failure leaves previous Current valid");
    SnapshotStorage write;write.directory=root.getChildFile("failed-copy");write.directory.createDirectory();write.createWriter=[](const juce::File&)->std::unique_ptr<SnapshotWriter>{return {};};
    expect(generateSectionProposal(request,token,write).error==ErrorCode::writeFailed&&write.directory.findChildFiles(juce::File::findFiles,true).isEmpty(),"write failure cleans only owned partial copy");
    JobControl active;auto future=std::async(std::launch::async,[&]{return generateSectionProposal(request,active);});
    while(active.progress<.12f&&future.wait_for(std::chrono::milliseconds(1))!=std::future_status::ready){}
    expect(active.progress>=.12f,"cancellation reaches active full measurement");auto start=juce::Time::getMillisecondCounterHiRes();active.cancel=true;auto cancelled=future.get();
    expect(!cancelled.value&&cancelled.error==ErrorCode::cancelled&&juce::Time::getMillisecondCounterHiRes()-start<1000,"active assessment cancellation under one second without candidate adoption");
    resources.getDynamicObject()->setProperty("cancel_ms",juce::Time::getMillisecondCounterHiRes()-start);
    JobControl replacing;auto late=std::async(std::launch::async,[&]{return generateSectionProposal(request,replacing);});
    while(replacing.progress<.12f&&late.wait_for(std::chrono::milliseconds(1))!=std::future_status::ready){}
    const auto originalTime=request.currentContext.sourceFile.getLastModificationTime();expect(request.currentContext.sourceFile.setLastModificationTime(originalTime+juce::RelativeTime::seconds(2)),"change owned fixture stamp during assessment");
    const auto replaced=late.get();expect(!replaced.value&&(replaced.error==ErrorCode::sourceChanged
        ||(replaced.error==ErrorCode::invalidRequest&&replaced.message.contains("source identity changed"))),"source replacement during assessment rejects late private result");
    expect(request.currentContext.sourceFile.setLastModificationTime(originalTime),"restore owned test fixture timestamp");
    setLra(s,0);auto strict=generateSectionProposal(s.sectionProposalRequest(),token);
    expect(strict.value&&(*strict.value)->status()==SectionProposalStatus::noSolution&&(*strict.value)->passes()<=3,"strict zero LRA produces honest bounded failure rather than unset or relaxed target");
    if(strict.value)root.getChildFile("strict-zero.txt").replaceWithText((*strict.value)->text());
    for(int channels:{1,8})
    {
        auto unsupported=lraFixture(root.getChildFile("unsupported-"+juce::String(channels)),channels);requireLraReduction(unsupported);
        auto out=generateSectionProposal(unsupported.sectionProposalRequest(),token);expect(out.value&&(*out.value)->status()==SectionProposalStatus::unavailable&&(*out.value)->passes()==0,"mono/7.1 changing LRA method explicitly unsupported; existing measurements preserved");
    }
    auto shortProgramme=lraFixture(root.getChildFile("short"),2,48000,6);setLra(shortProgramme,0);
    auto shortResult=generateSectionProposal(shortProgramme.sectionProposalRequest(),token);expect(shortResult.value&&(*shortResult.value)->status()==SectionProposalStatus::unavailable&&(*shortResult.value)->passes()==0,"short unstable LRA is available descriptively but cannot drive changing proposal");
}
void stress()
{
    auto a=std::make_shared<proposal::ProgramAnalysis>();a->sampleRate=192000;a->frames=int64_t(a->sampleRate)*10800;a->identity="synthetic-three-hour-population";a->p10=-35;a->p95=-20;a->lra={15,measurement::Status::available};a->lraStable=true;
    for(int n=0;n<107985;++n){proposal::DiagnosticWindow w;w.begin=int64_t(n)*19200;w.end=w.begin+576000;w.included=true;w.lufs=(n%100<10?-20.:-35.);a->lraWindows.push_back(w);}
    const auto start=juce::Time::getMillisecondCounterHiRes();auto p=proposal::findLraPassages(a,10);const auto elapsed=juce::Time::getMillisecondCounterHiRes()-start;
#if defined(__aarch64__) || defined(__arm64__)
    const double budget=2000;
#else
    const double budget=5000;
#endif
    expect(p&&elapsed<=budget&&p->storageBytes()<64*1024*1024,"three-hour-sized synthetic population meets grouping time and storage budgets");
    resources.getDynamicObject()->setProperty("synthetic_population_ms",elapsed);resources.getDynamicObject()->setProperty("synthetic_population_bytes",juce::int64(p->storageBytes()));
}
void gatingCase(const juce::File& root)
{
    const int rate=48000,frames=rate*600;auto file=root.getChildFile("gating.wav");auto stream=file.createOutputStream();juce::WavAudioFormat format;
    std::unique_ptr<juce::AudioFormatWriter> writer(format.createWriterFor(stream.get(),rate,2,24,{},0));expect(bool(writer),"gating fixture writer");stream.release();juce::AudioBuffer<float> block(2,4096);bool ok=true;
    for(int first=0;first<frames;first+=4096)
    {
        const int count=std::min(4096,frames-first);for(int n=0;n<count;++n)
        {
            const int frame=first+n;const double seconds=double(frame)/rate;
            const double db=seconds<30?0:seconds<135?-10:seconds<195?-1:seconds<300?-10:-29.99;
            const float x=float(.1*std::pow(10.,db/20)*std::sin(2*juce::MathConstants<double>::pi*997*frame/rate));block.setSample(0,n,x);block.setSample(1,n,x);
        }
        ok&=writer->writeFromAudioSampleBuffer(block,0,count);
    }
    writer.reset();expect(ok,"bounded ten-minute synthetic gating fixture written");JobControl job;auto media=validateMedia(file,job);expect(media.succeeded(),"gating fixture validated");
    StandaloneSessionModel s;s.setValidatedSource(file,media.media);s.sourceDocument.layout=layoutFor(2);s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
    measuredComposition(s,MeasurementDomain::source);s.initializeCurrentStateFromSource();if(!s.currentState.hasMeasurement)measuredComposition(s,MeasurementDomain::currentState);s.synchronizeSections();
    setLra(s,int(std::floor((*s.currentState.analysis->detail.lra.value-3.2)*10)));auto result=generateSectionProposal(s.sectionProposalRequest(),job);
    expect(result.value.has_value(),"actual gating-change assessment completes");const auto p=*result.value;root.getChildFile("proposal.txt").replaceWithText(p->text());
    std::cout<<p->text()<<'\n';expect(p->passes()==3&&p->status()==SectionProposalStatus::noSolution,"gating-change case consumes bounded budget and reports measured failure");
    expect(p->text().contains("Earlier measured attempt retained"),"actual deeper attenuation worsens gating result; better measured attempt retained");
    expect(p->lraSelection()&&p->lraSelection()->selected.size()==2&&p->assessment()&&*p->assessment()->detail.lra.value<10,"retained complete plan preserves better earlier local envelope after failed refinement");
    expect(s.acceptSectionProposal(p)&&!s.canApplySectionProposal()&&!s.applySectionProposal(),"better but still failing measured attempt never enables Apply");
}
}
#include "LraCinemaChecks.h"
#include "LraInputChecks.h"
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI gui;juce::String rootPath,output,group,input;bool passed=false;juce::String failure;
    for(int n=1;n<argc;++n){juce::String a(argv[n]);if(a=="--root"&&n+1<argc)rootPath=argv[++n];else if(a=="--output"&&n+1<argc)output=argv[++n];else if(a=="--input"&&n+1<argc)input=argv[++n];else if(a=="--group"&&n+1<argc)group=argv[++n];else return 2;}
    try
    {
        juce::File root(rootPath);expect(rootPath.isNotEmpty()&&!root.exists()&&root.createDirectory().wasOk(),"fresh explicit M19 test directory");
        if(group=="core"){pureCases();stress();}
        else if(group=="stereo")workflow(root,2);
        else if(group=="surround")workflow(root,6);
        else if(group=="rates")for(int rate:{44100,96000,192000})workflow(root.getChildFile(juce::String(rate)),2,rate);
        else if(group=="gating")gatingCase(root);
        else if(group=="real")realProgramme(root,juce::File(input));
        else if(group=="inputs")inputCases(root);
        else if(group=="failures")failureCases(root);
        else throw std::invalid_argument("Missing/unknown explicit M19 group");passed=true;
    }
    catch(const std::exception& e){failure=e.what();std::cerr<<failure<<'\n';}
    auto result=juce::var(new juce::DynamicObject());auto* o=result.getDynamicObject();o->setProperty("suite","M19 LRA passages and measured candidates");o->setProperty("group",group);o->setProperty("passed",passed);o->setProperty("failure",failure);o->setProperty("checks",checks);o->setProperty("cases",cases);o->setProperty("resources",resources);o->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);o->setProperty("physical_devices_opened",false);o->setProperty("video_loaded",false);
    if(output.isNotEmpty())juce::File(output).replaceWithText(juce::JSON::toString(result,true));std::cout<<"M19 "<<(passed?"PASS":"FAIL")<<" "<<checks.size()<<" checks\n";return passed?0:1;
}

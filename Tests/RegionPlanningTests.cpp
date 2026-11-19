#include "CompositionTestFixtures.h"
#include "Core/Proposal/RegionPlanning.h"
#include "Core/Measurement/GatePopulation.h"
#include "Standalone/StandaloneRegionProposal.h"
#include <iostream>

namespace
{
using namespace levelscope::proposal;
using Plan=levelscope::processing::SectionGainPlan;
Plan unity(int seconds=60,int rate=48000)
{Plan p;p.frames=int64_t(seconds)*rate;p.sampleRate=rate;p.inputIdentity=juce::String::repeatedString("a",64);p.intervals={{0,p.frames,0}};return p;}
std::shared_ptr<ProgramAnalysis> synthetic(const Plan& p)
{
    auto a=std::make_shared<ProgramAnalysis>();a->identity="analysis";a->snapshotIdentity="snapshot";a->sampleRate=p.sampleRate;a->frames=p.frames;a->lraStable=true;a->p10=-30;a->p95=-20;a->lra={10,measurement::Status::available};
    for(int64_t end=p.sampleRate*3;end<=p.frames;end+=p.sampleRate/10)
    {DiagnosticWindow w;w.begin=end-p.sampleRate*3;w.end=end;w.lufs=-20;w.energy=std::pow(10.,(-20+.691)/10);w.included=true;a->lraWindows.push_back(w);}return a;
}
RegionConstraints selected(const Plan& p,RegionSpan h)
{
    auto settings=RegionPlanningSettings::defaults(RegionMode::selection,p.inputIdentity);settings.selection=h;
    auto passages=findLraPassages(synthetic(p),5);return resolveRegionConstraints(settings,{},*passages,p);
}
RegionAssignment assignment(RegionSpan support,std::vector<int64_t> cuts,std::vector<double> reductions)
{RegionAssignment a;a.geometry.support=support;a.geometry.cuts=std::move(cuts);a.reductions=std::move(reductions);return a;}
void geometryTests()
{
    auto p=unity();auto c=selected(p,{480000,1920000});ProposalPolicy policy;policy.durationThousandths=1000;
    expect(c.error.isEmpty(),"explicit selected root resolves without global permission");expect(!c.settings.allowGlobal,"selected mode never silently permits whole-programme trim");
    auto v=c.settings.encode();expect(RegionPlanningSettings::decode(v).digest()==c.settings.digest(),"independent exact settings round trip");
    auto changed=c.settings;changed.protectionMethod="rendered-freeze";expect(changed.validationError().isNotEmpty(),"unsupported mandatory rendered freeze rejected");
    const auto simple=constructRegionPlan(p,0,c,policy,{assignment(c.roots[0].hint,{}, {2})},0);
    expect(simple.error.isEmpty(),"fixed processor represents exact selected support edges");expect(simple.activity.affectedFrames==1440000&&simple.activity.splits[0]==0,"outer refinement does not spend an internal split");
    expect(simple.assignments[0].nominal.begin==482400&&simple.assignments[0].nominal.end==1917600,"support edges distinguished from plateau cuts");
    double maxDelta=0;int64_t changedSamples=0;bool exterior=true;
    for(int64_t f=0;f<p.frames;++f)
    {
        const auto gain=simple.envelope.gainAt(f);if(f<c.roots[0].hint.begin||f>=c.roots[0].hint.end)exterior&=gain==0;
        if(gain!=0)++changedSamples;maxDelta=std::max(maxDelta,-gain);
    }
    expect(exterior,"independent per-frame unchanged exterior oracle");
    expect(maxDelta==2&&changedSamples<simple.activity.affectedFrames,"conservative support charges equal ramp endpoints without changing processing");
    auto split=constructRegionPlan(p,0,c,policy,{assignment(c.roots[0].hint,{960000,1440000},{2,0,2})},0);
    expect(split.error.isEmpty()&&split.activity.splits[0]==2&&split.activity.support.size()==2,"spared middle consumes two original-root splits and two components");
    c.settings.maximumSplits=1;expect(constructRegionPlan(p,0,c,policy,split.assignments,0).error.isNotEmpty(),"child decisions cannot evade root split cap");c.settings.maximumSplits=2;
    auto unequal=constructRegionPlan(p,0,c,policy,{assignment(c.roots[0].hint,{960000,1440000},{2,1,2})},0);
    expect(unequal.error.isEmpty()&&unequal.activity.splits[0]==2&&unequal.activity.support.size()==1,"unequal adjacent parts remain one component but two splits");
    auto equal=constructRegionPlan(p,0,c,policy,{assignment(c.roots[0].hint,{960000},{2,2})},0);
    expect(equal.error.isEmpty()&&equal.identity==simple.identity&&equal.activity.splits[0]==0,"redundant equal division canonicalizes to identical audio");
    c.settings.protectedLocal={{1000000,1400000}};
    expect(constructRegionPlan(p,0,c,policy,split.assignments,0).error.isEmpty(),"unchanged middle can protect local material");
    expect(constructRegionPlan(p,0,c,policy,simple.assignments,0).error.isNotEmpty(),"cannot fill a protected gap");
    c.settings.protectedLocal.clear();c.roots[0].start={480001,480001};c.roots[0].end={1919999,1919999};
    auto pins=constructRegionPlan(p,0,c,policy,{assignment({480001,1919999},{},{2})},0);
    expect(pins.error.isEmpty()&&pins.activity.support.front().begin==480001&&pins.activity.support.back().end==1919999,"off-grid exact support pins are retained sample-exactly");
    expect(constructRegionPlan(p,0,c,policy,{assignment({480000,1919999},{},{2})},0).error.isNotEmpty(),"one-frame pin violation rejected");
    expect(constructRegionPlan(p,0,c,policy,pins.assignments,1).error.isNotEmpty(),"global permission independent of legal local geometry");
    auto limited=policy;limited.durationThousandths=200;expect(constructRegionPlan(p,0,c,limited,pins.assignments,0).error.isNotEmpty(),"programme denominator includes full ramps and respects duration cap");
    auto prior=p;prior.intervals={{0,6*48000,0},{6*48000,50*48000,-1},{50*48000,p.frames,0}};auto area=selected(prior,{10*48000,40*48000});auto overlaid=constructRegionPlan(prior,0,area,policy,{assignment(area.roots[0].hint,{}, {2})},0);
    expect(overlaid.error.isEmpty(),"new local decision preserves an existing baseline envelope");bool outside=true,nonboost=true;double additional=0;
    for(int64_t f=0;f<p.frames;++f){const auto d=prior.gainAt(f)-overlaid.envelope.gainAt(f);nonboost&=d>=0;additional=std::max(additional,d);if(f<10*48000||f>=40*48000)outside&=d==0;}
    expect(outside&&nonboost&&additional==2,"independent per-frame prior-envelope oracle: unchanged exterior, no boost and exact additional cap");
    auto crowded=p;crowded.intervals={{0,10*48000,0},{10*48000,p.frames,-1}};auto tight=selected(crowded,{10*48000,40*48000});tight.roots[0].start={10*48000,10*48000};
    expect(constructRegionPlan(crowded,0,tight,policy,{assignment(tight.roots[0].hint,{}, {2})},0).error.isNotEmpty(),"old transition at exact area edge cannot be truncated or silently expanded to fit a new ramp");
    auto shortRamp=p;shortRamp.intervals={{0,480000,0},{480000,480020,-1},{480020,p.frames,0}};
    expect(shortRamp.radius(1)==10&&shortRamp.radius(2)==10,"processor can have ten-frame transition radii; no universal fifty-ms assumption");
    auto protectedRamp=selected(shortRamp,{479990,1920000});protectedRamp.settings.protectedLocal={{479990,480030}};
    expect(constructRegionPlan(shortRamp,0,protectedRamp,policy,{assignment(protectedRamp.roots[0].hint,{}, {2})},0).error.isNotEmpty(),"short pre-existing transition protection fails closed under full old/new support accounting");
}
void controlsAndFormation()
{
    auto p=unity();auto c=selected(p,{480000,1920000});auto settings=c.settings;
    expect(RegionPlanningSettings::decode(juce::JSON::parse(juce::JSON::toString(settings.encode(),true))).digest()==settings.digest(),"settings survive actual JSON null and exact-coordinate encoding");
    for(int kind=0;kind<8;++kind)
    {
        auto bad=settings;switch(kind){case 0:bad.maximumSplits=5;break;case 1:bad.minimumIntervalMs=499;break;case 2:bad.minimumIntervalMs=10001;break;case 3:bad.separationMs=3001;break;case 4:bad.needThresholdTenths=4;break;case 5:bad.persistenceMs=299;break;case 6:bad.leftMarginMs=3050;break;case 7:bad.rightMarginMs=1;break;}
        expect(bad.validationError().isNotEmpty(),"out-of-range or off-step control rejected");
    }
    auto json=settings.encode();json.getDynamicObject()->setProperty("mandatoryFuturePermission",true);bool rejected=false;try{RegionPlanningSettings::decode(json);}catch(...){rejected=true;}expect(rejected,"unknown mandatory settings cannot widen permissions");
    auto a=synthetic(p);for(auto& w:a->lraWindows){const auto centre=(w.begin+w.end)/2;w.lufs=centre>=1200000&&centre<1440000?-30:-20;}
    auto form=formRegions(c,*a,5);auto repeated=formRegions(c,*a,5);
    expect(form.constructions<=256&&form.alternatives.size()==1&&form.alternatives[0].size()<=16&&form.nominations[0].size()<=12,"formation retains bounded roots/cuts/topologies");
    bool same=form.alternatives[0].size()==repeated.alternatives[0].size();for(size_t n=0;same&&n<form.alternatives[0].size();++n)same=form.alternatives[0][n].cuts==repeated.alternatives[0][n].cuts&&form.alternatives[0][n].support==repeated.alternatives[0][n].support;
    expect(same,"formation order is reproducible");
    expect(std::any_of(form.alternatives[0].begin(),form.alternatives[0].end(),[](const auto& g){return g.cuts.size()==2;}),"two-cut topology survives first-version geometry budget");
    expect(formRegions(c,*a,5,[]{return true;}).cancelled,"formation cancellation returns no partial authority");
    c.settings.allowSplitting=false;auto unsplit=formRegions(c,*a,5);expect(std::all_of(unsplit.alternatives[0].begin(),unsplit.alternatives[0].end(),[](const auto& g){return g.cuts.empty();}),"split-off preserves unsplit-only enumeration");
    auto passages=findLraPassages(a,5);settings.leftMarginMs=500;settings.rightMarginMs=1500;auto margin=resolveRegionConstraints(settings,{},*passages,p);
    expect(margin.roots[0].allowed.begin==456000&&margin.roots[0].allowed.end==1992000,"margins expand explicit permissions in exact source frames");
    settings.selection={0,p.frames};settings.leftMarginMs=settings.rightMarginMs=3000;margin=resolveRegionConstraints(settings,{},*passages,p);expect(margin.roots[0].allowed==RegionSpan{0,p.frames},"margins clip to programme without extending EOF");
    ProposalPolicy policy;const auto need=regionGainNeed(p,*a,{480000,1200000},policy,5);expect(need.requested==5.5&&need.headroom==3&&need.choices.back()==3,"planning reserve remains inside additional cap with disclosed uncovered need");
}
void cachedPopulationStress()
{
    auto p=unity(10800);auto a=std::make_shared<ProgramAnalysis>();a->snapshotIdentity="stress";a->sampleRate=p.sampleRate;a->frames=p.frames;a->p10=a->p95=-40.691;
    PredictionSidecar sidecar;sidecar.sampleRate=p.sampleRate;sidecar.frames=p.frames;sidecar.channels=2;sidecar.roles={measurement::Role::left,measurement::Role::right};sidecar.baselineSnapshot="stress";sidecar.baselinePlan=regionPlanIdentity(p,0);
    const int64_t step=p.sampleRate/10;sidecar.bins.reserve(108015);
    for(int64_t begin=0;begin<p.frames+p.sampleRate*3/2;begin+=step){PredictionBin b;b.begin=begin;b.end=begin+step;b.drain=begin>=p.frames;b.weightedEnergy=double(step)*1e-4;sidecar.bins.push_back(b);}
    for(int64_t end=p.sampleRate*4/10;end<=p.frames;end+=step){DiagnosticWindow w;w.begin=end-p.sampleRate*4/10;w.end=end;w.energy=1e-4;w.included=true;w.lufs=-40.691;a->momentary.push_back(w);}
    for(int64_t end=p.sampleRate*3;end<=p.frames;end+=step){DiagnosticWindow w;w.begin=end-p.sampleRate*3;w.end=end;w.energy=1e-4;w.included=true;w.lufs=-40.691;a->shortTerm.push_back(w);a->lraWindows.push_back(w);}
    const auto start=std::chrono::steady_clock::now();bool available=true;
    for(int n=0;n<64;++n)
    {
        RegionPlan plan;plan.envelope=p;plan.envelope.intervals={{0,p.sampleRate*100,0},{p.sampleRate*100,p.sampleRate*110,-.5-n*.01},{p.sampleRate*110,p.frames,0}};plan.globalDb=n*.01;
        const auto forecast=predictRegionPlan(*a,sidecar,p,0,plan);available&=forecast.integrated.has_value()&&forecast.lra.has_value();
    }
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();expect(available,"64 complete proxy plans on three-hour-sized populations (no media decode)");
#if defined(__arm64__) || defined(__aarch64__)
    expect(seconds<=10,"cached model target within 10 seconds on Sequoia arm64");
#endif
    const auto cancelStart=std::chrono::steady_clock::now();RegionPlan changed;changed.envelope=p;changed.globalDb=1;auto cancelled=predictRegionPlan(*a,sidecar,p,0,changed,[]{return true;});
    expect(cancelled.cancelled&&std::chrono::duration<double>(std::chrono::steady_clock::now()-cancelStart).count()<1,"large-population cancellation within one second");
    auto row=juce::var(new juce::DynamicObject());row.getDynamicObject()->setProperty("synthetic_population_seconds",10800);row.getDynamicObject()->setProperty("64_proxy_seconds",seconds);row.getDynamicObject()->setProperty("analysis_sidecar_bytes",int64_t(a->storageBytes()+sidecar.storageBytes()));cases.add(row);
}
void rankTests()
{
    RegionReviewScore a;a.complete=true;a.lraExcess=.299999882961;a.failedRules=2;a.otherSum=.5;a.otherMaximum=.5;a.id="a";
    auto b=a;b.lraExcess=.299999887254;b.failedRules=1;b.otherSum=b.otherMaximum=0;b.id="b";
    auto third=a;third.lraExcess=19.29;third.id="c";
    expect(selectRegionReview({a,b,third})==1,"anchored failed shortlist prefers otherwise-satisfied rules without claiming LRA pass");
    expect(selectRegionReview({a,b,third},0)==0,"zero diagnostic band retains exact best LRA excess");
    auto pass=b;pass.allPass=true;pass.failedRules=0;pass.lraExcess=0;pass.activity.affectedFrames=500;
    expect(selectRegionReview({b,pass})==1,"measured all-rule pass dominates every failure");
    auto absent=pass;absent.complete=false;expect(selectRegionReview({absent,b})==1,"missing required metrics excluded from numerical headline");
}
void modelTests(int rate,int channels)
{
    auto before=unity(4,rate);before.frames+=137;before.intervals[0].end=before.frames;auto roles=OfflineAnalysisService::measurementRoles(layoutFor(channels));
    measurement::OfflineMeasurement meter(before.sampleRate,roles);
    PredictionCollector collector(before.sampleRate,roles,before.frames,"request","snapshot",regionPlanIdentity(before,0));
    juce::AudioBuffer<float> audio(channels,4096);
    for(int64_t first=0;first<before.frames;first+=4096)
    {
        const int count=int(std::min<int64_t>(4096,before.frames-first));
        for(int ch=0;ch<channels;++ch)for(int n=0;n<count;++n)audio.setSample(ch,n,float(.02*(ch+1)*std::sin(2*juce::MathConstants<double>::pi*997*(first+n)/before.sampleRate)));
        meter.push(audio.getArrayOfReadPointers(),count);collector.push(audio.getArrayOfReadPointers(),count);
    }
    const auto measured=meter.finish();const auto sidecar=collector.finish();const auto analysis=analyzeProgramme(measured,"snapshot");
    expect(sidecar->bins.back().end==before.frames+before.sampleRate*3/2&&sidecar->bins.back().drain,"auxiliary EOF drain has explicit non-programme provenance");
    bool allPeaks=true;for(int ch=0;ch<channels;++ch){double maximum=0;for(const auto& bin:sidecar->bins)maximum=std::max(maximum,bin.sampleMax[size_t(ch)]);allPeaks&=maximum==measured.channelSamplePeaks[size_t(ch)];}expect(allPeaks,"auxiliary maxima retain every channel including excluded-from-loudness LFE");
    double maxError=0;for(const auto& w:measured.momentary)
    {
        double sum=0;for(const auto& b:sidecar->bins)if(b.begin>=w.endFrame-before.sampleRate*4/10&&b.end<=w.endFrame)sum+=b.weightedEnergy;
        maxError=std::max(maxError,std::abs(sum/(before.sampleRate*4/10)-w.energy));
    }
    expect(maxError<1e-12,"independent auxiliary K-filter/order and bins agree with qualified meter within summation test tolerance");
    RegionPlan identity;identity.envelope=before;identity.globalDb=0;
    const auto same=predictRegionPlan(*analysis,*sidecar,before,0,identity);
    expect(same.identical&&same.integrated==measured.integrated.value&&same.lra==measured.lra.value,"identity proxy reuses exact measured values without accumulated bin error");
    auto gain=identity;gain.globalDb=-6;const auto scaled=predictRegionPlan(*analysis,*sidecar,before,0,gain);
    expect(scaled.integrated&&std::abs(*scaled.integrated-(*measured.integrated.value-6))<1e-10,"ideal constant scaling forecast matches analytical IL shift");
    expect(!scaled.truePeakUpperBound&&scaled.limitations.contains("true peak unavailable"),"unknown true peak cannot masquerade as a measured constraint");
    auto cancelled=predictRegionPlan(*analysis,*sidecar,before,0,gain,[]{return true;});expect(cancelled.cancelled,"bounded model cancellation");
    auto bad=*sidecar;bad.baselineSnapshot="other";bool rejected=false;try{predictRegionPlan(*analysis,bad,before,0,gain);}catch(...){rejected=true;}expect(rejected,"mismatched prediction sidecar rejected");
}
StandaloneSessionModel programme(const juce::File& root,int channels,bool heterogeneous=false)
{
    expect(root.createDirectory().wasOk(),"owned test root");const int rate=48000,seconds=65;
    const auto file=root.getChildFile("programme.wav");expect(!file.existsAsFile(),"never overwrite original or earlier fixture");
    auto stream=file.createOutputStream();juce::WavAudioFormat format;std::unique_ptr<juce::AudioFormatWriter> writer(format.createWriterFor(stream.get(),rate,unsigned(channels),24,{},0));expect(bool(writer),"PCM24 fixture writer");stream.release();
    juce::AudioBuffer<float> block(channels,4096);bool written=true;
    for(int first=0;first<rate*seconds;first+=4096)
    {
        const int count=std::min(4096,rate*seconds-first);
        for(int ch=0;ch<channels;++ch)for(int n=0;n<count;++n)
        {const int f=first+n;const bool high=heterogeneous?((f>=rate*15&&f<rate*23)||(f>=rate*35&&f<rate*43)):(f>=rate*30&&f<rate*36);const double amplitude=channels==6&&ch==3?.12:.03*(1+.07*ch);block.setSample(ch,n,float(amplitude*(high?2.:1.)*std::sin(2*juce::MathConstants<double>::pi*997*f/rate)));}
        written&=writer->writeFromAudioSampleBuffer(block,0,count);
    }
    expect(written,"complete fixture");writer.reset();JobControl job;const auto media=validateMedia(file,job);expect(media.succeeded(),"fixture validated");StandaloneSessionModel s;expect(s.setValidatedSource(file,media.media),"fresh programme");s.sourceDocument.layout=layoutFor(channels);s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(Applicability::preProducedLinearProgramme);measuredComposition(s,MeasurementDomain::source);if(!s.currentState.hasMeasurement){s.initializeCurrentStateFromSource();measuredComposition(s,MeasurementDomain::currentState);}s.synchronizeSections();return s;
}
juce::File modifiedCore(const juce::File& from,const juce::File& to,const std::function<void(juce::var&)>& edit)
{
    juce::MemoryBlock bytes;expect(from.loadFileAsData(bytes)&&bytes.getSize()>=80,"load disposable session container");const auto* data=static_cast<const unsigned char*>(bytes.getData());uint64_t size=0;for(int n=0;n<8;++n)size|=uint64_t(data[8+n])<<(8*n);
    auto value=juce::JSON::parse(juce::String::fromUTF8(reinterpret_cast<const char*>(data+80),int(size)));edit(value);const auto text=juce::JSON::toString(value,true).toStdString();
    auto out=to.createOutputStream();expect(bool(out),"owned session mutation stream");out->write("LSSESS01",8);out->writeInt64(juce::int64(text.size()));const auto hash=juce::SHA256(text.data(),text.size()).toHexString();out->write(hash.toRawUTF8(),64);out->write(text.data(),text.size());out->write(data+80+size,bytes.getSize()-80-size);out->flush();expect(out->getStatus().wasOk(),"exact session core checksum");return to;
}
void failures(const juce::File& root)
{
    auto s=programme(root,2);auto requirement=s.jobRequirements;requirement.maximumLraTenths=49;juce::String error;s.setJobRequirements(requirement,error);
    auto settings=RegionPlanningSettings::defaults(RegionMode::automatic,sectionInputIdentity(s.measurementContext(MeasurementDomain::source)));s.setRegionPlanning(settings,error);
    const auto before=s.currentState.analysis;const auto revision=s.revision;auto request=s.sectionProposalRequest();
    JobControl early;early.cancel=true;auto cancelled=generateSectionProposal(request,early);expect(!cancelled.value&&cancelled.error==ErrorCode::cancelled,"pre-cancelled job has no partial result or Apply authority");
    JobControl live;std::atomic<bool> finished{false};std::thread monitor([&]{while(!finished.load()){if(live.progress.load()>=.15f){live.cancel=true;return;}std::this_thread::sleep_for(std::chrono::milliseconds(1));}});
    auto result=generateSectionProposal(request,live);finished=true;monitor.join();expect(!result.value&&result.error==ErrorCode::cancelled&&s.currentState.analysis==before&&s.revision==revision,"auxiliary pass cancellation preserves frozen session");
    JobControl job;SnapshotStorage lowSpace;lowSpace.directory=root;lowSpace.availableBytes=[](const juce::File&){return int64_t(0);};auto disk=generateSectionProposal(request,job,lowSpace);expect(!disk.value&&disk.error==ErrorCode::insufficientDiskSpace,"insufficient source-copy space leaves no proposal");
    auto wrong=request;wrong.currentContext.processing.staticGainDb=1;auto stale=generateSectionProposal(wrong,job);expect(stale.value&&(*stale.value)->status()==SectionProposalStatus::unavailable&&(*stale.value)->passes()==0,"mismatched measured plan rejected without a private pass");
    auto original=s.sourceDocument.file;auto moved=root.getChildFile("temporarily-missing.wav");expect(original.moveFileTo(moved),"move only owned source for replacement guard");auto changed=generateSectionProposal(request,job);const bool restored=moved.moveFileTo(original);expect(restored&&!changed.value&&changed.error==ErrorCode::sourceChanged,"missing or replaced source has no late proposal authority");
    auto unmeasurable=request;unmeasurable.regionSettings->protectionMethod="future-freeze";auto unsupported=generateSectionProposal(unmeasurable,job);expect(unsupported.value&&(*unsupported.value)->status()==SectionProposalStatus::unavailable,"unknown mandatory protection fails closed");
    auto ordinary=s;auto noLra=ordinary.jobRequirements;noLra.maximumLraTenths.reset();expect(ordinary.setJobRequirements(noLra,error),"remove optional LRA while retaining experimental settings");const auto scalar=draftGlobalTrim(ordinary.currentTrimRequest());
    expect(ordinary.acceptTrimProposal(scalar)&&ordinary.applyTrim()&&ordinary.lastAppliedExplanation&&!ordinary.lastAppliedExplanation->regionSettings,"ordinary global Apply history stays independent of inactive M20 settings");
    auto selectedSettings=RegionPlanningSettings::defaults(RegionMode::selectedSection,settings.inputIdentity);selectedSettings.sectionId=s.sectionDocument.map().sections.front().id;s.setRegionPlanning(selectedSettings,error);
    auto binding=s.sectionProposalRequest().binding();s.sectionDocument.split(selectedSettings.sectionId,48000*20,error);expect(binding!=s.sectionProposalRequest().binding(),"selected section boundary mutation stales proposal binding");
    auto saved=sf::capture(s);auto path=root.getChildFile("settings.levelscope");expect(sf::save(saved,path,job).succeeded(),"independent settings session save");
    auto old=modifiedCore(path,root.getChildFile("legacy8.levelscope"),[](juce::var& v){v.getDynamicObject()->setProperty("schema",8);v["programme"].getDynamicObject()->setProperty("schema",8);v.getDynamicObject()->removeProperty("regionPlanning");});
    auto read=sf::read(old,job);expect(read.document&&read.document->loadedSchema==8&&!read.document->regionPlanning,"schema8 migrates without automatic experimental opt-in");
    auto bad=modifiedCore(path,root.getChildFile("unknown-permission.levelscope"),[](juce::var& v){v["regionPlanning"].getDynamicObject()->setProperty("protection","future-freeze");});expect(!sf::read(bad,job).document,"unknown saved mandatory protection cannot silently become permissive");
    auto missing=modifiedCore(path,root.getChildFile("missing-selection.levelscope"),[](juce::var& v){v.getDynamicObject()->removeProperty("regionPlanning");});expect(!sf::read(missing,job).document,"schema9 requires explicit method selection field");
}

void boundedExternal(const juce::File& root,const juce::File& file,bool cinema)
{
    expect(root.createDirectory().wasOk(),"new bounded comparison root");JobControl job;auto media=validateMedia(file,job);expect(media.succeeded(),"existing bounded input verified without mutation");
    expect(media.media.sampleRate==48000&&media.media.channels==(cinema?6:2)&&media.media.frames==int64_t(48000)*(cinema?480:600),"exact existing excerpt or gate-counterexample dimensions, no full film");
    StandaloneSessionModel s;expect(s.setValidatedSource(file,media.media),"bounded input adopted");s.sourceDocument.layout=layoutFor(media.media.channels);s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(Applicability::preProducedLinearProgramme);measuredComposition(s,MeasurementDomain::source);s.synchronizeSections();
    auto requirements=s.jobRequirements;requirements.maximumLraTenths=cinema?150:67;requirements.provenance="M20 bounded development example, not a delivery specification";juce::String error;s.setJobRequirements(requirements,error);
    auto settings=RegionPlanningSettings::defaults(RegionMode::automatic,sectionInputIdentity(s.measurementContext(MeasurementDomain::source)));expect(s.setRegionPlanning(settings,error),"explicit experimental external comparison");
    const auto baseline=s.currentState.analysis;auto outcome=generateSectionProposal(s.sectionProposalRequest(),job);expect(outcome.value&&(*outcome.value)->passes()<=3,"bounded external assessment completes within pass budget");const auto& result=**outcome.value;
    expect(s.currentState.analysis==baseline&&s.sourceDocument.media.sha256==media.media.sha256,"external comparison leaves Current and original source unchanged");
    expect(result.regionReview()&&result.regionReview()->proxyEvaluations<=64,"complete proxy count remains bounded");
    std::vector<RegionReviewScore> scores;bool gateChange=false;for(const auto& attempt:result.regionReview()->attempts)
    {
        scores.push_back(regionMeasuredScore(attempt.evaluation,attempt.plan.activity,attempt.plan.identity));
        gateChange|=attempt.prediction.gateSensitive;
        if(attempt.evaluation.verdict==Verdict::pass)expect(attempt.measured[2].value&&*attempt.measured[2].value<=requirements.maximumLraTenths.value()/10.,"strict actual LRA threshold, diagnostic shortlist is no allowance");
    }
    const auto best=selectRegionReview(scores);expect(!best||result.regionReview()->retainedAttempt==*best,"headline is selected over the complete measured attempt set");
    expect(s.acceptSectionProposal(*outcome.value),"exact external result can be inspected");
    if(result.status()!=SectionProposalStatus::ready)expect(!s.canApplySectionProposal()&&!s.applySectionProposal(),"measured failure cannot Apply regardless of attractive forecast");
    root.getChildFile("proposal.txt").replaceWithText(result.text());
    if(!cinema)
    {
        // The bounded search may correctly avoid the bad deeper plan. Exercise
        // that known counterexample explicitly as a model test, not as a fourth
        // candidate assessment or an assertion about which proxies must win.
        auto analysis=analyzeProgramme(baseline->detail,baseline->identity);auto passages=findLraPassages(analysis,6.7);
        Plan before;before.inputIdentity=settings.inputIdentity;before.sampleRate=48000;before.frames=media.media.frames;before.intervals={{0,before.frames,0}};
        auto seed=selectLraTrims(*passages,before,s.proposalPolicy);auto deeper=refineLraTrims(seed,.3,s.proposalPolicy);expect(seed.available&&deeper.available,"known legal M19 deeper gate counterexample reconstructed");
        auto copy=prepareVerifiedSource(s.measurementContext(MeasurementDomain::source),job);expect(bool(copy.value),"independent model-only baseline copy");auto reader=(*copy.value)->createReader();PredictionCollector bins(48000,OfflineAnalysisService::measurementRoles(s.sourceDocument.layout),before.frames,"explicit-gate-model-test",baseline->identity,regionPlanIdentity(before,0));juce::AudioBuffer<float> block(2,4096);
        for(int64_t first=0;first<before.frames;first+=4096){int count=int(std::min<int64_t>(4096,before.frames-first));if(!reader->read(block.getArrayOfWritePointers(),2,first,count))throw std::runtime_error("model counterexample read");bins.push(block.getArrayOfReadPointers(),count);}auto sidecar=bins.finish();
        RegionPlan first;first.envelope=seed.plan;RegionPlan changed;changed.envelope=deeper.plan;
        const auto predictedSeed=predictRegionPlan(*analysis,*sidecar,before,0,first),predictedDeep=predictRegionPlan(*analysis,*sidecar,before,0,changed);
        expect(predictedDeep.lraEntered>0&&predictedSeed.lra&&predictedDeep.lra&&*predictedDeep.lra>*predictedSeed.lra,"full-population proxy reveals deeper attenuation admitting excluded windows and worsening LRA");
        auto record=juce::var(new juce::DynamicObject());record.getDynamicObject()->setProperty("model_only_counterexample",true);record.getDynamicObject()->setProperty("seed_predicted_lra",*predictedSeed.lra);record.getDynamicObject()->setProperty("deeper_predicted_lra",*predictedDeep.lra);record.getDynamicObject()->setProperty("newly_included_windows",int(predictedDeep.lraEntered));cases.add(record);
    }
    auto row=juce::var(new juce::DynamicObject());auto* o=row.getDynamicObject();o->setProperty("content_sha256",media.media.sha256);o->setProperty("frames",media.media.frames);o->setProperty("before_lra",*baseline->detail.lra.value);o->setProperty("before_integrated",*baseline->detail.integrated.value);o->setProperty("target_lra",requirements.maximumLraTenths.value()/10.);o->setProperty("outcome",result.outcomeTitle());o->setProperty("passes",result.passes());o->setProperty("proxies",int(result.regionReview()->proxyEvaluations));o->setProperty("gate_change_predicted",gateChange);o->setProperty("planning_seconds",result.regionReview()->planningSeconds);o->setProperty("planning_bytes",int64_t(result.regionReview()->planningBytes));
    if(result.assessment()){o->setProperty("retained_lra",*result.assessment()->detail.lra.value);o->setProperty("retained_integrated",*result.assessment()->detail.integrated.value);}cases.add(row);
}

#include "RegionAudioChecks.h"
void workflow(const juce::File& root,int channels,bool heterogeneous)
{
    auto s=programme(root,channels,heterogeneous);auto requirements=s.jobRequirements;requirements.maximumLraTenths=49;requirements.provenance="M20 synthetic evaluation requirement";juce::String error;expect(s.setJobRequirements(requirements,error),"LRA requirement without remeasurement");
    auto settings=RegionPlanningSettings::defaults(heterogeneous?RegionMode::selection:RegionMode::automatic,sectionInputIdentity(s.measurementContext(MeasurementDomain::source)));
    if(heterogeneous){settings.selection={int64_t(48000)*12,int64_t(48000)*46};settings.allowGlobal=true;auto policy=s.proposalPolicy;policy.durationThousandths=600;expect(s.setProposalPolicy(policy,error),"explicit heterogeneous test duration preference");}
    expect(s.setRegionPlanning(settings,error),"explicit M20 opt-in preserves measurements");const auto frozen=s;
    const auto source=s.source.analysis->identity,current=s.currentState.analysis->identity;const auto revision=s.revision;
    {
        auto analysis=analyzeProgramme(s.currentState.analysis->detail,current);auto passages=findLraPassages(analysis,4.9);Plan before;before.inputIdentity=settings.inputIdentity;before.sampleRate=48000;before.frames=s.sourceDocument.media.frames;before.intervals={{0,before.frames,0}};
        auto constraints=resolveRegionConstraints(settings,s.sectionDocument.map(),*passages,before);auto formation=formRegions(constraints,*analysis,4.9);juce::String text;
        for(size_t r=0;r<formation.alternatives.size();++r){for(auto cut:formation.nominations[r])text+="cut "+juce::String(cut.frame/48000.,6)+" strength "+juce::String(cut.strength,16)+" "+cut.reason+"\n";for(auto g:formation.alternatives[r]){text+="geometry "+juce::String(g.support.begin/48000.,6)+":"+juce::String(g.support.end/48000.,6)+" "+g.reason;for(auto f:g.cuts)text+=" "+juce::String(f/48000.,6);text+="\n";}}
        root.getChildFile("formation.txt").replaceWithText(text);
    }
    JobControl control;const auto result=generateSectionProposal(s.sectionProposalRequest(),control);expect(bool(result.value),"region service completed safely");const auto proposed=*result.value;
    root.getChildFile("proposal.txt").replaceWithText(proposed->text());
    expect(proposed->methodVersion()==regionSearchVersion&&proposed->passes()<=3,"M20 selected explicitly and bounded assessment count");
    expect(proposed->regionReview()&&proposed->regionReview()->proxyEvaluations<=64,"bounded proxy accounting exposed");
    expect(s.source.analysis->identity==source&&s.currentState.analysis->identity==current&&s.revision==revision,"planning does not alter original/current audio or measurements");
    expect(s.acceptSectionProposal(proposed),"immutable request completion accepted");
    auto stale=s;expect(stale.setCurrentGainDb(.1)&&!stale.canApplySectionProposal(),"gain edit rejects a stale region proposal");
    stale=s;auto edited=settings;edited.maximumSplits=1;expect(stale.setRegionPlanning(edited,error)&&!stale.canApplySectionProposal()&&stale.currentState.analysis->identity==current,"changed permissions stale proposals but reuse compatible measurement");
    stale=s;stale.sessionId="reopened-runtime";expect(!stale.acceptSectionProposal(proposed),"late completion rejected by new runtime identity");
    stale=s;stale.auditionBusy=[] {return true;};expect(!stale.applySectionProposal(),"Apply remains stopped-only");
    if(!heterogeneous)
    {
        auto dense=s;for(int n=1;n<10;++n)expect(dense.sectionDocument.split(n,n*48000,error),"redundant annotation split");
        expect(dense.sectionProposalRequest().binding()==s.sectionProposalRequest().binding(),"unselected markers and names do not alter automatic region binding");
    }
    s.setOutputRoute(channels==6?3:0,-1);expect(s.sectionProposalRequest().binding()==proposed->identity(),"partial monitor muting preserves all-channel region identity");
    expect(proposed->status()==SectionProposalStatus::ready,"representative fixture has a complete measured passing candidate");
    if(heterogeneous)
    {
        const auto& chosen=proposed->regionReview()->attempts[proposed->regionReview()->retainedAttempt].plan;
        expect(chosen.activity.support.size()==2&&chosen.activity.splits[0]==2,"passing broad-root plan spares the middle with two original-root decisions");
        auto unsplit=frozen;auto withoutSplits=settings;withoutSplits.allowSplitting=false;expect(unsplit.setRegionPlanning(withoutSplits,error),"explicit split-off comparison under identical requirements");JobControl compare;auto whole=generateSectionProposal(unsplit.sectionProposalRequest(),compare);expect(whole.value&&(*whole.value)->status()==SectionProposalStatus::noSolution,"split geometry finds a measured solution missed by bounded whole-root search");root.getChildFile("unsplit-comparison.txt").replaceWithText((*whole.value)->text());
    }
    expect(s.applySectionProposal()&&s.currentStateNeedsRemeasurement(),"explicit Apply invalidates Current; private result never substitutes");
    const auto actual=measuredComposition(s,MeasurementDomain::currentState);expect(s.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass,"independent Current passes selected requirements");
    expect(actual.detail.integrated.value==proposed->assessment()->detail.integrated.value&&actual.detail.lra.value==proposed->assessment()->detail.lra.value
        &&actual.detail.channelTruePeaks==proposed->assessment()->detail.channelTruePeaks&&actual.detail.truePeakEighthFrames==proposed->assessment()->detail.truePeakEighthFrames,"independent Current exactly matches complete candidate and every peak position");
    expect(s.source.analysis->identity==source,"Source preserved after Apply and remeasurement");
    expect(s.lastAppliedExplanation&&s.lastAppliedExplanation->regionSettings.has_value(),"historical explanation retains complete region controls");
    const auto bytes=encodeAppliedExplanation(*s.lastAppliedExplanation);const auto history=decodeAppliedExplanation(bytes,control);expect(history->regionSettings->digest()==settings.digest(),"versioned historical region definition round trip");
    auto longHistory=*s.lastAppliedExplanation;longHistory.narrative=juce::String::repeatedString("Complete bounded historical narrative.\n",3000);
    expect(longHistory.narrative.getNumBytesAsUTF8()>65536&&decodeAppliedExplanation(encodeAppliedExplanation(longHistory),control)->narrative==longHistory.narrative,"historical narrative beyond generic 64 KiB string cap retains exact content within total attachment budget");
    s.discardSectionProposal();const auto doc=sf::capture(s);const auto path=root.getChildFile("session.levelscope");expect(sf::save(doc,path,control).succeeded(),"schema9 coherent session saved");const auto reopened=sf::read(path,control);expect(reopened.succeeded()&&reopened.document->loadedSchema==9&&reopened.document->regionPlanning->digest()==settings.digest(),"schema9 settings recovered independently of caches");
    auto prepared=sf::verify(*reopened.document,control);expect(prepared.prepared&&prepared.prepared->model.regionPlanning->digest()==settings.digest()&&!prepared.prepared->model.sectionProposal,"session reopens settings and history without private Apply authority");
    const auto exportRoot=root.getChildFile("export");exportRoot.createDirectory();const auto rendered=ex::render(ex::capture(s),exportRoot,control);
    expect(rendered.value&&ex::artifactEvaluation(**rendered.value,ex::capture(s)).verdict==Verdict::pass,"finalized PCM24 independently measured against every selected requirement");
    expect((*rendered.value)->analysis().context.processing.staticGainDb==0,"artifact verification never reapplies processing gain");
    const auto noChange=generateSectionProposal(s.sectionProposalRequest(),control);expect(noChange.value&&(*noChange.value)->status()==SectionProposalStatus::noChange&&(*noChange.value)->passes()==0,"already passing M20 preserves valid measurements without another pass");
    auto row=juce::var(new juce::DynamicObject());auto* o=row.getDynamicObject();o->setProperty("channels",channels);o->setProperty("heterogeneous",heterogeneous);o->setProperty("before_lra",*s.source.analysis->detail.lra.value);o->setProperty("after_lra",*actual.detail.lra.value);o->setProperty("proxy_count",int(proposed->regionReview()->proxyEvaluations));o->setProperty("passes",proposed->passes());o->setProperty("planning_seconds",proposed->regionReview()->planningSeconds);o->setProperty("planning_bytes",int64_t(proposed->regionReview()->planningBytes));cases.add(row);
}
}
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI gui;juce::String output,group,root,input;bool passed=false;juce::String failure;
    for(int n=1;n<argc;++n){juce::String arg(argv[n]);if(arg=="--group"&&n+1<argc)group=argv[++n];else if(arg=="--output"&&n+1<argc)output=argv[++n];else if(arg=="--root"&&n+1<argc)root=argv[++n];else if(arg=="--input"&&n+1<argc)input=argv[++n];else return 2;}
    try{if(group=="core"){geometryTests();controlsAndFormation();rankTests();for(int rate:{44100,48000,96000,192000})for(int channels:{2,6})modelTests(rate,channels);}else if(group=="cinema"||group=="gating")boundedExternal(juce::File(root),juce::File(input),group=="cinema");else if(group=="audio")audioCases(juce::File(root));else if(group=="failures")failures(juce::File(root));else if(group=="population")cachedPopulationStress();else if(group=="stereo")workflow(juce::File(root),2,false);else if(group=="surround")workflow(juce::File(root),6,false);else if(group=="heterogeneous")workflow(juce::File(root),2,true);else throw std::runtime_error("Required explicit M20 group missing/unknown");passed=true;}
    catch(const std::exception& e){failure=e.what();}
    auto result=juce::var(new juce::DynamicObject());auto* o=result.getDynamicObject();o->setProperty("suite","M20 Logic development checks");o->setProperty("group",group);o->setProperty("passed",passed);o->setProperty("failure",failure);o->setProperty("checks",checks);o->setProperty("cases",cases);o->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);o->setProperty("devices_opened",false);o->setProperty("video_loaded",false);
    if(output.isNotEmpty())juce::File(output).replaceWithText(juce::JSON::toString(result,true));std::cout<<(passed?"PASS ":"FAIL ")<<group<<" "<<failure<<std::endl;return passed?0:1;
}

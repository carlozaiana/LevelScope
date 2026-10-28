#include "CompositionTestFixtures.h"
#include "Standalone/StandaloneSectionProposal.h"
#include <algorithm>
namespace
{
double oracleGain(const processing::SectionGainPlan& plan,int64_t frame)
{
    std::vector<processing::GainInterval> rows;
    for(const auto& i:plan.intervals)
        if(!rows.empty()&&rows.back().db==i.db)rows.back().end=i.end;else rows.push_back(i);
    for(size_t n=1;n<rows.size();++n)
    {
        const auto& a=rows[n-1];const auto& b=rows[n];
        const int64_t half=std::min({int64_t(plan.sampleRate/20),(a.end-a.begin)/2,(b.end-b.begin)/2});
        const auto start=b.begin-half,end=b.begin+half;
        if(frame==start)return a.db;
        if(frame==end-1)return b.db;
        if(frame>start&&frame<end-1)return a.db+(b.db-a.db)*(double(frame-start)/double(end-start-1));
    }
    for(const auto& row:rows)if(frame>=row.begin&&frame<row.end)return row.db;
    throw std::runtime_error("oracle frame outside programme");
}
juce::AudioBuffer<float> oracleAudio(const juce::AudioBuffer<float>& raw,const ProcessingPlan& plan,int64_t first=0)
{
    auto out=raw;const double global=std::pow(10.0,plan.staticGainDb/20.0);
    for(int n=0;n<out.getNumSamples();++n)
    {
        const auto db=plan.sectionGains?oracleGain(*plan.sectionGains,first+n):0;
        const auto coefficient=std::pow(10.0,db/20.0);
        for(int c=0;c<out.getNumChannels();++c)
        {
            const float section=db==0?raw.getSample(c,n):static_cast<float>(double(raw.getSample(c,n))*coefficient);
            out.setSample(c,n,static_cast<float>(double(section)*global));
        }
    }
    return out;
}

proposal::SectionMap mapFor(int rate,int64_t frames,int count)
{
    proposal::SectionMap m;m.inputIdentity=juce::String::repeatedString("a",64);m.sampleRate=rate;m.frames=frames;m.nextId=uint64_t(count+1);
    for(int n=0;n<count;++n)m.sections.push_back({uint64_t(n+1),frames*n/count,frames*(n+1)/count,"Section "+juce::String(n+1)});
    return m;
}
proposal::TrimRequest scalar(double loudness,double peak)
{
    proposal::TrimRequest r;r.profile=proposal::ebuR128ProgrammeQc();r.sourceInputIdentity=r.currentInputIdentity="same-programme";
    r.current.domain=proposal::SnapshotDomain::currentState;r.current.currentIdentityMatches=true;r.current.snapshotIdentity="current";
    r.current.measurementVersion=measurement::version;r.current.applicability=proposal::Applicability::preProducedLinearProgramme;
    r.current.integrated={loudness,measurement::Status::available};r.current.truePeak={peak,measurement::Status::available};
    r.source=r.current;r.source.domain=proposal::SnapshotDomain::source;r.source.snapshotIdentity="source";return r;
}
proposal::SectionPeakMap peakMap(const proposal::SectionMap& m,double peak,int changed=1)
{
    proposal::SectionPeakMap p;p.geometry=proposal::sectionGeometryKey(m);p.snapshotIdentity="current";p.channels=6;
    p.sections.resize(m.sections.size());
    for(size_t n=0;n<p.sections.size();++n)
    {
        p.sections[n].sample.resize(6);p.sections[n].truePeak.resize(6);
        if(n<size_t(changed))p.sections[n].truePeak[3]={std::pow(10.,peak/20.),m.sections[n].begin*8,false};
    }
    return p;
}
void policyCases()
{
    auto m=mapFor(48000,480000,100);auto r=scalar(-24,-.5);auto peaks=peakMap(m,-.5);
    auto s=proposal::selectSectionTrims(r,m,{},peaks);
    expect(s.available&&s.corrections.size()==1&&s.corrections[0].channel==3,"all-channel LFE peak selects one bounded section");
    expect(std::abs(s.corrections[0].reductionDb-2)<1e-12&&s.corrections[0].reserveDb==.5,"reserve included inside additional attenuation");
    r=scalar(-24,.8);peaks=peakMap(m,.8);s=proposal::selectSectionTrims(r,m,{},peaks);
    expect(s.available&&s.corrections[0].reductionDb==3&&std::abs(s.corrections[0].reserveDb-.2)<1e-12&&s.corrections[0].uncoveredDb==0,"partial reserve disclosed at three dB cap");
    r=scalar(-24,1.2);peaks=peakMap(m,1.2);s=proposal::selectSectionTrims(r,m,{},peaks);
    expect(s.available&&s.corrections[0].reductionDb==3&&s.corrections[0].reserveDb==0&&std::abs(s.corrections[0].uncoveredDb-.2)<1e-12,"uncovered deficit disclosed without widening ceiling");
    peaks=peakMap(m,1.2,8);expect(proposal::selectSectionTrims(r,m,{},peaks).available,"eight changed sections accepted within extent budget");
    peaks=peakMap(m,1.2,9);s=proposal::selectSectionTrims(r,m,{},peaks);expect(!s.available&&s.corrections.size()==9,"ninth section rejects full set without truncation");
    auto edge=mapFor(48000,10005,2);edge.sections[0].end=1334;edge.sections[1].begin=1334;
    expect(proposal::selectSectionTrims(r,edge,{},peakMap(edge,1.2)).available,"exact twenty percent including transition support accepted");
    --edge.frames;--edge.sections.back().end;
    expect(!proposal::selectSectionTrims(r,edge,{},peakMap(edge,1.2)).available,"one frame beyond twenty percent rejects without epsilon");
    processing::SectionGainPlan fractional;fractional.inputIdentity=m.inputIdentity;fractional.sampleRate=m.sampleRate;fractional.frames=m.frames;
    for(double gain:{-1.23456789012345,-.1,.1,3.12345678901234})
    {
        fractional.intervals={{0,m.frames,gain}};const auto precise=proposal::selectSectionTrims(r,m,fractional,peakMap(m,1.2));
        expect(precise.available&&gain-precise.plan.intervals.front().db<=3&&precise.plan.intervals.front().db<gain,"representable additional attenuation never exceeds exact three dB cap");
    }
    auto broad=mapFor(48000,480000,4);expect(!proposal::selectSectionTrims(r,broad,{},peakMap(broad,1.2)).available,"more than twenty percent local extent rejected");
    processing::SectionGainPlan before;before.inputIdentity=m.inputIdentity;before.sampleRate=48000;before.frames=10000;before.intervals={{0,1000,0},{1000,10000,0}};
    auto after=before;after.intervals[0].db=-1;
    expect(proposal::changedEnvelopeExtent(before,after)==1500,"transition support counted in local extent at EOF-side plateau");
    auto manual=s.plan;manual.inputIdentity=m.inputIdentity;manual.frames=m.frames;manual.sampleRate=m.sampleRate;manual.intervals={{0,m.frames,-11}};
    expect(!proposal::selectSectionTrims(r,m,manual,peakMap(m,1.2)).available,"existing manual cumulative minus twelve dB bound protected");
    peaks=peakMap(m,1.2);peaks.geometry+="stale";expect(!proposal::selectSectionTrims(r,m,{},peaks).available,"stale peak geometry rejected");
    peaks=peakMap(m,1.2);peaks.sections[0].truePeak[0].magnitude=std::numeric_limits<double>::quiet_NaN();expect(!proposal::selectSectionTrims(r,m,{},peaks).available,"nonfinite attribution rejected");
    auto many=mapFor(48000,480000,10000);auto manyPeaks=peakMap(many,1.2,10000);auto manySeed=proposal::selectSectionTrims(r,many,{},manyPeaks);
    const auto manyText=proposal::sectionSeedText(manySeed);
    expect(!manySeed.available&&manySeed.corrections.size()==10000&&manyText.contains("complete set was rejected"),"large over-budget failure retains full data and bounds explanation text without truncating a plan");
    expect(manyPeaks.storageBytes()+manySeed.corrections.capacity()*sizeof(proposal::SectionCorrection)
        +manySeed.plan.intervals.capacity()*sizeof(processing::GainInterval)+manyText.getNumBytesAsUTF8()<16*1024*1024,"combined retained attribution proposal and explanation under sixteen MiB at ten thousand sections");
    r=scalar(-23,-1);expect(proposal::draftGlobalTrim(r).status==proposal::TrimStatus::noChange,"exact passing thresholds no change");
    r.current.truePeak.value=std::nextafter(-1.,INFINITY);expect(proposal::draftGlobalTrim(r).status!=proposal::TrimStatus::noChange,"immediately outside peak ceiling is not pass");
    r.current.integrated.value=-22.8;r.current.truePeak.value=-2;expect(proposal::draftGlobalTrim(r).status==proposal::TrimStatus::noChange,"exact loudness upper boundary no change");
    r.current.integrated.value=std::nextafter(-22.8,INFINITY);expect(proposal::draftGlobalTrim(r).status!=proposal::TrimStatus::noChange,"immediately outside loudness allowance is not pass");
    r=scalar(-23,-2);r.current.domain=proposal::SnapshotDomain::privateCandidate;expect(proposal::draftGlobalTrim(r).status==proposal::TrimStatus::unavailable,"candidate cannot enter regular M3 Apply authority");
}
struct Event {size_t channel;double value;int64_t position;bool padding;};
struct Observer final:measurement::PeakObserver
{
    proposal::SectionPeakCollector collector;std::vector<Event> samples,peaks;
    Observer(const proposal::SectionMap& m,int ch):collector(m,ch,"observed"){}
    void sample(size_t c,double v,int64_t f)override{samples.push_back({c,v,f*8,false});collector.sample(c,v,f);}
    void truePeak(size_t c,double v,int64_t f,bool p)override{peaks.push_back({c,v,f,p});collector.truePeak(c,v,f,p);}
};
void collectorCases()
{
    for(int rate:{44100,48000,96000,192000})
    {
        const int frames=rate/2;auto m=mapFor(rate,frames,20);Observer observer(m,6);
        const auto roles=OfflineAnalysisService::measurementRoles(layoutFor(6));
        measurement::OfflineMeasurement plain(rate,roles),observed(rate,roles,&observer);
        juce::AudioBuffer<float> audio(6,frames);
        for(int c=0;c<6;++c)for(int n=0;n<frames;++n)audio.setSample(c,n,float(.1*std::sin(n*.237+c)));
        for(int c=0;c<6;++c){audio.setSample(c,0,.5f);audio.setSample(c,frames/2,1.f);audio.setSample(c,frames-1,-1.f);}
        plain.push(audio.getArrayOfReadPointers(),frames);observed.push(audio.getArrayOfReadPointers(),frames);
        const auto a=plain.finish(),b=observed.finish();auto p=observer.collector.finish();
        const auto sameWindows=[](const auto& x,const auto& y)
        {return x.size()==y.size()&&std::equal(x.begin(),x.end(),y.begin(),[](const auto& u,const auto& v){return u.endFrame==v.endFrame&&u.energy==v.energy;});};
        expect(a.momentaryMax.value==b.momentaryMax.value&&a.shortTermMax.value==b.shortTermMax.value&&a.integratedRelativeGate==b.integratedRelativeGate
            &&a.lraRelativeGate==b.lraRelativeGate&&a.p10==b.p10&&a.p95==b.p95&&a.integratedIncluded==b.integratedIncluded&&a.lraIncluded==b.lraIncluded
            &&a.lraStable==b.lraStable&&sameWindows(a.momentary,b.momentary)&&sameWindows(a.shortTerm,b.shortTerm)&&sameWindows(a.lraWindows,b.lraWindows),"observer preserves every full-resolution loudness population gate and stability result");
        expect(a.integrated.value==b.integrated.value&&a.truePeak.value==b.truePeak.value&&a.lra.value==b.lra.value
            &&a.channelTruePeaks==b.channelTruePeaks&&a.truePeakEighthFrames==b.truePeakEighthFrames&&a.samplePeakFrames==b.samplePeakFrames,"observer leaves ordinary metrics and extrema bit-exact");
        bool equal=true;
        for(size_t n=0;n<m.sections.size();++n)for(size_t c=0;c<6;++c)
        {
            proposal::LocatedPeak sample,peak;const auto& section=m.sections[n];
            const int64_t margin=int64_t(std::llround(rate*.05))+64;
            for(const auto& e:observer.samples)if(e.channel==c&&e.position>=section.begin*8&&e.position<section.end*8
                &&e.value>sample.magnitude)sample={e.value,e.position,false};
            for(const auto& e:observer.peaks)
            {
                const auto loc=std::clamp<int64_t>(e.position,0,frames*8-1);
                if(e.channel==c&&loc>=std::max<int64_t>(0,section.begin-margin)*8&&loc<std::min<int64_t>(frames,section.end+margin)*8
                    &&(e.value>peak.magnitude||(e.value==peak.magnitude&&e.position<peak.eighthFrame)))
                    peak={e.value,e.position,e.padding||e.position<0||e.position>=frames*8};
            }
            equal&=sample.magnitude==p->sections[n].sample[c].magnitude&&sample.eighthFrame==p->sections[n].sample[c].eighthFrame
                &&peak.magnitude==p->sections[n].truePeak[c].magnitude&&peak.eighthFrame==p->sections[n].truePeak[c].eighthFrame
                &&peak.padding==p->sections[n].truePeak[c].padding;
        }
        expect(equal,"all section/channel maxima match independent full-event oracle including boundaries and EOF drain");
    }
    const auto m=mapFor(48000,480000,10000);proposal::SectionPeakCollector collector(m,12,"large");
    collector.truePeak(11,2,-3,true);collector.truePeak(11,3,m.frames*8+1,true);
    const auto p=collector.finish();expect(p->storageBytes()<16*1024*1024,"10000-section twelve-channel retained attribution under sixteen MiB");
    expect(p->sections.front().truePeak[11].eighthFrame==-3&&p->sections.back().truePeak[11].eighthFrame==m.frames*8+1,"padding coordinates retained without extending programme duration");
}
StandaloneSessionModel proposalFixture(const juce::File& root,int rate,int channels,int mode,float crest=.8f)
{
    root.createDirectory();const int frames=rate*10;juce::AudioBuffer<float> audio(channels,frames);
    // Pin decoded PCM across libm implementations; tiny sine residuals must not
    // create different source identities on the two actual qualification Macs.
    for(int c=0;c<channels;++c)for(int n=0;n<frames;++n)
        audio.setSample(c,n,c==3?0.f:float(std::llround(32768*(channels==2?.06:.034)*std::sin(2*3.141592653589793*1000*n/rate)))/32768.f);
    audio.setSample(channels==6?3:0,rate*5+rate/4,crest);
    StandaloneSessionModel s;
    if(mode<2)
    {
        auto p=component(root,"programme",rate,audio,mode==1);
        if(mode==1)expect(s.setValidatedProgramme(*p.assembly),"adopt grouped fixture");
        else{expect(s.setValidatedSource(p.file,p.media),"adopt fixture");s.sourceDocument.layout=p.layout;}
    }
    else
    {
        auto me=audio;const int dialogueChannels=mode>=4?1:channels;
        juce::AudioBuffer<float> dialogue(dialogueChannels,frames);
        for(int c=0;c<dialogueChannels;++c)for(int n=0;n<frames;++n)
        {
            const int dest=mode>=4?2:c;const float half=audio.getSample(dest,n)*.5f;
            dialogue.setSample(c,n,half);me.setSample(dest,n,half);
        }
        ProgrammeComposition composition;composition.preset=channels==2?CompositionPreset::stereoStereo
            :mode>=4?CompositionPreset::surroundMonoCentre:CompositionPreset::surroundSurround;
        composition.me=component(root,"me",rate,me,mode%2==1);
        composition.dialogue=component(root,"dialogue",rate,dialogue,mode%2==1);
        expect(composition.confirmComplementaryBalancedAlignment()&&s.setValidatedComposition(composition),"adopt aligned unity composition fixture");
    }
    expect(s.initializeCurrentStateFromSource(),"initialize fixture Current");s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
    s.synchronizeSections();juce::String error;for(int n=1;n<20;++n)expect(s.sectionDocument.split(uint64_t(n),n*rate/2,error),"explicit half-second manual annotations");
    Capture raw(channels,frames);measuredComposition(s,MeasurementDomain::source,&raw);
    expect(root.getChildFile("independent-input.f32").replaceWithData(raw.data.getReadPointer(0),0),"create owned parity capture");
    {auto stream=root.getChildFile("independent-input.f32").createOutputStream();for(int ch=0;ch<channels;++ch)expect(stream->write(raw.data.getReadPointer(ch),size_t(frames)*4),"retain source samples for independent export quantizer oracle");}
    measuredComposition(s,MeasurementDomain::currentState);return s;
}
void workflow(const juce::File& root,int rate,int channels,int mode)
{
    auto s=proposalFixture(root,rate,channels,mode);
    if(rate==48000&&channels==6&&mode==5)
    {
        expect(s.beginSectionTrimDraft()&&s.setSectionDraftGain(1,-.1)&&s.applySectionTrimDraft(),"existing explicit manual envelope before automatic proposal");
        measuredComposition(s,MeasurementDomain::currentState);
    }
    const auto source=s.source.analysis,current=s.currentState.analysis;
    const auto revision=s.revision;auto request=s.sectionProposalRequest();
    if(rate==48000&&mode==0){JobControl save;expect(sf::save(sf::capture(s),root.getChildFile("baseline.levelscope"),save).succeeded(),"retain measured baseline UI fixture");}
    JobControl control;auto result=generateSectionProposal(request,control);
    if(!result.value)std::cerr<<result.message<<'\n';else std::cout<<(*result.value)->text()<<'\n';
    expect(result.value&&(*result.value)->status()==SectionProposalStatus::ready,"bounded high-crest fixture yields exactly assessed ready candidate");
    auto candidate=*result.value;
    if(rate==48000&&channels==6&&mode==5)
        expect(candidate->selection().plan.intervals[1].db==-.1,"automatic local proposal preserves unrelated pre-existing manual gain exactly");
    expect(candidate->passes()==3&&candidate->assessment()->domain==MeasurementDomain::privateCandidate,"one envelope uses exactly three private passes");
    expect(s.acceptSectionProposal(candidate)&&s.revision==revision&&s.currentState.analysis==current&&s.source.analysis==source,"generation changes no applied audio or measured authority");
    auto cached=s.sectionProposalRequest();JobControl repeatControl;auto repeat=generateSectionProposal(cached,repeatControl);
    expect(repeat.value&&(*repeat.value)->passes()==2&&(*repeat.value)->plan().sectionGains->audioKey()==candidate->plan().sectionGains->audioKey()
        &&(*repeat.value)->plan().staticGainDb==candidate->plan().staticGainDb,"compatible attribution reuse skips baseline and repeats identical plan");
    s.setOutputRoute(channels==6?3:0,-1);s.referenceVideo.startOnAudio=video::Time(-1001,30000);
    expect(s.canApplySectionProposal(),"muted peak channel and rational video offset do not stale audio proposal");
    auto changed=s;changed.setCurrentGainDb(.125);expect(!changed.acceptSectionProposal(candidate),"processing mutation rejects late proposal");
    changed=s;changed.setProgrammeApplicability(proposal::Applicability::unconfirmed);expect(!changed.acceptSectionProposal(candidate),"profile applicability invalidates proposal");
    changed=s;juce::String error;expect(changed.sectionDocument.moveEnd(1,rate/2+1,error)&&!changed.applySectionProposal(),"geometry edit rejects Apply without remapping manual gains");
    changed=s;expect(changed.sectionDocument.rename(1,"cosmetic",error)&&changed.canApplySectionProposal(),"cosmetic annotation names preserve proposal/audio identity");
    changed=s;changed.sessionId="new-runtime-after-reopen";expect(!changed.acceptSectionProposal(candidate),"late private completion cannot enter another runtime session with matching audio");
    changed=s;changed.auditionBusy=[] { return true; };
    expect(!changed.canApplySectionProposal()&&!changed.applySectionProposal()&&changed.currentState.analysis==current&&changed.revision==revision,"active audition blocks proposal Apply without changing audio or measured state");
    bool saveRefused=false;try{sf::capture(s);}catch(...){saveRefused=true;}expect(saveRefused,"transient private candidate cannot persist as Current authority");
    expect(s.applySectionProposal()&&s.revision==revision+1&&!s.currentState.hasMeasurement&&s.currentStateNeedsRemeasurement()&&s.source.analysis==source,"explicit Apply atomically changes plan once and leaves Current unverified");
    expect(s.job.message=="Proposal applied explicitly. Current State UNVERIFIED; re-measure independently.","Apply replaces obsolete private-completion status with current unverified status");
    const auto plan=s.currentStateDocument.processingPlan();s.measurementUnsuccessful(MeasurementDomain::currentState,ErrorCode::cancelled);
    expect(s.currentStateNeedsRemeasurement()&&s.currentStateDocument.processingPlan().sectionGains->audioKey()==plan.sectionGains->audioKey(),"unsuccessful verification retains applied plan and unverified state");
    const auto verification=measuredComposition(s,MeasurementDomain::currentState);
    expect(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass
        &&verification.measured.integratedLufs==candidate->assessment()->measured.integratedLufs
        &&verification.measured.truePeakDbtp==candidate->assessment()->measured.truePeakDbtp,"independent Current verification equals exact private candidate and passes selected rules");
    if(rate==48000&&mode==0)
    {
        auto faulted=s;expect(faulted.job.begin(JobStage::analysis,faulted.revision),"start independent verification fault injection");
        auto unexpected=verification;unexpected.context=faulted.measurementContext(MeasurementDomain::currentState);
        unexpected.measured.truePeakDbtp=-.5;unexpected.detail.truePeak.value=-.5;
        expect(faulted.job.finish(unexpected.context.ticket,faulted.revision,ErrorCode::none,"injected differing measured outcome")&&faulted.acceptMeasurement(unexpected),"accept explicitly injected independent measured disagreement");
        expect(faulted.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::fail&&faulted.currentStateDocument.staticGainDb==plan.staticGainDb
            &&faulted.currentStateDocument.sectionGains->audioKey()==plan.sectionGains->audioKey()&&!faulted.sectionProposal,"actual verification failure displayed without rollback retry or new correction");
        auto privateForgery=*candidate->assessment();privateForgery.context=faulted.measurementContext(MeasurementDomain::privateCandidate);
        expect(!faulted.acceptMeasurement(privateForgery),"private candidate domain cannot become Current through measurement acceptance");
        auto unavailable=request;auto missingSnapshot=std::make_shared<AnalysisSnapshot>(*request.current);
        missingSnapshot->measured.integratedLufs.reset();unavailable.current=missingSnapshot;JobControl missing;
        auto rejected=generateSectionProposal(unavailable,missing);expect(rejected.value&&(*rejected.value)->status()==SectionProposalStatus::unavailable&&(*rejected.value)->passes()==0,"missing required metric rejects generation before measurement");
        auto whole=request;whole.annotations=mapFor(rate,rate*10,1);whole.annotations.inputIdentity=request.annotations.inputIdentity;
        JobControl extentControl;auto extent=generateSectionProposal(whole,extentControl);
        expect(extent.value&&(*extent.value)->status()==SectionProposalStatus::noSolution&&(*extent.value)->passes()==1&&(*extent.value)->text().contains("above 20%"),"whole-programme local attenuation rejected honestly after attribution without candidate pass claim");
        auto manual=request;manual.pendingManualDraft=true;JobControl manualControl;auto blocked=generateSectionProposal(manual,manualControl);
        expect(blocked.value&&(*blocked.value)->status()==SectionProposalStatus::unavailable&&(*blocked.value)->passes()==0,"pending manual draft is not overwritten");
        JobControl active;auto activeJob=std::async(std::launch::async,[&]{return generateSectionProposal(request,active);});
        while(active.progress.load()<.15f&&activeJob.wait_for(std::chrono::milliseconds(1))!=std::future_status::ready)std::this_thread::sleep_for(std::chrono::milliseconds(1));
        expect(active.progress.load()>=.15f&&activeJob.wait_for(std::chrono::milliseconds(0))!=std::future_status::ready,"cancel test reaches active private measurement");
        const auto start=juce::Time::getMillisecondCounterHiRes();active.cancel=true;const auto cancelledResult=activeJob.get();
        expect(cancelledResult.error==ErrorCode::cancelled&&!cancelledResult.value&&juce::Time::getMillisecondCounterHiRes()-start<=1000,"active private measurement cancellation acknowledged within one second with no draft publication");

    }
    const auto valid=s.currentState.analysis;const auto finalRevision=s.revision;JobControl noopControl;
    auto noop=generateSectionProposal(s.sectionProposalRequest(),noopControl);
    expect(noop.value&&(*noop.value)->status()==SectionProposalStatus::noChange&&(*noop.value)->passes()==0,"already passing returns genuine no-op with no private passes");
    expect(s.acceptSectionProposal(*noop.value)&&s.applySectionProposal()&&s.revision==finalRevision&&s.currentState.analysis==valid,"no-op preserves gain, revision and valid measurement");
    expect(!s.hasActionableSectionProposal(),"already-pass explanation does not create an unapplied audio draft or block saving");
    const auto noOpDocument=sf::capture(s);expect(noOpDocument.current.staticGainDb==plan.staticGainDb,"no-op session capture preserves exact settings without private authority");
    s.discardSectionProposal();const auto doc=sf::capture(s);auto file=root.getChildFile("applied.levelscope");JobControl io;
    expect(sf::save(doc,file,io).succeeded(),"applied proposal saves as existing exact M12 plan");
    expect(root.getChildFile("measurement-only.json").replaceWithText(juce::JSON::toString(ex::report(ex::capture(s),nullptr),true)),"retain measured Source and applied Current report without artifact authority");
    auto read=sf::read(file,io);expect(read.document.has_value(),"read applied session");auto prepared=sf::verify(*read.document,io);
    expect(prepared.prepared&&prepared.prepared->model.currentStateDocument.sectionGains->audioKey()==plan.sectionGains->audioKey()
        &&!prepared.prepared->model.sectionProposal,"safe reopen restores exact processing without transient proposal authority");
    if(rate==48000)
    {
        auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));OutputRouting muted;
        s.selectAuditionDomain(MeasurementDomain::currentState);auto auditionRequest=s.auditionRequest();audition.synchronize(auditionRequest,muted);
        expect(audition.connect(rate).isEmpty(),"proposal injected audition device only");
        std::unique_ptr<ProgrammeReader> sourceReader;
        // Open the exact original programme through the established verified-copy boundary.
        JobControl sourceJob;auto copy=prepareVerifiedSource(s.measurementContext(MeasurementDomain::source),sourceJob);
        expect(copy.value.has_value(),"verified source for independent sample oracle");sourceReader=(*copy.value)->createReader();
        auto captureFile=root.getChildFile("audition-capture.f32").createOutputStream();
        expect(bool(captureFile),"owned cross-platform capture file opens");
        for(int64_t first:{int64_t(0),int64_t(rate*5-1000),int64_t(rate*10-4096)})
        {
            expect(audition.seek(first,auditionRequest,muted).isEmpty(),"proposal exact source-coordinate seek");Capture captured(channels,4096,first);audition.setCaptureObserver(&captured);
            expect(audition.play(auditionRequest,muted,false,true).isEmpty(),"proposal playback only via injected muted device");awaitPlaying(audition,*injected);juce::AudioBuffer<float> out(maxChannels,256);
            for(int n=0;n<16;++n){injected->tick(out);audition.poll();std::this_thread::sleep_for(std::chrono::milliseconds(1));}audition.stop();audition.setCaptureObserver(nullptr);
            juce::AudioBuffer<float> raw(channels,4096);expect(sourceReader->read(raw.getArrayOfWritePointers(),channels,first,4096),"read original oracle span");
            expect(captured.valid&&captured.next==first+4096&&sameSamples(captured.data,oracleAudio(raw,plan,first))&&audition.underruns()==0,"applied proposal pre-monitor samples equal independent section/global oracle at seek and EOF");
            for(int ch=0;ch<channels;++ch)expect(captureFile->write(captured.data.getReadPointer(ch),4096*4),"retain bounded pre-monitor capture for cross-platform comparison");
        }
        captureFile->flush();expect(captureFile->getStatus().wasOk(),"bounded cross-platform capture finalized");captureFile.reset();
        audition.disconnect();
        const auto exportRequest=ex::capture(s);JobControl exporting;auto rendered=ex::render(exportRequest,root,exporting);
        expect(rendered.value&&ex::artifactEvaluation(**rendered.value,exportRequest).verdict==proposal::Verdict::pass,"actual proposal PCM24 artifact independently remeasured and passes");
        expect(ex::publish(**rendered.value,exportRequest,root.getChildFile("export"),ex::approvalIdentity(**rendered.value,exportRequest),exporting).value.has_value(),"explicit local verified export includes actual report/hash");
    }
    if(rate==48000&&mode==0)
    {
        auto globalOnly=s;expect(globalOnly.setCurrentGainDb(plan.staticGainDb-.5),"explicit global-only baseline below target");measuredComposition(globalOnly,MeasurementDomain::currentState);
        const auto fixed=globalOnly.currentStateDocument.sectionGains->audioKey();JobControl globalControl;auto global=generateSectionProposal(globalOnly.sectionProposalRequest(),globalControl);
        expect(global.value&&(*global.value)->status()==SectionProposalStatus::globalOnly&&(*global.value)->passes()==0,"feasible global-only decision precedes section work with zero private passes");
        expect(globalOnly.acceptSectionProposal(*global.value)&&globalOnly.applySectionProposal()&&globalOnly.currentStateDocument.sectionGains->audioKey()==fixed,"global-only Apply preserves all existing section gains");
        measuredComposition(globalOnly,MeasurementDomain::currentState);expect(globalOnly.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"global-only independent measured verification passes");
    }
    auto stale=request;stale.currentContext.processing.staticGainDb=.1;JobControl badControl;auto bad=generateSectionProposal(stale,badControl);
    expect(bad.value&&(*bad.value)->status()==SectionProposalStatus::unavailable&&(*bad.value)->passes()==0,"mismatched measured plan rejected before audio work");
    JobControl cancelled;cancelled.cancel=true;expect(generateSectionProposal(request,cancelled).error==ErrorCode::cancelled,"pre-cancellation preserves previous session");
    auto caseRow=juce::var(new juce::DynamicObject());auto* cr=caseRow.getDynamicObject();cr->setProperty("rate",rate);cr->setProperty("channels",channels);cr->setProperty("mode",mode);
    cr->setProperty("private_passes",candidate->passes());cr->setProperty("gain",plan.staticGainDb);cr->setProperty("integrated",*verification.measured.integratedLufs);cr->setProperty("true_peak",*verification.measured.truePeakDbtp);cases.add(caseRow);
    if(rate==48000)cr->setProperty("capture_sha256",juce::SHA256(root.getChildFile("audition-capture.f32")).toHexString());
    SnapshotStorage noSpace;noSpace.availableBytes=[](const juce::File&){return int64_t(0);};JobControl disk;
    expect(generateSectionProposal(request,disk,noSpace).error==ErrorCode::insufficientDiskSpace,"verified-copy disk exhaustion fails without replacing session");
    if(rate==48000&&mode==0)
    {
        SnapshotStorage unwritable;unwritable.directory=root.getChildFile("failed-private-copy");
        expect(unwritable.directory.createDirectory().wasOk(),"create owned injected-write-failure area");
        unwritable.createWriter=[](const juce::File&)->std::unique_ptr<SnapshotWriter>{return {};};JobControl writeFailure;
        const auto failed=generateSectionProposal(request,writeFailure,unwritable);
        expect(failed.error==ErrorCode::writeFailed&&!failed.value&&unwritable.directory.findChildFiles(juce::File::findFiles,true).isEmpty(),"private-copy write failure publishes no candidate and cleans owned partial files");
    }
}
void assessmentBranches(const juce::File& root)
{
    auto failure=proposalFixture(root.getChildFile("insufficient-local-cap"),48000,6,0,2.f);
    const auto current=failure.currentState.analysis;const auto revision=failure.revision;JobControl control;
    auto failed=generateSectionProposal(failure.sectionProposalRequest(),control);
    expect(failed.value&&(*failed.value)->status()==SectionProposalStatus::noSolution&&(*failed.value)->passes()==2
        &&(*failed.value)->assessment()&&(*failed.value)->selection().corrections.front().uncoveredDb>0,
        "capped envelope measured once and fails honestly with uncovered deficit and no second candidate");
    expect(failure.acceptSectionProposal(*failed.value)&&!failure.canApplySectionProposal()&&!failure.applySectionProposal()
        &&failure.currentState.analysis==current&&failure.revision==revision,"measured private failure cannot change Current or enable Apply");
    auto feasible=proposalFixture(root.getChildFile("identical-final-pass"),48000,6,0);
    expect(feasible.setCurrentGainDb(-23-*feasible.currentState.integratedLufs),"explicit baseline at intended loudness with excess peak");
    measuredComposition(feasible,MeasurementDomain::currentState);const auto gain=feasible.currentStateDocument.staticGainDb;JobControl same;
    auto result=generateSectionProposal(feasible.sectionProposalRequest(),same);
    expect(result.value&&(*result.value)->status()==SectionProposalStatus::ready&&(*result.value)->passes()==2
        &&(*result.value)->plan().staticGainDb==gain,"already-acceptable exact local candidate skips identical final measurement without centering gain");
}
}
// Resource helpers only; no historical runner invoked.
#include "SectionProposalLongTests.h"
#include "SectionProposalCinemaTests.h"
#include "SectionProposalIntegrationTests.h"
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI gui;juce::String rootPath,output,longSession,cinema,videoDirectory,sessionPath,exchange;bool passed=false,smoke=false;juce::String failure;
    for(int n=1;n<argc;++n){juce::String arg(argv[n]);if(arg=="--smoke")smoke=true;else if(arg=="--root"&&n+1<argc)rootPath=argv[++n];else if(arg=="--video"&&n+1<argc)videoDirectory=argv[++n];else if(arg=="--session"&&n+1<argc)sessionPath=argv[++n];else if(arg=="--exchange"&&n+1<argc)exchange=argv[++n];else if(arg=="--cinema"&&n+1<argc)cinema=argv[++n];else if(arg=="--long-session"&&n+1<argc)longSession=argv[++n];else if(arg=="--output"&&n+1<argc)output=argv[++n];else return 2;}
    try
    {
        const juce::File root(rootPath);expect(rootPath.isNotEmpty()&&!root.exists()&&root.createDirectory().wasOk(),"fresh M13 fixture root required");
        if(videoDirectory.isNotEmpty())proposalVideo(juce::File(sessionPath),juce::File(videoDirectory));else if(exchange.isNotEmpty())proposalExchange(juce::File(exchange));else if(cinema.isNotEmpty())cinemaEvaluation(root,juce::File(cinema));else if(longSession.isNotEmpty())longSectionProposal(root,juce::File(longSession));else {
        policyCases();collectorCases();assessmentBranches(root.getChildFile("assessment-branches"));
        for(int rate:{44100,48000,96000,192000})for(int channels:{2,6})for(int mode=0;mode<(channels==2?4:6);++mode)
        {
            if(smoke&&(rate!=48000||mode!=0))continue;
            workflow(root.getChildFile(juce::String(rate)+"-"+juce::String(channels)+"-"+juce::String(mode)),rate,channels,mode);
        }
        }
        passed=true;
    }
    catch(const std::exception& e){failure=e.what();std::cerr<<failure<<'\n';}
    auto v=juce::var(new juce::DynamicObject());auto* r=v.getDynamicObject();r->setProperty("suite","M13 section-led proposals");r->setProperty("passed",passed);r->setProperty("failure",failure);
    r->setProperty("checks",checks);r->setProperty("cases",cases);r->setProperty("phase",videoDirectory.isNotEmpty()?"video":exchange.isNotEmpty()?"exchange":cinema.isNotEmpty()?"cinema":longSession.isNotEmpty()?"long":smoke?"development-smoke":"audio");r->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);r->setProperty("os",juce::SystemStats::getOperatingSystemName());r->setProperty("long_workflow",longEvidence);r->setProperty("physical_devices_opened",false);r->setProperty("video_loaded",videoDirectory.isNotEmpty());r->setProperty("video_cases",videoEvidence);r->setProperty("exchange_cases",exchangeEvidence);
#if defined(__arm64__) || defined(__aarch64__)
    r->setProperty("architecture","arm64");
#else
    r->setProperty("architecture","x86_64");
#endif
    if(output.isNotEmpty())juce::File(output).replaceWithText(juce::JSON::toString(v,true));std::cout<<"M13 "<<(passed?"PASS":"FAIL")<<" "<<checks.size()<<" checks\n";return passed?0:1;
}

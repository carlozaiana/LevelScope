#include "CompositionTestFixtures.h"
#include "Core/Proposal/EventPlanning.h"
#include "Standalone/StandaloneRegionProposal.h"
#include <thread>
#include <type_traits>
using namespace levelscope::proposal;
namespace {
using Plan=levelscope::processing::SectionGainPlan;
Plan unity(int seconds=60,int rate=48000){Plan p;p.frames=int64_t(seconds)*rate;p.sampleRate=rate;p.inputIdentity=juce::String::repeatedString("a",64);p.intervals={{0,p.frames,0}};return p;}
std::shared_ptr<ProgramAnalysis> analysis(const Plan& p)
{
 auto a=std::make_shared<ProgramAnalysis>();a->identity="analysis";a->snapshotIdentity="snapshot";a->sampleRate=p.sampleRate;a->frames=p.frames;a->lraStable=true;a->p10=-30;a->p95=-20;
 for(int64_t end=3*p.sampleRate;end<=p.frames;end+=p.sampleRate/10){DiagnosticWindow w;w.begin=end-3*p.sampleRate;w.end=end;w.lufs=-20;w.energy=.001;w.included=true;a->lraWindows.push_back(w);}return a;
}
RegionConstraints area(const Plan& p,RegionSpan h)
{auto s=RegionPlanningSettings::eventDefaults(RegionMode::selection,p.inputIdentity);s.selection=h;auto passages=findLraPassages(analysis(p),5);return resolveRegionConstraints(s,{},*passages,p);}
RegionAssignment assignment(RegionSpan span,std::vector<int64_t> cuts={},std::vector<double> gains={2}){RegionAssignment a;a.geometry.support=span;a.geometry.cuts=cuts;a.reductions=gains;return a;}
EventEvidence evidence(const RegionConstraints& c,const ProgramAnalysis& a)
{EventEvidence e;e.finished=true;e.identity="evidence";e.constraintsIdentity=c.identity;e.currentSnapshot=a.snapshotIdentity;e.definitionDigest=c.settings.eventFeatures.digest();e.frames=c.frames;e.sampleRate=c.sampleRate;e.incompleteRoots={false};EventGroup g;g.id="event";g.root=0;g.onset={15*c.sampleRate,16*c.sampleRate};g.decay={24*c.sampleRate,25*c.sampleRate};g.lastAttack={20*c.sampleRate,21*c.sampleRate};g.corroborated=g.needAssociated=g.endingSupported=true;e.groups={g};e.needCuts.resize(1);return e;}
void core()
{
 static_assert(!std::is_default_constructible_v<RegionAssessmentReceipt>);static_assert(!std::is_copy_constructible_v<RegionAssessmentReceipt>);
 auto p=unity();auto a=analysis(p);auto c=area(p,{10*48000,30*48000});auto e=evidence(c,*a);ProposalPolicy policy;policy.durationThousandths=600;
 auto legacy=RegionPlanningSettings::defaults(RegionMode::automatic,p.inputIdentity);expect(legacy.schema==1&&!legacy.preferCompleteEvent&&!legacy.allowAutomaticMargins,"legacy defaults stay legacy");
 expect(RegionPlanningSettings::decode(legacy.encode()).digest()==legacy.digest(),"legacy exact round trip");
 expect(RegionPlanningSettings::decode(c.settings.encode()).digest()==c.settings.digest(),"v2 full independent definition round trip");
 auto bad=c.settings;bad.compareAlternatives=false;expect(bad.validationError().isNotEmpty(),"new mode cannot stop at first partial pass");
 bad=RegionPlanningSettings::eventDefaults(RegionMode::automatic,p.inputIdentity);bad.leftMarginMs=1000;expect(bad.validationError().isNotEmpty(),"context does not authorize automatic margin");bad.allowAutomaticMargins=true;expect(bad.validationError().isEmpty(),"explicit automatic margin permission");
 auto json=c.settings.encode();json["eventFeatures"].getDynamicObject()->setProperty("hopMs",7);bool rejected=false;try{RegionPlanningSettings::decode(json);}catch(...){rejected=true;}expect(rejected,"unknown feature hypothesis fails closed");
 auto complete=constructRegionPlan(p,0,c,policy,{assignment({12*48000,28*48000})},0);expect(complete.error.isEmpty(),"complete envelope legal");auto full=assessEventCoverage(p,complete,c,e,*a,5);expect(full.complete(),"correction before onset maintained through supported decay");
 auto onTime=constructRegionPlan(p,0,c,policy,{assignment({e.groups[0].onset.begin-4799,28*48000})},0);expect(onTime.error.isEmpty()&&assessEventCoverage(p,onTime,c,e,*a,5).complete(),"exact final ramp sample at earliest onset is complete");auto oneLate=constructRegionPlan(p,0,c,policy,{assignment({e.groups[0].onset.begin-4798,28*48000})},0);expect(oneLate.error.isEmpty()&&assessEventCoverage(p,oneLate,c,e,*a,5).category==EventCoverageClass::knownPartial,"one-frame late full correction is partial without epsilon allowance");
 auto late=constructRegionPlan(p,0,c,policy,{assignment({17*48000,28*48000})},0);auto partial=assessEventCoverage(p,late,c,e,*a,5);expect(partial.category==EventCoverageClass::knownPartial,"ramp wholly after first hit is partial");
 auto early=constructRegionPlan(p,0,c,policy,{assignment({12*48000,23*48000})},0);expect(assessEventCoverage(p,early,c,e,*a,5).category==EventCoverageClass::knownPartial,"return before last decay rejected");
 auto inside=constructRegionPlan(p,0,c,policy,{assignment({12*48000,28*48000},{18*48000,22*48000},{2,1,2})},0);expect(inside.error.isEmpty()&&assessEventCoverage(p,inside,c,e,*a,5).complete(),"positive internal schedule can cover continuing event");
 auto hole=constructRegionPlan(p,0,c,policy,{assignment({12*48000,28*48000},{18*48000,22*48000},{2,0,2})},0);expect(hole.error.isEmpty()&&assessEventCoverage(p,hole,c,e,*a,5).category==EventCoverageClass::knownPartial,"zero hole in continuing event is partial");
 auto unknown=e;unknown.groups[0].incomplete=true;expect(assessEventCoverage(p,late,c,unknown,*a,5).category==EventCoverageClass::partialAndUncertain,"known missed entrance and incomplete evidence both retained");
 expect(assessEventCoverage(p,complete,c,unknown,*a,5).category==EventCoverageClass::uncertain,"uncertainty earns no complete credit");
 auto unrelated=*a;for(auto& w:unrelated.lraWindows)w.lufs=-30;expect(assessEventCoverage(p,complete,c,e,unrelated,5).category==EventCoverageClass::notApplicable,"no need association earns no event credit");
 auto stale=e;stale.constraintsIdentity="other";expect(!assessEventCoverage(p,complete,c,stale,*a,5).complete(),"stale constraints reject event evidence");
 RegionReviewScore shortScore,longScore;shortScore.complete=shortScore.allPass=longScore.complete=longScore.allPass=true;shortScore.id="short";longScore.id="long";shortScore.activity=late.activity;longScore.activity=complete.activity;
 partial.signature={1,2};expect(selectCompleteEventPass({shortScore,longScore},{partial,full})==1,"complete pass beats shorter partial across signatures");
 auto second=full;second.signature={1,2};auto cheaper=longScore;cheaper.id="different";cheaper.activity.affectedFrames-=1;expect(selectCompleteEventPass({longScore,cheaper},{full,second})==1,"duration ranks supported complete alternatives without event count reward");
 auto emptyCoverage=full;emptyCoverage.signature.clear();expect(!selectCompleteEventPass({longScore},{emptyCoverage}),"empty signature cannot earn complete-event reward");expect(selectCompleteEventPass({longScore,shortScore},{full,partial})==0,"permuting candidate order preserves complete recommendation");longScore.allPass=false;expect(!selectCompleteEventPass({longScore},{full}),"coverage never promotes measured failure");longScore.allPass=true;
 auto mixed=e;auto incidental=e.groups[0];incidental.id="incidental";incidental.needAssociated=false;incidental.incomplete=true;mixed.groups.push_back(incidental);expect(assessEventCoverage(p,complete,c,mixed,*a,5).category==EventCoverageClass::uncertain,"uncertain incidental group prevents complete credit");
 auto zero=constructRegionPlan(p,0,c,policy,{},0);expect(assessEventCoverage(p,zero,c,e,*a,5).category==EventCoverageClass::notApplicable,"zero envelope earns no complete credit");
 // Compare coverage with every actual sample across old automation, including
 // adjacent old boundaries that shorten the processor's nominal 100 ms ramps.
 auto automated=p;automated.intervals={{0,15*48000-480,-1},{15*48000-480,15*48000+480,-.5},{15*48000+480,25*48000,-1.5},{25*48000,25*48000+960,-.5},{25*48000+960,p.frames,-1}};
 auto oldConstraints=area(automated,{10*48000,30*48000});auto oldEvidence=evidence(oldConstraints,*a);
 auto crossing=constructRegionPlan(automated,0,oldConstraints,policy,{assignment({12*48000,28*48000})},0);
 expect(crossing.error.isEmpty()&&automated.radius(1)<2400,"old close boundaries retain shortened processor ramps");
 bool fullyCovered=true;for(int64_t f=oldEvidence.groups[0].onset.begin;f<oldEvidence.groups[0].decay.end;++f)
 {const double b=automated.gainAt(f);fullyCovered&=b-crossing.envelope.gainAt(f)>=b-(b-2.);}
 expect(assessEventCoverage(automated,crossing,oldConstraints,oldEvidence,*a,5).complete()==fullyCovered,"coverage agrees with exhaustive exact additional-envelope samples through existing ramps");
 auto oldUnchanged=constructRegionPlan(automated,0,oldConstraints,policy,{},0);
 expect(assessEventCoverage(automated,oldUnchanged,oldConstraints,oldEvidence,*a,5).category==EventCoverageClass::notApplicable,"pre-existing attenuation alone earns no new complete-event credit");
 // A compact, valid processor plan returns during the supported ending. Its
 // shortened ramp must not evade the decay-end obligation.
 auto shortReturn=complete;shortReturn.envelope.intervals={{0,12*48000,0},{12*48000,25*48000-480,-2},{25*48000-480,25*48000+480,-1},{25*48000+480,p.frames,0}};
 expect(shortReturn.envelope.validationError().isEmpty()&&shortReturn.envelope.radius(2)<2400,"shortened return fixture uses a valid canonical processor");
 expect(assessEventCoverage(p,shortReturn,c,e,*a,5).category==EventCoverageClass::knownPartial,"shortened ramp returning inside the decay is partial");
 auto protectedArea=c;protectedArea.settings.protectedLocal={{14*48000,16*48000}};expect(constructRegionPlan(p,0,protectedArea,policy,complete.assignments,0).error.isNotEmpty(),"complete preference cannot cross protection");
 LraPassages passages;passages.analysis=a;passages.identity="ownership";
 for(int n=0;n<10;++n)passages.high.push_back({int64_t(1+n*5)*48000,int64_t(4+n*5)*48000,10,0,double(100-n),-20});
 auto automatic=RegionPlanningSettings::eventDefaults(RegionMode::automatic,p.inputIdentity);automatic.allowAutomaticMargins=true;automatic.leftMarginMs=automatic.rightMarginMs=3000;
 auto owned=resolveRegionConstraints(automatic,{},passages,p);expect(owned.roots.size()==8&&owned.omittedRoots==2,"ownership formed against all roots before eight-root pruning");
 for(const auto& root:owned.roots){const auto n=root.hint.begin/48000;expect(root.allowed.end<=root.hint.end+48000,"outward permission clips at fixed midpoint even beside omitted root");(void)n;}
 const auto frozen=owned.roots[0].allowed;owned.roots.erase(owned.roots.begin()+1);expect(owned.roots[0].allowed==frozen,"dropping neighbour cannot borrow its ownership");
 EventAssessmentChoice original{shortScore,partial,false,true,true},whole{longScore,full,false,false,true},other{shortScore,partial,false,false,true};
 expect(chooseEventAssessment({original,whole,other},true,false)==1,"eligible complete first");whole.screening=false;expect(chooseEventAssessment({original,whole,other},true,false)==0,"original first when complete forecast fails");original.measured=true;expect(chooseEventAssessment({original,whole,other},false,false)==1,"partial pass cannot suppress unmeasured complete challenger");original.measured=false;expect(chooseEventAssessment({original,whole,other},false,true)==0,"complete pass reserves original comparison");original.measured=whole.measured=other.measured=true;expect(!chooseEventAssessment({original,whole,other},false,true),"exhausted pool has no fabricated assessment");
}
std::shared_ptr<const EventEvidence> extract(int rate,int hop,int chunk,int mode,int seconds=10)
{
 auto p=unity(seconds,rate);auto a=analysis(p);auto c=area(p,{3*rate,(seconds-3)*rate});c.settings.eventFeatures.hopMs=hop;c.identity="feature-constraints";
 auto formed=formRegions(c,*a,5);const int channels=(mode==3||mode==6||mode==7)?6:2;auto roles=OfflineAnalysisService::measurementRoles(layoutFor(channels));if(mode==7)std::reverse(roles.begin(),roles.end());EventFeatureCollector collector(c,formed,*a,roles,5,"binding","source",regionPlanIdentity(p,0));
 juce::AudioBuffer<float> source(channels,chunk),current(channels,chunk);
 for(int64_t first=0;first<p.frames;first+=chunk){const int n=int(std::min<int64_t>(chunk,p.frames-first));for(int ch=0;ch<channels;++ch)for(int k=0;k<n;++k){double t=double(first+k)/rate;double gain=t<3?0:t<6?.15:t<6.7?.15*std::exp(-(t-6)*15):0;double v=gain*std::sin(2*juce::MathConstants<double>::pi*60*t);if(mode==1&&ch==1)v=-v;if(mode==2)v=0;if(mode==3&&ch!=3)v=0;if(mode==4)v=.1;if(mode==9)v=t>=2&&std::fmod(t-2,1.4)<.15?.15*std::sin(2*juce::MathConstants<double>::pi*60*t):0;if(mode==10&&t>=5)v=std::numeric_limits<double>::quiet_NaN();if(mode==11&&t>=3)v=.15*std::sin(2*juce::MathConstants<double>::pi*60*t);if(mode==8)v=.02*std::sin(2*juce::MathConstants<double>::pi*60*t);if(mode==5)v*=20;if((mode==6&&ch==3)||(mode==7&&ch==2))v=0;source.setSample(ch,k,float(v));current.setSample(ch,k,float(v*(mode==8&&t>=4&&t<6?1.:.5)));}collector.push(source.getArrayOfReadPointers(),current.getArrayOfReadPointers(),n);}return collector.finish();
}
void features()
{
 for(int rate:{44100,48000,96000,192000})for(int hop:{5,10,20})
 {
  auto a=extract(rate,hop,4096,0),b=extract(rate,hop,137,0);expect(a->identity==b->identity,"feature filters/cadence invariant to chunking");expect(a->highWaterBytes<=10*1024*1024,"feature storage within cap");
  auto g=std::find_if(a->groups.begin(),a->groups.end(),[](auto& x){return x.needAssociated&&x.corroborated&&x.endingSupported&&!x.incomplete;});
  expect(g!=a->groups.end(),"known synthetic attack and decay supported");expect(g->firstNeedWindow>=0&&g->needWindowCount>0&&g->onsetPeakFrame>=g->onset.begin&&g->onsetPeakFrame<g->onset.end,"raw attack witness and exact baseline-window association retained");expect(g->onset.begin<=3*rate&&g->onset.end>=3*rate&&g->decay.end>6*rate,"onset bracket includes first sample and return follows decay");
  auto row=juce::var(new juce::DynamicObject());row.getDynamicObject()->setProperty("rate",rate);row.getDynamicObject()->setProperty("hop_ms",hop);row.getDynamicObject()->setProperty("onset_begin",g->onset.begin);row.getDynamicObject()->setProperty("decay_end",g->decay.end);row.getDynamicObject()->setProperty("feature_bytes",int64_t(a->highWaterBytes));cases.add(row);
 }
 auto phase=extract(48000,10,4096,1),normal=extract(48000,10,4096,0);expect(phase->identity==normal->identity,"antiphase stereo preserves per-channel event evidence");expect(extract(48000,10,4096,2)->groups.empty(),"silence does not fabricate events");expect(extract(48000,10,4096,3)->groups.empty(),"LFE-only activity cannot create correction event");
 auto over=extract(48000,10,4096,5);expect(!over->groups.empty(),"finite over-full-scale composition remains analyzable without attenuation");auto surround=extract(48000,10,4096,6),reordered=extract(48000,10,4096,7);expect(!surround->groups.empty()&&surround->groups.size()==reordered->groups.size()&&surround->groups[0].onset==reordered->groups[0].onset&&surround->groups[0].decay==reordered->groups[0].decay,"ordered 5.1 role handling preserves event extents after explicit permutation");
 auto overflow=extract(48000,5,4096,9,120);expect(overflow->groups.size()==64&&overflow->omittedGroups>0&&overflow->incompleteRoots[0],"summary overflow marks root incomplete without raising storage cap");auto openEnd=extract(48000,10,4096,11);expect(std::any_of(openEnd->groups.begin(),openEnd->groups.end(),[](const auto& g){return !g.endingSupported&&g.incomplete;}),"EOF never invents a supported decay");bool nonfinite=false;try{extract(48000,10,4096,10);}catch(...){nonfinite=true;}expect(nonfinite,"non-finite Source evidence fails safely");
 auto processingOnly=extract(48000,10,4096,8);expect(std::any_of(processingOnly->cells.begin(),processingOnly->cells.end(),[](const auto& c){return c.begin>=4*48000&&(c.flags&2);}),"Current-only rise is labelled processing context");expect(std::none_of(processingOnly->groups.begin(),processingOnly->groups.end(),[](const auto& g){return g.onset.begin>48000&&g.corroborated;}),"existing gain automation cannot create Source event authority");
 auto dc=extract(48000,10,4096,4);expect(std::none_of(dc->groups.begin(),dc->groups.end(),[](const auto& g){return g.onset.begin>24000&&g.corroborated;}),"constant DC is not a later musical event");
 auto p=unity(10800);auto a=analysis(p);RegionConstraints c;c.frames=p.frames;c.sampleRate=p.sampleRate;c.settings=RegionPlanningSettings::eventDefaults(RegionMode::automatic,p.inputIdentity);c.settings.eventFeatures.hopMs=5;c.identity="stress";RegionFormation formation;formation.nominations.resize(8);
 for(int n=0;n<8;++n){PlanningRoot root;root.hint={int64_t(100+n*1200)*48000,int64_t(1100+n*1200)*48000};root.allowed=root.hint;c.roots.push_back(root);for(int k=0;k<12;++k)formation.nominations[n].push_back({root.hint.begin+int64_t(k+1)*70*48000,1,"bounded nomination"});}
 EventFeatureCollector stress(c,formation,*a,{measurement::Role::left,measurement::Role::right},5,"stress","source",regionPlanIdentity(p,0));expect(stress.reservedBytes()<=10*1024*1024,"eight-root three-hour-sized 5ms reservation fits without PCM decode");
 auto row=juce::var(new juce::DynamicObject());row.getDynamicObject()->setProperty("three_hour_feature_reserved",int64_t(stress.reservedBytes()));cases.add(row);
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

void workflow(const juce::File& root,int channels,bool heterogeneous,bool fallback=false)
{
 auto s=programme(root,channels,heterogeneous);auto requirements=s.jobRequirements;requirements.maximumLraTenths=49;requirements.provenance="Event workflow test";juce::String error;expect(s.setJobRequirements(requirements,error),"test requirement");
 auto settings=RegionPlanningSettings::eventDefaults(RegionMode::automatic,sectionInputIdentity(s.measurementContext(MeasurementDomain::source)));settings.allowAutomaticMargins=true;settings.leftMarginMs=settings.rightMarginMs=1000;if(fallback){settings=RegionPlanningSettings::eventDefaults(RegionMode::selection,settings.inputIdentity);settings.selection={int64_t(48000)*302/10,int64_t(48000)*36};settings.allowGlobal=true;}expect(s.setRegionPlanning(settings,error),"explicit event preference");auto policy=s.proposalPolicy;policy.durationThousandths=600;expect(s.setProposalPolicy(policy,error),"explicit 60 percent experiment, default unchanged");
 const auto original=s;const auto request=s.sectionProposalRequest();JobControl job;auto outcome=generateSectionProposal(request,job);if(!outcome.value)std::cerr<<outcome.message<<std::endl;expect(bool(outcome.value),"event service complete");auto proposed=*outcome.value;root.getChildFile("proposal.txt").replaceWithText(proposed->text());
 expect(proposed->methodVersion()==eventSearchVersion&&proposed->passes()<=3,"versioned bounded service");auto review=proposed->regionReview();expect(review&&review->eventEvidence&&review->proxyEvaluations<=64&&review->planningBytes<=64*1024*1024,"bounded auxiliary/formation/search evidence");expect(!proposed->assessment(),"private full populations not retained");
 expect(s.job.begin(JobStage::proposal,s.revision),"ordinary UI job begins after request capture");auto ticket=s.job.ticket;expect(s.job.finish(ticket,s.revision,ErrorCode::none,"test complete"),"ordinary completion ticket accepted");expect(s.acceptSectionProposal(proposed),"matching immutable completion");expect(s.revision==original.revision&&s.currentState.analysis==original.currentState.analysis,"planning preserves Current");
 size_t passing=0;for(auto& a:review->attempts){expect(a.receipt&&a.receipt->measured().momentary.empty()&&a.receipt->measured().shortTerm.empty()&&a.receipt->measured().lraWindows.empty(),"compact receipt holds no window populations");if(a.evaluation.verdict==Verdict::pass)++passing;}
 expect(passing>0,"representative candidate has measured passing alternative");
 auto i=std::find_if(review->attempts.begin(),review->attempts.end(),[](auto& a){return a.evaluation.verdict==Verdict::pass;});
 if(fallback)expect(proposed->status()==SectionProposalStatus::reviewRequired&&!review->recommendedAttempt,"passing partial-only set has no recommendation");if(proposed->status()==SectionProposalStatus::reviewRequired)expect(!s.canApplySectionProposal()&&!s.applySectionProposal(),"inspection-only review cannot Apply");
 expect(!s.selectMeasuredRegionAlternative("wrong-result",i->attemptIdentity)&&!s.selectMeasuredRegionAlternative(review->resultIdentity,"wrong-attempt"),"foreign result and fabricated attempt rejected");
 expect(s.selectMeasuredRegionAlternative(review->resultIdentity,i->attemptIdentity),"explicit owned measured alternative selection");proposed=s.sectionProposal;expect(s.revision==original.revision&&s.currentState.analysis==original.currentState.analysis,"selection changes no audio or measurement");
 auto stale=s;stale.setCurrentGainDb(.1);expect(!stale.canApplySectionProposal(),"processing invalidates selected authority");auto changed=request;changed.currentContext.sessionId="new-runtime";expect(!selectMeasuredEventAlternative(changed,proposed,review->resultIdentity,i->attemptIdentity).value,"New/Open identity rejects late selection");changed=request;changed.pendingManualDraft=true;expect(validateEventRegionProposalForApply(changed,*proposed).isNotEmpty(),"manual draft rejects old authority");
 changed=request;changed.currentContext.fileLayout.roles[0]=juce::AudioChannelSet::centre;expect(validateEventRegionProposalForApply(changed,*proposed).isNotEmpty(),"changed ordered roles reject receipt");changed=request;changed.regionSettings->eventFeatures.hopMs=5;expect(validateEventRegionProposalForApply(changed,*proposed).isNotEmpty(),"feature definition change invalidates proposal, not measurement");
 changed=request;changed.currentContext.profile.limits.maxLoudnessRangeLu=.1;expect(validateEventRegionProposalForApply(changed,*proposed).isNotEmpty(),"changed requirements cannot reuse old verdict authority");
 s.setOutputRoute(0,-1);expect(s.canApplySectionProposal(),"monitor muting preserves all-channel receipt");expect(s.applySectionProposal()&&s.currentStateNeedsRemeasurement(),"explicit Apply makes Current unverified");
 const auto measured=measuredComposition(s,MeasurementDomain::currentState);const auto& expected=proposed->eventAssessment()->measured();expect(measured.detail.integrated.value==expected.integrated.value&&measured.detail.lra.value==expected.lra.value&&measured.detail.channelTruePeaks==expected.channelTruePeaks&&measured.detail.channelSamplePeaks==expected.channelSamplePeaks&&measured.detail.truePeakEighthFrames==expected.truePeakEighthFrames&&measured.detail.samplePeakFrames==expected.samplePeakFrames,"independent Current exact scalar and every-channel peak/location agreement");
 expect(s.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass,"independent Current passes");
 auto history=decodeAppliedExplanation(encodeAppliedExplanation(*s.lastAppliedExplanation),job);expect(history&&history->method==eventSearchVersion&&history->regionSettings->digest()==settings.digest()&&history->narrative.contains("Selected passing alternative"),"attachment3 keeps selection explanation and independent definition");
 const auto path=root.getChildFile("session.levelscope");expect(sf::save(sf::capture(s),path,job).succeeded(),"schema10 session save");auto read=sf::read(path,job);expect(read.document&&read.document->loadedSchema==10,"explicit document10/programme9 pairing");auto opened=sf::verify(*read.document,job);expect(opened.prepared&&!opened.prepared->model.sectionProposal&&opened.prepared->model.lastAppliedExplanation,"reopen retains history not private Apply authority");
 auto oldPath=modifiedCore(path,root.getChildFile("bad-schema9.levelscope"),[](juce::var& v){v.getDynamicObject()->setProperty("schema",9);});expect(!sf::read(oldPath,job).succeeded(),"schema9 cannot claim new event settings");
 auto destination=root.getChildFile("export");destination.createDirectory();auto rendered=ex::render(ex::capture(s),destination,job);expect(rendered.value&&ex::artifactEvaluation(**rendered.value,ex::capture(s)).verdict==Verdict::pass&&(*rendered.value)->analysis().context.processing.staticGainDb==0,"finalized PCM export independently measured without second gain");
 const auto noop=generateSectionProposal(s.sectionProposalRequest(),job);expect(noop.value&&(*noop.value)->status()==SectionProposalStatus::noChange&&(*noop.value)->passes()==0,"already passing all requirements remains no-op");
 JobControl early;early.cancel=true;expect(!generateSectionProposal(request,early).value,"pre-cancelled proposal cannot authorize");JobControl cancel;std::atomic<bool> done{false};std::thread interrupter([&]{while(!done){if(cancel.progress>=.15f){cancel.cancel=true;break;}std::this_thread::sleep_for(std::chrono::milliseconds(1));}});auto cancelled=generateSectionProposal(request,cancel);done=true;interrupter.join();expect(!cancelled.value&&cancelled.error==ErrorCode::cancelled,"cancel auxiliary collection yields no partial authority");
 if(channels==2&&!heterogeneous&&!fallback)
 {
  JobControl late;std::atomic<bool> ended{false};std::thread stop([&]{while(!ended){if(late.progress>=.8f){late.cancel=true;break;}std::this_thread::sleep_for(std::chrono::milliseconds(1));}});auto stopped=generateSectionProposal(request,late);ended=true;stop.join();expect(!stopped.value&&stopped.error==ErrorCode::cancelled,"cancel third assessment discards earlier passing Apply authority");
  auto missing=request;const auto file=original.sourceDocument.file;auto moved=root.getChildFile("owned-temporarily-missing.wav");expect(file.moveFileTo(moved),"move only owned fixture");const auto staleError=validateEventRegionProposalForApply(missing,*proposed);const bool restored=moved.moveFileTo(file);expect(restored&&staleError.isNotEmpty(),"missing source cannot use old measured receipt");
 }
 auto row=juce::var(new juce::DynamicObject());auto* o=row.getDynamicObject();o->setProperty("channels",channels);o->setProperty("passes",proposed->passes());o->setProperty("proxies",int(review->proxyEvaluations));o->setProperty("feature_high_water",int64_t(review->featureHighWaterBytes));o->setProperty("planning_bytes",int64_t(review->planningBytes));o->setProperty("planning_seconds",review->planningSeconds);o->setProperty("extraction_seconds",review->sidecarSeconds);o->setProperty("outcome",proposed->outcomeTitle());cases.add(row);
}
void population()
{
    auto p=unity(10800);auto a=std::make_shared<ProgramAnalysis>();a->snapshotIdentity="stress";a->sampleRate=p.sampleRate;a->frames=p.frames;a->p10=a->p95=-40.691;
    PredictionSidecar sidecar;sidecar.sampleRate=p.sampleRate;sidecar.frames=p.frames;sidecar.channels=2;sidecar.roles={measurement::Role::left,measurement::Role::right};sidecar.baselineSnapshot="stress";sidecar.baselinePlan=regionPlanIdentity(p,0);
    const int64_t step=p.sampleRate/10;sidecar.bins.reserve(108015);
    for(int64_t begin=0;begin<p.frames+p.sampleRate*3/2;begin+=step){PredictionBin b;b.begin=begin;b.end=begin+step;b.drain=begin>=p.frames;b.weightedEnergy=double(step)*1e-4;sidecar.bins.push_back(b);}
    for(int64_t end=p.sampleRate*4/10;end<=p.frames;end+=step){DiagnosticWindow w;w.begin=end-p.sampleRate*4/10;w.end=end;w.energy=1e-4;w.included=true;w.lufs=-40.691;a->momentary.push_back(w);}
    for(int64_t end=p.sampleRate*3;end<=p.frames;end+=step){DiagnosticWindow w;w.begin=end-p.sampleRate*3;w.end=end;w.energy=1e-4;w.included=true;w.lufs=-40.691;a->shortTerm.push_back(w);a->lraWindows.push_back(w);}
    RegionConstraints c;c.frames=p.frames;c.sampleRate=p.sampleRate;c.settings=RegionPlanningSettings::eventDefaults(RegionMode::automatic,p.inputIdentity);c.settings.eventFeatures.hopMs=5;c.identity="combined-resource";RegionFormation formation;formation.nominations.resize(8);
    for(int n=0;n<8;++n){PlanningRoot root;root.hint={int64_t(100+n*1200)*48000,int64_t(1100+n*1200)*48000};root.allowed=root.hint;c.roots.push_back(root);for(int k=0;k<12;++k)formation.nominations[n].push_back({root.hint.begin+int64_t(k+1)*70*48000,1,"bounded nomination"});}
    EventFeatureCollector featureReservation(c,formation,*a,{measurement::Role::left,measurement::Role::right},5,"combined","source",regionPlanIdentity(p,0));
    const size_t scratch=sidecar.bins.size()*2*sizeof(double)+(a->momentary.size()+a->shortTerm.size()+a->lraWindows.size())*sizeof(measurement::Window)+a->lraWindows.size()*2*sizeof(double)+1024*1024;
    const size_t concurrent=a->storageBytes()+sidecar.storageBytes()+featureReservation.reservedBytes()+scratch;
    expect(concurrent<64*1024*1024,"three-hour full populations, proxy scratch and worst 5ms feature reservation fit concurrently");
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
    auto row=juce::var(new juce::DynamicObject());row.getDynamicObject()->setProperty("synthetic_population_seconds",10800);row.getDynamicObject()->setProperty("64_proxy_seconds",seconds);row.getDynamicObject()->setProperty("analysis_sidecar_bytes",int64_t(a->storageBytes()+sidecar.storageBytes()));row.getDynamicObject()->setProperty("concurrent_reserved_and_scratch_bytes",int64_t(concurrent));row.getDynamicObject()->setProperty("maximum_feature_reservation",int64_t(featureReservation.reservedBytes()));cases.add(row);
}

void external(const juce::File& root,const juce::File& file,int experiment)
{
 expect(root.createDirectory().wasOk(),"new read-only media comparison root");JobControl job;auto media=validateMedia(file,job);expect(media.succeeded()&&media.media.frames<=int64_t(media.media.sampleRate)*600,"bounded existing programme, no full film");
 StandaloneSessionModel s;expect(s.setValidatedSource(file,media.media),"adopt verified existing input");s.sourceDocument.layout=layoutFor(media.media.channels);s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(Applicability::preProducedLinearProgramme);measuredComposition(s,MeasurementDomain::source);s.synchronizeSections();
 auto requirements=s.jobRequirements;requirements.maximumLraTenths=media.media.channels==6?150:media.media.frames==28800000?67:180;requirements.provenance="Predeclared development comparison, not delivery default";juce::String error;expect(s.setJobRequirements(requirements,error),"external requirement");
 auto settings=experiment==0?RegionPlanningSettings::defaults(RegionMode::automatic,sectionInputIdentity(s.measurementContext(MeasurementDomain::source))):RegionPlanningSettings::eventDefaults(RegionMode::automatic,sectionInputIdentity(s.measurementContext(MeasurementDomain::source)));
 if(experiment>=2){settings.allowAutomaticMargins=true;settings.leftMarginMs=settings.rightMarginMs=1000;}expect(s.setRegionPlanning(settings,error),"explicit version and margins");auto policy=s.proposalPolicy;policy.durationThousandths=experiment==3?400:experiment>=4?600:200;expect(s.setProposalPolicy(policy,error),"explicit duration comparison");
 auto before=s.currentState.analysis;auto result=generateSectionProposal(s.sectionProposalRequest(),job);if(!result.value)std::cerr<<result.message<<std::endl;expect(bool(result.value),"bounded external comparison completes safely");auto proposed=*result.value;root.getChildFile("proposal.txt").replaceWithText(proposed->text());expect(s.currentState.analysis==before&&!s.currentStateNeedsRemeasurement(),"external comparison leaves Current untouched");auto review=proposed->regionReview();expect(review&&proposed->passes()<=3&&review->proxyEvaluations<=64&&review->planningBytes<=64*1024*1024,"external shared budgets");
 if(proposed->status()==SectionProposalStatus::ready)expect(validateRegionProposalForApply(s.sectionProposalRequest(),*proposed).isEmpty(),"any ready external result has exact measured authority");
 auto row=juce::var(new juce::DynamicObject());auto* o=row.getDynamicObject();o->setProperty("content_sha256",media.media.sha256);o->setProperty("frames",media.media.frames);o->setProperty("experiment",experiment);o->setProperty("outcome",proposed->outcomeTitle());o->setProperty("passes",proposed->passes());o->setProperty("proxies",int(review->proxyEvaluations));o->setProperty("planning_seconds",review->planningSeconds);o->setProperty("extraction_seconds",review->sidecarSeconds);o->setProperty("planning_bytes",int64_t(review->planningBytes));o->setProperty("feature_high_water",int64_t(review->featureHighWaterBytes));cases.add(row);
}

}
int main(int argc,char** argv)
{
 juce::ScopedJuceInitialiser_GUI gui;juce::String group,root,output,input;int experiment=0;bool passed=false;juce::String failure;
 for(int n=1;n<argc;++n){juce::String a(argv[n]);if(a=="--group"&&n+1<argc)group=argv[++n];else if(a=="--root"&&n+1<argc)root=argv[++n];else if(a=="--output"&&n+1<argc)output=argv[++n];else if(a=="--input"&&n+1<argc)input=argv[++n];else if(a=="--experiment"&&n+1<argc)experiment=juce::String(argv[++n]).getIntValue();else return 2;}
 try{if(group=="core")core();else if(group=="features")features();else if(group=="stereo")workflow(juce::File(root),2,false);else if(group=="surround")workflow(juce::File(root),6,false);else if(group=="heterogeneous")workflow(juce::File(root),2,true);else if(group=="fallback")workflow(juce::File(root),2,false,true);else if(group=="population")population();else if(group=="external")external(juce::File(root),juce::File(input),experiment);else throw std::runtime_error("Explicit required event test group missing");passed=true;}catch(const std::exception& e){failure=e.what();}
 auto v=juce::var(new juce::DynamicObject());auto* o=v.getDynamicObject();o->setProperty("suite","M20 complete-event Logic checks");o->setProperty("group",group);o->setProperty("passed",passed);o->setProperty("failure",failure);o->setProperty("checks",checks);o->setProperty("cases",cases);o->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);o->setProperty("devices_opened",false);o->setProperty("combined_qualified",false);if(output.isNotEmpty())juce::File(output).replaceWithText(juce::JSON::toString(v,true));std::cout<<(passed?"PASS ":"FAIL ")<<group<<" "<<failure<<std::endl;return passed?0:1;
}

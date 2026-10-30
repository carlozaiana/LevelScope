#pragma once
#include "Standalone/StandaloneProposalPreferences.h"
#include <sys/resource.h>
namespace
{
juce::StringArray policyGroups;
juce::var policyResources(new juce::DynamicObject());
void m16PolicyValidation()
{
    using P=proposal::ProposalPolicy;P p;
    expect(!p.regionLimit&&p.durationThousandths==200&&p.additionalTenths==30&&p.cumulativeTenths==120&&p.reserveTenths==5,"M16 approved defaults: unrestricted count, preservation duration, independent gain budgets");
    expect(P::decode(p.encode()).digest()==p.digest(),"policy round trip deterministic digest");
    for(int n:{1,65,10000}){auto q=p;q.regionLimit=n;expect(q.validationError().isEmpty()&&q.digest()!=p.digest(),"optional count within structural capacity incl above64");}
    for(int n:{0,-1,10001}){auto q=p;q.regionLimit=n;expect(q.validationError().isNotEmpty(),"invalid explicit count cannot mean unrestricted");}
    for(int field=0;field<4;++field)
    {
        const int maxima[]={1000,60,120,10};
        for(int value:{-1,0,maxima[field],maxima[field]+1})
        {
            auto q=p;q.reserveTenths=0;int* fields[]={&q.durationThousandths,&q.additionalTenths,&q.cumulativeTenths,&q.reserveTenths};*fields[field]=value;
            expect(q.validationError().isEmpty()==(value>=0&&value<=maxima[field]),"exact policy control boundaries");
        }
    }
    auto conflict=p;conflict.additionalTenths=0;expect(conflict.validationError().isNotEmpty(),"reserve must fit additional cap without auto clamp");
    for(int mode=0;mode<5;++mode)
    {
        auto v=p.encode();auto* o=v.getDynamicObject();
        if(mode==0)o->setProperty("version","future");if(mode==1)o->setProperty("extra",0);
        if(mode==2)o->setProperty("additional_tenths_db",3.5);if(mode==3)o->removeProperty("region_limit");
        if(mode==4)o->setProperty("duration_thousandths",1001);
        auto recovered=P::decode(v);expect(recovered.validationError().isNotEmpty(),"unknown/missing/invalid policy unavailable, no default pass");expect(P::decode(recovered.encode()).digest()==recovered.digest(),"unavailable policy identity stable across repeated saves");
    }
    for(auto value:{"NaN","inf","-1","3.01","1e1",".5","1,0",""})expect(!ProposalPreferencesPanel::exactUnits(value,true),"UI rejects ambiguous/nonfinite/overprecision without rounding");
    expect(ProposalPreferencesPanel::exactUnits("20.1",true)==201&&ProposalPreferencesPanel::exactUnits("10000",false)==10000,"UI parses exact units");
    policyGroups.add("policy-validation");
}
void m16SupportAndLimits()
{
    for(int rate:{44100,48000,96000,192000})
    {
        processing::SectionGainPlan a;a.inputIdentity=juce::String::repeatedString("a",64);a.sampleRate=rate;a.frames=rate;
        a.intervals={{0,rate/5,0},{rate/5,rate*2/5,-1},{rate*2/5,rate,0}};
        auto b=a;b.intervals[1].db=-2;
        const auto spans=proposal::changedEnvelopeSupport(a,b);bool contains=true;int64_t count=0;
        // Independent sample oracle: conservative support must contain all unequal
        // samples, plus precisely the two equal ramp endpoints in this fixture.
        for(int64_t frame=0;frame<a.frames;++frame)
        {
            bool in=false;for(auto span:spans)in|=frame>=span.first&&frame<span.second;
            const bool differs=oracleGain(a,frame)!=oracleGain(b,frame);contains&=!differs||in;if(differs)++count;
        }
        expect(contains&&spans.size()==1&&proposal::changedEnvelopeExtent(a,b)==count+2,"support oracle contains all actual changes and disclosed conservative endpoints at every rate");
        auto split=b;split.intervals={{0,rate/10,0},{rate/10,rate/5,0},{rate/5,rate*3/10,-2},{rate*3/10,rate*2/5,-2},{rate*2/5,rate,0}};
        expect(proposal::changedEnvelopeSupport(a,split)==spans&&proposal::changedEnvelopeSupport(b,split).empty(),"redundant markers leave canonical processing/count/support identical");
    }
    auto m=mapFor(48000,4800000,1000);auto r=scalar(-24,1.2);auto peaks=peakMap(m,1.2,0);
    for(int n=0;n<130;n+=2)peaks.sections[size_t(n)].truePeak[3]={std::pow(10.,1.2/20.),m.sections[size_t(n)].begin*8,false};
    proposal::ProposalPolicy p;p.durationThousandths=1000;
    auto result=proposal::selectSectionTrims(r,m,{},peaks,p);
    // Rows are 0.1 s: two ramps touch across unchanged plateau, so choose wider rows.
    m=mapFor(48000,48000000,1000);peaks=peakMap(m,1.2,0);
    for(int n=0;n<130;n+=2)peaks.sections[size_t(n)].truePeak[3]={std::pow(10.,1.2/20.),m.sections[size_t(n)].begin*8,false};
    result=proposal::selectSectionTrims(r,m,{},peaks,p);
    expect(result.available&&result.processingRegions==65&&result.corrections.size()==65,"unrestricted default permits 65 distinct processing regions");
    p.regionLimit=64;auto limited=proposal::selectSectionTrims(r,m,{},peaks,p);
    expect(!limited.available&&limited.corrections.size()==65&&limited.reason.contains("Complete set retained"),"optional64 count rejects complete set without truncation");
    p.regionLimit=65;expect(proposal::selectSectionTrims(r,m,{},peaks,p).available,"exact optional count boundary passes");
    auto broad=mapFor(48000,4800000,5);broad.sections[3].end=4297419;broad.sections[4].begin=4297419;
    auto wide=proposal::selectSectionTrims(r,broad,{},peakMap(broad,1.2,4));
    expect(!wide.available&&wide.reason.contains("preservation preference")&&wide.affectedFrames>broad.frames*89/100,"cinema-like four-row 89.53 percent remains blocked by duration default, not count");
    p={};p.durationThousandths=1000;expect(proposal::selectSectionTrims(r,broad,{},peakMap(broad,1.2,4),p).available,"explicit duration relaxation permits assessment but is not a measured verdict");
    processing::SectionGainPlan existing;existing.inputIdentity=m.inputIdentity;existing.sampleRate=m.sampleRate;existing.frames=m.frames;existing.intervals={{0,m.frames,-2.5}};
    p.cumulativeTenths=30;auto capped=proposal::selectSectionTrims(r,m,existing,peakMap(m,1.2),p);
    expect(capped.available&&capped.corrections.front().resultingDb==-3&&capped.corrections.front().reductionDb==.5&&capped.corrections.front().uncoveredDb>0,"cumulative cap reduces headroom independently from per-request cap");
    existing.intervals[0].db=std::nextafter(-3.,-INFINITY);expect(!proposal::selectSectionTrims(r,m,existing,peakMap(m,1.2),p).available,"immediately over cumulative baseline rejects local plan without modifying settings");
    p.cumulativeTenths=120;p.additionalTenths=60;p.reserveTenths=10;existing.intervals[0].db=2.1234567890123;
    auto full=proposal::selectSectionTrims(r,m,existing,peakMap(m,1.2),p);
    expect(full.available&&full.corrections.front().reserveDb>.999999999999&&full.corrections.front().reductionDb<=6,"positive manual gain and full precision reserve retain distinct meanings");
    p.additionalTenths=0;p.reserveTenths=0;expect(!proposal::selectSectionTrims(r,m,existing,peakMap(m,1.2),p).available,"zero additional disables local proposal");
    policyGroups.add("support-and-gain-limits");
}
void m16Capacity()
{
    struct rusage before{},after{};getrusage(RUSAGE_SELF,&before);
    auto map=mapFor(48000,4800000,10000);auto r=scalar(-24,1.2);auto peaks=peakMap(map,1.2,10000);proposal::ProposalPolicy p;p.durationThousandths=1000;
    const auto start=juce::Time::getMillisecondCounterHiRes();auto seed=proposal::selectSectionTrims(r,map,{},peaks,p);
    const auto elapsed=(juce::Time::getMillisecondCounterHiRes()-start)/1000.;
    struct rusage selected{};getrusage(RUSAGE_SELF,&selected);
    expect(seed.available&&seed.corrections.size()==10000&&seed.processingRegions==1,"10000 changed annotations canonicalize to one connected region without truncation");
    const auto textStart=juce::Time::getMillisecondCounterHiRes();const auto text=proposal::sectionSeedText(seed);
    const auto textSeconds=(juce::Time::getMillisecondCounterHiRes()-textStart)/1000.;
    policyResources.getDynamicObject()->setProperty("complete_text_seconds",textSeconds);expect(textSeconds<=1,"complete 10000-row text fits cancellation boundary budget");
    ProposalInspectionPanel view(text);
    expect(view.pageCount()>1&&view.allText()==text&&text.contains("Section 10000")&&!text.contains("omitted")&&!text.contains("rejected"),"every accepted correction accessible in exact paginated inspection, no rejection-only cutoff");
    getrusage(RUSAGE_SELF,&after);const auto bytes=peaks.storageBytes()+seed.corrections.capacity()*sizeof(proposal::SectionCorrection)+seed.plan.intervals.capacity()*sizeof(processing::GainInterval)+text.getNumBytesAsUTF8();
    policyResources.getDynamicObject()->setProperty("selection_seconds",elapsed);policyResources.getDynamicObject()->setProperty("selection_peak_rss_delta",double(selected.ru_maxrss-before.ru_maxrss));policyResources.getDynamicObject()->setProperty("inspection_peak_rss_delta",double(after.ru_maxrss-before.ru_maxrss));policyResources.getDynamicObject()->setProperty("retained_plan_attribution_text_bytes",double(bytes));
    std::cerr << "M16 capacity: seconds=" << elapsed << " selection RSS=" << selected.ru_maxrss-before.ru_maxrss << " inspection RSS=" << after.ru_maxrss-before.ru_maxrss << " retained=" << bytes << '\n';
    expect(elapsed<=.5&&selected.ru_maxrss-before.ru_maxrss<=32*1024*1024&&bytes<16*1024*1024&&after.ru_maxrss-before.ru_maxrss<=256*1024*1024,"structural-capacity support time/transient/retained memory budgets");
    auto tooMany=mapFor(48000,4800000,10001);expect(!proposal::selectSectionTrims(r,tooMany,{},peakMap(tooMany,1.2,10001),p).available,"10001 interval structure rejected despite unrestricted preference");
    auto doc=StandaloneSessionModel();auto report=ex::capture(doc);report.proposalInspection=text;report.proposalPolicy=p;
    const auto reportStart=juce::Time::getMillisecondCounterHiRes();
    const auto json=ex::reportJson(ex::report(report)),html=ex::reportHtml(ex::report(report));
    const auto reportSeconds=(juce::Time::getMillisecondCounterHiRes()-reportStart)/1000.;expect(reportSeconds<=5,"complete large JSON/HTML generation within five seconds");policyResources.getDynamicObject()->setProperty("large_report_seconds",reportSeconds);
    expect(json.contains("Section 10000")&&html.contains("Section 10000")&&json.getNumBytesAsUTF8()+html.getNumBytesAsUTF8()<=16*1024*1024,"complete large proposal available in bounded JSON/HTML reports");
    policyResources.getDynamicObject()->setProperty("ten_thousand_selection_seconds",elapsed);
    policyResources.getDynamicObject()->setProperty("peak_rss_delta_bytes",double(after.ru_maxrss-before.ru_maxrss));
    policyResources.getDynamicObject()->setProperty("retained_bytes",double(bytes));policyResources.getDynamicObject()->setProperty("report_pair_bytes",double(json.getNumBytesAsUTF8()+html.getNumBytesAsUTF8()));
    auto escapedReport=report;escapedReport.programmeId=juce::String::fromUTF8("漢é🙂");escapedReport.proposalInspection=juce::String::fromUTF8("<script>\"&'漢é🙂</script>");escapedReport.programmeId=escapedReport.proposalInspection;const auto escapedHtml=ex::reportHtml(ex::report(escapedReport));
    std::cerr << "escaping probes=" << !escapedHtml.contains("<script>") << escapedHtml.contains("&lt;script&gt;") << escapedHtml.contains("&amp;") << escapedHtml.contains("&#39;") << escapedHtml.contains(juce::String::fromUTF8("漢")) << '\n';
    expect(!escapedHtml.contains("<script>")&&escapedHtml.contains("&lt;script&gt;")&&escapedHtml.contains("&amp;")&&escapedHtml.contains("&#39;")&&escapedHtml.contains(juce::String::fromUTF8("漢")),"linear report escaping preserves UTF8 and HTML safety");
    juce::String unicode=juce::String::repeatedString("x",11999)+juce::String::fromUTF8("漢é🙂")+"end";ProposalInspectionPanel utfPages(unicode);expect(utfPages.allText()==unicode,"page boundaries preserve UTF8 even without newline");
    policyGroups.add("capacity-and-complete-inspection");
}
void m16UiControls()
{
    StandaloneSessionModel s;const auto revision=s.revision;const auto doc=s.documentRevision;
    ProposalPreferencesPanel panel(s,[]{});
    auto button=[&](const juce::String& name)->juce::Button*{for(auto* c:panel.getChildren())if(auto* b=dynamic_cast<juce::Button*>(c))if(b->getButtonText()==name)return b;return nullptr;};
    auto field=[&](const juce::String& prefix)->juce::TextEditor*{for(auto* c:panel.getChildren())if(auto* e=dynamic_cast<juce::TextEditor*>(c))if(e->getName().startsWith(prefix))return e;return nullptr;};
    auto click=[&](const char* name){auto* b=button(name);expect(b!=nullptr,"policy UI exposes named action");b->triggerClick();juce::MessageManager::getInstance()->runDispatchLoopUntil(20);};
    auto* duration=field("Newly affected");expect(duration!=nullptr,"policy UI exposes duration units");duration->setText("89.6");click("Cancel edits");expect(s.documentRevision==doc&&s.proposalPolicy.durationThousandths==200,"UI Cancel preserves active preferences");
    duration->setText("100.1");click("Use these limits");expect(s.documentRevision==doc,"UI rejects outside duration without clamp");
    duration->setText("100");click("Use these limits");expect(s.proposalPolicy.durationThousandths==1000&&s.revision==revision&&s.documentRevision==doc+1,"UI explicit confirmation changes preferences only");
    click("No region-count restriction");auto* count=field("Optional maximum");expect(count&&count->isEnabled(),"optional explicit count field enabled");count->setText("65");click("Use these limits");expect(s.proposalPolicy.regionLimit==65,"UI supports count greater than64");
    click("Draft defaults");expect(s.proposalPolicy.regionLimit==65,"default draft requires confirmation");click("Use these limits");expect(!s.proposalPolicy.regionLimit&&s.proposalPolicy.durationThousandths==200,"UI confirmed defaults remove count restriction");
    struct Heartbeat:juce::Timer{std::vector<double> ticks;double last=juce::Time::getMillisecondCounterHiRes();void timerCallback()override{auto now=juce::Time::getMillisecondCounterHiRes();ticks.push_back(now-last);last=now;}} heartbeat;
    juce::MemoryOutputStream largeStream;for(int n=0;n<10000;++n)largeStream << "Changed annotation " << n << "; explicit proposal details; no omitted rows.\n";
    const auto large=largeStream.toUTF8();
    heartbeat.startTimer(5);juce::MessageManager::getInstance()->runDispatchLoopUntil(40);
    ProposalInspectionPanel inspect(large);for(int n=0;n<inspect.pageCount();++n){inspect.nextPage();juce::MessageManager::getInstance()->runDispatchLoopUntil(5);}
    juce::MessageManager::getInstance()->runDispatchLoopUntil(100);heartbeat.stopTimer();auto samples=heartbeat.ticks;std::sort(samples.begin(),samples.end());
    expect(!samples.empty()&&samples[size_t(.95*(samples.size()-1))]<=100,"10000 row paginated UI heartbeat p95 within100ms");
    policyResources.getDynamicObject()->setProperty("inspection_heartbeat_p95_ms",samples[size_t(.95*(samples.size()-1))]);policyResources.getDynamicObject()->setProperty("inspection_heartbeat_max_ms",samples.back());
    policyGroups.add("ui-controls-and-responsiveness");
}
void rewritePolicyCore(const juce::File& from,const juce::File& to,const std::function<void(juce::var&)>& edit)
{
    juce::MemoryBlock b;expect(from.loadFileAsData(b),"read owned session");auto* data=static_cast<const uint8_t*>(b.getData());uint64_t length=0;for(int i=0;i<8;++i)length|=uint64_t(data[8+i])<<(8*i);
    auto v=juce::JSON::parse(juce::String::fromUTF8(reinterpret_cast<const char*>(data+80),int(length)));edit(v);
    const auto json=juce::JSON::toString(v,true).toStdString();std::string output="LSSESS01";
    for(int i=0;i<8;++i)output.push_back(char(uint64_t(json.size())>>(8*i)));output+=juce::SHA256(json.data(),json.size()).toHexString().toStdString();output+=json;
    output.append(reinterpret_cast<const char*>(data+80+length),b.getSize()-80-size_t(length));expect(to.replaceWithData(output.data(),output.size()),"write independent core mutation fixture");
}
void m16IdentitySessions(const juce::File& root)
{
    auto s=proposalFixture(root.getChildFile("source"),48000,2,0);JobControl job;
    const auto source=s.source.analysis,current=s.currentState.analysis;const auto audio=OfflineAnalysisService::identityFor(s.measurementContext(MeasurementDomain::currentState),MeasurementDomain::currentState,s.sourceDocument.media.sha256);
    JobControl cancelled;cancelled.cancel=true;const auto cancelStart=juce::Time::getMillisecondCounterHiRes();auto stopped=generateSectionProposal(s.sectionProposalRequest(),cancelled);
    const auto cancelSeconds=(juce::Time::getMillisecondCounterHiRes()-cancelStart)/1000.;expect(stopped.error==ErrorCode::cancelled&&!stopped.value&&cancelSeconds<=1,"cancelled policy request publishes nothing within one second");policyResources.getDynamicObject()->setProperty("cancel_seconds",cancelSeconds);
    auto generated=generateSectionProposal(s.sectionProposalRequest(),job);expect(generated.value&&(*generated.value)->status()==SectionProposalStatus::ready,"fresh exact candidate under default policy");
    expect(s.acceptSectionProposal(*generated.value),"accept bound candidate");const auto revision=s.revision;const auto doc=s.documentRevision;
    auto oldGlobal=proposal::draftGlobalTrim(s.currentTrimRequest());
    auto p=s.proposalPolicy;p.durationThousandths=1000;juce::String error;
    expect(s.setProposalPolicy(p,error)&&s.revision==revision&&s.documentRevision==doc+1&&s.source.analysis==source&&s.currentState.analysis==current,"policy edit only document revision, measurements unchanged");
    expect(!s.acceptTrimProposal(oldGlobal)&&!s.canApplySectionProposal()&&!s.acceptSectionProposal(*generated.value)&&s.sectionProposalDisplayText().contains("STALE"),"old policy result and late completion cannot Apply or replace current result");
    expect(s.setProposalPolicy(p,error)&&s.documentRevision==doc+1,"identical preference confirmation is metadata no-op");
    const auto saved=sf::capture(s);auto file=root.getChildFile("preferences.levelscope");expect(sf::save(saved,file,job).succeeded(),"stale transient result does not block coherent settings save");
    p.additionalTenths=40;expect(s.setProposalPolicy(p,error)&&sf::saved(s,saved,file)&&s.documentDirty(),"policy edit during save remains unsaved");
    auto read=sf::read(file,job);expect(read.document&&read.document->proposalPolicy.digest()==saved.proposalPolicy.digest(),"exact saved policy snapshot, no mixed revisions");
    auto verified=sf::verify(*read.document,job);StandaloneSessionModel restored;expect(verified.document&&sf::install(*verified.document,restored,error),"policy session verifies and reopens safely");
    expect(restored.source.hasMeasurement&&restored.currentState.hasMeasurement&&!restored.sectionProposal&&OfflineAnalysisService::identityFor(restored.measurementContext(MeasurementDomain::currentState),MeasurementDomain::currentState,restored.sourceDocument.media.sha256)==audio,"validated caches reused and exact applied audio retained, proposals not cached");
    for(int route:restored.sourceDocument.routing.fileToDevice)expect(route==-1,"reopened routes muted");
    for(int schema=1;schema<=5;++schema)
    {
        auto old=root.getChildFile("schema-"+juce::String(schema)+".levelscope");rewritePolicyCore(file,old,[schema](juce::var& v){v.getDynamicObject()->setProperty("schema",schema);v["programme"].getDynamicObject()->setProperty("schema",schema);v.getDynamicObject()->removeProperty("proposalPolicy");v.getDynamicObject()->removeProperty("lastAppliedProposalPolicy");v.getDynamicObject()->removeProperty("lastAppliedProposalBinding");});
        const auto hash=juce::SHA256(old).toHexString();auto older=sf::read(old,job);
        expect(older.document&&older.document->proposalPolicy.digest()==proposal::ProposalPolicy{}.digest()&&older.document->notices.joinIntoString(" ").contains("default proposal preferences")&&juce::SHA256(old).toHexString()==hash,"schemas1-5 adopt new defaults in memory; no legacy executable mode or file rewrite");
    }
    auto invalid=root.getChildFile("invalid-policy.levelscope");rewritePolicyCore(file,invalid,[](juce::var& v){v["proposalPolicy"].getDynamicObject()->setProperty("version","future");});
    auto missingPolicy=root.getChildFile("missing-policy.levelscope");rewritePolicyCore(file,missingPolicy,[](juce::var& v){v.getDynamicObject()->removeProperty("proposalPolicy");});auto missingRead=sf::read(missingPolicy,job);expect(missingRead.document&&missingRead.document->proposalPolicy.validationError().isNotEmpty(),"missing independently identifiable policy recovers core with proposals blocked");
    auto bad=sf::read(invalid,job);expect(bad.document&&bad.document->proposalPolicy.validationError().isNotEmpty()&&bad.document->current.staticGainDb==saved.current.staticGainDb,"invalid independently identified policy recovers core settings with unavailable policy");
    auto wrong=sf::verify(*bad.document,job);StandaloneSessionModel recovered;expect(wrong.document&&sf::install(*wrong.document,recovered,error),"unavailable policy recovery retains session");
    auto blocked=generateSectionProposal(recovered.sectionProposalRequest(),job);expect(blocked.value&&(*blocked.value)->status()==SectionProposalStatus::unavailable&&(*blocked.value)->passes()==0,"unknown policy cannot generate or pass, even with valid measurements");
    s.discardSectionProposal();auto repeated=generateSectionProposal(s.sectionProposalRequest(),job);expect(repeated.value&&(*repeated.value)->status()==SectionProposalStatus::ready&&(*repeated.value)->passes()<=2,"changed preferences reuse compatible baseline attribution, exact candidate remeasured");
    expect(s.acceptSectionProposal(*repeated.value)&&s.applySectionProposal()&&s.currentStateNeedsRemeasurement()&&s.lastAppliedProposalPolicy->digest()==p.digest(),"explicit Apply retains generating policy, requires independent verification");
    measuredComposition(s,MeasurementDomain::currentState);expect(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"independent Current measurement verifies actual applied outcome");
    auto report=ex::report(ex::capture(s));expect(report["last_applied_proposal_policy_historical"]["additional_tenths_db"]==juce::var(40),"report applied policy distinct from active preferences");
    const auto exportBefore=ex::capture(s);auto artifact=ex::render(exportBefore,root,job);expect(artifact.value.has_value(),"render artifact from independently verified Current under applied policy");
    const auto approval=ex::approvalIdentity(**artifact.value,exportBefore);
    p.additionalTenths=30;expect(s.setProposalPolicy(p,error),"later preference edit");
    const auto exportAfter=ex::capture(s);expect((*artifact.value)->matchesAudio(exportAfter)&&ex::audioIdentity(exportAfter)==ex::audioIdentity(exportBefore)&&ex::approvalIdentity(**artifact.value,exportAfter)!=approval,"preference edit reuses exact rendered audio but invalidates report approval");
    auto publication=ex::publish(**artifact.value,exportAfter,root.getChildFile("stale-approval-must-not-publish"),approval,job);expect(!publication.value&&!root.getChildFile("stale-approval-must-not-publish").exists(),"old metadata approval cannot publish changed policy report");
    const auto appliedFile=root.getChildFile("applied-policy.levelscope");expect(sf::save(sf::capture(s),appliedFile,job).succeeded(),"applied policy receipt saved with exact envelope");auto appliedRead=sf::read(appliedFile,job);expect(appliedRead.document&&appliedRead.document->lastAppliedProposalPolicy&&appliedRead.document->lastAppliedProposalPolicy->additionalTenths==40&&appliedRead.document->proposalPolicy.additionalTenths==30,"roundtrip distinguishes applied historical policy from active policy");
    auto guarded=s;guarded.auditionBusy=[]{return true;};p.durationThousandths=500;expect(!guarded.setProposalPolicy(p,error),"active playback guard prevents policy edits without stopping audio implicitly");p.durationThousandths=1000;report=ex::report(ex::capture(s));expect(report["proposal_preferences"]["additional_tenths_db"]==juce::var(30)&&report["last_applied_proposal_policy_historical"]["additional_tenths_db"]==juce::var(40),"policy edits cannot relabel earlier application");
    const auto passed=s.currentState.analysis;const auto rev=s.revision;auto noChange=generateSectionProposal(s.sectionProposalRequest(),job);expect(noChange.value&&(*noChange.value)->status()==SectionProposalStatus::noChange&&(*noChange.value)->passes()==0&&s.currentState.analysis==passed&&s.revision==rev,"already passing remains true no-op under new policy");
    policyGroups.add("identity-session-report");
}
void m16LargeMeasuredProposal(const juce::File& root)
{
    auto directory=root.getChildFile("large-measured-proposal");directory.createDirectory();constexpr int rate=48000,frames=rate*10;
    juce::AudioBuffer<float> audio(2,frames);
    for(int c=0;c<2;++c)for(int n=0;n<frames;++n)audio.setSample(c,n,float(std::llround(32768*.06*std::sin(2*3.141592653589793*1000*n/rate)))/32768.f);
    for(int n=0;n<130;n+=2)audio.setSample(0,n*480+240,.8f);
    const auto file=directory.getChildFile("sixty-five-peaks.wav");floatFile(file,rate,audio);JobControl job;const auto media=validateMedia(file,job);expect(media.succeeded(),"large proposal fixture verified");
    StandaloneSessionModel s;expect(s.setValidatedSource(file,media.media),"large proposal source imported");s.sourceDocument.layout=layoutFor(2);s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);expect(s.initializeCurrentStateFromSource(),"large proposal Current initialized");
    measuredComposition(s,MeasurementDomain::source);measuredComposition(s,MeasurementDomain::currentState);s.synchronizeSections();auto map=mapFor(rate,frames,1000);map.inputIdentity=s.sectionDocument.map().inputIdentity;
    expect(s.sectionDocument.restore(map,map.inputIdentity),"explicit synthetic marker-dense geometry");s.documentChanged();
    expect(sf::save(sf::capture(s),directory.getChildFile("large-proposal-input.levelscope"),job).succeeded(),"save reproducible silent UI input with exact measured baseline");
    const auto before=s.currentState.analysis;auto generated=generateSectionProposal(s.sectionProposalRequest(),job);
    directory.getChildFile("large-proposal-result.txt").replaceWithText(generated.value?(*generated.value)->text():generated.message);
    if(generated.value)std::cerr << "Large candidate: status=" << int((*generated.value)->status()) << " rows=" << (*generated.value)->selection().corrections.size() << " passes=" << (*generated.value)->passes() << "\n";
    expect(generated.value&&(*generated.value)->status()==SectionProposalStatus::ready&&(*generated.value)->selection().corrections.size()>64&&(*generated.value)->passes()<=3,"actual measured large candidate with more than64 changed annotation rows passes in bounded private passes");
    const auto lastLabel="Section "+juce::String((*generated.value)->selection().corrections.back().id)+" ";
    policyResources.getDynamicObject()->setProperty("large_candidate_changed_rows",int((*generated.value)->selection().corrections.size()));policyResources.getDynamicObject()->setProperty("large_candidate_regions",int((*generated.value)->selection().processingRegions));policyResources.getDynamicObject()->setProperty("large_candidate_affected_frames",juce::String((*generated.value)->selection().affectedFrames));
    expect((*generated.value)->text().startsWith("READY FOR EXPLICIT APPLY")&&(*generated.value)->text().contains(lastLabel)&&(*generated.value)->text().contains("66720 frames / 1.390000 seconds / 13.900000% of programme"),"large accepted outcome leads summary; complete rows and affected frames/seconds/percentage remain inspectable");
    expect(s.acceptSectionProposal(*generated.value)&&s.currentState.analysis==before,"large candidate assessment never becomes Current automatically");
    const auto report=ex::report(ex::capture(s));expect(report["proposal_inspection_not_current_measurement"].toString().contains(lastLabel),"full accepted proposal appears in reporting without replacing Current measurements");
    expect(directory.getChildFile("large-proposal-report.html").replaceWithText(ex::reportHtml(report))&&directory.getChildFile("large-proposal-report.json").replaceWithText(ex::reportJson(report)),"write complete accepted-proposal inspection reports");
    policyGroups.add("large-measured-proposal");
}
void m16MonoPolicy(const juce::File& root)
{
    for(auto role:{juce::AudioChannelSet::left,juce::AudioChannelSet::right,juce::AudioChannelSet::centre})
    {
        const int rate=48000,frames=rate*4;auto directory=root.getChildFile("mono-"+juce::String(int(role)));directory.createDirectory();
        juce::AudioBuffer<float> audio(1,frames);for(int n=0;n<frames;++n)audio.setSample(0,n,float(.085*std::sin(2*3.141592653589793*1000*n/rate)));audio.setSample(0,rate*2+rate/100,.8f);
        const auto file=directory.getChildFile("mono.wav");floatFile(file,rate,audio);JobControl job;const auto media=validateMedia(file,job);StandaloneSessionModel s;expect(media.succeeded()&&s.setValidatedSource(file,media.media),"mono preference fixture verified");
        s.sourceDocument.layout=layoutFor(1);s.sourceDocument.layout.roles[0]=role;s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);expect(s.initializeCurrentStateFromSource(),"mono Current initialized");
        measuredComposition(s,MeasurementDomain::source);measuredComposition(s,MeasurementDomain::currentState);s.synchronizeSections();auto map=mapFor(rate,frames,100);map.inputIdentity=s.sectionDocument.map().inputIdentity;expect(s.sectionDocument.restore(map,map.inputIdentity),"mono explicit sections");s.documentChanged();
        auto candidate=generateSectionProposal(s.sectionProposalRequest(),job);expect(candidate.value&&(*candidate.value)->status()==SectionProposalStatus::ready&&(*candidate.value)->passes()<=3,"confirmed mono L/R/C policy path exact candidate passes");
        expect(s.acceptSectionProposal(*candidate.value)&&s.applySectionProposal()&&s.currentStateNeedsRemeasurement(),"mono proposal Apply is explicit and unverified");measuredComposition(s,MeasurementDomain::currentState);expect(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass&&s.sourceDocument.layout.count==1,"independent mono Current verification with no channel duplication");
    }
    policyGroups.add("mono-policy-path");
}
void m16Preferences(const juce::File& root)
{m16PolicyValidation();m16SupportAndLimits();m16Capacity();m16UiControls();m16IdentitySessions(root);m16LargeMeasuredProposal(root);m16MonoPolicy(root);}
}

#include "CompositionTestFixtures.h"
#include "Standalone/StandaloneJobRequirements.h"
#include <algorithm>
#include <fstream>

namespace
{
using proposal::JobRequirements;
using proposal::Verdict;
juce::StringArray groups;
juce::var resources(new juce::DynamicObject());
proposal::MeasuredProgramme input(proposal::SnapshotDomain domain=proposal::SnapshotDomain::currentState)
{
    proposal::MeasuredProgramme m;m.domain=domain;m.snapshotIdentity="explicit-synthetic-boundary";m.currentIdentityMatches=true;
    m.measurementVersion=measurement::version;m.applicability=proposal::Applicability::preProducedLinearProgramme;
    m.integrated={-23,measurement::Status::available};m.truePeak={-1,measurement::Status::available};m.lra={10,measurement::Status::available};m.lraStable=true;return m;
}
template<class F> void rejects(F f,const char* label){bool failed=false;try{f();}catch(...){failed=true;}expect(failed,label);}
void definitionTests(const juce::File& root)
{
    JobRequirements r;const auto base=proposal::ebuR128ProgrammeQc();
    expect(base.digest()==r.baseDigest&&!base.limits.hasMaxLoudnessRange,"immutable EBU digest and absent universal LRA ceiling");
    expect(r.validationError().isEmpty()&&!r.maximumLraTenths,"default is genuinely unset");
    for(int ceiling:{0,1,100,999,1000}){r.maximumLraTenths=ceiling;expect(JobRequirements::fromJson(juce::JSON::toString(r.encode())).digest()==r.digest(),"exact named definition JSON round trip");}
    for(const auto* s:{"0","0.0","100.0","17.4"})expect(JobRequirements::parseCeiling(s).has_value(),"explicit legal tenths");
    for(const auto* s:{"","-0","-1","NaN","inf","1e1",".5","1,0","10.01","100.1","1001"})expect(!JobRequirements::parseCeiling(s),"invalid value never rounds/clamps");
    for(int n=0;n<10;++n)
    {
        auto bad=r;switch(n){case 0:bad.schema=2;break;case 1:bad.id="unknown";break;case 2:bad.baseId="invented";break;case 3:bad.baseDigest="bad";break;case 4:bad.method="dialogue-gated";break;case 5:bad.applicability="live";break;case 6:bad.maximumLraTenths=-1;break;case 7:bad.maximumLraTenths=1001;break;case 8:bad.name=juce::String::repeatedString("x",81);break;case 9:bad.provenance=juce::String::repeatedString("x",2049);break;}
        expect(bad.validationError().isNotEmpty(),"invalid requirement definition rejected");rejects([&]{JobRequirements::decode(bad.encode());},"invalid serialized definition cannot acquire authority");
    }
    r.name=juce::String::fromUTF8("Mix italiano — <script>alert(1)</script>");r.provenance="User supplied <delivery> & reference";
    auto bad=r.encode();bad.getDynamicObject()->setProperty("maximum_lra_tenths",1.0);rejects([&]{JobRequirements::decode(bad);},"integer units reject floating storage");
    bad=r.encode();bad.getDynamicObject()->setProperty("authoritative",true);rejects([&]{JobRequirements::decode(bad);},"no external authority flag");
    bad=r.encode();bad.getDynamicObject()->setProperty("name","tampered");rejects([&]{JobRequirements::decode(bad);},"digest tampering rejected");
    auto text=juce::JSON::toString(r.encode(),true);rejects([&]{JobRequirements::fromJson(text.replaceFirstOccurrenceOf("{","{\"name\":\"duplicate\","));},"duplicate fields rejected");
    rejects([&]{JobRequirements::fromJson(text+"{}");},"trailing object rejected");
    rejects([&]{JobRequirements::fromJson(text+"garbage");},"trailing non-JSON rejected");
    rejects([&]{JobRequirements::fromJson(juce::String::repeatedString("[",1000));},"deep untrusted JSON rejected before general parser");
    rejects([&]{JobRequirements::fromJson(juce::String::repeatedString(" ",65537));},"oversized JSON rejected");
    const auto file=root.getChildFile("portable.levelscope-target");juce::String error;
    expect(saveRequirementsPreset(r,file,error),"portable preset saved atomically");const auto original=juce::SHA256(file).toHexString();
    for(const char* stage:{"open","write","replace"})
        expect(!saveRequirementsPreset(r,file,error,[&](const char* name){return std::strcmp(name,stage)==0;})&&juce::SHA256(file).toHexString()==original,"write failure preserves previous preset");
    const auto loaded=loadRequirementsPreset(file,error);expect(loaded&&loaded->digest()==r.digest(),"portable preset reloaded independently");
    expect(root.findChildFiles(juce::File::findFiles,false).size()==1,"owned partial preset files cleaned");
    auto duplicate=r;duplicate.id=juce::Uuid().toString();expect(duplicate.name==r.name&&duplicate.digest()!=r.digest(),"duplicate display names do not alias identity");
    groups.add("definitions-presets-failures");
}
void ruleTests()
{
    JobRequirements r;r.maximumLraTenths=100;const auto p=r.applyTo(proposal::ebuR128ProgrammeQc());
    for(auto domain:{proposal::SnapshotDomain::source,proposal::SnapshotDomain::currentState,proposal::SnapshotDomain::renderedOutput,proposal::SnapshotDomain::privateCandidate})
    {
        auto m=input(domain);const auto evaluate=[&]{return proposal::evaluateMeasuredProgramme(p,m);};
        for(double value:{std::nextafter(10.,0.),10.,std::nextafter(10.,20.)})
        {m.lra.value=value;expect(evaluate().verdict==(value<=10?Verdict::pass:Verdict::fail),"full precision inclusive LRA ceiling at adjacent values");expect(juce::String(value,6)=="10.000000","display rounding cannot determine verdict");}
        for(auto value:{-1.,std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
        {m=input(domain);m.lra.value=value;expect(evaluate().verdict==Verdict::indeterminate,"negative/non-finite required LRA never passes");}
        for(auto status:{measurement::Status::silence,measurement::Status::belowGate,measurement::Status::insufficientDuration})
        {m=input(domain);m.lra.status=status;expect(evaluate().verdict==Verdict::indeterminate,"unavailable statuses retain no-pass semantics");}
        m=input(domain);m.lra.value.reset();expect(evaluate().verdict==Verdict::indeterminate,"missing required LRA cannot pass");
        m=input(domain);m.lraStable=false;expect(evaluate().verdict==Verdict::indeterminate,"short-programme available LRA remains uninterpretable for rule");
        expect(proposal::evaluateMeasuredProgramme(JobRequirements{}.applyTo(proposal::ebuR128ProgrammeQc()),m).verdict==Verdict::pass,"unset LRA ignores optional unavailable/unstable metric");
        m=input(domain);m.currentIdentityMatches=false;expect(evaluate().verdict==Verdict::notEvaluated,"stale measured input rejected");
    }
    auto m=input(proposal::SnapshotDomain::predicted);expect(proposal::evaluateMeasuredProgramme(p,m).verdict==Verdict::notEvaluated,"prediction cannot satisfy selected rules");
    r.maximumLraTenths=0;m=input();m.lra.value=0;expect(proposal::evaluateMeasuredProgramme(r.applyTo(proposal::ebuR128ProgrammeQc()),m).verdict==Verdict::pass,"zero is a strict available ceiling");
    m.lra.value=std::nextafter(0.,1.);expect(proposal::evaluateMeasuredProgramme(r.applyTo(proposal::ebuR128ProgrammeQc()),m).verdict==Verdict::fail,"positive value immediately above zero fails");
    r.maximumLraTenths=100;auto current=input(),artifact=input(proposal::SnapshotDomain::renderedOutput);artifact.lra.value=std::nextafter(10.,20.);
    const auto ce=proposal::evaluateMeasuredProgramme(p,current),ae=proposal::evaluateMeasuredProgramme(p,artifact);
    expect(!ex::rulesPermitPublication(ce,ae),"artifact-only LRA boundary failure blocks production publication gate");
    artifact.lra.value=10;expect(ex::rulesPermitPublication(ce,proposal::evaluateMeasuredProgramme(p,artifact)),"both domains pass the same requirements gate");
    r.name="Different requirement provenance";expect(!ex::rulesPermitPublication(ce,proposal::evaluateMeasuredProgramme(r.applyTo(proposal::ebuR128ProgrammeQc()),artifact)),"differing requirements identity cannot share export approval");
    groups.add("measured-rule-boundaries");
}
void proposalTests()
{
    proposal::TrimRequest r;r.source=input(proposal::SnapshotDomain::source);r.current=input();r.sourceInputIdentity=r.currentInputIdentity="same-content-roles-range";
    JobRequirements job;job.maximumLraTenths=100;r.profile=job.applyTo(proposal::ebuR128ProgrammeQc());r.currentGainDb=-3.123456789;
    const auto evaluated=proposal::evaluateMeasuredProgramme(r.profile,r.current);
    expect(evaluated.profileDigest==job.baseDigest&&evaluated.jobRequirementsDigest==job.digest()&&evaluated.rulesDigest==r.profile.digest(),"base and job identities are separately retained in combined evaluation");
    auto p=proposal::draftGlobalTrim(r);expect(p.status==proposal::TrimStatus::noChange&&p.totalGainDb==r.currentGainDb&&p.deltaDb==0,"all selected rules preserve exact existing gain");
    for(int n=0;n<3;++n)
    {r.current=input();if(n==0)r.current.lra.value=11;else if(n==1)r.current.lraStable=false;else r.current.integrated.value=-24;
        p=proposal::draftGlobalTrim(r);expect(p.status==proposal::TrimStatus::unavailable&&p.explanation.contains("not supported yet")&&p.explanation.contains("not a claim that no solution exists"),"unmet LRA or other rule has truthful capability limitation");}
    job.maximumLraTenths.reset();r.profile=job.applyTo(proposal::ebuR128ProgrammeQc());r.current.truePeak.value=-5;
    expect(proposal::draftGlobalTrim(r).status==proposal::TrimStatus::draft,"plain EBU automatic global correction preserved");
    r.currentInputIdentity="different";expect(proposal::draftGlobalTrim(r).status==proposal::TrimStatus::unavailable,"individually valid mismatched inputs rejected");
    groups.add("proposal-capability-noop");
}
StandaloneSessionModel model(const juce::File& folder,int channels,int rate,int frames)
{
    expect(folder.createDirectory().wasOk(),"owned programme fixture folder");auto file=folder.getChildFile("signal.wav");
    juce::WavAudioFormat format;auto stream=file.createOutputStream();expect(bool(stream),"fixture stream");
    std::unique_ptr<juce::AudioFormatWriter> writer(format.createWriterFor(stream.get(),rate,channels,24,{},0));expect(bool(writer),"PCM fixture writer");stream.release();
    juce::AudioBuffer<float> block(channels,4096);bool written=true;
    for(int first=0;first<frames;first+=4096){const int n=std::min(4096,frames-first);for(int c=0;c<channels;++c)for(int k=0;k<n;++k)
        block.setSample(c,k,float((first+k<frames/2?.015:.045)*std::sin((first+k)*2*3.141592653589793*(997+c*101)/rate)));written&=writer->writeFromAudioSampleBuffer(block,0,n);}
    writer.reset();expect(written,"whole fixture written");JobControl job;auto v=validateMedia(file,job);expect(v.succeeded(),"fixture content verified");
    StandaloneSessionModel s;expect(s.setValidatedSource(file,v.media),"adopt source");s.sourceDocument.layout=layoutFor(channels);s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
    measuredComposition(s,MeasurementDomain::source);expect(s.initializeCurrentStateFromSource(),"initialize current");expect(s.setCurrentGainDb(-23-*s.source.integratedLufs),"explicit normalization gain");measuredComposition(s,MeasurementDomain::currentState);return s;
}
juce::File rewriteCore(const juce::File& from,const juce::File& to,const std::function<void(juce::var&)>& change)
{
    juce::MemoryBlock bytes;expect(from.loadFileAsData(bytes)&&bytes.getSize()>=80,"load owned session for mutation fixture");const auto* data=static_cast<const unsigned char*>(bytes.getData());uint64_t size=0;
    for(int i=0;i<8;++i)size|=uint64_t(data[8+i])<<(i*8);
    auto json=juce::JSON::parse(juce::String::fromUTF8(reinterpret_cast<const char*>(data+80),int(size)));change(json);const auto text=juce::JSON::toString(json,true).toStdString();
    auto out=to.createOutputStream();expect(bool(out),"mutated test document stream");out->write("LSSESS01",8);out->writeInt64(juce::int64(text.size()));const auto hash=juce::SHA256(text.data(),text.size()).toHexString();out->write(hash.toRawUTF8(),64);out->write(text.data(),text.size());out->write(data+80+size,bytes.getSize()-80-size);out->flush();expect(out->getStatus().wasOk(),"mutated core checksum finalized");return to;
}
void sessionTests(const juce::File& root,StandaloneSessionModel& s)
{
    root.createDirectory();JobControl control;juce::String error;JobRequirements r;r.name="Named delivery job";r.maximumLraTenths=1000;r.provenance="User-specified test contract";r.id=juce::Uuid().toString();
    const auto source=s.source.analysis,current=s.currentState.analysis;const auto audio=s.revision,document=s.documentRevision;const auto identity=ex::audioIdentity(ex::capture(s));
    expect(s.setJobRequirements(r,error)&&s.revision==audio&&s.documentRevision==document+1&&s.source.analysis==source&&s.currentState.analysis==current,"requirement edit preserves exact measurement and audio identities");
    expect(ex::audioIdentity(ex::capture(s))==identity&&s.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass,"immediate full selected-rule reevaluation without worker");
    const auto revision=s.documentRevision;expect(s.setJobRequirements(r,error)&&s.documentRevision==revision,"identical requirement confirmation is no-op");
    const auto noChange=generateSectionProposal(s.sectionProposalRequest(),control);expect(noChange.value&&(*noChange.value)->status()==SectionProposalStatus::noChange&&(*noChange.value)->passes()==0,"all selected rules section proposal no-op uses zero private passes");
    const auto requestBefore=s.sectionProposalRequest();r.maximumLraTenths=0;++r.revision;expect(s.setJobRequirements(r,error),"select strict explicit requirement");
    expect(s.evaluation(MeasurementDomain::currentState).verdict==Verdict::fail&&!s.acceptSectionProposal(*noChange.value),"old completion rejected and strict rule fails without remeasurement");
    expect(ex::readiness(ex::capture(s)).isNotEmpty(),"LRA failure blocks audio export");
    auto report=ex::report(ex::capture(s));expect(!report.hasProperty("artifact")&&report["current"]["evaluation"]["verdict"].toString()=="FAIL","measurement-only report exposes actual failure without artifact");
    auto unsupported=generateSectionProposal(s.sectionProposalRequest(),control);expect(unsupported.value&&(*unsupported.value)->status()==SectionProposalStatus::unavailable&&(*unsupported.value)->passes()==0,"selected unmet requirements generate no private correction passes");
    r.maximumLraTenths=1000;++r.revision;s.setJobRequirements(r,error);
    auto captured=sf::capture(s);const auto file=root.getChildFile("named.levelscope");expect(sf::save(captured,file,control).succeeded(),"schema 7 coherent session save");
    r.name="Edits during save";++r.revision;s.setJobRequirements(r,error);sf::saved(s,captured,file);expect(s.documentRevision!=s.savedDocumentRevision,"changes after capture remain unsaved");
    const auto reopenStart=juce::Time::getMillisecondCounterHiRes();
    auto read=sf::read(file,control);expect(read.document&&read.document->loadedSchema==7,"schema 7 decoded");auto verified=sf::verify(*read.document,control);StandaloneSessionModel reopened;
    expect(verified.document&&sf::install(*verified.document,reopened,error),"verified named session reopen");
    expect(reopened.jobRequirements.digest()==captured.jobRequirements.digest()&&reopened.currentState.hasMeasurement&&reopened.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass,"independent embedded definition and compatible caches reconstruct evaluation");
    const auto reopenMs=juce::Time::getMillisecondCounterHiRes()-reopenStart;expect(reopenMs<=5000,"short-session verify/reopen within five seconds");resources.getDynamicObject()->setProperty("short_session_reopen_ms",reopenMs);
    expect(reopened.sessionId!=s.sessionId&&!reopened.acceptMeasurement(*current),"old-runtime completion cannot modify newly opened session");
    auto broken=rewriteCore(file,root.getChildFile("invalid-requirements.levelscope"),[](juce::var& core){core["jobRequirements"].getDynamicObject()->setProperty("method","unsupported");});
    auto recover=sf::read(broken,control);expect(recover.document&&recover.document->jobRequirements.validationError().isNotEmpty(),"damaged independent requirements recover core without fallback");
    auto resolved=sf::verify(*recover.document,control);StandaloneSessionModel partial;expect(resolved.document&&sf::install(*resolved.document,partial,error)&&partial.currentState.hasMeasurement,"invalid rule preserves valid settings/media/cache");
    expect(partial.evaluation(MeasurementDomain::currentState).verdict==Verdict::notEvaluated&&ex::readiness(ex::capture(partial)).isNotEmpty(),"recovered invalid requirement cannot pass or publish");
    // Empty core has no unsupported newer processing; prove every old schema defaults unset.
    StandaloneSessionModel empty;auto emptyFile=root.getChildFile("empty.levelscope");expect(sf::save(sf::capture(empty),emptyFile,control).succeeded(),"empty migration fixture");
    for(int schema=1;schema<=6;++schema)
    {
        auto old=rewriteCore(emptyFile,root.getChildFile("schema"+juce::String(schema)+".levelscope"),[&](juce::var& core){core.getDynamicObject()->setProperty("schema",schema);core["programme"].getDynamicObject()->setProperty("schema",schema);core.getDynamicObject()->removeProperty("jobRequirements");});
        auto legacy=sf::read(old,control);expect(legacy.document&&!legacy.document->jobRequirements.maximumLraTenths&&legacy.document->jobRequirements.validationError().isEmpty(),"schemas 1 through 6 migrate with unset LRA");
    }
    auto cacheBroken=root.getChildFile("damaged-cache.levelscope");juce::MemoryBlock cacheData;file.loadFileAsData(cacheData);static_cast<char*>(cacheData.getData())[cacheData.getSize()-1]^=1;cacheBroken.replaceWithData(cacheData.getData(),cacheData.getSize());
    const auto cacheRecovered=sf::read(cacheBroken,control);expect(cacheRecovered.document&&cacheRecovered.document->jobRequirements.digest()==captured.jobRequirements.digest()&&!cacheRecovered.document->currentCache,"damaged independent cache preserves complete requirements and core settings");
    auto malformed=file.getSiblingFile("bad-core.levelscope");malformed.replaceWithText("invalid core");expect(!sf::read(malformed,control).document,"invalid core still fails safely");
    resources.getDynamicObject()->setProperty("session_bytes",juce::int64(file.getSize()));groups.add("session-identity-reuse-migration");
}
void exportTests(const juce::File& root,StandaloneSessionModel& s)
{
    root.createDirectory();juce::String error;JobRequirements r;r.maximumLraTenths=1000;r.name="<script>not executed</script>";s.setJobRequirements(r,error);
    auto request=ex::capture(s);JobControl control;auto result=ex::render(request,root,control);if(!result.value)std::cerr<<result.message<<'\n';expect(bool(result.value),"real PCM render and reopened artifact measurement with required LRA");auto candidate=*result.value;
    expect(ex::artifactEvaluation(*candidate,request).verdict==Verdict::pass&&candidate->analysis().context.processing.staticGainDb==0,"actual artifact meets all rules without a second processing gain");
    const auto approval=ex::approvalIdentity(*candidate,request),hash=candidate->hash();r.name="Changed provenance";++r.revision;s.setJobRequirements(r,error);auto fresh=ex::capture(s);
    expect(candidate->matchesAudio(fresh)&&candidate->hash()==hash&&ex::approvalIdentity(*candidate,fresh)!=approval,"requirement edit reuses exact artifact while refreshing approval");
    expect(!ex::publish(*candidate,fresh,root.getChildFile("stale"),approval,control).value,"stale approval cannot publish");
    auto report=ex::report(request,candidate.get());const auto html=ex::reportHtml(report);
    expect(report["schema"].toString()=="levelscope-measurement-report-v5"&&report["profile"]["digest"].toString()==r.baseDigest,"report versions addition while retaining base EBU digest");
    expect(html.contains("&lt;script&gt;")&&!html.contains("<script>not executed")&&html.contains("required unstable or unavailable LRA cannot pass"),"HTML escapes provenance and distinguishes required LRA");
    expect(bool(ex::publish(*candidate,fresh,root.getChildFile("accepted"),ex::approvalIdentity(*candidate,fresh),control).value),"explicit fresh approval publishes verified PCM/report");
    groups.addIfNotAlreadyThere("real-export-report-gates");
}

void manualAndInputTests(const juce::File& root,StandaloneSessionModel& s)
{
    juce::String error;JobRequirements r;r.maximumLraTenths=1000;s.setJobRequirements(r,error);
    const auto source=s.source.analysis;const auto oldRevision=s.revision;
    expect(s.beginSectionTrimDraft()&&s.setSectionDraftGain(0,-1.123456789),"manual section draft remains available under an LRA requirement");
    expect(s.currentState.hasMeasurement&&s.revision==oldRevision,"unapplied manual draft leaves measured authority unchanged");
    expect(s.applySectionTrimDraft()&&s.currentStateNeedsRemeasurement()&&!s.currentState.hasMeasurement&&s.source.analysis==source,"explicit manual Apply invalidates only affected Current authority");
    JobControl cancel;cancel.cancel=true;const auto context=s.measurementContext(MeasurementDomain::currentState);
    auto abandoned=OfflineAnalysisService().analyze(context,MeasurementDomain::currentState,cancel);
    expect(!abandoned.value&&abandoned.error==ErrorCode::cancelled&&s.currentStateNeedsRemeasurement()&&s.currentStateDocument.sectionGains->intervals.front().db==-1.123456789,"cancelled verification retains exact applied gain and UNVERIFIED state");
    measuredComposition(s,MeasurementDomain::currentState);
    expect(s.evaluation(MeasurementDomain::currentState).verdict==Verdict::fail,"actual measured manual result displays failure against complete rule set");
    expect(s.setCurrentGainDb(s.currentStateDocument.staticGainDb-23-*s.currentState.integratedLufs),"explicit manual global correction");measuredComposition(s,MeasurementDomain::currentState);
    expect(s.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass,"independent Current verification meets all selected rules");
    auto before=s.currentState.analysis;const auto evaluation=s.evaluation(MeasurementDomain::currentState).identity;
    expect(s.setOutputRoute(0,-1),"partial monitoring route changed while stopped");
    s.referenceVideo.startOnAudio=video::Time::parse("-1001/30000");
    expect(s.currentState.analysis==before&&s.evaluation(MeasurementDomain::currentState).identity==evaluation,"monitor route and rational video-offset metadata preserve evaluation without video execution");
    auto policy=s.proposalPolicy;policy.durationThousandths=1000;expect(s.setProposalPolicy(policy,error)&&s.evaluation(MeasurementDomain::currentState).identity==evaluation,"proposal preferences remain separate from selected job requirements");
    s.referenceVideo={};
    // Exercise programme assemblies/compositions with short existing sample paths;
    // stable LRA qualification uses the separate real 65-second fixtures.
    auto fixture=fixtureComposition(root.getChildFile("composed"),48000,CompositionPreset::surroundMonoCentre,true,48000);
    auto composed=compositionModel(fixture.recipe);measuredComposition(composed,MeasurementDomain::source);measuredComposition(composed,MeasurementDomain::currentState);
    expect(composed.setJobRequirements(r,error)&&composed.evaluation(MeasurementDomain::currentState).verdict!=Verdict::pass,"required unstable LRA on composed input cannot pass");
    StandaloneSessionModel assembled;expect(assembled.setValidatedProgramme(*fixture.recipe.me.assembly)&&assembled.initializeCurrentStateFromSource(),"existing six-mono assembly remains available");assembled.setSelectedTargetProfileIndex(1);assembled.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);measuredComposition(assembled,MeasurementDomain::source);measuredComposition(assembled,MeasurementDomain::currentState);assembled.setJobRequirements(r,error);
    expect(assembled.evaluation(MeasurementDomain::currentState).verdict!=Verdict::pass,"assembled-channel workflow evaluates selected LRA without changing programme assembly");
    groups.add("manual-processing-input-identity");
}

void uiAndResources(StandaloneSessionModel& s)
{
    JobRequirementsPanel panel(s,[]{});panel.setSize(1100,88);panel.resized();
    juce::TextEditor* ceiling=nullptr;juce::TextButton* use=nullptr;juce::ToggleButton* enabled=nullptr;
    for(auto* c:panel.getChildren())
    {if(auto* f=dynamic_cast<juce::TextEditor*>(c);f&&f->getName().contains("Draft maximum"))ceiling=f;if(auto* b=dynamic_cast<juce::TextButton*>(c);b&&b->getButtonText()=="Use EBU requirements")use=b;if(auto* t=dynamic_cast<juce::ToggleButton*>(c))enabled=t;}
    expect(ceiling&&use&&enabled,"inline controls expose named fields");enabled->setToggleState(true,juce::dontSendNotification);ceiling->setText("",false);auto before=s.jobRequirements.digest();use->onClick();expect(s.jobRequirements.digest()==before,"empty required ceiling cannot silently become zero");
    ceiling->setText("0.0",false);use->onClick();expect(s.jobRequirements.maximumLraTenths==0,"explicit zero is retained");
    const auto confirmed=s.jobRequirements.digest();const auto revision=s.documentRevision;
    ceiling->onTextChange();use->onClick();bool draftNotice=false;
    for(auto* c:panel.getChildren())if(auto* label=dynamic_cast<juce::Label*>(c))draftNotice|=label->getText().startsWith("Draft only");
    expect(!draftNotice&&s.jobRequirements.digest()==confirmed&&s.documentRevision==revision,"identical confirmation clears draft notice without editing active definition");
    enabled->setToggleState(false,juce::dontSendNotification);use->onClick();expect(!s.jobRequirements.maximumLraTenths,"explicit unset restores ordinary automatic workflow");
    std::vector<double> times;for(int n=0;n<200;++n){auto t=juce::Time::getMillisecondCounterHiRes();s.evaluation(MeasurementDomain::source);s.evaluation(MeasurementDomain::currentState);times.push_back(juce::Time::getMillisecondCounterHiRes()-t);}
    std::sort(times.begin(),times.end());expect(times[190]<=10,"cached reevaluation p95 within 10 ms");resources.getDynamicObject()->setProperty("cached_evaluation_p95_ms",times[190]);
    times.clear();for(int n=0;n<100;++n){auto t=juce::Time::getMillisecondCounterHiRes();panel.refresh();panel.resized();times.push_back(juce::Time::getMillisecondCounterHiRes()-t);}std::sort(times.begin(),times.end());expect(times[95]<=100,"inline feedback p95 within 100 ms");resources.getDynamicObject()->setProperty("inline_feedback_p95_ms",times[95]);
    expect(juce::JSON::toString(s.jobRequirements.encode()).getNumBytesAsUTF8()<1048576,"retained requirements definition bounded");
    s.auditionBusy=[]{return true;};juce::String error;JobRequirements change;change.maximumLraTenths=1;expect(!s.setJobRequirements(change,error),"playing guard blocks requirement changes");s.auditionBusy={};
    groups.add("inline-ui-resources");
}
}
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI init;juce::String rootPath,output;bool passed=false;juce::String failure;
    for(int n=1;n<argc;++n){juce::String arg(argv[n]);if(arg=="--root"&&n+1<argc)rootPath=argv[++n];else if(arg=="--output"&&n+1<argc)output=argv[++n];else return 2;}
    try
    {
        const juce::File root(rootPath);expect(rootPath.isNotEmpty()&&!root.exists()&&root.createDirectory().wasOk(),"new M17 owned test root required");auto presets=root.getChildFile("presets");presets.createDirectory();definitionTests(presets);ruleTests();proposalTests();
        auto stereo=model(root.getChildFile("stereo"),2,48000,65*48000);sessionTests(root.getChildFile("sessions"),stereo);exportTests(root.getChildFile("stereo-export"),stereo);uiAndResources(stereo);manualAndInputTests(root,stereo);
        for(int channels:{1,6}){auto s=model(root.getChildFile("channels-"+juce::String(channels)),channels,48000,65*48000);exportTests(root.getChildFile("export-"+juce::String(channels)),s);}
        for(int offset:{-1,0,1}){auto s=model(root.getChildFile("duration-"+juce::String(offset)),1,44100,60*44100+offset);expect(s.currentState.analysis->detail.lraStable==(offset>=0),"actual under/exactly/above 60s stability boundary");}
        groups.add("actual-layout-duration-boundaries");passed=true;
    }
    catch(const std::exception& e){failure=e.what();std::cerr<<failure<<'\n';}
    auto result=juce::var(new juce::DynamicObject());auto* r=result.getDynamicObject();r->setProperty("suite","M17 job LRA requirements");r->setProperty("passed",passed);r->setProperty("failure",failure);r->setProperty("checks",checks);r->setProperty("groups",groups);r->setProperty("resources",resources);r->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);r->setProperty("video_loaded",false);r->setProperty("physical_devices_opened",false);
    if(output.isNotEmpty())juce::File(output).replaceWithText(juce::JSON::toString(result,false));std::cout<<"M17 "<<(passed?"PASS":"FAIL")<<" "<<checks.size()<<" checks\n";return passed?0:1;
}

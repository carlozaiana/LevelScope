#include "CompositionTestFixtures.h"
#include "Standalone/StandaloneSectionsPanel.h"
#include <algorithm>
namespace
{
// Independent slow oracle: locate explicit plateaus/ramps by scanning, without
// calling the production canonicalizer, radius, evaluator or processor.
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
void envelopeCases()
{
    for(int rate:{44100,48000,96000,192000})
    {
        processing::SectionGainPlan p;p.sampleRate=rate;p.frames=rate;p.inputIdentity=juce::String::repeatedString("a",64);
        p.intervals={{0,2,-12},{2,5,12},{5,rate/2,-3.305074118944},{rate/2,rate,0}};
        expect(p.validationError().isEmpty(),"valid minimum/short/ordinary transitions at all rates");
        const auto canonical=p.canonical();bool all=true;
        for(int64_t n=0;n<rate;++n)all&=canonical.gainAt(n)==oracleGain(p,n);
        expect(all,"every envelope sample equals independent boundary oracle");
        auto invalid=p;invalid.intervals[0].end=1;invalid.intervals[1].begin=1;
        expect(invalid.validationError().isNotEmpty(),"nonzero unrepresentable transition rejected");
        invalid=p;invalid.intervals[1].db=std::numeric_limits<double>::quiet_NaN();expect(invalid.validationError().isNotEmpty(),"nonfinite section gain rejected");
        invalid=p;invalid.intervals[1].db=std::nextafter(12.,INFINITY);expect(invalid.validationError().isNotEmpty(),"immediately outside +12 rejected without tolerance");
        invalid=p;invalid.intervals[1].db=std::nextafter(-12.,-INFINITY);expect(invalid.validationError().isNotEmpty(),"immediately outside -12 rejected");
        invalid=p;invalid.intervals[1].begin++;expect(invalid.validationError().isNotEmpty(),"gap rejected");
        invalid=p;invalid.schema=2;expect(invalid.validationError().isNotEmpty(),"unknown envelope version rejected");
        invalid=p;invalid.intervals.back().end--;expect(invalid.validationError().isNotEmpty(),"incomplete programme coverage rejected");
        auto merged=p;merged.intervals={{0,2,1.25},{2,rate,1.25}};auto one=merged;one.intervals={{0,rate,1.25}};
        expect(merged.audioKey()==one.audioKey(),"equal-gain boundaries excluded from audible identity");
    }
    processing::SectionGainPlan many;many.sampleRate=48000;many.frames=100000;many.inputIdentity=juce::String::repeatedString("b",64);
    for(int n=0;n<10000;++n)many.intervals.push_back({n*10,(n+1)*10,n%2?12.:-12.});
    const auto start=juce::Time::getMillisecondCounterHiRes();expect(many.validationError().isEmpty(),"10000 intervals supported");
    expect(sizeof(many)+many.intervals.capacity()*sizeof(processing::GainInterval)<4*1024*1024,"plan memory below 4 MiB at 10000 intervals");
    double checksum=0;for(int n=0;n<100000;++n)checksum+=many.gainAt(n);
    expect(std::isfinite(checksum)&&juce::Time::getMillisecondCounterHiRes()-start<1000,"bounded 10000-interval lookup responsiveness");
    ProcessingPlan plan;plan.sectionGains=many;plan.engineVersion="section-gain-v1";plan.schemaVersion=2;
    juce::AudioBuffer<float> samples(12,1000);for(int ch=0;ch<12;++ch)juce::FloatVectorOperations::fill(samples.getWritePointer(ch),.5f,1000);
    const auto expected=oracleAudio(samples,plan);
    expect(ProgrammeGainProcessor(plan).process(samples.getArrayOfWritePointers(),12,1000,0)&&sameSamples(samples,expected),"channel-driven evaluator retains 12-channel extension capacity without qualifying layouts");
    expect(!ProgrammeGainProcessor(plan).process(samples.getArrayOfWritePointers(),12,1000,99999),"processing refuses range overrun");
}
void draftThree(StandaloneSessionModel& s)
{
    s.synchronizeSections();const auto frames=s.sourceDocument.media.frames;juce::String error;
    expect(s.sectionDocument.split(1,frames/3,error)&&s.sectionDocument.split(2,frames*2/3,error),"manual exact-frame sections");s.documentChanged();
    expect(s.beginSectionTrimDraft()&&s.setSectionDraftGain(0,-6.02059991328)&&s.setSectionDraftGain(1,3.125)&&s.setSectionDraftGain(2,0),"manual draft contains exact full-precision gains");
    expect(s.applySectionTrimDraft(),"explicit section Apply");
}
StandaloneSessionModel inputModel(const juce::File& root,int rate,int channels,int mode,juce::AudioBuffer<float>& raw)
{
    const int frames=rate*4;StandaloneSessionModel s;
    if(mode>=2)
    {
        auto preset=channels==2?CompositionPreset::stereoStereo:mode==2?CompositionPreset::surroundSurround:CompositionPreset::surroundMonoCentre;
        auto f=fixtureComposition(root,rate,preset,mode==3,frames);raw=f.expected;s=compositionModel(f.recipe);
    }
    else
    {
        root.createDirectory();raw.setSize(channels,frames);
        for(int ch=0;ch<channels;++ch)for(int n=0;n<frames;++n)raw.setSample(ch,n,float(.003*(ch+1)*std::sin((n+ch*7)*.11)));
        auto c=component(root,"mix",rate,raw,mode==1);
        if(mode==1)expect(s.setValidatedProgramme(*c.assembly),"adopt aligned mono group");
        else{expect(s.setValidatedSource(c.file,c.media),"adopt interleaved source");s.sourceDocument.layout=c.layout;}
        expect(s.initializeCurrentStateFromSource(),"initialize Current");s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
    }
    return s;
}
void playbackParity(StandaloneSessionModel& s,const juce::AudioBuffer<float>& expected,bool seekCases)
{
    const int channels=expected.getNumChannels(),frames=expected.getNumSamples();
    auto device=std::make_unique<Device>();auto* d=device.get();StandaloneAudition audio(std::move(device));OutputRouting muted;
    s.selectAuditionDomain(MeasurementDomain::currentState);auto r=s.auditionRequest();
    audio.synchronize(r,muted);expect(audio.connect(int(r.context.media.sampleRate)).isEmpty(),"audio-only injected device");
    Capture played(channels,frames);audio.setCaptureObserver(&played);expect(audio.play(r,muted,true,true).isEmpty(),"section audition starts explicitly");awaitPlaying(audio,*d);drain(audio,*d,997);
    expect(played.valid&&played.next==frames&&sameSamples(played.data,expected)&&audio.underruns()==0,"all-channel offline and pre-monitor samples identical with every route muted and continuous EOF");
    const auto* copy=audio.preparedCopy();audio.setCaptureObserver(nullptr);
    if(seekCases)for(int64_t first:{int64_t(1),int64_t(frames/3-13),int64_t(frames/3+23),int64_t(frames-300)})
    {
        expect(audio.seek(first,r,muted).isEmpty(),"exact audio-only seek around transition/EOF");Capture captured(channels,frames-int(first),first);audio.setCaptureObserver(&captured);
        expect(audio.play(r,muted,false,true).isEmpty(),"restart from exact audio frame");awaitPlaying(audio,*d);drain(audio,*d,127);
        bool same=captured.valid&&captured.next==frames;
        for(int ch=0;ch<channels;++ch)same&=std::memcmp(captured.data.getReadPointer(ch),expected.getReadPointer(ch,int(first)),size_t(frames-first)*4)==0;
        expect(same&&audio.preparedCopy()==copy&&audio.positionFrames()==frames,"seek uses unchanged verified source copy and final source coordinates");audio.setCaptureObserver(nullptr);
    }
    audio.disconnect();
}
void matrixCase(const juce::File& root,int rate,int channels,int mode,bool thorough)
{
    juce::AudioBuffer<float> raw;auto s=inputModel(root,rate,channels,mode,raw);
    auto original=measuredComposition(s,MeasurementDomain::source);measuredComposition(s,MeasurementDomain::currentState);
    const auto savedSource=s.source.analysis;draftThree(s);expect(s.source.analysis==savedSource&&!s.currentState.hasMeasurement&&s.currentStateNeedsRemeasurement(),"Apply preserves Source and discards previous Current measurement");
    expect(s.setCurrentGainDb(-3.305074118944),"full precision global trim after sections");
    Capture offline(channels,raw.getNumSamples());auto current=measuredComposition(s,MeasurementDomain::currentState,&offline);
    const auto expected=oracleAudio(raw,s.currentStateDocument.processingPlan());
    expect(offline.valid&&sameSamples(offline.data,expected),"offline section plus global samples equal independent composition/section/trim oracle");
    playbackParity(s,expected,thorough);
    recordSamples(offline,rate,s.currentStateDocument.staticGainDb,"M12 matrix input-"+juce::String(mode),current.detail);
    if(thorough)
    {
        const auto planKey=s.currentStateDocument.sectionGains->audioKey();const auto draft=proposal::draftGlobalTrim(s.currentTrimRequest());
        expect(draft.status==proposal::TrimStatus::draft&&s.acceptTrimProposal(draft)&&s.applyTrim(),"M3 proposes/applies only final global trim above fixed sections");
        expect(s.currentStateDocument.sectionGains->audioKey()==planKey&&!s.currentState.hasMeasurement,"M3 preserves section plan and invalidates measured Current");
        measuredComposition(s,MeasurementDomain::currentState);expect(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"remeasured actual section-processed Current passes selected rules");
        const auto revision=s.revision;const auto analysis=s.currentState.analysis;
        const auto noChange=proposal::draftGlobalTrim(s.currentTrimRequest());expect(noChange.status==proposal::TrimStatus::noChange&&s.acceptTrimProposal(noChange)&&s.applyTrim()&&s.revision==revision&&s.currentState.analysis==analysis,"already-pass global trim retains exact plan and measurement");
        JobControl control;Capture rendered(channels,raw.getNumSamples());auto request=ex::capture(s);auto artifact=ex::render(request,root,control,{},&rendered);
        if(!artifact.value)std::cerr<<artifact.message<<'\n';expect(artifact.value.has_value(),"section-processed render and actual artifact verification");
        expect(sameSamples(rendered.data,oracleAudio(raw,s.currentStateDocument.processingPlan())),"pre-quantization rendering uses same samples without dither in processing");
        const auto rep=ex::report(request,artifact.value->get());expect(rep["current_processing_plan"]["applied_intervals"].size()==3&&rep["current_processing_plan"]["transitions"].size()==2,"report distinguishes exact applied plan and transition frames");
        auto oracleFile=root.getChildFile("independent-input.f32").createOutputStream();expect(oracleFile!=nullptr,"independent raw oracle capture stream");
        for(int ch=0;ch<channels;++ch)expect(oracleFile->write(raw.getReadPointer(ch),size_t(raw.getNumSamples())*4),"write planar independent input for Python quantizer oracle");
        oracleFile->flush();expect(oracleFile->getStatus().wasOk(),"independent input oracle finalized");oracleFile.reset();
        expect((*artifact.value)->analysis().context.processing.sectionGains==std::nullopt&&(*artifact.value)->analysis().context.processing.staticGainDb==0,"finalized artifact measured without section gain or second global gain");
        expect(ex::publish(**artifact.value,request,root.getChildFile("export"),ex::approvalIdentity(**artifact.value,request),control).value.has_value(),"strict Current and actual artifact gates permit explicit local publication");
        auto document=sf::capture(s);const auto file=root.getChildFile("section.levelscope");expect(sf::save(document,file,control).succeeded(),"schema4 exact applied plan save");auto read=sf::read(file,control);
        expect(read.document&&read.document->loadedSchema==4&&read.document->current.sectionGains->audioKey()==planKey,"schema4 full-precision plan round trip");
        auto verified=sf::verify(*read.document,control);StandaloneSessionModel reopened;juce::String message;
        expect(verified.prepared&&sf::adopt(*verified.document,*verified.prepared,reopened,message)&&reopened.currentState.hasMeasurement&&reopened.source.hasMeasurement,"reopen restores valid caches after per-member verification");
        expect(reopened.sessionId!=s.sessionId&&reopened.sourceDocument.routing.fileToDevice[0]==-1,"reopen new runtime with all routes muted");
        auto late=current;expect(!reopened.acceptMeasurement(late),"old runtime result cannot publish into reopened session");
    }
}
void modelCases(const juce::File& root)
{
    juce::AudioBuffer<float> raw;auto s=inputModel(root,48000,2,0,raw);measuredComposition(s,MeasurementDomain::source);measuredComposition(s,MeasurementDomain::currentState);
    s.synchronizeSections();const auto revision=s.revision;const auto measured=s.currentState.analysis;
    expect(s.beginSectionTrimDraft()&&s.applySectionTrimDraft()&&s.revision==revision&&s.currentState.analysis==measured&&!s.currentStateDocument.sectionGains,"all-zero section draft preserves legacy identity and measurement");
    expect(s.beginSectionTrimDraft()&&s.setSectionDraftGain(0,1.234567890123456),"draft edits full precision");
    bool blocked=false;try{sf::capture(s);}catch(const std::exception&){blocked=true;}expect(blocked,"saving cannot silently drop unapplied draft");
    juce::String error;expect(s.sectionDocument.rename(1,"Renamed",error)&&s.sectionDraftFresh(),"cosmetic name preserves draft identity");
    expect(s.sectionDocument.split(1,48000,error)&&!s.sectionDraftFresh()&&!s.applySectionTrimDraft(),"geometry edit stales draft without mutating audio");
    expect(s.revision==revision&&s.currentState.analysis==measured,"annotation edits retain passing or failing measured results unchanged");s.discardSectionTrimDraft();
    expect(s.beginSectionTrimDraft()&&s.setSectionDraftGain(0,-2)&&s.applySectionTrimDraft(),"new explicit draft can apply");
    const auto active=s.currentStateDocument.sectionGains->audioKey();const auto appliedRevision=s.revision;
    expect(s.sectionDocument.moveEnd(1,60000,error)&&s.currentStateDocument.sectionGains->audioKey()==active&&s.revision==appliedRevision,"moving annotations never moves applied gain boundaries");
    s.measurementUnsuccessful(MeasurementDomain::currentState,ErrorCode::cancelled);
    expect(s.currentStateNeedsRemeasurement()&&!s.currentState.hasMeasurement&&s.currentStateDocument.sectionGains->audioKey()==active,"cancelled verification retains applied gains and unverified state");
    bool busy=true;s.auditionBusy=[&]{return busy;};expect(!s.beginSectionTrimDraft()&&!s.setCurrentGainDb(2),"active audition rejects new processing edits");busy=false;
    expect(s.beginSectionTrimDraft()&&s.setSectionDraftGain(0,0)&&s.setSectionDraftGain(1,0)&&s.applySectionTrimDraft()&&!s.currentStateDocument.sectionGains&&s.revision>appliedRevision,"explicit reset to unity is an audio change and retains no stale pass");
    for(int channels:{1,5,8})
    {
        auto request=fixture(root.getChildFile("legacy-"+juce::String(channels)+".wav"),48000,channels,48000);
        StandaloneSessionModel old;expect(old.setValidatedSource(request.context.sourceFile,request.context.media),"legacy source import");old.sourceDocument.layout=request.context.fileLayout;expect(old.initializeCurrentStateFromSource(),"legacy Current initialized");
        expect(!old.beginSectionTrimDraft(),"mono/5.0/7.1 section processing explicitly unavailable");
        expect(OfflineAnalysisService::validateRequest(old.measurementContext(MeasurementDomain::currentState),MeasurementDomain::currentState).isEmpty(),"legacy mono/5.0/7.1 static gain remains supported");
    }
}
void formatCases(const juce::File& root)
{
    root.createDirectory();
    for(const char* extension:{"wav","aiff","flac"})for(int bits:{16,24})
    {
        auto r=fixture(root.getChildFile(juce::String(extension)+"-"+juce::String(bits)+"."+extension),48000,2,195013,0,bits);
        StandaloneSessionModel s;expect(s.setValidatedSource(r.context.sourceFile,r.context.media),"format imported");s.sourceDocument.layout=r.context.fileLayout;
        expect(s.initializeCurrentStateFromSource(),"format Current initialized");draftThree(s);expect(s.setCurrentGainDb(2.125),"format global trim");
        juce::AudioFormatManager formats;registerSourceFormats(formats);std::unique_ptr<juce::AudioFormatReader> reader(formats.createReaderFor(r.context.sourceFile));
        juce::AudioBuffer<float> raw(2,int(r.context.media.frames));expect(reader&&reader->read(raw.getArrayOfWritePointers(),2,0,raw.getNumSamples()),"independent raw decoder read");
        const auto expected=oracleAudio(raw,s.currentStateDocument.processingPlan());Capture actual(2,raw.getNumSamples());measuredComposition(s,MeasurementDomain::currentState,&actual);
        expect(sameSamples(expected,actual.data),"WAV AIFF FLAC section offline sample parity");playbackParity(s,expected,true);
    }
}
void identityAndRecovery(const juce::File& root)
{
    juce::AudioBuffer<float> raw;auto s=inputModel(root,48000,6,3,raw);measuredComposition(s,MeasurementDomain::source);draftThree(s);measuredComposition(s,MeasurementDomain::currentState);
    const auto good=s.measurementContext(MeasurementDomain::currentState);const auto originalAnalysis=s.currentState.analysis;const auto originalSource=s.source.analysis;
    const auto id=OfflineAnalysisService::identityFor(good,MeasurementDomain::currentState,good.media.sha256);
    for(int n=0;n<8;++n)
    {
        auto bad=good;
        if(n==0)bad.processing.sectionGains->inputIdentity=juce::String::repeatedString("c",64);
        if(n==1)std::swap(bad.fileLayout.roles[0],bad.fileLayout.roles[3]);
        if(n==2)bad.processing.sectionGains->frames--;
        if(n==3)bad.processing.sectionGains->sampleRate=44100;
        if(n==4)bad.processing.engineVersion="unknown-section-v2";
        if(n==5)bad.processing.sectionGains->intervals[0].db=INFINITY;
        if(n==6)bad.processing.schemaVersion=1;
        if(n==7)bad.processing.graphJson="unknown";
        expect(OfflineAnalysisService::validateRequest(bad,MeasurementDomain::currentState).isNotEmpty(),"mismatched roles input frames rate version or nonfinite plan rejected");
        expect(evaluateSnapshot(originalAnalysis.get(),bad,MeasurementDomain::currentState,s.getSelectedTargetProfile(),s.programmeApplicability).verdict!=proposal::Verdict::pass,"invalid plan cannot retain previous measured pass");
    }
    auto changed=good;changed.processing.sectionGains->intervals[0].db=std::nextafter(changed.processing.sectionGains->intervals[0].db,INFINITY);
    expect(OfflineAnalysisService::identityFor(changed,MeasurementDomain::currentState,changed.media.sha256)!=id,"one representable gain step changes measurement identity before display rounding");
    expect(s.beginSectionTrimDraft()&&s.applySectionTrimDraft()&&s.currentState.analysis==originalAnalysis,"identical applied section plan is an exact no-op");
    s.setProgrammeApplicability(proposal::Applicability::unconfirmed);s.setSelectedTargetProfileIndex(0);
    expect(s.currentState.analysis==originalAnalysis&&s.source.analysis==originalSource&&!s.currentStateNeedsRemeasurement(),"profile/applicability changes reuse section-processed measurements");
    s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
    s.setOutputRoute(3,-1);s.referenceVideo.startOnAudio=video::Time(-1001,30000);
    expect(id==OfflineAnalysisService::identityFor(s.measurementContext(MeasurementDomain::currentState),MeasurementDomain::currentState,good.media.sha256),"monitor mute and serialized rational offset preserve audio identity without loading video");
    expect(s.beginSectionTrimDraft()&&s.setSectionDraftGain(0,12)&&s.applySectionTrimDraft(),"applied loud section plan");
    expect(s.setCurrentGainDb(24),"explicit excessive global gain retained");Capture tooLoud(6,raw.getNumSamples());measuredComposition(s,MeasurementDomain::currentState,&tooLoud);
    expect(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::fail&&s.currentStateDocument.staticGainDb==24&&tooLoud.data.getMagnitude(0,raw.getNumSamples())>1,"finite over-full-scale samples measured without clipping; actual failure displayed and gain retained");
    expect(ex::readiness(ex::capture(s)).isNotEmpty(),"failed measured Current cannot render for publication");
    const auto report=ex::report(ex::capture(s),nullptr);expect(!report.hasProperty("artifact")&&report["current"]["evaluation"]["verdict"].toString()=="FAIL","measurement-only failure report invents no artifact verdict/hash");
    expect(s.setCurrentGainDb(-5),"restore explicit test gain");measuredComposition(s,MeasurementDomain::currentState);
    expect(s.setCurrentGainDb(-4.5)&&s.job.begin(JobStage::analysis,s.revision),"start explicit verification of changed section-processed Current");
    JobControl interrupted;
    struct CancelAtSamples final:ProcessedAudioObserver
    {
        JobControl& control;bool called=false;explicit CancelAtSamples(JobControl& c):control(c){}
        void observe(const float* const*,int,int,int64_t)noexcept override{called=true;control.cancel=true;}
    } cancelAtSamples(interrupted);
    const auto cancelledContext=s.measurementContext(MeasurementDomain::currentState);const auto retainedPlan=s.currentStateDocument.sectionGains->audioKey();
    auto cancelled=OfflineAnalysisService(&cancelAtSamples).analyze(cancelledContext,MeasurementDomain::currentState,interrupted);
    expect(cancelAtSamples.called&&cancelled.error==ErrorCode::cancelled&&!cancelled.value,"actual Current worker cancellation after section samples publishes no result");
    expect(!s.job.finish(cancelledContext.ticket,s.revision,cancelled.error,cancelled.message)&&s.job.state==JobState::cancelled,"cancelled Current job finishes without success");s.measurementUnsuccessful(MeasurementDomain::currentState,cancelled.error);
    expect(s.currentStateNeedsRemeasurement()&&!s.currentState.hasMeasurement&&s.currentStateDocument.staticGainDb==-4.5&&s.currentStateDocument.sectionGains->audioKey()==retainedPlan&&s.source.analysis==originalSource,"cancelled section verification preserves applied processing and Source while Current remains unverified");
    measuredComposition(s,MeasurementDomain::currentState);
    JobControl control;const auto doc=sf::capture(s);const auto file=root.getChildFile("recovery.levelscope");expect(sf::save(doc,file,control).succeeded(),"save section recovery fixture");
    auto changedDoc=doc;changedDoc.current.sectionGains->intervals[0].db=13;expect(!sf::save(changedDoc,root.getChildFile("invalid.levelscope"),control).succeeded(),"invalid section core rejected safely");
    changedDoc=doc;changedDoc.current.sectionGainAnnotations.sections[0].end--;changedDoc.current.sectionGainAnnotations.sections[1].begin--;
    expect(!sf::save(changedDoc,root.getChildFile("invalid-provenance.levelscope"),control).succeeded(),"applied boundaries must belong to their original annotation provenance");
    expect(s.beginSectionTrimDraft()&&sf::saved(s,doc,file)&&s.documentDirty(),"draft begun during background save remains unsaved and cannot trigger a pending close");
    s.discardSectionTrimDraft();expect(!s.documentDirty(),"discarding the only pending draft restores saved status without changing audio");
    const auto newer=s.currentStateDocument.staticGainDb+1;expect(s.setCurrentGainDb(newer),"edit after coherent save capture");expect(sf::saved(s,doc,file)&&s.documentDirty(),"older snapshot completion keeps newer edits unsaved");
    auto other=s;other.newDocument();expect(!sf::saved(other,doc,file),"late save cannot affect a new runtime session");
    juce::MemoryBlock data;expect(file.loadFileAsData(data),"read disposable recovery file");auto* bytes=static_cast<uint8_t*>(data.getData());bytes[data.getSize()-1]^=1;
    const auto damaged=root.getChildFile("damaged-cache.levelscope");expect(damaged.replaceWithData(data.getData(),data.getSize()),"create independently damaged cache fixture");auto read=sf::read(damaged,control);
    expect(read.document&&read.document->current.sectionGains->audioKey()==doc.current.sectionGains->audioKey()&&!read.document->currentCache&&read.document->sourceCache,"damaged Current cache discarded with exact applied plan and Source cache recovered");
    expect(!read.document->notices.isEmpty(),"cache recovery explains required remeasurement");
    const auto member=doc.source.composition->files().front().file;const auto held=member.getSiblingFile(member.getFileName()+".held");expect(member.moveFileTo(held),"temporarily move owned test member");
    struct Restore{juce::File from,to;~Restore(){from.moveFileTo(to);}}restore{held,member};
    auto incomplete=sf::verify(*read.document,control);expect(incomplete.prepared&&incomplete.prepared->model.sourceContentInvalid&&!incomplete.prepared->model.currentState.hasMeasurement&&incomplete.prepared->model.currentStateDocument.sectionGains,"missing group member retains processing settings and blocks partial audio/current pass");
}
void editorPrecision(const juce::File& root)
{
    juce::AudioBuffer<float> raw;auto s=inputModel(root,48000,2,0,raw);expect(s.beginSectionTrimDraft()&&s.setSectionDraftGain(0,1e-20),"tiny finite draft gain accepted without display quantization");
    StandaloneSectionsPanel panel(s);panel.synchronize();juce::TextEditor* editor=nullptr;
    for(auto* child:panel.getChildren())if(child->getName()=="Section draft gain dB")editor=dynamic_cast<juce::TextEditor*>(child);
    expect(editor&&std::stod(editor->getText().toStdString())==1e-20,"editable gain text round-trips full double precision including scientific notation");panel.suspend();
}
}
// Audio-only sustained workflow; historical helpers are not runners.
#include "SectionProcessingLongTests.h"
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI gui;juce::String rootPath,output,longSession;bool smoke=false,passed=false;juce::String failure;
    for(int n=1;n<argc;++n){const juce::String arg(argv[n]);if(arg=="--smoke")smoke=true;else if(arg=="--root"&&n+1<argc)rootPath=argv[++n];else if(arg=="--output"&&n+1<argc)output=argv[++n];else if(arg=="--long-session"&&n+1<argc)longSession=argv[++n];else return 2;}
    try
    {
        expect(rootPath.isNotEmpty(),"explicit fresh evidence root required");const juce::File root(rootPath);expect(!root.exists()&&root.createDirectory().wasOk(),"fresh M12 owned fixtures");
        if(longSession.isNotEmpty())longSections(root,juce::File(longSession));
        else { envelopeCases();modelCases(root.getChildFile("model"));
        if(!smoke){formatCases(root.getChildFile("formats"));identityAndRecovery(root.getChildFile("recovery"));editorPrecision(root.getChildFile("editor"));}
        if(smoke)matrixCase(root.getChildFile("smoke"),48000,6,3,true);
        else for(int rate:{44100,48000,96000,192000})for(int channels:{2,6})for(int mode=0;mode<4;++mode)
            matrixCase(root.getChildFile(juce::String(rate)+"-"+juce::String(channels)+"-"+juce::String(mode)),rate,channels,mode,rate==48000);
        }
        passed=true;
    }
    catch(const std::exception& e){failure=e.what();std::cerr<<failure<<'\n';}
    auto v=juce::var(new juce::DynamicObject());auto* r=v.getDynamicObject();r->setProperty("suite","M12 section processing");r->setProperty("passed",passed);r->setProperty("failure",failure);r->setProperty("checks",checks);r->setProperty("cases",cases);r->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);r->setProperty("physical_devices_opened",false);r->setProperty("video_loaded",false);r->setProperty("phase",longSession.isNotEmpty()?"long":smoke?"development-smoke":"audio");r->setProperty("long_workflow",longEvidence);
#if defined(__arm64__) || defined(__aarch64__)
    r->setProperty("architecture","arm64");
#else
    r->setProperty("architecture","x86_64");
#endif
    if(output.isNotEmpty())juce::File(output).replaceWithText(juce::JSON::toString(v,true));std::cout<<"M12 "<<(passed?"PASS":"FAIL")<<" "<<checks.size()<<" checks\n";return passed?0:1;
}

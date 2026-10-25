#include "AuditionTestSupport.h"
#include "Standalone/SessionDocument.h"
#include "Standalone/StandaloneDiagnostics.h"
#include <numeric>
#include <sys/stat.h>
#include <unistd.h>
void beginVideoReaderAudit();
int videoReaderAuditCount(int);
namespace sf=levelscope::standalone::sessionfile;
namespace
{
struct Fixture { ProgrammeAssembly group; juce::File interleaved; int rate=0,frames=0; };
void writePcm(const juce::File& file,int rate,const juce::AudioBuffer<float>& samples,int bits)
{
    juce::WavAudioFormat wav;juce::AiffAudioFormat aiff;juce::FlacAudioFormat flac;
    juce::AudioFormat* format=file.hasFileExtension("flac")?static_cast<juce::AudioFormat*>(&flac):file.hasFileExtension("aiff")?static_cast<juce::AudioFormat*>(&aiff):static_cast<juce::AudioFormat*>(&wav);
    file.deleteFile();auto stream=file.createOutputStream();expect(stream!=nullptr,"M9 fixture stream");
    std::unique_ptr<juce::AudioFormatWriter> writer(format->createWriterFor(stream.get(),rate,samples.getNumChannels(),bits,{},0));
    expect(writer!=nullptr,"M9 fixture encoder");stream.release();expect(writer->writeFromAudioSampleBuffer(samples,0,samples.getNumSamples()),"M9 fixture payload");
}
Fixture makeFixture(const juce::File& root,int channels,int rate,int frames,bool mixed=false,int peak=-1,bool silent=false)
{
    expect(root.createDirectory().wasOk(),"M9 owned fixture folder");Fixture f;f.rate=rate;f.frames=frames;
    juce::AudioBuffer<float> data(channels,frames);
    for(int ch=0;ch<channels;++ch)for(int n=0;n<frames;++n)
    {
        // Exactly PCM16-representable oracle values, independent of decoder/assembler.
        const float x=silent?0:static_cast<float>(std::llround(500*(ch+1)*std::sin((n+ch*7)*(.011+ch*.003))))/32768.f;
        data.setSample(ch,n,x);
    }
    if(peak>=0)data.setSample(peak,frames/2,30000.f/32768);
    f.interleaved=root.getChildFile("interleaved.wav");writePcm(f.interleaved,rate,data,24);
    std::vector<juce::File> files;
    for(int ch=0;ch<channels;++ch)
    {
        const char* ext=mixed?(ch%3==0?".wav":ch%3==1?".aiff":".flac"):".wav";
        auto file=root.getChildFile("mono-"+juce::String(ch)+ext);juce::AudioBuffer<float> mono(1,frames);mono.copyFrom(0,0,data,ch,0,frames);writePcm(file,rate,mono,mixed&&ch%2?16:24);files.push_back(file);
    }
    JobControl control;auto inspected=inspectProgrammeFiles(files,control);if(inspected.error!=ErrorCode::none)std::cerr<<inspected.message<<'\n';
    expect(inspected.error==ErrorCode::none,"all mono members inspected without device");f.group=std::move(inspected.group);
    const auto roles=assemblyRoles(channels);for(int ch=0;ch<channels;++ch)f.group.members[ch].role=roles[ch];
    expect(f.group.orderByRoles()&&f.group.confirmAlignment(),"explicit roles and exact aligned-mix confirmation");return f;
}
StandaloneSessionModel model(const ProgrammeAssembly& group)
{
    StandaloneSessionModel s;expect(s.setValidatedProgramme(group),"atomic programme adoption");expect(s.initializeCurrentStateFromSource(),"group Current initialized");
    s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);return s;
}
void measure(StandaloneSessionModel& s,MeasurementDomain domain,ProcessedAudioObserver* observer=nullptr)
{
    expect(s.job.begin(JobStage::analysis,s.revision),"group measurement begin");auto context=s.measurementContext(domain);JobControl c;
    auto r=OfflineAnalysisService(observer).analyze(context,domain,c);if(!r.value)std::cerr<<r.message<<'\n';expect(r.value.has_value(),"whole-group measurement succeeds");
    expect(s.job.finish(context.ticket,s.revision,r.error,r.message)&&s.acceptMeasurement(std::move(*r.value)),"group measurement current ticket accepted");
}
bool equalMetrics(const measurement::Result& a,const measurement::Result& b)
{
    if(a.sampleRate!=b.sampleRate||a.frames!=b.frames||a.integratedIncluded!=b.integratedIncluded||a.lraIncluded!=b.lraIncluded||a.lraStable!=b.lraStable
        ||a.integratedRelativeGate!=b.integratedRelativeGate||a.lraRelativeGate!=b.lraRelativeGate||a.p10!=b.p10||a.p95!=b.p95
        ||a.channelSamplePeaks!=b.channelSamplePeaks||a.channelTruePeaks!=b.channelTruePeaks||a.samplePeakFrames!=b.samplePeakFrames||a.truePeakEighthFrames!=b.truePeakEighthFrames)return false;
    auto am=std::array{a.integrated,a.lra,a.momentaryMax,a.shortTermMax,a.samplePeak,a.truePeak};auto bm=std::array{b.integrated,b.lra,b.momentaryMax,b.shortTermMax,b.samplePeak,b.truePeak};
    for(size_t n=0;n<am.size();++n)if(am[n].status!=bm[n].status||am[n].value!=bm[n].value)return false;
    for(int kind=0;kind<3;++kind)
    {
        const auto& aw=kind==0?a.momentary:kind==1?a.shortTerm:a.lraWindows;const auto& bw=kind==0?b.momentary:kind==1?b.shortTerm:b.lraWindows;
        if(aw.size()!=bw.size())return false;for(size_t n=0;n<aw.size();++n)if(aw[n].endFrame!=bw[n].endFrame||aw[n].energy!=bw[n].energy)return false;
    }
    return true;
}
ServiceContext interleavedContext(const Fixture& f,const ServiceContext& group)
{
    auto c=group;JobControl job;auto v=validateMedia(f.interleaved,job);expect(v.succeeded(),"independent interleaved oracle validates");c.assembly.reset();c.sourceFile=f.interleaved;c.media=v.media;return c;
}
void recordGroupSamples(const Capture& capture,int rate,double gain,const juce::String& kind,const measurement::Result& result)
{
    recordSamples(capture,rate,gain,kind,result);auto* row=cases.getReference(cases.size()-1).getDynamicObject();
    juce::Array<juce::var> sampleFrames,trueFrames,roles;
    for(auto frame:result.samplePeakFrames)sampleFrames.add(static_cast<juce::int64>(frame));
    for(auto frame:result.truePeakEighthFrames)trueFrames.add(static_cast<juce::int64>(frame));
    for(auto role:assemblyRoles(capture.data.getNumChannels()))roles.add(static_cast<int>(role));
    row->setProperty("sample_peak_frames",sampleFrames);row->setProperty("true_peak_eighth_frames",trueFrames);row->setProperty("ordered_roles",roles);
    row->setProperty("integrated_included",static_cast<juce::int64>(result.integratedIncluded));row->setProperty("lra_included",static_cast<juce::int64>(result.lraIncluded));
}
void permutationsAndIdentity(const juce::File& root)
{
    for(int channels:{2,6})
    {
        auto f=makeFixture(root.getChildFile(juce::String(channels)),channels,48000,129);
        std::vector<int> indices(channels);std::iota(indices.begin(),indices.end(),0);int count=0;
        do
        {
            auto group=f.group;group.members.clear();group.alignmentConfirmation.clear();for(auto i:indices)group.members.push_back(f.group.members[i]);
            expect(group.orderByRoles()&&group.definitionIdentity()==f.group.definitionIdentity(),"picker permutations resolve explicit roles identically");
            expect(group.validate().isNotEmpty()&&group.confirmAlignment(),"ordering alone never confirms alignment");
            std::vector<juce::File> files;for(const auto& m:group.members)files.push_back(m.file);
            auto reader=ProgrammeReader::open(files,group.summary());auto oracle=ProgrammeReader::open({f.interleaved},group.summary());
            juce::AudioBuffer<float> a(channels,129),b(channels,129);
            expect(reader&&oracle&&reader->read(a.getArrayOfWritePointers(),channels,0,129)&&oracle->read(b.getArrayOfWritePointers(),channels,0,129),"complete group block read");
            bool same=true;for(int ch=0;ch<channels;++ch)same&=std::memcmp(a.getReadPointer(ch),b.getReadPointer(ch),129*sizeof(float))==0;
            expect(same,"every permutation produces canonical independently interleaved samples");++count;
        }while(std::next_permutation(indices.begin(),indices.end()));
        expect(count==(channels==2?2:720),"all picker permutations covered");
        for(int which=0;which<5;++which)
        {
            auto g=f.group;if(which==0)g.members[0].media.sha256=juce::String::repeatedString("0",64);
            if(which==1)std::swap(g.members[0].role,g.members[1].role);if(which==2)g.members[0].offsetFrames=1;
            if(which==3)g.members[0].media.frames++;if(which==4)g.members[0].media.sampleRate=44100;
            expect(g.validate().isNotEmpty(),"member content mapping alignment or format change invalidates confirmation");
        }
        auto c=model(f.group).measurementContext(MeasurementDomain::source);c.assembly->alignmentConfirmation.clear();expect(!c.isValid(),"unconfirmed group cannot be measured evaluated or auditioned");
        auto g=f.group;g.members[1].role=g.members[0].role;expect(!g.orderByRoles(),"duplicate role rejected");
        g=f.group;g.members[0].id="renamed-slot";expect(g.definitionIdentity()==f.group.definitionIdentity()&&g.validate().isEmpty(),"non-audible stable slot identifiers do not alter programme identity");
        std::vector<juce::File> files;for(const auto& m:f.group.members)files.push_back(m.file);
        files[1]=files[0];JobControl j;expect(inspectProgrammeFiles(files,j).error!=ErrorCode::none,"duplicate physical member rejected");
        auto link=root.getChildFile("alias.wav");expect(f.group.members[0].file.createSymbolicLink(link,true),"member symlink fixture");files[1]=link;expect(inspectProgrammeFiles(files,j).error!=ErrorCode::none,"symlink member alias rejected");link.deleteFile();
        expect(::link(f.group.members[0].file.getFullPathName().toRawUTF8(),link.getFullPathName().toRawUTF8())==0,"member hardlink fixture");expect(inspectProgrammeFiles(files,j).error!=ErrorCode::none,"hardlink member alias rejected");link.deleteFile();
        expect(f.group.members[0].file.copyFileTo(link),"distinct identical-content member");expect(inspectProgrammeFiles(files,j).error==ErrorCode::none,"distinct files with equal hashes remain legitimate channels");link.deleteFile();
    }
}
void parity(const juce::File& root)
{
    for(int channels:{2,6})for(int rate:{44100,48000,96000,192000})
    {
        auto f=makeFixture(root.getChildFile(juce::String(channels)+"-"+juce::String(rate)),channels,rate,rate/2+17,true);
        auto s=model(f.group);auto device=std::make_unique<Device>();auto* d=device.get();StandaloneAudition a(std::move(device));auto map=routes(channels);
        auto initial=s.auditionRequest();a.synchronize(initial,map);expect(a.prepare(initial).isEmpty(),"group preparation begins device closed");waitPreparation(a);const auto* prepared=a.preparedCopy();expect(prepared&&d->opens==0&&prepared->files().size()==static_cast<size_t>(channels),"no device or interleaved temporary required");
        int64_t bytes=0;for(auto file:prepared->files())bytes+=file.getSize();expect(bytes==s.sourceDocument.media.bytes,"owned payload is exactly aggregate original bytes");
        expect(a.connect(rate).isEmpty(),"explicit injected device only");a.setMonitorGainDb(-50);
        for(double gain:{0.,6.020599913279624,-6.020599913279624,-3.305074118944,-60.,24.})
        {
            s.setCurrentGainDb(gain);s.selectAuditionDomain(MeasurementDomain::currentState);auto r=s.auditionRequest();a.synchronize(r,map);
            Capture offline(channels,f.frames),oracle(channels,f.frames);measure(s,MeasurementDomain::currentState,&offline);JobControl control;
            auto truth=OfflineAnalysisService(&oracle).analyze(interleavedContext(f,r.context),r.domain,control);
            expect(truth.value&&identical(offline,oracle)&&equalMetrics(s.currentState.analysis->detail,truth.value->detail),"all-channel samples and complete metrics exactly equal independent interleaved oracle");
            expect(s.currentState.analysis->identity!=truth.value->identity,"equivalent samples retain distinct programme provenance");
            Capture heard(channels,f.frames);a.setCaptureObserver(&heard);expect(a.play(r,map,true).isEmpty(),"explicit group Current playback");awaitPlaying(a,*d);drain(a,*d,1024);
            expect(identical(heard,offline)&&a.preparedCopy()==prepared,"all-channel pre-monitor capture bit-exact; gain changes reuse group");recordGroupSamples(heard,rate,gain,"M9 mono group Current",s.currentState.analysis->detail);
            const int start=f.frames-37;Capture seek(channels,37,start);a.setCaptureObserver(&seek);expect(a.seek(start,r,map).isEmpty(),"near EOF group seek");
            expect(a.play(r,map,false).isEmpty(),"seek playback explicit");awaitPlaying(a,*d);drain(a,*d,256);
            bool same=seek.valid&&seek.next==f.frames;for(int ch=0;ch<channels;++ch)same&=std::memcmp(seek.data.getReadPointer(ch),offline.data.getReadPointer(ch,start),37*4)==0;
            expect(same&&a.preparedCopy()==prepared,"near-EOF mixed FLAC group seek exact with reused copies");
        }
        s.selectAuditionDomain(MeasurementDomain::source);auto source=s.auditionRequest();a.synchronize(source,map);Capture original(channels,f.frames),played(channels,f.frames);measure(s,MeasurementDomain::source,&original);a.setCaptureObserver(&played);expect(a.play(source,map,true).isEmpty(),"group Source audition");awaitPlaying(a,*d);drain(a,*d,997);expect(identical(original,played),"Source zero processing capture exact");recordGroupSamples(played,rate,0,"M9 mono group Source",s.source.analysis->detail);a.setCaptureObserver(nullptr);a.disconnect();
    }
}
void extendedFormats(const juce::File& root)
{
    expect(root.createDirectory().wasOk(),"extended format folder");const int rate=48000,frames=24001;std::vector<juce::File> files;juce::AudioBuffer<float> expected(6,frames);
    for(int ch=0;ch<6;++ch)
    {
        juce::AudioBuffer<float> mono(1,frames);for(int n=0;n<frames;++n){auto value=static_cast<float>(std::llround(700*(ch+1)*std::sin(n*.013)))/32768.f;mono.setSample(0,n,value);expected.setSample(ch,n,value);}
        auto file=root.getChildFile("channel-"+juce::String(ch)+(ch%3==2?".aiff":".wav"));
        if(ch%3==0)writePcm(file,rate,mono,32); // IEEE float WAV
        else
        {
            auto out=file.createOutputStream();expect(out!=nullptr,"independent PCM32 writer");
            if(ch%3==1)
            {
                out->write("RIFF",4);out->writeInt(frames*4+36);out->write("WAVEfmt ",8);out->writeInt(16);out->writeShort(1);out->writeShort(1);out->writeInt(rate);out->writeInt(rate*4);out->writeShort(4);out->writeShort(32);out->write("data",4);out->writeInt(frames*4);
            }
            else
            {
                out->write("FORM",4);out->writeIntBigEndian(frames*4+46);out->write("AIFFCOMM",8);out->writeIntBigEndian(18);out->writeShortBigEndian(1);out->writeIntBigEndian(frames);out->writeShortBigEndian(32);const int exponent=static_cast<int>(std::floor(std::log2(rate)));out->writeShortBigEndian(static_cast<short>(16383+exponent));out->writeInt64BigEndian(static_cast<juce::int64>(static_cast<uint64_t>(rate)<<(63-exponent)));out->write("SSND",4);out->writeIntBigEndian(frames*4+8);out->writeIntBigEndian(0);out->writeIntBigEndian(0);
            }
            for(int n=0;n<frames;++n){const auto value=static_cast<int32_t>(std::llround(mono.getSample(0,n)*2147483648.));if(ch%3==1)out->writeInt(value);else out->writeIntBigEndian(value);}out->flush();expect(out->getStatus().wasOk(),"independent PCM32 payload complete");
        }
        files.push_back(file);
    }
    JobControl control;auto result=inspectProgrammeFiles(files,control);expect(result.error==ErrorCode::none,"mixed float WAV PCM32 WAV/AIFF mono import");auto roles=assemblyRoles(6);for(int ch=0;ch<6;++ch)result.group.members[ch].role=roles[ch];expect(result.group.confirmAlignment(),"explicit extended-format alignment");auto s=model(result.group);Capture capture(6,frames);measure(s,MeasurementDomain::source,&capture);bool same=capture.valid;for(int ch=0;ch<6;++ch)same&=std::memcmp(capture.data.getReadPointer(ch),expected.getReadPointer(ch),frames*4)==0;expect(same,"mixed 32-bit containers decode exact canonical samples");
    auto unsupported=root.getChildFile("unqualified8.wav");juce::AudioBuffer<float> mono(1,frames);mono.clear();writePcm(unsupported,rate,mono,8);files[0]=unsupported;expect(inspectProgrammeFiles(files,control).error==ErrorCode::unsupportedMedia,"unqualified 8-bit assembly rejected explicitly");
}

void sessions(const juce::File& root)
{
    for(int channels:{2,6})
    {
        auto f=makeFixture(root.getChildFile(juce::String(channels)),channels,48000,192001,true);auto s=model(f.group);s.setCurrentGainDb(-3.305074118944);measure(s,MeasurementDomain::source);measure(s,MeasurementDomain::currentState);
        s.synchronizeSections();juce::String notice;expect(s.sectionDocument.split(1,111111,notice),"group section split");s.documentChanged();s.referenceVideo.startOnAudio=video::Time(-1001,30000);
        const auto confirmation=f.group.alignmentConfirmation;auto saved=sf::capture(s);auto file=root.getChildFile("group.levelscope");JobControl control;
        expect(sf::save(saved,file,control).succeeded(),"group session save");auto read=sf::read(file,control);expect(read.document&&read.document->sourceCache&&read.document->currentCache,"group cache parse retains full descriptor");
        auto verified=sf::verify(*read.document,control);expect(verified.document&&verified.document->audioResolved,"group all-member reopening");StandaloneSessionModel opened;expect(sf::adopt(*verified.document,*verified.prepared,opened,notice),"group atomic adoption");
        expect(opened.source.hasMeasurement&&opened.currentState.hasMeasurement&&opened.sourceDocument.assembly->alignmentConfirmation==confirmation,"unchanged reopening retains confirmation and compatible measurements");
        expect(opened.currentStateDocument.staticGainDb==s.currentStateDocument.staticGainDb&&opened.sectionDocument.map().digest()==s.sectionDocument.map().digest()&&opened.referenceVideo.startOnAudio==s.referenceVideo.startOnAudio,"gain bits sections rational video alignment preserved");
        expect(opened.sessionId!=s.sessionId&&opened.documentId==s.documentId&&!opened.trimDraft,"runtime opening distinct; predictions not restored");
        for(auto route:opened.sourceDocument.routing.fileToDevice)expect(route==-1,"reopen routes muted");
        auto late=*s.currentState.analysis;expect(!opened.acceptMeasurement(late),"old runtime measurement rejected");
        auto group=*read.document->source.assembly;std::vector<juce::File> moved;
        for(size_t i=0;i<group.members.size();++i){auto dest=group.members[i].file.getSiblingFile("moved-"+group.members[i].file.getFileName());expect(group.members[i].file.moveFileTo(dest),"relocate member fixture");moved.push_back(dest);}
        auto missing=sf::verify(*read.document,control);expect(missing.document&&!missing.document->audioResolved,"missing members make whole group incomplete");expect(sf::install(*missing.document,opened,notice)&&opened.sourceContentInvalid&&!opened.currentState.hasMeasurement&&!opened.source.hasMeasurement,"unresolved session cannot inherit measured pass");
        expect(opened.auditionRequest().validationError().isNotEmpty(),"incomplete group cannot audition substitute channels");
        for(size_t i=0;i<moved.size();++i)
        {
            auto d=sf::capture(opened);expect(sf::relinkAudioMember(d,group.members[i].id,moved[i]),"explicit per-slot relink");auto r=sf::verify(d,control);expect(r.document&&r.document->verifiedAudioMembers.contains(group.members[i].id),"selected relink content verified");expect(sf::install(*r.document,opened,notice),"partial relink settings install safely");
            expect(opened.sourceContentInvalid==(i+1<moved.size()),"no partial group becomes available");
        }
        expect(opened.source.hasMeasurement&&opened.currentState.hasMeasurement&&opened.sourceDocument.assembly->alignmentConfirmation==confirmation,"complete identical-content relink retains original confirmation and compatible caches");
        auto d=sf::capture(opened);auto changed=f.group;changed.members[0].file=moved[0];changed.members[0].media.sha256=juce::String::repeatedString("0",64);expect(changed.validate().isNotEmpty(),"changed member never inherits old alignment confirmation");
        auto wrong=makeFixture(root.getChildFile("different"+juce::String(channels)),channels,48000,192001,false,0);
        expect(sf::relinkAudioMember(d,group.members[0].id,wrong.group.members[0].file),"changed-content relink requested");auto rejected=sf::verify(d,control);expect(rejected.document&&!rejected.document->audioResolved&&!rejected.document->verifiedAudioMembers.contains(group.members[0].id),"same-rate same-length changed content fails relink");
        auto consistent=sf::capture(opened);sf::FileHooks hooks;hooks.reached=[&](const char* stage){if(juce::String(stage)=="replace")opened.setCurrentGainDb(-1.23456789012345);};
        expect(sf::save(consistent,file,control,hooks).succeeded()&&sf::saved(opened,consistent,file)&&opened.documentDirty(),"edits during group save remain unsaved");auto coherent=sf::read(file,control);expect(coherent.document&&coherent.document->current.staticGainDb==consistent.current.staticGainDb,"saved group revision coherent");
        for(const auto& member:consistent.source.assembly->members){auto alias=root.getChildFile("member-alias.levelscope");expect(member.file.createSymbolicLink(alias,true),"member destination alias");expect(!sf::save(consistent,alias,control).succeeded(),"save cannot overwrite any group member");alias.deleteFile();}
        juce::MemoryBlock b;expect(file.loadFileAsData(b),"group document bytes");static_cast<char*>(b.getData())[b.getSize()-1]^=1;auto damaged=root.getChildFile("damaged.levelscope");expect(damaged.replaceWithData(b.getData(),b.getSize()),"group damaged cache fixture");auto recovered=sf::read(damaged,control);expect(recovered.document&&recovered.document->source.assembly&&recovered.document->sourceCache&&!recovered.document->currentCache,"independently damaged cache preserves assembly and settings");
    }
}
struct FaultWriter final:SnapshotWriter
{
    std::unique_ptr<juce::FileOutputStream> stream;juce::String fault;JobControl* cancel=nullptr;
    FaultWriter(const juce::File& f,juce::String mode,JobControl* c):stream(f.createOutputStream()),fault(mode),cancel(c){}
    bool write(const void* bytes,size_t count) override{bool ok=stream&&stream->write(bytes,count);if(cancel)cancel->cancel=true;return ok&&fault!="write";}
    bool flush() override{if(stream)stream->flush();return stream&&stream->getStatus().wasOk()&&fault!="flush";}
};
void failures(const juce::File& root)
{
    auto f=makeFixture(root.getChildFile("inputs"),6,48000,48001);auto s=model(f.group);const auto before=s.sessionId;auto context=s.measurementContext(MeasurementDomain::source);
    auto temp=root.getChildFile("snapshots");expect(temp.createDirectory().wasOk(),"owned snapshot failure folder");SnapshotStorage storage;storage.directory=temp;
    storage.availableBytes=[&](auto&){return context.media.bytes-1;};JobControl low;expect(prepareVerifiedSource(context,low,storage).error==ErrorCode::insufficientDiskSpace,"aggregate free-space preflight");storage.availableBytes={};
    for(int member:{0,2,5})for(const char* mode:{"create","write","flush","cancel"})
    {
        JobControl c;int next=0;storage.createWriter=[&](const juce::File& file)->std::unique_ptr<SnapshotWriter>
        {const bool hit=next++==member;if(hit&&juce::String(mode)=="create")return{};return std::make_unique<FaultWriter>(file,hit?mode:"",hit&&juce::String(mode)=="cancel"?&c:nullptr);};
        auto result=prepareVerifiedSource(context,c,storage);expect(!result.value,"first middle last snapshot failure refuses entire group");expect(temp.findChildFiles(juce::File::findFiles,false).isEmpty(),"failed group cleans every owned partial copy");expect(s.sessionId==before&&s.sourceDocument.assembly->definitionIdentity()==f.group.definitionIdentity(),"failed preparation preserves previous live session");
    }
    storage.createWriter={};JobControl cancelled;cancelled.cancel=true;expect(prepareVerifiedSource(context,cancelled,storage).error==ErrorCode::cancelled,"preparation cancellation before copy");
    for(int member:{0,2,5})
    {
        auto path=f.group.members[member].file;auto absent=path.getSiblingFile("away.wav");expect(path.moveFileTo(absent),"missing-member injection");JobControl c;
        expect(!prepareVerifiedSource(context,c,storage).value,"missing member cannot become silence");expect(!s.setValidatedProgramme(f.group)&&s.sessionId==before,"incomplete adoption preserves previous session");expect(absent.moveFileTo(path),"restore failure fixture");
    }
    auto group=f.group;group.members.pop_back();expect(group.validate().isNotEmpty(),"unsupported member count");
    std::vector<juce::File> files;for(const auto& m:f.group.members)files.push_back(m.file);files[2]=f.interleaved;JobControl c;expect(inspectProgrammeFiles(files,c).error==ErrorCode::unsupportedChannels,"non-mono member rejected");
    auto shortFile=root.getChildFile("one-frame-short.wav");fixture(shortFile,48000,1,48000);files[2]=shortFile;expect(inspectProgrammeFiles(files,c).error==ErrorCode::invalidMetadata,"one-sample mismatch rejected without truncation");
    auto rateFile=root.getChildFile("other-rate.wav");fixture(rateFile,44100,1,48001);files[2]=rateFile;expect(inspectProgrammeFiles(files,c).error==ErrorCode::invalidMetadata,"rate mismatch rejected without resampling");
    auto shifted=root.getChildFile("shifted.wav");fixture(shifted,48000,1,48001);files[2]=shifted;auto compatible=inspectProgrammeFiles(files,c);expect(compatible.error==ErrorCode::none&&compatible.group.alignmentConfirmation.isEmpty(),"equal-length different signal accepted only as unconfirmed candidate; no automatic alignment claim");
}
void lifecycleAndVerification(const juce::File& root)
{
    auto f=makeFixture(root,6,48000,192001);auto s=model(f.group);measure(s,MeasurementDomain::source);measure(s,MeasurementDomain::currentState);
    const auto originalGain=s.currentStateDocument.staticGainDb;const auto originalRevision=s.revision;auto original=s.currentState.analysis;
    auto draft=proposal::draftGlobalTrim(s.currentTrimRequest());expect(draft.status==proposal::TrimStatus::draft&&s.acceptTrimProposal(draft),"group deterministic feasible trim draft");
    expect(s.currentStateDocument.staticGainDb==originalGain&&s.revision==originalRevision&&s.currentState.analysis==original,"group draft alone makes no change");
    expect(s.applyTrim()&&!s.currentState.hasMeasurement&&s.currentStateNeedsRemeasurement(),"explicit Apply alone changes group Current and clears old result");
    const auto gain=s.currentStateDocument.staticGainDb;s.measurementUnsuccessful(MeasurementDomain::currentState,ErrorCode::cancelled);expect(s.currentStateDocument.staticGainDb==gain&&!s.currentState.hasMeasurement,"failed verification leaves group unverified without rollback");
    measure(s,MeasurementDomain::currentState);expect(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"group applied trim passes only after actual measured verification");
    const auto current=s.currentState.analysis;const auto revision=s.revision;auto noChange=proposal::draftGlobalTrim(s.currentTrimRequest());expect(noChange.status==proposal::TrimStatus::noChange&&s.acceptTrimProposal(noChange)&&s.applyTrim(),"already-passing group accepts no-change outcome");
    expect(s.revision==revision&&s.currentState.analysis==current&&s.currentStateDocument.staticGainDb==gain,"already-passing group retains gain revision and measurement");
    auto different=makeFixture(root.getChildFile("other"),6,48000,192001,false,0);auto other=model(different.group);other.sessionId=s.sessionId;measure(other,MeasurementDomain::source);
    expect(trimRequest(other.source.analysis.get(),s.currentState.analysis.get(),s.measurementContext(MeasurementDomain::source),s.measurementContext(MeasurementDomain::currentState),s.getSelectedTargetProfile(),s.programmeApplicability).sourceInputIdentity!=s.currentTrimRequest().currentInputIdentity,"individually valid mismatched groups cannot form matching proposal inputs");
    auto mismatch=*s.source.analysis;mismatch.context.assembly->alignmentConfirmation.clear();expect(!measuredProgrammeForSnapshot(&mismatch,s.measurementContext(MeasurementDomain::source),MeasurementDomain::source,s.programmeApplicability).currentIdentityMatches,"lost alignment confirmation rejects saved measured authority");
    auto path=f.group.members[0].file;juce::MemoryBlock originalBytes;expect(path.loadFileAsData(originalBytes),"member mutation backup");const auto time=path.getLastModificationTime();
    struct Mutator final:ProcessedAudioObserver{juce::File path;juce::MemoryBlock data;bool changed=false;void observe(const float*const*,int,int,int64_t)noexcept override{if(changed)return;changed=true;static_cast<char*>(data.getData())[data.getSize()-1]^=1;path.replaceWithData(data.getData(),data.getSize());path.setLastModificationTime(juce::Time::getCurrentTime()+juce::RelativeTime::seconds(2));}} mutator;mutator.path=path;mutator.data=originalBytes;JobControl control;
    auto changed=OfflineAnalysisService(&mutator).analyze(s.measurementContext(MeasurementDomain::source),MeasurementDomain::source,control);expect(!changed.value&&changed.error==ErrorCode::sourceChanged,"member changed during measurement rejects obsolete result");
    expect(path.replaceWithData(originalBytes.getData(),originalBytes.getSize())&&path.setLastModificationTime(time),"restore owned mutation fixture");
    for(int hit:{0,2,5})
    {
        auto member=f.group.members[hit].file;juce::MemoryBlock bytes;expect(member.loadFileAsData(bytes),"short-read fixture backup");auto stamp=member.getLastModificationTime();SnapshotStorage storage;int next=0;
        storage.createWriter=[&](const juce::File& destination)->std::unique_ptr<SnapshotWriter>{if(next++==hit)::truncate(member.getFullPathName().toRawUTF8(),32);return std::make_unique<FaultWriter>(destination,"",nullptr);};JobControl c;auto prepared=prepareVerifiedSource(s.measurementContext(MeasurementDomain::source),c,storage);expect(!prepared.value,"first middle last short input read rejects entire group");expect(member.replaceWithData(bytes.getData(),bytes.getSize())&&member.setLastModificationTime(stamp),"restore short-read fixture");
    }
    std::vector<juce::File> files;for(const auto& m:f.group.members)files.push_back(m.file);auto nonfinite=root.getChildFile("nonfinite.wav");juce::AudioBuffer<float> invalid(1,192001);invalid.clear();invalid.setSample(0,31,std::numeric_limits<float>::quiet_NaN());writePcm(nonfinite,48000,invalid,32);files[2]=nonfinite;JobControl invalidControl;
    expect(inspectProgrammeFiles(files,invalidControl).error==ErrorCode::nonFiniteAudio,"non-finite member fails whole import without silence replacement");
    int notifications=0;auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));StandaloneAudioPanel panel(s,audition,[&]{++notifications;});auto away=path.getSiblingFile("temporarily-missing.wav");expect(path.moveFileTo(away),"idle missing-member fixture");panel.synchronize();expect(s.sourceContentInvalid&&!s.source.hasMeasurement&&!s.currentState.hasMeasurement&&notifications>0&&injected->opens==0,"idle missing member immediately clears visible measured authority without opening a device");expect(away.moveFileTo(path),"restore idle missing-member fixture");expect(s.sourceContentInvalid,"member return requires explicit complete relink verification");
}

void mutedPeakAndWorkflow(const juce::File& root)
{
    for(int high:{3,5})
    {
        auto f=makeFixture(root.getChildFile(juce::String(high)),6,48000,192001,true,high);auto s=model(f.group);measure(s,MeasurementDomain::source);measure(s,MeasurementDomain::currentState);
        auto source=s.source.analysis,current=s.currentState.analysis;auto verdict=s.evaluation(MeasurementDomain::currentState);auto draft=proposal::draftGlobalTrim(s.currentTrimRequest());
        expect(current->detail.channelTruePeaks[high]==*std::max_element(current->detail.channelTruePeaks.begin(),current->detail.channelTruePeaks.end()),"unique highest peak on LFE or surround");
        for(int c=0;c<6;++c)s.setOutputRoute(c,c==high?-1:c);s.referenceVideo.startOnAudio=video::Time(7,3);
        expect(source==s.source.analysis&&current==s.currentState.analysis&&s.evaluation(MeasurementDomain::currentState).snapshotIdentity==verdict.snapshotIdentity,"muting highest peak and video offsets preserve measurements/evaluation");
        auto other=proposal::draftGlobalTrim(s.currentTrimRequest());expect(draft.totalGainDb==other.totalGainDb,"muted highest peak still constrains M3 proposal");
        auto device=std::make_unique<Device>();auto* d=device.get();StandaloneAudition a(std::move(device));auto request=s.auditionRequest();a.synchronize(request,s.sourceDocument.routing);a.connect(48000);a.setMonitorGainDb(-40);Capture heard(6,f.frames),offline(6,f.frames);JobControl c;
        expect(OfflineAnalysisService(&offline).analyze(request.context,request.domain,c).value.has_value(),"muted-channel offline reference");a.setCaptureObserver(&heard);expect(a.play(request,s.sourceDocument.routing,true).isEmpty(),"explicit partial monitoring");awaitPlaying(a,*d);drain(a,*d,1024);expect(identical(heard,offline),"muted LFE/surround remain present pre-monitor");a.setCaptureObserver(nullptr);a.disconnect();
        s.setSelectedTargetProfileIndex(0);expect(s.currentState.analysis==current,"profile switch reuses group measurements");s.setSelectedTargetProfileIndex(1);
        s.setCurrentGainDb(-1.23456789012345);expect(s.source.analysis==source&&!s.currentState.hasMeasurement&&s.currentStateNeedsRemeasurement(),"processing invalidates only Current group measurement");s.measurementUnsuccessful(MeasurementDomain::currentState,ErrorCode::cancelled);expect(s.currentStateDocument.staticGainDb==-1.23456789012345&&!s.currentState.hasMeasurement,"cancelled group verification retains gain unverified");
    }
}
}
#include "ProgrammeQualificationSupport.h"
#include "ProgrammeExtendedOracle.h"
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI init;juce::String output,video,prores,exchangeWrite,exchangeRead,legacy,longInput,generate,prepareInput,destination,realRequest;int rate=48000,seconds=10800;bool only=false,extended=false;
    for(int n=1;n<argc;++n){auto arg=juce::String(argv[n]);if(arg=="--only")only=true;else if(arg=="--extended-oracle")extended=true;else if(n+1<argc){auto value=juce::String(argv[++n]);if(arg=="--output")output=value;else if(arg=="--samples")sampleDirectory=juce::File(value);else if(arg=="--video")video=value;else if(arg=="--prores")prores=value;else if(arg=="--exchange-out")exchangeWrite=value;else if(arg=="--exchange-in")exchangeRead=value;else if(arg=="--m8-session")legacy=value;else if(arg=="--long-programme")longInput=value;else if(arg=="--generate-long")generate=value;else if(arg=="--confirmed-request")realRequest=value;else if(arg=="--prepare-session")prepareInput=value;else if(arg=="--destination")destination=value;else if(arg=="--rate")rate=value.getIntValue();else if(arg=="--seconds")seconds=value.getIntValue();}}
    if(sampleDirectory!=juce::File())sampleDirectory.createDirectory();Temp temp;bool passed=false;juce::String failure;
    try{if(!only){permutationsAndIdentity(temp.directory.getChildFile("permutations"));failures(temp.directory.getChildFile("failures"));extendedFormats(temp.directory.getChildFile("formats"));sessions(temp.directory.getChildFile("sessions"));lifecycleAndVerification(temp.directory.getChildFile("lifecycle"));parity(temp.directory.getChildFile("parity"));mutedPeakAndWorkflow(temp.directory.getChildFile("peaks"));}
        if(video.isNotEmpty()&&prores.isNotEmpty())nativeReferences(temp.directory.getChildFile("video"),juce::File(video),juce::File(prores));
        if(exchangeWrite.isNotEmpty())exchangeOut(juce::File(exchangeWrite));if(exchangeRead.isNotEmpty())exchangeIn(juce::File(exchangeRead));
        if(legacy.isNotEmpty())
        {
            const juce::File input(legacy);
            if(input.isDirectory())
                for(const char* name:{"stereo.levelscope","5.1.levelscope","7.1.levelscope"})
                    migrateM8(input.getChildFile(name),temp.directory.getChildFile(name));
            else migrateM8(input,temp.directory.getChildFile("migrated.levelscope"));
        }
        if(extended)extendedOracles(temp.directory.getChildFile("extended-oracles"));
        if(realRequest.isNotEmpty())confirmedRealProgrammes(juce::File(realRequest),juce::File(output).getParentDirectory());
        if(prepareInput.isNotEmpty())prepareProgramme(juce::File(prepareInput),juce::File(destination));
        if(generate.isNotEmpty())generateLong(juce::File(generate),rate,seconds);
        if(longInput.isNotEmpty())longProgramme(juce::File(longInput),juce::File(output).getSiblingFile("long-programme.levelscope"));passed=true;}
    catch(const std::exception& e){failure=e.what();std::cerr<<failure<<'\n';}
    auto v=juce::var(new juce::DynamicObject());auto* o=v.getDynamicObject();o->setProperty("suite","M9 separate-file programme");o->setProperty("passed",passed);o->setProperty("failure",failure);o->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);o->setProperty("os",juce::SystemStats::getOperatingSystemName());o->setProperty("checks",checks);o->setProperty("cases",cases);o->setProperty("reference_cases",referenceCases);o->setProperty("long_programme",longProgrammeEvidence);o->setProperty("real_programmes",realProgrammeEvidence);o->setProperty("m8_migrations",migrationEvidence);o->setProperty("physical_devices_opened",false);o->setProperty("extended_oracles",extendedOracleEvidence);
#if defined(__arm64__) || defined(__aarch64__)
    o->setProperty("architecture","arm64");
#else
    o->setProperty("architecture","x86_64");
#endif
    const auto json=juce::JSON::toString(v,true);if(output.isNotEmpty())juce::File(output).replaceWithText(json);std::cout<<"M9 "<<(passed?"PASS":"FAIL")<<" checks="<<checks.size()<<" cases="<<cases.size()<<'\n';return passed?0:1;
}

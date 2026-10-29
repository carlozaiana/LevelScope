#include "CompositionTestFixtures.h"
#include "Standalone/StandaloneSectionProposal.h"
#include "Standalone/Pcm24Wave.h"
#include <algorithm>
#include <array>
#include <map>
namespace pcm=levelscope::standalone::pcm24;
namespace
{
juce::var containerEvidence,cinemaEvidence;
void le(std::vector<uint8_t>& b,uint64_t value,int n){for(int i=0;i<n;++i)b.push_back(uint8_t(value>>(i*8)));}
void four(std::vector<uint8_t>& b,const char* id){b.insert(b.end(),id,id+4);}
void setle(std::vector<uint8_t>& b,size_t at,uint64_t value,int n){for(int i=0;i<n;++i)b.at(at+size_t(i))=uint8_t(value>>(i*8));}
uint64_t getle(const uint8_t* b,int n){uint64_t v=0;for(int i=0;i<n;++i)v|=uint64_t(b[i])<<(i*8);return v;}
struct RfFixture{std::vector<uint8_t> bytes;size_t fmt=0,data=0,end=0;int frames=0,channels=0;};
RfFixture rfFixture(int frames=101,int channels=1,bool extended=true,bool zeroCount=true,bool staleOuter=true)
{
    RfFixture x;x.frames=frames;x.channels=channels;auto& b=x.bytes;
    four(b,"RF64");le(b,0xffffffffu,4);four(b,"WAVE");four(b,"ds64");le(b,extended?92:28,4);
    le(b,0,8);le(b,uint64_t(frames)*unsigned(channels)*3,8);le(b,zeroCount?0:uint64_t(frames),8);le(b,0,4);
    if(extended)b.resize(b.size()+64,0);
    four(b,"bext");le(b,4,4);le(b,0,4);
    four(b,"fmt ");le(b,40,4);x.fmt=b.size();le(b,0xfffe,2);le(b,channels,2);le(b,48000,4);le(b,48000*channels*3,4);
    le(b,channels*3,2);le(b,24,2);le(b,22,2);le(b,24,2);le(b,channels==1?4:channels==2?3:0x3f,4);
    const uint8_t guid[]={1,0,0,0,0,0,16,0,128,0,0,170,0,56,155,113};b.insert(b.end(),std::begin(guid),std::end(guid));
    four(b,"data");le(b,0xffffffffu,4);x.data=b.size();
    for(int n=0;n<frames;++n)for(int ch=0;ch<channels;++ch)le(b,uint32_t((n%2?1:-1)*(100000+ch*1000+n)),3);
    if((frames*channels*3)&1)b.push_back(0);x.end=b.size();
    four(b,"regn");le(b,92,4);b.resize(b.size()+92,0);
    four(b,"umid");le(b,24,4);b.resize(b.size()+24,0);
    four(b,"DGDA");le(b,101,4);b.resize(b.size()+101,0x7f);b.push_back(0);
    setle(b,20,staleOuter?x.end+4:b.size()-8,8);return x;
}
juce::File putFixture(const juce::File& root,const juce::String& name,const std::vector<uint8_t>& b)
{const auto f=root.getChildFile(name);expect(!f.exists()&&f.replaceWithData(b.data(),b.size()),"owned RF64 fixture written");return f;}
std::unique_ptr<juce::AudioFormatReader> reader(const juce::File& f)
{juce::AudioFormatManager formats;registerSourceFormats(formats);return std::unique_ptr<juce::AudioFormatReader>(formats.createReaderFor(f));}
void containerCases(const juce::File& root)
{
    root.createDirectory();JobControl job;pcm::Format format;juce::String error,notice;
    for(int channels:{1,2,6})for(int frames:{100,101})for(bool extended:{false,true})
    {
        auto x=rfFixture(frames,channels,extended);const auto f=putFixture(root,juce::String(channels)+"-"+juce::String(frames)+"-"+juce::String(int(extended))+".wav",x.bytes);
        expect(pcm::inspectRf64Input(f,format,error,job,notice)&&format.frames==frames&&format.dataBytes==uint64_t(frames*channels*3),"extended ds64 zero count and metadata import uses exact audio extent");
        expect(notice.contains("derived")&&notice.contains("outer size"),"input compatibility explains derived count and stale outer size");
        expect(!pcm::inspect(f,format,error,job),"input compatibility never relaxes exported artifact validation");
        auto v=validateMedia(f,job);if(!v.succeeded())std::cerr<<v.message<<'\n';expect(v.succeeded()&&v.media.frames==frames,"real import validates complete compatible PCM input");
        auto r=reader(f);juce::AudioBuffer<float> audio(channels,frames+2);audio.clear();
        expect(r&&r->lengthInSamples==frames&&r->read(audio.getArrayOfWritePointers(),channels,0,frames+2),"decoder honors PCM extent and ignores trailing metadata");
        bool equal=true;for(int c=0;c<channels;++c)for(int n=0;n<frames;++n)equal&=audio.getSample(c,n)==float((n%2?1:-1)*(100000+c*1000+n))/8388608.f;
        for(int c=0;c<channels;++c)equal&=audio.getSample(c,frames)==0&&audio.getSample(c,frames+1)==0;
        expect(equal,"every short compatible RF64 decoded sample matches independent PCM bytes");
    }
    const auto base=rfFixture();std::vector<std::pair<juce::String,std::function<void(std::vector<uint8_t>&)>>> bad={
        {"count",[](auto& b){setle(b,36,102,8);}}, {"table",[](auto& b){setle(b,44,1,4);}},
        {"extension",[](auto& b){b[48]=1;}}, {"outer-overflow",[](auto& b){setle(b,20,UINT64_MAX,8);}},
        {"outer-before-audio",[&](auto& b){setle(b,20,base.end-9,8);}},
        {"data-overflow",[](auto& b){setle(b,28,UINT64_MAX,8);}}, {"partial-frame",[](auto& b){setle(b,28,302,8);}},
        {"alignment",[&](auto& b){setle(b,base.fmt+12,2,2);}}, {"byte-rate",[&](auto& b){setle(b,base.fmt+8,1,4);}},
        {"mask",[&](auto& b){setle(b,base.fmt+20,8,4);}}, {"bits",[&](auto& b){setle(b,base.fmt+14,16,2);}},
        {"float",[&](auto& b){setle(b,base.fmt+24,3,4);}}, {"rate",[&](auto& b){setle(b,base.fmt+4,8000,4);}},
        {"sentinel",[&](auto& b){setle(b,base.data-4,303,4);}}, {"truncated",[](auto& b){b.pop_back();}},
        {"truncated-header",[](auto& b){b.push_back(1);}},
        {"duplicate-data",[](auto& b){four(b,"data");le(b,0xffffffffu,4);}},
        {"duplicate-fmt",[](auto& b){four(b,"fmt ");le(b,40,4);b.resize(b.size()+40,0);}},
        {"metadata-overrun",[&](auto& b){setle(b,base.end+4,0xfffffffe,4);}},
        {"unknown-size",[&](auto& b){setle(b,base.end+4,0xffffffffu,4);}},
        {"fact-contradiction",[](auto& b){four(b,"fact");le(b,4,4);le(b,102,4);}},
        {"bw64",[](auto& b){b[0]='B';b[1]='W';}}
    };
    for(const auto& item:bad)
    {auto bytes=base.bytes;item.second(bytes);const auto f=putFixture(root,"invalid-"+item.first+".wav",bytes);expect(!pcm::inspectRf64Input(f,format,error,job,notice)&&error.isNotEmpty(),"malformed or unsupported RF64 rejected with reason");}
    JobControl cancelled;cancelled.cancel=true;const auto cancelledFile=root.getChildFile("1-101-1.wav");
    expect(!pcm::inspectRf64Input(cancelledFile,format,error,cancelled,notice),"RF64 inspection cancellation honored");

    // Real sparse logical payload around the RIFF limit; bounded reads only.
    juce::Array<juce::var> boundaries;
    const int64_t edge=(0xfffffffeLL-60)/3;
    for(int64_t frames:{edge-2,edge-1,edge,edge+1,edge+2,edge+4097})
    {
        expect(pcm::describe(192000,1,frames,format,error,4),"mono 64-bit boundary format within product duration");
        const uint64_t data=uint64_t(frames)*3,pad=data&1;
        expect(format.rf64==(data+pad+60>0xfffffffeULL)&&format.fileBytes==format.dataOffset+data+pad,"RIFF threshold counts header and odd padding independently");
        const auto f=root.getChildFile("sparse-"+juce::String(frames)+".wav");auto out=f.createOutputStream();const auto header=pcm::header(format);
        expect(out&&out->write(header.data(),header.size())&&out->setPosition(int64_t(format.fileBytes-1))&&out->writeByte(0),"construct sparse complete logical payload without full-duration claim");
        const auto tail=format.dataOffset+data-3;expect(out->setPosition(int64_t(tail))&&out->writeByte(0)&&out->writeByte(0)&&out->writeByte(16),"known sample beyond four GiB");out->flush();expect(out->getStatus().wasOk(),"sparse fixture finalized");out.reset();
        pcm::Format inspected;expect(pcm::inspect(f,inspected,error,job)&&inspected.frames==frames,"strict inspection accepts valid complete sparse container");
        auto r=reader(f);juce::AudioBuffer<float> sample(1,1);expect(r&&r->lengthInSamples==frames&&r->read(sample.getArrayOfWritePointers(),1,frames-1,1)&&sample.getSample(0,0)==.125f,"bounded seek beyond four GiB decodes expected last PCM sample");
        auto row=juce::var(new juce::DynamicObject());row.getDynamicObject()->setProperty("frames",juce::String(frames));row.getDynamicObject()->setProperty("logical_bytes",juce::String(format.fileBytes));row.getDynamicObject()->setProperty("rf64",format.rf64);boundaries.add(row);
        r.reset();expect(f.deleteFile(),"remove only owned sparse fixture after bounded checks");
    }
    containerEvidence=juce::var(new juce::DynamicObject());containerEvidence.getDynamicObject()->setProperty("boundaries",boundaries);
    containerEvidence.getDynamicObject()->setProperty("full_duration_qualified",false);
}
double independentGain(const processing::SectionGainPlan& p,int64_t frame)
{
    std::vector<processing::GainInterval> rows;for(auto i:p.intervals)if(!rows.empty()&&rows.back().db==i.db)rows.back().end=i.end;else rows.push_back(i);
    for(size_t n=1;n<rows.size();++n)
    {const auto& a=rows[n-1];const auto& b=rows[n];const auto half=std::min({int64_t(p.sampleRate/20),(a.end-a.begin)/2,(b.end-b.begin)/2});
     const auto begin=b.begin-half,end=b.begin+half;if(frame==begin)return a.db;if(frame==end-1)return b.db;if(frame>begin&&frame<end-1)return a.db+(b.db-a.db)*(double(frame-begin)/double(end-begin-1));}
    for(auto i:rows)if(frame>=i.begin&&frame<i.end)return i.db;throw std::runtime_error("independent sample outside plan");
}
juce::AudioBuffer<float> monoOracle(const juce::AudioBuffer<float>& raw,const ProcessingPlan& plan,int64_t first=0)
{
    auto result=raw;for(int n=0;n<raw.getNumSamples();++n)
    {const auto db=plan.sectionGains?independentGain(*plan.sectionGains,first+n):0;const float section=db==0?raw.getSample(0,n):float(double(raw.getSample(0,n))*std::pow(10.,db/20));result.setSample(0,n,float(double(section)*std::pow(10.,plan.staticGainDb/20)));}return result;
}
StandaloneSessionModel monoModel(const juce::File& root,int rate,ChannelRole role,juce::AudioBuffer<float>& raw)
{
    root.createDirectory();const int frames=rate*10+1;raw.setSize(1,frames);
    for(int n=0;n<frames;++n)raw.setSample(0,n,float(std::llround(32768*.08*std::sin(2*3.141592653589793*1000*n/rate)))/32768.f);
    raw.setSample(0,rate*5+rate/4,.8f);const auto file=root.getChildFile("mono.wav");floatFile(file,rate,raw);
    JobControl job;const auto v=validateMedia(file,job);StandaloneSessionModel s;expect(v.succeeded()&&s.setValidatedSource(file,v.media),"single mono file imported");
    s.setChannelRole(0,role);expect(s.confirmLayout()&&s.initializeCurrentStateFromSource(),"explicit mono role and Current initialized");
    s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
    s.synchronizeSections();juce::String error;for(int n=1;n<20;++n)expect(s.sectionDocument.split(uint64_t(n),n*rate/2,error),"manual mono annotations without audio changes");return s;
}
void monoPlayback(StandaloneSessionModel& s,const juce::AudioBuffer<float>& expected)
{
    auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));OutputRouting muted;
    s.selectAuditionDomain(MeasurementDomain::currentState);const auto request=s.auditionRequest();audition.synchronize(request,muted);
    expect(audition.connect(int(request.context.media.sampleRate)).isEmpty(),"mono injected device only");Capture capture(1,expected.getNumSamples());audition.setCaptureObserver(&capture);
    expect(audition.play(request,muted,true,true).isEmpty(),"explicit mono Current audition");awaitPlaying(audition,*injected);drain(audition,*injected,997);
    expect(capture.valid&&sameSamples(capture.data,expected)&&audition.underruns()==0,"muted mono remains present in exact pre-monitor capture through EOF");
    const auto* copy=audition.preparedCopy();audition.setCaptureObserver(nullptr);
    for(int64_t first:{int64_t(1),int64_t(request.context.media.sampleRate*5-10),int64_t(expected.getNumSamples()-1)})
    {
        expect(audition.seek(first,request,muted).isEmpty(),"mono absolute-frame seek");Capture part(1,expected.getNumSamples()-int(first),first);audition.setCaptureObserver(&part);
        expect(audition.play(request,muted,false,true).isEmpty(),"mono explicit seek restart");awaitPlaying(audition,*injected);drain(audition,*injected,127);
        expect(part.valid&&std::memcmp(part.data.getReadPointer(0),expected.getReadPointer(0,int(first)),size_t(part.data.getNumSamples())*4)==0&&audition.preparedCopy()==copy,"mono seek uses unchanged copy and exact section coefficients");audition.setCaptureObserver(nullptr);
    }
    audition.disconnect();
}
void monoWorkflow(const juce::File& root,int rate,ChannelRole role)
{
    juce::AudioBuffer<float> raw;auto s=monoModel(root,rate,role,raw);measuredComposition(s,MeasurementDomain::source);measuredComposition(s,MeasurementDomain::currentState);
    const auto source=s.source.analysis;JobControl job;
    if(rate==48000)expect(sf::save(sf::capture(s),root.getChildFile("baseline.levelscope"),job).succeeded(),"retain mono baseline for silent UI");
    auto request=s.sectionProposalRequest();auto result=generateSectionProposal(request,job);
    if(result.value)std::cout<<(*result.value)->text()<<'\n';else std::cerr<<result.message<<'\n';
    expect(result.value&&(*result.value)->status()==SectionProposalStatus::ready&&(*result.value)->passes()<=3,"existing bounded method finds measured mono candidate within pass cap");
    expect(s.acceptSectionProposal(*result.value)&&s.applySectionProposal(),"explicit mono section proposal Apply");
    expect(s.source.analysis==source&&!s.currentState.hasMeasurement&&s.currentStateNeedsRemeasurement(),"mono Apply preserves Source and makes Current unverified");
    Capture measured(1,raw.getNumSamples());measuredComposition(s,MeasurementDomain::currentState,&measured);const auto expected=monoOracle(raw,s.currentStateDocument.processingPlan());
    expect(sameSamples(measured.data,expected)&&s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"independent Current mono measurement and sample oracle pass");
    const auto revision=s.revision;const auto current=s.currentState.analysis;auto noChange=generateSectionProposal(s.sectionProposalRequest(),job);
    expect(noChange.value&&(*noChange.value)->passes()==0&&s.acceptSectionProposal(*noChange.value)&&s.revision==revision&&s.currentState.analysis==current,"already passing mono makes no change or private scan");
    s.discardSectionProposal();
    if(rate==48000)monoPlayback(s,expected);
    Capture rendered(1,raw.getNumSamples());const auto exportRequest=ex::capture(s);auto artifact=ex::render(exportRequest,root,job,{},&rendered);
    if(!artifact.value)std::cerr<<artifact.message<<'\n';expect(artifact.value&&sameSamples(rendered.data,expected),"mono verified export shares exact processing samples");
    const auto a=*artifact.value;const uint32_t mask=role==juce::AudioChannelSet::left?1u:role==juce::AudioChannelSet::right?2u:4u;
    expect(a->spec().format.channels==1&&a->spec().format.frames==raw.getNumSamples()&&a->spec().format.mask==mask&&a->spec().outputToInput==std::vector<int>{0},"export remains mono with exact confirmed role and explicit mapping");
    expect(a->file().getSize()==68+int64_t(raw.getNumSamples())*3+1&&a->analysis().context.processing.staticGainDb==0&&!a->analysis().context.processing.sectionGains,"odd PCM24 padding is not audio and artifact has zero extra processing");
    expect(ex::publish(*a,exportRequest,root.getChildFile("export"),ex::approvalIdentity(*a,exportRequest),job).value.has_value(),"explicit mono publication requires measured Current and artifact passes");
    expect(root.getChildFile("independent-input.f32").replaceWithData(raw.getReadPointer(0),size_t(raw.getNumSamples())*4),"retain mono input for independent final PCM oracle");
    auto doc=sf::capture(s);const auto file=root.getChildFile("applied.levelscope");expect(sf::save(doc,file,job).succeeded(),"save applied mono schema4 plan");auto read=sf::read(file,job);
    expect(read.document&&read.document->current.sectionGains&&read.document->source.layout.roles[0]==role,"mono session preserves exact plan and role");
    auto reopened=sf::verify(*read.document,job);expect(reopened.prepared&&reopened.prepared->model.currentState.hasMeasurement&&reopened.prepared->model.sourceDocument.routing.fileToDevice[0]==-1,"mono same-build caches revalidated and reopening muted");
    const auto bound=s.currentStateDocument.sectionGains->audioKey();expect(s.setCurrentGainDb(s.currentStateDocument.staticGainDb+.125),"explicit full-precision mono gain change");
    s.measurementUnsuccessful(MeasurementDomain::currentState,ErrorCode::cancelled);expect(!s.currentState.hasMeasurement&&s.currentStateNeedsRemeasurement()&&s.currentStateDocument.sectionGains->audioKey()==bound,"cancelled mono verification keeps plan unverified without rollback");
    expect(!ex::render(ex::capture(s),root,job).value,"unverified mono cannot render for publication");
    auto row=juce::var(new juce::DynamicObject());auto* r=row.getDynamicObject();r->setProperty("rate",rate);r->setProperty("role",int(role));r->setProperty("frames",raw.getNumSamples());r->setProperty("mask",int(mask));r->setProperty("private_passes",(*result.value)->passes());r->setProperty("artifact_sha256",a->hash());cases.add(row);
}
void manualAndFormatCases(const juce::File& root)
{
    root.createDirectory();JobControl job;
    for(auto extension:{"wav","aiff","flac"})
    {
        auto request=fixture(root.getChildFile(juce::String("format.")+extension),48000,1,192001);
        StandaloneSessionModel s;expect(s.setValidatedSource(request.context.sourceFile,request.context.media),"mono format imported without assembly or composition");s.sourceDocument.layout=request.context.fileLayout;
        expect(s.initializeCurrentStateFromSource(),"mono format initializes Current");s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
        Capture source(1,192001);measuredComposition(s,MeasurementDomain::source,&source);measuredComposition(s,MeasurementDomain::currentState);
        s.synchronizeSections();juce::String error;expect(s.sectionDocument.split(1,96000,error)&&s.beginSectionTrimDraft()&&s.setSectionDraftGain(0,1.234567890123)&&s.setSectionDraftGain(1,-6.02059991328),"mono manual positive negative full-precision draft");
        const auto before=s.currentState.analysis;const auto revision=s.revision;expect(s.currentState.analysis==before&&s.revision==revision,"unapplied mono draft preserves measurements");
        expect(s.applySectionTrimDraft()&&!s.currentState.hasMeasurement,"manual mono Apply invalidates Current");Capture applied(1,192001);measuredComposition(s,MeasurementDomain::currentState,&applied);
        expect(sameSamples(applied.data,monoOracle(source.data,s.currentStateDocument.processingPlan())),"all mono WAV AIFF FLAC manual samples equal independent oracle");
        const auto identity=s.currentState.analysis->identity,plan=s.currentStateDocument.sectionGains->audioKey();const auto audioRevision=s.revision;
        expect(s.sectionDocument.moveEnd(1,100000,error)&&s.currentState.analysis->identity==identity&&s.currentStateDocument.sectionGains->audioKey()==plan&&s.revision==audioRevision,"mono annotation boundary changes leave applied envelope and truth unchanged");
        const auto report=ex::report(ex::capture(s),nullptr);expect(!report.hasProperty("artifact"),"mono measurement-only report has no invented artifact");
        auto bad=s.measurementContext(MeasurementDomain::currentState);bad.fileLayout.roles[0]=juce::AudioChannelSet::LFE;expect(OfflineAnalysisService::validateRequest(bad,MeasurementDomain::currentState).isNotEmpty(),"isolated LFE rejected as mono section programme");
        bad=s.measurementContext(MeasurementDomain::currentState);bad.fileLayout.roles[0]=juce::AudioChannelSet::left;expect(OfflineAnalysisService::validateRequest(bad,MeasurementDomain::currentState).isNotEmpty(),"equal-weight changed role cannot reuse bound section plan");
        const auto old=s.currentState.analysis;s.setChannelRole(0,juce::AudioChannelSet::left);expect(!s.source.hasMeasurement&&!s.currentState.hasMeasurement,"mono role change invalidates Source and Current");
        expect(!s.acceptMeasurement(*old),"stale mono measurement rejected after role change");
    }
    // Independent mono-vs-dual-mono level relationship exposes accidental duplication.
    const int rate=48000,frames=rate*4;juce::AudioBuffer<float> mono(1,frames),dual(2,frames);
    for(int n=0;n<frames;++n){const auto x=float(.1*std::sin(2*3.141592653589793*997*n/rate));mono.setSample(0,n,x);dual.setSample(0,n,x);dual.setSample(1,n,x);}
    measurement::OfflineMeasurement one(rate,{measurement::Role::centre}),two(rate,{measurement::Role::left,measurement::Role::right});
    one.push(mono.getArrayOfReadPointers(),frames);two.push(dual.getArrayOfReadPointers(),frames);const auto a=one.finish(),b=two.finish();
    expect(a.integrated.value&&b.integrated.value&&std::abs(*b.integrated.value-*a.integrated.value-10*std::log10(2.))<1e-10,"one-channel loudness differs from explicit duplicated stereo by calculated energy sum");
    expect(std::abs(*a.integrated.value+23.01)<.05,"independent 997 Hz mono level has no hidden compensation");
}
void rawMonoPcm(const juce::File& root)
{
    root.createDirectory();juce::Array<juce::var> oracles;bool changed=false;
    for(int rate:{44100,48000,96000,192000})for(uint32_t mask:{1u,2u,4u})for(int frames:{1,2,101,102})
    {
        const auto key=juce::String(rate)+"-"+juce::String(mask)+"-"+juce::String(frames);const auto file=root.getChildFile(key+".wav");pcm::Format f;juce::String error;JobControl job;
        expect(pcm::describe(rate,1,frames,f,error,mask),"mono PCM24 role/rate/odd-even format");juce::AudioBuffer<float> raw(1,frames);
        for(int n=0;n<frames;++n)raw.setSample(0,n,float((n%13-6)*.0123456789));
        {pcm::Writer writer(file,f,{},job);writer.write(raw.getArrayOfReadPointers(),{0},frames,0x1a2b3c4d5e6f7788ULL,job);writer.finish(job);}
        pcm::Format parsed;expect(pcm::inspect(file,parsed,error,job)&&parsed.frames==frames&&parsed.mask==mask&&parsed.fileBytes==68+uint64_t(frames)*3+uint64_t(frames&1),"strict mono artifact sizes and padding retain exact frames");
        auto r=reader(file);juce::AudioBuffer<float> decoded(1,frames+1);decoded.clear();expect(r&&r->lengthInSamples==frames&&r->read(decoded.getArrayOfWritePointers(),1,0,frames+1)&&decoded.getSample(0,frames)==0,"decoder never interprets mono padding as a frame");
        for(int n=0;n<frames;++n)changed|=decoded.getSample(0,n)!=raw.getSample(0,n);
        expect(root.getChildFile(key+".f32").replaceWithData(raw.getReadPointer(0),size_t(frames)*4),"retain short raw mono quantizer oracle");
        auto v=juce::var(new juce::DynamicObject());auto* o=v.getDynamicObject();o->setProperty("file",file.getFileName());o->setProperty("capture",key+".f32");o->setProperty("rate",rate);o->setProperty("mask",int(mask));o->setProperty("frames",frames);o->setProperty("seed_hex","1a2b3c4d5e6f7788");oracles.add(v);
    }
    expect(changed,"zero processing gain dithered mono is not a bit-identical copy");
    pcm::Format odd;juce::String reason;JobControl job;expect(pcm::describe(48000,1,1,odd,reason,4),"odd mono failure fixture");
    const auto broken=root.getChildFile("failed-padding.wav");bool failed=false;
    try { pcm::Hooks hook;hook.fail=[](const char* point){return std::strcmp(point,"audio-padding")==0;};pcm::Writer writer(broken,odd,hook,job);float sample=.1f;const float* channels[]={&sample};writer.write(channels,{0},1,1,job);writer.finish(job); }
    catch(const std::exception&){failed=true;}
    pcm::Format inspected;expect(failed&&!pcm::inspect(broken,inspected,reason,job),"failed mono pad write cannot become a verified artifact");
    expect(root.getChildFile("raw-oracles.json").replaceWithText(juce::JSON::toString(oracles,true)),"write independent mono quantizer cases");
}
void cinemaCompare(const juce::File& original,const juce::File& copy)
{
    JobControl job;const auto beforeSize=original.getSize();const auto beforeTime=original.getLastModificationTime();
    auto imported=validateMedia(original,job);if(!imported.succeeded())std::cerr<<imported.message<<'\n';
    expect(imported.succeeded()&&imported.media.sha256=="1fb7f93a33d600bf0a62d63e642e5703f03aec8e6920f8a6e510cf65f61fdeaf","direct original cinema import validates content identity without rewrap");
    expect(imported.media.channels==6&&imported.media.sampleRate==48000&&imported.media.frames==269875200&&imported.media.bitsPerSample==24,"original cinema duration and PCM format unchanged");
    expect(pcm::hashFile(copy,job)=="2b00ba5ae88ef06b5d949adf7745144abf8b38d4cb97f78e8abe94d5916b9af6","preserved M13 comparison copy identity verified");
    auto a=reader(original),b=reader(copy);expect(a&&b&&a->lengthInSamples==b->lengthInSamples&&a->getChannelLayout()==b->getChannelLayout(),"original and M13 copy decoder duration and ordered roles agree");
    juce::Array<juce::var> spans;const int64_t crossing=(0x100000000LL-16384)/18;
    for(int64_t first:{int64_t(0),int64_t(1),int64_t(123456789),crossing-2048,crossing+2048,int64_t(269875200-4096)})
    {
        const int count=4096;juce::AudioBuffer<float> x(6,count),y(6,count);
        expect(a->read(x.getArrayOfWritePointers(),6,first,count)&&b->read(y.getArrayOfWritePointers(),6,first,count)&&sameSamples(x,y),"bounded original cinema decoded samples match PCM-identical copy including beyond four GiB");
        auto row=juce::var(new juce::DynamicObject());row.getDynamicObject()->setProperty("first_frame",juce::String(first));row.getDynamicObject()->setProperty("frames",count);row.getDynamicObject()->setProperty("channels",6);spans.add(row);
    }
    expect(beforeSize==original.getSize()&&beforeTime==original.getLastModificationTime(),"original cinema remains untouched");
    cinemaEvidence=juce::var(new juce::DynamicObject());auto* e=cinemaEvidence.getDynamicObject();e->setProperty("original_sha256",imported.media.sha256);e->setProperty("frames",juce::String(imported.media.frames));e->setProperty("rate",48000);e->setProperty("channels",6);e->setProperty("input_notice",imported.media.layoutProvenance);e->setProperty("bounded_reads",spans);e->setProperty("full_proposal_export_repeated",false);
}
}

#include "MonoLongTests.h"
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI gui;juce::String path,output,cinema,copy;bool containers=false,smoke=false,longRun=false,passed=false;juce::String failure;
    for(int n=1;n<argc;++n){juce::String a(argv[n]);if(a=="--root"&&n+1<argc)path=argv[++n];else if(a=="--output"&&n+1<argc)output=argv[++n];else if(a=="--containers")containers=true;else if(a=="--smoke")smoke=true;else if(a=="--long")longRun=true;else if(a=="--cinema"&&n+1<argc)cinema=argv[++n];else if(a=="--copy"&&n+1<argc)copy=argv[++n];else return 2;}
    try
    {
        const juce::File root(path);expect(path.isNotEmpty()&&!root.exists()&&root.createDirectory().wasOk(),"fresh M14 evidence root required");
        if(longRun)monoLong(root);else if(cinema.isNotEmpty())cinemaCompare(juce::File(cinema),juce::File(copy));else{
        containerCases(root.getChildFile("containers"));rawMonoPcm(root.getChildFile("raw-pcm"));
        if(!containers&&!smoke)manualAndFormatCases(root.getChildFile("manual-formats"));
        if(!containers)for(int rate:{44100,48000,96000,192000})for(auto role:{juce::AudioChannelSet::left,juce::AudioChannelSet::right,juce::AudioChannelSet::centre})
        {if(smoke&&(rate!=48000||role!=juce::AudioChannelSet::centre))continue;monoWorkflow(root.getChildFile(juce::String(rate)+"-"+juce::String(int(role))),rate,role);}
        }
        passed=true;
    }
    catch(const std::exception& e){failure=e.what();std::cerr<<failure<<'\n';}
    auto report=juce::var(new juce::DynamicObject());auto* r=report.getDynamicObject();r->setProperty("suite","M14 mono and bounded RF64 input");r->setProperty("passed",passed);r->setProperty("failure",failure);r->setProperty("checks",checks);r->setProperty("cases",cases);r->setProperty("containers",containerEvidence);r->setProperty("cinema",cinemaEvidence);r->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);r->setProperty("physical_devices_opened",false);r->setProperty("video_loaded",false);
    r->setProperty("long_workflow",longEvidence);r->setProperty("phase",longRun?"long":cinema.isNotEmpty()?"cinema":containers?"containers":smoke?"development-smoke":"short");
#if defined(__arm64__) || defined(__aarch64__)
    r->setProperty("architecture","arm64");
#else
    r->setProperty("architecture","x86_64");
#endif
    if(output.isNotEmpty())juce::File(output).replaceWithText(juce::JSON::toString(report,true));std::cout<<"M14 "<<(passed?"PASS":"FAIL")<<" "<<checks.size()<<" checks\n";return passed?0:1;
}

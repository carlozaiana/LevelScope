#include "AuditionTestSupport.h"
#include "Standalone/StandaloneExport.h"
#include "Standalone/SessionDocument.h"
#include "Standalone/VerifiedSource.h"
#include "Standalone/StandaloneDiagnostics.h"
#include <fstream>
namespace ex=levelscope::standalone::exporting;
namespace sf=levelscope::standalone::sessionfile;
namespace
{
void floatFile(const juce::File& file,int rate,const juce::AudioBuffer<float>& data)
{
    expect(!file.exists(),"oracle fixture never overwrites existing media");auto out=file.createOutputStream();expect(out!=nullptr,"float oracle stream");
    const int bytes=data.getNumChannels()*data.getNumSamples()*4;
    out->write("RIFF",4);out->writeInt(36+bytes);out->write("WAVEfmt ",8);out->writeInt(16);out->writeShort(3);out->writeShort(static_cast<short>(data.getNumChannels()));
    out->writeInt(rate);out->writeInt(rate*data.getNumChannels()*4);out->writeShort(static_cast<short>(data.getNumChannels()*4));out->writeShort(32);out->write("data",4);out->writeInt(bytes);
    for(int n=0;n<data.getNumSamples();++n)for(int c=0;c<data.getNumChannels();++c){uint32_t bits;const float x=data.getSample(c,n);std::memcpy(&bits,&x,4);out->writeInt(static_cast<int>(bits));}
    out->flush();expect(out->getStatus().wasOk(),"float oracle finalized");
}
ProgrammeComponent component(const juce::File& root,const juce::String& name,int rate,const juce::AudioBuffer<float>& data,bool grouped)
{
    ProgrammeComponent p;p.id=name;p.purpose=name=="me"?"music-effects":"dialogue";p.layout=layoutFor(data.getNumChannels());
    JobControl job;
    if(grouped&&data.getNumChannels()>1)
    {
        ProgrammeAssembly group;const auto roles=assemblyRoles(data.getNumChannels());
        for(int c=0;c<data.getNumChannels();++c)
        {
            juce::AudioBuffer<float> mono(1,data.getNumSamples());mono.copyFrom(0,0,data,c,0,data.getNumSamples());
            auto file=root.getChildFile(name+"-"+juce::String(c)+".wav");floatFile(file,rate,mono);auto v=validateMedia(file,job);expect(v.succeeded(),"group component member verified");
            group.members.push_back({"audio-"+juce::String(c),file,v.media,roles[static_cast<size_t>(c)],0});
        }
        expect(group.confirmAlignment(),"test generator explicitly confirms component channel alignment");p.assembly=group;p.media=group.summary();p.layout=group.layout();
    }
    else
    {p.file=root.getChildFile(name+".wav");floatFile(p.file,rate,data);auto v=validateMedia(p.file,job);expect(v.succeeded(),"interleaved component verified");p.media=v.media;}
    if(data.getNumChannels()==1){p.layout.count=1;p.layout.roles[0]=juce::AudioChannelSet::centre;p.layout.confirmed=true;}
    return p;
}
struct CompositionFixture
{
    ProgrammeComposition recipe;
    juce::File oracle;
    juce::AudioBuffer<float> expected;
};
CompositionFixture fixtureComposition(const juce::File& root,int rate,CompositionPreset preset,bool grouped,int frames)
{
    expect(root.createDirectory().wasOk(),"composition fixture directory");const int channels=preset==CompositionPreset::stereoStereo?2:6;
    const int dialogueChannels=preset==CompositionPreset::surroundMonoCentre?1:channels;
    juce::AudioBuffer<float> me(channels,frames),dialogue(dialogueChannels,frames);
    for(int c=0;c<channels;++c)for(int n=0;n<frames;++n)me.setSample(c,n,static_cast<float>(std::llround((200+c*80)*std::sin((n+c*11)*2*3.14159265358979323846/97)))/32768.f);
    for(int c=0;c<dialogueChannels;++c)for(int n=0;n<frames;++n)dialogue.setSample(c,n,static_cast<float>(std::llround((100+c*50)*std::sin((n+c*7)*2*3.14159265358979323846/53)))/32768.f);
    CompositionFixture f;f.recipe.preset=preset;f.recipe.me=component(root,"me",rate,me,grouped);f.recipe.dialogue=component(root,"dialogue",rate,dialogue,grouped);
    expect(f.recipe.confirmComplementaryBalancedAlignment(),"generator confirms complementary balanced aligned composition");
    f.expected.setSize(channels,frames);
    // Independent preset oracle: no production mapping or summing helper.
    for(int c=0;c<channels;++c)for(int n=0;n<frames;++n)
    {
        const bool absent=dialogueChannels==1&&c!=2;const int dc=dialogueChannels==1?0:c;
        const float result=absent?me.getSample(c,n):static_cast<float>(static_cast<double>(me.getSample(c,n))+static_cast<double>(dialogue.getSample(dc,n)));
        f.expected.setSample(c,n,result);
    }
    f.oracle=root.getChildFile("independent-composed-oracle.wav");floatFile(f.oracle,rate,f.expected);return f;
}
StandaloneSessionModel compositionModel(const ProgrammeComposition& recipe)
{
    StandaloneSessionModel s;expect(s.setValidatedComposition(recipe),"adopt explicit composition atomically");expect(s.initializeCurrentStateFromSource(),"initialize composed Current");
    s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);return s;
}
AnalysisSnapshot measuredComposition(StandaloneSessionModel& s,MeasurementDomain domain,ProcessedAudioObserver* capture=nullptr)
{
    expect(s.job.begin(JobStage::analysis,s.revision),"begin composition measurement");auto c=s.measurementContext(domain);c.ticket=s.job.ticket;JobControl job;
    auto result=OfflineAnalysisService(capture).analyze(c,domain,job);if(!result.value)std::cerr<<result.message<<'\n';expect(result.value.has_value(),"composition measurement completes");
    expect(s.job.finish(c.ticket,s.revision,result.error,result.message)&&s.acceptMeasurement(*result.value),"accept exact composition measurement");return *result.value;
}
bool sameSamples(const juce::AudioBuffer<float>& a,const juce::AudioBuffer<float>& b)
{
    if(a.getNumChannels()!=b.getNumChannels()||a.getNumSamples()!=b.getNumSamples())return false;
    for(int c=0;c<a.getNumChannels();++c)if(std::memcmp(a.getReadPointer(c),b.getReadPointer(c),static_cast<size_t>(a.getNumSamples())*4))return false;return true;
}
void coreCase(const juce::File& root,int rate,CompositionPreset preset,bool grouped,bool smoke)
{
    auto f=fixtureComposition(root,rate,preset,grouped,rate*4);auto s=compositionModel(f.recipe);const int channels=f.expected.getNumChannels(),frames=f.expected.getNumSamples();
    Capture sourceCapture(channels,frames);auto source=measuredComposition(s,MeasurementDomain::source,&sourceCapture);expect(sameSamples(f.expected,sourceCapture.data),"Source samples match independent composed oracle exactly");
    JobControl job;auto valid=validateMedia(f.oracle,job);expect(valid.succeeded(),"independent oracle imported");auto oracle=source.context;oracle.composition.reset();oracle.sourceFile=f.oracle;oracle.media=valid.media;
    auto analysis=OfflineAnalysisService().analyze(oracle,MeasurementDomain::source,job);expect(analysis.value.has_value(),"independent composed file measured");
    expect(source.measured.integratedLufs==analysis.value->measured.integratedLufs&&source.measured.loudnessRangeLu==analysis.value->measured.loudnessRangeLu&&source.detail.channelTruePeaks==analysis.value->detail.channelTruePeaks&&source.detail.truePeakEighthFrames==analysis.value->detail.truePeakEighthFrames,"composition metrics and per-channel true peaks equal independently composed file");
    auto diagnostic=proposal::analyzeProgramme(source.detail,source.identity),oracleDiagnostic=proposal::analyzeProgramme(analysis.value->detail,analysis.value->identity);
    expect(diagnostic&&oracleDiagnostic&&diagnostic->integratedRank==oracleDiagnostic->integratedRank&&diagnostic->excludedIntegrated==oracleDiagnostic->excludedIntegrated&&diagnostic->lraOrder==oracleDiagnostic->lraOrder&&diagnostic->frames==frames,"composition full-resolution diagnostic ranks match independent file and source duration");
    bool windows=true;for(size_t n=0;n<diagnostic->momentary.size();++n){const auto& a=diagnostic->momentary[n];const auto& b=oracleDiagnostic->momentary[n];windows&=a.begin==b.begin&&a.end==b.end&&a.energy==b.energy&&a.included==b.included&&a.padding==b.padding;}
    expect(windows,"composition diagnostic positions energy gate inclusion and padding match oracle");
    s.synchronizeSections();juce::String sectionError;expect(s.sectionDocument.split(1,frames/2,sectionError),"composed sections use unchanged audio-source frames");
    const auto map=s.sectionDocument.map();const auto sections=proposal::analyzeSections(*diagnostic,map,map.inputIdentity),oracleSections=proposal::analyzeSections(*oracleDiagnostic,map,map.inputIdentity);
    expect(sections&&oracleSections&&sections->sections.size()==2&&sections->sections[0].momentary.included==oracleSections->sections[0].momentary.included,"composition section summaries use same measured windows as independent file");
    auto prepared=prepareVerifiedSource(source.context,job);expect(prepared.value.has_value(),"composition private copies prepared");auto reader=(*prepared.value)->createReader();expect(reader!=nullptr,"shared composition reader available");
    for(int block:{1,127,4096,8193})
    {
        const int count=std::min(block,frames-11);juce::AudioBuffer<float> data(channels,count);expect(reader->read(data.getArrayOfWritePointers(),channels,11,count),"random access bounded composition read");
        bool same=true;for(int c=0;c<channels;++c)same&=std::memcmp(data.getReadPointer(c),f.expected.getReadPointer(c,11),count*4)==0;expect(same,"seek and chunk independent composition samples");
    }
    for(double gain:smoke?std::vector<double>{0}:std::vector<double>{0,3.125,-6.02059991328})
    {
        expect(s.setCurrentGainDb(gain),"set full precision global gain");Capture offline(channels,frames);const auto measured=measuredComposition(s,MeasurementDomain::currentState,&offline);
        juce::AudioBuffer<float> expected(f.expected);const double coefficient=std::pow(10.0,gain/20.0);
        for(int c=0;c<channels;++c)for(int n=0;n<frames;++n)expected.setSample(c,n,static_cast<float>(static_cast<double>(expected.getSample(c,n))*coefficient));
        expect(sameSamples(offline.data,expected),"all-channel global trim follows sum with exact existing precision");
        auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));
        s.selectAuditionDomain(MeasurementDomain::currentState);auto request=s.auditionRequest();Capture capture(channels,frames);audition.setCaptureObserver(&capture);
        expect(audition.connect(rate).isEmpty(),"injected device connection only");
        OutputRouting muted;audition.synchronize(request,muted);expect(audition.prepare(request).isEmpty(),"composed audition prepare with injected device");waitPreparation(audition);
        expect(audition.play(request,muted,true,true).isEmpty(),"explicit inaudible injected playback");awaitPlaying(audition,*injected);drain(audition,*injected,1024);
        expect(identical(offline,capture),"all channels remain in pre-monitor capture when every route muted");audition.disconnect();
        recordSamples(offline,rate,gain,"composition-"+juce::String(static_cast<int>(preset))+(grouped?"-grouped":"-interleaved"),measured.detail);
    }
    expect(s.setCurrentGainDb(-23-*source.measured.integratedLufs),"explicit test gain for export");measuredComposition(s,MeasurementDomain::currentState);
    expect(ex::readiness(ex::capture(s)).isEmpty(),"composed Current passes selected export rules");
    auto rendered=ex::render(ex::capture(s),root,job);if(!rendered.value)std::cerr<<rendered.message<<'\n';expect(rendered.value.has_value(),"composed render and actual artifact verification");
    const auto report=ex::report(ex::capture(s),rendered.value->get());expect(report.hasProperty("composition")&&report["composition"]["components"].size()==2,"report separately identifies accessible original components");
    expect(!(*rendered.value)->analysis().context.composition&&(*rendered.value)->analysis().context.processing.staticGainDb==0,"artifact is measured in its own domain without composition or second gain");
    const auto published=ex::publish(**rendered.value,ex::capture(s),root.getChildFile("export"),ex::approvalIdentity(**rendered.value,ex::capture(s)),job);expect(published.value.has_value(),"explicit complete composed export publication");
    auto doc=sf::capture(s);const auto file=root.getChildFile("composition.levelscope");expect(sf::save(doc,file,job).succeeded(),"coherent composed session save");
    auto read=sf::read(file,job);if(!read.document)std::cerr<<read.message<<'\n';expect(read.document.has_value()&&read.document->source.composition.has_value(),"schema-3 composition recovered");
    auto verified=sf::verify(*read.document,job);if(!verified.document)std::cerr<<verified.message<<'\n';expect(verified.document&&verified.document->audioResolved,"all reopened components verified");
    StandaloneSessionModel reopened;juce::String error;expect(sf::adopt(*verified.document,*verified.prepared,reopened,error),"safe composition session adoption");
    expect(reopened.source.hasMeasurement&&reopened.currentState.hasMeasurement&&reopened.currentStateDocument.staticGainDb==s.currentStateDocument.staticGainDb,"exact Current plan and valid compatible measured caches restored");
    expect(reopened.sourceDocument.composition->confirmation==f.recipe.confirmation,"exact composition confirmation persists");
}
}
#include "CompositionSafetyTests.h"
#include "CompositionEdgeTests.h"
#include "CompositionLongTests.h"
#include "CompositionIntegrationTests.h"
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI gui;bool smoke=false,passed=false;juce::String rootPath,output,message,phase="core",video,prores,exchange;
    for(int i=1;i<argc;++i){const juce::String a(argv[i]);if(a=="--smoke")smoke=true;else if(a=="--phase"&&i+1<argc)phase=argv[++i];else if(a=="--video"&&i+1<argc)video=argv[++i];else if(a=="--prores"&&i+1<argc)prores=argv[++i];else if(a=="--exchange"&&i+1<argc)exchange=argv[++i];else if(a=="--root"&&i+1<argc)rootPath=argv[++i];else if(a=="--output"&&i+1<argc)output=argv[++i];else if(a=="--samples"&&i+1<argc)sampleDirectory=juce::File(argv[++i]);else{std::cerr<<"Unknown/incomplete argument\n";return 2;}}
    try
    {
        expect(rootPath.isNotEmpty(),"explicit owned test root required");juce::File root(rootPath);expect(!root.exists()&&root.createDirectory().wasOk(),"new evidence directory required");
        if(sampleDirectory!=juce::File())expect(sampleDirectory.createDirectory().wasOk(),"capture directory ready");
        if(smoke)coreCase(root.getChildFile("smoke"),48000,CompositionPreset::surroundMonoCentre,false,true);
        else if(phase=="long")longComposition(root);
        else if(phase=="video"){expect(video.isNotEmpty()&&prores.isNotEmpty(),"explicit native fixture paths required");compositionVideo(root,juce::File(video),juce::File(prores));}
        else if(phase=="exchange"){expect(exchange.isNotEmpty(),"explicit foreign composition path required");compositionExchange(juce::File(exchange));}
        else if(phase=="core"||phase=="safety")
        {
            compositionSafety(root.getChildFile("safety"));compositionOverload(root.getChildFile("overload"));compositionEdges(root.getChildFile("edges"));compositionFormats(root.getChildFile("formats"));compositionMutedPeaks(root.getChildFile("muted-peaks"));
            if(phase=="core")for(int rate:{44100,48000,96000,192000})for(int p=0;p<3;++p)for(bool grouped:{false,true})
                coreCase(root.getChildFile(juce::String(rate)+"-"+juce::String(p)+(grouped?"-groups":"-files")),rate,static_cast<CompositionPreset>(p),grouped,false);
        }
        else throw std::runtime_error("Unknown phase");
        passed=true;
    }
    catch(const std::exception& e){message=e.what();std::cerr<<message<<'\n';}
    auto result=juce::var(new juce::DynamicObject());auto* r=result.getDynamicObject();r->setProperty("suite","M11 aligned composition");r->setProperty("passed",passed);r->setProperty("phase",smoke?"development-smoke":phase);r->setProperty("message",message);r->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);r->setProperty("checks",checks);r->setProperty("cases",cases);r->setProperty("physical_devices_opened",false);r->setProperty("long_composition",longEvidence);r->setProperty("video_cases",compositionVideoCases);
#if defined(__aarch64__)
    r->setProperty("architecture","arm64");
#else
    r->setProperty("architecture","x86_64");
#endif
    r->setProperty("video_selections",compositionVideoEvidence);
    r->setProperty("exchange_cases",compositionExchangeEvidence);
    if(output.isNotEmpty())juce::File(output).replaceWithText(juce::JSON::toString(result,true));std::cout<<checks.size()<<" assertions; "<<(passed?"PASS":"FAIL")<<'\n';return passed?0:1;
}

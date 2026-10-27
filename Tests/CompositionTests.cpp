#include "CompositionTestFixtures.h"
namespace
{
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

#include "AuditionTestSupport.h"
#include "Standalone/StandaloneExport.h"
#include "Standalone/Pcm24Wave.h"
#include <fstream>
#include <fcntl.h>
#include <unistd.h>
#include "Standalone/SessionDocument.h"
namespace ex=levelscope::standalone::exporting;
namespace pcm=levelscope::standalone::pcm24;
#include "ExportNumericTests.h"
namespace
{
void numericTests()
{
    int32_t q=0;
    expect(pcm::quantize(0,0,q)&&q==0,"PCM24 zero without test dither");
    expect(pcm::quantize(0,.5,q)&&q==0,"ties to even zero");
    expect(pcm::quantize(1.f/8388608.f,.5,q)&&q==2,"ties to even positive");
    expect(pcm::quantize(-1.f/8388608.f,-.5,q)&&q==-2,"ties to even negative");
    expect(!pcm::quantize(1,0,q)&&!pcm::quantize(-1,-.75,q),"no hidden full-scale clipping");
    expect(!pcm::quantize(std::numeric_limits<float>::quiet_NaN(),0,q),"non-finite output rejected");
    double sum=0,energy=0,cross=0;bool bounded=true;
    for(uint64_t n=0;n<100000;++n)
    {const auto a=pcm::dither(123,n,0),b=pcm::dither(123,n,1);bounded=bounded&&a>-1&&a<1;sum+=a;energy+=a*a;cross+=a*b;}
    expect(bounded,"bounded TPDF");
    expect(std::abs(sum/100000)<.005&&std::abs(energy/100000-1.0/6)<.005&&std::abs(cross/100000)<.005,"TPDF zero mean variance and channel decorrelation");
    pcm::Format f;juce::String error;
    expect(pcm::describe(48000,6,10800LL*48000,f,error)&&f.rf64&&f.dataBytes==9331200000ULL,"three-hour 5.1 RF64 payload");
    expect(pcm::header(f).size()==104,"bounded RF64 header");
    expect(pcm::describe(48000,2,48000,f,error)&&!f.rf64&&pcm::header(f).size()==68,"ordinary extensible RIFF");
    expect(!pcm::describe(48000,8,48000,f,error),"7.1 export deliberately unavailable");
    for(int rate:{44100,48000,96000,192000})for(int channels:{2,6})
    {
        const uint64_t frames=uint64_t(rate)*10800,payload=frames*channels*3;
        expect(pcm::describe(rate,channels,static_cast<int64_t>(frames),f,error)
            &&f.dataBytes==payload&&f.rf64==(payload+60>UINT64_C(0xfffffffe))
            &&f.fileBytes==payload+(f.rf64?104:68)&&pcm::header(f).size()==f.dataOffset,"three-hour size arithmetic through 192kHz, without claiming untested long-run rate performance");
    }
}
void measure(StandaloneSessionModel& s,MeasurementDomain domain,ProcessedAudioObserver* observer=nullptr)
{
    expect(s.job.begin(JobStage::analysis,s.revision),"begin test measurement");auto c=s.measurementContext(domain);c.ticket=s.job.ticket;JobControl control;
    auto a=OfflineAnalysisService(observer).analyze(c,domain,control);expect(a.value.has_value(),"test source/current measured");
    expect(s.job.finish(c.ticket,s.revision,a.error,a.message)&&s.acceptMeasurement(std::move(*a.value)),"accept exact measured test state");
}
StandaloneSessionModel measured(const juce::File& root,int channels,int rate)
{
    auto fixtureRequest=fixture(root.getChildFile("source.wav"),rate,channels,rate*4);
    StandaloneSessionModel s;expect(s.setValidatedSource(fixtureRequest.context.sourceFile,fixtureRequest.context.media),"adopt test programme");
    s.sourceDocument.layout=fixtureRequest.context.fileLayout;s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
    measure(s,MeasurementDomain::source);expect(s.initializeCurrentStateFromSource(),"initialize Current");
    expect(s.setCurrentGainDb(-23-*s.source.integratedLufs),"set exact normalizing trim");measure(s,MeasurementDomain::currentState);return s;
}
void workflowTests()
{
    Temp temp;auto s=measured(temp.directory,2,48000);auto r=ex::capture(s);
    expect(ex::readiness(r).isEmpty(),"valid measured Current ready for render");
    JobControl control;auto result=ex::render(r,temp.directory,control);
    if(!result.value)std::cerr<<result.message<<'\n';expect(result.value.has_value(),"render finalized PCM and verify");auto candidate=*result.value;
    expect(ex::artifactEvaluation(*candidate,r).verdict==proposal::Verdict::pass,"artifact has independent passing evaluation");
    expect(candidate->analysis().context.processing.staticGainDb==0&&candidate->analysis().domain==MeasurementDomain::renderedOutput,"encoded artifact measured without additional gain");
    const auto leased=candidate->verifiedSource();pcm::Hooks reuseSpace;reuseSpace.availableBytes=[&](const juce::File& dir){return dir==temp.directory?int64_t(1024)*1024*1024:0;};
    auto repeated=ex::render(r,temp.directory,control,reuseSpace,nullptr,leased);expect(repeated.value&&(*repeated.value)->verifiedSource()==leased&&(*repeated.value)->hash()==candidate->hash(),"repeated render reuses verified source, seed and exact PCM bytes without another copy");
    repeated.value.reset();
    const auto model=ex::report(r,candidate.get());const auto plain=ex::report(r);
    expect(model.hasProperty("artifact")&&!plain.hasProperty("artifact"),"measurement-only report never invents artifact");
    expect(ex::reportHtml(plain).contains("MEASUREMENT-ONLY"),"plain report clearly labelled");
    auto dest=temp.directory.getChildFile("export");auto published=ex::publish(*candidate,r,dest,ex::approvalIdentity(*candidate,r),control);
    if(!published.value)std::cerr<<published.message<<'\n';expect(published.value.has_value(),"explicit publication succeeds");
    expect(dest.getChildFile("programme.wav").existsAsFile()&&dest.getChildFile("manifest.json").existsAsFile(),"complete audio/report bundle");
    expect(!ex::publish(*candidate,r,dest,"",control).value,"cannot republish or overwrite");
    s.setCurrentGainDb(-55);auto unavailable=ex::capture(s);expect(ex::readiness(unavailable).isNotEmpty(),"unverified Current cannot render");
    auto report=ex::measurementReport(unavailable,temp.directory.getChildFile("measurement"),control);
    expect(report.value.has_value(),"measurement report available despite blocked export");
    expect(!report.value->getChildFile("programme.wav").exists(),"measurement-only publishes no audio");
    s.setCurrentGainDb(-10);measure(s,MeasurementDomain::currentState);auto failing=ex::capture(s);expect(ex::readiness(failing).isNotEmpty(),"measured failing Current cannot render");
    auto failed=ex::measurementReport(failing,temp.directory.getChildFile("measured-failure"),control);
    expect(failed.value&&juce::JSON::parse(failed.value->getChildFile("report.json"))["current"]["evaluation"]["verdict"].toString()=="FAIL","measurement-only report retains actual measured failure before rendering");
}
}
#include "ExportFailureTests.h"
#include "ExportLongTests.h"
#include "ExportVideoTests.h"
#include "ExportCrashTests.h"
#include "ExportMutedPeakTests.h"
#include "ExportExchangeTests.h"
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI gui;bool passed=false;juce::String message,output,artifacts,longSession,longDestination,video,prores,clockFixture,crashFolder,crashStage,exchange,rf64Audition;
    for(int i=1;i<argc;++i){const juce::String arg(argv[i]);if(i+1<argc&&arg=="--output")output=argv[++i];else if(i+1<argc&&arg=="--artifacts")artifacts=argv[++i];else if(i+1<argc&&arg=="--long-session")longSession=argv[++i];else if(i+1<argc&&arg=="--long-destination")longDestination=argv[++i];else if(i+1<argc&&arg=="--video")video=argv[++i];else if(i+1<argc&&arg=="--prores")prores=argv[++i];else if(i+1<argc&&arg=="--prepare-clock")clockFixture=argv[++i];else if(i+1<argc&&arg=="--crash-folder")crashFolder=argv[++i];else if(i+1<argc&&arg=="--rf64-audition")rf64Audition=argv[++i];else if(i+1<argc&&arg=="--exchange")exchange=argv[++i];else if(i+1<argc&&arg=="--crash-stage")crashStage=argv[++i];else if(!arg.startsWith("--"))output=arg;}
    if(crashFolder.isNotEmpty())crashCase(juce::File(crashFolder),crashStage);
    Temp temporary;const auto root=artifacts.isNotEmpty()?juce::File(artifacts):temporary.directory.getChildFile("artifacts");
    try{expect(!root.exists()&&root.createDirectory().wasOk(),"fresh evidence artifact directory");if(rf64Audition.isNotEmpty())auditionLargeArtifact(juce::File(rf64Audition));else if(exchange.isNotEmpty())exchangeExports(juce::File(exchange));else if(clockFixture.isNotEmpty())prepareExportClock(juce::File(clockFixture),juce::File(longDestination));else if(longSession.isNotEmpty())longExport(juce::File(longSession),juce::File(longDestination));else{numericTests();rawEncoderMatrix(root);containerFailures(root);workflowTests();failureTests();crashTests();exportMutedPeaks();mappingAndParity(root);if(video.isNotEmpty()&&prores.isNotEmpty())exportVideo(root.getChildFile("verified-6-48000/programme.wav"),juce::File(video),juce::File(prores));}passed=true;}catch(const std::exception& e){message=e.what();std::cerr<<e.what()<<'\n';}
    auto* result=new juce::DynamicObject();result->setProperty("suite","M10 verified PCM export and reporting");result->setProperty("passed",passed);result->setProperty("message",message);result->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);result->setProperty("checks",checks);result->setProperty("cases",cases);result->setProperty("long_export",longEvidence);result->setProperty("video_cases",exportVideoCases);result->setProperty("physical_devices_opened",false);
#if defined(__aarch64__)
    result->setProperty("architecture","arm64");
#else
    result->setProperty("architecture","x86_64");
#endif
    const auto json=ex::reportJson(juce::var(result));if(output.isNotEmpty())juce::File(output).replaceWithText(json);else std::cout<<json;
    return passed?0:1;
}

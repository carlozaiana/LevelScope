// Close the original long driver's identity-only comparison gap without
// repeating its proposal, Current, audition or export workflow.
#define main unusedWaveformQualificationMain
#include "WaveformTests.cpp"
#undef main
#include "Standalone/SessionDocument.h"

int main(int argc,char** argv)
{
    if(argc!=3)return 2;
    juce::ScopedJuceInitialiser_GUI gui;
    const juce::File session(argv[1]),output(argv[2]);
    auto report=juce::var(new juce::DynamicObject());auto* r=report.getDynamicObject();
    bool passed=false;juce::String failure;const auto begin=juce::Time::getMillisecondCounterHiRes();
    try
    {
        expect(!output.exists(),"fresh exact-long comparison receipt");
        JobControl job;
        auto read=sessionfile::read(session,job);
        expect(read.document&&read.document->sourceCache,"completed long session retains its observed Source cache");
        auto verified=sessionfile::verify(*read.document,job);
        expect(verified.prepared&&verified.prepared->model.source.hasMeasurement,"source bytes and saved observed cache reverified on the same host");
        auto& model=verified.prepared->model;
        const auto expected=model.source.analysis;
        expect(expected&&expected->detail.frames==2073600000LL&&expected->detail.sampleRate==192000,"actual three-hour 192 kHz Source result selected");
        expect(model.job.begin(JobStage::analysis,model.revision),"new independent plain Source job");
        auto context=model.measurementContext(MeasurementDomain::source);context.ticket=model.job.ticket;
        const auto measuredAt=juce::Time::getMillisecondCounterHiRes();
        auto measured=OfflineAnalysisService().analyze(context,MeasurementDomain::source,job);
        r->setProperty("plain_measurement_ms",juce::Time::getMillisecondCounterHiRes()-measuredAt);
        expect(measured.value&&model.job.finish(context.ticket,model.revision,measured.error,measured.message),"independent plain Source measurement completes");
        // Existing full oracle: exact metric values/status, all gates/counts,
        // percentiles, coverage, channel peak values/positions and every window.
        sameMeasurement(*expected,*measured.value);
        r->setProperty("source_sha256",expected->measuredContentHash);
        r->setProperty("measurement_identity",expected->identity);
        r->setProperty("momentary_windows",int(expected->detail.momentary.size()));
        r->setProperty("short_term_windows",int(expected->detail.shortTerm.size()));
        r->setProperty("lra_windows",int(expected->detail.lraWindows.size()));
        r->setProperty("session_sha256",juce::SHA256(session).toHexString());
        passed=true;
    }
    catch(const std::exception& e){failure=e.what();}
    r->setProperty("passed",passed);r->setProperty("failure",failure);r->setProperty("checks",checks);
    r->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);r->setProperty("os",juce::SystemStats::getOperatingSystemName());
    r->setProperty("total_ms",juce::Time::getMillisecondCounterHiRes()-begin);
    r->setProperty("scope","One fresh plain Source pass compared exactly with the verified observed Source cache from the completed three-hour workflow; existing fixture reused; no new proposal/export or physical output.");
    output.replaceWithText(juce::JSON::toString(report,true));std::cout<<juce::JSON::toString(report,true)<<std::endl;
    return passed?0:1;
}

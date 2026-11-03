// Supplemental synthetic display-density check; this does not decode long PCM.
// Built with -fno-access-control only in this diagnostic executable. Seed sealed
// finest-bin extrema, then use the unchanged production tree/reduction/painting.
#define main unusedDiagnosticsQualificationMain
#include "DiagnosticsTests.cpp"
#undef main
#include "Core/Waveform/WaveformOverview.h"

namespace {
waveform::Extrema densityBin(const waveform::Overview& wave,size_t n,int channel,int domain)
{
    const auto begin=int64_t(n)*wave.binFrames,end=std::min(wave.spec.frames,begin+wave.binFrames);
    return {-.05f*float((n+3*channel+domain)%17+1),.04f*float((2*n+channel+domain)%23+1),
            begin+std::min<int64_t>(channel*7+5,end-begin-1),begin+std::min<int64_t>(channel*13+9,end-begin-1)};
}
std::shared_ptr<waveform::Overview> densityOverview(int channels,int domain,const std::shared_ptr<waveform::Budget>& budget)
{
    auto wave=waveform::Overview::create({10800LL*192000,192000,channels},budget);
    expect(wave!=nullptr,"synthetic long overview allocation");
    while(wave->consumed<wave->spec.frames)
    {
        for(int c=0;c<channels;++c)wave->pending[size_t(c)]=densityBin(*wave,wave->sealed,c,domain);
        wave->consumed=std::min(wave->spec.frames,wave->consumed+wave->binFrames);
        wave->seal();
    }
    expect(wave->finish()&&wave->coveredFrames()==wave->spec.frames,"synthetic sealed bins cover the whole three-hour interval");
    return wave;
}
void densityOracle(const waveform::Overview& wave,int domain)
{
    for(int columns:{1,1200,4096})
    {
        const auto data=wave.reduce(0,wave.spec.frames,columns);
        expect(data.size()==size_t(columns),"long multichannel display columns bounded");
        bool exact=true;
        for(const auto& col:data)for(int c=0;c<wave.spec.channels;++c)
        {
            auto truth=densityBin(wave,size_t(col.begin/wave.binFrames),c,domain);
            for(auto n=col.begin/wave.binFrames+1;n<(col.end+wave.binFrames-1)/wave.binFrames;++n)
            {
                const auto x=densityBin(wave,size_t(n),c,domain);
                if(x.minimum<truth.minimum){truth.minimum=x.minimum;truth.minimumFrame=x.minimumFrame;}
                if(x.maximum>truth.maximum){truth.maximum=x.maximum;truth.maximumFrame=x.maximumFrame;}
            }
            const auto& actual=col.channels[size_t(c)];
            exact &= actual.minimum==truth.minimum&&actual.maximum==truth.maximum
                &&actual.minimumFrame==truth.minimumFrame&&actual.maximumFrame==truth.maximumFrame;
        }
        expect(exact,"every reduced synthetic channel preserves both extrema and earliest exact positions");
    }
}
juce::var densityCase(int channels)
{
    auto budget=std::make_shared<waveform::Budget>();
    auto source=densityOverview(channels,0,budget),current=densityOverview(channels,1,budget);
    expect(budget->used()<=waveform::combinedBudgetBytes,"both synthetic overview domains stay within 128 MiB");
    densityOracle(*source,0);densityOracle(*current,1);
    auto measured=longFixture();
    while(measured.channelSamplePeaks.size()<size_t(channels))
    {
        measured.channelSamplePeaks.push_back(.5);measured.channelTruePeaks.push_back(.51);
        measured.samplePeakFrames.push_back(measured.frames-2);measured.truePeakEighthFrames.push_back((measured.frames-2)*8+3);
    }
    const auto analysis=analyzeProgramme(measured,"synthetic-long-density-no-PCM-claim");
    MeasurementTimeline timeline;timeline.setSize(1200,760);timeline.setAnalysis(analysis);
    juce::StringArray roles{"L","R","C","LFE","Ls","Rs"};if(channels==8){roles.add("Lrs");roles.add("Rrs");}
    std::vector<double> preparation,paint;
    for(int n=0;n<100;++n)
    {
        const double begin=n%4?double(n*7919%1000)/1000*measured.frames:0;
        const double end=n%4?std::min(double(measured.frames),begin+192000.*(n%2?30:1800)):double(measured.frames);
        auto t=now();timeline.setWaveform(n%2?source:current,true,roles);timeline.setViewRange(begin,end);preparation.push_back(now()-t);
        t=now();@autoreleasepool{auto image=timeline.createComponentSnapshot(timeline.getLocalBounds());}paint.push_back(now()-t);
        expect(timeline.displayedWaveform().size()<=4096,"synthetic viewport retains bounded waveform columns");
    }
    expect(timeline.inspectSeconds(10799.9)&&timeline.inspectionText().contains("frames ["),"synthetic long multichannel EOF remains inspectable");
    auto result=juce::var(new juce::DynamicObject());auto* r=result.getDynamicObject();
    r->setProperty("channels",channels);r->setProperty("frames",juce::String(measured.frames));r->setProperty("sample_rate",192000);
    r->setProperty("combined_overview_bytes",juce::String(budget->used()));r->setProperty("finest_interval_frames",juce::String(source->binFrames));
    r->setProperty("prepare_p95_ms",p95(preparation));r->setProperty("paint_p95_ms",p95(paint));
    return result;
}
}
int main(int argc,char** argv)
{
    if(argc!=2)return 2;juce::ScopedJuceInitialiser_GUI gui;juce::Array<juce::var> rows;bool passed=false;juce::String failure;
    try{for(int channels:{6,8}){const auto row=densityCase(channels);rows.add(row);
        expect(double(row["prepare_p95_ms"])<=50&&double(row["paint_p95_ms"])<=100,"synthetic long multichannel preparation and paint meet agreed p95 budgets");}passed=true;}
    catch(const std::exception& e){failure=e.what();}
    auto result=juce::var(new juce::DynamicObject());auto* r=result.getDynamicObject();
    r->setProperty("passed",passed);r->setProperty("failure",failure);r->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);
    r->setProperty("os",juce::SystemStats::getOperatingSystemName());r->setProperty("checks",checks);r->setProperty("cases",rows);
    r->setProperty("scope","Synthetic three-hour 6/8-channel sealed extrema and measurement-sized data; real production reduction and offscreen painting. No long multichannel PCM or hardware playback claim.");
    juce::File(argv[1]).replaceWithText(juce::JSON::toString(result,true));std::cout<<juce::JSON::toString(result,true)<<std::endl;
    return passed?0:1;
}

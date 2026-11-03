#include "CompositionTestFixtures.h"
#include "Standalone/StandaloneWaveform.h"
#include "Standalone/MeasurementTimeline.h"
#include <future>

namespace
{
double observedMs=0,plainMs=0;
void sameMeasurement(const AnalysisSnapshot& a,const AnalysisSnapshot& b)
{
    expect(a.identity==b.identity && a.measuredContentHash==b.measuredContentHash,"preview preserves measurement identity");
    const auto& x=a.detail;const auto& y=b.detail;
    auto metric=[](const measurement::Metric& p,const measurement::Metric& q){return p.value==q.value&&p.status==q.status;};
    expect(metric(x.integrated,y.integrated)&&metric(x.lra,y.lra)&&metric(x.momentaryMax,y.momentaryMax)
        &&metric(x.shortTermMax,y.shortTermMax)&&metric(x.samplePeak,y.samplePeak)&&metric(x.truePeak,y.truePeak),"preview preserves every metric and availability bit exactly");
    expect(x.frames==y.frames&&x.sampleRate==y.sampleRate&&x.lraStable==y.lraStable
        &&x.integratedRelativeGate==y.integratedRelativeGate&&x.lraRelativeGate==y.lraRelativeGate
        &&x.p10==y.p10&&x.p95==y.p95&&x.integratedIncluded==y.integratedIncluded&&x.lraIncluded==y.lraIncluded,"preview preserves gates, populations, percentiles and coverage");
    expect(x.channelSamplePeaks==y.channelSamplePeaks&&x.channelTruePeaks==y.channelTruePeaks
        &&x.samplePeakFrames==y.samplePeakFrames&&x.truePeakEighthFrames==y.truePeakEighthFrames,"preview preserves per-channel peak values and positions");
    auto windows=[](const auto& p,const auto& q){if(p.size()!=q.size())return false;for(size_t n=0;n<p.size();++n)if(p[n].endFrame!=q[n].endFrame||p[n].energy!=q[n].energy)return false;return true;};
    expect(windows(x.momentary,y.momentary)&&windows(x.shortTerm,y.shortTerm)&&windows(x.lraWindows,y.lraWindows),"preview preserves every M, S and LRA measured window exactly");
}
void exactPreview(const LoudnessProgress& progress,const WaveformBinding& binding,const measurement::Result& result)
{
    auto data=progress.view(binding);expect(data&&data->coverage()==result.frames,"measured loudness preview covers actual EOF");
    for(bool s:{false,true})
    {
        bool equal=true;for(const auto& w:(s?result.shortTerm:result.momentary))
        {auto point=data->inspect(s,w.endFrame,data->coverage());equal=equal&&point&&point->endFrame==w.endFrame&&point->energy==w.energy;}
        expect(equal,"all progressive M/S window ends and energies exactly equal the accepted measurement");
    }
}
void oracle(const waveform::Overview& data,const juce::AudioBuffer<float>& audio)
{
    expect(data.complete()&&data.spec.frames==audio.getNumSamples()&&data.spec.channels==audio.getNumChannels(),"waveform covers exactly the actual programme and ordered channels");
    for(int width:{1,17,4096})
    {
        const auto columns=data.reduce(0,data.spec.frames,width);bool exact=true;
        for(const auto& column:columns)for(int c=0;c<data.spec.channels;++c)
        {
            waveform::Extrema expected{audio.getSample(c,int(column.begin)),audio.getSample(c,int(column.begin)),column.begin,column.begin};
            for(auto n=column.begin+1;n<column.end;++n)
            {
                const float value=audio.getSample(c,int(n));
                if(value<expected.minimum){expected.minimum=value;expected.minimumFrame=n;}
                if(value>expected.maximum){expected.maximum=value;expected.maximumFrame=n;}
            }
            const auto& e=column.channels[size_t(c)];
            exact=exact&&column.coveredEnd==column.end&&e.minimum==expected.minimum&&e.maximum==expected.maximum
                &&e.minimumFrame==expected.minimumFrame&&e.maximumFrame==expected.maximumFrame;
        }
        expect(exact,"display reduction preserves exact extrema and earliest sample positions from independent PCM oracle");
    }
}
StandaloneSessionModel model(const juce::File& root,int rate,int channels,int kind,juce::AudioBuffer<float>& raw)
{
    const int frames=rate*3+101;StandaloneSessionModel s;
    expect(root.createDirectory().wasOk(),"fixture directory");
    if(kind==2)
    {
        auto f=fixtureComposition(root,rate,channels==2?CompositionPreset::stereoStereo:CompositionPreset::surroundMonoCentre,true,frames);
        raw=f.expected;expect(s.setValidatedComposition(f.recipe),"confirmed composed programme");
    }
    else
    {
        raw.setSize(channels,frames);
        for(int c=0;c<channels;++c)for(int n=0;n<frames;++n)
        {
            float value=float(.15*std::sin(n*.017));
            if(c==1)value=-value;else if(c>1)value=float((c+1)*.02*std::sin(n*(c+1)*.011));
            if(n==frames-1)value= c%2 ? -1.25f:1.5f;
            raw.setSample(c,n,value);
        }
        auto p=component(root,"mix",rate,raw,kind==1);
        if(p.assembly)expect(s.setValidatedProgramme(*p.assembly),"confirmed aligned channel programme");
        else {expect(s.setValidatedSource(p.file,p.media),"confirmed interleaved programme");s.sourceDocument.layout=p.layout;}
    }
    return s;
}
void runCase(const juce::File& root,int rate,int channels,int kind)
{
    juce::AudioBuffer<float> raw;auto s=model(root,rate,channels,kind,raw);WaveformCache cache;
    expect(s.job.begin(JobStage::analysis,s.revision),"begin initial analysis");
    auto context=s.measurementContext(MeasurementDomain::source);auto binding=WaveformBinding::from(context,MeasurementDomain::source);
    auto preview=cache.begin(binding);LoudnessProgress sourceWindows(binding);Capture capture(channels,raw.getNumSamples());JobControl control;
    auto started=juce::Time::getMillisecondCounterHiRes();
    auto result=OfflineAnalysisService(&capture,preview.get(),&sourceWindows).analyze(context,MeasurementDomain::source,control);
    observedMs+=juce::Time::getMillisecondCounterHiRes()-started;
    expect(result.value.has_value()&&capture.valid,"Source analysis with waveform completes");
    exactPreview(sourceWindows,binding,result.value->detail);oracle(*preview->overview(),raw);expect(sameSamples(raw,capture.data),"Source waveform input is confirmed programme samples");
    started=juce::Time::getMillisecondCounterHiRes();
    auto plain=OfflineAnalysisService().analyze(context,MeasurementDomain::source,control);plainMs+=juce::Time::getMillisecondCounterHiRes()-started;expect(plain.value.has_value(),"observer-off Source analysis completes");sameMeasurement(*result.value,*plain.value);
    expect(s.job.finish(context.ticket,s.revision,result.error,result.message)&&s.acceptMeasurement(*result.value),"Source measurement accepted before waveform promotion");
    expect(cache.accept(binding,WaveformBinding::from(s.measurementContext(MeasurementDomain::source),MeasurementDomain::source)),"matching complete waveform promoted");
    expect(s.currentState.analysis&&cache.shareInitialCurrent(*s.source.analysis,*s.currentState.analysis),"tested automatic Current initialization shares Source waveform");
    auto currentBinding=WaveformBinding::from(s.measurementContext(MeasurementDomain::currentState),MeasurementDomain::currentState);
    expect(cache.view(binding).data==cache.view(currentBinding).data,"identical Source and Current own one overview allocation");
    auto wrongOrigin=*s.currentState.analysis;wrongOrigin.originJobId+=1;
    expect(!cache.shareInitialCurrent(*s.source.analysis,wrongOrigin),"incorrect Current measurement origin cannot authorize sharing");
    expect(s.setCurrentGainDb(2.375),"apply global gain");
    if(channels==2||channels==6)
    {
        s.synchronizeSections();juce::String error;expect(s.sectionDocument.split(1,raw.getNumSamples()/2,error),"create explicit section boundary");
        expect(s.beginSectionTrimDraft()&&s.setSectionDraftGain(0,-6)&&s.setSectionDraftGain(1,3)&&s.applySectionTrimDraft(),"apply exact gain transitions");
    }
    currentBinding=WaveformBinding::from(s.measurementContext(MeasurementDomain::currentState),MeasurementDomain::currentState);
    cache.synchronize(WaveformBinding::from(s.measurementContext(MeasurementDomain::source),MeasurementDomain::source),currentBinding,true,true);
    expect(!cache.view(currentBinding).data&&cache.view(binding).data&&s.currentStateDocument.needsRemeasurement,"processing edit invalidates only Current waveform and measurement");
    expect(s.job.begin(JobStage::analysis,s.revision),"begin independent Current analysis");
    context=s.measurementContext(MeasurementDomain::currentState);binding=WaveformBinding::from(context,MeasurementDomain::currentState);
    preview=cache.begin(binding);LoudnessProgress currentWindows(binding);Capture processed(channels,raw.getNumSamples());
    started=juce::Time::getMillisecondCounterHiRes();
    result=OfflineAnalysisService(&processed,preview.get(),&currentWindows).analyze(context,MeasurementDomain::currentState,control);observedMs+=juce::Time::getMillisecondCounterHiRes()-started;
    expect(result.value.has_value()&&processed.valid,"Current waveform observes the shared applied processor including transitions");exactPreview(currentWindows,binding,result.value->detail);oracle(*preview->overview(),processed.data);
    started=juce::Time::getMillisecondCounterHiRes();
    plain=OfflineAnalysisService().analyze(context,MeasurementDomain::currentState,control);plainMs+=juce::Time::getMillisecondCounterHiRes()-started;expect(plain.value.has_value(),"observer-off Current completes");sameMeasurement(*result.value,*plain.value);
    expect(s.job.finish(context.ticket,s.revision,result.error,result.message)&&s.acceptMeasurement(*result.value)&&cache.accept(binding,binding),"independently measured Current waveform accepted");
    const auto sourceSnapshot=s.source.analysis,currentSnapshot=s.currentState.analysis;
    expect(!cache.shareInitialCurrent(*sourceSnapshot,*currentSnapshot),"processed Current never shares Source extrema");
    expect(s.job.begin(JobStage::waveform,s.revision),"begin replaceable waveform-only job");
    context=s.measurementContext(MeasurementDomain::currentState);binding=WaveformBinding::from(context,MeasurementDomain::currentState);preview=cache.begin(binding);
    auto rebuilt=OfflineAnalysisService(nullptr,preview.get()).buildWaveform(context,MeasurementDomain::currentState,control);
    expect(rebuilt.value.has_value()&&s.job.finish(context.ticket,s.revision,rebuilt.error,rebuilt.message)&&cache.accept(binding,binding),"waveform-only rebuild accepted");
    oracle(**rebuilt.value,processed.data);
    expect(s.source.analysis==sourceSnapshot&&s.currentState.analysis==currentSnapshot&&!s.currentStateDocument.needsRemeasurement,"Build waveform preserves accepted snapshots, origins and verification");
    if(rate==48000&&channels==2&&kind==0)
    {
        expect(s.job.begin(JobStage::waveform,s.revision),"begin cancellable waveform-only rebuild");context=s.measurementContext(MeasurementDomain::currentState);
        binding=WaveformBinding::from(context,MeasurementDomain::currentState);preview=cache.begin(binding);JobControl cancelled;cancelled.cancel=true;
        auto stopped=OfflineAnalysisService(nullptr,preview.get()).buildWaveform(context,MeasurementDomain::currentState,cancelled);
        expect(!stopped.value&&stopped.error==ErrorCode::cancelled&&!cache.accept(binding,binding),"cancelled waveform-only rebuild cannot publish");cache.reject(binding);
        expect(s.source.analysis==sourceSnapshot&&s.currentState.analysis==currentSnapshot&&!s.currentStateDocument.needsRemeasurement,"cancelled waveform rebuild preserves accepted measurements");
    }
    expect(cache.retainedBytes()<=waveform::combinedBudgetBytes,"both live domains and temporarily held worker references remain within 128 MiB");
}
struct CancelObserver final : ProcessedAudioObserver
{
    JobControl& control;WaveformProgress& progress;bool sawPartial=false;
    CancelObserver(JobControl& c,WaveformProgress& p):control(c),progress(p){}
    void observe(const float* const*,int,int,int64_t first) noexcept override
    {
        if(first>0){auto data=progress.overview();sawPartial=data&&data->coveredFrames()>0&&!data->complete();control.cancel=true;}
    }
};
void lifecycle(const juce::File& root)
{
    juce::AudioBuffer<float> raw;auto s=model(root,48000,2,0,raw);auto c=s.measurementContext(MeasurementDomain::source);
    c.ticket={17,s.revision};auto binding=WaveformBinding::from(c,MeasurementDomain::source);WaveformCache cache;
    auto progress=cache.begin(binding);JobControl control;CancelObserver cancel(control,*progress);
    const auto result=OfflineAnalysisService(&cancel,progress.get()).analyze(c,MeasurementDomain::source,control);
    expect(cancel.sawPartial&&!result.value&&result.error==ErrorCode::cancelled,"partial actual PCM is visible during analysis and cancellation returns no final measurement");
    expect(!progress->overview()&&!cache.accept(binding,binding),"cancelled partial coverage cannot be promoted");cache.reject(binding);
    auto pool=std::make_shared<waveform::Budget>(waveform::displayReserveBytes);WaveformProgress unavailable(pool);JobControl fresh;
    auto missing=OfflineAnalysisService(nullptr,&unavailable).analyze(c,MeasurementDomain::source,fresh);
    expect(missing.value.has_value()&&unavailable.state()==WaveformProgress::State::unavailable,"waveform memory exhaustion does not fail measurement");
    auto plain=OfflineAnalysisService().analyze(c,MeasurementDomain::source,fresh);expect(plain.value.has_value(),"plain lifecycle measurement");sameMeasurement(*missing.value,*plain.value);
    auto fill=[&](std::shared_ptr<WaveformProgress> p){auto o=p->begin({raw.getNumSamples(),48000,2});expect(o&&o->append(raw.getArrayOfReadPointers(),2,raw.getNumSamples(),0),"synthetic completed publication");p->finish(o);};
    progress=cache.begin(binding);fill(progress);auto newer=binding;++newer.ticket.id;
    expect(!cache.accept(binding,newer),"late job cannot be accepted for another ticket");
    auto next=cache.begin(newer);cache.reject(binding);fill(next);expect(cache.accept(newer,newer),"rejecting old completion leaves new worker intact");
    auto changed=newer;changed.runtime="reopened-runtime";cache.synchronize(changed,changed,true,false);
    expect(!cache.view(newer).data,"reopen runtime retires cached waveform");
    progress=cache.begin(binding);fill(progress);changed=binding;++changed.ticket.sessionRevision;
    cache.synchronize(changed,changed,true,false);expect(progress->state()==WaveformProgress::State::abandoned&&!cache.accept(binding,binding),"revision change abandons active publication");
    progress=cache.begin(binding);fill(progress);cache.synchronize(binding,binding,false,false);
    expect(!cache.view(binding).data,"changed or missing source clears waveform coverage");
    auto identity=binding;auto modified=c;std::swap(modified.fileLayout.roles[0],modified.fileLayout.roles[1]);
    expect(!identity.sameAudio(WaveformBinding::from(modified,MeasurementDomain::source)),"ordered roles are part of waveform identity");
}
void timeline(const juce::File& root)
{
    auto budget=std::make_shared<waveform::Budget>();auto data=waveform::Overview::create({48000*4+13,48000,2},budget);
    juce::AudioBuffer<float> raw(2,int(data->spec.frames));raw.clear();
    for(int n=64;n<raw.getNumSamples();++n){raw.setSample(0,n,float((.25+.15*std::sin(n*.0001))*std::sin(n*.025)));raw.setSample(1,n,float((.3+.2*std::sin(n*.00007))*std::sin(n*.031)));}
    raw.setSample(0,7,-.625f);raw.setSample(0,41,.75f);raw.setSample(1,99,-1.25f);
    expect(data->append(raw.getArrayOfReadPointers(),2,4096,0),"first actual PCM block");
    MeasurementTimeline ui;ui.setSize(1200,650);ui.setWaveform(data,false,{"L","R"});
    expect(ui.inspectSeconds(0)&&ui.inspectionText().contains("[0, 64)")&&ui.inspectionText().contains("frame 7")&&ui.inspectionText().contains("frame 41"),"deepest inspection identifies represented interval and both exact sample positions");
    expect(ui.inspectionText().contains("ANALYSING")&&ui.inspectionText().contains("measurements pending"),"partial waveform is not a completed measurement");
    expect(ui.inspectSeconds(2)&&ui.inspectionText().contains("Pending coverage"),"unread region is pending, not silence");
    auto image=ui.createComponentSnapshot(ui.getLocalBounds());auto output=root.getChildFile("partial.png").createOutputStream();expect(output&&juce::PNGImageFormat().writeImageToStream(image,*output),"render partial waveform inspection artifact");output.reset();
    const auto token=data->publication();const float* remainder[]={raw.getReadPointer(0,4096),raw.getReadPointer(1,4096)};
    expect(data->append(remainder,2,raw.getNumSamples()-4096,4096)&&data->finish(),"complete progressive overview");
    auto frozen=data->reduce(0,data->spec.frames,1,token);expect(frozen[0].coveredEnd==4096&&data->coveredFrames(token)==4096,"prepared coverage snapshot stays consistent after worker advances");
    measurement::OfflineMeasurement meter(48000,{measurement::Role::left,measurement::Role::right});meter.push(raw.getArrayOfReadPointers(),raw.getNumSamples());
    const auto measured=proposal::analyzeProgramme(meter.finish(),"actual-independent-measurement");ui.setAnalysis(measured);ui.setWaveform(data,true,{"L","R"});
    expect(ui.inspectSeconds(.7)&&ui.inspectionText().contains("Measured loudness")&&!ui.inspectionText().contains("ANALYSING"),"final inspection uses independent measured loudness windows");
    expect(ui.displayedPlot().momentary.size()>0,"measured M curve remains present with waveform");
    ui.setWaveform(data,false,{"L","R"});
    expect(ui.inspectSeconds(.7)&&ui.inspectionText().contains("Measured loudness")&&!ui.inspectionText().contains("measurement pending"),"waveform rebuild keeps existing measured loudness available for inspection");
    ui.setWaveform(data,true,{"L","R"});
    ui.setViewRange(0,1);expect(ui.displayedWaveform().size()==1&&ui.displayedWaveform()[0].end==64,"overview stops at its actual bin resolution");
    ui.setViewRange(0,double(data->spec.frames));ui.setSize(10000,650);expect(ui.displayedWaveform().size()<=4096,"oversized display remains bounded");ui.setSize(1200,650);
    image=ui.createComponentSnapshot(ui.getLocalBounds());output=root.getChildFile("complete.png").createOutputStream();expect(output&&juce::PNGImageFormat().writeImageToStream(image,*output),"render complete waveform with measured loudness");output.reset();
    double worst=0;for(int n=0;n<100;++n){const auto start=juce::Time::getMillisecondCounterHiRes();ui.setViewRange(n*100,48000+n*100);image=ui.createComponentSnapshot(ui.getLocalBounds());worst=std::max(worst,juce::Time::getMillisecondCounterHiRes()-start);}
    auto* row=new juce::DynamicObject();row->setProperty("viewport_prepare_and_software_paint_max_ms",worst);cases.add(juce::var(row));
    expect(worst<100,"short development timeline preparation and software paint stays responsive");
    ui.setWaveform({},false,{});ui.setAnalysis({});expect(ui.displayedWaveform().empty()&&!ui.inspectSeconds(0),"clearing source removes all inspection data");
}
}
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI gui;juce::String rootPath,output;bool passed=false;juce::String failure;
    for(int n=1;n<argc;++n){const juce::String arg(argv[n]);if(arg=="--root"&&n+1<argc)rootPath=argv[++n];else if(arg=="--output"&&n+1<argc)output=argv[++n];else return 2;}
    try
    {
        expect(rootPath.isNotEmpty(),"explicit owned fixture root required");const juce::File root(rootPath);expect(!root.exists()&&root.createDirectory().wasOk(),"new owned fixture root");
        for(int rate:{44100,48000,96000,192000})for(int channels:{1,2,5,6,8})runCase(root.getChildFile(juce::String(rate)+"-"+juce::String(channels)),rate,channels,0);
        runCase(root.getChildFile("assembly"),48000,6,1);runCase(root.getChildFile("composition"),48000,6,2);
        lifecycle(root.getChildFile("lifecycle"));timeline(root);
        auto* timing=new juce::DynamicObject();timing->setProperty("observed_analysis_total_ms",observedMs);timing->setProperty("plain_analysis_total_ms",plainMs);timing->setProperty("observed_to_plain_ratio",observedMs/plainMs);cases.add(juce::var(timing));
        expect(observedMs<=plainMs*1.15,"aggregate short waveform observation overhead <=15 percent");passed=true;
    }
    catch(const std::exception& e){failure=e.what();std::cerr<<failure<<'\n';}
    auto* report=new juce::DynamicObject();report->setProperty("passed",passed);report->setProperty("failure",failure);
    report->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);report->setProperty("checks",juce::var(checks));report->setProperty("performance",juce::var(cases));
    if(output.isNotEmpty())juce::File(output).replaceWithText(juce::JSON::toString(juce::var(report),true));
    std::cout<<"M18 waveform "<<(passed?"PASS":"FAIL")<<" "<<checks.size()<<" checks\n";return passed?0:1;
}

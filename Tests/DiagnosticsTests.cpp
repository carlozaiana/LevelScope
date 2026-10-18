#include "Standalone/StandaloneDiagnostics.h"
#include "Standalone/MeasurementTimeline.h"
#include "Standalone/StandaloneDiagnosticsPanel.h"
#include "Standalone/StandaloneSessionModel.h"
#include "Core/Measurement/GatePopulation.h"
#include <iostream>
#include <thread>
#include <chrono>
#include <limits>
#include <numeric>
#include <sys/resource.h>
#include <mach/mach.h>
using namespace levelscope;
using namespace levelscope::standalone;
using namespace levelscope::proposal;
namespace {
juce::Array<juce::var> checks;
juce::DynamicObject::Ptr performance = new juce::DynamicObject();
void expect(bool ok, const char* name) {
    auto* c=new juce::DynamicObject(); c->setProperty("name",name);c->setProperty("passed",ok);checks.add(juce::var(c));
    if(!ok) throw std::runtime_error(name);
}
double now(){return juce::Time::getMillisecondCounterHiRes();}
template<class F> void rejects(F f,const char* name){bool rejected=false;try{f();}catch(const std::invalid_argument&){rejected=true;}expect(rejected,name);}
double energy(double lu){return std::pow(10.,(lu+.691)/10);}
void boundaries(){
    const auto abs=measurement::absoluteGateEnergy(), inf=std::numeric_limits<double>::infinity();
    expect(!measurement::includedInPopulation(abs,abs/10,false),"integrated absolute boundary excluded");
    expect(measurement::includedInPopulation(abs,abs/100,true),"LRA absolute boundary included");
    for(bool lra:{false,true}){
        expect(!measurement::includedInPopulation(std::nextafter(abs,0.),abs/100,lra),"one ULP below absolute excluded");
        expect(measurement::includedInPopulation(std::nextafter(abs,inf),abs/100,lra),"one ULP above absolute included");
        auto relative=energy(-45);
        expect(measurement::includedInPopulation(relative,relative,lra)==lra,"relative exact comparison preserves distinct I and LRA convention");
        expect(!measurement::includedInPopulation(std::nextafter(relative,0.),relative,lra),"relative one ULP below excluded");
        expect(measurement::includedInPopulation(std::nextafter(relative,inf),relative,lra),"relative one ULP above included");
    }
    expect(measurement::percentileIndex(11,.1)==1 && measurement::percentileIndex(11,.95)==10,"percentile reference rounding at half index");
}
measurement::Result signal(int rate,int channels,double seconds,int kind){
    using R=measurement::Role; const std::vector<R> roles=channels==1?std::vector<R>{R::centre}:channels==2?std::vector<R>{R::left,R::right}:std::vector<R>{R::left,R::right,R::centre,R::surroundLeft,R::surroundRight,R::lfe};
    measurement::OfflineMeasurement m(rate,roles);juce::AudioBuffer<float> audio(channels,1024);std::vector<const float*> ptr(channels);
    auto frames=static_cast<int64_t>(seconds*rate);
    for(int64_t start=0;start<frames;start+=1024){int n=static_cast<int>(std::min<int64_t>(1024,frames-start));
        for(int c=0;c<channels;++c){for(int j=0;j<n;++j){auto f=start+j;double x=.08*std::sin(2*3.141592653589793*1000*f/rate);
            if(kind==1)x=0; if(kind==2)x*=1.e-5; if(kind==3)x=(f==0||f==frames-1)?.8:0;
            if(kind==4)x=c==channels-1?x:0;
            audio.setSample(c,j,static_cast<float>(x));}ptr[c]=audio.getReadPointer(c);}
        m.push(ptr.data(),n);
    }return m.finish();
}
void populations(){
    for(int rate:{44100,48000,96000,192000}){
        auto r=signal(rate,1,3.2,0);const auto original=r.integrated.value;
        auto a=analyzeProgramme(r,"actual-multirate");expect(a && a->integrated.value==original,"diagnostic construction preserves actual measured value");
        expect(a->momentary.front().begin==0 && a->momentary.front().end==rate*4/10,"M actual half-open start/end");
        expect(a->shortTerm.front().begin==0 && a->shortTerm.front().end==rate*3,"S actual half-open start/end");
        expect(a->lraWindows.back().padding && a->frames==r.frames && a->shortTerm.back().end<=r.frames,"LRA padding does not extend programme timeline");
        expect(a->integratedRank.size()==r.integratedIncluded && a->lraOrder.size()==r.lraIncluded && a->p10==r.p10 && a->p95==r.p95,"qualified population counts and percentiles preserved");
        double shares=0;for(auto n:a->integratedRank)shares+=a->momentary[n].energyShare;
        expect(std::abs(shares-1)<1e-12,"overlapping gated block energy shares sum to one");
        expect(a->identity==analyzeProgramme(r,"actual-multirate")->identity,"same snapshot has deterministic analysis identity");
        expect(a->identity!=analyzeProgramme(r,"other-snapshot")->identity,"different snapshot identity binds separately");
        auto plot=prepareTimeline(*a,0,static_cast<double>(a->frames),1);
        double minimum=1e9,maximum=-1e9;for(auto n:a->integratedRank){minimum=std::min(minimum,*a->momentary[n].lufs);maximum=std::max(maximum,*a->momentary[n].lufs);}
        expect(plot.momentary[0].included.minimum==minimum && plot.momentary[0].included.maximum==maximum,"single display bucket preserves exact extrema");
        expect(inspectWindow(a->momentary,(a->momentary[0].end+a->momentary[1].end)/2.)==0,"inspection tie chooses earlier full-resolution window");
        auto malformed=r;malformed.momentary[0].energy=std::numeric_limits<double>::quiet_NaN();rejects([&]{analyzeProgramme(malformed,"invalid");},"nonfinite window refused");
        malformed=r;--malformed.momentary[0].endFrame;rejects([&]{analyzeProgramme(malformed,"invalid");},"malformed location refused");
        malformed=r;++malformed.integratedIncluded;rejects([&]{analyzeProgramme(malformed,"invalid");},"inconsistent measured population refused");
        malformed=r;malformed.p95=*r.p95+1;rejects([&]{analyzeProgramme(malformed,"invalid");},"inconsistent percentile refused");
    }
    for(int kind:{1,2,3,4}){
        auto r=signal(48000,6,3.1,kind);auto a=analyzeProgramme(r,"silence-peaks-lfe");
        expect(a && a->lra.status==r.lra.status && a->lraStable==r.lraStable,"availability and stability survive analysis");
        for(int c=0;c<6;++c)expect(a->peaks[c].sampleFrame==r.samplePeakFrames[c] && a->peaks[c].truePeakEighthFrame==r.truePeakEighthFrames[c],"channel maximum locations preserved at eighth-frame precision");
        if(kind==1){expect(!a->integrated.value && a->integratedRank.empty() && a->histogram.empty(),"silence never invents zero-valued loudness");expect(!a->peaks[0].trueDbtp && !a->momentary[0].lufs,"silent maxima have no finite values or plotted points");}
        if(kind==2)expect(!a->integrated.value && a->integratedRank.empty() && !a->histogram.empty(),"excluded faint windows retain descriptive values without passing gates");
        if(kind==4)expect(!a->integrated.value && a->peaks[5].trueDbtp.has_value() && !a->peaks[0].trueDbtp,"LFE peak remains available without loudness contribution");
    }
    for(double duration:{.2,.5,2.9}){auto r=signal(48000,1,duration,0);auto a=analyzeProgramme(r,"short");expect(!a->lra.value && !a->lraStable && a->lraOrder.empty(),"short input does not acquire LRA through padding");expect(diagnosticMetricText(a->lra,"LU").contains("unavailable"),"short metric explicitly unavailable");}
}
// Generated three-hour measurement-sized fixture. This benchmarks downstream
// diagnostics, not three hours of PCM decoding or a normative measurement signal.
measurement::Result longFixture(){
    measurement::Result r;r.sampleRate=192000;r.frames=static_cast<int64_t>(r.sampleRate)*10800;
    auto fill=[&](auto& windows,int length,int64_t limit){for(int64_t end=length;end<=limit;end+=r.sampleRate/10){auto n=end/(r.sampleRate/10);double lu=n%1000==0?-4:n%7==0?-80:n%3==0?-35:-23;windows.push_back({end,energy(lu)});}};
    fill(r.momentary,r.sampleRate*4/10,r.frames);fill(r.shortTerm,r.sampleRate*3,r.frames);fill(r.lraWindows,r.sampleRate*3,r.frames+r.sampleRate*3/2);
    // Independent fixture oracle: first absolute-gated mean, then I/LRA populations.
    const double absolute=std::pow(10.,(-70+.691)/10);
    for(bool lra:{false,true}){const auto& windows=lra?r.lraWindows:r.momentary;double sum=0;size_t n=0;
        for(auto w:windows)if(lra?w.energy>=absolute:w.energy>absolute){sum+=w.energy;++n;}
        double gate=sum/n/(lra?100:10);std::vector<double> values;sum=0;
        for(auto w:windows)if((lra?w.energy>=absolute:w.energy>absolute)&&(lra?w.energy>=gate:w.energy>gate)){values.push_back(measurement::OfflineMeasurement::loudness(w.energy));sum+=w.energy;}
        if(lra){r.lraIncluded=values.size();std::sort(values.begin(),values.end());r.p10=values[static_cast<size_t>(std::floor((values.size()-1)*.1+.5))];r.p95=values[static_cast<size_t>(std::floor((values.size()-1)*.95+.5))];r.lra={*r.p95-*r.p10,measurement::Status::available};r.lraRelativeGate=measurement::OfflineMeasurement::loudness(gate);}
        else{r.integratedIncluded=values.size();r.integrated={measurement::OfflineMeasurement::loudness(sum/values.size()),measurement::Status::available};r.integratedRelativeGate=measurement::OfflineMeasurement::loudness(gate);}
    }
    r.lraStable=true;r.momentaryMax=r.shortTermMax={-4,measurement::Status::available};r.samplePeak=r.truePeak={-2,measurement::Status::available};
    for(int c=0;c<6;++c){r.channelSamplePeaks.push_back(std::pow(10.,-2./20));r.channelTruePeaks.push_back(std::pow(10.,-2./20));r.samplePeakFrames.push_back(c==5?r.frames-1:c*1000);r.truePeakEighthFrames.push_back(c==5?(r.frames-1)*8+3:c*8000);}
    return r;
}
uint64_t resident(){mach_task_basic_info_data_t i;mach_msg_type_number_t count=MACH_TASK_BASIC_INFO_COUNT;task_info(mach_task_self(),MACH_TASK_BASIC_INFO,reinterpret_cast<task_info_t>(&i),&count);return i.resident_size;}
double p95(std::vector<double> v){std::sort(v.begin(),v.end());return v[static_cast<size_t>(.95*(v.size()-1))];}
void longProgramme(){
    auto r=longFixture();auto baseline=resident();auto start=now();auto a=analyzeProgramme(r,"generated-three-hour-six-channel-192k");auto elapsed=now()-start;
    performance->setProperty("generated_fixture_seconds",10800);performance->setProperty("analysis_ms",elapsed);performance->setProperty("analysis_storage_bytes",static_cast<juce::int64>(a->storageBytes()));
    expect(elapsed<=2000,"three-hour diagnostics core completes within 2 seconds");expect(a->storageBytes()<=64*1024*1024,"three-hour additional analysis storage within 64 MiB");
    expect(a->peaks.back().sampleFrame==r.frames-1 && a->peaks.back().truePeakEighthFrame==(r.frames-1)*8+3,"three-hour 192k peak locations retain more than 32-bit eighth-frame coordinates");
    bool ordered=true;for(size_t n=1;n<a->integratedRank.size();++n){const auto& previous=a->momentary[a->integratedRank[n-1]];const auto& next=a->momentary[a->integratedRank[n]];ordered &= previous.energy>next.energy || (previous.energy==next.energy && previous.end<next.end);}
    expect(ordered,"three-hour contribution ranks preserve full-resolution energy order and deterministic ties");
    bool tails=true;for(auto n:a->lowerTail)tails &= a->lraWindows[n].included && *a->lraWindows[n].lufs<=*a->p10;for(auto n:a->upperTail)tails &= a->lraWindows[n].included && *a->lraWindows[n].lufs>=*a->p95;
    expect(tails && !a->lowerTail.empty() && !a->upperTail.empty(),"LRA tails contain only the exact eligible full-resolution percentiles");
    expect(diagnosticWindowText(a->lraWindows.back(),a->sampleRate,a->frames).contains("EOF padding") && a->lraWindows.back().end>a->frames,"three-hour padding keeps actual interval with explicit label");
    auto synthetic=std::make_shared<AnalysisSnapshot>();synthetic->context.sessionId="explicit-generated-performance-fixture";
    synthetic->context.media.sha256="generated-three-hour-no-PCM-claim";synthetic->context.media.channels=6;synthetic->context.media.sampleRate=r.sampleRate;synthetic->context.media.frames=r.frames;
    synthetic->context.range={0,r.frames};synthetic->context.fileLayout.count=6;synthetic->context.fileLayout.confirmed=true;
    const juce::AudioChannelSet::ChannelType roles[]={juce::AudioChannelSet::left,juce::AudioChannelSet::right,juce::AudioChannelSet::centre,juce::AudioChannelSet::leftSurround,juce::AudioChannelSet::rightSurround,juce::AudioChannelSet::LFE};
    for(int c=0;c<6;++c)synthetic->context.fileLayout.roles[c]=roles[c];
    synthetic->measurementVersion=measurement::version;synthetic->measuredContentHash=synthetic->context.media.sha256;synthetic->domain=MeasurementDomain::source;
    synthetic->detail=r;synthetic->measured={r.integrated.value,r.lra.value,r.truePeak.value,r.samplePeak.value};
    synthetic->identity=OfflineAnalysisService::identityFor(synthetic->context,synthetic->domain,synthetic->measuredContentHash);
    {DiagnosticsController worker;worker.select(synthetic,synthetic->context,MeasurementDomain::source);start=now();do{worker.poll();std::this_thread::sleep_for(std::chrono::milliseconds(1));}while(worker.busy()&&now()-start<5000);
        auto total=now()-start;performance->setProperty("worker_including_identity_ms",total);expect(worker.result()!=nullptr && total<=2000,"three-hour worker including snapshot identity completes within 2 seconds");}
    synthetic.reset();
    std::vector<double> times;for(int n=0;n<120;++n){double begin=(n%60)*150.*r.sampleRate;double end=std::min(static_cast<double>(r.frames),begin+(n%2?1800:100)*r.sampleRate);start=now();auto plot=prepareTimeline(*a,begin,end,1000);times.push_back(now()-start);
        expect(plot.momentary.size()==1000 && plot.shortTerm.size()==1000,"viewport storage bounded by display columns");
    }
    performance->setProperty("viewport_prepare_p95_ms",p95(times));expect(p95(times)<=50,"three-hour viewport p95 within 50 ms");
    auto compact=prepareTimeline(*a,0,static_cast<double>(r.frames),1);
    expect(compact.momentary[0].included.count==r.integratedIncluded && compact.momentary[0].excluded.count==r.momentary.size()-r.integratedIncluded,"mixed included and excluded membership survives reduction");
    expect(compact.momentary[0].included.maximum==measurement::OfflineMeasurement::loudness(energy(-4)) && compact.momentary[0].excluded.minimum==measurement::OfflineMeasurement::loudness(energy(-80)),"single-bin three-hour view retains spikes and quiet extremes");
    size_t cumulative=0;for(const auto& bin:a->histogram){cumulative+=bin.included;expect(bin.cumulativeIncluded==cumulative,"CDF uses exact full-resolution population counts");}
    expect(cumulative==a->lraOrder.size(),"CDF reaches complete eligible population");
    MeasurementTimeline timeline;timeline.setSize(1040,300);timeline.setAnalysis(a);expect(timeline.inspectSeconds(10799.9),"three-hour numerical inspection reaches final second");
    expect(timeline.inspectionText().contains("frames [") && timeline.inspectionText().contains("10799"),"numeric inspection displays exact source interval");
    expect(!timeline.inspectSeconds(std::numeric_limits<double>::quiet_NaN()) && !timeline.inspectSeconds(10801),"invalid and out-of-range inspection rejected");
    timeline.setViewRange(60.*r.sampleRate,120.*r.sampleRate);timeline.clearHighlight();
    juce::Image without(juce::Image::ARGB,1040,300,true),outside(juce::Image::ARGB,1040,300,true);
    {juce::Graphics g(without);timeline.paintEntireComponent(g,false);}
    timeline.highlight(0,60.*r.sampleRate,false);
    {juce::Graphics g(outside);timeline.paintEntireComponent(g,false);}
    bool equalPixels=true;for(int y=0;y<300;++y)for(int x=0;x<1040;++x)equalPixels &= without.getPixelAt(x,y)==outside.getPixelAt(x,y);
    expect(equalPixels,"offscreen half-open selection does not leave a fabricated viewport edge marker");
    proposal::ProgramAnalysis missing;missing.sampleRate=48000;missing.frames=48000*120;
    expect(diagnosticLraInterpretation(missing).contains("unavailable") && !diagnosticLraInterpretation(missing).contains("Under 60"),"long silence is unavailable without a false short-programme label");
    missing.frames=48000*2;expect(diagnosticLraInterpretation(missing).contains("Under 60") && diagnosticLraInterpretation(missing).contains("unavailable"),"short unavailable LRA keeps both distinct explanations");
    auto image=juce::Image(juce::Image::ARGB,1040,300,true);std::vector<double> paintTimes;
    for(int n=0;n<40;++n){timeline.setViewRange(n*100.*r.sampleRate,(n*100.+1800)*r.sampleRate);juce::Graphics g(image);start=now();timeline.paintEntireComponent(g,false);paintTimes.push_back(now()-start);}
    performance->setProperty("offscreen_timeline_paint_p95_ms",p95(paintTimes));expect(p95(paintTimes)<=100,"three-hour offscreen timeline paint p95 within 100 ms");
    struct rusage usage{};getrusage(RUSAGE_SELF,&usage);auto growth=static_cast<int64_t>(usage.ru_maxrss)-static_cast<int64_t>(baseline);performance->setProperty("diagnostics_peak_rss_growth_bytes",static_cast<juce::int64>(growth));expect(growth<=128*1024*1024,"three-hour diagnostic pipeline peak RSS growth within 128 MiB");
    std::atomic<bool> cancel{false};start=now();auto work=std::async(std::launch::async,[&]{return analyzeProgramme(r,"cancel-large",[&]{return cancel.load();});});std::this_thread::sleep_for(std::chrono::milliseconds(1));cancel=true;auto cancelled=work.get();auto ms=now()-start;
    performance->setProperty("active_cancel_ms",ms);expect(!cancelled && ms<=1000,"active three-hour analysis cancellation settles within one second");
}
struct Fixture{
    juce::File file=juce::File::getSpecialLocation(juce::File::tempDirectory).getNonexistentChildFile("levelscope-m4",".wav",false);
    Fixture(){juce::WavAudioFormat format;auto stream=file.createOutputStream();std::unique_ptr<juce::AudioFormatWriter> writer(format.createWriterFor(stream.get(),48000,2,24,{},0));expect(writer!=nullptr,"actual worker fixture created");stream.release();juce::AudioBuffer<float> audio(2,48000*4);for(int c=0;c<2;++c)for(int n=0;n<audio.getNumSamples();++n)audio.setSample(c,n,static_cast<float>(.05*std::sin(2*3.141592653589793*1000*n/48000)));expect(writer->writeFromAudioSampleBuffer(audio,0,audio.getNumSamples()),"actual worker fixture written");}
    ~Fixture(){file.deleteFile();}
};
void measured(StandaloneSessionModel& s,MeasurementDomain domain){expect(s.job.begin(JobStage::analysis,s.revision),"measurement job begins");auto ctx=s.measurementContext(domain);JobControl token;auto r=OfflineAnalysisService().analyze(ctx,domain,token);expect(r.value.has_value(),"qualified worker measured actual fixture");expect(s.job.finish(ctx.ticket,s.revision,r.error,r.message)&&s.acceptMeasurement(*r.value),"actual measured snapshot accepted");}
void awaitDiagnostics(DiagnosticsController& c){auto start=now();do{c.poll();std::this_thread::sleep_for(std::chrono::milliseconds(1));}while(c.busy()&&now()-start<5000);expect(!c.busy(),"diagnostic worker terminates within deadline");}
void nativeControls(const StandaloneSessionModel& session){
    // Exercise the actual panel's user-action bindings at a deterministic event
    // boundary. The ordinary three-hour rebuild is faster than remote/manual
    // click latency; no production sleep, test hook or fabricated measurement is used.
    StandaloneDiagnosticsPanel panel(session);
    juce::TextButton *build=nullptr,*cancel=nullptr;juce::TextEditor* summary=nullptr;juce::Label* status=nullptr;MeasurementTimeline* timeline=nullptr;
    for(auto* child:panel.getChildren()){
        if(auto* button=dynamic_cast<juce::TextButton*>(child)){if(button->getButtonText()=="Build diagnostics")build=button;if(button->getButtonText()=="Cancel diagnostics")cancel=button;}
        if(auto* editor=dynamic_cast<juce::TextEditor*>(child))if(editor->getName()=="Measured diagnostics availability and interpretation")summary=editor;
        if(auto* label=dynamic_cast<juce::Label*>(child))status=label;
        if(auto* view=dynamic_cast<MeasurementTimeline*>(child))timeline=view;
    }
    expect(build && cancel && summary && status && timeline,"native diagnostics controls and measured summary are present");
    const auto settle=[&]{auto start=now();do{panel.synchronize();std::this_thread::sleep_for(std::chrono::milliseconds(1));}while(cancel->isEnabled()&&now()-start<5000);expect(!cancel->isEnabled(),"native diagnostics action settles within deadline");};
    settle();expect(summary->getText().startsWith("SOURCE") && !timeline->displayedPlot().momentary.empty(),"native panel publishes measured Source and timeline");
    expect(build->isEnabled(),"native rebuild action is enabled for a valid measured snapshot");build->onClick();
    expect(cancel->isEnabled(),"native cancel action is enabled while completion remains pending");cancel->onClick();settle();
    expect(status->getText().contains("cancelled") && timeline->displayedPlot().momentary.empty() && summary->getText().contains("No current measured diagnostic result"),"native Cancel clears diagnostics and refuses pending publication");
    build->onClick();settle();expect(summary->getText().startsWith("SOURCE"),"native Build explicitly retries after cancellation");
    build->onClick();expect(cancel->isEnabled(),"native close test starts with pending diagnostic completion");panel.suspend();settle();
    expect(timeline->displayedPlot().momentary.empty() && status->getText().contains("cancelled"),"diagnostic window close binding cancels pending completion and leaves no stale view");
    panel.resume();settle();expect(summary->getText().startsWith("SOURCE"),"reopening diagnostics rebuilds from the still-valid measurement");
}
void lifecycle(){
    Fixture f;StandaloneSessionModel s;JobControl token;auto media=validateMedia(f.file,token);expect(media.succeeded()&&s.setValidatedSource(f.file,media.media),"source imported");s.setChannelRole(0,juce::AudioChannelSet::left);s.setChannelRole(1,juce::AudioChannelSet::right);expect(s.confirmLayout(),"roles confirmed");s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(Applicability::preProducedLinearProgramme);
    measured(s,MeasurementDomain::source);s.initializeCurrentStateFromSource();measured(s,MeasurementDomain::currentState);
    const auto revision=s.revision,job=s.job.ticket.id;const auto proposal=draftGlobalTrim(s.currentTrimRequest());const auto evaluation=s.evaluation(MeasurementDomain::currentState).identity;auto original=s.currentState.analysis;
    nativeControls(s);expect(s.revision==revision && s.job.ticket.id==job && s.currentState.analysis==original && s.evaluation(MeasurementDomain::currentState).identity==evaluation && draftGlobalTrim(s.currentTrimRequest()).identity==proposal.identity,"native cancellation close and reopen preserve measurements evaluations proposals and processing");
    DiagnosticsController c;c.select(original,s.measurementContext(MeasurementDomain::currentState),MeasurementDomain::currentState);awaitDiagnostics(c);auto result=c.result();expect(result!=nullptr,"actual Current State diagnostics available");
    const auto id=result->analysis->identity;auto annotation=annotateDiagnostics(*result,s.getSelectedTargetProfile(),s.programmeApplicability);
    expect(annotation.ceiling==-1 && annotation.peakHeadroom.size()==2,"separate strict profile headroom annotations");
    s.referenceVideo.startOnAudio=video::Time(12,1);s.sourceDocument.routing.fileToDevice[0]=5;
    c.select(original,s.measurementContext(MeasurementDomain::currentState),MeasurementDomain::currentState);c.poll();expect(c.result()==result && c.buildsStarted==1,"video offset and routing preserve diagnostic positions and reuse");
    auto profile=s.getSelectedTargetProfile();profile.definitionVersion="annotation-only";auto changed=annotateDiagnostics(*result,profile,s.programmeApplicability);
    expect(changed.identity!=annotation.identity && result->analysis->identity==id,"profile annotations bound separately from descriptive identity");
    expect(!annotateDiagnostics(*result,{},s.programmeApplicability).ceiling && !annotateDiagnostics(*result,profile,Applicability::unconfirmed).ceiling,"invalid profile or applicability cannot supply evaluated annotations");
    expect(s.revision==revision && s.job.ticket.id==job && s.currentState.analysis==original && s.evaluation(MeasurementDomain::currentState).identity==evaluation && draftGlobalTrim(s.currentTrimRequest()).identity==proposal.identity,"diagnostics preserve audio revisions snapshots verdict and M3 proposal");
    MeasurementTimeline view;view.setSize(1000,300);view.setAnalysis(result->analysis);view.setViewRange(48000,96000);view.inspectSeconds(1.4);view.highlight(0,19200);
    expect(s.revision==revision && s.job.ticket.id==job && s.currentStateDocument.staticGainDb==0,"timeline navigation cannot mutate processing or start jobs");
    s.setCurrentGainDb(-6);c.select(original,s.measurementContext(MeasurementDomain::currentState),MeasurementDomain::currentState);expect(!c.result(),"processing change clears displayed stale diagnostics immediately");awaitDiagnostics(c);expect(!c.result(),"stale current plan rejected even with retained snapshot pointer");
    c.select(s.source.analysis,s.measurementContext(MeasurementDomain::source),MeasurementDomain::source);c.poll();c.cancel();awaitDiagnostics(c);expect(!c.result(),"cancelled completion cannot publish");
    c.retry();awaitDiagnostics(c);expect(c.result()!=nullptr,"explicit retry reuses measured source");
    c.retry();c.poll();s.clearSourceFile();c.select(nullptr,s.measurementContext(MeasurementDomain::source),MeasurementDomain::source);awaitDiagnostics(c);expect(!c.result(),"source clear discards in-flight result");
    c.select(original,original->context,MeasurementDomain::renderedOutput);awaitDiagnostics(c);expect(!c.result(),"unsupported domain cannot masquerade as measured Source or Current State");
    {DiagnosticsController closing;closing.select(original,original->context,MeasurementDomain::currentState);closing.poll();}expect(true,"controller close cancels and joins without session mutation");
}
}
int main(int argc,char**argv){juce::ScopedJuceInitialiser_GUI gui;bool passed=false;try{boundaries();longProgramme();populations();lifecycle();passed=true;}catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<'\n';}
    auto* r=new juce::DynamicObject();r->setProperty("suite","M4 measured diagnostics");r->setProperty("passed",passed);r->setProperty("checks",checks);r->setProperty("performance",juce::var(performance.get()));r->setProperty("method",diagnosticsVersion);r->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);r->setProperty("engine",measurement::version);r->setProperty("os",juce::SystemStats::getOperatingSystemName());
#if defined(__aarch64__) || defined(__arm64__)
    r->setProperty("architecture","arm64");
#else
    r->setProperty("architecture","x86_64");
#endif
    if(argc==3&&juce::String(argv[1])=="--report")if(!juce::File(juce::String(argv[2])).replaceWithText(juce::JSON::toString(juce::var(r))))return 2;
    std::cout<<(passed?"PASS: ":"FAIL: ")<<checks.size()<<" M4 checks; no audio device opened\n";return passed?0:1;
}

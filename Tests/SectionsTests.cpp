#include "Standalone/StandaloneSectionsPanel.h"
#include "Core/Measurement/GatePopulation.h"
#include <iostream>
#include <thread>
#include <chrono>
#include <limits>
#include <numeric>
using namespace levelscope;
using namespace levelscope::standalone;
using namespace levelscope::proposal;
namespace {
juce::Array<juce::var> checks;
juce::Array<juce::var> measuredSectionValues;
juce::DynamicObject::Ptr performance = new juce::DynamicObject();
void expect(bool ok,const char* name){auto* c=new juce::DynamicObject();c->setProperty("name",name);c->setProperty("passed",ok);checks.add(juce::var(c));if(!ok)throw std::runtime_error(name);}
double now(){return juce::Time::getMillisecondCounterHiRes();}
template<class F>void rejects(F f,const char* name){bool rejected=false;try{f();}catch(const std::invalid_argument&){rejected=true;}expect(rejected,name);}
double p95(std::vector<double> v){std::sort(v.begin(),v.end());return v[static_cast<size_t>(.95*(v.size()-1))];}
double energy(double lu){return std::pow(10.,(lu+.691)/10);}
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

void boundaries(){
    for(int rate:{44100,48000,96000,192000}){
        SectionDocument d;const int64_t frames=static_cast<int64_t>(rate)*10800;juce::String error;
        expect(d.bind("content-roles-range",rate,frames)&&d.map().valid(),"whole programme initializes an exact full-coverage partition");
        const auto original=d.map().digest();expect(!d.bind("content-roles-range",rate,frames)&&d.map().digest()==original,"same source reuses section map");
        for(int64_t bad:{int64_t(-1),int64_t(0),frames,frames+1})expect(!d.split(1,bad,error)&&d.map().digest()==original,"invalid endpoint split is transactional");
        expect(d.split(1,1,error)&&d.map().sections[0].id==1&&d.map().sections[1].id==2,"one-frame split preserves left identity and creates stable right identity");
        expect(d.moveEnd(1,frames-1,error)&&d.map().sections[1].begin==frames-1,"move boundary preserves coverage and both IDs through final sample");
        const auto before=d.map().digest();expect(d.moveEnd(1,frames-1,error)&&d.map().digest()==before,"unchanged boundary is revision-free no-op");
        expect(!d.moveEnd(2,frames-1,error)&&!d.moveEnd(1,frames,error)&&d.map().digest()==before,"programme end and empty neighbour cannot move");
        expect(d.rename(2,"Final section",error)&&d.map().sections[1].id==2,"rename preserves stable section ID");
        const auto named=d.map().digest();expect(d.rename(2,"Final section",error)&&d.map().digest()==named,"unchanged name is revision-free no-op");
        expect(!d.rename(2," ",error)&&!d.rename(2,juce::String::repeatedString("a",129),error)&&d.map().digest()==named,"invalid names cannot partially mutate document");
        expect(d.mergeNext(1,error)&&d.map().sections.size()==1&&d.map().sections[0].id==1&&d.map().sections[0].end==frames,"merge retains left ID name and exact full coverage");
        expect(!d.mergeNext(1,error)&&!d.split(99,2,error),"missing neighbour and removed IDs rejected");
        expect(d.split(1,frames-1,error)&&d.map().sections[1].id==3,"removed ID is never reused");
        auto bad=d.map();bad.sections[1].begin++;expect(!bad.valid(),"gap rejected");bad=d.map();bad.sections[1].begin--;expect(!bad.valid(),"overlap rejected");
        bad=d.map();bad.sections[1].id=1;expect(!bad.valid(),"duplicate section identity rejected");bad=d.map();bad.schemaVersion=2;expect(!bad.valid(),"unsupported section schema rejected");
        bad=d.map();bad.sections[1].end=bad.sections[1].begin;expect(!bad.valid(),"empty section rejected");
        expect(d.bind("different-content",rate,frames)&&d.map().sections.size()==1,"source replacement creates a new bound partition");
        expect(d.bind({},0,0)&&!d.map().valid()&&d.map().sections.empty(),"source clear removes bound annotations");
    }
}
// Independent membership oracle scans every interval; production uses indexed ranges.
void checkSubset(const WindowSubset& subset,const std::vector<size_t>& indices){
    expect(subset.size()==indices.size(),"membership count matches independent interval oracle");
    for(size_t n=0;n<indices.size();++n)expect(subset.index(n)==indices[n],"inspection index preserves every original full-resolution window");
    expect(!subset.index(indices.size()),"inspection beyond subset is unavailable");
}
void checkSummary(const std::vector<DiagnosticWindow>& windows,const ProgrammeSection& s,int64_t frames,const SectionWindowSummary& got,bool gated,bool lra,const ProgramAnalysis& a){
    std::vector<size_t> contained,crossing,padding;std::vector<double> finite;size_t included=0,missing=0,lower=0,upper=0;double share=0;
    for(size_t n=0;n<windows.size();++n){const auto& w=windows[n];if(w.end<=s.begin||w.begin>=s.end)continue;
        if(w.end>frames){padding.push_back(n);continue;}
        if(w.begin<s.begin||w.end>s.end){crossing.push_back(n);continue;}
        contained.push_back(n);if(w.lufs)finite.push_back(*w.lufs);else ++missing;
        if(gated&&w.included){++included;if(!lra)share+=w.energyShare;}
        if(lra&&w.included&&w.lufs){lower+=a.p10&&*w.lufs<=*a.p10;upper+=a.p95&&*w.lufs>=*a.p95;}
    }
    checkSubset(got.contained,contained);checkSubset(got.crossing,crossing);checkSubset(got.padding,padding);
    expect(got.finiteCount==finite.size()&&got.missingCount==missing,"finite versus missing measurements remain distinct");
    expect(got.included==included&&got.excluded==(gated?contained.size()-included:0)&&got.lowerTail==lower&&got.upperTail==upper,"programme gate and inclusive tail membership preserved without local regating");
    expect(got.programmeBlockEnergyShare==(gated&&!lra&&!a.integratedRank.empty()?std::optional<double>(share):std::nullopt),"contained overlapping block shares use the full-resolution programme denominator or remain unavailable");
    if(finite.empty())expect(!got.minimum&&!got.median&&!got.maximum,"unavailable sections never fabricate zero values");
    else{std::sort(finite.begin(),finite.end());const auto n=finite.size();const auto median=n%2?finite[n/2]:(finite[n/2-1]+finite[n/2])/2;
        expect(got.minimum==finite.front()&&got.maximum==finite.back()&&std::abs(*got.median-median)<1.e-12,"sampled min median max match independently sorted values");}
}
void populations(){
    // Explicit diagnostic fixture: ties, silence, excluded finite data, and EOF drain.
    ProgramAnalysis a;a.identity="canonical-diagnostic-fixture";a.snapshotIdentity="measured-fixture";a.sampleRate=10;a.frames=80;a.p10=-40;a.p95=-10;
    // Use a qualified source time base while retaining readable ten-frame cadence.
    a.sampleRate=48000;const int64_t scale=4800;a.frames*=scale;
    for(int n=4;n<=80;++n){DiagnosticWindow w;w.begin=(n-4)*scale;w.end=n*scale;w.energy=n%7?energy(n%2?-10:-40):0;
        if(w.energy)w.lufs=measurement::OfflineMeasurement::loudness(w.energy);w.included=n%3!=0&&w.lufs.has_value();w.energyShare=w.included?.01:0;a.momentary.push_back(w);}
    for(int n=30;n<=95;++n){DiagnosticWindow w;w.begin=(n-30)*scale;w.end=n*scale;w.padding=n>80;w.energy=n%7?energy(n%2?-10:-40):0;if(w.energy)w.lufs=measurement::OfflineMeasurement::loudness(w.energy);w.included=n%3!=0&&w.lufs.has_value();a.lraWindows.push_back(w);if(!w.padding)a.shortTerm.push_back(w);}
    for(size_t n=0;n<a.momentary.size();++n)if(a.momentary[n].included)a.integratedRank.push_back(n);
    a.peaks.push_back({0,scale*40,scale*40*8,0.,0.});a.peaks.push_back({1,0,-1,-3.,-2.});
    SectionDocument d;d.bind("input",a.sampleRate,a.frames);juce::String error;
    for(auto f:{int64_t(1),scale*4,scale*40,a.frames-1}){const auto id=d.map().sections.back().id;expect(d.split(id,f,error),"irregular and one-frame boundary fixture created");}
    auto result=analyzeSections(a,d.map(),"input");expect(result&&result->sections.size()==d.map().sections.size(),"all sections summarized");
    for(const auto& s:result->sections){checkSummary(a.momentary,s.section,a.frames,s.momentary,true,false,a);checkSummary(a.shortTerm,s.section,a.frames,s.shortTerm,false,false,a);checkSummary(a.lraWindows,s.section,a.frames,s.lraWindows,true,true,a);}
    expect(result->sections[2].programmeMaxima.empty()&&result->sections[3].programmeMaxima[0].sample&&result->sections[3].programmeMaxima[0].truePeak,"peak at a shared boundary belongs only to following half-open section");
    expect(result->sections[0].programmeMaxima.size()==1&&!result->sections[0].programmeMaxima[0].truePeak,"out-of-source interpolated true peak is not assigned to a section");
    expect(result->identity==analyzeSections(a,d.map(),"input")->identity,"same measured analysis and map produce deterministic identity");
    auto changed=a;changed.identity="other-measured-current";expect(result->identity!=analyzeSections(changed,d.map(),"input")->identity,"measured snapshot identity binds independently");
    rejects([&]{analyzeSections(a,d.map(),"wrong-content-roles-range");},"mismatched section input rejected");
    auto bad=d.map();bad.sections.front().end++;rejects([&]{analyzeSections(a,bad,"input");},"malformed partition rejected before aggregation");
    changed=a;changed.frames--;rejects([&]{analyzeSections(changed,d.map(),"input");},"mismatched programme range rejected");
    expect(!analyzeSections(a,d.map(),"input",[]{return true;}),"pre-cancelled aggregation cannot publish");
    auto reduced=prepareTimeline(a,0,a.frames,1);(void)reduced;expect(analyzeSections(a,d.map(),"input")->identity==result->identity,"plot reduction cannot change detailed section analysis");
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
void measuredPopulations(){
    for(int rate:{44100,48000,96000,192000})for(int kind:{0,1,2}){
        auto r=signal(rate,1,kind==0?2.:4.,kind);auto a=analyzeProgramme(r,"actual-section-measurement");SectionDocument d;d.bind("actual-input",rate,r.frames);juce::String error;d.split(1,r.frames/2+1,error);
        auto result=analyzeSections(*a,d.map(),"actual-input");
        for(const auto& row:result->sections){checkSummary(a->momentary,row.section,a->frames,row.momentary,true,false,*a);checkSummary(a->shortTerm,row.section,a->frames,row.shortTerm,false,false,*a);checkSummary(a->lraWindows,row.section,a->frames,row.lraWindows,true,true,*a);
            auto* value=new juce::DynamicObject();value->setProperty("sample_rate",rate);value->setProperty("signal_kind",kind);value->setProperty("begin",static_cast<juce::int64>(row.section.begin));value->setProperty("end",static_cast<juce::int64>(row.section.end));
            for(const auto& item:std::array<std::pair<const char*,const SectionWindowSummary*>,3>{{{"m",&row.momentary},{"s",&row.shortTerm},{"lra",&row.lraWindows}}}){
                auto* population=new juce::DynamicObject();const auto& v=*item.second;
                population->setProperty("contained",static_cast<juce::int64>(v.contained.size()));population->setProperty("crossing",static_cast<juce::int64>(v.crossing.size()));population->setProperty("padding",static_cast<juce::int64>(v.padding.size()));population->setProperty("included",static_cast<juce::int64>(v.included));population->setProperty("missing",static_cast<juce::int64>(v.missingCount));
                population->setProperty("minimum",v.minimum?juce::var(*v.minimum):juce::var());population->setProperty("median",v.median?juce::var(*v.median):juce::var());population->setProperty("maximum",v.maximum?juce::var(*v.maximum):juce::var());population->setProperty("share",v.programmeBlockEnergyShare?juce::var(*v.programmeBlockEnergyShare):juce::var());value->setProperty(item.first,juce::var(population));
            }measuredSectionValues.add(juce::var(value));}
        if(kind==0)expect(!result->sections[0].shortTerm.median&&a->lra.status==measurement::Status::insufficientDuration,"short real programmes retain unavailable S/LRA without fabricated section metrics");
        if(kind==1)expect(!result->sections[0].momentary.median&&!result->sections[0].momentary.programmeBlockEnergyShare,"digital silence has neither zero loudness nor a fabricated zero contribution denominator");
        if(kind==2)expect(result->sections[0].momentary.median.has_value()&&!result->sections[0].momentary.programmeBlockEnergyShare&&result->sections[0].momentary.included==0,"finite below-gate windows retain descriptive values but unavailable programme contribution denominator");
    }
}
void longProgramme(){
    auto r=longFixture();auto a=analyzeProgramme(r,"three-hour-section-benchmark");SectionDocument d;d.bind("three-hour-input",r.sampleRate,r.frames);juce::String error;
    for(int n=1;n<1000;++n)if(!d.split(d.map().sections.back().id,r.frames*n/1000,error))throw std::runtime_error(error.toStdString());
    auto start=now();auto result=analyzeSections(*a,d.map(),"three-hour-input");const auto ms=now()-start;
    // Includes output, annotation map and worst-case median workspace; no PCM/window copies retained.
    const auto storage=result->storageBytes()+d.map().sections.capacity()*sizeof(ProgrammeSection)+a->momentary.size()*sizeof(double)+128000;
    performance->setProperty("generated_fixture_seconds",10800);performance->setProperty("section_count",1000);performance->setProperty("aggregation_ms",ms);performance->setProperty("additional_section_storage_bytes",static_cast<juce::int64>(storage));
    expect(ms<=2000&&result->sections.size()==1000,"three-hour 1000-section aggregation within two seconds");expect(storage<=16*1024*1024,"section output map and median workspace within 16 MiB");
    const auto& last=result->sections.back();expect(last.section.end==r.frames&&last.programmeMaxima.back().truePeak,"three-hour final peak retains exact eighth-frame identity");
    std::vector<double> uiTimes;MeasurementTimeline view;view.setSize(1040,300);view.setAnalysis(a);std::vector<int64_t> boundaries;for(const auto& s:d.map().sections)boundaries.push_back(s.end);view.setSectionBoundaries(boundaries);
    juce::Image image(juce::Image::ARGB,1040,300,true);
    for(int n=0;n<80;++n){const auto& s=result->sections[static_cast<size_t>(n*12)];start=now();expect(d.moveEnd(s.section.id,s.section.end-1,error)&&d.moveEnd(s.section.id,s.section.end,error),"long-programme boundary edits preserve neighbouring sections");view.highlight(s.section.begin,s.section.end);view.setViewRange(s.section.begin,std::min(r.frames,s.section.end+192000*120));juce::Graphics g(image);view.paintEntireComponent(g,false);uiTimes.push_back(now()-start);}
    performance->setProperty("interaction_p95_ms",p95(uiTimes));expect(p95(uiTimes)<=100,"three-hour section selection pan and paint p95 within 100 ms");
    std::atomic<bool> entered{false},cancel{false};start=now();auto future=std::async(std::launch::async,[&]{return analyzeSections(*a,d.map(),"three-hour-input",[&]{entered=true;return cancel.load();});});
    while(!entered.load())std::this_thread::yield();cancel=true;auto stopped=future.get();const auto cancelMs=now()-start;
    performance->setProperty("active_cancel_ms",cancelMs);expect(!stopped&&cancelMs<=1000,"active section aggregation cancels within one second");
}
struct Fixture{
    juce::File file=juce::File::getSpecialLocation(juce::File::tempDirectory).getNonexistentChildFile("levelscope-m5",".wav",false);
    Fixture(){juce::WavAudioFormat format;auto stream=file.createOutputStream();std::unique_ptr<juce::AudioFormatWriter> writer(format.createWriterFor(stream.get(),48000,2,24,{},0));expect(writer!=nullptr,"actual section worker fixture created");stream.release();juce::AudioBuffer<float> audio(2,48000*8);for(int c=0;c<2;++c)for(int n=0;n<audio.getNumSamples();++n)audio.setSample(c,n,static_cast<float>(.05*std::sin(2*3.141592653589793*1000*n/48000)));expect(writer->writeFromAudioSampleBuffer(audio,0,audio.getNumSamples()),"actual section worker fixture written");}
    ~Fixture(){file.deleteFile();}
};
void measuredState(StandaloneSessionModel& s,MeasurementDomain domain){expect(s.job.begin(JobStage::analysis,s.revision),"measurement begins");auto ctx=s.measurementContext(domain);JobControl token;auto r=OfflineAnalysisService().analyze(ctx,domain,token);expect(r.value.has_value(),"qualified worker measures actual section input");expect(s.job.finish(ctx.ticket,s.revision,r.error,r.message)&&s.acceptMeasurement(*r.value),"measured section input accepted");}
template<class C>void settle(C& c){auto start=now();do{c.poll();std::this_thread::sleep_for(std::chrono::milliseconds(1));}while(c.busy()&&now()-start<5000);expect(!c.busy(),"worker settles within deadline");}
template<class T>T* control(juce::Component& parent,const juce::String& name){for(auto* c:parent.getChildren())if(auto* p=dynamic_cast<T*>(c))if(p->getName()==name)return p;return nullptr;}
juce::TextButton* button(juce::Component& p,const juce::String& text){for(auto* c:p.getChildren())if(auto* b=dynamic_cast<juce::TextButton*>(c))if(b->getButtonText()==text)return b;return nullptr;}
void nativePanel(StandaloneSessionModel& session){
    StandaloneSectionsPanel p(session);p.setSize(1000,820);auto* split=button(p,"Split at frame");auto* move=button(p,"Move end");auto* merge=button(p,"Merge next");auto* rename=button(p,"Rename");auto* build=button(p,"Build sections");auto* cancel=button(p,"Cancel sections");
    auto* frame=control<juce::TextEditor>(p,"Section boundary source frame");auto* name=control<juce::TextEditor>(p,"Section name");auto* summary=control<juce::TextEditor>(p,"Measured section summary");auto* list=control<juce::ListBox>(p,"Programme sections");auto* domain=control<juce::ComboBox>(p,"Section measured state");
    expect(split&&move&&merge&&rename&&build&&cancel&&frame&&name&&summary&&list&&domain,"native section controls are present");
    for(auto* child:p.getChildren())expect(p.getLocalBounds().contains(child->getBounds())&&child->getHeight()>0,"all section controls fit the minimum High Sierra window size");
    auto wait=[&]{auto start=now();do{p.synchronize();std::this_thread::sleep_for(std::chrono::milliseconds(1));}while(cancel->isEnabled()&&now()-start<5000);expect(!cancel->isEnabled(),"native section worker settles");};wait();
    expect(summary->getText().contains("SOURCE")&&list->getListBoxModel()->getNumRows()==1,"native panel initializes whole programme measured Source");
    frame->setText("192000",false);split->onClick();wait();expect(list->getListBoxModel()->getNumRows()==2&&summary->getText().contains("[0, 192000) frames"),"native split updates exact frame interval");
    expect(session.sectionDocument.map().sections.size()==2&&session.documentDirty(),"section edits belong to persistent session and mark it unsaved");
    frame->setText("192001",false);move->onClick();wait();expect(summary->getText().contains("[0, 192001) frames"),"native boundary move retains sample precision");
    name->setText("Opening",false);rename->onClick();wait();expect(summary->getText().startsWith("Opening"),"native rename is applied explicitly");
    frame->setText("nan",false);split->onClick();expect(list->getListBoxModel()->getNumRows()==2,"invalid numeric boundary cannot mutate sections");
    domain->setSelectedId(2,juce::dontSendNotification);domain->onChange();wait();expect(summary->getText().contains("CURRENT STATE")&&list->getListBoxModel()->getNumRows()==2,"Source and Current reuse the same section map");
    merge->onClick();wait();expect(list->getListBoxModel()->getNumRows()==1&&summary->getText().startsWith("Opening"),"native merge retains left section name");
    build->onClick();expect(cancel->isEnabled(),"native build exposes pending completion");cancel->onClick();wait();expect(summary->getText().contains("No current measured section result"),"native cancel rejects stale section publication");
    build->onClick();wait();expect(summary->getText().contains("CURRENT STATE"),"native retry restores measured section result");
    build->onClick();p.suspend();wait();expect(summary->getText().contains("No current measured section result"),"closing section window cancels pending completion");p.resume();wait();expect(summary->getText().contains("CURRENT STATE")&&list->getListBoxModel()->getNumRows()==1,"reopen preserves annotation document and rebuilds valid results");
}
void lifecycle(){
    Fixture f;StandaloneSessionModel s;JobControl token;auto media=validateMedia(f.file,token);expect(media.succeeded()&&s.setValidatedSource(f.file,media.media),"source imports for section workflow");s.setChannelRole(0,juce::AudioChannelSet::left);s.setChannelRole(1,juce::AudioChannelSet::right);expect(s.confirmLayout(),"source roles confirmed");s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(Applicability::preProducedLinearProgramme);
    measuredState(s,MeasurementDomain::source);s.initializeCurrentStateFromSource();measuredState(s,MeasurementDomain::currentState);
    const auto revision=s.revision,job=s.job.ticket.id;const auto oldSource=s.source.analysis,oldCurrent=s.currentState.analysis;const auto eval=s.evaluation(MeasurementDomain::currentState).identity,proposal=draftGlobalTrim(s.currentTrimRequest()).identity;
    nativePanel(s);expect(s.revision==revision&&s.job.ticket.id==job&&s.source.analysis==oldSource&&s.currentState.analysis==oldCurrent&&s.evaluation(MeasurementDomain::currentState).identity==eval&&draftGlobalTrim(s.currentTrimRequest()).identity==proposal,"native section edits navigation cancellation and close preserve all audio measurement verdict and proposal identities");
    const auto ctx=s.measurementContext(MeasurementDomain::currentState);const auto input=sectionInputIdentity(ctx);SectionDocument d;d.bind(input,48000,media.media.frames);juce::String error;d.split(1,192000,error);
    DiagnosticsController dc;dc.select(s.currentState.analysis,ctx,MeasurementDomain::currentState);settle(dc);expect(dc.result()!=nullptr,"valid measured diagnostic input available");SectionsController c;c.select(dc.result(),d.map());settle(c);const auto result=c.result();expect(result!=nullptr,"valid measured sections available");
    s.referenceVideo.startOnAudio=video::Time(12,1);s.sourceDocument.routing.fileToDevice[0]=5;
    expect(sectionInputIdentity(s.measurementContext(MeasurementDomain::currentState))==input,"video offset and monitor routing cannot move section coordinates or identity");
    c.select(dc.result(),d.map());c.poll();expect(c.result()==result&&c.buildsStarted==1,"same section and measurement reused without rebuilding");
    s.setProgrammeApplicability(Applicability::unconfirmed);expect(sectionInputIdentity(s.measurementContext(MeasurementDomain::currentState))==input,"profile applicability cannot change source section identity");
    const auto oldMap=d.map().digest();d.rename(1,"Changed annotation",error);c.select(dc.result(),d.map());expect(!c.result(),"section edit immediately invalidates only section result");settle(c);expect(c.result()&&c.result()->identity!=result->identity&&d.map().digest()!=oldMap,"section revisions bind new derived result");
    expect(s.currentState.analysis==oldCurrent&&s.source.analysis==oldSource&&s.revision==revision,"section worker never invalidates audio snapshots or audio revision");
    auto wrong=d.map();wrong.inputIdentity="different-source";c.select(dc.result(),wrong);settle(c);expect(!c.result(),"individually valid mismatched map and measurement rejected");
    auto forged=std::make_shared<MeasuredDiagnostics>(*dc.result());forged->measuredBinding.currentIdentityMatches=false;c.select(forged,d.map());settle(c);expect(!c.result(),"stale measured binding cannot populate sections");
    s.setCurrentGainDb(-6);dc.select(oldCurrent,s.measurementContext(MeasurementDomain::currentState),MeasurementDomain::currentState);c.select(dc.result(),d.map());expect(!c.result()&&sectionInputIdentity(s.measurementContext(MeasurementDomain::currentState))==input,"processing invalidates Current summaries while retaining source section map");settle(dc);expect(!dc.result(),"stale Current measurement rejected by authoritative bridge");
    measuredState(s,MeasurementDomain::currentState);dc.select(s.currentState.analysis,s.measurementContext(MeasurementDomain::currentState),MeasurementDomain::currentState);settle(dc);c.select(dc.result(),d.map());settle(c);expect(c.result()&&c.result()->analysisIdentity!=result->analysisIdentity,"remeasurement replaces section truth after processing change");
    c.retry();c.poll();c.cancel();settle(c);expect(!c.result(),"cancelled section work cannot publish after completion");c.retry();settle(c);expect(c.result()!=nullptr,"explicit retry preserves section map and measurements");
    auto changed=ctx;std::swap(changed.fileLayout.roles[0],changed.fileLayout.roles[1]);expect(sectionInputIdentity(changed)!=input,"channel interpretation binds section input");changed=ctx;changed.range.end--;expect(sectionInputIdentity(changed).isEmpty(),"unsupported partial programme rejected");
    s.clearSourceFile();dc.select(nullptr,s.measurementContext(MeasurementDomain::source),MeasurementDomain::source);c.select(nullptr,{});settle(c);expect(!c.result(),"source clear discards section results");
}
}
int main(int argc,char**argv){juce::ScopedJuceInitialiser_GUI gui;bool passed=false;try{boundaries();populations();measuredPopulations();longProgramme();lifecycle();passed=true;}catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<'\n';}
    auto* r=new juce::DynamicObject();r->setProperty("suite","M5 programme sections");r->setProperty("passed",passed);r->setProperty("checks",checks);r->setProperty("performance",juce::var(performance.get()));r->setProperty("method",sectionsVersion);r->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);r->setProperty("engine",measurement::version);r->setProperty("os",juce::SystemStats::getOperatingSystemName());
#if defined(__aarch64__) || defined(__arm64__)
    r->setProperty("architecture","arm64");
#else
    r->setProperty("architecture","x86_64");
#endif
    r->setProperty("measured_section_values",measuredSectionValues);
    if(argc==3&&juce::String(argv[1])=="--report")if(!juce::File(juce::String(argv[2])).replaceWithText(juce::JSON::toString(juce::var(r))))return 2;
    std::cout<<(passed?"PASS: ":"FAIL: ")<<checks.size()<<" M5 checks; no audio device opened\n";return passed?0:1;
}

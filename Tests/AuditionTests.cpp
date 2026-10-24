#include "Standalone/StandaloneAudition.h"
#include "Standalone/StandaloneSessionModel.h"
#include "Standalone/StandaloneAudioPanel.h"
#include <iostream>
#include <thread>
#include <chrono>
#include <cstring>
#include <limits>
#include <juce_cryptography/juce_cryptography.h>
#if JUCE_MAC
#include <mach/mach.h>
#endif
using namespace levelscope;
using namespace levelscope::standalone;
namespace
{
juce::Array<juce::var> checks, cases;
juce::File sampleDirectory;
juce::var longEvidence;
void expect (bool ok, const char* name)
{
    auto* row = new juce::DynamicObject(); row->setProperty ("name", name); row->setProperty ("passed", ok);
    checks.add (juce::var (row)); if (! ok) throw std::runtime_error (name);
}
struct Temp
{
    juce::File directory = juce::File::getSpecialLocation (juce::File::tempDirectory).getChildFile ("levelscope-m6-test-" + juce::Uuid().toString());
    Temp() { expect (directory.createDirectory().wasOk(), "create owned fixture directory"); }
    ~Temp() { directory.deleteRecursively(); }
};
struct Capture final : ProcessedAudioObserver
{
    juce::AudioBuffer<float> data;
    int64_t begin, next;
    bool valid = true;
    Capture (int channels, int frames, int64_t first = 0) : data (channels, frames), begin (first), next (first) { data.clear(); }
    void observe (const float* const* input, int channels, int frames, int64_t first) noexcept override
    {
        if (first != next || channels != data.getNumChannels() || first - begin + frames > data.getNumSamples()) { valid = false; return; }
        for (int c = 0; c < channels; ++c) std::memcpy (data.getWritePointer (c, static_cast<int> (first - begin)), input[c], static_cast<size_t> (frames) * sizeof (float));
        next += frames;
    }
};
struct Device final : juce::AudioIODevice
{
    Device() : AudioIODevice ("M6 injected outputs", "Injected") {}
    juce::StringArray getOutputChannelNames() override { juce::StringArray n; for (int i=0;i<maxChannels;++i)n.add(juce::String(i+1));return n; }
    juce::StringArray getInputChannelNames() override { return {}; }
    juce::Array<double> getAvailableSampleRates() override { return { 44100,48000,96000,192000 }; }
    juce::Array<int> getAvailableBufferSizes() override { return { 256,1024 }; }
    int getDefaultBufferSize() override { return 256; }
    juce::String open (const juce::BigInteger& input, const juce::BigInteger& output, double r, int) override
    { expect (input.isZero(), "no input device opened"); ++opens; outputs=output;rate=r;opened=!rejectOpen;return rejectOpen?"injected open failure":""; }
    void close() override { stop();opened=false; }
    bool isOpen() override { return opened; }
    void start (juce::AudioIODeviceCallback* c) override { callback=c;c->audioDeviceAboutToStart(this); }
    void stop() override { if(callback)callback->audioDeviceStopped();callback=nullptr; }
    bool isPlaying() override { return callback!=nullptr; }
    juce::String getLastError() override { return {}; }
    int getCurrentBufferSizeSamples() override { return 256; }
    double getCurrentSampleRate() override { return rate; }
    int getCurrentBitDepth() override { return 32; }
    juce::BigInteger getActiveOutputChannels() const override { return outputs; }
    juce::BigInteger getActiveInputChannels() const override { return {}; }
    int getOutputLatencyInSamples() override { return 0; }
    int getInputLatencyInSamples() override { return 0; }
    void tick (juce::AudioBuffer<float>& out)
    { if(callback)callback->audioDeviceIOCallbackWithContext(nullptr,0,out.getArrayOfWritePointers(),out.getNumChannels(),out.getNumSamples(),{}); }
    bool opened=false,rejectOpen=false;int opens=0;double rate=48000;
    juce::BigInteger outputs;juce::AudioIODeviceCallback* callback=nullptr;
};
void waitPreparation (StandaloneAudition& a)
{
    const auto deadline=juce::Time::getMillisecondCounterHiRes()+10000;
    while(a.preparing()&&juce::Time::getMillisecondCounterHiRes()<deadline){a.poll();std::this_thread::sleep_for(std::chrono::milliseconds(1));}
    expect(!a.preparing(),"preparation/cancellation completes");
}
void awaitPlaying (StandaloneAudition& a, Device& d)
{
    const auto deadline=juce::Time::getMillisecondCounterHiRes()+10000;
    while(a.active()&&!d.isPlaying()&&juce::Time::getMillisecondCounterHiRes()<deadline){a.poll();std::this_thread::sleep_for(std::chrono::milliseconds(1));}
    if(!d.isPlaying())std::cerr<<a.status()<<'\n';
    expect(d.isPlaying(),"verified preparation and prefill start injected callback");
}
OutputRouting routes (int n) { OutputRouting r;for(int i=0;i<n;++i)r.fileToDevice[i]=maxChannels-1-i;return r; }
ChannelLayout layoutFor (int n)
{
    ChannelLayout l;l.count=n;l.confirmed=true;
    auto roles=n==1?juce::AudioChannelSet::mono():n==2?juce::AudioChannelSet::stereo():n==5?juce::AudioChannelSet::create5point0():n==6?juce::AudioChannelSet::create5point1():n==8?juce::AudioChannelSet::create7point1():juce::AudioChannelSet::create7point1point4();
    for(int c=0;c<n;++c)l.roles[c]=roles.getTypeOfChannel(c);
    return l;
}
AuditionRequest fixture (const juce::File& file, int rate, int channels, int frames, int signal=0, int bits=24)
{
    juce::WavAudioFormat wav;juce::AiffAudioFormat aiff;juce::FlacAudioFormat flac;
    juce::AudioFormat* format=file.hasFileExtension("flac")?static_cast<juce::AudioFormat*>(&flac):file.hasFileExtension("aiff")?static_cast<juce::AudioFormat*>(&aiff):static_cast<juce::AudioFormat*>(&wav);
    auto stream=file.createOutputStream();expect(stream!=nullptr,"fixture stream");
    std::unique_ptr<juce::AudioFormatWriter> writer(format->createWriterFor(stream.get(),rate,channels,bits,{},0));
    expect(writer!=nullptr,"fixture writer");stream.release();juce::AudioBuffer<float> block(channels,1024);
    for(int first=0;first<frames;first+=1024){const auto n=std::min(1024,frames-first);
        for(int c=0;c<channels;++c)for(int i=0;i<n;++i){const int pos=first+i;
            const double x=signal==1?.9:signal==2?0:signal==3?(c==3?.5*std::sin(pos*.1):0)
                :signal==4?(pos==0?1.25:pos==1?static_cast<double>(std::numeric_limits<float>::denorm_min()):pos==2?-0.0:pos==frames-1?-1.5:0)
                :signal==5?(pos<rate*2?.0001:.01)*std::sin(pos*.09)
                :.03*(c+1)*std::sin(pos*(c+1)*.017);
            block.setSample(c,i,static_cast<float>(x));}
        expect(writer->writeFromAudioSampleBuffer(block,0,n),"fixture payload");}
    writer.reset();JobControl control;const auto validated=validateMedia(file,control);expect(validated.succeeded(),"fixture validated");
    auto r=AuditionRequest::source(file,validated.media,layoutFor(channels));r.context.sessionId="m6-test-session";r.context.buildVersion=LEVELSCOPE_MEASUREMENT_BUILD;
    r.context.platformAndFloatingPointMode="float32-decode;float64-gain;fp-contract-off";r.currentInitialized=true;r.context.ticket.sessionRevision=1;return r;
}
void drain (StandaloneAudition& a, Device& d, int blockSize)
{
    juce::AudioBuffer<float> block(maxChannels,blockSize);
    const auto deadline=juce::Time::getMillisecondCounterHiRes()+20000;
    while(a.active()&&juce::Time::getMillisecondCounterHiRes()<deadline){a.poll();d.tick(block);std::this_thread::sleep_for(std::chrono::milliseconds(1));}
    if(a.getState()!=AuditionState::finished)std::cerr<<a.status()<<'\n';
    expect(a.getState()==AuditionState::finished&&!d.isPlaying(),"complete captured playback reaches audio EOF");
}
bool identical (const Capture& a,const Capture& b)
{
    if(!a.valid||!b.valid||a.next!=b.next||a.begin!=b.begin||a.data.getNumChannels()!=b.data.getNumChannels())return false;
    for(int c=0;c<a.data.getNumChannels();++c)if(std::memcmp(a.data.getReadPointer(c),b.data.getReadPointer(c),static_cast<size_t>(a.next-a.begin)*sizeof(float)))return false;
    return true;
}
void recordSamples(const Capture& capture,int rate,double gain,const juce::String& kind,const measurement::Result& metrics)
{
    auto* row=new juce::DynamicObject();row->setProperty("rate",rate);row->setProperty("channels",capture.data.getNumChannels());
    row->setProperty("gain_db",gain);row->setProperty("frames",capture.next-capture.begin);row->setProperty("kind",kind);
    auto* values=new juce::DynamicObject();
    auto metric=[&](const char* name,const measurement::Metric& m)
    {auto* entry=new juce::DynamicObject();entry->setProperty("status",static_cast<int>(m.status));entry->setProperty("value",m.value?juce::var(*m.value):juce::var());values->setProperty(name,juce::var(entry));};
    metric("integrated",metrics.integrated);metric("lra",metrics.lra);metric("true_peak",metrics.truePeak);
    metric("sample_peak",metrics.samplePeak);metric("momentary",metrics.momentaryMax);metric("short_term",metrics.shortTermMax);
    row->setProperty("metrics",juce::var(values));row->setProperty("lra_stable",metrics.lraStable);
    if(sampleDirectory!=juce::File())
    {
        const auto file=sampleDirectory.getChildFile(juce::String(cases.size())+".f32");auto output=file.createOutputStream();
        expect(output!=nullptr,"create owned cross-platform sample capture");
        for(int c=0;c<capture.data.getNumChannels();++c)
            expect(output->write(capture.data.getReadPointer(c),static_cast<size_t>(capture.next-capture.begin)*sizeof(float)),"write complete planar float capture");
        output->flush();expect(output->getStatus().wasOk(),"sample capture flushed");output.reset();
        juce::FileInputStream stream(file);row->setProperty("sha256",juce::SHA256(stream).toHexString());row->setProperty("file",file.getFileName());
    }
    cases.add(juce::var(row));
}
void parityMatrix (const juce::File& dir)
{
    for(int rate:{44100,48000,96000,192000})for(int channels:{1,2,5,6})
    {
        const char* extension=channels==1?".aiff":channels==2?".flac":".wav";
        auto r=fixture(dir.getChildFile(juce::String(rate)+"-"+juce::String(channels)+extension),rate,channels,rate/2+17);
        auto device=std::make_unique<Device>();auto* d=device.get();auto a=std::make_unique<StandaloneAudition>(std::move(device));
        const auto map=routes(channels);a->sourceChanged(r.context.media);a->synchronize(r,map);
        expect(a->prepare(r).isEmpty()&&!d->opened,"prepare with closed device never opens hardware");waitPreparation(*a);
        const auto* copy=a->preparedCopy();expect(copy&&d->opens==0,"private source is verified before device use");
        expect(a->connect(rate).isEmpty(),"explicit injected device connection");a->setMonitorGainDb(-40);
        for(double gain:{0.,6.020599913279624,-6.020599913279624,-3.305074118944,-60.,24.})
        {
            r.domain=MeasurementDomain::currentState;r.context.processing.staticGainDb=gain;++r.context.ticket.sessionRevision;
            a->synchronize(r,map);Capture offline(channels,static_cast<int>(r.context.media.frames));JobControl control;
            const auto measured=OfflineAnalysisService(&offline).analyze(r.context,r.domain,control);
            expect(measured.value.has_value(),"ordinary offline Current measurement succeeds");
            Capture played(channels,static_cast<int>(r.context.media.frames));expect(a->setCaptureObserver(&played),"install bounded pre-monitor capture while stopped");
            expect(a->play(r,map,true).isEmpty(),"explicit Current play accepts exact revision");awaitPlaying(*a,*d);drain(*a,*d,channels%2?997:256);
            expect(identical(offline,played)&&played.next==r.context.media.frames,"actual consumed pre-monitor stream bit-identical to offline samples");
            expect(a->preparedCopy()==copy,"gain/Current revision changes reuse verified bytes");
            measurement::OfflineMeasurement meter(rate,OfflineAnalysisService::measurementRoles(r.context.fileLayout));meter.push(played.data.getArrayOfReadPointers(),played.data.getNumSamples());
            const auto truth=meter.finish();const auto& expected=measured.value->detail;
            expect(truth.integrated.value==expected.integrated.value&&truth.lra.value==expected.lra.value&&truth.truePeak.value==expected.truePeak.value
                &&truth.samplePeak.value==expected.samplePeak.value&&truth.momentaryMax.value==expected.momentaryMax.value&&truth.shortTermMax.value==expected.shortTermMax.value,
                "remeasured captured samples reproduce I LRA M S SP TP without prediction");
            recordSamples(played,rate,gain,"Current matrix",truth);
        }
        r.domain=MeasurementDomain::source;r.context.processing.staticGainDb=0;a->synchronize(r,map);
        Capture source(channels,static_cast<int>(r.context.media.frames));a->setCaptureObserver(&source);a->play(r,map,true);awaitPlaying(*a,*d);drain(*a,*d,1024);
        juce::AudioFormatManager formats;registerSourceFormats(formats);std::unique_ptr<juce::AudioFormatReader> reader(formats.createReaderFor(r.context.sourceFile));
        juce::AudioBuffer<float> decoded(channels,source.data.getNumSamples());reader->read(decoded.getArrayOfWritePointers(),channels,0,decoded.getNumSamples());
        bool exact=true;for(int c=0;c<channels;++c)exact=exact&&!std::memcmp(decoded.getReadPointer(c),source.data.getReadPointer(c),static_cast<size_t>(decoded.getNumSamples())*sizeof(float));
        expect(exact,"Source identity gain is bit-identical to decoded file");
        measurement::OfflineMeasurement sourceMeter(rate,OfflineAnalysisService::measurementRoles(r.context.fileLayout));
        sourceMeter.push(source.data.getArrayOfReadPointers(),source.data.getNumSamples());
        recordSamples(source,rate,0,"Source matrix",sourceMeter.finish());
        a->setCaptureObserver(nullptr);a->disconnect();
    }
}
void exceptionalSignals(const juce::File& dir)
{
    for(int signal:{2,3,4,5})for(double gain:{0.,6.020599913279624,-17.384714111})
    {
        const int channels=signal==3?6:2,frames=signal==4?19:192001;
        auto r=fixture(dir.getChildFile("signal-"+juce::String(signal)+"-"+juce::String(gain,9)+".wav"),48000,channels,frames,signal,32);
        r.domain=MeasurementDomain::currentState;r.context.processing.staticGainDb=gain;
        auto device=std::make_unique<Device>();auto* d=device.get();auto a=std::make_unique<StandaloneAudition>(std::move(device));
        const auto map=routes(channels);a->synchronize(r,map);a->connect(48000);a->setMonitorGainDb(-40);
        Capture expected(channels,frames),played(channels,frames);JobControl control;
        const auto measured=OfflineAnalysisService(&expected).analyze(r.context,r.domain,control);
        expect(measured.value.has_value(),"silent LFE-only short impulse and gate-crossing offline cases succeed");
        a->setCaptureObserver(&played);a->play(r,map,true);awaitPlaying(*a,*d);drain(*a,*d,1023);
        expect(identical(expected,played),"exceptional samples including over-full-scale float remain exact");
        measurement::OfflineMeasurement meter(48000,OfflineAnalysisService::measurementRoles(r.context.fileLayout));
        meter.push(played.data.getArrayOfReadPointers(),frames);const auto actual=meter.finish();const auto& truth=measured.value->detail;
        expect(actual.integrated.value==truth.integrated.value&&actual.integrated.status==truth.integrated.status
            &&actual.lra.value==truth.lra.value&&actual.lra.status==truth.lra.status&&actual.truePeak.value==truth.truePeak.value
            &&actual.samplePeak.value==truth.samplePeak.value,"captured silence and gate-crossing metrics/statuses reproduce offline truth");
        recordSamples(played,48000,gain,"signal-"+juce::String(signal),actual);
    }
    for(const auto* ext:{".wav",".aiff",".flac"})for(int bits:{16,24})
    {
        auto r=fixture(dir.getChildFile("format-"+juce::String(bits)+ext),48000,2,15001,0,bits);
        r.domain=MeasurementDomain::currentState;r.context.processing.staticGainDb=-3.305074118944;
        auto device=std::make_unique<Device>();auto* d=device.get();auto a=std::make_unique<StandaloneAudition>(std::move(device));a->synchronize(r,routes(2));a->connect(48000);
        Capture expected(2,15001),played(2,15001);JobControl c;
        expect(OfflineAnalysisService(&expected).analyze(r.context,r.domain,c).value.has_value(),"PCM format offline processing");
        a->setCaptureObserver(&played);a->play(r,routes(2));awaitPlaying(*a,*d);drain(*a,*d,127);
        expect(identical(expected,played),"WAV AIFF FLAC PCM16/24 audition sample parity");
    }
}
void bufferAndOracle()
{
    juce::AudioBuffer<float> extreme(12,1);extreme.clear();extreme.setSample(11,0,std::numeric_limits<float>::max());
    expect(!StaticGainProcessor(24).process(extreme.getArrayOfWritePointers(),12,1),"non-finite processing overflow rejected without limiting");
    for(double db:{0.,20.,-20.,-3.305074118944,-60.,24.})
    {
        juce::AudioBuffer<float> data(maxChannels,17);auto original=data;
        for(int c=0;c<maxChannels;++c)for(int i=0;i<17;++i)data.setSample(c,i,(i-8)*(c+1)*.013f);
        original=data;expect(StaticGainProcessor(db).process(data.getArrayOfWritePointers(),maxChannels,17),"generic 12-channel processing accepts finite samples");
        bool reference=true;
        for(int c=0;c<maxChannels;++c)for(int i=0;i<17;++i){const float expected=static_cast<float>(static_cast<long double>(original.getSample(c,i))*std::pow(10.L,static_cast<long double>(db)/20.L));
            const float actual=data.getSample(c,i);reference=reference&&(actual==expected||actual==std::nextafter(expected,INFINITY)||actual==std::nextafter(expected,-INFINITY));}
        expect(reference,"independent long-double reference agrees within one float ULP");
    }
    auto fifo=std::make_unique<AuditionBuffer>();const auto map=routes(12);fifo->reset(12,map);fifo->gain.store(1);
    juce::AudioBuffer<float> input(12,4096),output(12,4096);for(int c=0;c<12;++c)juce::FloatVectorOperations::fill(input.getWritePointer(c),.2f,4096);
    fifo->push(input.getArrayOfReadPointers(),4096);fifo->render(output.getArrayOfWritePointers(),12,4096);
    const auto before=fifo->playedFrames();input.setSample(11,4095,2.f);fifo->push(input.getArrayOfReadPointers(),4096);
    fifo->render(output.getArrayOfWritePointers(),12,4096);
    expect(fifo->overloaded()&&fifo->playedFrames()==before&&output.getMagnitude(0,4096)==0,"late overload silences entire callback and consumes no source frames");
    expect(fifo->faultChannel()==11&&fifo->faultFrame()==before+4095,"overload records exact original channel and source offset");
    fifo->gain.store(.1f);fifo->reset(12,map);Capture capture(12,4096);fifo->setObserver(&capture,0);
    fifo->push(input.getArrayOfReadPointers(),4096);fifo->render(output.getArrayOfWritePointers(),12,4096);
    expect(!fifo->overloaded()&&capture.data.getSample(11,4095)==2.f,"restart clears old ramp and preserves over-full-scale pre-monitor sample");
    const float previous=output.getSample(0,4094);fifo->gain.store(.05f);fifo->push(input.getArrayOfReadPointers(),4096);fifo->setObserver(nullptr,0);fifo->render(output.getArrayOfWritePointers(),12,4096);
    expect(!fifo->overloaded()&&std::abs(output.getSample(0,0)-previous)<.001f&&output.getSample(0,4094)<previous,"live monitor volume changes remain smoothly ramped");
}
struct FailedWriter final : SnapshotWriter
{
    FailedWriter(const juce::File& f,bool atFlush,JobControl* c=nullptr):stream(f.createOutputStream()),flushFailure(atFlush),cancel(c){}
    bool write(const void* p,size_t n) override {const bool wrote=stream&&stream->write(p,n);if(cancel)cancel->cancel.store(true);return wrote&&(flushFailure||cancel);}
    bool flush() override {if(stream)stream->flush();return !flushFailure;}
    std::unique_ptr<juce::FileOutputStream> stream;bool flushFailure;JobControl* cancel;
};
struct SlowWriter final : SnapshotWriter
{
    explicit SlowWriter(const juce::File& file):stream(file.createOutputStream()){}
    bool write(const void* p,size_t n) override
    {std::this_thread::sleep_for(std::chrono::milliseconds(5));return stream&&stream->write(p,n);}
    bool flush() override{if(!stream)return false;stream->flush();return stream->getStatus().wasOk();}
    std::unique_ptr<juce::FileOutputStream> stream;
};
void asynchronousPreparation(const juce::File& dir)
{
    auto r=fixture(dir.getChildFile("async.wav"),48000,2,192000);
    StandaloneSessionModel session;session.setValidatedSource(r.context.sourceFile,r.context.media);session.sourceDocument.layout=r.context.fileLayout;
    session.initializeCurrentStateFromSource();session.setCurrentGainDb(-3.305074118944);
    const auto revision=session.revision;const auto hash=session.sourceDocument.media.sha256;
    const auto copies=dir.getChildFile("async-copies");copies.createDirectory();SnapshotStorage storage;storage.directory=copies;
    storage.availableBytes=[](const juce::File&){return 0;};
    auto device=std::make_unique<Device>();auto* d=device.get();auto a=std::make_unique<StandaloneAudition>(std::move(device));
    a->synchronize(r,routes(2));a->setSnapshotStorage(storage);
    expect(a->prepare(r).isEmpty(),"disk failure is reported asynchronously without blocking UI");waitPreparation(*a);
    expect(a->getState()==AuditionState::failed&&a->status().contains("Insufficient")&&!a->active()
        &&session.revision==revision&&session.sourceDocument.media.sha256==hash&&session.currentStateDocument.staticGainDb==-3.305074118944
        &&d->opens==0,"failed preparation retains valid session and gain, stopped with device closed");
    storage.availableBytes={};storage.createWriter=[](const juce::File& f){return std::make_unique<SlowWriter>(f);};a->setSnapshotStorage(storage);
    for(int action=0;action<3;++action)
    {
        a->synchronize(r,routes(2));expect(a->prepare(r).isEmpty(),"start cancellable source preparation");
        expect(a->connect(48000).isNotEmpty()&&d->opens==0,"pending preparation cannot unexpectedly open device");
        const auto began=juce::Time::getMillisecondCounterHiRes();
        if(action==0)a->stop();
        else if(action==1){auto other=r;other.context.sessionId="replacement";a->synchronize(other,routes(2));}
        else a->sourceChanged({});
        waitPreparation(*a);
        expect(juce::Time::getMillisecondCounterHiRes()-began<=2000&&!a->preparedCopy()&&!a->active()
            &&copies.getNumberOfChildFiles(juce::File::findFiles)==0,"cancel replacement and clear discard preparation within two seconds and remove owned files");
    }
    a->synchronize(r,routes(2));a->prepare(r);a.reset();
    expect(copies.getNumberOfChildFiles(juce::File::findFiles)==0,"closing during preparation joins and cleans partial copy");
    // Stale keys are checked before creating a copy or starting any callback.
    auto rejected=std::make_unique<StandaloneAudition>();rejected->synchronize(r,routes(2));
    for(int mismatch=0;mismatch<6;++mismatch)
    {
        auto bad=r;
        if(mismatch==0)bad.context.media.sha256="different";
        if(mismatch==1)std::swap(bad.context.fileLayout.roles[0],bad.context.fileLayout.roles[1]);
        if(mismatch==2)bad.context.range.begin=1;
        if(mismatch==3)bad.domain=MeasurementDomain::renderedOutput;
        if(mismatch==4)bad.context.processing.staticGainDb=std::numeric_limits<double>::quiet_NaN();
        if(mismatch==5)bad.context.processing.graphJson="unsupported";
        expect(rejected->prepare(bad).isNotEmpty()&&!rejected->preparing(),"wrong hash roles range domain or plan rejected before preparation");
    }
    for(double invalid:{-60.0000001,24.0000001,std::numeric_limits<double>::infinity()})
    {auto bad=r;bad.domain=MeasurementDomain::currentState;bad.context.processing.staticGainDb=invalid;expect(bad.validationError().isNotEmpty(),"invalid Current gain rejected without rounding allowance");}
}
void storageFailures(const juce::File& dir)
{
    auto r=fixture(dir.getChildFile("storage.wav"),48000,2,10000);const auto copies=dir.getChildFile("copies");copies.createDirectory();SnapshotStorage storage;storage.directory=copies;
    for(int failure=0;failure<4;++failure)
    {
        JobControl c;storage.availableBytes=[failure](const juce::File&){return failure==0?0:INT64_MAX;};
        storage.createWriter=[failure,&c](const juce::File& f)->std::unique_ptr<SnapshotWriter>{return failure==1?nullptr:std::make_unique<FailedWriter>(f,failure==3);};
        const auto result=prepareVerifiedSource(r.context,c,storage);
        expect(!result.value&&result.error==(failure==0?ErrorCode::insufficientDiskSpace:ErrorCode::writeFailed),"insufficient space/create/write/flush failures are explicit");
        expect(copies.getNumberOfChildFiles(juce::File::findFiles)==0,"storage failure removes owned partial copy");
    }
    JobControl cancel;storage.availableBytes={};storage.createWriter=[&cancel](const juce::File& f){return std::make_unique<FailedWriter>(f,false,&cancel);};
    auto result=prepareVerifiedSource(r.context,cancel,storage);expect(!result.value&&result.error==ErrorCode::cancelled&&copies.getNumberOfChildFiles(juce::File::findFiles)==0,"copy cancellation removes partial snapshot");
    storage.createWriter={};juce::MemoryBlock bytes;expect(r.context.sourceFile.loadFileAsData(bytes),"mutation fixture read");
    static_cast<char*>(bytes.getData())[bytes.getSize()-3]^=0x20;r.context.sourceFile.replaceWithData(bytes.getData(),bytes.getSize());r.context.sourceFile.setLastModificationTime(r.context.media.modified);
    expect(r.context.media.matchesFileStamp(r.context.sourceFile),"same-size timestamp-restored mutation defeats stamp-only check");
    JobControl changed;result=prepareVerifiedSource(r.context,changed,storage);expect(!result.value&&result.error==ErrorCode::sourceChanged,"content hash rejects same-stamp mutation");
    expect(copies.getNumberOfChildFiles(juce::File::findFiles)==0,"hash rejection cleans owned copy");
}
void lifecycle(const juce::File& dir)
{
    auto r=fixture(dir.getChildFile("overload.wav"),48000,2,48000,1);r.domain=MeasurementDomain::currentState;r.context.processing.staticGainDb=24;
    auto device=std::make_unique<Device>();auto* d=device.get();auto a=std::make_unique<StandaloneAudition>(std::move(device));auto map=routes(2);
    a->sourceChanged(r.context.media);a->synchronize(r,map);
    expect(a->play(r,map).isNotEmpty()&&d->opens==0,"Play with closed device never opens one implicitly");
    d->rejectOpen=true;expect(a->connect(48000).isNotEmpty()&&!d->opened&&!a->active(),"device open failure remains stopped");
    d->rejectOpen=false;a->connect(48000);a->setMonitorGainDb(0);
    a->play(r,map);awaitPlaying(*a,*d);juce::AudioBuffer<float> out(12,256);d->tick(out);a->poll();
    expect(a->getState()==AuditionState::failed&&out.getMagnitude(0,256)==0&&a->positionFrames()==0,"monitor fault silences entire first callback and preserves presented position");
    const auto* copy=a->preparedCopy();a->setMonitorGainDb(-40);a->poll();expect(!a->active()&&!d->isPlaying(),"monitor adjustment after fault never auto-resumes");
    Capture recovered(2,48000);a->setCaptureObserver(&recovered);a->play(r,map);awaitPlaying(*a,*d);drain(*a,*d,256);
    expect(a->preparedCopy()==copy&&recovered.data.getMagnitude(0,48000)>1,"explicit recovery reuses copy and preserves over-full-scale processing");
    a->setCaptureObserver(nullptr);a->play(r,map,true);awaitPlaying(*a,*d);d->tick(out);
    std::this_thread::sleep_for(std::chrono::milliseconds(10));const auto presented=a->positionFrames();
    a->setMonitorGainDb(0);d->tick(out);a->poll();
    expect(presented==256&&a->positionFrames()==presented&&out.getMagnitude(0,256)==0&&a->getState()==AuditionState::failed,
        "overload after delivered audio preserves exact presented position and suppresses the entire new callback");
    a->setMonitorGainDb(-40);a->play(r,map);awaitPlaying(*a,*d);d->tick(out);a->poll();
    expect(a->getState()==AuditionState::playing&&out.getMagnitude(0,256)>0,"lower volume and explicit resume cannot reuse previous louder ramp");a->stop();
    a->setCaptureObserver(nullptr);a->seek(12345,r,map);Capture seeked(2,48000-12345,12345);a->setCaptureObserver(&seeked);a->play(r,map);awaitPlaying(*a,*d);drain(*a,*d,997);
    expect(seeked.valid&&seeked.next==48000&&a->preparedCopy()==copy,"seek/restart retains exact source coordinates and verified copy");
    auto stale=r;stale.context.processing.staticGainDb=0;++r.context.ticket.sessionRevision;a->synchronize(r,map);
    expect(a->play(stale,map).isNotEmpty()&&!a->active(),"old processing request rejected even if independently valid");
    a->setCaptureObserver(nullptr);a->play(r,map,true);awaitPlaying(*a,*d);auto changed=r;std::swap(changed.context.fileLayout.roles[0],changed.context.fileLayout.roles[1]);a->synchronize(changed,map);
    expect(!a->active()&&!d->isPlaying(),"unexpected role mutation stops queued playback");
    a->synchronize(r,map);a->play(r,map,true);awaitPlaying(*a,*d);auto remap=map;remap.fileToDevice[0]=5;a->synchronize(r,remap);
    expect(!a->active()&&a->boundIdentity()==r.identity(),"routing mutation stops output without changing processing identity");
    a->play(r,remap,true);awaitPlaying(*a,*d);d->close();a->poll();expect(a->getState()==AuditionState::failed&&!a->active(),"device loss stops without automatic reopen");
    a->connect(48000);a->play(r,remap,true);awaitPlaying(*a,*d);d->rate=44100;d->callback->audioDeviceAboutToStart(d);d->tick(out);a->poll();
    expect(a->getState()==AuditionState::failed&&out.getMagnitude(0,256)==0,"runtime rate change stops silently without resampling");
    auto wide=fixture(dir.getChildFile("twelve.wav"),48000,12,4096);a->synchronize(wide,routes(12));a->connect(48000);a->play(wide,routes(12));awaitPlaying(*a,*d);drain(*a,*d,256);
    wide.domain=MeasurementDomain::currentState;a->synchronize(wide,routes(12));expect(a->play(wide,routes(12)).isNotEmpty(),"12-channel Source remains available but unsupported Current never falls back");
    wide.context.fileLayout=layoutFor(8);wide.context.media.channels=8;expect(wide.validationError().isNotEmpty(),"future 7.1 Current remains explicitly unsupported");
}
void modelGuards(const juce::File& dir)
{
    auto r=fixture(dir.getChildFile("session.wav"),48000,2,1000);StandaloneSessionModel s;
    expect(s.setValidatedSource(r.context.sourceFile,r.context.media),"session imported");s.sourceDocument.layout=r.context.fileLayout;s.initializeCurrentStateFromSource();
    s.selectAuditionDomain(MeasurementDomain::currentState);const auto revision=s.revision;const auto key=s.auditionRequest().identity();bool busy=true;s.auditionBusy=[&]{return busy;};
    expect(!s.setCurrentGainDb(6)&&!s.selectAuditionDomain(MeasurementDomain::source)&&!s.setOutputRoute(0,1),"model rejects processing selection and routing changes while active");
    s.setChannelRole(0,juce::AudioChannelSet::centre);s.clearCurrentState();expect(s.revision==revision&&s.auditionRequest().identity()==key,"blocked roles/Current clear preserve source and processing identity");
    s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);s.setSelectedTargetProfileIndex(1);
    expect(s.auditionRequest().identity()==key,"profile and applicability stay outside audition identity");
    busy=false;expect(s.setCurrentGainDb(6)&&s.currentStateNeedsRemeasurement()&&!s.currentState.hasMeasurement,"stopped gain edits keep Current explicitly unverified");
    expect(s.auditionRequest().validationError().isEmpty(),"unmeasured valid Current can be auditioned");
    expect(s.selectAuditionDomain(MeasurementDomain::source)&&s.auditionRequest().context.processing.staticGainDb==0&&s.currentStateDocument.staticGainDb==6,"Source selection preserves applied Current gain");
}
void panelControls(const juce::File& dir)
{
    auto r=fixture(dir.getChildFile("panel.wav"),48000,2,48000);StandaloneSessionModel model;
    model.setValidatedSource(r.context.sourceFile,r.context.media);model.sourceDocument.layout=r.context.fileLayout;model.initializeCurrentStateFromSource();
    auto device=std::make_unique<Device>();auto* d=device.get();auto audio=std::make_unique<StandaloneAudition>(std::move(device));
    std::unique_ptr<StandaloneAudioPanel> panel;
    panel=std::make_unique<StandaloneAudioPanel>(model,*audio,[&]{if(panel)panel->refresh();});panel->setSize(720,400);
    juce::ComboBox* domain=nullptr;juce::Slider* monitor=nullptr;juce::Label* status=nullptr;
    for(auto* child:panel->getChildren())
    {
        if(auto* box=dynamic_cast<juce::ComboBox*>(child);box&&box->getItemText(0)=="Source")domain=box;
        if(auto* slider=dynamic_cast<juce::Slider*>(child))monitor=slider;
        if(auto* label=dynamic_cast<juce::Label*>(child);label&&label->getText().contains("Processing"))status=label;
    }
    expect(domain&&monitor&&status&&d->opens==0,"product controls exist without physical output or playback");
    domain->setSelectedId(2,juce::sendNotificationSync);
    expect(model.auditionDomain==MeasurementDomain::currentState&&status->getText().contains("UNVERIFIED"),"Current selection explicitly labels unavailable measurements UNVERIFIED");
    model.job.begin(JobStage::analysis,model.revision);
    JobControl control;const auto measured=OfflineAnalysisService().analyze(model.measurementContext(MeasurementDomain::currentState),MeasurementDomain::currentState,control);
    model.job.finish(model.job.ticket,model.revision,measured.error,measured.message);
    expect(measured.value&&model.acceptMeasurement(*measured.value),"real Current measurement accepted");panel->refresh();
    expect(status->getText().contains("Measurements current")&&status->getText().contains("does not measure"),"current measurements and playback limits are separately visible");
    const auto identity=model.currentState.analysis->identity;const auto revision=model.revision;
    SnapshotStorage unavailable;unavailable.availableBytes=[](const juce::File&){return 0;};audio->setSnapshotStorage(unavailable);
    audio->prepare(model.auditionRequest());waitPreparation(*audio);panel->refresh();
    expect(audio->getState()==AuditionState::failed&&model.currentState.hasMeasurement&&model.currentState.analysis->identity==identity
        &&model.revision==revision&&d->opens==0,"product disk failure preserves prior valid measured session while stopped");audio->setSnapshotStorage({});
    model.setOutputRoute(0,0);model.setOutputRoute(1,1);audio->synchronize(model.auditionRequest(),model.sourceDocument.routing);
    audio->connect(48000);audio->play(model.auditionRequest(),model.sourceDocument.routing);awaitPlaying(*audio,*d);panel->refresh();
    expect(!domain->isEnabled()&&monitor->isEnabled(),"playing disables state selection but leaves listening volume enabled");
    monitor->setValue(-20,juce::sendNotificationSync);
    expect(model.revision==revision&&model.currentState.analysis->identity==identity&&model.currentState.hasMeasurement
        &&audio->boundIdentity()==model.auditionRequest().identity(),"live monitor change preserves processing measurement and audition identity");
    expect(!model.setCurrentGainDb(3)&&model.setCurrentGainDb(0),"active processing edit rejected and identical gain remains a no-op");
    audio->stop();panel->refresh();model.setCurrentGainDb(-3.305074118944);panel->refresh();
    expect(domain->isEnabled()&&status->getText().contains("UNVERIFIED")&&status->getText().contains("-3.305074118944"),"stopped gain edit shows full precision and stale Current status");
    auto stale=*measured.value;expect(!model.acceptMeasurement(stale),"previous measured result cannot become current after gain edit");
    domain->setSelectedId(1,juce::sendNotificationSync);
    expect(status->getText().startsWith("Source")&&model.currentStateDocument.staticGainDb==-3.305074118944,"Source control uses identity processing and retains applied Current trim");
    r.context.sourceFile.setLastModificationTime(r.context.media.modified+juce::RelativeTime::seconds(5));
    expect(audio->prepare(model.auditionRequest()).isNotEmpty(),"changed original stamp rejects even cached playback");panel->refresh();
    expect(model.sourceContentInvalid&&!model.currentState.hasMeasurement&&model.currentStateDocument.staticGainDb==-3.305074118944,
        "audition source-integrity failure invalidates measured verdicts while retaining applied gain");
    audio->disconnect();
}
uint64_t residentBytes()
{
#if JUCE_MAC
    mach_task_basic_info_data_t info{};mach_msg_type_number_t count=MACH_TASK_BASIC_INFO_COUNT;
    if(task_info(mach_task_self(),MACH_TASK_BASIC_INFO,reinterpret_cast<task_info_t>(&info),&count)==KERN_SUCCESS)return info.resident_size;
#endif
    return 0;
}
void longProgramme(const juce::File& file)
{
    JobControl c;const auto media=validateMedia(file,c);expect(media.succeeded(),"long programme fully validated");
    expect(media.media.channels==2&&media.media.sampleRate==48000&&media.media.frames==48000LL*10800,"actual three-hour stereo 48k input");
    auto request=AuditionRequest::source(file,media.media,layoutFor(2));request.currentInitialized=true;
    auto d=std::make_unique<Device>();auto* device=d.get();auto audio=std::make_unique<StandaloneAudition>(std::move(d));
    StandaloneSessionModel model;model.setValidatedSource(file,media.media);model.sourceDocument.layout=request.context.fileLayout;
    model.initializeCurrentStateFromSource();StandaloneAudioPanel panel(model,*audio,[]{});panel.setSize(720,400);
    request=model.auditionRequest();audio->synchronize(request,model.sourceDocument.routing);
    const uint64_t baseline=residentBytes();std::atomic<uint64_t> peak{baseline};std::atomic<bool> finished{false};
    std::thread sample([&]{while(!finished.load()){const auto rss=residentBytes();peak.store(std::max(peak.load(),rss));std::this_thread::sleep_for(std::chrono::milliseconds(5));}});
    struct Join {std::atomic<bool>& done;std::thread& thread;~Join(){done.store(true);thread.join();}} join{finished,sample};
    std::vector<double> refreshMs;juce::Array<juce::var> cancels;
    for(float stage:{0.02f,0.35f,0.65f})
    {
        expect(audio->prepare(request).isEmpty(),"long preparation begins with closed device");
        const auto deadline=juce::Time::getMillisecondCounterHiRes()+180000;
        while(audio->preparing()&&audio->preparationProgress()<stage&&juce::Time::getMillisecondCounterHiRes()<deadline)
        {audio->poll();std::this_thread::sleep_for(std::chrono::milliseconds(1));}
        expect(audio->preparing(),"cancellation exercised during copy hash and payload stages");
        const auto began=juce::Time::getMillisecondCounterHiRes();audio->stop();waitPreparation(*audio);
        const double elapsed=juce::Time::getMillisecondCounterHiRes()-began;cancels.add(elapsed);
        expect(elapsed<=2000&&!audio->preparedCopy(),"long preparation cancellation is bounded and cannot publish");
    }
    const auto began=juce::Time::getMillisecondCounterHiRes();audio->prepare(request);float progress=0;
    const auto deadline=began+300000;
    while(audio->preparing()&&juce::Time::getMillisecondCounterHiRes()<deadline)
    {
        const auto tick=juce::Time::getMillisecondCounterHiRes();panel.refresh();audio->poll();
        refreshMs.push_back(juce::Time::getMillisecondCounterHiRes()-tick);
        expect(audio->preparationProgress()>=progress,"source preparation progress monotonic");progress=audio->preparationProgress();
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    const double preparationMs=juce::Time::getMillisecondCounterHiRes()-began;
    expect(!audio->preparing()&&audio->preparedCopy()&&!device->opened,"long verified copy completes with device closed");
    const auto* copy=audio->preparedCopy();const auto copyPath=copy->file();
    expect(copyPath.getSize()==media.media.bytes,"temporary disk copy is exact source size");
    expect(audio->connect(48000).isEmpty(),"explicit long-test injected device only");audio->setMonitorGainDb(-40);
    juce::AudioFormatManager formats;registerSourceFormats(formats);std::unique_ptr<juce::AudioFormatReader> decoder(formats.createReaderFor(file));
    for(int64_t frame:{0LL,123456789LL,48000LL*10799})
    {
        request.domain=frame==0?MeasurementDomain::source:MeasurementDomain::currentState;request.currentInitialized=true;
        request.context.processing.staticGainDb=frame==0?0:-3.305074118944;request.context.ticket.sessionRevision++;
        audio->synchronize(request,routes(2));audio->seek(frame,request,routes(2));
        Capture captured(2,48000,frame);audio->setCaptureObserver(&captured);audio->play(request,routes(2));awaitPlaying(*audio,*device);
        juce::AudioBuffer<float> output(12,256);for(int n=0;n<187;++n){device->tick(output);std::this_thread::sleep_for(std::chrono::milliseconds(1));}
        audio->stop();audio->setCaptureObserver(nullptr);
        const int count=static_cast<int>(captured.next-frame);expect(captured.valid&&count>0,"long random seek captures valid source interval");
        juce::AudioBuffer<float> expected(2,count);decoder->read(expected.getArrayOfWritePointers(),2,frame,count);
        StaticGainProcessor(request.context.processing.staticGainDb).process(expected.getArrayOfWritePointers(),2,count);
        bool exact=true;for(int ch=0;ch<2;++ch)exact=exact&&!std::memcmp(expected.getReadPointer(ch),captured.data.getReadPointer(ch),count*sizeof(float));
        expect(exact&&audio->preparedCopy()==copy,"three-hour seeks and Source/Current gain switches reuse copy with exact samples");
    }
    std::sort(refreshMs.begin(),refreshMs.end());const double p95=refreshMs[static_cast<size_t>(.95*(refreshMs.size()-1))];
    const auto growth=peak.load()-baseline+sizeof(AuditionBuffer);expect(baseline>0&&growth<=16*1024*1024,"M6 long-run additional resident workspace within 16 MiB including FIFO");
    expect(p95<=100,"three-hour panel refresh p95 within 100 ms during preparation");
    auto* record=new juce::DynamicObject();record->setProperty("source_sha256",media.media.sha256);record->setProperty("frames",media.media.frames);
    record->setProperty("preparation_ms",preparationMs);record->setProperty("disk_bytes",media.media.bytes);record->setProperty("cancellation_ms",cancels);
    record->setProperty("refresh_p95_ms",p95);record->setProperty("additional_peak_resident_bytes",static_cast<juce::int64>(growth));longEvidence=juce::var(record);
    audio->disconnect();audio->sourceChanged({});expect(!copyPath.exists(),"source teardown removes owned three-hour copy");
}
}
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI init;bool passed=false;
    try{
        for(int i=1;i+1<argc;++i)if(juce::String(argv[i])=="--samples")
        {sampleDirectory=juce::File(argv[i+1]);expect(!sampleDirectory.exists()&&sampleDirectory.createDirectory().wasOk(),"fresh sample evidence directory");}
        juce::File longSource;for(int i=1;i+1<argc;++i)if(juce::String(argv[i])=="--long-source")longSource=juce::File(argv[i+1]);
        if(longSource!=juce::File())longProgramme(longSource);
        else{Temp temp;bufferAndOracle();storageFailures(temp.directory);asynchronousPreparation(temp.directory);modelGuards(temp.directory);panelControls(temp.directory);lifecycle(temp.directory);parityMatrix(temp.directory);exceptionalSignals(temp.directory);}
        passed=true;}
    catch(const std::exception& e){std::cerr<<"M6 FAIL: "<<e.what()<<'\n';}
    auto* report=new juce::DynamicObject();report->setProperty("suite","M6 verified audition");report->setProperty("passed",passed);report->setProperty("checks",checks);report->setProperty("sample_cases",cases);report->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);
    report->setProperty("long_programme",longEvidence);
#if defined(__aarch64__) || defined(__arm64__)
    report->setProperty("architecture","arm64");
#else
    report->setProperty("architecture","x86_64");
#endif
    for(int i=1;i+1<argc;++i)if(juce::String(argv[i])=="--report")juce::File(argv[i+1]).replaceWithText(juce::JSON::toString(juce::var(report),true));
    std::cout<<"M6 audition: "<<checks.size()<<" checks, "<<(passed?"PASS":"FAIL")<<'\n';return passed?0:1;
}

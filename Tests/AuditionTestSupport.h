#pragma once
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
    juce::StringArray getOutputChannelNames() override { juce::StringArray n; for (int i=0;i<outputCount;++i)n.add(juce::String(i+1));return n; }
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
    bool opened=false,rejectOpen=false;int opens=0,outputCount=maxChannels;double rate=48000;
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
}

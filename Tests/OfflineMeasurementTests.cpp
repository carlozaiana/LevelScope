#include "Standalone/StandaloneMeasurementService.h"
#include "Standalone/StandaloneSessionModel.h"
#include <juce_cryptography/juce_cryptography.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <thread>
#include <complex>
#include <future>
#include <limits>

using namespace levelscope;
using namespace levelscope::standalone;
using R = measurement::Role;
namespace
{
int checks = 0;
constexpr double pi = 3.14159265358979323846;
void expect (bool ok, const char* message) { ++checks; if (! ok) throw std::runtime_error (message); }
void near (std::optional<double> value, double target, double tolerance, const char* label)
{
    if (! value || std::abs (*value - target) > tolerance)
    {
        std::cerr << label << ": " << (value ? std::to_string (*value) : "missing") << " expected " << target << " +/- " << tolerance << '\n';
        expect (false, label);
    }
    ++checks;
}
template <typename F> measurement::Result signal (int rate, int64_t frames, const std::vector<R>& roles, F sample, int blockSize = 1024)
{
    measurement::OfflineMeasurement engine (rate, roles);
    juce::AudioBuffer<float> buffer (static_cast<int> (roles.size()), blockSize);
    for (int64_t position = 0; position < frames;)
    {
        const auto n = static_cast<int> (std::min<int64_t> (blockSize, frames-position));
        for (int c = 0; c < buffer.getNumChannels(); ++c)
            for (int i = 0; i < n; ++i) buffer.setSample (c, i, static_cast<float> (sample (c, position+i)));
        engine.push (buffer.getArrayOfReadPointers(), n); position += n;
    }
    return engine.finish();
}
void numericChecks()
{
    // Independent response oracle: directly evaluate the published 48 kHz
    // transfer polynomials, not the implementation's other-rate design equations.
    const auto referenceGain = [] (double frequency)
    {
        const auto z = std::exp (std::complex<double> (0, -2*pi*frequency/48000));
        return std::abs ((1.53512485958697-2.69169618940638*z+1.19839281085285*z*z)
            /(1.0-1.69065929318241*z+.73248077421585*z*z)
            *(1.0-2.0*z+z*z)/(1.0-1.99004745483398*z+.99007225036621*z*z));
    };
    for (int rate : {44100,48000,96000,192000})
    {
        for (double frequency : {25.,100.,500.,1000.,2000.,10000.})
        {
            const double amplitude = std::pow(10.,(-23+.691)/20)/referenceGain(frequency);
            auto calibration = signal(rate,rate*4,{R::left,R::right},[=](int,int64_t n){return amplitude*std::sin(2*pi*frequency*n/rate);});
            near(calibration.integrated.value,-23,.1,"multirate published-filter response");
        }
        const auto tone = [=] (int, int64_t n) { return std::pow (10.0,-23.0/20) * std::sin (2*pi*1000*n/rate); };
        auto result = signal (rate, rate*4, {R::left,R::right}, tone);
        near (result.integrated.value, -23, .1, "stereo 1 kHz calibration I");
        near (result.momentaryMax.value, -23, .1, "stereo calibration M");
        near (result.shortTermMax.value, -23, .1, "stereo calibration S");
        // Four-second file: the documented 1.5 s LRA drain occupies enough of
        // the population to affect P10. Analytic rectangular-window reference.
        near (result.lra.value, 10*std::log10(3.0/1.8), .01, "short-file LRA drain");
        near (result.samplePeak.value, -23, .01, "sample peak");
        near (result.truePeak.value, -23, .2, "true peak");
        expect (! result.lraStable, "short program LRA marked unstable");
        const auto chunked = signal (rate, rate*4, {R::left,R::right}, tone, 997);
        expect (result.integrated.value == chunked.integrated.value && result.lra.value == chunked.lra.value
            && result.truePeak.value == chunked.truePeak.value, "decoder chunk independence");
        for (int64_t n : {int64_t (rate*4/10-1), int64_t (rate*4/10), int64_t (rate*4/10+1), int64_t (rate*3-1), int64_t (rate*3)})
        {
            auto edge = signal (rate, n, {R::left}, tone);
            expect (edge.integrated.value.has_value() == (n >= rate*4/10), "400 ms boundary");
            expect (edge.lra.value.has_value() == (n >= rate*3), "3 s boundary");
        }
        // Published Tech 3341 TP cases 15-19, synthesized independently.
        for (auto test : { std::array<double,4>{4,0,.5,-6}, {4,45,.5,-6}, {6,60,.5,-6}, {8,67.5,.5,-6}, {4,45,1.41,3} })
        {
            auto tp = signal (rate, rate, {R::left,R::right}, [=] (int, int64_t n)
            {
                const double fade = std::min (1.0, std::min (static_cast<double> (n), static_cast<double> (rate-1-n))/(rate*.01));
                return test[2]*fade*std::sin (2*pi*n/test[0] + test[1]*pi/180);
            });
            expect (tp.truePeak.value && *tp.truePeak.value >= test[3]-.4 && *tp.truePeak.value <= test[3]+.2, "published true-peak tolerance");
        }
    }
    auto silence = signal (48000, 48000*4, {R::left}, [] (int, int64_t) {return 0;});
    expect (! silence.integrated.value && ! silence.lra.value && ! silence.truePeak.value, "silence is not numeric zero dB/LUFS");
    auto lfe = signal (48000, 48000*4, {R::left,R::right,R::centre,R::lfe,R::surroundLeft,R::surroundRight}, [] (int c, int64_t n)
    { return c == 3 ? .5*std::sin (2*pi*100*n/48000) : 0; });
    expect (! lfe.integrated.value && ! lfe.lra.value, "LFE absent from loudness");
    near (lfe.samplePeak.value,-6.0206,.0001,"LFE included in peak");
    // Tech 3341 case 6: channel gains and summing, independent signal definition.
    auto surround = signal (48000, 48000*4, {R::left,R::right,R::centre,R::surroundLeft,R::surroundRight}, [] (int c, int64_t n)
    { const double db[] {-28,-28,-24,-30,-30}; return std::pow (10.0, db[c]/20)*std::sin (2*pi*1000*n/48000); });
    near (surround.integrated.value,-23,.1,"5.0 weighted summing");
    // Tech 3342 cases 1-4: independently synthesized levels and durations.
    for (const auto& levels : {std::vector<double>{-20,-30}, {-20,-15}, {-40,-20}, {-50,-35,-20,-35,-50}})
    {
        auto lra = signal (48000, 48000*20*static_cast<int64_t> (levels.size()), {R::left,R::right}, [&] (int, int64_t n)
        { return std::pow (10.0,levels[static_cast<size_t> (n/(48000*20))]/20)*std::sin (2*pi*1000*n/48000); });
        const double expected = levels.size() == 5 ? 15 : std::abs (levels[0]-levels[1]);
        near (lra.lra.value,expected,1,"published synthetic LRA case");
    }
    bool threw = false;
    try { measurement::OfflineMeasurement bad (32000,{R::centre}); } catch (...) { threw = true; }
    expect (threw,"unqualified rate rejected");
    threw=false;
    try { signal(48000,1,{R::centre},[](int,int64_t){return std::numeric_limits<float>::quiet_NaN();}); } catch (...) {threw=true;}
    expect(threw,"nonfinite samples rejected");
    for (double level : {-70.1,-69.9})
    {
        auto gate=signal(48000,48000*4,{R::left,R::right},[=](int,int64_t n){return std::pow(10.,level/20)*std::sin(2*pi*1000*n/48000);});
        expect(gate.integrated.value.has_value()==(level > -70),"absolute gate sides");
    }
    const std::vector<R> ordered {R::left,R::right,R::centre,R::lfe,R::surroundLeft,R::surroundRight};
    const std::vector<R> reversed(ordered.rbegin(),ordered.rend());
    auto perChannel=[](int c,int64_t n){return (c+1)*.01*std::sin(2*pi*(100+c*200)*n/48000);};
    auto forward=signal(48000,48000*4,ordered,perChannel);
    auto backward=signal(48000,48000*4,reversed,[&](int c,int64_t n){return perChannel(5-c,n);});
    near(backward.integrated.value,*forward.integrated.value,1e-10,"role/sample permutation");
    for(int c=0;c<6;++c) expect(forward.channelSamplePeaks[c]==backward.channelSamplePeaks[5-c],"file-indexed peak map");
    auto startImpulse=signal(48000,48000,{R::centre},[](int,int64_t n){return n==0?.5:0.;});
    auto endImpulse=signal(48000,48000,{R::centre},[](int,int64_t n){return n==47999?.5:0.;});
    expect(startImpulse.truePeak.value==endImpulse.truePeak.value,"true-peak FIR EOF drain");

    // Tech 3341 cases 20-23: synthesize at 4*fs, independently low-pass and
    // decimate with all four offsets. No supplied EBU samples are copied.
    constexpr int highRate=192000, length=highRate;
    std::vector<double> high(length), filtered(length);
    double phase=0;
    for(int i=0;i<length;++i)
    {
        const bool event=i>=96000 && i<96016;
        const double fade=std::min(1.,std::min(double(i),double(length-1-i))/1920);
        high[i]=(event?1.:.5)*std::sin(phase)*fade;
        phase+=2*pi/(event?16:24);
    }
    std::array<double,257> fir{};double normalization=0;
    for(int k=0;k<257;++k)
    {const double x=k-128;fir[k]=(x==0?.25:std::sin(pi*x/4)/(pi*x))*(.42+.5*std::cos(pi*x/128)+.08*std::cos(2*pi*x/128));normalization+=fir[k];}
    for(int i=0;i<length;++i) for(int k=0;k<257;++k)
    {const int j=i+k-128;if(j>=0 && j<length)filtered[i]+=high[j]*fir[k]/normalization;}
    for(int offset=0;offset<4;++offset)
    {
        auto peak=signal(48000,48000,{R::left,R::right},[&](int,int64_t n){return filtered[static_cast<size_t>(n*4+offset)];});
        expect(peak.truePeak.value && *peak.truePeak.value>=-.4 && *peak.truePeak.value<=.2,"synthetic TP case 20-23");
    }
}

void serviceChecks()
{
    const auto file = juce::File::getSpecialLocation (juce::File::tempDirectory).getChildFile ("levelscope-test-"+juce::Uuid().toString()+".wav");
    struct Cleanup {juce::File f; ~Cleanup(){f.deleteFile();}} cleanup {file};
    {
        auto stream = file.createOutputStream(); juce::WavAudioFormat wav;
        std::unique_ptr<juce::AudioFormatWriter> writer (wav.createWriterFor (stream.get(),48000,2,24,{},0));
        expect (writer != nullptr,"fixture encoder"); stream.release();
        juce::AudioBuffer<float> data (2,48000*4);
        for (int c=0;c<2;++c) for(int n=0;n<data.getNumSamples();++n) data.setSample(c,n,static_cast<float>(.1*std::sin(2*pi*1000*n/48000)));
        expect (writer->writeFromAudioSampleBuffer (data,0,data.getNumSamples()),"write fixture");
    }
    JobControl control;
    const auto media = validateMedia (file,control); expect (media.succeeded(),"validate service fixture");
    ServiceContext context; context.sourceFile=file; context.media=media.media;
    context.fileLayout.count=2; context.fileLayout.roles[0]=juce::AudioChannelSet::right; context.fileLayout.roles[1]=juce::AudioChannelSet::left;
    context.fileLayout.confirmed=true; context.range={0,media.media.frames}; context.sessionId="test-session"; context.ticket={1,1};
    OfflineAnalysisService service;
    auto source=service.analyze(context,MeasurementDomain::source,control);
    expect (source.value.has_value(), source.message.toRawUTF8());
    for (int format : {1,2})
    {
        const auto alternate=file.getSiblingFile(juce::Uuid().toString()+".bin");
        Cleanup alternateCleanup{alternate};
        juce::AiffAudioFormat aiff;juce::FlacAudioFormat flac;
        juce::AudioFormat* codec=format==1?static_cast<juce::AudioFormat*>(&aiff):static_cast<juce::AudioFormat*>(&flac);
        auto stream=alternate.createOutputStream();
        std::unique_ptr<juce::AudioFormatWriter> writer(codec->createWriterFor(stream.get(),48000,2,24,{},0));
        expect(writer!=nullptr,"alternate format writer");stream.release();
        juce::AudioFormatManager formats;registerSourceFormats(formats);
        std::unique_ptr<juce::AudioFormatReader> reader(formats.createReaderFor(file));
        expect(writer->writeFromAudioReader(*reader,0,reader->lengthInSamples),"alternate format encode");writer.reset();
        auto alternateValidation=validateMedia(alternate,control);
        expect(alternateValidation.succeeded(),"alternate content validated without suffix");
        auto alternateContext=context;alternateContext.sourceFile=alternate;alternateContext.media=alternateValidation.media;
        auto measured=service.analyze(alternateContext,MeasurementDomain::source,control);
        expect(measured.value.has_value(),measured.message.toRawUTF8());
        near(measured.value->measured.integratedLufs,*source.value->measured.integratedLufs,1e-8,"AIFF/FLAC snapshot content probing");
    }
    auto neutral=service.analyze(context,MeasurementDomain::currentState,control);
    expect (neutral.value && source.value->measured.integratedLufs == neutral.value->measured.integratedLufs,"identity current equals source");
    context.processing.staticGainDb=-6;
    auto current=service.analyze(context,MeasurementDomain::currentState,control);
    expect (current.value.has_value(),current.message.toRawUTF8());
    near (current.value->measured.integratedLufs,*source.value->measured.integratedLufs-6,1e-5,"actual current gain");
    auto rendered=service.verifyTemporaryRender(context,control);
    expect (rendered.value.has_value(),rendered.message.toRawUTF8());
    near (rendered.value->measured.integratedLufs,*current.value->measured.integratedLufs,.0001,"PCM24 re-opened verification");
    expect (rendered.value->measuredContentHash != media.media.sha256,"render identity differs");
    expect (juce::SHA256(file).toHexString()==media.media.sha256,"original media preserved");

    // This checks identity binding only; no authoritative delivery profile exists.
    auto certificateContext=context;
    certificateContext.profile.id="test-only"; certificateContext.profileVersion="test-v1";
    certificateContext.profile.limits.valuesAreAuthoritative=true;
    RenderArtifact artifact; artifact.sha256=rendered.value->measuredContentHash; artifact.verifiedOutput=*rendered.value;
    artifact.verifiedOutput.context=certificateContext;
    ExportApproval approval{certificateContext.sessionId,artifact.verifiedOutput.identity,artifact.sha256,"test-only","test-v1",1};
    expect(approval.matches(certificateContext,artifact),"approval matches exact artifact identity");
    auto changed=certificateContext; changed.processing.staticGainDb=-5;
    expect(!approval.matches(changed,artifact),"plan edit rejects approval even with unchanged revision");
    changed=certificateContext; ++changed.ticket.sessionRevision;
    expect(!approval.matches(changed,artifact),"revision edit rejects approval");
    changed=certificateContext; changed.profileVersion="test-v2";
    expect(!approval.matches(changed,artifact),"profile edition rejects approval");
    changed=certificateContext;changed.profile.jobRequirementsDigest="different-job-requirements";
    expect(!approval.matches(changed,artifact),"job definition identity rejects legacy artifact approval");
    artifact.sha256="changed";
    expect(!approval.matches(certificateContext,artifact),"artifact edit rejects approval");
    expect(!service.analyze(context,MeasurementDomain::renderedOutput,control).value,"cannot relabel source as rendered output");
    changed=context; changed.range.end--;
    expect(!service.analyze(changed,MeasurementDomain::currentState,control).value,"segment cannot claim full-program metrics");
    changed=context; changed.processing.graphJson="{}";
    expect(!service.analyze(changed,MeasurementDomain::currentState,control).value,"unsupported graph rejected");
    changed=context; changed.processing.staticGainDb=24;
    expect(!service.verifyTemporaryRender(changed,control).value,"PCM clipping rejected");

    StandaloneSessionModel session;
    expect(session.setValidatedSource(file,media.media),"session fixture source");
    session.setChannelRole(0,juce::AudioChannelSet::right);session.setChannelRole(1,juce::AudioChannelSet::left);session.confirmLayout();
    session.initializeCurrentStateFromSource();
    expect(session.job.begin(JobStage::analysis,session.revision),"analysis job begins");
    const auto request=session.measurementContext(MeasurementDomain::currentState);
    auto measured=service.analyze(request,MeasurementDomain::currentState,control);
    expect(measured.value.has_value(),"session worker result");
    expect(session.job.finish(request.ticket,session.revision,measured.error,measured.message),"analysis finishes");
    expect(session.acceptMeasurement(*measured.value) && !session.currentStateNeedsRemeasurement(),"publish current measurement");
    session.sourceDocument.routing.fileToDevice[0]=11;
    expect(session.currentState.hasMeasurement && session.measurementContext(MeasurementDomain::currentState).processing.staticGainDb==0,"routing independent of current truth");
    expect(session.setCurrentGainDb(-6) && !session.currentState.hasMeasurement && session.currentStateNeedsRemeasurement(),"gain invalidates measured current");
    expect(!session.acceptMeasurement(*measured.value),"stale result cannot republish");
    expect(!session.getWorkflowReadiness().canExport(),"measurement alone cannot enable export");
    session.setChannelRole(0,juce::AudioChannelSet::left);
    expect(!session.hasCurrentStateInitialized(),"roles invalidate current plan");

    JobControl inFlight;
    auto future=std::async(std::launch::async,[&]{return service.analyze(context,MeasurementDomain::currentState,inFlight);});
    while(inFlight.progress.load()<.5f && future.wait_for(std::chrono::milliseconds(1))!=std::future_status::ready) {}
    inFlight.cancel.store(true);
    expect(future.get().error==ErrorCode::cancelled,"in-flight measurement cancellation");
    control.cancel.store(true);
    expect (service.analyze(context,MeasurementDomain::currentState,control).error==ErrorCode::cancelled,"cancelled job has no result");
    control.cancel.store(false);
    const auto stamp=file.getLastModificationTime();
    { auto stream=file.createOutputStream(); stream->setPosition(file.getSize()-4); stream->writeInt(0x12345678); }
    file.setLastModificationTime(stamp);
    expect (media.media.matchesFileStamp(file),"same-size same-time mutation fixture");
    expect (service.analyze(context,MeasurementDomain::currentState,control).error==ErrorCode::sourceChanged,"hash rejects stamp-preserving mutation");
}

juce::var metricValue (const measurement::Metric& metric)
{
    return metric.value ? juce::var (*metric.value) : juce::var();
}

juce::var referenceMetrics (const measurement::Result& result)
{
    auto* object = new juce::DynamicObject;
    object->setProperty ("I", metricValue (result.integrated));
    object->setProperty ("LRA", metricValue (result.lra));
    object->setProperty ("TP", metricValue (result.truePeak));
    object->setProperty ("SP", metricValue (result.samplePeak));
    object->setProperty ("M", metricValue (result.momentaryMax));
    object->setProperty ("S", metricValue (result.shortTermMax));
    return juce::var (object);
}

void writeUiFixtures (const juce::File& directory)
{
    if (directory.exists()) throw std::runtime_error ("Fixture directory already exists; choose a new path");
    if (! directory.createDirectory()) throw std::runtime_error ("Cannot create fixture directory");
    juce::WavAudioFormat wav; juce::AiffAudioFormat aiff; juce::FlacAudioFormat flac;
    // WAVE speaker-mask roles, all representable by JUCE's WAV writer. This is
    // deliberately outside M1's measurement layouts but valid for source import.
    auto extended = juce::AudioChannelSet::create5point1();
    for (auto role : {juce::AudioChannelSet::leftSurroundSide, juce::AudioChannelSet::rightSurroundSide,
                      juce::AudioChannelSet::topFrontLeft, juce::AudioChannelSet::topFrontRight,
                      juce::AudioChannelSet::topRearLeft, juce::AudioChannelSet::topRearRight})
        extended.addChannel (role);
    struct Fixture { const char* name; juce::AudioFormat* format; int seconds; juce::AudioChannelSet layout; };
    const std::vector<Fixture> fixtures {
        { "stereo.wav", &wav, 8, juce::AudioChannelSet::stereo() },
        { "stereo.aiff", &aiff, 8, juce::AudioChannelSet::stereo() },
        { "stereo.flac", &flac, 8, juce::AudioChannelSet::stereo() },
        { "long-stereo.wav", &wav, 600, juce::AudioChannelSet::stereo() },
        { "unsupported-12ch.wav", &wav, 4, extended }
    };
    for (const auto& fixture : fixtures)
    {
        const auto file = directory.getChildFile (fixture.name);
        auto stream = file.createOutputStream();
        std::unique_ptr<juce::AudioFormatWriter> writer (fixture.format->createWriterFor (
            stream.get(), 48000, fixture.layout, 16, {}, 0));
        if (! writer) throw std::runtime_error ("Cannot create UI fixture writer");
        stream.release();
        juce::AudioBuffer<float> samples (fixture.layout.size(), 48000);
        for (int c = 0; c < samples.getNumChannels(); ++c)
            for (int n = 0; n < samples.getNumSamples(); ++n)
                samples.setSample (c, n, static_cast<float> (std::pow (10., -23./20) * std::sin (2*pi*1000*n/48000)));
        for (int second = 0; second < fixture.seconds; ++second)
            if (! writer->writeFromAudioSampleBuffer (samples, 0, samples.getNumSamples()))
                throw std::runtime_error ("UI fixture encode failure");
        writer.reset();
        std::cout << fixture.name << ' ' << juce::SHA256 (file).toHexString() << '\n';
    }
}

void measureReference (const juce::File& file, bool curves, int64_t checkpointFrames,
                       const juce::String& explicitRoles)
{
    juce::AudioFormatManager formats; registerSourceFormats(formats);
    std::unique_ptr<juce::AudioFormatReader> reader(formats.createReaderFor(file));
    if (!reader) throw std::runtime_error("Cannot read reference");
    std::vector<R> roles = reader->numChannels==1 ? std::vector<R>{R::centre} : reader->numChannels==2 ? std::vector<R>{R::left,R::right}
        : reader->numChannels==5 ? std::vector<R>{R::left,R::right,R::centre,R::surroundLeft,R::surroundRight}
        : std::vector<R>{R::left,R::right,R::centre,R::lfe,R::surroundLeft,R::surroundRight};
    if (explicitRoles.isNotEmpty())
    {
        roles.clear();
        for (const auto& name : juce::StringArray::fromTokens (explicitRoles, ",", ""))
        {
            if (name == "L") roles.push_back (R::left);
            else if (name == "R") roles.push_back (R::right);
            else if (name == "C") roles.push_back (R::centre);
            else if (name == "Ls") roles.push_back (R::surroundLeft);
            else if (name == "Rs") roles.push_back (R::surroundRight);
            else if (name == "LFE") roles.push_back (R::lfe);
            else throw std::runtime_error ("Unknown reference channel role");
        }
    }
    if(roles.size()!=reader->numChannels) throw std::runtime_error("Reference channel layout needs explicit support");
    if (checkpointFrames < 0 || checkpointFrames > reader->lengthInSamples)
        throw std::runtime_error ("Invalid checkpoint interval");
    measurement::OfflineMeasurement engine(static_cast<int>(reader->sampleRate),roles);
    juce::AudioBuffer<float> data(static_cast<int>(roles.size()),4096);
    juce::Array<juce::var> checkpoints;
    for(juce::int64 pos=0;pos<reader->lengthInSamples;)
    { int n=static_cast<int>(std::min<juce::int64>(4096,reader->lengthInSamples-pos));
      if (checkpointFrames > 0) n = static_cast<int> (std::min<int64_t> (n, checkpointFrames - pos % checkpointFrames));
      if(!reader->read(data.getArrayOfWritePointers(),data.getNumChannels(),pos,n)) throw std::runtime_error("reference decode failure");
      engine.push(data.getArrayOfReadPointers(),n);pos+=n;
      if (checkpointFrames > 0 && pos % checkpointFrames == 0)
      {
          // Finish a copy: checkpoint observation must not drain/reset live state.
          auto snapshot = engine;
          auto checkpoint = referenceMetrics (snapshot.finish());
          checkpoint.getDynamicObject()->setProperty ("end_frame", pos);
          checkpoints.add (checkpoint);
      }
    }
    const auto& r=engine.finish();
    auto output = referenceMetrics (r);
    auto* object = output.getDynamicObject();
    object->setProperty ("sample_rate", r.sampleRate);
    object->setProperty ("frames", static_cast<juce::int64> (r.frames));
    object->setProperty ("channels", static_cast<int> (reader->numChannels));
    object->setProperty ("bits_per_sample", static_cast<int> (reader->bitsPerSample));
    object->setProperty ("roles", explicitRoles);
    object->setProperty ("engine_version", measurement::version);
    object->setProperty ("build_identity", LEVELSCOPE_MEASUREMENT_BUILD);
#if JUCE_ARM
    object->setProperty ("execution_arch", "arm64");
#else
    object->setProperty ("execution_arch", "x86_64");
#endif
    if (checkpointFrames > 0) object->setProperty ("checkpoints", checkpoints);
    if (curves)
    {
        const auto windows = [] (const std::vector<measurement::Window>& values)
        {
            juce::Array<juce::var> array;
            for (const auto& window : values)
            {
                auto* item = new juce::DynamicObject;
                item->setProperty ("end_frame", static_cast<juce::int64> (window.endFrame));
                item->setProperty ("value", window.energy > 0 ? juce::var (measurement::OfflineMeasurement::loudness (window.energy)) : juce::var());
                array.add (juce::var (item));
            }
            return array;
        };
        object->setProperty ("M_curve", windows (r.momentary));
        object->setProperty ("S_curve", windows (r.shortTerm));
    }
    std::cout << juce::JSON::toString (output, true, 15) << '\n';
}
}
int main(int argc,char**argv)
{
    const bool reporting=argc==3&&std::string(argv[1])=="--report";
    const auto report=[&](bool passed,const juce::String& failure)
    {
        if(!reporting)return;
        auto result=juce::var(new juce::DynamicObject());auto* r=result.getDynamicObject();
        r->setProperty("suite","Offline measurement and service identity checks");r->setProperty("passed",passed);r->setProperty("failure",failure);
        r->setProperty("checks",checks);r->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);r->setProperty("physical_devices_opened",false);r->setProperty("video_loaded",false);
        if(!juce::File(argv[2]).replaceWithText(juce::JSON::toString(result,false)))throw std::runtime_error("Cannot write measurement test receipt");
    };
    try
    {
        if (argc == 3 && std::string (argv[1]) == "--write-ui-fixtures")
        {
            writeUiFixtures (juce::File (argv[2])); return 0;
        }
        if (argc >= 3 && std::string (argv[1]) == "--file")
        {
            bool curves = false; int64_t checkpointFrames = 0; juce::String roles;
            for (int i = 3; i < argc; ++i)
            {
                const std::string option (argv[i]);
                if (option == "--curves") curves = true;
                else if (option == "--roles" && i+1 < argc) roles = argv[++i];
                else if (option == "--checkpoint-frames" && i+1 < argc)
                {
                    const std::string text (argv[++i]); size_t parsed = 0;
                    checkpointFrames = std::stoll (text, &parsed);
                    if (parsed != text.size() || checkpointFrames <= 0) throw std::runtime_error ("Invalid checkpoint interval");
                }
                else throw std::runtime_error ("Unknown/incomplete reference option");
            }
            measureReference (juce::File (argv[2]), curves, checkpointFrames, roles); return 0;
        }
        if (argc != 1 && !reporting) throw std::runtime_error ("Usage: --file path [--roles L,R,...] [--curves] [--checkpoint-frames N]");
        numericChecks();serviceChecks();report(true,{});std::cout<<checks<<" measurement checks passed\n";return 0;
    }
    catch(const std::exception& error){report(false,error.what());std::cerr<<"FAIL: "<<error.what()<<'\n';return 1;}
}

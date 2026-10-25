#include "AuditionTestSupport.h"
#include "Standalone/StandaloneMeasurementLayout.h"
#include "Standalone/StandaloneDiagnostics.h"
#include <numeric>
#include "Standalone/MeasurementTimeline.h"
#include "Core/Proposal/SectionAnalysis.h"
#include <sys/resource.h>
void beginVideoReaderAudit();
int videoReaderAuditCount(int);
using C = juce::AudioChannelSet;
using R = measurement::Role;
namespace
{
juce::Array<juce::var> videoCases;
void writePcm (const juce::File& file, const juce::AudioBuffer<float>& data, int rate,
               int bits = 24, bool floating = false, std::optional<uint32_t> mask = 0x63f)
{
    auto out = file.createOutputStream(); expect (out != nullptr, "create PCM fixture");
    const int channels = data.getNumChannels(), bytes = bits / 8;
    const uint32_t payload = static_cast<uint32_t> (data.getNumSamples() * channels * bytes);
    const bool extended = mask.has_value();
    out->write ("RIFF",4); out->writeInt (payload + (extended ? 60 : 36)); out->write ("WAVEfmt ",8);
    out->writeInt (extended ? 40 : 16); out->writeShort (static_cast<short> (extended ? 65534 : floating ? 3 : 1));
    out->writeShort (static_cast<short> (channels)); out->writeInt (rate); out->writeInt (rate * channels * bytes);
    out->writeShort (static_cast<short> (channels * bytes)); out->writeShort (static_cast<short> (bits));
    if (extended)
    {
        out->writeShort (22); out->writeShort (static_cast<short> (bits)); out->writeInt (static_cast<int> (*mask));
        const uint8_t guid[] { static_cast<uint8_t> (floating ? 3 : 1),0,0,0,0,0,16,0,128,0,0,170,0,56,155,113 };
        out->write (guid,16);
    }
    out->write ("data",4); out->writeInt (static_cast<int> (payload));
    for (int i=0;i<data.getNumSamples();++i) for (int ch=0;ch<channels;++ch)
    {
        if (floating) out->writeFloat (data.getSample(ch,i));
        else
        {
            const double scale=std::ldexp(1.,bits-1);
            const auto value=static_cast<uint32_t>(static_cast<int32_t>(std::llround(std::max(-scale,std::min(scale-1.,data.getSample(ch,i)*scale)))));
            for(int b=0;b<bytes;++b)out->writeByte(static_cast<char>(value>>(8*b)));
        }
    }
    out->flush();expect(out->getStatus().wasOk(),"complete PCM fixture");
}
void writeCompressedOrAiff (const juce::File& file,const juce::AudioBuffer<float>& data,int rate,int bits)
{
    // JUCE's AIFF writer offers 16/24 only, while its reader supports signed
    // PCM32. Generate that container independently instead of dropping the case.
    if(file.hasFileExtension("aiff")&&bits==32)
    {
        auto out=file.createOutputStream();expect(out!=nullptr,"create PCM32 AIFF fixture");
        const int payload=data.getNumChannels()*data.getNumSamples()*4;
        out->write("FORM",4);out->writeIntBigEndian(payload+46);out->write("AIFFCOMM",8);out->writeIntBigEndian(18);
        out->writeShortBigEndian(static_cast<short>(data.getNumChannels()));out->writeIntBigEndian(data.getNumSamples());out->writeShortBigEndian(32);
        const int exponent=static_cast<int>(std::floor(std::log2(rate)));
        out->writeShortBigEndian(static_cast<short>(16383+exponent));
        const auto mantissa=static_cast<uint64_t>(rate)<<(63-exponent);out->writeInt64BigEndian(static_cast<juce::int64>(mantissa));
        out->write("SSND",4);out->writeIntBigEndian(payload+8);out->writeIntBigEndian(0);out->writeIntBigEndian(0);
        for(int i=0;i<data.getNumSamples();++i)for(int ch=0;ch<data.getNumChannels();++ch)
            out->writeIntBigEndian(static_cast<int32_t>(std::llround(std::max(-2147483648.,std::min(2147483647.,data.getSample(ch,i)*2147483648.)))));
        out->flush();expect(out->getStatus().wasOk(),"complete independent PCM32 AIFF fixture");return;
    }
    juce::FlacAudioFormat flac;juce::AiffAudioFormat aiff;
    juce::AudioFormat& format=file.hasFileExtension("flac")?static_cast<juce::AudioFormat&>(flac):static_cast<juce::AudioFormat&>(aiff);
    auto stream=file.createOutputStream();expect(stream!=nullptr,"create format fixture");
    std::unique_ptr<juce::AudioFormatWriter> writer(format.createWriterFor(stream.get(),rate,data.getNumChannels(),bits,{},0));
    expect(writer!=nullptr,"format fixture writer");stream.release();expect(writer->writeFromAudioSampleBuffer(data,0,data.getNumSamples()),"format fixture samples");
}
juce::AudioBuffer<float> signal (int rate,int channels,int frames)
{
    juce::AudioBuffer<float> b(channels,frames);
    for(int ch=0;ch<channels;++ch)for(int i=0;i<frames;++i)
        b.setSample(ch,i,static_cast<float>(.008*(ch+1)*std::sin(2*3.14159265358979323846*(251+113*ch)*i/rate)));
    for(int ch=0;ch<channels;++ch)b.setSample(ch,frames-1-ch*17,.2f+.05f*ch);
    return b;
}
AuditionRequest validated (const juce::File& file,bool manuallyAssign=false)
{
    JobControl control;const auto m=validateMedia(file,control);
    if(!m.succeeded())std::cerr<<m.message<<'\n';expect(m.succeeded(),"M7 file validates without device");
    auto l=m.media.suggestedLayout;if(manuallyAssign){l=systemILayout();l.manuallyInterpreted=true;}
    l.confirmed=true;auto r=AuditionRequest::source(file,m.media,l);
    r.context.sessionId="m7-qualification";r.context.buildVersion=LEVELSCOPE_MEASUREMENT_BUILD;
    r.context.platformAndFloatingPointMode="float32-decode;float64-gain;fp-contract-off";
    r.context.ticket.sessionRevision=1;r.currentInitialized=true;return r;
}
measurement::Result measuredBuffer (const juce::AudioBuffer<float>& b,int rate,const ChannelLayout& l)
{
    measurement::OfflineMeasurement meter(rate,OfflineAnalysisService::measurementRoles(l));meter.push(b.getArrayOfReadPointers(),b.getNumSamples());return meter.finish();
}
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
    for(int c=0;c<8;++c){r.channelSamplePeaks.push_back(std::pow(10.,-2./20));r.channelTruePeaks.push_back(std::pow(10.,-2./20));r.samplePeakFrames.push_back(c==7?r.frames-1:c*1000);r.truePeakEighthFrames.push_back(c==7?(r.frames-1)*8+3:c*8000);}
    return r;
}

void layoutChecks()
{
    expect(proposal::ebuR128ProgrammeQc().digest()=="ef6f5cc071a8c9d5a347da0dd61fcc46e9f24da0aed787bc9f3ea48cb3db7e66","authoritative EBU profile digest exactly matches published M2");
    auto l=systemILayout();std::sort(l.roles.begin(),l.roles.begin()+8);int n=0;
    do {expect(isSystemI(l)&&measurement::OfflineMeasurement::supports(48000,OfflineAnalysisService::measurementRoles(l)),"explicit permutation preserves eight-role capability");++n;}while(std::next_permutation(l.roles.begin(),l.roles.begin()+8));
    expect(n==40320,"all 40320 ordered role permutations tested");
    for(auto invalid:{C::leftSurround,C::topFrontLeft,C::LFE2,C::leftCentre,C::unknown,C::left})
    {l=systemILayout();l.roles[4]=invalid;expect(!isSystemI(l)&&!measurement::OfflineMeasurement::supports(48000,OfflineAnalysisService::measurementRoles(l)),"ambiguous duplicate height extra-LFE wide or missing role rejected");}
    expect(systemIRoleForAlias("  Left\t Rear   Surround ")==C::leftSurroundRear&&systemIRoleForAlias("M-090")==C::rightSurroundSide,"bounded whitespace and ITU alias normalization");
    expect(systemIRoleForAlias("Ls")==C::unknown&&systemIRoleForAlias("rear-left")==C::unknown&&systemIRoleForAlias("Surround")==C::unknown,"ambiguous and fuzzy aliases rejected");
    for(const auto& label:std::initializer_list<const char*>{"L","FL","Left","Front Left","M+030"})expect(systemIRoleForAlias(label)==C::left,"left aliases");
    expect(juce::String(measurement::version)=="bs1770-5-ebu2023-m7.1","expanded measurement method explicitly versioned");
}
void weightingChecks()
{
    for(int rate:{44100,48000,96000,192000})
    {
        juce::AudioBuffer<float> one(1,rate*4),b(8,rate*4);for(int i=0;i<one.getNumSamples();++i)one.setSample(0,i,.05f*static_cast<float>(std::sin(2*3.14159265358979323846*1000*i/rate)));
        measurement::OfflineMeasurement mono(rate,{R::centre});mono.push(one.getArrayOfReadPointers(),one.getNumSamples());const double reference=*mono.finish().integrated.value;
        const double weights[]{1,1,1,0,1,1,1.41,1.41};
        for(int ch=0;ch<8;++ch)
        {
            b.clear();b.copyFrom(ch,0,one,0,0,one.getNumSamples());const auto result=measuredBuffer(b,rate,systemILayout());
            if(ch==3)expect(!result.integrated.value&&result.truePeak.value.has_value(),"LFE excluded from loudness but included in peaks");
            else expect(result.integrated.value&&std::abs(*result.integrated.value-reference-10*std::log10(weights[ch]))<1e-9,"independent per-role coefficient oracle within 1e-9 LU");
            expect(result.channelSamplePeaks[ch]>0&&result.samplePeakFrames[ch]>=0,"every channel retains peak and source location");
        }
        for(int ch=0;ch<8;++ch)b.copyFrom(ch,0,one,0,0,one.getNumSamples());
        const auto all=measuredBuffer(b,rate,systemILayout());expect(std::abs(*all.integrated.value-reference-10*std::log10(7.82))<1e-9,"independent seven-main coefficient sum 7.82");
    }
}

void peakAndPermutationChecks()
{
    for(int owner=0;owner<8;++owner)
    {
        juce::AudioBuffer<float> b(8,4800);b.clear();
        for(int i=0;i<4800;++i)b.setSample(owner,i,static_cast<float>(.8*std::sin(2*3.14159265358979323846*.25*i+3.14159265358979323846/4)));
        const auto r=measuredBuffer(b,48000,systemILayout());
        expect(r.truePeak.value&&r.samplePeak.value&&*r.truePeak.value>*r.samplePeak.value+.5,"each original channel owns measurable intersample peak including LFE and rear");
        expect(std::max_element(r.channelTruePeaks.begin(),r.channelTruePeaks.end())-r.channelTruePeaks.begin()==owner,"eight independent true-peak owner cases");
        b.clear();b.setSample(owner,0,.4f);b.setSample(owner,4799,.6f);const auto impulse=measuredBuffer(b,48000,systemILayout());
        expect(impulse.samplePeakFrames[owner]==4799&&impulse.channelTruePeaks[owner]>=.6f,"first/last boundary impulse retained on every channel");
    }
    auto b=signal(48000,8,192017);const auto l=systemILayout();const auto original=measuredBuffer(b,48000,l);
    for(const auto& order:std::vector<std::array<int,8>>{{7,6,5,4,3,2,1,0},{3,1,6,5,0,7,4,2},{0,1,2,3,6,7,4,5}})
    {
        auto permuted=l;juce::AudioBuffer<float> pcm(8,b.getNumSamples());
        for(int ch=0;ch<8;++ch){permuted.roles[ch]=l.roles[order[ch]];pcm.copyFrom(ch,0,b,order[ch],0,b.getNumSamples());}
        const auto result=measuredBuffer(pcm,48000,permuted);
        expect(std::abs(*result.integrated.value-*original.integrated.value)<=1e-9&&result.truePeak.value==original.truePeak.value,"role plus PCM permutation preserves aggregate programme values");
        bool reordered=true;for(int ch=0;ch<8;++ch)reordered&=result.channelTruePeaks[ch]==original.channelTruePeaks[order[ch]]&&result.truePeakEighthFrames[ch]==original.truePeakEighthFrames[order[ch]];
        expect(reordered,"per-channel peak values and locations follow original file order exactly");
    }
    juce::AudioBuffer<float> silence(8,19200);silence.clear();const auto empty=measuredBuffer(silence,48000,l);
    expect(!empty.integrated.value&&!empty.truePeak.value&&!empty.lra.value&&empty.channelTruePeaks.size()==8,"eight-channel silence preserves unavailable metrics and all-channel storage");
    silence.setSample(7,0,std::numeric_limits<float>::quiet_NaN());bool rejected=false;try{measuredBuffer(silence,48000,l);}catch(const std::invalid_argument&){rejected=true;}expect(rejected,"non-finite eighth channel rejected");
}

void appendFlacComments(const juce::File& file,const juce::StringArray& comments)
{
    juce::MemoryBlock original;expect(file.loadFileAsData(original),"read fixture FLAC metadata");
    auto* data=static_cast<uint8_t*>(original.getData());size_t pos=4,last=4;
    do{last=pos;const size_t length=(size_t(data[pos+1])<<16)|(size_t(data[pos+2])<<8)|data[pos+3];pos+=4+length;}while((data[last]&128)==0);
    data[last]&=127;juce::MemoryOutputStream block;block.writeInt(0);block.writeInt(comments.size());
    for(const auto& comment:comments){block.writeInt(static_cast<int>(comment.getNumBytesAsUTF8()));block.write(comment.toRawUTF8(),comment.getNumBytesAsUTF8());}
    juce::MemoryOutputStream out;out.write(data,pos);out.writeByte(static_cast<char>(132));const auto size=block.getDataSize();out.writeByte(static_cast<char>(size>>16));out.writeByte(static_cast<char>(size>>8));out.writeByte(static_cast<char>(size));out.write(block.getData(),size);out.write(data+pos,original.getSize()-pos);
    expect(file.replaceWithData(out.getData(),out.getDataSize()),"add explicit FLAC mask metadata");
}
void metadataAndLegacy(const juce::File& dir)
{
    auto b=signal(48000,8,24001);
    for(const auto mask:{0u,0x63fu,0xffu,0x3fu})
    {
        auto file=dir.getChildFile("mask-"+juce::String(mask)+".wav");writePcm(file,b,48000,24,false,mask);JobControl c;const auto m=validateMedia(file,c);
        if(mask==0x3f)expect(!m.succeeded(),"eight-channel mask with six bits rejected");
        else expect(m.succeeded()&&(isSystemI(m.media.suggestedLayout)==(mask==0x63f)),"WAVE exact mask alone supplies System I ordered suggestion");
    }
    auto plain=dir.getChildFile("unlabelled.wav");writePcm(plain,b,48000,24,false,{});JobControl c;auto m=validateMedia(plain,c);expect(m.succeeded()&&!m.media.suggestedLayout.isValid(),"plain eight-channel WAV count does not infer roles");
    for(const auto text:{"0x63f","0X0000063F","0xff","0x3f","0x0","0xQQ"})
    {
        const auto file=dir.getChildFile("override-"+juce::String(text)+".flac");writeCompressedOrAiff(file,b,48000,24);appendFlacComments(file,{"waveformatextensible_channel_mask="+juce::String(text)});JobControl job;const auto result=validateMedia(file,job);
        if(juce::String(text)=="0x3f"||juce::String(text)=="0xQQ")expect(!result.succeeded(),"malformed or channel-inconsistent FLAC override rejected");
        else expect(result.succeeded()&&(isSystemI(result.media.suggestedLayout)==juce::String(text).containsIgnoreCase("63f")),"native FLAC override respected without count fallback");
    }
    auto duplicate=dir.getChildFile("duplicate.flac");writeCompressedOrAiff(duplicate,b,48000,24);appendFlacComments(duplicate,{"WAVEFORMATEXTENSIBLE_CHANNEL_MASK=0x63f","WAVEFORMATEXTENSIBLE_CHANNEL_MASK=0x63f"});JobControl dc;expect(!validateMedia(duplicate,dc).succeeded(),"duplicate mask fails closed");
    for(const char* ext:{"wav","flac"})
    {
        auto six=signal(48000,6,192000);six.clear();for(int i=0;i<six.getNumSamples();++i)six.setSample(4,i,.05f*static_cast<float>(std::sin(i*.13)));
        const auto f=dir.getChildFile("legacy-51."+juce::String(ext));
        if(juce::String(ext)=="wav")writePcm(f,six,48000,24,false,0x3f);else writeCompressedOrAiff(f,six,48000,24);
        auto r=validated(f);expect(r.context.fileLayout.roles[4]==C::leftSurround&&r.context.fileLayout.roles[5]==C::rightSurround,"concrete conventional WAV/FLAC back-labelled 5.1 keeps legacy decoder surrounds");
        expect(OfflineAnalysisService::measurementRoles(r.context.fileLayout)[4]==R::surroundLeft,"legacy 5.1 maps to unchanged 1.41 role");
        JobControl job;auto actual=OfflineAnalysisService().analyze(r.context,r.domain,job);expect(actual.value.has_value(),"legacy 5.1 still measures");
        juce::AudioFormatManager fm;registerSourceFormats(fm);std::unique_ptr<juce::AudioFormatReader> reader(fm.createReaderFor(f));juce::AudioBuffer<float> decoded(6,192000);reader->read(decoded.getArrayOfWritePointers(),6,0,192000);
        const auto legacy=measuredBuffer(decoded,48000,layoutFor(6));expect(actual.value->detail.integrated.value==legacy.integrated.value&&actual.value->detail.truePeak.value==legacy.truePeak.value,"legacy 5.1 numerical results exactly match established M1 interpretation");
        measurement::OfflineMeasurement mono(48000,{R::centre});const float* monoData[]{decoded.getReadPointer(4)};mono.push(monoData,192000);
        expect(std::abs(*actual.value->detail.integrated.value-*mono.finish().integrated.value-10*std::log10(1.41))<1e-9,"WAV and FLAC 5.1 surround retains measured 1.41 coefficient");
    }
}

void acceptMeasured (StandaloneSessionModel& session,MeasurementDomain domain)
{
    session.job.state=JobState::succeeded;
    JobControl job;auto result=OfflineAnalysisService().analyze(session.measurementContext(domain),domain,job);
    expect(result.value.has_value()&&session.acceptMeasurement(std::move(*result.value)),"device-closed worker result accepted by session");
}
void mutedPeakWorkflow(const juce::File& dir)
{
    for(int owner:{4,3}) // file-order Lrs, then LFE
    {
        const int frames=48000*4;
        auto b=signal(48000,8,frames);
        for(int ch=0;ch<8;++ch)for(int i=0;i<frames;++i)b.setSample(ch,i,.008f*static_cast<float>(std::sin(i*(.071+ch*.013))));
        const int peakFrame=frames/2+19;b.setSample(owner,peakFrame,.7f);
        for(int iteration=0;iteration<4;++iteration)
        {
            const auto result=measuredBuffer(b,48000,systemILayout());
            const double desiredPeak=*result.integrated.value+22.1;
            b.setSample(owner,peakFrame,static_cast<float>(b.getSample(owner,peakFrame)*std::pow(10.,(desiredPeak-*result.truePeak.value)/20.)));
        }
        const auto file=dir.getChildFile("muted-peak-"+juce::String(owner)+".wav");writePcm(file,b,48000,32,true);
        auto r=validated(file);StandaloneSessionModel session;
        expect(session.setValidatedSource(file,r.context.media)&&session.confirmLayout(),"System I roles explicitly confirmed in session");
        session.initializeCurrentStateFromSource();session.setSelectedTargetProfileIndex(1);session.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
        auto injected=std::make_unique<Device>();auto* device=injected.get();device->rejectOpen=true;
        auto audio=std::make_unique<StandaloneAudition>(std::move(injected));
        acceptMeasured(session,MeasurementDomain::source);acceptMeasured(session,MeasurementDomain::currentState);
        const auto source=session.source.analysis,current=session.currentState.analysis;
        const auto& d=current->detail;expect(std::max_element(d.channelTruePeaks.begin(),d.channelTruePeaks.end())-d.channelTruePeaks.begin()==owner,"unique highest measured true peak belongs to requested rear/LFE channel");
        const auto evaluation=session.evaluation(MeasurementDomain::currentState);
        const auto proposal=proposal::draftGlobalTrim(session.currentTrimRequest());
        expect(proposal.status==proposal::TrimStatus::draft&&proposal.deltaDb&&proposal.predictedTruePeak
            &&std::abs(*proposal.predictedTruePeak+1)<1e-12&&*proposal.deltaDb<-23-*d.integrated.value,"real measured rear/LFE peak constrains trim within QC allowance");
        for(int ch=0;ch<8;++ch)session.setOutputRoute(ch,ch);session.setOutputRoute(owner,-1);
        expect(session.source.analysis==source&&session.currentState.analysis==current
            &&session.evaluation(MeasurementDomain::currentState).identity==evaluation.identity
            &&proposal::draftGlobalTrim(session.currentTrimRequest()).identity==proposal.identity,"muting highest peak changes neither measurement nor M2 verdict nor M3 proposal");
        const auto id=current->identity;session.referenceVideo.startOnAudio=video::Time(-1001,30000);
        expect(OfflineAnalysisService::identityFor(session.measurementContext(MeasurementDomain::currentState),MeasurementDomain::currentState,session.sourceDocument.media.sha256)==id,"video offset preserves eight-channel audio measurement identity");
        auto analysis=proposal::analyzeProgramme(d,id);expect(analysis&&analysis->peaks.size()==8,"diagnostics retains all eight peak locations including muted rear/LFE");
        auto stale=*current;stale.measurementVersion="bs1770-5-ebu2023-m1.1";
        expect(evaluateSnapshot(&stale,current->context,MeasurementDomain::currentState,session.getSelectedTargetProfile(),session.programmeApplicability).verdict!=proposal::Verdict::pass,"old measurement capability cannot satisfy expanded profile evaluation");
        auto wrong=current->context;std::swap(wrong.fileLayout.roles[4],wrong.fileLayout.roles[6]);
        expect(evaluateSnapshot(current.get(),wrong,MeasurementDomain::currentState,session.getSelectedTargetProfile(),session.programmeApplicability).verdict!=proposal::Verdict::pass,"ordered side/rear identity mismatch rejected");
        expect(proposal::draftGlobalTrim(trimRequest(source.get(),current.get(),source->context,wrong,session.getSelectedTargetProfile(),session.programmeApplicability)).status==proposal::TrimStatus::unavailable,"M3 mismatched ordered inputs rejected");
        expect(session.acceptTrimProposal(proposal)&&session.applyTrim()&&session.currentStateNeedsRemeasurement(),"explicit trim application invalidates Current while preserving full gain");
        session.measurementUnsuccessful(MeasurementDomain::currentState,ErrorCode::cancelled);
        expect(session.currentStateNeedsRemeasurement()&&session.currentStateDocument.staticGainDb==*proposal.totalGainDb&&!session.currentState.hasMeasurement,"cancelled verification retains applied gain and UNVERIFIED state");
        acceptMeasured(session,MeasurementDomain::currentState);
        expect(session.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"actual remeasurement evaluates corrected eight-channel result");
        const auto revision=session.revision;const auto accepted=session.currentState.analysis;
        const auto noChange=proposal::draftGlobalTrim(session.currentTrimRequest());
        expect(noChange.status==proposal::TrimStatus::noChange&&session.acceptTrimProposal(noChange)&&session.applyTrim()
            &&session.revision==revision&&session.currentState.analysis==accepted,"already passing 7.1 is a true no-op");
        expect(device->opens==0,"entire import measure evaluate propose apply verify workflow opens zero devices");
        // Physical output policy is exercised separately with an injected backend.
        device->rejectOpen=false;r=session.auditionRequest();r.domain=MeasurementDomain::currentState;r.context=session.measurementContext(r.domain);r.currentInitialized=true;
        auto map=session.sourceDocument.routing;audio->synchronize(r,map);expect(audio->connect(48000).isEmpty(),"explicit injected backend opened only for routing evidence");audio->setMonitorGainDb(-30);
        Capture offline(8,frames),played(8,frames);JobControl control;auto measured=OfflineAnalysisService(&offline).analyze(r.context,r.domain,control);expect(measured.value.has_value(),"muted fixture offline capture");
        audio->setCaptureObserver(&played);expect(audio->play(r,map,true).isEmpty(),"partial monitoring accepts explicit mute");awaitPlaying(*audio,*device);drain(*audio,*device,997);
        expect(identical(offline,played)&&played.data.getMagnitude(owner,0,frames)>0,"muted rear/LFE remains present bit-exact in every pre-monitor processing capture");
        recordSamples(played,48000,r.context.processing.staticGainDb,"muted highest peak "+juce::String(owner),measured.value->detail);audio->disconnect();
    }
}
void routingCapture()
{
    for(int owner=0;owner<8;++owner)
    {
        auto fifo=std::make_unique<AuditionBuffer>();const auto map=routes(8);fifo->gain.store(1);fifo->reset(8,map);
        juce::AudioBuffer<float> input(8,4096),out(12,4096);input.clear();juce::FloatVectorOperations::fill(input.getWritePointer(owner),2.f,4096);
        fifo->push(input.getArrayOfReadPointers(),4096);fifo->render(out.getArrayOfWritePointers(),12,4096);
        expect(fifo->overloaded()&&fifo->faultChannel()==owner&&fifo->playedFrames()==0&&out.getMagnitude(0,4096)==0,"each mapped 7.1 channel overload silences entire callback without advancing presentation");
        fifo->gain.store(.1f);fifo->reset(8,map);Capture capture(8,4096);fifo->setObserver(&capture,0);fifo->push(input.getArrayOfReadPointers(),4096);fifo->render(out.getArrayOfWritePointers(),12,4096);
        expect(!fifo->overloaded()&&capture.data.getSample(owner,4095)==2.f&&out.getMagnitude(0,4096)>0,"explicit lower-monitor restart resets ramp while preserving over-full-scale processing samples");
    }

    for(int outputs:{2,8,12})
    {
        auto fifo=std::make_unique<AuditionBuffer>();OutputRouting map;
        for(int ch=0;ch<std::min(outputs,8);++ch)map.fileToDevice[ch]=outputs-1-ch;
        expect(map.validate(8,outputs).isEmpty(),"explicit partial/full noncontiguous output mapping valid");
        fifo->gain.store(.25f);fifo->reset(8,map);Capture capture(8,4096);fifo->setObserver(&capture,0);
        juce::AudioBuffer<float> input(8,4096),out(outputs,4096);
        for(int ch=0;ch<8;++ch)juce::FloatVectorOperations::fill(input.getWritePointer(ch),.01f*(ch+1),4096);
        fifo->push(input.getArrayOfReadPointers(),4096);fifo->render(out.getArrayOfWritePointers(),outputs,4096);
        bool mapped=true,exact=true;
        for(int ch=0;ch<8;++ch){exact&=!std::memcmp(input.getReadPointer(ch),capture.data.getReadPointer(ch),4096*sizeof(float));if(map.fileToDevice[ch]>=0)mapped&=std::abs(out.getSample(map.fileToDevice[ch],4095)-input.getSample(ch,4095)*.25*(1-std::pow(.998,4096)))<1e-6;}
        expect(mapped&&exact,"injected output map receives only assigned channels; pre-monitor capture retains eight");
        if(outputs==2)expect(monitoringCoverage(systemILayout(),map).contains("Partial monitoring: 2 of 8"),"two outputs explicitly labelled partial monitoring without downmix");
        auto duplicate=map;duplicate.fileToDevice[1]=duplicate.fileToDevice[0];expect(duplicate.validate(8,outputs).isNotEmpty(),"duplicate output rejected");
        auto invalid=map;invalid.fileToDevice[0]=outputs;expect(invalid.validate(8,outputs).isNotEmpty(),"out-of-range route never silently reassigned");
    }
}


double milliseconds(){return juce::Time::getMillisecondCounterHiRes();}
uint64_t residentM7(){mach_task_basic_info_data_t info{};mach_msg_type_number_t count=MACH_TASK_BASIC_INFO_COUNT;if(task_info(mach_task_self(),MACH_TASK_BASIC_INFO,reinterpret_cast<task_info_t>(&info),&count)!=KERN_SUCCESS)throw std::runtime_error("resident memory query failed");return info.resident_size;}
double percentile95(std::vector<double> values){std::sort(values.begin(),values.end());return values[static_cast<size_t>(.95*(values.size()-1))];}
juce::var downstreamChecks(const measurement::Result& r,const juce::String& identity)
{
    using namespace proposal;
    const auto baseline=residentM7();auto began=milliseconds();auto a=analyzeProgramme(r,identity);const double core=milliseconds()-began;
    expect(a&&a->peaks.size()==8&&core<=2000&&a->storageBytes()<=64*1024*1024,"eight-channel three-hour diagnostic core and storage budgets");
    expect(a->frames==r.frames&&a->peaks.back().sampleFrame==r.frames-1,"late channel-eight peak remains at actual programme EOF");
    auto snapshot=std::make_shared<AnalysisSnapshot>();snapshot->context.sessionId="m7-downstream-performance";
    snapshot->context.media.sha256=identity;snapshot->context.media.channels=8;snapshot->context.media.sampleRate=r.sampleRate;snapshot->context.media.frames=r.frames;
    snapshot->context.range={0,r.frames};snapshot->context.fileLayout=systemILayout();snapshot->context.fileLayout.confirmed=true;
    snapshot->measurementVersion=measurement::version;snapshot->measuredContentHash=identity;snapshot->domain=MeasurementDomain::source;snapshot->detail=r;
    snapshot->measured={r.integrated.value,r.lra.value,r.truePeak.value,r.samplePeak.value};snapshot->identity=OfflineAnalysisService::identityFor(snapshot->context,snapshot->domain,identity);
    double workerMs=0;
    {DiagnosticsController worker;began=milliseconds();worker.select(snapshot,snapshot->context,snapshot->domain);do{worker.poll();std::this_thread::sleep_for(std::chrono::milliseconds(1));}while(worker.busy()&&milliseconds()-began<5000);workerMs=milliseconds()-began;expect(worker.result()!=nullptr&&workerMs<=2000,"eight-channel diagnostic worker including full snapshot binding within two seconds");}
    snapshot.reset();

    std::vector<double> plots,paints,edits;MeasurementTimeline view;view.setSize(1040,300);view.setAnalysis(a);
    juce::Image image(juce::Image::ARGB,1040,300,true);
    for(int n=0;n<80;++n)
    {
        const double begin=(n%60)*150.*r.sampleRate,end=std::min(double(r.frames),begin+1800*r.sampleRate);
        began=milliseconds();const auto plot=prepareTimeline(*a,begin,end,1000);plots.push_back(milliseconds()-began);expect(plot.momentary.size()==1000,"long diagnostic plot bounded by viewport width");
        view.setViewRange(begin,end);juce::Graphics g(image);began=milliseconds();view.paintEntireComponent(g,false);paints.push_back(milliseconds()-began);
    }
    expect(percentile95(plots)<=50&&percentile95(paints)<=100,"eight-channel long viewport and paint responsiveness");
    expect(view.inspectSeconds(10799.9),"eight-channel final-second numerical inspection");
    const auto reduced=prepareTimeline(*a,0,double(r.frames),1);
    expect(reduced.momentary[0].included.count==r.integratedIncluded&&reduced.momentary[0].excluded.count==r.momentary.size()-r.integratedIncluded,"long plot gate inclusion remains exact after reduction");
    SectionDocument document;document.bind(identity,r.sampleRate,r.frames);juce::String error;
    for(int n=1;n<1000;++n)expect(document.split(document.map().sections.back().id,r.frames*n/1000,error),"create bounded 1000-section annotation benchmark");
    began=milliseconds();auto sections=analyzeSections(*a,document.map(),identity);const auto sectionMs=milliseconds()-began;
    const auto workspace=sections->storageBytes()+document.map().sections.capacity()*sizeof(ProgrammeSection)+a->momentary.size()*sizeof(double)+128000;
    expect(sectionMs<=2000&&workspace<=16*1024*1024&&sections->sections.size()==1000,"eight-channel section work within two seconds and 16 MiB");
    expect(!sections->sections.back().programmeMaxima.empty(),"late programme peak belongs to final section");
    for(int n=0;n<80;++n){const auto& section=sections->sections[n*12];began=milliseconds();expect(document.moveEnd(section.section.id,section.section.end-1,error)&&document.moveEnd(section.section.id,section.section.end,error),"section edit remains reversible");view.highlight(section.section.begin,section.section.end);edits.push_back(milliseconds()-began);}
    expect(percentile95(edits)<=100,"eight-channel section edit responsiveness");
    std::atomic<bool> entered{false},cancel{false};began=milliseconds();auto task=std::async(std::launch::async,[&]{return analyzeProgramme(r,"cancel",[&]{entered=true;return cancel.load();});});while(!entered.load())std::this_thread::yield();cancel=true;expect(!task.get()&&milliseconds()-began<=1000,"active long eight-channel diagnostic cancellation within one second");
    entered=false;cancel=false;began=milliseconds();auto task2=std::async(std::launch::async,[&]{return analyzeSections(*a,document.map(),identity,[&]{entered=true;return cancel.load();});});while(!entered.load())std::this_thread::yield();cancel=true;expect(!task2.get()&&milliseconds()-began<=1000,"active long eight-channel section cancellation within one second");
    struct rusage usage{};getrusage(RUSAGE_SELF,&usage);const auto extra=std::max<int64_t>(0,static_cast<int64_t>(usage.ru_maxrss)-static_cast<int64_t>(baseline));expect(extra<=128*1024*1024,"diagnostic pipeline additional peak RSS within 128 MiB");
    auto* out=new juce::DynamicObject();out->setProperty("core_ms",core);out->setProperty("worker_ms",workerMs);out->setProperty("storage_bytes",static_cast<juce::int64>(a->storageBytes()));out->setProperty("plot_p95_ms",percentile95(plots));out->setProperty("paint_p95_ms",percentile95(paints));out->setProperty("edit_p95_ms",percentile95(edits));out->setProperty("section_ms",sectionMs);out->setProperty("section_workspace_bytes",static_cast<juce::int64>(workspace));out->setProperty("additional_peak_rss_bytes",static_cast<juce::int64>(extra));return juce::var(out);
}
void makeLong(const juce::File& file,int rate,int seconds)
{
    expect(!file.exists()&&file.getParentDirectory().createDirectory().wasOk(),"fresh generated long fixture path");
    juce::FlacAudioFormat format;auto stream=file.createOutputStream();std::unique_ptr<juce::AudioFormatWriter> writer(format.createWriterFor(stream.get(),rate,8,16,{},0));expect(writer!=nullptr,"native eight-channel FLAC long fixture writer");stream.release();
    const int64_t total=static_cast<int64_t>(rate)*seconds;juce::AudioBuffer<float> block(8,4096);
    for(int64_t first=0;first<total;first+=4096)
    {
        const int count=static_cast<int>(std::min<int64_t>(4096,total-first));
        for(int ch=0;ch<8;++ch)for(int i=0;i<count;++i)
        {
            const int64_t frame=first+i;const int period=64+32*ch;const int phase=static_cast<int>(frame%period);
            const int triangle=phase<period/2?phase:period-phase;
            const int level=1+static_cast<int>((frame/(rate*60LL))%3);
            float x=static_cast<float>((triangle-period/4)*(ch+1)*level)/32768.f;
            if(frame==total-1-(7-ch)*17)x=.25f+.04f*ch;
            block.setSample(ch,i,x);
        }
        if(!writer->writeFromAudioSampleBuffer(block,0,count))throw std::runtime_error("long fixture write failure");
    }
    writer.reset();std::cout<<"Generated "<<file.getFullPathName()<<" frames="<<total<<" bytes="<<file.getSize()<<std::endl;
}
void longProgrammeM7(const juce::File& file,bool measureProgramme=true)
{
    const auto began=milliseconds();auto request=validated(file);const auto importMs=milliseconds()-began;
    expect(request.context.media.channels==8&&request.context.media.sampleRate==48000&&request.context.media.frames==48000LL*10800,"actual three-hour eight-channel 48 kHz native FLAC fully validated");
    auto* evidence=new juce::DynamicObject();longEvidence=juce::var(evidence);evidence->setProperty("source_sha256",request.context.media.sha256);evidence->setProperty("encoded_bytes",request.context.media.bytes);evidence->setProperty("decoded_frames",request.context.media.frames);evidence->setProperty("import_ms",importMs);
    auto device=std::make_unique<Device>();auto* fake=device.get();auto audio=std::make_unique<StandaloneAudition>(std::move(device));auto map=routes(8);audio->synchronize(request,map);
    const auto baseline=residentM7();const auto preparation=milliseconds();expect(audio->prepare(request).isEmpty(),"prepare real three-hour eight-channel file while closed");
    uint64_t peak=baseline;std::vector<double> pollTimes;
    while(audio->preparing()&&milliseconds()-preparation<600000){const auto start=milliseconds();audio->poll();pollTimes.push_back(milliseconds()-start);peak=std::max(peak,residentM7());std::this_thread::sleep_for(std::chrono::milliseconds(20));}
    expect(audio->preparedCopy()&&!audio->preparing()&&fake->opens==0,"three-hour source copy prepared without device");
    const auto* copy=audio->preparedCopy();evidence->setProperty("prepare_ms",milliseconds()-preparation);evidence->setProperty("prepare_poll_p95_ms",percentile95(pollTimes));evidence->setProperty("additional_owned_workspace_rss_bytes",static_cast<juce::int64>(peak-baseline+sizeof(AuditionBuffer)));
    expect(peak-baseline+sizeof(AuditionBuffer)<=16*1024*1024&&percentile95(pollTimes)<=100,"long preparation workspace and responsiveness budgets");
    if(measureProgramme) for(double gain:{0.,-3.305074118944})
    {
        request.domain=gain==0?MeasurementDomain::source:MeasurementDomain::currentState;request.context.processing.staticGainDb=gain;++request.context.ticket.sessionRevision;
        const auto start=milliseconds();JobControl job;
        auto measurement=std::async(std::launch::async,[&]{return OfflineAnalysisService().analyze(request.context,request.domain,job);});
        while(measurement.wait_for(std::chrono::seconds(15))!=std::future_status::ready)
            std::cout<<"Three-hour "<<(gain==0?"Source":"Current")<<" progress="<<job.progress.load()<<" elapsed_s="<<(milliseconds()-start)/1000<<std::endl;
        auto result=measurement.get();
        expect(result.value.has_value()&&result.value->detail.frames==48000LL*10800,"complete real three-hour eight-channel measurement");
        const auto& r=result.value->detail;const auto key=gain==0?"source":"current";
        auto* run=new juce::DynamicObject();run->setProperty("elapsed_ms",milliseconds()-start);run->setProperty("identity",result.value->identity);run->setProperty("gain_db",gain);run->setProperty("integrated_lufs",*r.integrated.value);run->setProperty("true_peak_dbtp",*r.truePeak.value);run->setProperty("lra_lu",r.lra.value?juce::var(*r.lra.value):juce::var());
        juce::Array<juce::var> peaks;for(int ch=0;ch<8;++ch){expect(r.samplePeakFrames[ch]==r.frames-1-(7-ch)*17,"all eight late peak locations retained");auto* p=new juce::DynamicObject();p->setProperty("channel",ch);p->setProperty("sample_peak",r.channelSamplePeaks[ch]);p->setProperty("true_peak",r.channelTruePeaks[ch]);p->setProperty("sample_frame",static_cast<juce::int64>(r.samplePeakFrames[ch]));p->setProperty("true_peak_eighth_frame",static_cast<juce::int64>(r.truePeakEighthFrames[ch]));peaks.add(juce::var(p));}run->setProperty("peaks",peaks);run->setProperty("diagnostics",downstreamChecks(r,result.value->identity));evidence->setProperty(key,juce::var(run));
        std::cout<<"Three-hour "<<key<<" measurement complete in "<<(milliseconds()-start)/1000<<" seconds"<<std::endl;
    }
    expect(fake->opens==0,"both full long measurements complete device-closed");
    expect(audio->connect(48000).isEmpty(),"explicit injected device for long seek captures");audio->setMonitorGainDb(-40);
    juce::AudioFormatManager fm;registerSourceFormats(fm);std::unique_ptr<juce::AudioFormatReader> reader(fm.createReaderFor(file));
    juce::Array<juce::var> seekCases;
    for(double gain:{0.,-3.305074118944})
    for(int64_t frame:{0LL,123456789LL,48000LL*10799,48000LL*10800-1})
    {
        request.domain=gain==0?MeasurementDomain::source:MeasurementDomain::currentState;request.context.processing.staticGainDb=gain;++request.context.ticket.sessionRevision;
        audio->synchronize(request,map);expect(audio->seek(frame,request,map).isEmpty(),"long seek uses exact audio-source frame");
        const int count=static_cast<int>(std::min<int64_t>(256,request.context.media.frames-frame));Capture captured(8,count,frame);audio->setCaptureObserver(&captured);
        const auto seekStart=milliseconds();expect(audio->play(request,map).isEmpty(),"long seek play");awaitPlaying(*audio,*fake);const auto seekMs=milliseconds()-seekStart;
        expect(seekMs<=2000,"prepared three-hour FLAC seek including final sample buffers within two seconds");
        juce::AudioBuffer<float> output(12,count);fake->tick(output);audio->stop();
        juce::AudioBuffer<float> expected(8,count);reader->read(expected.getArrayOfWritePointers(),8,frame,count);StaticGainProcessor(request.context.processing.staticGainDb).process(expected.getArrayOfWritePointers(),8,count);
        bool same=captured.valid&&captured.next==frame+count;for(int ch=0;ch<8;++ch)same&=!std::memcmp(expected.getReadPointer(ch),captured.data.getReadPointer(ch),count*sizeof(float));
        expect(same&&audio->preparedCopy()==copy,"random and EOF seek samples match on all eight channels with reused copy");
        auto* row=new juce::DynamicObject();row->setProperty("frame",static_cast<juce::int64>(frame));row->setProperty("gain_db",gain);row->setProperty("buffer_ms",seekMs);row->setProperty("all_eight_exact",same);seekCases.add(juce::var(row));
        audio->setCaptureObserver(nullptr);
    }
    evidence->setProperty("seek_cases",seekCases);
    expect(audio->seek(request.context.media.frames-1,request,map).isEmpty(),"final-sample cancellation seek");
    const auto cancelStart=milliseconds();expect(audio->play(request,map).isEmpty(),"start final-sample prefill before cancellation");audio->stop();
    expect(milliseconds()-cancelStart<=2000&&!audio->active()&&!fake->isPlaying()&&audio->preparedCopy()==copy,"final-sample FLAC cancellation bounded, stopped and verified copy retained");
    audio->disconnect();
}


void videoIsolation(const juce::File& dir,const juce::File& movie)
{
    const auto file=dir.getChildFile("video-eight.wav");const auto pcm=signal(48000,8,96000);writePcm(file,pcm,48000,24);
    auto r=validated(file);auto injected=std::make_unique<Device>();auto* device=injected.get();device->outputCount=2;
    auto audio=std::make_unique<StandaloneAudition>(std::move(injected));OutputRouting map;map.fileToDevice[0]=0;map.fileToDevice[1]=1;
    expect(audio->connect(48000).isEmpty(),"explicit two-output injected device for video isolation");audio->setMonitorGainDb(-40);
    beginVideoReaderAudit();video::Player player;player.attach(movie);uint64_t generation=1;
    const auto picture=[&](video::Time time)
    {
        video::Player::View view;const auto deadline=milliseconds()+30000;
        do{view=player.update(time,generation);if(view.state==video::Player::State::failed)throw std::runtime_error(view.message.toStdString());if(view.picture&&view.picture->generation==generation)break;std::this_thread::sleep_for(std::chrono::milliseconds(1));}while(milliseconds()<deadline);
        expect(view.picture.has_value()&&view.asset&&view.queued<=8,"native video presents bounded generation-correct picture");
        expect(video::frameAt(view.asset->frames,time)==view.picture->index,"native video frame follows exact rational audio coordinate and signed offset");return view;
    };
    auto initial=picture({0,1});expect((initial.asset->embeddedAudioTracks>0)==(movie.getFileName()!="vfr-exact.mov"),"reference embedded-audio presence matches fixture definition (VFR fixture is video-only)");
    for(double gain:{0.,-3.305074118944})
    {
        r.domain=gain==0?MeasurementDomain::source:MeasurementDomain::currentState;r.context.processing.staticGainDb=gain;++r.context.ticket.sessionRevision;
        const auto id=r.identity();audio->synchronize(r,map);
        for(int64_t frame:{0LL,72001LL,95999LL})
        {
            const auto offset=frame==0?video::Time(-1001,30000):video::Time(1001,30000);++generation;
            picture(video::Time(frame,48000)-offset);
            expect(audio->seek(frame,r,map).isEmpty(),"native-video linked audition seek accepts original source coordinate");
            const int count=static_cast<int>(std::min<int64_t>(256,96000-frame));Capture capture(8,count,frame);audio->setCaptureObserver(&capture);expect(audio->play(r,map).isEmpty(),"Source/Current video isolation play");awaitPlaying(*audio,*device);
            juce::AudioBuffer<float> output(2,count);device->tick(output);audio->stop();
            juce::AudioFormatManager formats;registerSourceFormats(formats);std::unique_ptr<juce::AudioFormatReader> reader(formats.createReaderFor(file));juce::AudioBuffer<float> expected(8,count);reader->read(expected.getArrayOfWritePointers(),8,frame,count);StaticGainProcessor(gain).process(expected.getArrayOfWritePointers(),8,count);
            bool exact=capture.valid&&capture.next==frame+count;for(int ch=0;ch<8;++ch)exact&=!std::memcmp(expected.getReadPointer(ch),capture.data.getReadPointer(ch),count*sizeof(float));
            expect(exact&&r.identity()==id,"eight selected channels remain exact with H.264/ProRes and embedded audio excluded");
        }
    }
    player.clear();audio->disconnect();expect(videoReaderAuditCount(0)>0&&videoReaderAuditCount(1)==0,"runtime reader audit proves zero embedded-audio decode outputs");
    auto* item=new juce::DynamicObject();item->setProperty("file",movie.getFileName());item->setProperty("sha256",initial.asset->sha256);item->setProperty("embedded_audio_tracks",initial.asset->embeddedAudioTracks);item->setProperty("observed_video_outputs",videoReaderAuditCount(0));item->setProperty("observed_audio_outputs",videoReaderAuditCount(1));item->setProperty("selected_channels",8);item->setProperty("source_and_current_exact_capture",true);videoCases.add(juce::var(item));
}

void panelM7(const juce::File& dir)
{
    const auto file=dir.getChildFile("panel-71.wav");writePcm(file,signal(48000,8,96000),48000,24);auto r=validated(file);
    StandaloneSessionModel session;session.setValidatedSource(file,r.context.media);
    auto device=std::make_unique<Device>();device->outputCount=2;auto audio=std::make_unique<StandaloneAudition>(std::move(device));
    StandaloneAudioPanel panel(session,*audio,[]{});panel.setSize(700,680);panel.resized();panel.refresh();
    std::vector<juce::ComboBox*> roles;juce::ComboBox* preset=nullptr;int labels=0;
    const std::function<void(juce::Component&)> visit=[&](juce::Component& component)
    {
        if(auto* box=dynamic_cast<juce::ComboBox*>(&component))
        {if(box->getItemText(0)=="Unassigned")roles.push_back(box);if(box->getItemText(0).startsWith("L R C LFE Lrs"))preset=box;}
        if(auto* label=dynamic_cast<juce::Label*>(&component))if(label->getText().startsWith("Channel "))++labels;
        for(int i=0;i<component.getNumChildComponents();++i)visit(*component.getChildComponent(i));
    };visit(panel);
    expect(roles.size()==12&&labels==12&&preset&&preset->isEnabled(),"reusable twelve-channel UI uses Channel labels and enables explicit 7.1 presets");
    expect(!session.sourceDocument.layout.confirmed,"format interpretation still requires user role confirmation");
    preset->setSelectedId(2,juce::sendNotificationSync);panel.refresh();
    expect(session.sourceDocument.layout.roles[4]==C::leftSurroundSide&&session.sourceDocument.layout.roles[6]==C::leftSurroundRear
        &&session.sourceDocument.layout.manuallyInterpreted&&!session.sourceDocument.layout.confirmed,"named preset changes interpretation only with explicit unconfirmed manual provenance");
    preset->setSelectedId(1,juce::sendNotificationSync);panel.refresh();
    expect(roles[4]->isTextEditable(),"System I role fields accept bounded explicit aliases");
    roles[4]->setText("BL",juce::sendNotificationSync);expect(session.sourceDocument.layout.roles[4]==C::leftSurroundRear,"typed back-left alias explicitly identifies System I rear");
    roles[4]->setText("Ls",juce::sendNotificationSync);expect(!session.confirmLayout(),"ambiguous generic typed surround cannot confirm a 7.1 layout");
    roles[4]->setText("Lrs",juce::sendNotificationSync);expect(session.confirmLayout(),"canonical rear alias can be confirmed explicitly");
    session.initializeCurrentStateFromSource();session.selectAuditionDomain(MeasurementDomain::currentState);panel.refresh();
    bool unverified=false,partial=false;session.setOutputRoute(0,0);session.setOutputRoute(1,1);panel.refresh();
    const std::function<void(juce::Component&)> status=[&](juce::Component& component)
    {if(auto* label=dynamic_cast<juce::Label*>(&component)){unverified|=label->getText().contains("Current State")&&label->getText().contains("UNVERIFIED");partial|=label->getText().contains("Partial monitoring: 2 of 8");}for(int i=0;i<component.getNumChildComponents();++i)status(*component.getChildComponent(i));};status(panel);
    expect(unverified&&partial&&!audio->connected(),"eight-channel UI shows UNVERIFIED and explicit partial monitoring while device closed");
}

void createUiFixtures(const juce::File& dir)
{
    expect(!dir.exists()&&dir.createDirectory().wasOk(),"fresh M7 UI fixture directory");
    auto b=signal(48000,8,48000*8);writePcm(dir.getChildFile("system-i-mask.wav"),b,48000,24);
    writePcm(dir.getChildFile("system-i-unlabelled.wav"),b,48000,24,false,{});
    writeCompressedOrAiff(dir.getChildFile("system-i-default.flac"),b,48000,24);
    writeCompressedOrAiff(dir.getChildFile("system-i-manual.aiff"),b,48000,24);
    writePcm(dir.getChildFile("unsupported-sdds.wav"),b,48000,24,false,0xff);
    writePcm(dir.getChildFile("invalid-mask.wav"),b,48000,24,false,0x3f);
    auto six=signal(48000,6,48000*8);writePcm(dir.getChildFile("legacy-51.wav"),six,48000,24,false,0x3f);writeCompressedOrAiff(dir.getChildFile("legacy-51.flac"),six,48000,24);
}

void parityAndFormats(const juce::File& dir)
{
    for(int rate:{44100,48000,96000,192000})for(int format=0;format<9;++format)
    {
        const int bits=format<3?16+8*format:format==3?32:format<7?16+8*(format-4):16+8*(format-7);
        const auto file=dir.getChildFile(juce::String(rate)+"-"+juce::String(format)+(format<4?".wav":format<7?".aiff":".flac"));
        const auto b=signal(rate,8,rate/2+17);
        if(format<4)writePcm(file,b,rate,bits,format==3);else writeCompressedOrAiff(file,b,rate,bits);
        auto r=validated(file,format>=4&&format<7);expect(isSystemI(r.context.fileLayout),"format-aware or explicitly assigned System I layout");
        auto injected=std::make_unique<Device>();auto* device=injected.get();auto audio=std::make_unique<StandaloneAudition>(std::move(injected));
        auto map=routes(8);audio->synchronize(r,map);expect(audio->prepare(r).isEmpty(),"eight-channel prepare with device closed");waitPreparation(*audio);expect(device->opens==0&&audio->preparedCopy(),"offline preparation opens no audio device");const auto* copy=audio->preparedCopy();
        expect(audio->connect(rate).isEmpty(),"explicit injected device connects at source rate");audio->setMonitorGainDb(-40);
        for(double gain:{0.,-6.020599913279624,6.020599913279624,-3.305074118944,-60.,24.})
        {
            r.domain=MeasurementDomain::currentState;r.context.processing.staticGainDb=gain;++r.context.ticket.sessionRevision;audio->synchronize(r,map);
            Capture offline(8,b.getNumSamples()),playback(8,b.getNumSamples());JobControl job;auto result=OfflineAnalysisService(&offline).analyze(r.context,r.domain,job);expect(result.value.has_value(),"eight-channel offline Current measures");
            expect(audio->setCaptureObserver(&playback)&&audio->play(r,map,true).isEmpty(),"actual eight-channel audition starts explicitly");awaitPlaying(*audio,*device);drain(*audio,*device,997);
            expect(identical(offline,playback)&&audio->preparedCopy()==copy,"all eight pre-monitor channels bit-exact to offline and verified copy reused");
            recordSamples(playback,rate,gain,"M7 format "+juce::String(format),result.value->detail);
        }
        audio->disconnect();
    }
}
}
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI init;bool passed=false;
    try
    {
        for(int i=1;i+1<argc;++i)if(juce::String(argv[i])=="--samples"){sampleDirectory=juce::File(argv[i+1]);expect(!sampleDirectory.exists()&&sampleDirectory.createDirectory().wasOk(),"fresh M7 sample evidence directory");}
        juce::File longSource,seekSource,generate,uiDirectory;std::vector<juce::File> movies;int generateRate=48000,generateSeconds=10800;
        for(int i=1;i+1<argc;++i){const juce::String key(argv[i]);if(key=="--long-source")longSource=juce::File(argv[i+1]);if(key=="--long-seeks")seekSource=juce::File(argv[i+1]);if(key=="--ui-fixtures")uiDirectory=juce::File(argv[i+1]);if(key=="--video")movies.emplace_back(juce::String(argv[i+1]));if(key=="--generate")generate=juce::File(argv[i+1]);if(key=="--rate")generateRate=juce::String(argv[i+1]).getIntValue();if(key=="--seconds")generateSeconds=juce::String(argv[i+1]).getIntValue();}
        if(uiDirectory!=juce::File())createUiFixtures(uiDirectory);
        else if(generate!=juce::File())makeLong(generate,generateRate,generateSeconds);
        else if(longSource!=juce::File())longProgrammeM7(longSource);
        else if(seekSource!=juce::File())longProgrammeM7(seekSource,false);
        else if(!movies.empty()){Temp temp;for(const auto& movie:movies)videoIsolation(temp.directory,movie);}
        else{Temp temp;layoutChecks();weightingChecks();peakAndPermutationChecks();metadataAndLegacy(temp.directory);mutedPeakWorkflow(temp.directory);routingCapture();panelM7(temp.directory);parityAndFormats(temp.directory);longEvidence=downstreamChecks(longFixture(),"synthetic-three-hour-eight-channel-192k-not-PCM");}
        passed=true;
    }
    catch(const std::exception& e){std::cerr<<"M7 FAIL: "<<e.what()<<'\n';}
    auto* report=new juce::DynamicObject();report->setProperty("suite","M7 System I PCM 7.1");report->setProperty("passed",passed);report->setProperty("checks",checks);report->setProperty("sample_cases",cases);report->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);
    report->setProperty("profile_digest",proposal::ebuR128ProgrammeQc().digest());report->setProperty("layout_capability",systemICapability);report->setProperty("mapping_version",layoutMappingVersion);
    report->setProperty("long_programme",longEvidence);report->setProperty("video_cases",videoCases);
    report->setProperty("physical_71_listening","pending; not established by injected devices");
#if defined(__aarch64__) || defined(__arm64__)
    report->setProperty("architecture","arm64");
#else
    report->setProperty("architecture","x86_64");
#endif
    for(int i=1;i+1<argc;++i)if(juce::String(argv[i])=="--report")juce::File(argv[i+1]).replaceWithText(juce::JSON::toString(juce::var(report),true));
    std::cout<<"M7 multichannel: "<<checks.size()<<" checks, "<<(passed?"PASS":"FAIL")<<'\n';return passed?0:1;
}

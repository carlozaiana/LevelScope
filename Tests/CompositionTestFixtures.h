#pragma once
#include "AuditionTestSupport.h"
#include "Standalone/StandaloneExport.h"
#include "Standalone/SessionDocument.h"
#include "Standalone/VerifiedSource.h"
#include "Standalone/StandaloneDiagnostics.h"
#include <fstream>
namespace ex=levelscope::standalone::exporting;
namespace sf=levelscope::standalone::sessionfile;
namespace
{
void floatFile(const juce::File& file,int rate,const juce::AudioBuffer<float>& data)
{
    expect(!file.exists(),"oracle fixture never overwrites existing media");auto out=file.createOutputStream();expect(out!=nullptr,"float oracle stream");
    const int bytes=data.getNumChannels()*data.getNumSamples()*4;
    out->write("RIFF",4);out->writeInt(36+bytes);out->write("WAVEfmt ",8);out->writeInt(16);out->writeShort(3);out->writeShort(static_cast<short>(data.getNumChannels()));
    out->writeInt(rate);out->writeInt(rate*data.getNumChannels()*4);out->writeShort(static_cast<short>(data.getNumChannels()*4));out->writeShort(32);out->write("data",4);out->writeInt(bytes);
    for(int n=0;n<data.getNumSamples();++n)for(int c=0;c<data.getNumChannels();++c){uint32_t bits;const float x=data.getSample(c,n);std::memcpy(&bits,&x,4);out->writeInt(static_cast<int>(bits));}
    out->flush();expect(out->getStatus().wasOk(),"float oracle finalized");
}
ProgrammeComponent component(const juce::File& root,const juce::String& name,int rate,const juce::AudioBuffer<float>& data,bool grouped)
{
    ProgrammeComponent p;p.id=name;p.purpose=name=="me"?"music-effects":"dialogue";p.layout=layoutFor(data.getNumChannels());
    JobControl job;
    if(grouped&&data.getNumChannels()>1)
    {
        ProgrammeAssembly group;const auto roles=assemblyRoles(data.getNumChannels());
        for(int c=0;c<data.getNumChannels();++c)
        {
            juce::AudioBuffer<float> mono(1,data.getNumSamples());mono.copyFrom(0,0,data,c,0,data.getNumSamples());
            auto file=root.getChildFile(name+"-"+juce::String(c)+".wav");floatFile(file,rate,mono);auto v=validateMedia(file,job);expect(v.succeeded(),"group component member verified");
            group.members.push_back({"audio-"+juce::String(c),file,v.media,roles[static_cast<size_t>(c)],0});
        }
        expect(group.confirmAlignment(),"test generator explicitly confirms component channel alignment");p.assembly=group;p.media=group.summary();p.layout=group.layout();
    }
    else
    {p.file=root.getChildFile(name+".wav");floatFile(p.file,rate,data);auto v=validateMedia(p.file,job);expect(v.succeeded(),"interleaved component verified");p.media=v.media;}
    if(data.getNumChannels()==1){p.layout.count=1;p.layout.roles[0]=juce::AudioChannelSet::centre;p.layout.confirmed=true;}
    return p;
}
struct CompositionFixture
{
    ProgrammeComposition recipe;
    juce::File oracle;
    juce::AudioBuffer<float> expected;
};
CompositionFixture fixtureComposition(const juce::File& root,int rate,CompositionPreset preset,bool grouped,int frames)
{
    expect(root.createDirectory().wasOk(),"composition fixture directory");const int channels=preset==CompositionPreset::stereoStereo?2:6;
    const int dialogueChannels=preset==CompositionPreset::surroundMonoCentre?1:channels;
    juce::AudioBuffer<float> me(channels,frames),dialogue(dialogueChannels,frames);
    for(int c=0;c<channels;++c)for(int n=0;n<frames;++n)me.setSample(c,n,static_cast<float>(std::llround((200+c*80)*std::sin((n+c*11)*2*3.14159265358979323846/97)))/32768.f);
    for(int c=0;c<dialogueChannels;++c)for(int n=0;n<frames;++n)dialogue.setSample(c,n,static_cast<float>(std::llround((100+c*50)*std::sin((n+c*7)*2*3.14159265358979323846/53)))/32768.f);
    CompositionFixture f;f.recipe.preset=preset;f.recipe.me=component(root,"me",rate,me,grouped);f.recipe.dialogue=component(root,"dialogue",rate,dialogue,grouped);
    expect(f.recipe.confirmComplementaryBalancedAlignment(),"generator confirms complementary balanced aligned composition");
    f.expected.setSize(channels,frames);
    // Independent preset oracle: no production mapping or summing helper.
    for(int c=0;c<channels;++c)for(int n=0;n<frames;++n)
    {
        const bool absent=dialogueChannels==1&&c!=2;const int dc=dialogueChannels==1?0:c;
        const float result=absent?me.getSample(c,n):static_cast<float>(static_cast<double>(me.getSample(c,n))+static_cast<double>(dialogue.getSample(dc,n)));
        f.expected.setSample(c,n,result);
    }
    f.oracle=root.getChildFile("independent-composed-oracle.wav");floatFile(f.oracle,rate,f.expected);return f;
}
StandaloneSessionModel compositionModel(const ProgrammeComposition& recipe)
{
    StandaloneSessionModel s;expect(s.setValidatedComposition(recipe),"adopt explicit composition atomically");expect(s.initializeCurrentStateFromSource(),"initialize composed Current");
    s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);return s;
}
AnalysisSnapshot measuredComposition(StandaloneSessionModel& s,MeasurementDomain domain,ProcessedAudioObserver* capture=nullptr)
{
    expect(s.job.begin(JobStage::analysis,s.revision),"begin composition measurement");auto c=s.measurementContext(domain);c.ticket=s.job.ticket;JobControl job;
    auto result=OfflineAnalysisService(capture).analyze(c,domain,job);if(!result.value)std::cerr<<result.message<<'\n';expect(result.value.has_value(),"composition measurement completes");
    expect(s.job.finish(c.ticket,s.revision,result.error,result.message)&&s.acceptMeasurement(*result.value),"accept exact composition measurement");return *result.value;
}
bool sameSamples(const juce::AudioBuffer<float>& a,const juce::AudioBuffer<float>& b)
{
    if(a.getNumChannels()!=b.getNumChannels()||a.getNumSamples()!=b.getNumSamples())return false;
    for(int c=0;c<a.getNumChannels();++c)if(std::memcmp(a.getReadPointer(c),b.getReadPointer(c),static_cast<size_t>(a.getNumSamples())*4))return false;return true;
}
}

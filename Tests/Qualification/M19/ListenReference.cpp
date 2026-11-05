// User-requested repeat: unchanged 99.4 s reference, then four 8 s captures, monitor 0 dB.
#include "Standalone/StandaloneAudition.h"
#include <juce_cryptography/juce_cryptography.h>
#include <iostream>
#include <thread>
#include <chrono>
using namespace levelscope::standalone;
int main(int argc,char** argv)
{
    if(argc!=5||juce::String(argv[1])!="--allow-audio-saffire-zero-reference")return 2;
    juce::ScopedJuceInitialiser_GUI gui;const juce::String name(argv[2]);const juce::File manifestFile(argv[3]),report(argv[4]);
    if(name!="Saffire Pro 40 (Experimental)"||report.exists())return 2;
    StandaloneAudition audition;juce::Array<juce::var> clips;juce::var result(new juce::DynamicObject());auto* r=result.getDynamicObject();bool passed=false;
    try
    {
        auto require=[](bool ok,const juce::String& error){if(!ok)throw std::runtime_error(error.toStdString());};
        const auto manifest=juce::JSON::parse(manifestFile.loadFileAsString());require(manifest["clips"].isArray()&&manifest["clips"].getArray()->size()==5&&int(manifest["rate"])==48000&&int(manifest["monitor_db"])==0,"Exact coordinated sequence required");
        std::vector<MediaValidationResult> media;const auto& sequence=*manifest["clips"].getArray();
        for(int i=0;i<5;++i)
        {
            const auto& row=sequence[i];juce::File file(row["path"].toString());require(juce::SHA256(file).toHexString()==row["sha256"].toString(),"Capture hash changed");JobControl job;media.push_back(validateMedia(file,job));
            require(media.back().succeeded()&&media.back().media.sampleRate==48000&&media.back().media.channels==2&&media.back().media.frames==(i==0?4771200:384000),"Only coordinated reference and stereo captures allowed");
        }
        int type=-1;for(int n=0;n<audition.deviceTypes().size();++n)if(audition.scanOutputs(n).contains(name)){require(type==-1,"Ambiguous device");type=n;}
        require(type>=0&&audition.selectDevice(type,name).isEmpty()&&audition.outputNames().size()>=2,"Approved output pair unavailable");
        ChannelLayout roles;roles.count=2;roles.confirmed=true;roles.roles[0]=juce::AudioChannelSet::left;roles.roles[1]=juce::AudioChannelSet::right;OutputRouting routing;routing.fileToDevice[0]=0;routing.fileToDevice[1]=1;
        for(int i=0;i<5;++i)
        {
            const auto& row=sequence[i];const auto& m=media[size_t(i)].media;audition.sourceChanged(m);if(!audition.connected())require(audition.connect(48000).isEmpty(),"Cannot open approved 48 kHz output");audition.setMonitorGainDb(0);
            std::cout<<"PLAY "<<row["label"].toString()<<" | monitor 0 dB | "<<m.frames/48000.<<" s\n"<<std::flush;
            require(audition.play(juce::File(row["path"].toString()),m,roles,routing,true).isEmpty(),"Cannot start explicit capture");auto start=std::chrono::steady_clock::now();
            while(audition.active())
            {
                audition.poll();std::this_thread::sleep_for(std::chrono::milliseconds(5));
                require(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<m.frames/48000.+15,"Bounded playback timeout");
            }
            require(audition.getState()==AuditionState::finished&&audition.underruns()==0&&audition.positionFrames()==m.frames,"Playback failure, underrun or wrong EOF");
            auto receipt=juce::var(new juce::DynamicObject());auto* o=receipt.getDynamicObject();o->setProperty("label",row["label"]);o->setProperty("sha256",row["sha256"]);o->setProperty("frames",juce::int64(audition.positionFrames()));o->setProperty("underruns",juce::int64(audition.underruns()));clips.add(receipt);
            audition.stop();std::cout<<"STOP\n"<<std::flush;std::this_thread::sleep_for(std::chrono::seconds(i==0?6:3));
        }
        passed=true;
    }
    catch(const std::exception& e){r->setProperty("failure",e.what());std::cerr<<e.what()<<'\n';}
    audition.stop();audition.disconnect();r->setProperty("passed",passed);r->setProperty("device",name);r->setProperty("rate",48000);r->setProperty("monitor_db",0);r->setProperty("outputs","1/2");r->setProperty("left_right_only",true);r->setProperty("clips",clips);r->setProperty("device_closed",!audition.connected());r->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);r->setProperty("perceptual_result","Awaiting user observations; no full 5.1 listening claim");
    return report.replaceWithText(juce::JSON::toString(result))?(passed?0:1):2;
}

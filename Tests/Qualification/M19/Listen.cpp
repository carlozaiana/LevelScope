// Optional coordinated listening launcher; never invoked by the silent runner.
#include "Standalone/StandaloneAudition.h"
#include <juce_cryptography/juce_cryptography.h>
#include <iostream>
#include <thread>
#include <chrono>
using namespace levelscope::standalone;
int main(int argc,char** argv)
{
    // Exact opt-in, device, clip directory and fresh receipt are mandatory.
    if(argc!=5||juce::String(argv[1])!="--allow-audio-saffire-minus5")return 2;
    juce::ScopedJuceInitialiser_GUI gui;const juce::String deviceName(argv[2]);const juce::File directory(argv[3]),report(argv[4]);
    if(deviceName!="Saffire Pro 40 (Experimental)"||report.exists())return 2;
    juce::var record(new juce::DynamicObject());auto* o=record.getDynamicObject();juce::Array<juce::var> clips;bool passed=false;
    StandaloneAudition audition;
    try
    {
        auto require=[](bool condition,const juce::String& error){if(!condition)throw std::runtime_error(error.toStdString());};
        const auto manifest=juce::JSON::parse(directory.getChildFile("manifest.json").loadFileAsString());require(manifest.isObject()&&manifest["published_export"]==juce::var(false),"Explicit diagnostic-capture manifest required");
        require(manifest["captures"].isArray()&&manifest["captures"].getArray()->size()==4,"Exactly four identified captures required");
        int selected=-1;for(int type=0;type<audition.deviceTypes().size();++type)if(audition.scanOutputs(type).contains(deviceName)){require(selected==-1,"Ambiguous device");selected=type;}
        require(selected>=0,"Approved Saffire device unavailable");require(audition.selectDevice(selected,deviceName).isEmpty(),"Cannot select approved device");require(audition.outputNames().size()>=2,"Outputs 1/2 unavailable");
        for(auto region:{"onset","return"})for(bool after:{false,true})
        {
            const auto name=juce::String(region)+(after?"-candidate":"-baseline");auto file=directory.getChildFile(name+".wav");
            bool partial=false;if(directory.getChildFile(name+"-LEFT-RIGHT-ONLY.wav").existsAsFile()){partial=true;file=directory.getChildFile(name+"-LEFT-RIGHT-ONLY.wav");}
            const auto hash=juce::SHA256(file).toHexString();bool identified=false;
            for(const auto& row:*manifest["captures"].getArray())identified|=partial?(row["optional_partial_capture"].toString()==file.getFileName()&&row["partial_sha256"].toString()==hash):(row["file"].toString()==file.getFileName()&&row["sha256"].toString()==hash);
            require(identified,"Clip does not match diagnostic manifest");JobControl job;auto media=validateMedia(file,job);
            require(media.succeeded()&&media.media.channels==2&&media.media.sampleRate==48000&&media.media.frames==384000,"Only eight-second stereo 48 kHz captures allowed");
            ChannelLayout roles;roles.count=2;roles.roles[0]=juce::AudioChannelSet::left;roles.roles[1]=juce::AudioChannelSet::right;roles.confirmed=true;
            OutputRouting routing;routing.fileToDevice[0]=0;routing.fileToDevice[1]=1;
            audition.sourceChanged(media.media);if(!audition.connected())require(audition.connect(48000).isEmpty(),"Cannot open approved output at 48 kHz");audition.setMonitorGainDb(-5);
            std::cout<<"PLAY "<<file.getFileName()<<" | outputs 1/2 | monitor -5 dB\n"<<std::flush;
            require(audition.play(file,media.media,roles,routing,true).isEmpty(),"Cannot start explicit capture");const auto start=std::chrono::steady_clock::now();
            while(audition.active())
            {
                audition.poll();std::this_thread::sleep_for(std::chrono::milliseconds(5));
                require(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<15,"Playback did not finish within bounded timeout");
            }
            require(audition.getState()==AuditionState::finished&&audition.underruns()==0&&audition.positionFrames()==384000,"Playback failure, underrun or wrong EOF");
            auto row=juce::var(new juce::DynamicObject());row.getDynamicObject()->setProperty("file",file.getFileName());row.getDynamicObject()->setProperty("sha256",hash);row.getDynamicObject()->setProperty("frames",juce::int64(audition.positionFrames()));row.getDynamicObject()->setProperty("underruns",juce::int64(audition.underruns()));row.getDynamicObject()->setProperty("partial_left_right_only",partial);clips.add(row);
            audition.stop();std::cout<<"STOP\n"<<std::flush;std::this_thread::sleep_for(std::chrono::seconds(3));
        }
        passed=true;
    }
    catch(const std::exception& e){o->setProperty("failure",e.what());std::cerr<<e.what()<<'\n';}
    audition.stop();audition.disconnect();o->setProperty("passed",passed);o->setProperty("device",deviceName);o->setProperty("rate",48000);o->setProperty("monitor_db",-5);o->setProperty("output_pair","1/2");o->setProperty("clips",clips);o->setProperty("device_closed",!audition.connected());o->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);o->setProperty("perceptual_result","Awaiting user observation; transport success is not perceptual transparency");
    return report.replaceWithText(juce::JSON::toString(record))?(passed?0:1):2;
}

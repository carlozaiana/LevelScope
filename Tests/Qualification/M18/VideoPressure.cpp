#define main unusedShortAcceptanceMain
#include "../../Profiling/ShortAcceptance/ShortAcceptance.cpp"
#undef main
#include <mach/mach.h>
namespace {
uint64_t resident(){mach_task_basic_info_data_t info{};mach_msg_type_number_t n=MACH_TASK_BASIC_INFO_COUNT;return task_info(mach_task_self(),MACH_TASK_BASIC_INFO,reinterpret_cast<task_info_t>(&info),&n)==KERN_SUCCESS?info.resident_size:0;}
void settle(Main& ui){
    until("picture did not settle at stopped audio clock",[&]{ui.videoPanel.synchronize();const auto& v=ui.videoPanel.view;if(!v.asset)return false;
        const auto time=video::Time(ui.audition.positionFrames(),48000)-ui.session.referenceVideo.startOnAudio;
        const auto target=video::frameAt(v.asset->frames,time);
        return target? v.picture&&v.picture->index==*target:!v.picture;});
}
void qualifyVideo(Main& ui,const juce::File& saved,const juce::File& fixtures,const juce::File& out,juce::var& report){
    auto device=std::make_unique<Device>();auto* output=device.get();output->outputCount=6;ui.audition.device=std::move(device);
    click(ui.openSessionButton);choose(ui.sessionChooser.get(),saved);until("video fixture session did not reopen",[&]{return !ui.documentResult.valid();});
    auto& s=ui.session;expect(s.source.hasMeasurement&&s.currentState.hasMeasurement&&s.lastAppliedExplanation,"video uses accepted applied section session");
    const auto sourceIdentity=s.source.analysis->identity,currentIdentity=s.currentState.analysis->identity,history=s.lastAppliedExplanation->text();
    ui.inspectDomain.setSelectedId(2,juce::sendNotificationSync);
    click(ui.videoButton);expect(ui.videoWindow&&ui.videoWindow->isVisible(),"native reference window shown alongside timeline");
    const auto initial=resident();uint64_t peak=initial;juce::Array<juce::var> rows;
    for(auto name:{"h264-cfr.mov","h264-vfr.mov","422-25-1.mov","lt-30000-1001.mov"})for(int domain:{1,2}){
        click(ui.videoPanel.attach);choose(ui.videoPanel.chooser.get(),fixtures.getChildFile(name));
        until("reference asset did not load",[&]{ui.videoPanel.synchronize();return bool(ui.videoPanel.view.asset);});
        ui.auditionDomainBox.setSelectedId(domain,juce::sendNotificationSync);
        int activePaints=0,partialPaints=0;double longestPaint=0;
        // Actual main-window measurement holds the reference picture while the
        // progressive timeline changes. Normal UI forbids concurrent audition.
        click(ui.remeasureCurrentStateButton);
        until("analysis with reference picture did not finish",[&]{@autoreleasepool {
            if(s.job.isActive()){const auto t=Clock::now();auto a=ui.timeline.createComponentSnapshot(ui.timeline.getLocalBounds());auto b=ui.videoPanel.createComponentSnapshot(ui.videoPanel.getLocalBounds());longestPaint=std::max(longestPaint,since(t)*1000);++activePaints;if(!ui.timeline.waveformAccepted&&std::any_of(ui.timeline.displayedWaveform().begin(),ui.timeline.displayedWaveform().end(),[](const auto& c){return c.available();}))++partialPaints;}
            peak=std::max(peak,resident());return !s.job.isActive();
        }});
        expect(activePaints>0&&partialPaints>0&&s.currentState.hasMeasurement&&!ui.audition.active(),"native picture and progressive timeline remain responsive during analysis");
        expect(s.source.analysis->identity==sourceIdentity&&s.currentState.analysis->identity==currentIdentity&&s.lastAppliedExplanation->text()==history,"reference attachment and analysis preserve exact audio and explanation identities");
        selectOutput(ui,*output,6);
        for(const auto* offset:{"-1001/30000","7/3"}){
            ui.videoPanel.offset.setText(offset,false);click(ui.videoPanel.apply);
            for(int64_t first:{int64_t(0),int64_t(48000),int64_t(240000),int64_t(0),s.sourceDocument.media.frames}){
                ui.audition.stop();const auto req=s.auditionRequest();
                expect(ui.audition.seek(first,req,s.sourceDocument.routing).isEmpty(),"native stopped exact seek with signed reference offset");settle(ui);
                if(first<s.sourceDocument.media.frames){
                    const auto error=ui.audition.play(req,s.sourceDocument.routing,false,true);expect(error.isEmpty(),"explicit injected audition starts");until("injected device did not start with reference",[&]{ui.audition.poll();return output->isPlaying();});
                    juce::AudioBuffer<float> block(6,256);
                    for(int tick=0;tick<100;++tick){@autoreleasepool {output->tick(block);ui.audition.poll();ui.videoPanel.synchronize();const auto t=Clock::now();auto a=ui.timeline.createComponentSnapshot(ui.timeline.getLocalBounds());auto b=ui.videoPanel.createComponentSnapshot(ui.videoPanel.getLocalBounds());longestPaint=std::max(longestPaint,since(t)*1000);pump();peak=std::max(peak,resident());}}
                    ui.audition.stop();settle(ui);expect(ui.audition.underruns()==0,"no audio discontinuity during native picture and timeline painting");
                }
            }
        }
        click(ui.audioPanel.closeButton);auto row=juce::var(new juce::DynamicObject());auto* r=row.getDynamicObject();r->setProperty("file",name);r->setProperty("domain",domain==1?"Source":"Current");r->setProperty("sha256",ui.videoPanel.view.asset->sha256);r->setProperty("analysis_paints",activePaints);r->setProperty("partial_waveform_paints",partialPaints);r->setProperty("longest_snapshot_ms",longestPaint);r->setProperty("skipped_frames",juce::int64(ui.videoPanel.view.droppedFrames));r->setProperty("late_polls",juce::int64(ui.videoPanel.view.latePolls));r->setProperty("underruns",juce::int64(ui.audition.underruns()));rows.add(row);std::cout<<name<<" domain="<<domain<<" rss="<<resident()<<" peak="<<peak<<std::endl;
    }
    ui.timerCallback();snapshot(ui,out.getChildFile("timeline-after-video.png"));picture(ui.videoPanel.createComponentSnapshot(ui.videoPanel.getLocalBounds()),out.getChildFile("reference-after-video.png"));
    put(report,"video_cases",rows);put(report,"peak_rss_bytes",juce::String(peak));put(report,"peak_rss_growth",juce::String(peak-initial));put(report,"physical_devices_opened",false);put(report,"clock","injected callbacks, driver-assisted native windows");expect(peak<=768ULL*1024*1024,"native video total RSS within existing 768 MiB video budget");
}
}
int main(int argc,char** argv){
    if(argc!=4)return 2;const juce::File out(argv[3]);if(out.exists()||out.createDirectory().failed())return 2;juce::ScopedJuceInitialiser_GUI gui;juce::var report(new juce::DynamicObject());bool passed=false;
    try{juce::DocumentWindow window("M18 native video qualification",juce::Colours::darkgrey,juce::DocumentWindow::allButtons);auto content=std::make_unique<Main>();auto* ui=content.get();window.setUsingNativeTitleBar(true);window.setContentOwned(content.release(),true);window.centreWithSize(1240,1000);window.setVisible(true);window.toFront(true);pump();qualifyVideo(*ui,juce::File(argv[1]),juce::File(argv[2]),out,report);passed=true;}
    catch(const std::exception& e){put(report,"failure",e.what());std::cerr<<e.what()<<std::endl;}
    put(report,"passed",passed);put(report,"checks",checks);put(report,"build",LEVELSCOPE_MEASUREMENT_BUILD);out.getChildFile("result.json").replaceWithText(juce::JSON::toString(report,true));return passed?0:1;
}

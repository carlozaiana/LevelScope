#pragma once
#include "Standalone/StandaloneVideo.h"
void beginVideoReaderAudit();
int videoReaderAuditCount(int);
namespace
{
int compositionVideoCases=0;
juce::Array<juce::var> compositionVideoEvidence;
juce::Array<juce::var> compositionExchangeEvidence;
void compositionVideo(const juce::File& root,const juce::File& h264,const juce::File& prores)
{
    auto f=fixtureComposition(root,48000,CompositionPreset::surroundMonoCentre,true,48000*4);auto s=compositionModel(f.recipe);measuredComposition(s,MeasurementDomain::source);s.setCurrentGainDb(-3.305074118944);measuredComposition(s,MeasurementDomain::currentState);
    beginVideoReaderAudit();
    for(auto domain:{MeasurementDomain::source,MeasurementDomain::currentState})for(const auto& file:{h264.getChildFile("cfr-exact.mov"),h264.getChildFile("vfr-exact.mov"),prores.getChildFile("422-with-tone.mov")})
    {
        s.selectAuditionDomain(domain);const auto request=s.auditionRequest();const auto identity=request.identity();auto decoder=video::makeDecoder();const auto asset=decoder->load(file,0,[]{return false;});expect(asset!=nullptr,"composed audio native H264/ProRes reference opens");video::Player player;player.attach(file,asset->selectedStream);uint64_t generation=1;
        auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));OutputRouting muted;audition.synchronize(request,muted);expect(audition.connect(48000).isEmpty(),"composed video injected device only");
        for(auto offset:{video::Time(-1001,30000),video::Time(7,3)})for(int64_t first:{int64_t(0),int64_t(72000),int64_t(144000),int64_t(0),int64_t(192000)})
        {
            s.referenceVideo={file,asset->sha256,asset->selectedStream,offset};expect(audition.seek(first,request,muted).isEmpty(),"composition exact seek preserves audio programme coordinates");
            if(first<192000)
            {
                expect(audition.play(request,muted,false,true).isEmpty(),"composed video all-muted injected play");awaitPlaying(audition,*injected);juce::AudioBuffer<float> block(maxChannels,256);injected->tick(block);std::this_thread::sleep_for(std::chrono::milliseconds(10));audition.poll();audition.stop();
                expect(audition.positionFrames()==first+256&&audition.underruns()==0,"composed playback presented position follows uninterrupted authoritative audio callback");
            }
            const auto target=video::Time(audition.positionFrames(),48000)-offset;const auto wanted=video::frameAt(asset->frames,target);player.invalidate();++generation;bool matched=false;const auto start=compositionNow(),deadline=start+10000;video::Player::View view;
            while(compositionNow()<deadline){view=player.update(target,generation);if(view.asset&&((wanted&&view.picture&&view.picture->index==*wanted)||(!wanted&&!view.picture))){matched=true;break;}std::this_thread::sleep_for(std::chrono::milliseconds(2));}
            expect(matched&&s.auditionRequest().identity()==identity&&s.source.hasMeasurement&&s.currentState.hasMeasurement,"composed Source/Current rational offset seek restart EOF picture selection and measurement identity");
            auto evidence=juce::var(new juce::DynamicObject());auto* e=evidence.getDynamicObject();e->setProperty("file",file.getFileName());e->setProperty("sha256",asset->sha256);e->setProperty("domain",domain==MeasurementDomain::source?"Source":"Current State");e->setProperty("audio_frame",juce::String(audition.positionFrames()));e->setProperty("video_target",juce::String(target.text()));e->setProperty("offset",juce::String(offset.text()));e->setProperty("selected_frame",view.picture?juce::var(static_cast<juce::int64>(view.picture->index)):juce::var());e->setProperty("selection_ms",compositionNow()-start);e->setProperty("player_skipped_frames",static_cast<juce::int64>(view.droppedFrames));e->setProperty("late_polls",static_cast<juce::int64>(view.latePolls));e->setProperty("physical_paint_test",false);compositionVideoEvidence.add(evidence);
        }
        audition.disconnect();player.clear();++compositionVideoCases;
    }
    expect(compositionVideoCases==6&&videoReaderAuditCount(0)>0&&videoReaderAuditCount(1)==0,"six composed reference cases; embedded video audio remains undecoded");
}
void compositionExchange(const juce::File& from)
{
    for(int p=0;p<3;++p)for(bool grouped:{false,true})
    {
        const auto folder=from.getChildFile("48000-"+juce::String(p)+(grouped?"-groups":"-files"));JobControl job;auto read=sf::read(folder.getChildFile("composition.levelscope"),job);
        expect(read.document&&read.document->source.composition&&read.document->sourceCache,"foreign composition retains original measured provenance");const auto confirmation=read.document->source.composition->confirmation;const auto original=read.document->sourceCache;
        for(const auto& member:read.document->source.composition->files())expect(sf::relinkAudioMember(*read.document,member.id,folder.getChildFile(member.file.getFileName())),"explicit compound slot relink after cross-Mac exchange");
        auto verified=sf::verify(*read.document,job);expect(verified.prepared&&verified.document->audioResolved,"all exchanged component members verified");StandaloneSessionModel s;juce::String notice;expect(sf::adopt(*verified.document,*verified.prepared,s,notice),"safe composed exchange adoption");
        const auto context=s.measurementContext(MeasurementDomain::source);const bool compatible=original->context.buildVersion==context.buildVersion&&original->context.platformAndFloatingPointMode==context.platformAndFloatingPointMode;
        expect(original->context.buildVersion==context.buildVersion&&original->context.platformAndFloatingPointMode!=context.platformAndFloatingPointMode,"exchange requires same frozen build from the other actual platform");
        expect(s.source.hasMeasurement==compatible&&s.currentState.hasMeasurement==compatible&&s.sourceDocument.composition->confirmation==confirmation,"exchange keeps exact confirmation and conservative build/platform cache policy");
        for(int route:s.sourceDocument.routing.fileToDevice)expect(route==-1,"cross-platform reopen routes muted");
        measuredComposition(s,MeasurementDomain::source);measuredComposition(s,MeasurementDomain::currentState);expect(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"exchanged composed Current remeasured and evaluated on receiving Mac");
        auto evidence=juce::var(new juce::DynamicObject());auto* e=evidence.getDynamicObject();e->setProperty("preset",p);e->setProperty("representation",grouped?"grouped":"interleaved");e->setProperty("confirmation",confirmation);e->setProperty("origin_platform",original->context.platformAndFloatingPointMode);e->setProperty("receiving_platform",context.platformAndFloatingPointMode);e->setProperty("foreign_cache_reused",compatible);e->setProperty("remeasured_current_pass",true);compositionExchangeEvidence.add(evidence);
    }
}
}

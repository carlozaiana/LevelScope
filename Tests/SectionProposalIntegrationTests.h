#pragma once
#include "Standalone/StandaloneVideo.h"
void beginVideoReaderAudit();
int videoReaderAuditCount(int);
namespace
{
juce::Array<juce::var> videoEvidence,exchangeEvidence;
StandaloneSessionModel openApplied(const juce::File& file)
{
    JobControl job;auto read=sf::read(file,job);expect(read.document.has_value(),"read applied M13 fixture session");auto verified=sf::verify(*read.document,job);
    expect(verified.prepared&&verified.document->audioResolved,"verify applied M13 fixture members");StandaloneSessionModel s;juce::String notice;
    expect(sf::adopt(*verified.document,*verified.prepared,s,notice)&&s.currentStateDocument.sectionGains,"adopt exact saved section/global plan");return s;
}
void proposalVideo(const juce::File& session,const juce::File& fixtures)
{
    auto s=openApplied(session);expect(s.source.hasMeasurement&&s.currentState.hasMeasurement,"same-platform measured fixture caches retained");
    beginVideoReaderAudit();
    for(auto domain:{MeasurementDomain::source,MeasurementDomain::currentState})for(const auto* name:{"h264-cfr.mov","h264-vfr.mov","422-25-1.mov","lt-30000-1001.mov"})
    {
        const auto file=fixtures.getChildFile(name);s.selectAuditionDomain(domain);const auto request=s.auditionRequest();const auto identity=request.identity();
        auto decoder=video::makeDecoder();const auto asset=decoder->load(file,0,[]{return false;});expect(asset!=nullptr,"M13 native H264/ProRes reference opens");
        video::Player player;player.attach(file,asset->selectedStream);uint64_t generation=1;
        auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));OutputRouting muted;audition.synchronize(request,muted);expect(audition.connect(48000).isEmpty(),"section/global video injected device only");
        const int64_t frames=request.context.media.frames;
        for(auto offset:{video::Time(-1001,30000),video::Time(7,3)})for(int64_t first:{int64_t(0),int64_t(240000-128),int64_t(288000),int64_t(0),frames})
        {
            s.referenceVideo={file,asset->sha256,asset->selectedStream,offset};expect(audition.seek(first,request,muted).isEmpty(),"exact seek across applied transition preserves programme coordinates");
            if(first<frames)
            {
                expect(audition.play(request,muted,false,true).isEmpty(),"all-muted section/global injected play");awaitPlaying(audition,*injected);juce::AudioBuffer<float> block(maxChannels,256);injected->tick(block);std::this_thread::sleep_for(std::chrono::milliseconds(10));audition.poll();audition.stop();
                expect(audition.positionFrames()==first+256&&audition.underruns()==0,"presented position follows uninterrupted authoritative audio");
            }
            const auto target=video::Time(audition.positionFrames(),48000)-offset;const auto wanted=video::frameAt(asset->frames,target);player.invalidate();++generation;bool matched=false;
            const auto start=compositionNow();video::Player::View view;
            while(compositionNow()-start<10000){view=player.update(target,generation);if(view.asset&&((wanted&&view.picture&&view.picture->index==*wanted)||(!wanted&&!view.picture))){matched=true;break;}std::this_thread::sleep_for(std::chrono::milliseconds(2));}
            expect(matched&&s.auditionRequest().identity()==identity&&s.source.hasMeasurement&&s.currentState.hasMeasurement,"Source/Current rational offsets resume restart EOF and exact picture selection preserve audio identity");
            auto evidence=juce::var(new juce::DynamicObject());auto* e=evidence.getDynamicObject();e->setProperty("file",name);e->setProperty("sha256",asset->sha256);e->setProperty("domain",domain==MeasurementDomain::source?"Source":"Current State");
            e->setProperty("audio_frame",juce::String(audition.positionFrames()));e->setProperty("target",juce::String(target.text()));e->setProperty("offset",juce::String(offset.text()));e->setProperty("selection_ms",compositionNow()-start);e->setProperty("skipped_frames",juce::int64(view.droppedFrames));e->setProperty("late_polls",juce::int64(view.latePolls));e->setProperty("physical_paint_test",false);videoEvidence.add(evidence);
        }
        audition.disconnect();player.clear();
    }
    expect(videoEvidence.size()==80&&videoReaderAuditCount(0)>0&&videoReaderAuditCount(1)==0,"eight composed section/global video cases; embedded audio remains undecoded");
}
void proposalExchange(const juce::File& from)
{
    for(int channels:{2,6})for(int mode:{0,1,3})
    {
        const auto folder=from.getChildFile("48000-"+juce::String(channels)+"-"+juce::String(mode));JobControl job;auto read=sf::read(folder.getChildFile("applied.levelscope"),job);
        expect(read.document&&read.document->current.sectionGains&&read.document->currentCache,"foreign session has exact applied plan and provenance");
        const auto original=read.document->currentCache;const auto plan=read.document->current.sectionGains->audioKey();const auto gain=read.document->current.staticGainDb;
        auto members=read.document->source.composition?read.document->source.composition->files():read.document->source.assembly?read.document->source.assembly->members:std::vector<ProgrammeMember>{{"audio",read.document->source.file,read.document->source.media,read.document->source.layout.roles[0],0}};
        for(const auto& member:members)expect(read.document->source.composition||read.document->source.assembly
            ?sf::relinkAudioMember(*read.document,member.id,folder.getChildFile(member.file.getFileName()))
            :sf::relinkAudio(*read.document,folder.getChildFile(member.file.getFileName())),"foreign exact member relink");
        auto verified=sf::verify(*read.document,job);expect(verified.prepared&&verified.document->audioResolved,"all foreign members verified atomically");StandaloneSessionModel s;juce::String notice;expect(sf::adopt(*verified.document,*verified.prepared,s,notice),"safe cross-Mac applied session adoption");
        const auto context=s.measurementContext(MeasurementDomain::currentState);
        expect(original->context.buildVersion==context.buildVersion&&original->context.platformAndFloatingPointMode!=context.platformAndFloatingPointMode,"exchange from other actual platform with same frozen build");
        expect(!s.currentState.hasMeasurement&&!s.source.hasMeasurement&&!s.sectionProposal&&s.currentStateDocument.sectionGains->audioKey()==plan&&s.currentStateDocument.staticGainDb==gain,"foreign caches discarded with exact applied plan retained and no private authority");
        for(int route:s.sourceDocument.routing.fileToDevice)expect(route==-1,"exchanged session routes muted");
        measuredComposition(s,MeasurementDomain::source);measuredComposition(s,MeasurementDomain::currentState);expect(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"exchanged applied candidate independently remeasured and passes");
        auto evidence=juce::var(new juce::DynamicObject());auto* e=evidence.getDynamicObject();e->setProperty("case",folder.getFileName());e->setProperty("origin_platform",original->context.platformAndFloatingPointMode);e->setProperty("receiver",context.platformAndFloatingPointMode);e->setProperty("gain",gain);e->setProperty("cache_reused",false);exchangeEvidence.add(evidence);
    }
}
}

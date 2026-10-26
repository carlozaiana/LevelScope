#pragma once
#include "Standalone/StandaloneVideo.h"
void beginVideoReaderAudit();
int videoReaderAuditCount(int);
namespace
{
int exportVideoCases=0;
void exportVideo(const juce::File& audio,const juce::File& videoFolder,const juce::File& proresFolder)
{
    JobControl control;const auto validation=validateMedia(audio,control);expect(validation.succeeded(),"published PCM reimport for reference playback");
    StandaloneSessionModel s;expect(s.setValidatedSource(audio,validation.media)&&s.confirmLayout()&&s.initializeCurrentStateFromSource(),"exported PCM has explicit reimported roles");
    s.setCurrentGainDb(-3.305074118944);std::vector<juce::File> files{videoFolder.getChildFile("cfr-exact.mov"),videoFolder.getChildFile("vfr-exact.mov")};for(const char* name:{"proxy-with-tone.mov","lt-with-tone.mov","422-with-tone.mov"})files.push_back(proresFolder.getChildFile(name));
    beginVideoReaderAudit();
    for(auto domain:{MeasurementDomain::source,MeasurementDomain::currentState})for(const auto& file:files)
    {
        s.selectAuditionDomain(domain);const auto identity=s.auditionRequest().identity();auto decoder=video::makeDecoder();const auto asset=decoder->load(file,0,[]{return false;});expect(asset!=nullptr,"export native reference opens");
        video::Player player;player.attach(file,asset->selectedStream);uint64_t generation=1;
        for(auto offset:{video::Time(-1001,30000),video::Time(7,3)})for(auto time:{video::Time(0,1),video::Time(3,2),video::Time(3,1),video::Time(0,1),video::Time(validation.media.frames,48000)})
        {
            s.referenceVideo={file,asset->sha256,asset->selectedStream,offset};const auto target=time-offset;const auto wanted=video::frameAt(asset->frames,target);player.invalidate();++generation;bool matched=false;const auto deadline=exportNow()+10000;
            while(exportNow()<deadline){const auto view=player.update(target,generation);if(view.asset&&((wanted&&view.picture&&view.picture->index==*wanted)||(!wanted&&!view.picture))){matched=true;break;}std::this_thread::sleep_for(std::chrono::milliseconds(2));}
            expect(matched&&s.auditionRequest().identity()==identity,"export reimport Source/Current rational seek/offset/restart/EOF follows audio coordinates");
        }
        player.clear();++exportVideoCases;
    }
    expect(exportVideoCases==10&&videoReaderAuditCount(0)>0&&videoReaderAuditCount(1)==0,"ten export reference cases; embedded audio remains undecoded");
}
void prepareExportClock(const juce::File& fixtureFolder,const juce::File& destination)
{
    // Explicit stereo test programme from the two known-aligned generated clock files.
    // This is fixture construction, not implicit channel dropping in import/playback.
    JobControl control;auto inspected=inspectProgrammeFiles({fixtureFolder.getChildFile("mono-0.flac"),fixtureFolder.getChildFile("mono-1.flac")},control);
    expect(inspected.error==ErrorCode::none,"stereo qualification clock files inspected");inspected.group.members[0].role=juce::AudioChannelSet::left;inspected.group.members[1].role=juce::AudioChannelSet::right;
    expect(inspected.group.confirmAlignment(),"known generator alignment explicitly confirmed");StandaloneSessionModel s;
    expect(s.setValidatedProgramme(inspected.group)&&s.initializeCurrentStateFromSource(),"stereo export clock session");s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
    measure(s,MeasurementDomain::source);s.setCurrentGainDb(-23-*s.source.integratedLufs);measure(s,MeasurementDomain::currentState);
    const auto request=ex::capture(s);auto rendered=ex::render(request,destination.getParentDirectory(),control);expect(rendered.value.has_value(),"clock PCM rendered and independently remeasured");
    expect(ex::publish(**rendered.value,request,destination,ex::approvalIdentity(**rendered.value,request),control).value.has_value(),"clock actual PCM passes and is published for reimported video tests");
}
}

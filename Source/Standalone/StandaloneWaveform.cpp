#include "StandaloneWaveform.h"
#include "StandaloneMeasurementService.h"

namespace levelscope::standalone
{
WaveformProgress::WaveformProgress(std::shared_ptr<waveform::Budget> b):budget(std::move(b)){}
std::shared_ptr<waveform::Overview> WaveformProgress::begin(waveform::Spec spec) noexcept
{
    try
    {
        if(state()!=State::waiting)return {};
        auto next=waveform::Overview::create(spec,budget);
        std::lock_guard<std::mutex> lock(mutex);
        if(state()==State::abandoned)return {};
        data=next;phase.store(next?State::collecting:State::unavailable);return next;
    }
    catch(...) {unavailable();return {};}
}
std::shared_ptr<const waveform::Overview> WaveformProgress::overview() const
{std::lock_guard<std::mutex> lock(mutex);return data;}
void WaveformProgress::finish(const std::shared_ptr<waveform::Overview>& value) noexcept
{
    if(!value || !value->finish()){unavailable();return;}
    State expected=State::collecting;phase.compare_exchange_strong(expected,State::complete);
}
void WaveformProgress::unavailable() noexcept
{std::lock_guard<std::mutex> lock(mutex);data.reset();if(state()!=State::abandoned)phase.store(State::unavailable);}
void WaveformProgress::abandon() noexcept
{std::lock_guard<std::mutex> lock(mutex);phase.store(State::abandoned);data.reset();}
WaveformBinding WaveformBinding::from(const ServiceContext& c,MeasurementDomain d)
{return {c.sessionId,OfflineAnalysisService::identityFor(c,d,c.media.sha256),c.ticket,d};}
bool WaveformBinding::sameAudio(const WaveformBinding& b) const
{return !runtime.isEmpty() && runtime==b.runtime && identity==b.identity && domain==b.domain;}
bool WaveformBinding::sameJob(const WaveformBinding& b) const
{return sameAudio(b)&&ticket.id==b.ticket.id&&ticket.sessionRevision==b.ticket.sessionRevision;}
std::shared_ptr<WaveformProgress> WaveformCache::begin(const WaveformBinding& binding)
{
    if(binding.domain!=MeasurementDomain::source&&binding.domain!=MeasurementDomain::currentState)return {};
    if(progress)progress->abandon();
    slots[binding.domain==MeasurementDomain::source?0:1]={};active=binding;
    progress=std::make_shared<WaveformProgress>(budget);return progress;
}
bool WaveformCache::accept(const WaveformBinding& completed,const WaveformBinding& expected)
{
    if(!progress || !active.sameJob(completed) || !completed.sameJob(expected)
        || progress->state()!=WaveformProgress::State::complete)return false;
    auto data=progress->overview();if(!data||!data->complete())return false;
    slots[completed.domain==MeasurementDomain::source?0:1]={completed,std::move(data)};progress.reset();return true;
}
void WaveformCache::reject(const WaveformBinding& completed)
{if(progress&&active.sameJob(completed)){progress->abandon();progress.reset();}}
void WaveformCache::synchronize(const WaveformBinding& source,const WaveformBinding& current,bool valid,bool initialized)
{
    if(!valid||!slots[0].binding.sameAudio(source))slots[0]={};
    if(!valid||!initialized||!slots[1].binding.sameAudio(current))slots[1]={};
    if(progress&&(!valid||!active.sameJob(active.domain==MeasurementDomain::source?source:current)
        ||(active.domain==MeasurementDomain::currentState&&!initialized)))
    {progress->abandon();progress.reset();}
}
bool WaveformCache::shareInitialCurrent(const AnalysisSnapshot& source,const AnalysisSnapshot& current)
{
    const auto originRuntime=source.originSessionId.isEmpty()?source.context.sessionId:source.originSessionId;
    const auto originJob=source.originSessionId.isEmpty()?source.context.ticket.id:source.originJobId;
    if(source.domain!=MeasurementDomain::source||current.domain!=MeasurementDomain::currentState
        || current.identicalSourceMeasurementIdentity!=source.identity
        || current.originSessionId!=originRuntime || current.originJobId!=originJob
        || current.context.sessionId!=source.context.sessionId || current.detail.frames!=source.detail.frames
        || current.detail.sampleRate!=source.detail.sampleRate || current.context.fileLayout.count!=source.context.fileLayout.count
        || OfflineAnalysisService::validateRequest(current.context,MeasurementDomain::source).isNotEmpty()
        || OfflineAnalysisService::identityFor(current.context,MeasurementDomain::source,current.measuredContentHash)!=source.identity
        || current.identity!=OfflineAnalysisService::identityFor(current.context,current.domain,current.measuredContentHash)
        || !slots[0].data || slots[0].binding.identity!=source.identity || slots[0].binding.runtime!=source.context.sessionId)return false;
    slots[1]={WaveformBinding::from(current.context,MeasurementDomain::currentState),slots[0].data};return true;
}
WaveformView WaveformCache::view(const WaveformBinding& binding) const
{
    if(progress&&active.sameJob(binding))return {progress->overview(),false,progress->state()};
    if(binding.domain!=MeasurementDomain::source&&binding.domain!=MeasurementDomain::currentState)return {};
    const auto& slot=slots[binding.domain==MeasurementDomain::source?0:1];
    if(slot.data&&slot.binding.sameAudio(binding))return {slot.data,true,WaveformProgress::State::complete};
    return {};
}
}

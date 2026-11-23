#include "StandaloneAudition.h"
#include <cmath>
#include <chrono>
#include <stdexcept>
#if JUCE_MAC
#include <mach/mach_time.h>
#include <CoreAudio/CoreAudio.h>
#endif

namespace levelscope::standalone
{
namespace
{
#if JUCE_MAC
// Read-only lookup. HAL output timestamps already include scheduling/buffering;
// only the device + stream hardware delay is added to that timestamp.
int hardwareDelay (const juce::AudioIODevice& device)
{
    if (device.getTypeName() != "CoreAudio") return -1;
    AudioObjectPropertyAddress address { kAudioHardwarePropertyDevices, kAudioObjectPropertyScopeGlobal, 0 };
    UInt32 size = 0;
    if (AudioObjectGetPropertyDataSize (kAudioObjectSystemObject, &address, 0, nullptr, &size)) return -1;
    std::vector<AudioDeviceID> ids (size / sizeof (AudioDeviceID));
    if (AudioObjectGetPropertyData (kAudioObjectSystemObject, &address, 0, nullptr, &size, ids.data())) return -1;
    int result = -1;
    for (auto id : ids)
    {
        address = { kAudioObjectPropertyName, kAudioObjectPropertyScopeGlobal, 0 };
        CFStringRef name = nullptr; size = sizeof name;
        if (AudioObjectGetPropertyData (id, &address, 0, nullptr, &size, &name) || name == nullptr) continue;
        const auto matches = juce::String::fromCFString (name) == device.getName(); CFRelease (name);
        if (!matches) continue;
        address = { kAudioDevicePropertyStreams, kAudioDevicePropertyScopeOutput, 0 };
        if (AudioObjectGetPropertyDataSize (id, &address, 0, nullptr, &size) || size == 0) continue;
        std::vector<AudioStreamID> streams (size / sizeof (AudioStreamID));
        if (AudioObjectGetPropertyData (id, &address, 0, nullptr, &size, streams.data())) return -1;
        address = { kAudioDevicePropertyLatency, kAudioDevicePropertyScopeOutput, 0 };
        UInt32 delay = 0; size = sizeof delay;
        if (AudioObjectGetPropertyData (id, &address, 0, nullptr, &size, &delay)) return -1;
        UInt32 longest = 0;
        for (auto stream : streams)
        {
            UInt32 streamDelay = 0; size = sizeof streamDelay;
            address = { kAudioStreamPropertyLatency, kAudioObjectPropertyScopeGlobal, 0 };
            if (AudioObjectGetPropertyData (stream, &address, 0, nullptr, &size, &streamDelay)) return -1;
            longest = std::max (longest, streamDelay);
        }
        if (result != -1 || static_cast<uint64_t> (delay) + longest > INT_MAX) return -1;
        result = static_cast<int> (delay + longest);
    }
    return result;
}
#endif
}

uint64_t StandaloneAudition::hostTimeNs() noexcept
{
#if JUCE_MAC
    static const auto scale = [] { mach_timebase_info_data_t s; mach_timebase_info (&s); return s; }();
    return static_cast<uint64_t> ((__uint128_t) mach_absolute_time() * scale.numer / scale.denom);
#else
    return static_cast<uint64_t> (std::chrono::duration_cast<std::chrono::nanoseconds>
        (std::chrono::steady_clock::now().time_since_epoch()).count());
#endif
}

int64_t StandaloneAudition::positionFrames() const
{
    if (state != AuditionState::playing || sourceRate <= 0) return heldFrame;
    return std::min (sourceFrames, presentationClock.position (hostTimeNs(), sourceRate, generation, startFrame));
}

void StandaloneAudition::sourceChanged (const MediaInfo& media)
{
    stop();
    overloadWarning.clear();
    preparedSource.reset(); preparedRuntime.clear(); preparedAssets.reset(); expectedAssets.reset(); expectedIdentity.clear();
    fault = ErrorCode::none;
    sourceRate = static_cast<int64_t> (media.sampleRate);
    sourceFrames = media.frames; sourceHash = media.sha256;
    heldFrame = startFrame = 0;
}
// [BEGIN LS-STANDALONE-AUDITION-IMPLEMENTATION]
void AuditionBuffer::reset (int channels, const OutputRouting& map)
{
    channelCount = juce::jlimit (0, maxChannels, channels);
    routing = map;
    read.store (0); written.store (0); underruns.store (0);
    endOfFile.store (false);
    lastGain = 0; // Every explicit restart uses a fresh ramp, including overload recovery.
    monitorOverload.store (false); overloadChannel.store (-1); overloadFrame.store (0);
}

int AuditionBuffer::writable() const noexcept
{
    return capacity - static_cast<int> (written.load (std::memory_order_relaxed) - read.load (std::memory_order_acquire));
}

void AuditionBuffer::push (const float* const* input, int frames) noexcept
{
    const auto start = written.load (std::memory_order_relaxed);
    for (int channel = 0; channel < channelCount; ++channel)
        for (int frame = 0; frame < frames; ++frame)
            samples[channel][(start + static_cast<uint64_t> (frame)) % capacity] = input[channel][frame];
    written.store (start + static_cast<uint64_t> (frames), std::memory_order_release);
}

void AuditionBuffer::render (float* const* output, int outputCount, int frames) noexcept
{
    for (int channel = 0; channel < outputCount; ++channel)
        if (output[channel] != nullptr) juce::FloatVectorOperations::clear (output[channel], frames);
    if (monitorOverload.load()) return;
    const auto start = read.load (std::memory_order_relaxed);
    const auto available = written.load (std::memory_order_acquire) - start;
    const auto count = static_cast<int> (juce::jmin<uint64_t> (static_cast<uint64_t> (frames), available));
    if (count < frames && ! endOfFile.load()) underruns.fetch_add (1, std::memory_order_relaxed);
    const auto target = gain.load (std::memory_order_relaxed);
    // Validate the ENTIRE routed callback first, without consuming or changing
    // ramp/clock state. The already-cleared outputs stay entirely silent on fault.
    float ramp = lastGain;
    for (int frame = 0; frame < count; ++frame)
    {
        ramp += (target - ramp) * 0.002f;
        for (int channel = 0; channel < channelCount; ++channel)
        {
            const auto destination = routing.fileToDevice[channel];
            if (destination < 0 || destination >= outputCount || output[destination] == nullptr) continue;
            const float value = samples[channel][(start + static_cast<uint64_t> (frame)) % capacity] * ramp;
            if (! std::isfinite (value) || std::abs (value) > 1.0f)
            {
                overloadChannel.store (channel); overloadFrame.store (start + static_cast<uint64_t> (frame));
                monitorOverload.store (true); return;
            }
        }
    }
    if (observer && count > 0)
    {
        std::array<const float*, maxChannels> pointers {};
        int offset = 0;
        while (offset < count)
        {
            const auto index = static_cast<int> ((start + static_cast<uint64_t> (offset)) % capacity);
            const int n = std::min (count - offset, capacity - index);
            for (int ch = 0; ch < channelCount; ++ch) pointers[ch] = samples[ch].data() + index;
            observer->observe (pointers.data(), channelCount, n, captureStart + static_cast<int64_t> (start) + offset);
            offset += n;
        }
    }
    for (int frame = 0; frame < count; ++frame)
    {
        lastGain += (target - lastGain) * 0.002f;
        for (int channel = 0; channel < channelCount; ++channel)
        {
            const auto destination = routing.fileToDevice[channel];
            if (destination >= 0 && destination < outputCount && output[destination] != nullptr)
                output[destination][frame] = samples[channel][(start + static_cast<uint64_t> (frame)) % capacity] * lastGain;
        }
    }
    read.store (start + static_cast<uint64_t> (count), std::memory_order_release);
}

StandaloneAudition::~StandaloneAudition()
{
    disconnect();
    if (preparation.valid()) { preparation.wait(); preparation.get(); }
}

const juce::OwnedArray<juce::AudioIODeviceType>& StandaloneAudition::deviceTypes()
{
    if (! catalog) catalog = std::make_unique<juce::AudioDeviceManager>();
    return catalog->getAvailableDeviceTypes(); // Never initialise/open through the manager.
}

juce::StringArray StandaloneAudition::scanOutputs (int index)
{
    const auto& types = deviceTypes();
    if (! juce::isPositiveAndBelow (index, types.size())) return {};
    types[index]->scanForDevices();
    return types[index]->getDeviceNames (false);
}

juce::String StandaloneAudition::selectDevice (int index, const juce::String& name)
{
    disconnect();
    device.reset();
    const auto& types = deviceTypes();
    if (! juce::isPositiveAndBelow (index, types.size()) || name.isEmpty()) return "Select an output device";
    device.reset (types[index]->createDevice (name, {}));
    message = device ? "Selected " + name + "; device closed" : "Device is unavailable";
    return device ? juce::String() : message;
}

juce::StringArray StandaloneAudition::outputConfiguration() const
{ return device ? device->getOutputChannelNames() : juce::StringArray{}; }

juce::StringArray StandaloneAudition::outputNames() const
{
    if (! device) return {};
    auto names = device->getOutputChannelNames();
    while (names.size() > maxChannels) names.remove (names.size() - 1);
    return names;
}

void StandaloneAudition::selectOutputDevice(std::unique_ptr<juce::AudioIODevice> selected)
{
    disconnect();
    device = std::move(selected);
    setMonitorGainDb(-30);
}

juce::String StandaloneAudition::connect (double sourceRate)
{
    if (preparation.valid() && playAfterPreparation) return "Cancel the pending Play before opening a device";
    // An explicit device connection is independent of device-free warmup.
    if (!preparation.valid() && !readyResources) stop();
    if (! device || ! std::isfinite (sourceRate) || sourceRate <= 0) return "Select a device and validated source";
    device->close();
    juce::BigInteger outputs;
    outputs.setRange (0, outputNames().size(), true);
    if (outputs.isZero()) return "Device has no available outputs";
    if (! device->getAvailableSampleRates().contains (sourceRate))
        return "Incompatible programme rate: device does not offer " + juce::String(sourceRate, 0) + " Hz; resampling is unavailable";
    auto error = device->open ({}, outputs, sourceRate, device->getDefaultBufferSize());
    if (error.isNotEmpty() || ! device->isOpen() || std::abs (device->getCurrentSampleRate() - sourceRate) > 0.01
        || device->getActiveOutputChannels() != outputs)
    {
        device->close();
        state = AuditionState::failed;
        message = error.isNotEmpty() ? "Output is busy or unavailable: " + error : "Device did not accept the requested sample rate and outputs";
        return message;
    }
    rate = sourceRate;
    latencyFrames = std::max (0, device->getOutputLatencyInSamples());
#if JUCE_MAC
    hardwareLatencyFrames = hardwareDelay (*device);
#endif
    (void) hostTimeNs(); // initialize timebase on the message thread, never first-use in RT
    expectedOutputs = outputs.countNumberOfSetBits();
    state = preparation.valid() ? (preparationControl->isCancelled() ? AuditionState::cancelling : AuditionState::preparing) : AuditionState::ready;
    message = "Device open at " + juce::String (rate, 0) + " Hz; stopped";
    return {};
}

void StandaloneAudition::stopWorker()
{
    // Remove/quiesce callback BEFORE changing any callback-owned buffer/configuration.
    if (device) device->stop();
    cancel.store (true);
    if (decoder.joinable()) decoder.join();
}

void StandaloneAudition::stop()
{
    heldFrame = positionFrames();
    if (preparationControl) preparationControl->cancel.store (true);
    playAfterPreparation = false;
    resumeLatestSeek = false;
    readyResources.reset();
    stopWorker();
    ++generation;
    presentationClock.reset();
    state = preparation.valid() ? AuditionState::cancelling : connected() ? AuditionState::ready : AuditionState::disconnected;
    message = preparation.valid() ? "Cancelling source preparation; stopped" : connected() ? "Stopped; Play resumes, Restart begins at zero" : "Device closed";
}

void StandaloneAudition::disconnect()
{
    stop();
    if (device) device->close();
    rate = 0;
    state = preparation.valid() ? AuditionState::cancelling : AuditionState::disconnected;
    message = preparation.valid() ? "Cancelling preparation; device closed" : "Device closed";
}

bool StandaloneAudition::playbackReady() const
{
    return readyResources && readyResources->identity == expectedIdentity
        && readyResources->first == heldFrame && readyResources->source->matches(pendingRequest.context.media)
        && pendingRequest.validationError().isEmpty() && pendingRequest.context.sourceFilesMatch();
}

juce::String StandaloneAudition::preparationStage() const
{
    if(state == AuditionState::cancelling) return "Cancelling playback preparation";
    return decoderPreparation && decoderPreparation->load() ? "Preparing playback decoder and buffer" : "Preparing verified audio for playback";
}

void StandaloneAudition::synchronize (const AuditionRequest& request, const OutputRouting& routing)
{
    const auto key = request.identity();
    const bool changed = key != expectedIdentity || !sameVerifiedAssetOwner(request.context.verifiedAssets,expectedAssets);
    const bool invalidSource = key.isNotEmpty() && !request.context.sourceFilesMatch();
    if(changed)overloadWarning.clear();
    if (((active() || preparing() || readyResources) && (changed || invalidSource))
        || (active() && (((state == AuditionState::playing || state == AuditionState::buffering) && routing.fileToDevice != playingRouting.fileToDevice)
            || (state == AuditionState::preparing && playAfterPreparation && routing.fileToDevice != pendingRouting.fileToDevice))))
    {
        stop(); message = "Audition stopped: source, roles, processing or routing changed";
    }
    expectedIdentity = key;
    expectedAssets = request.context.verifiedAssets;
    if (preparedSource && (preparedRuntime != request.context.sessionId || !preparedSource->matches(request.context.media)
        || !sameVerifiedAssetOwner(preparedAssets,request.context.verifiedAssets)
        || (request.context.verifiedAssets && request.context.verifiedAssets->expired())))
    { readyResources.reset(); preparedSource.reset(); preparedRuntime.clear(); preparedAssets.reset(); }
    if(invalidSource) { fault=ErrorCode::sourceChanged;message="Source changed; import again. Playback stopped."; }
}

juce::String StandaloneAudition::startPreparation (const AuditionRequest& request)
{
    const auto error = request.validationError();
    if (error.isNotEmpty()) return error;
    if (request.identity() != expectedIdentity || !sameVerifiedAssetOwner(request.context.verifiedAssets,expectedAssets))
        return "Stale audition request; select the current state again";
    if (! request.context.sourceFilesMatch())
    { fault = ErrorCode::sourceChanged; return "Source changed; import again"; }
    if (preparation.valid())
        return state==AuditionState::preparing && !preparationControl->isCancelled()
            && pendingRequest.identity()==request.identity() ? juce::String() : "Wait for playback preparation cancellation";
    fault = ErrorCode::none;
    pendingRequest = request;
    if (playbackReady()) return {};
    readyResources.reset();
    const auto& media=request.context.media;
    if(sourceHash != media.sha256)heldFrame=0;
    sourceHash=media.sha256;sourceRate=static_cast<int64_t>(media.sampleRate);sourceFrames=media.frames;
    startFrame=heldFrame;
    preparationControl = std::make_shared<JobControl>();
    preparationControl->feedback = std::make_shared<JobFeedback>();
    decoderPreparation = std::make_shared<std::atomic<bool>>(false);
    preparationFailure.clear();
    preparationGeneration = generation;
    auto retained=preparedRuntime==request.context.sessionId && sameVerifiedAssetOwner(preparedAssets,request.context.verifiedAssets)
        && preparedSource && preparedSource->matches(media) ? preparedSource : nullptr;
    preparation = std::async (std::launch::async, [request, begin=heldFrame, retained, storage=snapshotStorage, control=preparationControl, warming=decoderPreparation]
        () -> ServiceResult<std::unique_ptr<PlaybackResources>>
    {
        try
        {
            if(control->isCancelled())return {{},ErrorCode::cancelled,"Playback preparation cancelled"};
            auto verified=retained ? VerifiedSourceResult{retained,ErrorCode::none,{}} : prepareVerifiedSource(request.context,*control,storage);
            if(!verified.value)return {{},verified.error,verified.message};
            if(control->isCancelled())return {{},ErrorCode::cancelled,"Playback preparation cancelled"};
            warming->store(true);
            auto ready=std::make_unique<PlaybackResources>();ready->source=*verified.value;
            ready->identity=request.identity();ready->first=begin;ready->reader=ready->source->createReader();
            const auto& media=request.context.media;
            if(!ready->reader || ready->reader->sampleRate!=media.sampleRate || ready->reader->lengthInSamples!=media.frames
                || static_cast<int>(ready->reader->numChannels)!=media.channels)
                return {{},ErrorCode::decodeFailed,"Verified source cannot be reopened consistently"};
            ready->processing=std::make_unique<ProgrammeGainProcessor>(request.context.processing);
            constexpr int blockSize=2048;
            juce::AudioBuffer<float> block(media.channels,blockSize);
            // Bound FLAC tail seek pre-roll; discarded frames never enter the
            // processor, FIFO, observers or presentation clock.
            if(ready->reader->hasFlac && begin>0)
                for(int64_t first=std::max<int64_t>(0,begin-65536);first<begin;)
                {
                    if(control->isCancelled())return {{},ErrorCode::cancelled,"Playback preparation cancelled"};
                    const int count=int(std::min<int64_t>(blockSize,begin-first));
                    if(!ready->reader->read(block.getArrayOfWritePointers(),media.channels,first,count))
                        return {{},ErrorCode::decodeFailed,"Verified source seek pre-roll failed"};
                    first+=count;
                }
            const int frames=int(std::min<int64_t>(4096,std::max<int64_t>(0,media.frames-begin)));
            ready->prefix.setSize(media.channels,frames);
            for(int first=0;first<frames;first+=blockSize)
            {
                if(control->isCancelled())return {{},ErrorCode::cancelled,"Playback preparation cancelled"};
                const int count=std::min(blockSize,frames-first);
                if(!ready->reader->read(block.getArrayOfWritePointers(),media.channels,begin+first,count)
                    || !ready->processing->process(block.getArrayOfWritePointers(),media.channels,count,begin+first))
                    return {{},ErrorCode::decodeFailed,"Playback prefix could not be decoded and processed"};
                for(int c=0;c<media.channels;++c)ready->prefix.copyFrom(c,first,block,c,0,count);
            }
            if(control->isCancelled())return {{},ErrorCode::cancelled,"Playback preparation cancelled"};
            if(!request.context.sourceFilesMatch() || !ready->source->matches(media))
                return {{},ErrorCode::sourceChanged,"Source changed during playback preparation; import again"};
            return {std::move(ready),ErrorCode::none,"Playback prepared; stopped"};
        }
        catch(const std::exception& e){return {{},ErrorCode::decodeFailed,e.what()};}
    });
    state = AuditionState::preparing;
    message = "Preparing playback in background; stopped";
    return {};
}

juce::String StandaloneAudition::prepare (const AuditionRequest& request)
{
    if (active()) return "Stop before preparing a different audition";
    return startPreparation(request);
}

juce::String StandaloneAudition::play (const juce::File& file, const MediaInfo& media,
                                     const ChannelLayout& layout, const OutputRouting& routing, bool restart, bool allowAllMuted)
{
    const auto request = AuditionRequest::source (file, media, layout);
    synchronize (request, routing);
    return play (request, routing, restart, allowAllMuted);
}

juce::String StandaloneAudition::play (const AuditionRequest& request, const OutputRouting& routing, bool restart, bool allowAllMuted)
{
    const auto error = request.validationError();
    if (error.isNotEmpty()) return error;
    if (request.identity() != expectedIdentity || !sameVerifiedAssetOwner(request.context.verifiedAssets,expectedAssets))return "Stale audition request";
    if (! connected()) return "Open the selected device first";
    const auto& media = request.context.media;
    if (std::abs (media.sampleRate - device->getCurrentSampleRate()) > 0.01)
        return "Device and source sample rates must match";
    const auto routingError = routing.validate (media.channels, outputNames().size(), allowAllMuted);
    if (routingError.isNotEmpty()) return routingError;
    if(state==AuditionState::playing || state==AuditionState::buffering)stop();
    const int64_t wanted=restart || state==AuditionState::finished || heldFrame>=media.frames ? 0 : heldFrame;
    if(preparing() && (state!=AuditionState::preparing || preparationControl->isCancelled()
        || pendingRequest.identity()!=request.identity() || wanted!=startFrame))
        return "Wait for playback preparation cancellation, then press Play again";
    heldFrame=startFrame=wanted;
    pendingRouting = routing;
    pendingOutputs = outputNames();
    const auto preparationError = startPreparation (request);
    if (preparationError.isNotEmpty()) return preparationError;
    if (preparation.valid())
    { playAfterPreparation = true; message="Play requested - preparing playback; Stop or Cancel prevents playback"; }
    else startDecoder (request, routing);
    return {};
}

void StandaloneAudition::startDecoder (const AuditionRequest& request, const OutputRouting& routing)
{
    if(!playbackReady())
    {state=AuditionState::failed;message="Playback preparation changed; press Play to prepare again";return;}
    overloadWarning.clear();
    ++generation;presentationClock.reset();
    playingRequest = request; playingRouting = routing;
    const auto media = request.context.media;
    startFrame=heldFrame;
    fifo.reset (media.channels, routing);
    fifo.setObserver (capture, startFrame);
    fifo.push(readyResources->prefix.getArrayOfReadPointers(),readyResources->prefix.getNumSamples());
    cancel.store (false); decodeFailed.store (false); deviceStopped.store (false); workerError.clear(); workerFailure = ErrorCode::decodeFailed;
    decoder = std::thread ([this, request, resources=std::move(readyResources), media] ()
    {
        try
        {
            constexpr int blockSize = 2048;
            juce::AudioBuffer<float> block (media.channels, blockSize);
            juce::int64 position = resources->first+resources->prefix.getNumSamples();
            resources->prefix.setSize(0,0);
            while (position < media.frames && ! cancel.load())
            {
                const auto count = static_cast<int> (juce::jmin<juce::int64> (blockSize, media.frames - position));
                if (fifo.writable() < count) { std::this_thread::sleep_for (std::chrono::milliseconds (2)); continue; }
                if (! request.context.sourceFilesMatch() || !resources->source->matches(media))
                { workerFailure = ErrorCode::sourceChanged; throw std::runtime_error ("Source changed; import again"); }
                if (! resources->reader->read (block.getArrayOfWritePointers(), media.channels, position, count))
                    throw std::runtime_error ("Verified source decode failed");
                if (! resources->processing->process (block.getArrayOfWritePointers(), media.channels, count, position))
                    throw std::runtime_error ("Non-finite processed audition sample");
                fifo.push (block.getArrayOfReadPointers(), count);
                position += count;
            }
        }
        catch (const std::exception& e) { workerError = e.what(); decodeFailed.store (true); }
        fifo.endOfFile.store (true);
    });
    state = AuditionState::buffering;
    message = "Buffering " + juce::String (request.domain == MeasurementDomain::source ? "Source" : "Current State");
}

juce::String StandaloneAudition::seek (int64_t frame, const juce::File& file, const MediaInfo& media,
                                     const ChannelLayout& layout, const OutputRouting& routing)
{
    const auto request = AuditionRequest::source (file, media, layout);
    synchronize (request, routing);
    return seek (frame, request, routing);
}

juce::String StandaloneAudition::seek (int64_t frame, const AuditionRequest& request, const OutputRouting& routing)
{
    if (request.identity().isEmpty() || request.identity() != expectedIdentity
        || !sameVerifiedAssetOwner(request.context.verifiedAssets,expectedAssets)) return "Invalid or stale audition request";
    const auto& media = request.context.media;
    if (frame < 0 || frame > media.frames) return "Seek is outside the selected audio";
    const bool resume = active();
    const auto outputsBeforeSeek = pendingOutputs;
    stop();
    sourceRate = static_cast<int64_t> (media.sampleRate); sourceFrames = media.frames; sourceHash = media.sha256;
    heldFrame = startFrame = frame;
    if (frame == media.frames) { state = AuditionState::finished; message = "Audio EOF"; return {}; }
    if (!resume) return {};
    if (!connected() || outputsBeforeSeek!=outputNames() || std::abs(media.sampleRate-device->getCurrentSampleRate())>0.01
        || routing.validate(media.channels,outputNames().size()).isNotEmpty())
    {state=AuditionState::failed;message="Output configuration changed during seek; stopped";if(preparation.valid())preparationFailure=message;return message;}
    if (preparation.valid())
    {
        // Repeated live seeks replace this request; never overlap readers or
        // wait for a decoder on the message thread. The cancelled result cannot
        // publish its old prefix because Stop advanced the generation.
        pendingRequest=request;pendingRouting=routing;pendingOutputs=outputsBeforeSeek;
        resumeLatestSeek=true;latestSeekGeneration=generation;playAfterPreparation=true;
        state=AuditionState::preparing;message="Seeking - preparing latest position; Stop or Cancel prevents playback";
        return {};
    }
    return play(request,routing);
}

void StandaloneAudition::poll (bool picturePrerolled)
{
    if (preparation.valid())
    {
        if (playAfterPreparation && ! connected())
        {
            stop();
            preparationFailure = "Device closed during source preparation; stopped";
        }
        if (preparation.wait_for (std::chrono::seconds (0)) != std::future_status::ready) return;
        auto result = preparation.get();
        if (resumeLatestSeek)
        {
            resumeLatestSeek=false;playAfterPreparation=false;
            state=connected()?AuditionState::ready:AuditionState::disconnected;
            if(latestSeekGeneration!=generation || pendingRequest.identity()!=expectedIdentity
                || !sameVerifiedAssetOwner(pendingRequest.context.verifiedAssets,expectedAssets)
                || pendingRequest.validationError().isNotEmpty() || !pendingRequest.context.sourceFilesMatch()
                || !connected() || pendingOutputs!=outputNames()
                || std::abs(pendingRequest.context.media.sampleRate-device->getCurrentSampleRate())>0.01)
            {state=AuditionState::failed;message="Seek request or output changed; stopped. Press Play to retry.";return;}
            const auto error=play(pendingRequest,pendingRouting);
            if(error.isNotEmpty()){state=AuditionState::failed;message=error+"; stopped";}
            return; // The latest prefix and picture must both preroll first.
        }
        const bool accepted = preparationGeneration == generation && ! preparationControl->isCancelled()
            && pendingRequest.identity() == expectedIdentity;
        const bool atEnd=state==AuditionState::finished;
        state = atEnd ? AuditionState::finished : connected() ? AuditionState::ready : AuditionState::disconnected;
        if (! accepted)
        {
            playAfterPreparation = false;
            message = preparationFailure.isEmpty() ? (atEnd?"Audio EOF":"Source preparation cancelled; stopped") : preparationFailure;
            if (preparationFailure.isNotEmpty()) state = AuditionState::failed;
            return;
        }
        if (! result.value) { playAfterPreparation = false; fault = result.error; state = AuditionState::failed; message = result.message + "; stopped; session retained"; return; }
        readyResources = std::move(*result.value);
        preparedSource = readyResources->source; preparedRuntime = pendingRequest.context.sessionId; preparedAssets=pendingRequest.context.verifiedAssets;
        if (!playbackReady()) { playAfterPreparation=false;readyResources.reset();fault=ErrorCode::sourceChanged;state=AuditionState::failed;message="Source changed before playback preparation completed; stopped";return; }
        if (! playAfterPreparation) { message = "Playback ready; stopped. Press Play to listen."; return; }
        playAfterPreparation = false;
        if (!connected() || pendingOutputs!=outputNames() || std::abs(pendingRequest.context.media.sampleRate-device->getCurrentSampleRate())>0.01
            || pendingRouting.validate(pendingRequest.context.media.channels,outputNames().size(),true).isNotEmpty())
        { state=AuditionState::failed;message="Output configuration changed during preparation; stopped";return; }
        startDecoder (pendingRequest, pendingRouting);
        return; // The next UI tick observes the new transport generation before picture preroll authorizes output.
    }
    if (! active()) return;
    if((state==AuditionState::buffering || state==AuditionState::playing)
        && (playingRequest.validationError().isNotEmpty() || !playingRequest.context.sourceFilesMatch()
            || !preparedSource || !preparedSource->matches(playingRequest.context.media)))
    {
        stop();fault=ErrorCode::sourceChanged;state=AuditionState::failed;
        message="Source or session changed before playback; stopped. Import or reopen again.";return;
    }
    if (state == AuditionState::buffering && picturePrerolled && ! decodeFailed.load() && connected()
        && (fifo.writable() <= AuditionBuffer::capacity - 4096 || fifo.endOfFile.load()))
    {
        device->start (this);
        state = AuditionState::playing;
        message = "Auditioning " + juce::String (playingRequest.domain == MeasurementDomain::source ? "Source" : "Current State") + "; playback is not measurement or compliance";
    }
    const bool failed = fifo.overloaded() || decodeFailed.load() || deviceStopped.load() || ! connected()
        || (state == AuditionState::playing && ! device->isPlaying());
    if (failed || (fifo.endOfFile.load() && fifo.drained() && positionFrames() >= sourceFrames))
    {
        if(fifo.overloaded())
        {
            const auto channel=fifo.faultChannel();
            const auto frame=startFrame+static_cast<int64_t>(fifo.faultFrame());
            const auto role=channel>=0&&channel<playingRequest.context.fileLayout.count
                ?juce::AudioChannelSet::getChannelTypeName(playingRequest.context.fileLayout.roles[size_t(channel)]):juce::String("Unknown role");
            const auto output=channel>=0&&channel<maxChannels?playingRouting.fileToDevice[size_t(channel)]+1:0;
            overloadWarning="PLAYBACK STOPPED - Monitor overload: "+role+" (channel "+juce::String(channel+1)
                +", output "+juce::String(output)+"), "+juce::String(double(frame)/sourceRate,6)+" s (frame "+juce::String(frame)+").\n"
                +"Lower monitor volume in Audio setup, then press Play again. Processing gain is unchanged.";
        }
        if (decodeFailed.load()) fault = workerFailure;
        heldFrame = failed ? positionFrames() : sourceFrames;
        stopWorker();
        ++generation;
        state = failed ? AuditionState::failed : AuditionState::finished;
        message = fifo.overloaded() ? "Monitor overload at channel " + juce::String (fifo.faultChannel() + 1)
            + ", source frame " + juce::String (startFrame + static_cast<int64_t> (fifo.faultFrame()))
            + ". Reduce monitor volume and explicitly Play/Restart; processing gain retained."
            : failed ? (workerError.isNotEmpty() ? workerError : "Audition stopped: device unavailable") : "Audition finished";
    }
}

void StandaloneAudition::audioDeviceIOCallbackWithContext (const float* const*, int,
    float* const* outputs, int outputCount, int frames, const juce::AudioIODeviceCallbackContext& context)
{
    const auto callbackTime = hostTimeNs();
    if (outputCount < expectedOutputs) deviceStopped.store (true, std::memory_order_relaxed);
    if (fifo.overloaded() || decodeFailed.load (std::memory_order_relaxed) || deviceStopped.load (std::memory_order_relaxed))
    {
        for (int channel = 0; channel < outputCount; ++channel)
            if (outputs[channel]) juce::FloatVectorOperations::clear (outputs[channel], frames);
        return;
    }
    const auto before = fifo.playedFrames();
    fifo.render (outputs, outputCount, frames);
    if (fifo.overloaded()) return; // Suppressed callback never advances the presentation clock.
    const bool useHardware = hardwareLatencyFrames >= 0 && context.hostTimeNs != nullptr && *context.hostTimeNs != 0;
    const auto clockOrigin = useHardware ? *context.hostTimeNs : callbackTime;
    const auto delay = useHardware ? hardwareLatencyFrames : latencyFrames;
    const auto presentation = clockOrigin + static_cast<uint64_t> (static_cast<int64_t> (delay) * 1000000000 / sourceRate);
    hardwareClock.store (useHardware, std::memory_order_relaxed);
    presentationClock.submit ({ generation, presentation, startFrame + static_cast<int64_t> (before),
                                static_cast<int> (fifo.playedFrames() - before) });
}
// [END LS-STANDALONE-AUDITION-IMPLEMENTATION]
}

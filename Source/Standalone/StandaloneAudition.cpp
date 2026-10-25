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
    preparedSource.reset(); expectedIdentity.clear();
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

juce::StringArray StandaloneAudition::outputNames() const
{
    if (! device) return {};
    auto names = device->getOutputChannelNames();
    while (names.size() > maxChannels) names.remove (names.size() - 1);
    return names;
}

juce::String StandaloneAudition::connect (double sourceRate)
{
    if (preparation.valid()) return "Wait for source preparation or cancellation before opening a device";
    stop();
    if (! device || ! std::isfinite (sourceRate) || sourceRate <= 0) return "Select a device and validated source";
    device->close();
    juce::BigInteger outputs;
    outputs.setRange (0, outputNames().size(), true);
    if (outputs.isZero()) return "Device has no available outputs";
    if (! device->getAvailableSampleRates().contains (sourceRate))
        return "Device does not offer the source sample rate; resampling is unavailable";
    auto error = device->open ({}, outputs, sourceRate, device->getDefaultBufferSize());
    if (error.isNotEmpty() || ! device->isOpen() || std::abs (device->getCurrentSampleRate() - sourceRate) > 0.01
        || device->getActiveOutputChannels() != outputs)
    {
        device->close();
        state = AuditionState::failed;
        message = error.isNotEmpty() ? error : "Device did not accept the requested sample rate and outputs";
        return message;
    }
    rate = sourceRate;
    latencyFrames = std::max (0, device->getOutputLatencyInSamples());
#if JUCE_MAC
    hardwareLatencyFrames = hardwareDelay (*device);
#endif
    (void) hostTimeNs(); // initialize timebase on the message thread, never first-use in RT
    expectedOutputs = outputs.countNumberOfSetBits();
    state = AuditionState::ready;
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

void StandaloneAudition::synchronize (const AuditionRequest& request, const OutputRouting& routing)
{
    const auto key = request.identity();
    if (active() && (key != expectedIdentity
        || ((state == AuditionState::playing || state == AuditionState::buffering) && routing.fileToDevice != playingRouting.fileToDevice)
        || (state == AuditionState::preparing && playAfterPreparation && routing.fileToDevice != pendingRouting.fileToDevice)))
    {
        stop(); message = "Audition stopped: source, roles, processing or routing changed";
    }
    expectedIdentity = key;
    if (preparedSource && ! preparedSource->matches (request.context.media)) preparedSource.reset();
}

juce::String StandaloneAudition::startPreparation (const AuditionRequest& request)
{
    const auto error = request.validationError();
    if (error.isNotEmpty()) return error;
    if (request.identity() != expectedIdentity) return "Stale audition request; select the current state again";
    if (preparation.valid()) return "Wait for source preparation cancellation";
    if (! request.context.sourceFilesMatch())
    { fault = ErrorCode::sourceChanged; return "Source changed; import again"; }
    fault = ErrorCode::none;
    pendingRequest = request;
    if (preparedSource && preparedSource->matches (request.context.media)) return {};
    preparationControl = std::make_shared<JobControl>();
    preparationFailure.clear();
    preparationGeneration = generation;
    preparation = std::async (std::launch::async, [request, storage = snapshotStorage, control = preparationControl]
    { return prepareVerifiedSource (request.context, *control, storage); });
    state = AuditionState::preparing;
    message = "Preparing verified source copy (copy / hash / payload validation)";
    return {};
}

juce::String StandaloneAudition::prepare (const AuditionRequest& request)
{
    if (active()) return "Stop before preparing a different audition";
    stop();
    const auto error = startPreparation (request);
    if (error.isEmpty() && ! preparation.valid()) message = "Verified source ready; device and playback unchanged";
    return error;
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
    if (request.identity() != expectedIdentity) return "Stale audition request";
    if (preparation.valid()) return "Wait for preparation or cancel it";
    const bool atEnd = state == AuditionState::finished;
    stop();
    if (! connected()) return "Open the selected device first";
    const auto& media = request.context.media;
    if (std::abs (media.sampleRate - device->getCurrentSampleRate()) > 0.01)
        return "Device and source sample rates must match";
    const auto routingError = routing.validate (media.channels, outputNames().size(), allowAllMuted);
    if (routingError.isNotEmpty()) return routingError;
    if (sourceHash != media.sha256) { heldFrame = startFrame = 0; if (preparedSource && ! preparedSource->matches (media)) preparedSource.reset(); }
    sourceHash = media.sha256; sourceRate = static_cast<int64_t> (media.sampleRate); sourceFrames = media.frames;
    if (restart || atEnd || heldFrame >= sourceFrames) heldFrame = 0;
    startFrame = heldFrame;
    pendingRouting = routing;
    const auto preparationError = startPreparation (request);
    if (preparationError.isNotEmpty()) return preparationError;
    if (preparation.valid()) playAfterPreparation = true;
    else startDecoder (request, routing);
    return {};
}

void StandaloneAudition::startDecoder (const AuditionRequest& request, const OutputRouting& routing)
{
    playingRequest = request; playingRouting = routing;
    const auto media = request.context.media;
    fifo.reset (media.channels, routing);
    fifo.setObserver (capture, startFrame);
    cancel.store (false); decodeFailed.store (false); deviceStopped.store (false); workerError.clear(); workerFailure = ErrorCode::decodeFailed;
    decoder = std::thread ([this, request, copy = preparedSource, media, begin = startFrame] ()
    {
        try
        {
            auto source = copy->createReader();
            if (! source || source->sampleRate != media.sampleRate || source->lengthInSamples != media.frames
                || static_cast<int> (source->numChannels) != media.channels)
                throw std::runtime_error ("Verified source cannot be reopened consistently");
            constexpr int blockSize = 2048;
            juce::AudioBuffer<float> block (media.channels, blockSize);
            // JUCE 7's FLAC reader can scan for tens of seconds when a fresh
            // reader seeks into the final partial block. Seek before that block
            // and decode a bounded prefix instead. Discarded PCM never reaches
            // processing, the FIFO, capture observers or the presentation clock.
            // 65536 exceeds FLAC's maximum block size; cancellation stays checked
            // between normal 2048-frame reads. Other formats retain their path.
            if (source->hasFlac && begin > 0)
                for (juce::int64 prefix = std::max<juce::int64> (0, begin - 65536);
                     prefix < begin && ! cancel.load();)
                {
                    const int count = static_cast<int> (std::min<juce::int64> (blockSize, begin - prefix));
                    if (! source->read (block.getArrayOfWritePointers(), media.channels, prefix, count))
                        throw std::runtime_error ("Verified source seek pre-roll failed");
                    prefix += count;
                }
            const StaticGainProcessor processing (request.context.processing.staticGainDb);
            juce::int64 position = begin;
            while (position < media.frames && ! cancel.load())
            {
                const auto count = static_cast<int> (juce::jmin<juce::int64> (blockSize, media.frames - position));
                if (fifo.writable() < count) { std::this_thread::sleep_for (std::chrono::milliseconds (2)); continue; }
                if (! request.context.sourceFilesMatch())
                { workerFailure = ErrorCode::sourceChanged; throw std::runtime_error ("Source changed; import again"); }
                if (! source->read (block.getArrayOfWritePointers(), media.channels, position, count))
                    throw std::runtime_error ("Verified source decode failed");
                if (! processing.process (block.getArrayOfWritePointers(), media.channels, count))
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
    if (request.identity().isEmpty() || request.identity() != expectedIdentity) return "Invalid or stale audition request";
    const auto& media = request.context.media;
    if (frame < 0 || frame > media.frames) return "Seek is outside the selected audio";
    if (preparation.valid()) return "Wait for preparation or cancel it before seeking";
    const bool resume = active();
    stop();
    sourceRate = static_cast<int64_t> (media.sampleRate); sourceFrames = media.frames; sourceHash = media.sha256;
    heldFrame = startFrame = frame;
    if (frame == media.frames) { state = AuditionState::finished; message = "Audio EOF"; return {}; }
    return resume ? play (request, routing) : juce::String();
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
        const bool accepted = preparationGeneration == generation && ! preparationControl->isCancelled()
            && pendingRequest.identity() == expectedIdentity;
        state = connected() ? AuditionState::ready : AuditionState::disconnected;
        if (! accepted)
        {
            playAfterPreparation = false;
            message = preparationFailure.isEmpty() ? "Source preparation cancelled; stopped" : preparationFailure;
            if (preparationFailure.isNotEmpty()) state = AuditionState::failed;
            return;
        }
        if (! result.value) { playAfterPreparation = false; fault = result.error; state = AuditionState::failed; message = result.message + "; stopped; session retained"; return; }
        preparedSource = *result.value;
        if (! playAfterPreparation) { message = "Verified source ready; playback stopped"; return; }
        playAfterPreparation = false;
        if (! connected()) { state = AuditionState::failed; message = "Device closed during preparation; stopped"; return; }
        startDecoder (pendingRequest, pendingRouting);
    }
    if (! active()) return;
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

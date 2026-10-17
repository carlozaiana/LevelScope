#include "StandaloneAudition.h"
#include <cmath>
#include <chrono>
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
    lastGain = 0;
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
    const auto start = read.load (std::memory_order_relaxed);
    const auto available = written.load (std::memory_order_acquire) - start;
    const auto count = static_cast<int> (juce::jmin<uint64_t> (static_cast<uint64_t> (frames), available));
    if (count < frames && ! endOfFile.load()) underruns.fetch_add (1, std::memory_order_relaxed);
    const auto target = gain.load (std::memory_order_relaxed);
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

StandaloneAudition::~StandaloneAudition() { disconnect(); }

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
    stopWorker();
    ++generation;
    presentationClock.reset();
    state = connected() ? AuditionState::ready : AuditionState::disconnected;
    message = connected() ? "Stopped; Play resumes, Restart begins at zero" : "Device closed";
}

void StandaloneAudition::disconnect()
{
    stop();
    if (device) device->close();
    rate = 0;
    state = AuditionState::disconnected;
    message = "Device closed";
}

juce::String StandaloneAudition::play (const juce::File& file, const MediaInfo& media,
                                     const ChannelLayout& layout, const OutputRouting& routing, bool restart, bool allowAllMuted)
{
    const bool atEnd = state == AuditionState::finished;
    stop();
    if (! connected()) return "Open the selected device first";
    if (media.channels != layout.count || ! layout.confirmed || ! layout.isValid())
        return "Confirm valid file-channel roles first";
    if (! media.matchesFileStamp (file)) return "Source changed; import again";
    if (std::abs (media.sampleRate - device->getCurrentSampleRate()) > 0.01)
        return "Device and source sample rates must match";
    const auto routingError = routing.validate (media.channels, outputNames().size(), allowAllMuted);
    if (routingError.isNotEmpty()) return routingError;
    if (media.sampleRate != std::floor (media.sampleRate)) return "Transport requires an integer source sample rate";
    if (sourceHash != media.sha256) sourceChanged (media);
    sourceRate = static_cast<int64_t> (media.sampleRate); sourceFrames = media.frames;
    if (restart || atEnd || heldFrame >= sourceFrames) heldFrame = 0;
    startFrame = heldFrame;
    juce::AudioFormatManager formats;
    registerSourceFormats (formats);
    std::unique_ptr<juce::AudioFormatReader> reader (formats.createReaderFor (file.createInputStream()));
    if (! reader || reader->sampleRate != media.sampleRate || reader->lengthInSamples != media.frames
        || static_cast<int> (reader->numChannels) != media.channels) return "Source cannot be reopened consistently";

    fifo.reset (media.channels, routing);
    cancel.store (false); decodeFailed.store (false); deviceStopped.store (false);
    decoder = std::thread ([this, source = std::move (reader), media, file] () mutable
    {
        constexpr int blockSize = 2048;
        juce::AudioBuffer<float> block (media.channels, blockSize);
        juce::int64 position = startFrame;
        while (position < media.frames && ! cancel.load())
        {
            const auto count = static_cast<int> (juce::jmin<juce::int64> (blockSize, media.frames - position));
            if (fifo.writable() < count)
            {
                std::this_thread::sleep_for (std::chrono::milliseconds (2));
                continue;
            }
            if (! media.matchesFileStamp (file)
                || ! source->read (block.getArrayOfWritePointers(), media.channels, position, count))
            {
                decodeFailed.store (true);
                break;
            }
            bool finite = true;
            for (int channel = 0; channel < media.channels; ++channel)
                for (int frame = 0; frame < count; ++frame)
                    finite = finite && std::isfinite (block.getSample (channel, frame));
            if (! finite) { decodeFailed.store (true); break; }
            fifo.push (block.getArrayOfReadPointers(), count);
            position += count;
        }
        fifo.endOfFile.store (true);
    });
    state = AuditionState::buffering;
    message = "Buffering source";
    return {};
}

juce::String StandaloneAudition::seek (int64_t frame, const juce::File& file, const MediaInfo& media,
                                     const ChannelLayout& layout, const OutputRouting& routing)
{
    if (frame < 0 || frame > media.frames) return "Seek is outside the selected audio";
    const bool resume = active();
    stop();
    sourceRate = static_cast<int64_t> (media.sampleRate); sourceFrames = media.frames; sourceHash = media.sha256;
    heldFrame = startFrame = frame;
    if (frame == media.frames) { state = AuditionState::finished; message = "Audio EOF"; return {}; }
    return resume ? play (file, media, layout, routing) : juce::String();
}

void StandaloneAudition::poll (bool picturePrerolled)
{
    if (! active()) return;
    if (state == AuditionState::buffering && picturePrerolled && ! decodeFailed.load() && connected()
        && (fifo.writable() <= AuditionBuffer::capacity - 4096 || fifo.endOfFile.load()))
    {
        device->start (this);
        state = AuditionState::playing;
        message = "Auditioning original source (monitor gain only)";
    }
    const bool failed = decodeFailed.load() || deviceStopped.load() || ! connected()
        || (state == AuditionState::playing && ! device->isPlaying());
    if (failed || (fifo.endOfFile.load() && fifo.drained() && positionFrames() >= sourceFrames))
    {
        heldFrame = failed ? positionFrames() : sourceFrames;
        stopWorker();
        ++generation;
        state = failed ? AuditionState::failed : AuditionState::finished;
        message = failed ? "Audition stopped: source read failed or device stopped" : "Audition finished";
    }
}

void StandaloneAudition::audioDeviceIOCallbackWithContext (const float* const*, int,
    float* const* outputs, int outputCount, int frames, const juce::AudioIODeviceCallbackContext& context)
{
    const auto callbackTime = hostTimeNs();
    if (outputCount < expectedOutputs) deviceStopped.store (true, std::memory_order_relaxed);
    if (decodeFailed.load (std::memory_order_relaxed) || deviceStopped.load (std::memory_order_relaxed))
    {
        for (int channel = 0; channel < outputCount; ++channel)
            if (outputs[channel]) juce::FloatVectorOperations::clear (outputs[channel], frames);
        return;
    }
    const auto before = fifo.playedFrames();
    fifo.render (outputs, outputCount, frames);
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

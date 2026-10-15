#include "StandaloneAudition.h"
#include <cmath>
#include <chrono>

namespace levelscope::standalone
{
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
    stopWorker();
    state = connected() ? AuditionState::ready : AuditionState::disconnected;
    message = connected() ? "Stopped; Play starts at the beginning" : "Device closed";
}

void StandaloneAudition::disconnect()
{
    stopWorker();
    if (device) device->close();
    rate = 0;
    state = AuditionState::disconnected;
    message = "Device closed";
}

juce::String StandaloneAudition::play (const juce::File& file, const MediaInfo& media,
                                     const ChannelLayout& layout, const OutputRouting& routing)
{
    stop();
    if (! connected()) return "Open the selected device first";
    if (media.channels != layout.count || ! layout.confirmed || ! layout.isValid())
        return "Confirm valid file-channel roles first";
    if (! media.matchesFileStamp (file)) return "Source changed; import again";
    if (std::abs (media.sampleRate - device->getCurrentSampleRate()) > 0.01)
        return "Device and source sample rates must match";
    const auto routingError = routing.validate (media.channels, outputNames().size());
    if (routingError.isNotEmpty()) return routingError;
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
        juce::int64 position = 0;
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

void StandaloneAudition::poll()
{
    if (! active()) return;
    if (state == AuditionState::buffering && ! decodeFailed.load() && connected()
        && (fifo.writable() <= AuditionBuffer::capacity - 4096 || fifo.endOfFile.load()))
    {
        device->start (this);
        state = AuditionState::playing;
        message = "Auditioning original source (monitor gain only)";
    }
    const bool failed = decodeFailed.load() || deviceStopped.load() || ! connected()
        || (state == AuditionState::playing && ! device->isPlaying());
    if (failed || (fifo.endOfFile.load() && fifo.drained()))
    {
        stopWorker();
        state = failed ? AuditionState::failed : AuditionState::finished;
        message = failed ? "Audition stopped: source read failed or device stopped" : "Audition finished";
    }
}

void StandaloneAudition::audioDeviceIOCallbackWithContext (const float* const*, int,
    float* const* outputs, int outputCount, int frames, const juce::AudioIODeviceCallbackContext&)
{
    if (outputCount < expectedOutputs) deviceStopped.store (true, std::memory_order_relaxed);
    if (decodeFailed.load (std::memory_order_relaxed) || deviceStopped.load (std::memory_order_relaxed))
    {
        for (int channel = 0; channel < outputCount; ++channel)
            if (outputs[channel]) juce::FloatVectorOperations::clear (outputs[channel], frames);
        return;
    }
    fifo.render (outputs, outputCount, frames);
}
// [END LS-STANDALONE-AUDITION-IMPLEMENTATION]
}

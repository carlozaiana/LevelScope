#pragma once

#include "StandaloneMedia.h"
#include <juce_audio_devices/juce_audio_devices.h>
#include <thread>

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-AUDITION]
// Single decoder producer / single audio consumer. Reset/configure only while
// both are stopped. Callback performs no I/O, allocation, locks or notifications.
class AuditionBuffer
{
public:
    static constexpr int capacity = 32768;
    void reset (int channels, const OutputRouting&);
    int writable() const noexcept;
    void push (const float* const*, int frames) noexcept;
    void render (float* const* outputs, int outputCount, int frames) noexcept;
    bool drained() const noexcept { return read.load() == written.load(); }
    uint64_t playedFrames() const noexcept { return read.load(); }
    uint64_t underrunBlocks() const noexcept { return underruns.load(); }
    std::atomic<bool> endOfFile { false };
    std::atomic<float> gain { 0.03162278f }; // -30 dB monitor only
private:
    static_assert (std::atomic<uint64_t>::is_always_lock_free && std::atomic<float>::is_always_lock_free
                   && std::atomic<bool>::is_always_lock_free,
                   "Audition requires lock-free counters and gain");
    std::array<std::array<float, capacity>, maxChannels> samples;
    std::atomic<uint64_t> read { 0 }, written { 0 }, underruns { 0 };
    OutputRouting routing;
    int channelCount = 0;
    float lastGain = 0.0f;
};

enum class AuditionState { disconnected, ready, buffering, playing, finished, failed };

class StandaloneAudition final : private juce::AudioIODeviceCallback
{
public:
    StandaloneAudition() = default;
    explicit StandaloneAudition (std::unique_ptr<juce::AudioIODevice> selectedDevice)
        : device (std::move (selectedDevice)) {} // injectable device for headless checks
    ~StandaloneAudition() override;
    const juce::OwnedArray<juce::AudioIODeviceType>& deviceTypes();
    juce::StringArray scanOutputs (int typeIndex);
    juce::String selectDevice (int typeIndex, const juce::String& name);
    juce::StringArray outputNames() const;
    juce::String connect (double sourceRate);
    void disconnect();
    juce::String play (const juce::File&, const MediaInfo&, const ChannelLayout&, const OutputRouting&);
    void stop();
    void poll(); // message thread; finish/failed notifications and worker join
    void setMonitorGainDb (float db) { fifo.gain.store (juce::Decibels::decibelsToGain (juce::jlimit (-60.0f, 0.0f, db))); }
    bool connected() const { return device != nullptr && device->isOpen(); }
    AuditionState getState() const { return state; }
    bool active() const { return state == AuditionState::buffering || state == AuditionState::playing; }
    double positionSeconds() const { return rate > 0 ? static_cast<double> (fifo.playedFrames()) / rate : 0; }
    uint64_t underruns() const { return fifo.underrunBlocks(); }
    const juce::String& status() const { return message; }
private:
    // Lazy factory: injected/headless devices never construct the host manager.
    std::unique_ptr<juce::AudioDeviceManager> catalog;
    std::unique_ptr<juce::AudioIODevice> device;
    AuditionBuffer fifo;
    std::thread decoder;
    std::atomic<bool> cancel { false }, decodeFailed { false }, deviceStopped { false };
    AuditionState state = AuditionState::disconnected;
    juce::String message = "Device closed. Select an output explicitly.";
    double rate = 0;
    int expectedOutputs = 0;
    void stopWorker();
    void audioDeviceIOCallbackWithContext (const float* const*, int, float* const*, int, int,
                                          const juce::AudioIODeviceCallbackContext&) override;
    void audioDeviceAboutToStart (juce::AudioIODevice* opened) override
    { deviceStopped.store (std::abs (opened->getCurrentSampleRate() - rate) > 0.01); }
    void audioDeviceStopped() override { deviceStopped.store (true); }
    void audioDeviceError (const juce::String&) override { deviceStopped.store (true); }
};
// [END LS-STANDALONE-AUDITION]
}

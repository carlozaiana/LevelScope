#pragma once

#include "StandaloneMedia.h"
#include "AudioPresentationClock.h"
#include "StandaloneAuditionRequest.h"
#include "VerifiedSource.h"
#include <juce_audio_devices/juce_audio_devices.h>
#include <thread>
#include <future>

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
    bool overloaded() const noexcept { return monitorOverload.load(); }
    int faultChannel() const { return overloadChannel.load(); }
    uint64_t faultFrame() const { return overloadFrame.load(); }
    void setObserver (ProcessedAudioObserver* sink, int64_t sourceStart) { observer = sink; captureStart = sourceStart; }
private:
    static_assert (std::atomic<uint64_t>::is_always_lock_free && std::atomic<float>::is_always_lock_free
                   && std::atomic<bool>::is_always_lock_free,
                   "Audition requires lock-free counters and gain");
    std::array<std::array<float, capacity>, maxChannels> samples;
    std::atomic<uint64_t> read { 0 }, written { 0 }, underruns { 0 };
    OutputRouting routing;
    int channelCount = 0;
    float lastGain = 0.0f;
    std::atomic<bool> monitorOverload { false };
    std::atomic<int> overloadChannel { -1 };
    std::atomic<uint64_t> overloadFrame { 0 };
    ProcessedAudioObserver* observer = nullptr;
    int64_t captureStart = 0;
};

enum class AuditionState { disconnected, ready, preparing, cancelling, buffering, playing, finished, failed };

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
    void selectOutputDevice(std::unique_ptr<juce::AudioIODevice> selected);
    juce::StringArray outputNames() const;
    juce::StringArray outputConfiguration() const; // complete ordered metadata, also while disconnected
    juce::String connect (double sourceRate);
    void disconnect();
    // Authoritative session binding is refreshed before polling/start. Unexpected
    // changes stop active output; prepared jobs cannot publish for an older key.
    void synchronize (const AuditionRequest&, const OutputRouting&);
    juce::String prepare (const AuditionRequest&);
    juce::String play (const AuditionRequest&, const OutputRouting&, bool restart = false, bool allowAllMuted = false);
    juce::String seek (int64_t frame, const AuditionRequest&, const OutputRouting&);
    juce::String play (const juce::File&, const MediaInfo&, const ChannelLayout&, const OutputRouting&, bool restart = false, bool allowAllMuted = false);
    juce::String seek (int64_t frame, const juce::File&, const MediaInfo&, const ChannelLayout&, const OutputRouting&);
    void sourceChanged (const MediaInfo&);
    void stop();
    void poll (bool picturePrerolled = true); // message thread; finish/failed notifications and worker join
    void setMonitorGainDb (float db) { if (std::isfinite (db)) fifo.gain.store (juce::Decibels::decibelsToGain (juce::jlimit (-60.0f, 0.0f, db))); }
    bool connected() const { return device != nullptr && device->isOpen(); }
    double connectedRate() const { return connected() ? device->getCurrentSampleRate() : 0; }
    AuditionState getState() const { return state; }
    // Stop revokes transport intent before an asynchronous preparation finishes
    // cleaning up. Its future lifetime is guarded separately by preparing().
    bool active() const { return (state == AuditionState::preparing && playAfterPreparation)
                             || state == AuditionState::buffering || state == AuditionState::playing; }
    bool preparing() const { return preparation.valid(); }
    bool backgroundPreparing() const { return preparing() && !playAfterPreparation && state == AuditionState::preparing; }
    bool playbackReady() const;
    juce::String preparationStage() const;
    float preparationProgress() const { return decoderPreparation && decoderPreparation->load() ? 0.95f
        : preparationControl ? juce::jlimit(0.0f,0.9f,preparationControl->progress.load()*1.8f) : 0; }
    const juce::String& boundIdentity() const { return expectedIdentity; }
    const VerifiedSource* preparedCopy() const { return preparedSource.get(); }
    std::shared_ptr<const VerifiedSource> sharedPreparedCopy() const { return preparedSource; }
    bool setCaptureObserver (ProcessedAudioObserver* sink) { if (active()) return false; capture = sink; return true; }
    bool setSnapshotStorage (SnapshotStorage storage) { if (active() || preparing()) return false; snapshotStorage = std::move (storage); return true; }
    int64_t positionFrames() const;
    double positionSeconds() const { return sourceRate > 0 ? static_cast<double> (positionFrames()) / sourceRate : 0; }
    uint64_t transportGeneration() const { return generation; }
    int outputLatencyFrames() const { return latencyFrames; }
    bool hasHardwareClock() const { return hardwareClock.load(); }
    static uint64_t hostTimeNs() noexcept;
    uint64_t underruns() const { return fifo.underrunBlocks(); }
    const juce::String& status() const { return message; }
    const juce::String& monitorOverloadWarning() const { return overloadWarning; }
    ErrorCode failureCode() const { return fault; }
private:
    // Lazy factory: injected/headless devices never construct the host manager.
    std::unique_ptr<juce::AudioDeviceManager> catalog;
    std::unique_ptr<juce::AudioIODevice> device;
    AuditionBuffer fifo;
    std::thread decoder;
    std::atomic<bool> cancel { false }, decodeFailed { false }, deviceStopped { false };
    AuditionState state = AuditionState::disconnected;
    juce::String message = "Device closed. Select an output explicitly.";
    juce::String overloadWarning; // Message-thread latch; monitor adjustment alone does not restart playback.
    double rate = 0;
    int64_t sourceRate = 0, startFrame = 0, heldFrame = 0, sourceFrames = 0;
    juce::String sourceHash;
    uint64_t generation = 1;
    AudioPresentationClock presentationClock;
    int latencyFrames = 0, hardwareLatencyFrames = -1;
    std::atomic<bool> hardwareClock { false };
    int expectedOutputs = 0;
    // A single reader, processor and at most 4096 frames; never whole-programme
    // PCM. Worker prepares these without touching the callback or any device.
    struct PlaybackResources
    {
        std::shared_ptr<const VerifiedSource> source;
        std::unique_ptr<ProgrammeReader> reader;
        std::unique_ptr<ProgrammeGainProcessor> processing;
        juce::AudioBuffer<float> prefix;
        juce::String identity;
        int64_t first = 0;
    };
    std::shared_ptr<const VerifiedSource> preparedSource;
    juce::String preparedRuntime;
    std::optional<std::weak_ptr<SessionVerifiedAssets>> preparedAssets;
    std::optional<std::weak_ptr<SessionVerifiedAssets>> expectedAssets;
    std::unique_ptr<PlaybackResources> readyResources;
    std::future<ServiceResult<std::unique_ptr<PlaybackResources>>> preparation;
    std::shared_ptr<JobControl> preparationControl;
    std::shared_ptr<std::atomic<bool>> decoderPreparation;
    SnapshotStorage snapshotStorage;
    AuditionRequest pendingRequest, playingRequest;
    OutputRouting pendingRouting, playingRouting;
    juce::StringArray pendingOutputs;
    juce::String expectedIdentity, workerError, preparationFailure;
    ErrorCode fault = ErrorCode::none, workerFailure = ErrorCode::none;
    uint64_t preparationGeneration = 0;
    bool playAfterPreparation = false;
    ProcessedAudioObserver* capture = nullptr;
    juce::String startPreparation (const AuditionRequest&);
    void startDecoder (const AuditionRequest&, const OutputRouting&);
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

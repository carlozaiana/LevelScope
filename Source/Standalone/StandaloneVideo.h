#pragma once
#include "VideoTiming.h"
#include <juce_graphics/juce_graphics.h>
#include <atomic>
#include <condition_variable>
#include <deque>
#include <functional>
#include <mutex>
#include <thread>

namespace levelscope::standalone::video
{
struct Stream
{
    int id = 0, width = 0, height = 0, timeScale = 0;
    juce::String codec;
    Time nominalRate;
    bool supported = false;
};
struct Asset
{
    juce::File file;
    juce::String sha256, backend = "AVAssetReader / BGRA / JUCE image";
    int64_t bytes = 0, modifiedMs = 0;
    std::vector<Stream> streams;
    int selectedStream = 0, embeddedAudioTracks = 0;
    bool variableRate = false;
    std::optional<Timecode> timecode;
    juce::String timecodeStatus = "No timecode track";
    std::vector<FrameTime> frames;
    bool unchanged() const
    { return file.existsAsFile() && file.getSize() == bytes && file.getLastModificationTime().toMilliseconds() == modifiedMs; }
};
struct Picture { juce::Image image; size_t index = 0; uint64_t generation = 0, revision = 0; };
using Cancelled = std::function<bool()>;
class Decoder
{
public:
    virtual ~Decoder() = default;
    virtual std::shared_ptr<const Asset> load (const juce::File&, int streamId, const Cancelled&) = 0;
    virtual void seek (size_t index) = 0;
    virtual std::optional<Picture> next() = 0;
};
std::unique_ptr<Decoder> makeDecoder();

// Session reference state has no dependency on measurement revisions or routing.
struct Reference
{
    juce::File file;
    juce::String sha256;
    int streamId = 0;
    Time startOnAudio;
};

class Player
{
public:
    enum class State { empty, loading, ready, buffering, failed };
    struct View
    {
        State state = State::empty;
        juce::String message = "No reference picture";
        std::shared_ptr<const Asset> asset;
        std::optional<Picture> picture;
        uint64_t latePolls = 0, droppedFrames = 0;
        size_t queued = 0;
    };
    Player();
    ~Player();
    void attach (const juce::File&, int streamId = 0);
    void clear();
    // Message-thread only. Revision rejects replacement/offset completions; the
    // transport generation rejects pre-seek/pre-stop completions independently.
    void invalidate();
    View update (Time videoTime, uint64_t transportGeneration);
    bool prerollReady() const;
private:
    void run();
    mutable std::mutex mutex;
    std::condition_variable wake;
    std::thread worker;
    std::atomic<bool> closing { false };
    std::atomic<uint64_t> contentEpoch { 1 }, viewEpoch { 1 };
    juce::File requestedFile;
    int requestedStream = 0;
    uint64_t generation = 0;
    Time requestedTime;
    View view;
    std::optional<size_t> desired, lastDisplayed;
    std::deque<Picture> queue;
    static constexpr size_t maxQueued = 8;
};
}

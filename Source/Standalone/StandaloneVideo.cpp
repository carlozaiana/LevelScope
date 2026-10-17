#include "StandaloneVideo.h"

namespace levelscope::standalone::video
{
Player::Player() { worker = std::thread ([this] { run(); }); }
Player::~Player()
{
    closing.store (true); ++contentEpoch; ++viewEpoch; wake.notify_one();
    if (worker.joinable()) worker.join();
}
void Player::attach (const juce::File& file, int streamId)
{
    std::lock_guard<std::mutex> lock (mutex);
    ++contentEpoch; ++viewEpoch;
    requestedFile = file; requestedStream = streamId;
    queue.clear(); desired.reset(); lastDisplayed.reset(); view = {};
    view.state = State::loading; view.message = "Indexing reference picture; embedded audio excluded";
    wake.notify_one();
}
void Player::clear()
{
    std::lock_guard<std::mutex> lock (mutex);
    ++contentEpoch; ++viewEpoch; requestedFile = juce::File(); queue.clear();
    desired.reset(); lastDisplayed.reset(); view = {}; wake.notify_one();
}
void Player::invalidate()
{
    std::lock_guard<std::mutex> lock (mutex);
    ++viewEpoch; queue.clear(); view.picture.reset(); lastDisplayed.reset(); wake.notify_one();
}
Player::View Player::update (Time time, uint64_t transportGeneration)
{
    std::lock_guard<std::mutex> lock (mutex);
    requestedTime = time;
    if (generation != transportGeneration)
    {
        generation = transportGeneration; ++viewEpoch; queue.clear(); view.picture.reset(); lastDisplayed.reset();
    }
    if (view.asset && view.state != State::failed)
    {
        desired = frameAt (view.asset->frames, time);
        while (!queue.empty() && (!desired || queue.front().index < *desired)) queue.pop_front();
        view.picture.reset();
        if (!desired) { view.state = State::ready; view.message = "Outside reference-picture coverage"; lastDisplayed.reset(); }
        else
        {
            for (const auto& picture : queue)
                if (picture.index == *desired && picture.generation == generation && picture.revision == viewEpoch.load())
                    view.picture = picture;
            if (view.picture)
            {
                if (lastDisplayed && *desired > *lastDisplayed + 1) view.droppedFrames += *desired - *lastDisplayed - 1;
                lastDisplayed = desired;
                view.state = State::ready; view.message = "Reference picture - embedded audio excluded";
            }
            else { ++view.latePolls; view.state = State::buffering; view.message = "Picture buffering (audio clock remains authoritative)"; }
        }
    }
    view.queued = queue.size(); wake.notify_one(); return view;
}
bool Player::prerollReady() const
{
    std::lock_guard<std::mutex> lock (mutex);
    return view.state == State::empty || view.state == State::failed || view.state == State::ready;
}
void Player::run()
{
    std::unique_ptr<Decoder> decoder;
    uint64_t loadedEpoch = 0, decodedEpoch = 0;
    size_t nextIndex = 0;
    while (!closing.load())
    {
        uint64_t content, revision, ticket;
        juce::File file; int stream;
        std::optional<size_t> target;
        std::shared_ptr<const Asset> asset;
        {
            std::unique_lock<std::mutex> lock (mutex);
            wake.wait_for (lock, std::chrono::milliseconds (5), [this, &loadedEpoch]
                { return closing.load() || contentEpoch.load() != loadedEpoch; });
            if (closing.load()) break;
            content = contentEpoch.load(); revision = viewEpoch.load(); ticket = generation;
            file = requestedFile; stream = requestedStream; target = desired; asset = view.asset;
            if (content == loadedEpoch && (view.state == State::failed || !target || queue.size() >= maxQueued)) continue;
        }
        try
        {
            if (content != loadedEpoch)
            {
                decoder.reset(); loadedEpoch = content; decodedEpoch = 0;
                if (file == juce::File()) continue;
                decoder = makeDecoder();
                asset = decoder->load (file, stream, [this, content] { return closing.load() || contentEpoch.load() != content; });
                std::lock_guard<std::mutex> lock (mutex);
                if (contentEpoch.load() != content) continue;
                view.asset = asset;
                desired = frameAt (asset->frames, requestedTime);
                view.state = desired ? State::buffering : State::ready;
                continue;
            }
            if (!asset || !decoder || !target) continue;
            if (!asset->unchanged()) throw std::runtime_error ("Reference picture changed or is missing; attach it again");
            if (decodedEpoch != revision || nextIndex > *target + maxQueued || nextIndex < *target)
            {
                decoder->seek (*target); nextIndex = *target; decodedEpoch = revision;
            }
            if (nextIndex >= asset->frames.size()) continue;
            auto picture = decoder->next();
            if (!picture) throw std::runtime_error ("Video decoder ended before the indexed EOF");
            if (picture->index != nextIndex) throw std::runtime_error ("Decoder PTS does not match the indexed source frame");
            ++nextIndex;
            std::lock_guard<std::mutex> lock (mutex);
            if (contentEpoch.load() != content || viewEpoch.load() != revision || generation != ticket) continue;
            picture->generation = ticket; picture->revision = revision;
            if (desired && picture->index >= *desired && queue.size() < maxQueued) queue.push_back (std::move (*picture));
        }
        catch (const std::exception& error)
        {
            decoder.reset();
            std::lock_guard<std::mutex> lock (mutex);
            if (contentEpoch.load() != content) continue;
            queue.clear(); view.picture.reset(); view.state = State::failed;
            view.message = "Reference picture: " + juce::String (error.what());
        }
    }
}

#if !JUCE_MAC
std::unique_ptr<Decoder> makeDecoder() { throw std::runtime_error ("Reference video is currently supported on macOS only"); }
#endif
}

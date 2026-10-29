#pragma once
#include "StandaloneJob.h"
#include <juce_core/juce_core.h>
#include <functional>
#include <vector>

namespace levelscope::standalone::pcm24
{
inline constexpr const char* version = "pcm24-tpdf-splitmix64-counter-v1";
struct Format
{
    int rate = 0, channels = 0;
    int64_t frames = 0;
    uint32_t mask = 0;
    uint64_t dataBytes = 0, fileBytes = 0, dataOffset = 0;
    bool rf64 = false;
};
// Checked sizes include the entire RIFF header, not just its data chunk.
bool describe (int rate, int channels, int64_t frames, Format&, juce::String& error, uint32_t mask = 0);
std::vector<uint8_t> header (const Format&);
bool inspect (const juce::File&, Format&, juce::String& error, JobControl&);
// Bounded PCM24 RF64 input compatibility. Never use this in export verification.
// Notices describe tolerated input metadata; no source bytes are repaired/written.
bool inspectRf64Input (const juce::File&, Format&, juce::String& error, JobControl&, juce::String& notice);
uint64_t mix (uint64_t) noexcept;
double dither (uint64_t seed, uint64_t frame, int channel) noexcept;
bool quantize (float sample, double noiseInLsb, int32_t& code) noexcept;

struct Hooks
{
    std::function<bool (const char*)> fail;
    std::function<void (const char*)> reached;
    std::function<int64_t (const juce::File&)> availableBytes;
};
void checkpoint (const Hooks&, const char*, JobControl&);
juce::String hashFile (const juce::File&, JobControl&);

// Exclusive owned output; no destructor-time header finalization or implicit clipping.
class Writer
{
public:
    Writer (const juce::File&, const Format&, const Hooks&, JobControl&);
    ~Writer();
    void write (const float* const*, const std::vector<int>& outputToInput,
                int count, uint64_t seed, JobControl&);
    void finish (JobControl&);
    Writer (const Writer&) = delete;
    Writer& operator= (const Writer&) = delete;
private:
    void bytes (const void*, size_t);
    int fd = -1;
    Format format;
    Hooks hooks;
    uint64_t next = 0;
    std::vector<uint8_t> packed;
};
}

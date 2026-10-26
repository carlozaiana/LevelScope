#include "Pcm24Wave.h"
#include <juce_cryptography/juce_cryptography.h>
#include <array>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <fcntl.h>
#include <unistd.h>
#include <cerrno>

namespace levelscope::standalone::pcm24
{
namespace
{
void require (bool ok, const char* message) { if (! ok) throw std::runtime_error (message); }
void put (std::vector<uint8_t>& b, uint64_t value, int count)
{ for (int n = 0; n < count; ++n) b.push_back (static_cast<uint8_t> (value >> (n * 8))); }
void id (std::vector<uint8_t>& b, const char* text) { b.insert (b.end(), text, text + 4); }
uint64_t number (const uint8_t* b, int size)
{ uint64_t v = 0; for (int i = 0; i < size; ++i) v |= uint64_t (b[i]) << (8 * i); return v; }
constexpr uint8_t pcmGuid[] = { 1,0,0,0,0,0,16,0,128,0,0,170,0,56,155,113 };
}
bool describe (int rate, int channels, int64_t frames, Format& f, juce::String& error)
{
    error.clear(); f = {};
    if ((rate != 44100 && rate != 48000 && rate != 96000 && rate != 192000)
        || (channels != 2 && channels != 6) || frames <= 0 || frames > int64_t (rate) * 10800)
    { error = "PCM24 export supports whole stereo/5.1 programmes up to three hours at 44.1/48/96/192 kHz."; return false; }
    f.rate = rate; f.channels = channels; f.frames = frames; f.mask = channels == 2 ? 3u : 0x3fu;
    f.dataBytes = uint64_t (frames) * static_cast<unsigned> (channels) * 3u;
    f.rf64 = f.dataBytes + 60u > 0xfffffffeu;
    f.dataOffset = f.rf64 ? 104u : 68u; f.fileBytes = f.dataOffset + f.dataBytes;
    return true;
}
std::vector<uint8_t> header (const Format& f)
{
    Format expected; juce::String error;
    require (describe (f.rate, f.channels, f.frames, expected, error) && expected.fileBytes == f.fileBytes
             && expected.rf64 == f.rf64 && expected.mask == f.mask
             && expected.dataBytes == f.dataBytes && expected.dataOffset == f.dataOffset, "Invalid PCM24 format");
    std::vector<uint8_t> b;
    id (b, f.rf64 ? "RF64" : "RIFF"); put (b, f.rf64 ? 0xffffffffu : f.fileBytes - 8, 4); id (b, "WAVE");
    if (f.rf64)
    { id (b, "ds64"); put (b, 28, 4); put (b, f.fileBytes - 8, 8); put (b, f.dataBytes, 8); put (b, f.frames, 8); put (b, 0, 4); }
    id (b, "fmt "); put (b, 40, 4); put (b, 0xfffe, 2); put (b, f.channels, 2);
    put (b, f.rate, 4); put (b, uint64_t (f.rate) * f.channels * 3, 4); put (b, f.channels * 3, 2);
    put (b, 24, 2); put (b, 22, 2); put (b, 24, 2); put (b, f.mask, 4);
    b.insert (b.end(), std::begin (pcmGuid), std::end (pcmGuid));
    id (b, "data"); put (b, f.rf64 ? 0xffffffffu : f.dataBytes, 4); return b;
}
bool inspect (const juce::File& file, Format& f, juce::String& error, JobControl& control)
{
    f = {}; error.clear();
    try
    {
        juce::FileInputStream in (file); const auto total = in.getTotalLength();
        require (in.openedOk() && total >= 68, "Missing or truncated PCM24 WAVE header");
        uint8_t top[12]; require (in.read (top, 12) == 12, "Truncated RIFF header");
        const bool rf = std::memcmp (top, "RF64", 4) == 0;
        require ((rf || std::memcmp (top, "RIFF", 4) == 0) && std::memcmp (top + 8, "WAVE", 4) == 0, "Unsupported WAVE container");
        require (rf ? number (top + 4, 4) == 0xffffffffu : number (top + 4, 4) == uint64_t (total - 8), "RIFF size mismatch");
        bool fmt = false, ds = false, data = false;
        uint64_t size64 = 0, samples64 = 0;
        while (in.getPosition() < total)
        {
            if (control.isCancelled()) throw std::runtime_error ("Cancelled");
            uint8_t ch[8]; require (in.read (ch, 8) == 8, "Truncated chunk header");
            uint64_t size = number (ch + 4, 4), start = static_cast<uint64_t> (in.getPosition());
            if (rf && !ds) require (start == 20 && std::memcmp (ch, "ds64", 4) == 0, "RF64 requires first ds64 chunk");
            if (std::memcmp (ch, "ds64", 4) == 0)
            {
                require (rf && !ds && size == 28, "Unsupported or duplicate ds64 chunk");
                uint8_t b[28]; require (in.read (b, 28) == 28, "Truncated ds64");
                require (number (b, 8) == uint64_t (total - 8) && number (b + 24, 4) == 0, "Invalid ds64 size/table");
                size64 = number (b + 8, 8); samples64 = number (b + 16, 8); ds = true;
            }
            else if (std::memcmp (ch, "fmt ", 4) == 0)
            {
                require (!fmt && !data && size == 40, "Expected one extensible PCM24 format chunk");
                uint8_t b[40]; require (in.read (b, 40) == 40, "Truncated format chunk");
                f.channels = static_cast<int> (number (b + 2, 2)); f.rate = static_cast<int> (number (b + 4, 4));
                f.mask = static_cast<uint32_t> (number (b + 20, 4));
                require (number (b, 2) == 0xfffe && number (b + 14, 2) == 24 && number (b + 16, 2) == 22
                    && number (b + 18, 2) == 24 && std::memcmp (b + 24, pcmGuid, 16) == 0,
                    "Expected 24-valid-bit packed PCM subtype");
                require ((f.channels == 2 || f.channels == 6) && f.mask == (f.channels == 2 ? 3u : 0x3fu)
                    && number (b + 12, 2) == uint64_t (f.channels * 3)
                    && number (b + 8, 4) == uint64_t (f.rate) * f.channels * 3, "Invalid rate/alignment/channel mask");
                fmt = true;
            }
            else if (std::memcmp (ch, "data", 4) == 0)
            {
                require (fmt && !data && (!rf || (ds && size == 0xffffffffu)), "Invalid data chunk");
                if (rf) size = size64;
                require (size > 0 && size % uint64_t (f.channels * 3) == 0, "Incomplete PCM frame");
                f.frames = static_cast<int64_t> (size / uint64_t (f.channels * 3));
                require (!rf || samples64 == uint64_t (f.frames), "RF64 sample count mismatch");
                f.dataOffset = start; f.dataBytes = size; data = true;
            }
            require (start <= uint64_t (total) && size <= uint64_t (total) - start, "Chunk exceeds physical file");
            const auto next = start + size + (size & 1u);
            require (next <= uint64_t (total) && in.setPosition (static_cast<int64_t> (next)), "Invalid chunk padding/seek");
        }
        require (fmt && data && (!rf || ds), "Required PCM chunk missing");
        Format bounded; juce::String reason;
        require (describe (f.rate, f.channels, f.frames, bounded, reason), "Unsupported PCM24 rate/layout/duration");
        f.rf64 = rf; f.fileBytes = static_cast<uint64_t> (total); return true;
    }
    catch (const std::exception& e) { error = e.what(); return false; }
}
uint64_t mix (uint64_t v) noexcept
{
    v = (v ^ (v >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    v = (v ^ (v >> 27)) * UINT64_C(0x94d049bb133111eb); return v ^ (v >> 31);
}
double dither (uint64_t seed, uint64_t frame, int channel) noexcept
{
    // Counter-addressed independent channel streams: no block-size or thread state.
    const auto key = seed ^ mix (uint64_t (channel + 1) * UINT64_C(0x9e3779b97f4a7c15));
    const auto a = mix (key + frame * 2 + UINT64_C(0x9e3779b97f4a7c15));
    const auto b = mix (key + frame * 2 + 1 + UINT64_C(0x9e3779b97f4a7c15));
    return static_cast<double> (a >> 11) * 0x1p-53 - static_cast<double> (b >> 11) * 0x1p-53;
}
bool quantize (float sample, double noise, int32_t& code) noexcept
{
    if (!std::isfinite (sample) || sample < -1 || sample >= 1 || !std::isfinite (noise) || std::abs (noise) > 1) return false;
    const double value = static_cast<double> (sample) * 8388608.0 + noise;
    const double lower = std::floor (value), fraction = value - lower;
    const auto base = static_cast<int64_t> (lower);
    const auto rounded = base + (fraction > .5 || (fraction == .5 && (base & 1)) ? 1 : 0);
    if (rounded < -8388608 || rounded > 8388607) return false;
    code = static_cast<int32_t> (rounded); return true;
}
void checkpoint (const Hooks& hooks, const char* name, JobControl& control)
{
    if (hooks.reached) hooks.reached (name);
    require (!control.isCancelled(), "Cancelled");
    if (hooks.fail && hooks.fail (name)) throw std::runtime_error (std::string ("Injected failure: ") + name);
}
juce::String hashFile (const juce::File& file, JobControl& control)
{
    class Input final : public juce::InputStream
    {
    public:
        Input (const juce::File& f, JobControl& c)
            : raw(f), buffered(raw,65536), length(raw.getTotalLength()), job(c)
        { require(raw.openedOk() && length >= 0, "Cannot open file for hashing"); }
        int64_t getTotalLength() override { return length; }
        int64_t getPosition() override { return buffered.getPosition(); }
        bool setPosition (int64_t p) override { return buffered.setPosition(p); }
        bool isExhausted() override { return getPosition() >= length; }
        int read (void* dest, int n) override
        {
            require(!job.isCancelled(), "Cancelled");
            const auto wanted=static_cast<int>(std::min<int64_t>(n,length-getPosition()));
            const auto got=buffered.read(dest,wanted);
            require(got==wanted, "Short hash read");return got;
        }
        juce::FileInputStream raw;
        juce::BufferedInputStream buffered;
        const int64_t length;
        JobControl& job;
    } input(file,control);
    // SHA256 requests small blocks. Buffer them without repeated filesystem size
    // queries, while still rejecting a short read or a changed final length.
    const auto hash=juce::SHA256(input).toHexString();
    require(input.getPosition()==input.length && input.raw.getStatus().wasOk()
        && file.getSize()==input.length,"File length/read status changed during hashing");
    return hash;
}

Writer::Writer (const juce::File& file, const Format& f, const Hooks& h, JobControl& control) : format (f), hooks (h)
{
    const auto bytesToWrite = header (f); checkpoint (hooks, "audio-create", control);
    fd = ::open (file.getFullPathName().toRawUTF8(), O_WRONLY | O_CREAT | O_EXCL | O_NOFOLLOW, 0600);
    require (fd >= 0, "Cannot exclusively create PCM output");
    try { checkpoint (hooks, "audio-header", control); bytes (bytesToWrite.data(), bytesToWrite.size()); }
    catch (...) { ::close (fd); fd = -1; throw; }
    packed.reserve (4096 * static_cast<size_t> (f.channels) * 3);
}
Writer::~Writer() { if (fd >= 0) ::close (fd); }
void Writer::bytes (const void* data, size_t count)
{
    const auto* p = static_cast<const uint8_t*> (data);
    while (count)
    {
        const auto n = ::write (fd, p, count);
        if (n < 0 && errno == EINTR) continue;
        require (n > 0, "PCM output write failed (possibly insufficient disk space)");
        p += n; count -= static_cast<size_t> (n);
    }
}
void Writer::write (const float* const* audio, const std::vector<int>& map, int count, uint64_t seed, JobControl& control)
{
    require (fd >= 0 && count > 0 && count <= 4096 && map.size() == static_cast<size_t> (format.channels)
             && next + static_cast<unsigned> (count) <= uint64_t (format.frames), "Invalid PCM output span");
    checkpoint (hooks, "audio-write", control); packed.clear();
    for (int i = 0; i < count; ++i) for (int ch = 0; ch < format.channels; ++ch)
    {
        require (map[static_cast<size_t> (ch)] >= 0 && map[static_cast<size_t> (ch)] < format.channels, "Invalid output channel mapping");
        int32_t code = 0;
        require (quantize (audio[map[static_cast<size_t> (ch)]][i], dither (seed, next + static_cast<unsigned> (i), ch), code),
                 "PCM24 sample overflow/non-finite value; processing gain retained, no limiting or attenuation");
        put (packed, static_cast<uint32_t> (code), 3);
    }
    bytes (packed.data(), packed.size()); next += static_cast<unsigned> (count);
}
void Writer::finish (JobControl& control)
{
    require (fd >= 0 && next == uint64_t (format.frames), "Incomplete PCM output");
    checkpoint (hooks, "audio-flush", control); require (::fsync (fd) == 0, "PCM output flush failed");
    checkpoint (hooks, "audio-close", control); const int closing = fd; fd = -1;
    require (::close (closing) == 0, "PCM output close failed");
}
}

#include "StandaloneMedia.h"
#include "StandaloneMeasurementLayout.h"
#include <juce_cryptography/juce_cryptography.h>
#include <cmath>

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-MEDIA-VALIDATION]
void registerSourceFormats (juce::AudioFormatManager& formats)
{
    formats.registerFormat (new juce::WavAudioFormat(), true);
    formats.registerFormat (new juce::AiffAudioFormat(), false);
    formats.registerFormat (new juce::FlacAudioFormat(), false);
}

namespace
{
// Read only bounded metadata, never the PCM payload. JUCE's native FLAC reader
// does not expose Vorbis channel-mask overrides. Skipping a large unrelated
// comment is O(1), and cancellation is checked at each block/comment boundary.
juce::String flacChannelMask (const juce::File& file, JobControl& control, std::optional<uint32_t>& mask)
{
    juce::FileInputStream input (file);
    if (! input.setPosition (4)) return "Cannot read FLAC metadata";
    bool last = false;
    while (! last)
    {
        if (control.isCancelled()) return "Import cancelled";
        uint8_t header[4];
        if (input.read (header, 4) != 4) return "Truncated FLAC metadata header";
        last = (header[0] & 128) != 0;
        const auto length = (uint32_t (header[1]) << 16) | (uint32_t (header[2]) << 8) | header[3];
        const auto end = input.getPosition() + length;
        if (end > input.getTotalLength()) return "Truncated FLAC metadata block";
        if ((header[0] & 127) == 4)
        {
            const auto readLength = [&] (uint32_t& n)
            {
                if (end - input.getPosition() < 4) return false;
                n = static_cast<uint32_t> (input.readInt()); return true;
            };
            uint32_t vendor = 0, entries = 0;
            if (! readLength (vendor) || vendor > end - input.getPosition()
                || ! input.setPosition (input.getPosition() + vendor) || ! readLength (entries)
                || entries > (end - input.getPosition()) / 4) return "Malformed FLAC Vorbis comments";
            for (uint32_t i = 0; i < entries; ++i)
            {
                if (control.isCancelled()) return "Import cancelled";
                uint32_t size = 0;
                if (! readLength (size) || size > end - input.getPosition()) return "Malformed FLAC comment length";
                const auto next = input.getPosition() + size;
                char prefix[129] {};
                const int count = static_cast<int> (std::min<uint32_t> (128, size));
                if (input.read (prefix, count) != count) return "Truncated FLAC comment";
                const auto text = juce::String::fromUTF8 (prefix, count);
                if (text.upToFirstOccurrenceOf ("=", false, false).equalsIgnoreCase ("WAVEFORMATEXTENSIBLE_CHANNEL_MASK"))
                {
                    if (mask || size > 128 || text.getNumBytesAsUTF8() != size || ! text.contains ("="))
                        return "Duplicate or malformed FLAC channel mask";
                    const auto value = text.fromFirstOccurrenceOf ("=", false, false);
                    const auto digits = value.substring (2);
                    if (! value.startsWithIgnoreCase ("0x") || digits.isEmpty() || digits.length() > 8
                        || ! digits.containsOnly ("0123456789abcdefABCDEF")) return "Malformed FLAC channel mask";
                    mask = static_cast<uint32_t> (digits.getHexValue64());
                }
                if (! input.setPosition (next)) return "Cannot skip FLAC comment";
            }
            if (input.getPosition() != end) return "Malformed FLAC comment block size";
        }
        if (! input.setPosition (end)) return "Cannot skip FLAC metadata block";
    }
    return {};
}

juce::String interpretEightChannels (const juce::File& file, MediaInfo& media, JobControl& control)
{
    media.suggestedLayout = {}; media.suggestedLayout.count = 8;
    juce::FileInputStream container (file); char magic[4] {};
    container.read (magic, 4);
    if (std::memcmp (magic, "RF64", 4) == 0 || std::memcmp (magic, "BW64", 4) == 0)
        media.measurementRestriction = "M7 does not qualify RF64/BW64 container extensions; use supported ordinary PCM containers.";
    std::optional<uint32_t> mask;
    const bool flac = media.format.containsIgnoreCase ("flac");
    if (flac)
    {
        const auto error = flacChannelMask (file, control, mask);
        if (error.isNotEmpty()) return error;
        if (! mask)
        {
            media.suggestedLayout = systemILayout();
            media.layoutProvenance = "Native FLAC default: L R C LFE Lrs Rrs Lss Rss (RFC 9639). Confirm roles.";
            return {};
        }
    }
    else if (media.format.containsIgnoreCase ("wav") && media.tags.containsKey ("ChannelMask"))
        mask = static_cast<uint32_t> (media.tags["ChannelMask"].getLargeIntValue());
    if (mask && *mask != 0)
    {
        auto bits = *mask; int count = 0;
        while (bits) { count += static_cast<int> (bits & 1); bits >>= 1; }
        if (count != 8) return "Channel mask does not describe eight file channels";
        if (*mask == 0x63f)
        {
            media.suggestedLayout = systemILayout();
            media.layoutProvenance = (flac ? "FLAC override" : "WAVEFORMATEXTENSIBLE")
                + juce::String (" mask 0x63f: L R C LFE Lrs Rrs Lss Rss. Confirm roles.");
        }
        else media.layoutProvenance = "Unsupported channel mask 0x" + juce::String::toHexString (*mask)
            + ". Explicit manual interpretation required; no layout inferred.";
        return {};
    }
    media.layoutProvenance = "Unlabelled eight-channel PCM: assign roles or an explicit order preset, then confirm.";
    if (media.format.containsIgnoreCase ("aiff"))
    {
        juce::FileInputStream input (file); input.setPosition (12);
        while (input.getTotalLength() - input.getPosition() >= 8)
        {
            if (control.isCancelled()) return "Import cancelled";
            char id[4]; input.read (id, 4);
            const auto size = static_cast<uint32_t> (input.readIntBigEndian());
            if (std::memcmp (id, "CHAN", 4) == 0)
                media.layoutProvenance = "AIFF CHAN metadata is not interpreted. Explicit manual interpretation required.";
            const auto next = input.getPosition() + size + (size & 1);
            if (next > input.getTotalLength() || ! input.setPosition (next)) break;
        }
    }
    return {};
}

// JUCE PCM readers may report success after padding a short stream read. Track
// actual reads during decode so a truncated payload cannot validate as silence.
class CheckedInput final : public juce::InputStream
{
public:
    CheckedInput (const juce::File& file, JobControl& job) : source (file), control (job) {}
    juce::int64 getTotalLength() override { return source.getTotalLength(); }
    juce::int64 getPosition() override { return source.getPosition(); }
    bool setPosition (juce::int64 position) override { return source.setPosition (position); }
    bool isExhausted() override { return control.isCancelled() || source.isExhausted(); }
    int read (void* destination, int bytes) override
    {
        if (control.isCancelled()) return 0;
        const auto count = source.read (destination, bytes);
        if (checking && count != bytes) shortRead = true;
        return count;
    }
    bool checking = false, shortRead = false;
private:
    juce::FileInputStream source;
    JobControl& control;
};

class CancellableInput final : public juce::InputStream
{
public:
    CancellableInput (juce::FileInputStream& input, JobControl& job)
        : source (input, 65536), control (job), length (input.getTotalLength()) {}
    juce::int64 getTotalLength() override { return length; }
    juce::int64 getPosition() override { return source.getPosition(); }
    bool setPosition (juce::int64 position) override { return source.setPosition (position); }
    bool isExhausted() override { return control.isCancelled() || source.isExhausted(); }
    int read (void* destination, int bytes) override
    {
        if (control.isCancelled()) return 0;
        const auto count = source.read (destination, bytes);
        control.setProgress (0.8f + 0.2f * static_cast<float> (
            static_cast<double> (getPosition()) / juce::jmax<juce::int64> (1, getTotalLength())));
        return count;
    }
private:
    juce::BufferedInputStream source;
    JobControl& control;
    const juce::int64 length;
};

// FLAC readers also pad on missing frames. Verify the STREAMINFO PCM MD5 so
// decoder recovery/silence cannot hide corruption. JUCE 7 exposes no CRC status.
class DecodedPcmStream final : public juce::InputStream
{
public:
    DecodedPcmStream (juce::AudioFormatReader& input, JobControl& job)
        : reader (input), control (job), bytesPerSample ((static_cast<int> (input.bitsPerSample) + 7) / 8) {}
    juce::int64 getTotalLength() override { return reader.lengthInSamples * reader.numChannels * bytesPerSample; }
    juce::int64 getPosition() override { return position; }
    bool setPosition (juce::int64 next) override { return next == position; }
    bool isExhausted() override { return control.isCancelled() || position == getTotalLength() || failed; }
    int read (void* destination, int requested) override
    {
        int copied = 0;
        auto* output = static_cast<uint8_t*> (destination);
        while (copied < requested && ! isExhausted())
        {
            if (bufferOffset == bufferSize && ! refill()) break;
            const auto count = juce::jmin (requested - copied, bufferSize - bufferOffset);
            std::memcpy (output + copied, packed.data() + bufferOffset, static_cast<size_t> (count));
            copied += count; bufferOffset += count; position += count;
        }
        return copied;
    }
    bool failed = false;
private:
    bool refill()
    {
        const auto count = static_cast<int> (juce::jmin<juce::int64> (1024, reader.lengthInSamples - frame));
        std::array<int*, maxChannels> pointers;
        for (unsigned int channel = 0; channel < reader.numChannels; ++channel) pointers[channel] = pcm[channel].data();
        if (count <= 0 || ! reader.read (pointers.data(), static_cast<int> (reader.numChannels), frame, count, false))
        { failed = true; return false; }
        int out = 0;
        for (int sample = 0; sample < count; ++sample)
            for (unsigned int channel = 0; channel < reader.numChannels; ++channel)
            {
                const auto value = static_cast<uint32_t> (pcm[channel][sample] >> (32 - reader.bitsPerSample));
                for (int byte = 0; byte < bytesPerSample; ++byte)
                    packed[static_cast<size_t> (out++)] = static_cast<uint8_t> (value >> (8 * byte));
            }
        frame += count;
        bufferSize = out; bufferOffset = 0;
        return true;
    }
    juce::AudioFormatReader& reader;
    JobControl& control;
    int bytesPerSample, bufferOffset = 0, bufferSize = 0;
    juce::int64 frame = 0, position = 0;
    std::array<std::array<int, 1024>, maxChannels> pcm;
    std::array<uint8_t, maxChannels * 1024 * 4> packed;
};
}

MediaValidationResult validateMedia (const juce::File& file, JobControl& control)
{
    MediaValidationResult result;
    auto fail = [&] (ErrorCode code, const juce::String& message)
    {
        result.error = code;
        result.message = message;
        return result;
    };
    if (control.isCancelled()) return fail (ErrorCode::cancelled, "Import cancelled");
    if (! file.existsAsFile()) return fail (ErrorCode::missingFile, "Source file is missing");
    auto& media = result.media;
    media.bytes = file.getSize();
    media.modified = file.getLastModificationTime();

    juce::AudioFormatManager formats;
    registerSourceFormats (formats);
    // Stream overload probes contents, not the filename extension.
    auto checkedInput = std::make_unique<CheckedInput> (file, control);
    auto* checked = checkedInput.get();
    std::unique_ptr<juce::AudioFormatReader> reader (formats.createReaderFor (std::move (checkedInput)));
    if (control.isCancelled()) return fail (ErrorCode::cancelled, "Import cancelled");
    if (! reader) return fail (ErrorCode::unsupportedMedia, "Cannot decode source as WAV, AIFF or FLAC");
    media.sampleRate = reader->sampleRate;
    media.frames = reader->lengthInSamples;
    media.channels = static_cast<int> (reader->numChannels);
    media.bitsPerSample = static_cast<int> (reader->bitsPerSample);
    media.format = reader->getFormatName();
    checked->checking = ! media.format.containsIgnoreCase ("flac");
    media.tags = reader->metadataValues;
    if (! std::isfinite (media.sampleRate) || media.sampleRate <= 0 || media.frames <= 0)
        return fail (ErrorCode::invalidMetadata, "Source has invalid sample rate or no audio frames");
    if (media.channels < 1 || media.channels > maxChannels)
        return fail (ErrorCode::unsupportedChannels, "Source must contain 1 to 12 non-ambisonic channels");
    if (media.bitsPerSample < 1 || media.bitsPerSample > 32)
        return fail (ErrorCode::unsupportedMedia, "Unsupported source bit depth");

    juce::String flacChecksum;
    if (media.format.containsIgnoreCase ("flac"))
    {
        juce::FileInputStream header (file);
        std::array<uint8_t, 42> bytes {};
        if (header.read (bytes.data(), static_cast<int> (bytes.size())) != static_cast<int> (bytes.size())
            || std::memcmp (bytes.data(), "fLaC", 4) != 0 || (bytes[4] & 0x7f) != 0
            || bytes[5] != 0 || bytes[6] != 0 || bytes[7] != 34)
            return fail (ErrorCode::unsupportedMedia, "FLAC requires native STREAMINFO metadata");
        uint64_t frames = bytes[21] & 0x0f;
        for (int i = 22; i < 26; ++i) frames = (frames << 8) | bytes[i];
        if (frames == 0 || frames > 0xffffffffu || frames != static_cast<uint64_t> (media.frames))
            return fail (ErrorCode::unsupportedMedia, "FLAC frame count is unavailable or exceeds JUCE 7's supported range");
        flacChecksum = juce::String::toHexString (bytes.data() + 26, 16, 0);
        if (flacChecksum == juce::String::repeatedString ("0", 32))
            return fail (ErrorCode::unsupportedMedia, "FLAC requires a PCM checksum for integrity validation");
    }

    const auto layout = reader->getChannelLayout();
    if (layout.getAmbisonicOrder() >= 0)
        return fail (ErrorCode::unsupportedChannels, "Ambisonic sources are outside this scaffold's scope");
    media.suggestedLayout.count = media.channels;
    for (int channel = 0; channel < juce::jmin (layout.size(), media.channels); ++channel)
        media.suggestedLayout.roles[channel] = layout.getTypeOfChannel (channel);
    media.layoutProvenance = "Legacy decoder role suggestions; confirm file-channel interpretation.";
    if (media.channels == 8)
    {
        const auto error = interpretEightChannels (file, media, control);
        if (error.isNotEmpty()) return fail (control.isCancelled() ? ErrorCode::cancelled : ErrorCode::invalidMetadata, error);
    }

    constexpr int blockSize = 4096;
    juce::AudioBuffer<float> block (media.channels, blockSize);
    for (juce::int64 position = 0; position < media.frames;)
    {
        if (control.isCancelled()) return fail (ErrorCode::cancelled, "Import cancelled");
        const auto count = static_cast<int> (juce::jmin<juce::int64> (blockSize, media.frames - position));
        if (! reader->read (block.getArrayOfWritePointers(), media.channels, position, count))
            return fail (ErrorCode::decodeFailed, "Audio decoder reported a read failure");
        if (control.isCancelled()) return fail (ErrorCode::cancelled, "Import cancelled");
        if (checked->shortRead)
            return fail (ErrorCode::decodeFailed, "Truncated audio payload");
        for (int channel = 0; channel < media.channels; ++channel)
            for (int frame = 0; frame < count; ++frame)
                if (! std::isfinite (block.getSample (channel, frame)))
                    return fail (ErrorCode::nonFiniteAudio, "Source contains non-finite audio samples");
        position += count;
        control.setProgress (0.8f * static_cast<float> (static_cast<double> (position) / media.frames));
    }
    if (flacChecksum.isNotEmpty())
    {
        DecodedPcmStream pcm (*reader, control);
        const auto actual = juce::MD5 (pcm).toHexString();
        if (control.isCancelled()) return fail (ErrorCode::cancelled, "Import cancelled");
        if (pcm.failed || actual != flacChecksum)
            return fail (ErrorCode::decodeFailed, "FLAC audio checksum mismatch (corrupt or truncated source)");
    }
    auto input = file.createInputStream();
    if (! input) return fail (ErrorCode::missingFile, "Cannot reopen source for content hashing");
    CancellableInput stream (*input, control);
    media.sha256 = juce::SHA256 (stream).toHexString();
    if (control.isCancelled()) return fail (ErrorCode::cancelled, "Import cancelled");
    if (input->getStatus().failed() || stream.getPosition() != media.bytes)
        return fail (ErrorCode::decodeFailed, "Cannot read the complete source for hashing");
    if (! media.matchesFileStamp (file)) return fail (ErrorCode::sourceChanged, "Source changed during validation; import again");
    result.message = "Source validation succeeded";
    control.setProgress (1.0f);
    return result;
}
// [END LS-STANDALONE-MEDIA-VALIDATION]
}

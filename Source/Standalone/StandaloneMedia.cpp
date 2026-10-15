#include "StandaloneMedia.h"
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
// JUCE PCM readers may report success after padding a short stream read. Track
// actual reads during decode so a truncated payload cannot validate as silence.
class CheckedInput final : public juce::InputStream
{
public:
    CheckedInput (const juce::File& file, JobControl& job) : source (file), control (job) {}
    juce::int64 getTotalLength() override { return source.getTotalLength(); }
    juce::int64 getPosition() override { return source.getPosition(); }
    bool setPosition (juce::int64 position) override { return source.setPosition (position); }
    bool isExhausted() override { return control.cancel.load() || source.isExhausted(); }
    int read (void* destination, int bytes) override
    {
        if (control.cancel.load()) return 0;
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
    CancellableInput (juce::FileInputStream& input, JobControl& job) : source (input), control (job) {}
    juce::int64 getTotalLength() override { return source.getTotalLength(); }
    juce::int64 getPosition() override { return source.getPosition(); }
    bool setPosition (juce::int64 position) override { return source.setPosition (position); }
    bool isExhausted() override { return control.cancel.load() || source.isExhausted(); }
    int read (void* destination, int bytes) override
    {
        if (control.cancel.load()) return 0;
        const auto count = source.read (destination, bytes);
        control.progress.store (0.8f + 0.2f * static_cast<float> (
            static_cast<double> (getPosition()) / juce::jmax<juce::int64> (1, getTotalLength())));
        return count;
    }
private:
    juce::FileInputStream& source;
    JobControl& control;
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
    bool isExhausted() override { return control.cancel.load() || position == getTotalLength() || failed; }
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
    if (control.cancel.load()) return fail (ErrorCode::cancelled, "Import cancelled");
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
    if (control.cancel.load()) return fail (ErrorCode::cancelled, "Import cancelled");
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

    constexpr int blockSize = 4096;
    juce::AudioBuffer<float> block (media.channels, blockSize);
    for (juce::int64 position = 0; position < media.frames;)
    {
        if (control.cancel.load()) return fail (ErrorCode::cancelled, "Import cancelled");
        const auto count = static_cast<int> (juce::jmin<juce::int64> (blockSize, media.frames - position));
        if (! reader->read (block.getArrayOfWritePointers(), media.channels, position, count))
            return fail (ErrorCode::decodeFailed, "Audio decoder reported a read failure");
        if (control.cancel.load()) return fail (ErrorCode::cancelled, "Import cancelled");
        if (checked->shortRead)
            return fail (ErrorCode::decodeFailed, "Truncated audio payload");
        for (int channel = 0; channel < media.channels; ++channel)
            for (int frame = 0; frame < count; ++frame)
                if (! std::isfinite (block.getSample (channel, frame)))
                    return fail (ErrorCode::nonFiniteAudio, "Source contains non-finite audio samples");
        position += count;
        control.progress.store (0.8f * static_cast<float> (static_cast<double> (position) / media.frames));
    }
    if (flacChecksum.isNotEmpty())
    {
        DecodedPcmStream pcm (*reader, control);
        const auto actual = juce::MD5 (pcm).toHexString();
        if (control.cancel.load()) return fail (ErrorCode::cancelled, "Import cancelled");
        if (pcm.failed || actual != flacChecksum)
            return fail (ErrorCode::decodeFailed, "FLAC audio checksum mismatch (corrupt or truncated source)");
    }
    auto input = file.createInputStream();
    if (! input) return fail (ErrorCode::missingFile, "Cannot reopen source for content hashing");
    CancellableInput stream (*input, control);
    media.sha256 = juce::SHA256 (stream).toHexString();
    if (control.cancel.load()) return fail (ErrorCode::cancelled, "Import cancelled");
    if (input->getStatus().failed() || stream.getPosition() != media.bytes)
        return fail (ErrorCode::decodeFailed, "Cannot read the complete source for hashing");
    if (! media.matchesFileStamp (file)) return fail (ErrorCode::sourceChanged, "Source changed during validation; import again");
    result.message = "Source validation succeeded";
    control.progress.store (1.0f);
    return result;
}
// [END LS-STANDALONE-MEDIA-VALIDATION]
}

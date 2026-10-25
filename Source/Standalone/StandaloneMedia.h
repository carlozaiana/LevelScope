#pragma once

#include "StandaloneChannelLayout.h"
#include "StandaloneJob.h"
#include <juce_audio_formats/juce_audio_formats.h>

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-MEDIA-CONTRACT]
struct MediaInfo
{
    double sampleRate = 0.0;
    juce::int64 frames = 0, bytes = 0;
    int channels = 0, bitsPerSample = 0;
    juce::Time modified;
    juce::String format, sha256;
    juce::String layoutProvenance;
    juce::String measurementRestriction; // Broader Source audition may remain available.
    juce::StringPairArray tags;
    ChannelLayout suggestedLayout;

    double durationSeconds() const { return sampleRate > 0 ? static_cast<double> (frames) / sampleRate : 0.0; }
    bool matchesFileStamp (const juce::File& file) const
    {
        return file.existsAsFile() && file.getSize() == bytes && file.getLastModificationTime() == modified;
    }
};

struct MediaValidationResult
{
    MediaInfo media;
    ErrorCode error = ErrorCode::none;
    juce::String message;
    bool succeeded() const { return error == ErrorCode::none; }
};

// Deliberate initial matrix: WAV, AIFF, FLAC. No OS-dependent codecs.
void registerSourceFormats (juce::AudioFormatManager&);
MediaValidationResult validateMedia (const juce::File&, JobControl&);
// [END LS-STANDALONE-MEDIA-CONTRACT]
}

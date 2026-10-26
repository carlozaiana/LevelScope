#pragma once
#include "StandaloneMedia.h"
#include <memory>
#include <vector>

namespace levelscope::standalone
{
// Each member supplies one channel at programme sample zero. Locators and slot IDs
// are persistence concerns; neither is evidence of content or alignment.
struct ProgrammeMember
{
    juce::String id;
    juce::File file;
    MediaInfo media;
    ChannelRole role = juce::AudioChannelSet::unknown;
    juce::int64 offsetFrames = 0;
};
struct ProgrammeAssembly
{
    static constexpr int version = 2;
    std::vector<ProgrammeMember> members; // canonical programme order, not picker order
    juce::String alignmentConfirmation;
    juce::String definitionIdentity() const;
    juce::String validate (bool requireConfirmation = true) const;
    bool matchesFileStamps() const;
    ChannelLayout layout() const;
    MediaInfo summary() const;
    // Explicit user action only; never inferred from filenames, length or hashes.
    bool confirmAlignment();
    bool orderByRoles(); // explicit assignments only; no filename inference
};
struct ProgrammeInspection
{
    ProgrammeAssembly group;
    ErrorCode error = ErrorCode::none;
    juce::String message;
};
ProgrammeInspection inspectProgrammeFiles (const std::vector<juce::File>&, JobControl&);
bool samePhysicalFile (const juce::File&, const juce::File&);
std::vector<ChannelRole> assemblyRoles (int channels);

// Worker-only adapter. All reads stay inside the declared programme; a failed
// member read fails the entire block. Only compose() permits the explicit M11 sum.
class ProgrammeReader
{
public:
    static std::unique_ptr<ProgrammeReader> open (const std::vector<juce::File>&, const MediaInfo&);
    static std::unique_ptr<ProgrammeReader> compose (std::unique_ptr<ProgrammeReader> me,
        std::unique_ptr<ProgrammeReader> dialogue, const std::vector<int>& meMap, const std::vector<int>& dialogueMap);
    bool read (float* const*, int channels, juce::int64 start, int count);
    double sampleRate = 0;
    juce::int64 lengthInSamples = 0;
    unsigned int numChannels = 0;
    bool hasFlac = false;
private:
    std::vector<std::unique_ptr<juce::AudioFormatReader>> readers;
    std::unique_ptr<ProgrammeReader> meReader, dialogueReader;
    std::vector<int> meMapping, dialogueMapping;
    juce::AudioBuffer<float> meBuffer, dialogueBuffer;
};
}

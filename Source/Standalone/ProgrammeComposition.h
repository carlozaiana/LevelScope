#pragma once
#include "ProgrammeSource.h"
#include <optional>

namespace levelscope::standalone
{
// Original component media remains accessible; a component is not a rendered mix.
// Labels are document/report context. Only the audible definition enters identity.
struct ProgrammeComponent
{
    juce::String id, purpose, language, versionLabel;
    juce::File file;
    MediaInfo media;
    ChannelLayout layout;
    std::optional<ProgrammeAssembly> assembly;
    juce::int64 offsetFrames = 0;
    juce::String validate() const;
    juce::String audioIdentity() const;
    std::vector<ProgrammeMember> files() const;
    bool matchesFileStamps() const;
};

enum class CompositionPreset { stereoStereo = 0, surroundMonoCentre = 1, surroundSurround = 2 };
struct ProgrammeComposition
{
    static constexpr int schemaVersion = 1;
    static constexpr const char* engineVersion = "aligned-two-component-sum-v1";
    CompositionPreset preset = CompositionPreset::stereoStereo;
    ProgrammeComponent me, dialogue;
    juce::String confirmation;
    juce::String validate (bool requireConfirmation = true) const;
    juce::String definitionIdentity() const;
    bool confirmComplementaryBalancedAlignment();
    bool matchesFileStamps() const;
    ChannelLayout layout() const;
    MediaInfo summary() const;
    // Output channel -> input channel, -1 only for explicitly absent contribution.
    std::vector<int> map (bool dialogueComponent) const;
    std::vector<ProgrammeMember> files() const;
};
}

#pragma once
#include "ProgrammeSource.h"

namespace levelscope::standalone
{
struct ChannelRoleSuggestion
{
    ChannelLayout order;
    juce::String explanation;
};

class RememberedChannelRoles
{
public:
    explicit RememberedChannelRoles(juce::File directory = {});
    ChannelRoleSuggestion suggest(const MediaInfo&) const;
    ChannelRoleSuggestion suggest(const ProgrammeAssembly&) const;
    juce::Result remember(const ChannelLayout&) const;
private:
    juce::File directory;
    std::optional<ChannelLayout> read(int count) const;
    juce::File file(int count) const;
};
}

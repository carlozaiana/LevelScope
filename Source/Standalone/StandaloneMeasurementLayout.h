#pragma once

#include "StandaloneChannelLayout.h"
#include "Core/Measurement/OfflineMeasurement.h"

namespace levelscope::standalone
{
inline constexpr const char* layoutMappingVersion = "ordered-pcm-roles-m7.1";
inline constexpr const char* systemICapability = "pcm-7.1-bs2051-system-i-v1";

// These are FILE orders. Never sort or reorder PCM to match AudioChannelSet.
inline ChannelLayout systemILayout (bool sideBeforeRear = false)
{
    using C = juce::AudioChannelSet;
    ChannelLayout l; l.count = 8;
    l.roles = { C::left, C::right, C::centre, C::LFE,
        sideBeforeRear ? C::leftSurroundSide : C::leftSurroundRear,
        sideBeforeRear ? C::rightSurroundSide : C::rightSurroundRear,
        sideBeforeRear ? C::leftSurroundRear : C::leftSurroundSide,
        sideBeforeRear ? C::rightSurroundRear : C::rightSurroundSide };
    return l;
}

inline bool isSystemI (const ChannelLayout& l)
{
    if (l.count != 8 || ! l.isValid()) return false;
    const auto expected = systemILayout();
    for (int i = 0; i < l.count; ++i)
        if (std::find (expected.roles.begin(), expected.roles.begin() + 8, l.roles[i]) == expected.roles.begin() + 8) return false;
    return true;
}

inline ChannelRole systemIRoleForAlias (juce::String text)
{
    using C = juce::AudioChannelSet;
    text = juce::StringArray::fromTokens (text.trim(), " \t\r\n", "").joinIntoString (" ").toLowerCase();
    // StringArray retains empty tokens: collapse whitespace explicitly.
    while (text.contains ("  ")) text = text.replace ("  ", " ");
    text = text.replace (juce::String::charToString (0x2212), "-");
    const std::pair<ChannelRole, const char*> aliases[] {
        { C::left, "l|fl|left|front left|m+030" },
        { C::right, "r|fr|right|front right|m-030" },
        { C::centre, "c|fc|centre|center|front centre|front center|m+000" },
        { C::LFE, "lfe|lfe1|low frequency effects" },
        { C::leftSurroundSide, "lss|sl|side left|left side surround|m+090" },
        { C::rightSurroundSide, "rss|sr|side right|right side surround|m-090" },
        { C::leftSurroundRear, "lrs|bl|back left|rear left|left rear surround|m+135" },
        { C::rightSurroundRear, "rrs|br|back right|rear right|right rear surround|m-135" }
    };
    for (const auto& entry : aliases)
        if (juce::StringArray::fromTokens (entry.second, "|", "").contains (text)) return entry.first;
    return C::unknown;
}

inline juce::String monitoringCoverage (const ChannelLayout& l, const OutputRouting& routing)
{
    int assigned = 0;
    for (int i = 0; i < l.count; ++i) if (routing.fileToDevice[i] >= 0) ++assigned;
    if (assigned == 0) return "All monitoring routes muted. No downmix.";
    if (assigned < l.count) return "Partial monitoring: " + juce::String (assigned) + " of " + juce::String (l.count) + " channels. No downmix.";
    return "All channels explicitly routed. Physical speaker mapping is not verified.";
}
}

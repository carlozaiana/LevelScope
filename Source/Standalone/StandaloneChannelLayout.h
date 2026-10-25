#pragma once

#include <array>
#include <juce_audio_basics/juce_audio_basics.h>

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-CHANNEL-CONTRACT]
constexpr int maxChannels = 12;
using ChannelRole = juce::AudioChannelSet::ChannelType;
inline const std::array<ChannelRole, 25> supportedRoles {
    juce::AudioChannelSet::left, juce::AudioChannelSet::right,
    juce::AudioChannelSet::centre, juce::AudioChannelSet::LFE,
    juce::AudioChannelSet::leftSurround, juce::AudioChannelSet::rightSurround,
    juce::AudioChannelSet::leftSurroundSide, juce::AudioChannelSet::rightSurroundSide,
    juce::AudioChannelSet::leftSurroundRear, juce::AudioChannelSet::rightSurroundRear,
    juce::AudioChannelSet::centreSurround,
    juce::AudioChannelSet::topFrontLeft, juce::AudioChannelSet::topFrontRight,
    juce::AudioChannelSet::topRearLeft, juce::AudioChannelSet::topRearRight,
    juce::AudioChannelSet::topSideLeft, juce::AudioChannelSet::topSideRight,
    juce::AudioChannelSet::leftCentre, juce::AudioChannelSet::rightCentre,
    juce::AudioChannelSet::topMiddle, juce::AudioChannelSet::topFrontCentre,
    juce::AudioChannelSet::topRearCentre, juce::AudioChannelSet::LFE2,
    juce::AudioChannelSet::wideLeft, juce::AudioChannelSet::wideRight
};

struct ChannelLayout
{
    int count = 0;
    // Indexed by ORIGINAL FILE channel, never by AudioChannelSet's sorted order.
    std::array<ChannelRole, maxChannels> roles {};
    bool confirmed = false;
    bool manuallyInterpreted = false;

    bool isValid() const
    {
        if (count < 1 || count > maxChannels)
            return false;
        for (int i = 0; i < count; ++i)
        {
            if (std::find (supportedRoles.begin(), supportedRoles.end(), roles[i]) == supportedRoles.end())
                return false;
            for (int j = 0; j < i; ++j)
                if (roles[i] == roles[j])
                    return false;
        }
        return true;
    }
};

struct OutputRouting
{
    // Physical device channel indices, zero-based. -1 is deliberately muted.
    std::array<int, maxChannels> fileToDevice;
    OutputRouting() { fileToDevice.fill (-1); }

    juce::String validate (int fileChannels, int deviceChannels, bool allowAllMuted = false) const
    {
        if (fileChannels < 1 || fileChannels > maxChannels || deviceChannels < 1)
            return "Select validated media and an output device.";
        bool audible = false;
        for (int i = 0; i < fileChannels; ++i)
        {
            const auto output = fileToDevice[i];
            if (output < -1 || output >= deviceChannels)
                return "An output is outside the selected device's channel range.";
            if (output < 0)
                continue;
            audible = true;
            for (int j = 0; j < i; ++j)
                if (fileToDevice[j] == output)
                    return "Assign each output once; implicit summing is unavailable.";
        }
        return audible || allowAllMuted ? juce::String() : "All channels are muted. Assign an output explicitly.";
    }
};
// [END LS-STANDALONE-CHANNEL-CONTRACT]
}

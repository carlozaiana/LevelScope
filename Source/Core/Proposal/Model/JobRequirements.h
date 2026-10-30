#pragma once
#include "TargetProfile.h"
#include <optional>

namespace levelscope::proposal
{
// A value copied into each session/request. It never changes the catalogue profile.
struct JobRequirements
{
    static constexpr int maximumBytes = 65536;
    int schema = 1, revision = 1;
    juce::String id = "00000000000000000000000000000017";
    juce::String name = "Programme job requirements", provenance;
    juce::String baseId = "ebu-r128-2023-programme-linear-qc", baseVersion = "1.0.0";
    juce::String baseDigest = "ef6f5cc071a8c9d5a347da0dd61fcc46e9f24da0aed787bc9f3ea48cb3db7e66";
    juce::String method = "EBU Tech 3342-2023 whole-programme LRA";
    juce::String applicability = "pre-produced-whole-programme-linear-audio";
    std::optional<int> maximumLraTenths;
    juce::String unavailableReason; // independently recoverable invalid session field

    juce::String validationError() const;
    juce::String digest() const;
    juce::var encode() const;
    static JobRequirements decode(const juce::var&); // throws; no silent defaults
    static JobRequirements fromJson(const juce::String&);
    static JobRequirements unavailable(const juce::String&);
    static std::optional<int> parseCeiling(const juce::String&);
    TargetProfile applyTo(const TargetProfile&) const;
    juce::String text() const;
};
}

#pragma once

#include <JuceHeader.h>

class LevelScopeAudioProcessor;

//==============================================================================
// LoudnessStatsComponent
// Small UI panel showing:
// - Integrated (running, gated)
// - LRA (running)
// - Rolling LRA (30/60/120s selectable)
//==============================================================================

class LoudnessStatsComponent : public juce::Component,
                               private juce::Timer
{
public:
    explicit LoudnessStatsComponent (LevelScopeAudioProcessor& p);
    ~LoudnessStatsComponent() override = default;

    int getPreferredHeight() const noexcept { return 48; }

    void resized() override;

private:
    void timerCallback() override;

    LevelScopeAudioProcessor& processor;

    juce::Label integratedLabel;
    juce::Label lraLabel;
    juce::Label rollingLabel;

    juce::ComboBox rollingWindowBox;

    static juce::String formatLufs (float v);
    static juce::String formatLu (float v);
};
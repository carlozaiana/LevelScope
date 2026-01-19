#pragma once

#include <JuceHeader.h>

class LevelScopeAudioProcessor;

//==============================================================================
// LoudnessStatsComponent
// Small UI panel showing running Integrated (gated) and LRA.
//==============================================================================

class LoudnessStatsComponent : public juce::Component,
                               private juce::Timer
{
public:
    explicit LoudnessStatsComponent (LevelScopeAudioProcessor& p);
    ~LoudnessStatsComponent() override = default;

    void resized() override;

private:
    void timerCallback() override;

    LevelScopeAudioProcessor& processor;

    juce::Label integratedLabel;
    juce::Label lraLabel;

    static juce::String formatLufs (float v);
    static juce::String formatLu (float v);
};
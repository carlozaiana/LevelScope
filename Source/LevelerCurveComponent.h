#pragma once

#include <JuceHeader.h>

class LevelScopeAudioProcessor;

//==============================================================================
// Leveler mapping display:
// measured loudness (LUFS) -> applied gain (dB)
//==============================================================================

class LevelerCurveComponent : public juce::Component,
                              private juce::Timer
{
public:
    explicit LevelerCurveComponent (LevelScopeAudioProcessor& p);
    ~LevelerCurveComponent() override = default;

    void paint (juce::Graphics& g) override;

private:
    void timerCallback() override;

    LevelScopeAudioProcessor& processor;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (LevelerCurveComponent)
};
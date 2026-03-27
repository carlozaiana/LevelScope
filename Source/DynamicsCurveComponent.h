#pragma once

#include <JuceHeader.h>

class LevelScopeAudioProcessor;

//==============================================================================
// Simple low-cost 2D curve display for dynamics modules.
// Stage UI-CURVE-1: downward transfer curve
// Stage UI-CURVE-2: optional upward conceptual curve
//==============================================================================

class DynamicsCurveComponent : public juce::Component,
                               private juce::Timer
{
public:
    enum class CurveKind
    {
        downward,
        upwardConceptual
    };

    DynamicsCurveComponent (LevelScopeAudioProcessor& p, CurveKind k);
    ~DynamicsCurveComponent() override = default;

    void paint (juce::Graphics& g) override;

private:
    void timerCallback() override;

    LevelScopeAudioProcessor& processor;
    CurveKind kind;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (DynamicsCurveComponent)
};
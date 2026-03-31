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

    // [BEGIN UI-CURVE-LVLR-TARGET-DRAG-PUBLIC]
    void mouseMove (const juce::MouseEvent& e) override;
    void mouseExit (const juce::MouseEvent& e) override;
    void mouseDown (const juce::MouseEvent& e) override;
    void mouseDrag (const juce::MouseEvent& e) override;
    void mouseUp   (const juce::MouseEvent& e) override;
    // [END UI-CURVE-LVLR-TARGET-DRAG-PUBLIC]

private:
    void timerCallback() override;

    // [BEGIN UI-CURVE-LVLR-TARGET-DRAG-DECL]
    bool getTargetInteractionGeometry (juce::Rectangle<float>& plotOut,
                                       float& targetXOut,
                                       bool& editableOut) const;
    void updateMouseCursorForTarget (juce::Point<float> pos);

    bool targetDragging = false;
    bool targetGestureActive = false;
    // [END UI-CURVE-LVLR-TARGET-DRAG-DECL]

    LevelScopeAudioProcessor& processor;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (LevelerCurveComponent)
};
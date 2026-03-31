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

    // [BEGIN UI-CURVE-DOWN-DRAG-PUBLIC]
    void mouseMove (const juce::MouseEvent& e) override;
    void mouseExit (const juce::MouseEvent& e) override;
    void mouseDown (const juce::MouseEvent& e) override;
    void mouseDrag (const juce::MouseEvent& e) override;
    void mouseUp   (const juce::MouseEvent& e) override;
    // [END UI-CURVE-DOWN-DRAG-PUBLIC]

private:
    void timerCallback() override;

    // [BEGIN UI-CURVE-DOWN-DRAG-DECL]
    bool getDownwardInteractionGeometry (juce::Rectangle<float>& plotOut,
                                         float& t2XOut,
                                         float& t3XOut) const;
    void updateMouseCursorForThresholds (juce::Point<float> pos);
    void computeOrderedThresholdsWithPush (int changedIndex,
                                           float newValueLufs,
                                           float outVals[4]) const noexcept;
    void applyThresholdValuesDuringDrag (const float newVals[4]);
    void endAllThresholdGestures();

    bool thresholdDragging = false;
    int  activeThresholdIndex = -1; // 2 or 3
    int  hoverThresholdIndex  = -1; // 2 or 3, or -1

    std::array<bool, 4> thresholdGestureActive { { false, false, false, false } };
    static constexpr float kThreshMinGapLu = 0.1f;
    // [END UI-CURVE-DOWN-DRAG-DECL]

    LevelScopeAudioProcessor& processor;
    CurveKind kind;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (DynamicsCurveComponent)
};
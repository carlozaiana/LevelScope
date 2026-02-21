#pragma once

#include <JuceHeader.h>
// [BEGIN MTDM-PANEL-INCLUDE-MEMORY]
#include <memory>
// [END MTDM-PANEL-INCLUDE-MEMORY]
#include "PluginProcessor.h"
#include "LoudnessStatsComponent.h"
#include "VolumeHistoryComponent.h"

//==============================================================================
// [BEGIN MTDM-PANEL-DECL]
//==============================================================================
// Simple UI panel for MTDM testing (Stage UI-2)
//==============================================================================

class GrMeterComponent : public juce::Component
{
public:
    GrMeterComponent() = default;

    void setNameLabel (juce::String s) { name = std::move (s); }
    void setValuesDb (float currentDb, float holdDb) noexcept;

    void paint (juce::Graphics& g) override;

private:
    juce::String name;
    float current = 0.0f;
    float hold    = 0.0f;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (GrMeterComponent)
};

class MtdmControlPanel : public juce::Component,
                         private juce::Timer
{
public:
    explicit MtdmControlPanel (LevelScopeAudioProcessor& p);
    ~MtdmControlPanel() override;

    void paint (juce::Graphics& g) override;
    void resized() override;

private:
    void timerCallback() override;

    // [BEGIN MTDM-PANEL-THRESH-ORDER-DECL]
    void enforceThresholdOrderingFromSlider (int changedIndex);
    void endPushedThresholdGestures();

    bool thresholdCallbacksSuppressed = false;

    bool thresholdSliderDragging = false;
    int  thresholdSliderDraggingIndex = -1;

    // Gestures for pushed neighbors (not the actively dragged slider)
    std::array<bool, 4> pushedGestureActive { { false, false, false, false } };

    static constexpr float minGapLu = 0.1f;
    // [END MTDM-PANEL-THRESH-ORDER-DECL]

    void configureToggle (juce::ToggleButton& b, const juce::String& text);
    void configureSliderForParam (juce::Slider& s, const juce::String& paramID,
                                  juce::Slider::SliderStyle style,
                                  const juce::String& suffix);

    LevelScopeAudioProcessor& processor;
    juce::AudioProcessorValueTreeState& apvts;

    // Top row toggles
    juce::ToggleButton mtdmEnabledButton;
    juce::ToggleButton downEnabledButton;
    juce::ToggleButton limEnabledButton;

    // Threshold sliders
    juce::Label  t0Label, t1Label, t2Label, t3Label;
    juce::Slider t0Slider, t1Slider, t2Slider, t3Slider;

    // GR meters (read-only)
    GrMeterComponent limGrMeter;
    GrMeterComponent downGrMeter;

    // Attachments (must be kept alive)
    using ButtonAttachment = juce::AudioProcessorValueTreeState::ButtonAttachment;
    using SliderAttachment = juce::AudioProcessorValueTreeState::SliderAttachment;

    std::unique_ptr<ButtonAttachment> mtdmEnabledAtt;
    std::unique_ptr<ButtonAttachment> downEnabledAtt;
    std::unique_ptr<ButtonAttachment> limEnabledAtt;

    std::unique_ptr<SliderAttachment> t0Att;
    std::unique_ptr<SliderAttachment> t1Att;
    std::unique_ptr<SliderAttachment> t2Att;
    std::unique_ptr<SliderAttachment> t3Att;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MtdmControlPanel)
};
// [END MTDM-PANEL-DECL]

class LevelScopeAudioProcessorEditor : public juce::AudioProcessorEditor
{
public:
    explicit LevelScopeAudioProcessorEditor (LevelScopeAudioProcessor&);
    ~LevelScopeAudioProcessorEditor() override;

    void paint (juce::Graphics&) override;
    void resized() override;

private:
    LoudnessStatsComponent  statsComponent;
    VolumeHistoryComponent  historyComponent;

    // [BEGIN MTDM-PANEL-EDITOR-MEMBERS]
    MtdmControlPanel mtdmPanel;

    // Draggable horizontal splitter between history and panel
    juce::StretchableLayoutManager layoutHistoryAndPanel;
    juce::StretchableLayoutResizerBar layoutResizerBar;
    // [END MTDM-PANEL-EDITOR-MEMBERS]

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (LevelScopeAudioProcessorEditor)
};
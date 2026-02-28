#pragma once

#include <JuceHeader.h>
// [BEGIN MTDM-PANEL-INCLUDE-MEMORY]
#include <memory>
// [END MTDM-PANEL-INCLUDE-MEMORY]
#include "PluginProcessor.h"
#include "VolumeHistoryComponent.h"

// [BEGIN UI3A-MISSIONCONTROL-DECL]
//==============================================================================
// Mission Control Top Strip (preset + targets/current + MC routing + curve toggles)
//==============================================================================

class MissionControlComponent : public juce::Component,
                                private juce::Timer
{
public:
    MissionControlComponent (LevelScopeAudioProcessor& p, VolumeHistoryComponent& h);
    ~MissionControlComponent() override;

    int getPreferredHeight() const noexcept { return preferredHeightPx; }

    void paint (juce::Graphics& g) override;
    void resized() override;

private:
    void timerCallback() override;

    void loadTargetsFromState();
    void storeTargetToState (const juce::Identifier& key, double value);
    double getTargetFromState (const juce::Identifier& key, double defaultValue) const;

    void updateCurrentReadouts();

    void startSavePreset();
    void startLoadPreset();

    void refreshCurveToggleStatesFromHistory();

    // Small graphic: detector/apply channel squares (always visible)
    class RoutingGraphic : public juce::Component
    {
    public:
        RoutingGraphic (LevelScopeAudioProcessor& p);
        void paint (juce::Graphics& g) override;

    private:
        LevelScopeAudioProcessor& processor;

        juce::String channelLabelForType (juce::AudioChannelSet::ChannelType t) const;
        bool isDetectorChannelActive (juce::AudioChannelSet::ChannelType t,
                                      int mcPolicy, int detChoice, bool lfeInDet) const;
        bool isApplyChannelActive (juce::AudioChannelSet::ChannelType t,
                                   int mcPolicy, int applyChoice, bool lfeInApply) const;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (RoutingGraphic)
    };

    LevelScopeAudioProcessor& processor;
    VolumeHistoryComponent&   history;
    juce::AudioProcessorValueTreeState& apvts;

    // [BEGIN UI3A-MISSIONCONTROL-PREFERRED-HEIGHT]
    const int preferredHeightPx = 140;
    // [END UI3A-MISSIONCONTROL-PREFERRED-HEIGHT]

    // Presets
    juce::TextButton savePresetButton { "Save Preset..." };
    juce::TextButton loadPresetButton { "Load Preset..." };

    // [BEGIN UI3A-MISSIONCONTROL-TARGETS-CURRENT-LABELS]
    // Column headers
    juce::Label hdrILabel, hdrPeakLabel, hdrLraLabel;

    // Row headers
    juce::Label hdrTargetRow, hdrCurrentRow;

    // Targets (editable values)
    juce::Label targetILabel, targetPeakLabel, targetLraLabel;

    // Current (read-only values)
    juce::Label currentILabel, currentPeakLabel, currentLraLabel;
    // [END UI3A-MISSIONCONTROL-TARGETS-CURRENT-LABELS]

    // Policy controls (APVTS)
    juce::ComboBox mcPolicyBox, dialogDetBox, dialogApplyBox;
    juce::ToggleButton lfeDetToggle { "LFE det" };
    juce::ToggleButton lfeApplyToggle { "LFE apply" };

    using ComboAttachment  = juce::AudioProcessorValueTreeState::ComboBoxAttachment;
    using ButtonAttachment = juce::AudioProcessorValueTreeState::ButtonAttachment;

    std::unique_ptr<ComboAttachment>  mcPolicyAtt;
    std::unique_ptr<ComboAttachment>  dialogDetAtt;
    std::unique_ptr<ComboAttachment>  dialogApplyAtt;
    std::unique_ptr<ButtonAttachment> lfeDetAtt;
    std::unique_ptr<ButtonAttachment> lfeApplyAtt;

    RoutingGraphic routingGraphic;

    // [BEGIN UI3C3-MISSIONCONTROL-RLRA-WINDOW-MEMBERS]
    juce::ToggleButton toggleMomentary { "M" };
    juce::ToggleButton toggleShortTerm { "S" };
    juce::ToggleButton toggleGate      { "Gate" };
    juce::ToggleButton toggleRolling   { "rLRA" };

    juce::ComboBox     rollingWindowBox; // 30/60/120s selector (between rLRA and Follow)

    juce::ToggleButton toggleFollow    { "Follow" };
    // [END UI3C3-MISSIONCONTROL-RLRA-WINDOW-MEMBERS]

    // State keys (persisted inside APVTS state tree as properties)
    static const juce::Identifier kTargetI;
    static const juce::Identifier kTargetPeak;
    static const juce::Identifier kTargetLra;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MissionControlComponent)
};
// [END UI3A-MISSIONCONTROL-DECL]

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
    // [BEGIN UI3A-EDITOR-MEMBERS-TOPSTRIP]
    VolumeHistoryComponent  historyComponent;
    MissionControlComponent missionControl;
    // [END UI3A-EDITOR-MEMBERS-TOPSTRIP]

    // [BEGIN MTDM-PANEL-EDITOR-MEMBERS]
    MtdmControlPanel mtdmPanel;

    // Draggable horizontal splitter between history and panel
    juce::StretchableLayoutManager layoutHistoryAndPanel;
    juce::StretchableLayoutResizerBar layoutResizerBar;
    // [END MTDM-PANEL-EDITOR-MEMBERS]

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (LevelScopeAudioProcessorEditor)
};
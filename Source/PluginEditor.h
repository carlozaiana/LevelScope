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

// [BEGIN UI4A-MTDM-PANEL-CARDS-DECL]
//==============================================================================
// MTDM Control Panel (scrollable cards)
//==============================================================================

class MtdmCardComponent : public juce::Component
{
public:
    explicit MtdmCardComponent (juce::String titleText);
    ~MtdmCardComponent() override = default;

    void paint (juce::Graphics& g) override;
    void resized() override;

protected:
    juce::Rectangle<int> getContentArea() const;

    juce::Label title;

private:
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MtdmCardComponent)
};

//------------------------------------------------------------------------------
// Levelling placeholder (future module)
//------------------------------------------------------------------------------
class LevellingCard : public MtdmCardComponent
{
public:
    LevellingCard();
    void resized() override;

private:
    juce::Label info;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (LevellingCard)
};

//------------------------------------------------------------------------------
// Zones / Thresholds card (T0..T3 + MTDM Enabled) with ordering push behavior
//------------------------------------------------------------------------------
class MtdmZonesCard : public MtdmCardComponent
{
public:
    explicit MtdmZonesCard (LevelScopeAudioProcessor& p);
    ~MtdmZonesCard() override;

    void resized() override;

private:
    LevelScopeAudioProcessor& processor;
    juce::AudioProcessorValueTreeState& apvts;

    juce::ToggleButton mtdmEnabledButton { "MTDM Enabled" };

    juce::Label  t0Label, t1Label, t2Label, t3Label;
    juce::Slider t0Slider, t1Slider, t2Slider, t3Slider;

    using ButtonAttachment = juce::AudioProcessorValueTreeState::ButtonAttachment;
    using SliderAttachment = juce::AudioProcessorValueTreeState::SliderAttachment;

    std::unique_ptr<ButtonAttachment> mtdmEnabledAtt;
    std::unique_ptr<SliderAttachment> t0Att, t1Att, t2Att, t3Att;

    // Threshold ordering enforcement (same behavior as handle drag)
    void enforceThresholdOrderingFromSlider (int changedIndex);
    void endPushedThresholdGestures();

    bool callbacksSuppressed = false;
    bool thresholdSliderDragging = false;

    std::array<bool, 4> pushedGestureActive { { false, false, false, false } };
    static constexpr float minGapLu = 0.1f;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MtdmZonesCard)
};

//------------------------------------------------------------------------------
// Upward card (essentials)
//------------------------------------------------------------------------------
class MtdmUpwardCard : public MtdmCardComponent
{
public:
    explicit MtdmUpwardCard (LevelScopeAudioProcessor& p);
    ~MtdmUpwardCard() override = default;

    void resized() override;

private:
    juce::AudioProcessorValueTreeState& apvts;

    juce::Label modeLabel;
    juce::ComboBox modeBox;

    juce::Label amountLabel, maxBoostLabel, attackLabel, releaseLabel;
    juce::Slider amountSlider, maxBoostSlider, attackSlider, releaseSlider;

    using ComboAttachment  = juce::AudioProcessorValueTreeState::ComboBoxAttachment;
    using SliderAttachment = juce::AudioProcessorValueTreeState::SliderAttachment;

    std::unique_ptr<ComboAttachment>  modeAtt;
    std::unique_ptr<SliderAttachment> amountAtt, maxBoostAtt, attackAtt, releaseAtt;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MtdmUpwardCard)
};

//------------------------------------------------------------------------------
// Downward card (essentials)
//------------------------------------------------------------------------------
class MtdmDownwardCard : public MtdmCardComponent
{
public:
    explicit MtdmDownwardCard (LevelScopeAudioProcessor& p);
    ~MtdmDownwardCard() override = default;

    void resized() override;

private:
    juce::AudioProcessorValueTreeState& apvts;

    juce::ToggleButton enabledButton { "Downward Enabled" };

    juce::Label ratioLabel, attackLabel, releaseLabel, makeupLabel;
    juce::Slider ratioSlider, attackSlider, releaseSlider, makeupSlider;

    using ButtonAttachment = juce::AudioProcessorValueTreeState::ButtonAttachment;
    using SliderAttachment = juce::AudioProcessorValueTreeState::SliderAttachment;

    std::unique_ptr<ButtonAttachment> enabledAtt;
    std::unique_ptr<SliderAttachment> ratioAtt, attackAtt, releaseAtt, makeupAtt;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MtdmDownwardCard)
};

//------------------------------------------------------------------------------
// Limiter card (essentials)
//------------------------------------------------------------------------------
class MtdmLimiterCard : public MtdmCardComponent
{
public:
    explicit MtdmLimiterCard (LevelScopeAudioProcessor& p);
    ~MtdmLimiterCard() override = default;

    void resized() override;

private:
    juce::AudioProcessorValueTreeState& apvts;

    juce::ToggleButton enabledButton { "Limiter Enabled" };

    juce::Label ceilingLabel, lookLabel, osLabel, attackLabel, releaseLabel, driveLabel;
    juce::Slider ceilingSlider, lookSlider, attackSlider, releaseSlider, driveSlider;
    juce::ComboBox osBox;

    using ButtonAttachment = juce::AudioProcessorValueTreeState::ButtonAttachment;
    using ComboAttachment  = juce::AudioProcessorValueTreeState::ComboBoxAttachment;
    using SliderAttachment = juce::AudioProcessorValueTreeState::SliderAttachment;

    std::unique_ptr<ButtonAttachment> enabledAtt;
    std::unique_ptr<ComboAttachment>  osAtt;
    std::unique_ptr<SliderAttachment> ceilingAtt, lookAtt, attackAtt, releaseAtt, driveAtt;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MtdmLimiterCard)
};

//------------------------------------------------------------------------------
// Content component (stacked cards) hosted inside a Viewport
//------------------------------------------------------------------------------
class MtdmCardsContent : public juce::Component
{
public:
    explicit MtdmCardsContent (LevelScopeAudioProcessor& p);
    ~MtdmCardsContent() override = default;

    void resized() override;
    int getPreferredHeight() const noexcept;

private:
    LevellingCard    levelling;
    MtdmZonesCard    zones;
    MtdmUpwardCard   upward;
    MtdmDownwardCard downward;
    MtdmLimiterCard  limiter;

    // [BEGIN UI4A1-CARDS-RESIZABLE-LAYOUT-MEMBERS]
    juce::StretchableLayoutManager cardsLayout;

    // Horizontal resizer bars between cards (drag up/down)
    juce::StretchableLayoutResizerBar bar01 { &cardsLayout, 1, false };
    juce::StretchableLayoutResizerBar bar12 { &cardsLayout, 3, false };
    juce::StretchableLayoutResizerBar bar23 { &cardsLayout, 5, false };
    juce::StretchableLayoutResizerBar bar34 { &cardsLayout, 7, false };

    int contentPreferredHeightPx = 0;
    // [END UI4A1-CARDS-RESIZABLE-LAYOUT-MEMBERS]

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MtdmCardsContent)
};

//------------------------------------------------------------------------------
// Public panel component used by the editor (scrollable viewport)
//------------------------------------------------------------------------------
class MtdmControlPanel : public juce::Component
{
public:
    explicit MtdmControlPanel (LevelScopeAudioProcessor& p);
    ~MtdmControlPanel() override = default;

    void paint (juce::Graphics& g) override;
    void resized() override;

private:
    juce::Viewport viewport;
    MtdmCardsContent content;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MtdmControlPanel)
};
// [END UI4A-MTDM-PANEL-CARDS-DECL]

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
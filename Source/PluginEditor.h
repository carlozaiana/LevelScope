#pragma once

#include <JuceHeader.h>
// [BEGIN MTDM-PANEL-INCLUDE-MEMORY]
#include <memory>
// [END MTDM-PANEL-INCLUDE-MEMORY]
#include "PluginProcessor.h"
#include "VolumeHistoryComponent.h"
// [BEGIN UI-CURVE-INCLUDES]
#include "DynamicsCurveComponent.h"
#include "LevelerCurveComponent.h"
// [END UI-CURVE-INCLUDES]

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

    // [BEGIN UI3A-MISSIONCONTROL-SAVELOAD-MODES]
    void showSaveMenu();
    void showLoadMenu();

    void startSaveSettingsPreset();
    void startLoadSettingsPreset();

    void startSaveSnapshot();
    void startLoadSnapshot();
    // [END UI3A-MISSIONCONTROL-SAVELOAD-MODES]

    void refreshCurveToggleStatesFromHistory();
    void updatePresetLoadEnablement();

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

    // Presets / snapshots
    juce::TextButton savePresetButton { "Save..." };
    juce::TextButton loadPresetButton { "Load..." };

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
    bool loadPresetLocked = false;

    // [BEGIN UI3A-MISSIONCONTROL-FILECHOOSER-MEMBER]
    std::unique_ptr<juce::FileChooser> activeFileChooser;
    // [END UI3A-MISSIONCONTROL-FILECHOOSER-MEMBER]

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
// Leveler card
//------------------------------------------------------------------------------
class LevelerCard : public MtdmCardComponent
{
public:
    explicit LevelerCard (LevelScopeAudioProcessor& p);
    ~LevelerCard() override = default;

    void resized() override;

private:
    juce::AudioProcessorValueTreeState& apvts;

    juce::ToggleButton enabledButton { "Enabled" };

    juce::Label targetLabel, maxBoostLabel, maxCutLabel;
    juce::Label measLabel, modeLabel;
    juce::Label rateUpLabel, rateDownLabel;

    juce::Slider targetSlider, maxBoostSlider, maxCutSlider;
    juce::Slider rateUpSlider, rateDownSlider;

    juce::ComboBox measBox, modeBox;
    juce::ToggleButton learnButton { "Learn" };

    using ButtonAttachment = juce::AudioProcessorValueTreeState::ButtonAttachment;
    using ComboAttachment  = juce::AudioProcessorValueTreeState::ComboBoxAttachment;
    using SliderAttachment = juce::AudioProcessorValueTreeState::SliderAttachment;

    std::unique_ptr<ButtonAttachment> enabledAtt;
    std::unique_ptr<SliderAttachment> targetAtt, maxBoostAtt, maxCutAtt;
    std::unique_ptr<ComboAttachment>  measAtt, modeAtt;
    std::unique_ptr<ButtonAttachment> learnAtt;
    std::unique_ptr<SliderAttachment> rateUpAtt, rateDownAtt;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (LevelerCard)
};

//------------------------------------------------------------------------------
// Zones / Thresholds card (T0..T3 + MTDM Enabled) with ordering push behavior
//------------------------------------------------------------------------------
class MtdmZonesCard : public MtdmCardComponent,
                      private juce::Timer
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

    void timerCallback() override;
    bool latencyLocked = false;

    std::array<bool, 4> pushedGestureActive { { false, false, false, false } };
    static constexpr float minGapLu = 0.1f;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MtdmZonesCard)
};

// [BEGIN UI-ZONE-AUDITION-CARD-DECL]
//------------------------------------------------------------------------------
// Zone Audition card (NEW): 5 loudness zones, combinable (OR), empty selection = OFF
//------------------------------------------------------------------------------
class MtdmAuditionCard : public MtdmCardComponent
{
public:
    explicit MtdmAuditionCard (LevelScopeAudioProcessor& p);
    ~MtdmAuditionCard() override = default;

    void resized() override;

private:
    juce::AudioProcessorValueTreeState& apvts;

    // 5 combinable zone toggles
    juce::ToggleButton zBelowT0 { "Below T0" };
    juce::ToggleButton zT0T1    { "T0..T1" };
    juce::ToggleButton zT1T2    { "T1..T2" };
    juce::ToggleButton zT2T3    { "T2..T3" };
    juce::ToggleButton zAboveT3 { "Above T3" };

    // Optional helpers
    juce::TextButton clearButton { "Clear" };
    juce::TextButton allButton   { "All" };
    juce::Label      activeLabel;

    using ButtonAttachment = juce::AudioProcessorValueTreeState::ButtonAttachment;

    std::unique_ptr<ButtonAttachment> attBelowT0, attT0T1, attT1T2, attT2T3, attAboveT3;

    void setAllZoneAuditionToggles (bool newState);

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MtdmAuditionCard)
};
// [END UI-ZONE-AUDITION-CARD-DECL]

//------------------------------------------------------------------------------
// Upward card (essentials)
//------------------------------------------------------------------------------
class MtdmUpwardCard : public MtdmCardComponent,
                       private juce::Timer
{
public:
    explicit MtdmUpwardCard (LevelScopeAudioProcessor& p);
    ~MtdmUpwardCard() override = default;

    void paint (juce::Graphics& g) override;
    void resized() override;

    bool isAdvancedOpen() const noexcept { return showAdvanced; }
    void applyPersistedAdvancedOpenState (bool shouldOpen);

private:
    LevelScopeAudioProcessor& processor;
    juce::AudioProcessorValueTreeState& apvts;

    // [BEGIN UI4B1-UPWARD-ADVANCED-MEMBERS]
    // [BEGIN UI-UP-EN-BYP-HEADER-MEMBERS]
    juce::ToggleButton enabledButton { "Upward" };
    juce::ToggleButton bypassButton  { "Bypass" };

    juce::Label modeLabel;
    juce::ComboBox modeBox;

    juce::ToggleButton advancedToggle { "Advanced" };
    bool showAdvanced = false;
    // [END UI-UP-EN-BYP-HEADER-MEMBERS]

    // Essentials
    juce::Label amountLabel, maxBoostLabel, attackLabel, releaseLabel;
    juce::Slider amountSlider, maxBoostSlider, attackSlider, releaseSlider;

    // Advanced (curve + knees)
    juce::Label curveLabel, curveTypeLabel, lowKneeLabel, highKneeLabel, calTrimLabel;
    juce::Slider curveSlider, lowKneeSlider, highKneeSlider, calTrimSlider;
    juce::ComboBox curveTypeBox;

    // Advanced (spectral-only controls)
    juce::Label fftLabel, bandsLabel, minFreqLabel, maxFreqLabel;
    juce::ComboBox fftBox, bandsBox;
    juce::Slider minFreqSlider, maxFreqSlider;

    // [BEGIN UI-UP-ADD-BUTTONATTACHMENT-TYPEDEF]
    using ButtonAttachment = juce::AudioProcessorValueTreeState::ButtonAttachment;
    using ComboAttachment  = juce::AudioProcessorValueTreeState::ComboBoxAttachment;
    using SliderAttachment = juce::AudioProcessorValueTreeState::SliderAttachment;
    // [END UI-UP-ADD-BUTTONATTACHMENT-TYPEDEF]

    // [BEGIN UI-UP-EN-BYP-ATTACHMENTS]
    std::unique_ptr<ButtonAttachment> enabledAtt;
    std::unique_ptr<ButtonAttachment> bypassAtt;

    std::unique_ptr<ComboAttachment>  modeAtt;
    std::unique_ptr<SliderAttachment> amountAtt, maxBoostAtt, attackAtt, releaseAtt;
    // [END UI-UP-EN-BYP-ATTACHMENTS]

    std::unique_ptr<SliderAttachment> curveAtt, lowKneeAtt, highKneeAtt, calTrimAtt, minFreqAtt, maxFreqAtt;
    std::unique_ptr<ComboAttachment>  curveTypeAtt, fftAtt, bandsAtt;

    void updateAdvancedVisibility();
    void updateSpectralEnablement();
    void timerCallback() override;
    bool latencyLocked = false;
    // [END UI4B1-UPWARD-ADVANCED-MEMBERS]

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

    bool isAdvancedOpen() const noexcept { return showAdvanced; }
    void applyPersistedAdvancedOpenState (bool shouldOpen);

private:
    LevelScopeAudioProcessor& processor;
    juce::AudioProcessorValueTreeState& apvts;

    // [BEGIN UI4B3-DOWNWARD-ADVANCED-MEMBERS]
    // [BEGIN UI-DOWN-BYP-MEMBERS]
    juce::ToggleButton enabledButton { "Downward Enabled" };
    juce::ToggleButton bypassButton  { "Bypass" };

    juce::ToggleButton advancedToggle { "Advanced" };
    bool showAdvanced = false;
    // [END UI-DOWN-BYP-MEMBERS]

    // Essentials
    juce::Label ratioLabel, attackLabel, releaseLabel, makeupLabel;
    juce::Slider ratioSlider, attackSlider, releaseSlider, makeupSlider;

    // Advanced
    juce::Label kneeLabel;
    juce::Slider kneeSlider;

    void updateAdvancedVisibility();
    // [END UI4B3-DOWNWARD-ADVANCED-MEMBERS]

    using ButtonAttachment = juce::AudioProcessorValueTreeState::ButtonAttachment;
    using SliderAttachment = juce::AudioProcessorValueTreeState::SliderAttachment;

    // [BEGIN UI4B3-DOWNWARD-ADVANCED-ATTACHMENTS]
    // [BEGIN UI-DOWN-BYP-ATTACH]
    std::unique_ptr<ButtonAttachment> enabledAtt;
    std::unique_ptr<ButtonAttachment> bypassAtt;
    // [END UI-DOWN-BYP-ATTACH]
    std::unique_ptr<SliderAttachment> ratioAtt, attackAtt, releaseAtt, makeupAtt;
    std::unique_ptr<SliderAttachment> kneeAtt;
    // [END UI4B3-DOWNWARD-ADVANCED-ATTACHMENTS]

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MtdmDownwardCard)
};

//------------------------------------------------------------------------------
// Limiter card (essentials)
//------------------------------------------------------------------------------
class MtdmLimiterCard : public MtdmCardComponent,
                        private juce::Timer
{
public:
    explicit MtdmLimiterCard (LevelScopeAudioProcessor& p);
    ~MtdmLimiterCard() override = default;

    void paint (juce::Graphics& g) override;
    void resized() override;

    bool isAdvancedOpen() const noexcept { return showAdvanced; }
    void applyPersistedAdvancedOpenState (bool shouldOpen);

private:
    LevelScopeAudioProcessor& processor;
    juce::AudioProcessorValueTreeState& apvts;

    // [BEGIN UI4B2-LIMITER-ADVANCED-MEMBERS]
    // [BEGIN UI4B2-LIMITER-ENABLED-TEXT]
    // [BEGIN UI-LIM-BYP-MEMBERS]
    juce::ToggleButton enabledButton { "Limiter" };
    juce::ToggleButton bypassButton  { "Bypass" };

    juce::ToggleButton advancedToggle { "Advanced" };
    bool showAdvanced = false;
    // [END UI-LIM-BYP-MEMBERS]

    // Essentials
    juce::Label ceilingLabel, lookLabel, osLabel;
    juce::Slider ceilingSlider, lookSlider;
    juce::ComboBox osBox;

    // Advanced
    juce::Label attackLabel, releaseLabel, driveLabel;
    juce::Slider attackSlider, releaseSlider, driveSlider;

    void updateAdvancedVisibility();
    void updateLatencySensitiveEnablement();
    void timerCallback() override;
    bool latencyLocked = false;
    // [END UI4B2-LIMITER-ADVANCED-MEMBERS]

    using ButtonAttachment = juce::AudioProcessorValueTreeState::ButtonAttachment;
    using ComboAttachment  = juce::AudioProcessorValueTreeState::ComboBoxAttachment;
    using SliderAttachment = juce::AudioProcessorValueTreeState::SliderAttachment;

    // [BEGIN UI-LIM-BYP-ATTACH]
    std::unique_ptr<ButtonAttachment> enabledAtt;
    std::unique_ptr<ButtonAttachment> bypassAtt;
    // [END UI-LIM-BYP-ATTACH]
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

    // [BEGIN UI-CURVE-SHARED-STRIP-DECL]
    int getCurveStripWidthForContent (int totalContentWidth) const noexcept;
    // [END UI-CURVE-SHARED-STRIP-DECL]

    void applyPersistedUiStateFromProcessor();

    // [BEGIN UI4B1-UPWARD-AUTOEXPAND-DECL]
    void ensureUpwardHeightAtLeast (int px);
    // [END UI4B1-UPWARD-AUTOEXPAND-DECL]

    // [BEGIN UI4B2-LIMITER-AUTOEXPAND-DECL]
    void ensureLimiterHeightAtLeast (int px);
    // [END UI4B2-LIMITER-AUTOEXPAND-DECL]

    // [BEGIN UI4B3-CARDS-SHRINK-DECL]
    void setUpwardHeightPx  (int px);
    void setLimiterHeightPx (int px);
    // [END UI4B3-CARDS-SHRINK-DECL]

    // [BEGIN UI4B3-DOWNWARD-HEIGHT-DECL]
    void ensureDownwardHeightAtLeast (int px);
    void setDownwardHeightPx (int px);
    // [END UI4B3-DOWNWARD-HEIGHT-DECL]

private:
    LevelScopeAudioProcessor& processor;

    // [BEGIN UI4C-CONTENT-ADD-AUDITION-MEMBER]
    LevelerCard      levelling;
    MtdmZonesCard    zones;
    MtdmAuditionCard audition;
    MtdmUpwardCard   upward;
    // [END UI4C-CONTENT-ADD-AUDITION-MEMBER]
    MtdmDownwardCard downward;
    MtdmLimiterCard  limiter;

    // [BEGIN UI4C-CARDHEIGHTS-ADD-AUDITION]
    struct CardHeights
    {
        int levelling = 240;
        int zones     = 220;
        int audition  = 110;
        int upward    = 220;
        int downward  = 180;
        int limiter   = 260;
    };
    // [END UI4C-CARDHEIGHTS-ADD-AUDITION]

    CardHeights cardHeights;

    static constexpr int barHeightPx = 8;

    // [BEGIN UI-CURVE-SHARED-STRIP-CONSTANTS]
    static constexpr int curveStripDefaultWidthPx    = 280;
    static constexpr int curveStripMinWidthPx        = 180;
    static constexpr int curveStripMaxWidthPx        = 520;
    static constexpr int curveStripMinControlWidthPx = 420;

    int curveStripWidthPxUser = curveStripDefaultWidthPx;
    int dragStartCurveStripWidthPxUser = curveStripDefaultWidthPx;
    // [END UI-CURVE-SHARED-STRIP-CONSTANTS]

    static constexpr int minLevellingPx = 50;
    static constexpr int minZonesPx     = 34;
    // [BEGIN UI4C-MINHEIGHT-AUDITION]
    static constexpr int minAuditionPx  = 34;
    // [END UI4C-MINHEIGHT-AUDITION]
    static constexpr int minUpwardPx    = 34;
    static constexpr int minDownwardPx  = 26;
    static constexpr int minLimiterPx   = 26;

    // [BEGIN UI-CURVE-SHARED-STRIP-RESIZER-DECL]
    class CurveStripResizerBar : public juce::Component
    {
    public:
        explicit CurveStripResizerBar (MtdmCardsContent& ownerIn);

        void paint (juce::Graphics& g) override;
        void mouseDown (const juce::MouseEvent& e) override;
        void mouseDrag (const juce::MouseEvent& e) override;

    private:
        MtdmCardsContent& owner;
        int dragStartScreenX = 0;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (CurveStripResizerBar)
    };
    // [END UI-CURVE-SHARED-STRIP-RESIZER-DECL]

    class CardResizerBar : public juce::Component
    {
    public:
        // [BEGIN UI4C-BOUNDARY-ADD-AUDITION]
        enum class Boundary
        {
            levellingZones,
            zonesAudition,
            auditionUpward,
            upwardDownward,
            downwardLimiter,
            limiterTail
        };
        // [END UI4C-BOUNDARY-ADD-AUDITION]

        CardResizerBar (MtdmCardsContent& ownerIn, Boundary b);

        void paint (juce::Graphics& g) override;

        void mouseEnter (const juce::MouseEvent&) override;
        void mouseExit  (const juce::MouseEvent&) override;
        void mouseDown  (const juce::MouseEvent& e) override;
        void mouseDrag  (const juce::MouseEvent& e) override;

    private:
        MtdmCardsContent& owner;
        Boundary boundary;

        // [BEGIN UI4A3-RESIZER-SCREENPOS]
        int dragStartScreenY = 0; // stable even if the component moves during drag
        CardHeights dragStartHeights;
        // [END UI4A3-RESIZER-SCREENPOS]

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (CardResizerBar)
    };

    // [BEGIN UI4C-BARS-REPLACE]
    CardResizerBar bar01 { *this, CardResizerBar::Boundary::levellingZones };
    CardResizerBar bar12 { *this, CardResizerBar::Boundary::zonesAudition };
    CardResizerBar bar23 { *this, CardResizerBar::Boundary::auditionUpward };
    CardResizerBar bar34 { *this, CardResizerBar::Boundary::upwardDownward };
    CardResizerBar bar45 { *this, CardResizerBar::Boundary::downwardLimiter };
    CardResizerBar bar56 { *this, CardResizerBar::Boundary::limiterTail };

    // [BEGIN UI-CURVE-SHARED-STRIP-MEMBER]
    CurveStripResizerBar curveStripBar { *this };
    // [END UI-CURVE-SHARED-STRIP-MEMBER]
    // [END UI4C-BARS-REPLACE]

    void applyDragToBoundary (CardResizerBar::Boundary b, int dy);
    void syncPersistedCardHeightsToProcessor() const;

    int contentPreferredHeightPx = 0;
    // [END UI4A3-CARDS-ACCORDION-MEMBERS]

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

    void applyPersistedUiStateFromProcessor();

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

    void applyPersistedUiStateFromProcessor();

private:
    LevelScopeAudioProcessor& processor;

    // [BEGIN UI3A-EDITOR-MEMBERS-TOPSTRIP]
    VolumeHistoryComponent  historyComponent;
    MissionControlComponent missionControl;
    // [END UI3A-EDITOR-MEMBERS-TOPSTRIP]

    // [BEGIN MTDM-PANEL-EDITOR-MEMBERS]
    MtdmControlPanel mtdmPanel;

    // Draggable horizontal splitter between history and panel
    juce::StretchableLayoutManager layoutHistoryAndPanel;
    juce::StretchableLayoutResizerBar layoutResizerBar;

    // Enables JUCE tooltips for controls in the editor
    juce::TooltipWindow tooltipWindow;
    // [END MTDM-PANEL-EDITOR-MEMBERS]

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (LevelScopeAudioProcessorEditor)
};
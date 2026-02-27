#include "PluginEditor.h"
// [BEGIN MTDM-PANEL-INCLUDE-PARAMIDS]
#include "Core/Processing/Modules/MultiThresholdDynamicsParamIDs.h"
// [END MTDM-PANEL-INCLUDE-PARAMIDS]

// [BEGIN UI3A-MISSIONCONTROL-IMPL]
const juce::Identifier MissionControlComponent::kTargetI    ("ui.targetIntegratedLufs");
const juce::Identifier MissionControlComponent::kTargetPeak ("ui.targetMaxPeakDb");
const juce::Identifier MissionControlComponent::kTargetLra  ("ui.targetLraLu");

MissionControlComponent::RoutingGraphic::RoutingGraphic (LevelScopeAudioProcessor& p)
    : processor (p)
{
    setOpaque (false);
}

// [BEGIN UI3A-ROUTINGGRAPHIC-CHANNEL-LABEL-PORTABLE]
juce::String MissionControlComponent::RoutingGraphic::channelLabelForType (juce::AudioChannelSet::ChannelType t) const
{
    // Use JUCE's built-in abbreviated labels to stay compatible across JUCE versions
    // and all supported layouts (stereo .. 7.1.4).
    auto s = juce::AudioChannelSet::getAbbreviatedChannelTypeName (t);

    if (s.isNotEmpty())
        return s;

    // Fallback (shouldn't happen)
    return "?";
}
// [END UI3A-ROUTINGGRAPHIC-CHANNEL-LABEL-PORTABLE]

bool MissionControlComponent::RoutingGraphic::isDetectorChannelActive (juce::AudioChannelSet::ChannelType t,
                                                                      int mcPolicy, int detChoice, bool lfeInDet) const
{
    using CT = juce::AudioChannelSet::ChannelType;

    const bool isLfe = (t == CT::LFE);
    if (isLfe)
        return lfeInDet;

    // mcPolicy: 0 Linked, 1 Dialog-mask, 2 Unlinked  (per ParamIDs string array)
    if (mcPolicy == 1) // Dialog-mask
    {
        if (detChoice == 0) // C
            return t == CT::centre;
        else                // LCR
            return (t == CT::left || t == CT::centre || t == CT::right);
    }

    // Linked/Unlinked: detector uses all non-LFE channels
    return true;
}

bool MissionControlComponent::RoutingGraphic::isApplyChannelActive (juce::AudioChannelSet::ChannelType t,
                                                                   int mcPolicy, int applyChoice, bool lfeInApply) const
{
    using CT = juce::AudioChannelSet::ChannelType;

    const bool isLfe = (t == CT::LFE);
    if (isLfe)
        return lfeInApply;

    if (mcPolicy == 1) // Dialog-mask
    {
        if (applyChoice == 0) // C
            return t == CT::centre;
        else                  // LCR
            return (t == CT::left || t == CT::centre || t == CT::right);
    }

    // Linked/Unlinked: apply uses all non-LFE channels
    return true;
}

void MissionControlComponent::RoutingGraphic::paint (juce::Graphics& g)
{
    auto r = getLocalBounds().toFloat().reduced (6.0f, 6.0f);

    // [BEGIN UI3A-ROUTINGGRAPHIC-VISIBILITY]
    g.setColour (juce::Colours::white.withMultipliedAlpha (0.22f));
    // [END UI3A-ROUTINGGRAPHIC-VISIBILITY]
    g.drawRoundedRectangle (r, 6.0f, 1.0f);

    // [BEGIN UI3A-ROUTINGGRAPHIC-SMALL-AREA-FIX]
    r = r.reduced (4.0f, 4.0f);
    if (r.getWidth() < 60.0f || r.getHeight() < 18.0f)
        return;
    // [END UI3A-ROUTINGGRAPHIC-SMALL-AREA-FIX]

    auto& apvts = processor.getAPVTS();
    using namespace levelscope::mtdm::ParamIDs;

    const int mcPolicy = (int) std::lround (apvts.getRawParameterValue (mcPolicyChoice)->load());
    const int detChoice = (int) std::lround (apvts.getRawParameterValue (dialogDetectorChoice)->load());
    const int appChoice = (int) std::lround (apvts.getRawParameterValue (dialogApplyChoice)->load());

    const bool lfeDet = (apvts.getRawParameterValue (lfeInDetector)->load() >= 0.5f);
    const bool lfeApp = (apvts.getRawParameterValue (lfeInApply)->load() >= 0.5f);

    const auto layout = processor.getBusesLayout().getMainInputChannelSet();
    const int numCh = layout.size();

    const float rowH = r.getHeight() * 0.5f;
    const float cellW = (numCh > 0 ? r.getWidth() / (float) numCh : r.getWidth());

    auto drawRow = [&] (float y0, const juce::String& rowName, auto isActiveFn)
    {
        g.setColour (juce::Colours::white.withMultipliedAlpha (0.60f));
        g.setFont (12.0f);
        g.drawText (rowName, (int) r.getX() - 34, (int) y0, 32, (int) rowH, juce::Justification::centredRight);

        for (int ch = 0; ch < numCh; ++ch)
        {
            const auto t = layout.getTypeOfChannel (ch);
            const bool active = isActiveFn (t);

            juce::Rectangle<float> cell (r.getX() + cellW * (float) ch + 1.0f,
                                         y0 + 2.0f,
                                         cellW - 2.0f,
                                         rowH - 4.0f);

            g.setColour (active ? juce::Colours::red.withMultipliedAlpha (0.75f)
                                : juce::Colours::grey.withMultipliedAlpha (0.25f));
            g.fillRoundedRectangle (cell, 3.0f);

            g.setColour (juce::Colours::white.withMultipliedAlpha (0.80f));
            g.setFont (11.0f);
            g.drawFittedText (channelLabelForType (t), cell.toNearestInt(), juce::Justification::centred, 1);
        }
    };

    drawRow (r.getY(), "Det", [&] (juce::AudioChannelSet::ChannelType t)
    {
        return isDetectorChannelActive (t, mcPolicy, detChoice, lfeDet);
    });

    drawRow (r.getY() + rowH, "App", [&] (juce::AudioChannelSet::ChannelType t)
    {
        return isApplyChannelActive (t, mcPolicy, appChoice, lfeApp);
    });
}

MissionControlComponent::MissionControlComponent (LevelScopeAudioProcessor& p, VolumeHistoryComponent& h)
    : processor (p),
      history (h),
      apvts (p.getAPVTS()),
      routingGraphic (p)
{
    setOpaque (true);

    // Presets
    addAndMakeVisible (savePresetButton);
    addAndMakeVisible (loadPresetButton);

    savePresetButton.onClick = [this] { startSavePreset(); };
    loadPresetButton.onClick = [this] { startLoadPreset(); };

    auto setupLabelBox = [] (juce::Label& l)
    {
        l.setColour (juce::Label::textColourId, juce::Colours::white.withMultipliedAlpha (0.90f));
        l.setColour (juce::Label::backgroundColourId, juce::Colours::black.withMultipliedAlpha (0.15f));
        l.setColour (juce::Label::outlineColourId, juce::Colours::white.withMultipliedAlpha (0.10f));
        l.setJustificationType (juce::Justification::centred);
        l.setFont (juce::Font (13.0f));
    };

    // [BEGIN UI3A-MISSIONCONTROL-HEADERS-SETUP]
    auto setupHeader = [] (juce::Label& l)
    {
        l.setEditable (false);
        l.setJustificationType (juce::Justification::centred);
        l.setColour (juce::Label::textColourId, juce::Colours::white.withMultipliedAlpha (0.70f));
        l.setFont (juce::Font (11.0f));
    };

    setupHeader (hdrILabel);
    setupHeader (hdrPeakLabel);
    setupHeader (hdrLraLabel);

    hdrILabel.setText    ("Int (LUFS)", juce::dontSendNotification);
    hdrPeakLabel.setText ("Max Peak",   juce::dontSendNotification);
    hdrLraLabel.setText  ("LRA (LU)",   juce::dontSendNotification);

    addAndMakeVisible (hdrILabel);
    addAndMakeVisible (hdrPeakLabel);
    addAndMakeVisible (hdrLraLabel);

    setupHeader (hdrTargetRow);
    setupHeader (hdrCurrentRow);

    hdrTargetRow.setText  ("Target",  juce::dontSendNotification);
    hdrCurrentRow.setText ("Current", juce::dontSendNotification);

    hdrTargetRow.setJustificationType  (juce::Justification::centredLeft);
    hdrCurrentRow.setJustificationType (juce::Justification::centredLeft);

    addAndMakeVisible (hdrTargetRow);
    addAndMakeVisible (hdrCurrentRow);
    // [END UI3A-MISSIONCONTROL-HEADERS-SETUP]

    // Target labels (editable)
    setupLabelBox (targetILabel);
    setupLabelBox (targetPeakLabel);
    setupLabelBox (targetLraLabel);

    targetILabel.setEditable (true);
    targetPeakLabel.setEditable (true);
    targetLraLabel.setEditable (true);

    addAndMakeVisible (targetILabel);
    addAndMakeVisible (targetPeakLabel);
    addAndMakeVisible (targetLraLabel);

    // Current labels (read-only)
    setupLabelBox (currentILabel);
    setupLabelBox (currentPeakLabel);
    setupLabelBox (currentLraLabel);

    addAndMakeVisible (currentILabel);
    addAndMakeVisible (currentPeakLabel);
    addAndMakeVisible (currentLraLabel);

    // Policy controls
    using namespace levelscope::mtdm::ParamIDs;

    mcPolicyBox.addItemList (juce::StringArray { "Linked", "Dialog-mask", "Unlinked" }, 1);
    dialogDetBox.addItemList (juce::StringArray { "C", "LCR" }, 1);
    dialogApplyBox.addItemList (juce::StringArray { "C", "LCR" }, 1);

    addAndMakeVisible (mcPolicyBox);
    addAndMakeVisible (dialogDetBox);
    addAndMakeVisible (dialogApplyBox);
    addAndMakeVisible (lfeDetToggle);
    addAndMakeVisible (lfeApplyToggle);

    // [BEGIN UI3A-MISSIONCONTROL-POLICY-COLOURS]
    lfeDetToggle.setColour   (juce::ToggleButton::textColourId, juce::Colours::white.withMultipliedAlpha (0.90f));
    lfeApplyToggle.setColour (juce::ToggleButton::textColourId, juce::Colours::white.withMultipliedAlpha (0.90f));

    // Shorter labels so they fit in the right column
    lfeDetToggle.setButtonText   ("LFE Det");
    lfeApplyToggle.setButtonText ("LFE Apply");
    // [END UI3A-MISSIONCONTROL-POLICY-COLOURS]

    mcPolicyAtt   = std::make_unique<ComboAttachment>  (apvts, mcPolicyChoice, mcPolicyBox);
    dialogDetAtt  = std::make_unique<ComboAttachment>  (apvts, dialogDetectorChoice, dialogDetBox);
    dialogApplyAtt= std::make_unique<ComboAttachment>  (apvts, dialogApplyChoice, dialogApplyBox);
    lfeDetAtt     = std::make_unique<ButtonAttachment> (apvts, lfeInDetector, lfeDetToggle);
    lfeApplyAtt   = std::make_unique<ButtonAttachment> (apvts, lfeInApply,    lfeApplyToggle);

    addAndMakeVisible (routingGraphic);

    // Curve toggles row (bottom)
    auto setupToggle = [] (juce::ToggleButton& b)
    {
        b.setClickingTogglesState (true);
        b.setColour (juce::ToggleButton::textColourId, juce::Colours::white.withMultipliedAlpha (0.90f));
    };

    setupToggle (toggleMomentary);
    setupToggle (toggleShortTerm);
    setupToggle (toggleGate);
    setupToggle (toggleRolling);

    addAndMakeVisible (toggleMomentary);
    addAndMakeVisible (toggleShortTerm);
    addAndMakeVisible (toggleGate);
    addAndMakeVisible (toggleRolling);

    toggleMomentary.onClick = [this] { history.setShowMomentaryCurve (toggleMomentary.getToggleState()); };
    toggleShortTerm.onClick = [this] { history.setShowShortTermCurve (toggleShortTerm.getToggleState()); };
    toggleGate.onClick      = [this] { history.setShowGateCurve (toggleGate.getToggleState()); };
    toggleRolling.onClick   = [this] { history.setShowRollingLraLane (toggleRolling.getToggleState()); };

    // Target persistence
    loadTargetsFromState();

    auto hookTarget = [this] (juce::Label& lab, const juce::Identifier& key, double def)
    {
        lab.onTextChange = [this, &lab, key, def]
        {
            const double v = lab.getText().trim().getDoubleValue();
            storeTargetToState (key, v);
            // reformat
            const double vv = getTargetFromState (key, def);
            lab.setText (juce::String (vv, 1), juce::dontSendNotification);
        };
    };

    hookTarget (targetILabel,    kTargetI,    -23.0);
    hookTarget (targetPeakLabel, kTargetPeak, -1.0);
    hookTarget (targetLraLabel,  kTargetLra,   7.0);

    // Start polling current readouts + keep toggles synced
    startTimerHz (10);
}

MissionControlComponent::~MissionControlComponent()
{
    stopTimer();
}

void MissionControlComponent::paint (juce::Graphics& g)
{
    g.fillAll (juce::Colour::fromRGB (10, 18, 28));
    g.setColour (juce::Colours::white.withMultipliedAlpha (0.08f));
    g.drawRect (getLocalBounds());
}

double MissionControlComponent::getTargetFromState (const juce::Identifier& key, double defaultValue) const
{
    const auto& vt = apvts.state;
    if (vt.hasProperty (key))
        return (double) vt.getProperty (key);
    return defaultValue;
}

void MissionControlComponent::storeTargetToState (const juce::Identifier& key, double value)
{
    apvts.state.setProperty (key, value, nullptr);
}

void MissionControlComponent::loadTargetsFromState()
{
    const double ti = getTargetFromState (kTargetI, -23.0);
    const double tp = getTargetFromState (kTargetPeak, -1.0);
    const double tl = getTargetFromState (kTargetLra, 7.0);

    targetILabel.setText    (juce::String (ti, 1), juce::dontSendNotification);
    targetPeakLabel.setText (juce::String (tp, 1), juce::dontSendNotification);
    targetLraLabel.setText  (juce::String (tl, 1), juce::dontSendNotification);
}

void MissionControlComponent::updateCurrentReadouts()
{
    const float I = processor.getRunningIntegratedLufs();
    const float LRA = processor.getRunningLraLu();

    currentILabel.setText   ((I > -199.0f ? juce::String (I, 1) : "--"), juce::dontSendNotification);
    currentLraLabel.setText ((LRA > 0.0f ? juce::String (LRA, 1) : "--"), juce::dontSendNotification);

    // Max peak not implemented yet in analysis pipeline -> placeholder
    currentPeakLabel.setText ("--", juce::dontSendNotification);
}

void MissionControlComponent::refreshCurveToggleStatesFromHistory()
{
    toggleMomentary.setToggleState (history.getShowMomentaryCurve(), juce::dontSendNotification);
    toggleShortTerm.setToggleState (history.getShowShortTermCurve(), juce::dontSendNotification);
    toggleGate.setToggleState      (history.getShowGateCurve(), juce::dontSendNotification);
    toggleRolling.setToggleState   (history.getShowRollingLraLane(), juce::dontSendNotification);
}

void MissionControlComponent::timerCallback()
{
    updateCurrentReadouts();
    refreshCurveToggleStatesFromHistory();
    routingGraphic.repaint();
}

void MissionControlComponent::startSavePreset()
{
    juce::FileChooser chooser ("Save LevelScope preset",
                              juce::File::getSpecialLocation (juce::File::userDocumentsDirectory),
                              "*.lscpreset");

    chooser.launchAsync (juce::FileBrowserComponent::saveMode | juce::FileBrowserComponent::canSelectFiles,
                         [this] (const juce::FileChooser& fc)
                         {
                             auto f = fc.getResult();
                             if (f == juce::File())
                                 return;

                             juce::MemoryBlock mb;
                             processor.getStateInformation (mb);

                             (void) f.replaceWithData (mb.getData(), mb.getSize());
                         });
}

void MissionControlComponent::startLoadPreset()
{
    juce::FileChooser chooser ("Load LevelScope preset",
                              juce::File::getSpecialLocation (juce::File::userDocumentsDirectory),
                              "*.lscpreset");

    chooser.launchAsync (juce::FileBrowserComponent::openMode | juce::FileBrowserComponent::canSelectFiles,
                         [this] (const juce::FileChooser& fc)
                         {
                             auto f = fc.getResult();
                             if (f == juce::File() || ! f.existsAsFile())
                                 return;

                             juce::MemoryBlock mb;
                             if (! f.loadFileAsData (mb))
                                 return;

                             processor.setStateInformation (mb.getData(), (int) mb.getSize());
                             loadTargetsFromState();
                         });
}

void MissionControlComponent::resized()
{
    auto r = getLocalBounds().reduced (8);

    // Bottom row: curve toggles (minimal height, broad buttons)
    auto toggles = r.removeFromBottom (22);
    const int tW = 70;
    toggleMomentary.setBounds (toggles.removeFromLeft (tW));
    toggleShortTerm.setBounds (toggles.removeFromLeft (tW));
    toggleGate.setBounds      (toggles.removeFromLeft (tW));
    toggleRolling.setBounds   (toggles.removeFromLeft (tW));

    r.removeFromBottom (6);

    // Right: routing graphic + policy controls
    // [BEGIN UI3A-MISSIONCONTROL-RIGHT-COL-WIDTH-FIX]
    const int rightW = juce::jlimit (320, 560, r.getWidth() / 2);
    auto right = r.removeFromRight (rightW);
    // [END UI3A-MISSIONCONTROL-RIGHT-COL-WIDTH-FIX]
    // [BEGIN UI3A-MISSIONCONTROL-POLICY-ROW-LAYOUT]
    auto policyRow = right.removeFromTop (22);

    mcPolicyBox.setBounds    (policyRow.removeFromLeft (130));
    dialogDetBox.setBounds   (policyRow.removeFromLeft (60));
    dialogApplyBox.setBounds (policyRow.removeFromLeft (60));
    lfeDetToggle.setBounds   (policyRow.removeFromLeft (80));
    lfeApplyToggle.setBounds (policyRow.removeFromLeft (90));
    // [END UI3A-MISSIONCONTROL-POLICY-ROW-LAYOUT]

    right.removeFromTop (6);
    routingGraphic.setBounds (right);

    // Left: presets + targets/current aligned
    auto left = r;

    auto topRow = left.removeFromTop (22);
    savePresetButton.setBounds (topRow.removeFromLeft (120));
    topRow.removeFromLeft (6);
    loadPresetButton.setBounds (topRow.removeFromLeft (120));

    left.removeFromTop (8);

    // [BEGIN UI3A-MISSIONCONTROL-TARGETS-CURRENT-LAYOUT]
    const int rowHdrW = 62;
    const int colW    = 96;

    // Column headers
    auto hdrRow = left.removeFromTop (16);
    hdrRow.removeFromLeft (rowHdrW);
    hdrILabel.setBounds    (hdrRow.removeFromLeft (colW));
    hdrPeakLabel.setBounds (hdrRow.removeFromLeft (colW));
    hdrLraLabel.setBounds  (hdrRow.removeFromLeft (colW));

    left.removeFromTop (2);

    // Targets row
    auto targetRow = left.removeFromTop (26);
    hdrTargetRow.setBounds (targetRow.removeFromLeft (rowHdrW));
    targetILabel.setBounds    (targetRow.removeFromLeft (colW));
    targetPeakLabel.setBounds (targetRow.removeFromLeft (colW));
    targetLraLabel.setBounds  (targetRow.removeFromLeft (colW));

    // Current row
    auto currentRow = left.removeFromTop (26);
    hdrCurrentRow.setBounds (currentRow.removeFromLeft (rowHdrW));
    currentILabel.setBounds    (currentRow.removeFromLeft (colW));
    currentPeakLabel.setBounds (currentRow.removeFromLeft (colW));
    currentLraLabel.setBounds  (currentRow.removeFromLeft (colW));
    // [END UI3A-MISSIONCONTROL-TARGETS-CURRENT-LAYOUT]
}
// [END UI3A-MISSIONCONTROL-IMPL]

//==============================================================================
// [BEGIN MTDM-PANEL-IMPL]
//==============================================================================
// GR meter (simple read-only bar)
//==============================================================================

void GrMeterComponent::setValuesDb (float currentDb, float holdDb) noexcept
{
    current = juce::jmax (0.0f, currentDb);
    hold    = juce::jmax (0.0f, holdDb);
}

void GrMeterComponent::paint (juce::Graphics& g)
{
    auto r = getLocalBounds().toFloat();

    g.setColour (juce::Colour::fromRGB (10, 18, 28));
    g.fillRoundedRectangle (r, 4.0f);

    r = r.reduced (6.0f, 5.0f);
    if (r.getWidth() <= 2.0f || r.getHeight() <= 2.0f)
        return;

    // Meter scale
    constexpr float maxDb = 24.0f;
    const float curNorm  = juce::jlimit (0.0f, 1.0f, current / maxDb);
    const float holdNorm = juce::jlimit (0.0f, 1.0f, hold    / maxDb);

    // Background bar
    g.setColour (juce::Colours::black.withMultipliedAlpha (0.35f));
    g.fillRoundedRectangle (r, 3.0f);

    // Filled current
    auto filled = r.withWidth (r.getWidth() * curNorm);
    g.setColour (juce::Colours::limegreen.withMultipliedAlpha (0.80f));
    g.fillRoundedRectangle (filled, 3.0f);

    // Hold marker
    const float holdX = r.getX() + r.getWidth() * holdNorm;
    g.setColour (juce::Colours::white.withMultipliedAlpha (0.85f));
    g.drawLine (holdX, r.getY(), holdX, r.getBottom(), 1.2f);

    // Text
    g.setColour (juce::Colours::white.withMultipliedAlpha (0.90f));
    g.setFont (12.0f);

    const juce::String txt = name + "  GR " + juce::String (current, 1) + " dB";
    g.drawText (txt, getLocalBounds().reduced (8, 2), juce::Justification::centredLeft, true);
}

//==============================================================================
// MTDM control panel
//==============================================================================

MtdmControlPanel::MtdmControlPanel (LevelScopeAudioProcessor& p)
    : processor (p),
      apvts (p.getAPVTS())
{
    setOpaque (true);

    configureToggle (mtdmEnabledButton, "MTDM");
    configureToggle (downEnabledButton, "Down");
    configureToggle (limEnabledButton,  "Lim");

    addAndMakeVisible (mtdmEnabledButton);
    addAndMakeVisible (downEnabledButton);
    addAndMakeVisible (limEnabledButton);

    auto initLabel = [] (juce::Label& l, const juce::String& s)
    {
        l.setText (s, juce::dontSendNotification);
        l.setJustificationType (juce::Justification::centred);
        l.setColour (juce::Label::textColourId, juce::Colours::white.withMultipliedAlpha (0.85f));
        l.setFont (juce::Font (12.0f));
    };

    initLabel (t0Label, "T0");
    initLabel (t1Label, "T1");
    initLabel (t2Label, "T2");
    initLabel (t3Label, "T3");

    addAndMakeVisible (t0Label);
    addAndMakeVisible (t1Label);
    addAndMakeVisible (t2Label);
    addAndMakeVisible (t3Label);

    // Configure sliders from param ranges
    using namespace levelscope::mtdm::ParamIDs;

    configureSliderForParam (t0Slider, t0Lufs, juce::Slider::LinearVertical, " LUFS");
    configureSliderForParam (t1Slider, t1Lufs, juce::Slider::LinearVertical, " LUFS");
    configureSliderForParam (t2Slider, t2Lufs, juce::Slider::LinearVertical, " LUFS");
    configureSliderForParam (t3Slider, t3Lufs, juce::Slider::LinearVertical, " LUFS");

    addAndMakeVisible (t0Slider);
    addAndMakeVisible (t1Slider);
    addAndMakeVisible (t2Slider);
    addAndMakeVisible (t3Slider);

    // GR meters
    limGrMeter.setNameLabel ("Limiter");
    downGrMeter.setNameLabel ("Downward");

    addAndMakeVisible (limGrMeter);
    addAndMakeVisible (downGrMeter);

    // Attachments
    mtdmEnabledAtt = std::make_unique<ButtonAttachment> (apvts, enabled,       mtdmEnabledButton);
    downEnabledAtt = std::make_unique<ButtonAttachment> (apvts, downEnabled01, downEnabledButton);
    limEnabledAtt  = std::make_unique<ButtonAttachment> (apvts, limEnabled01,  limEnabledButton);

    t0Att = std::make_unique<SliderAttachment> (apvts, t0Lufs, t0Slider);
    t1Att = std::make_unique<SliderAttachment> (apvts, t1Lufs, t1Slider);
    t2Att = std::make_unique<SliderAttachment> (apvts, t2Lufs, t2Slider);
    t3Att = std::make_unique<SliderAttachment> (apvts, t3Lufs, t3Slider);

    // [BEGIN MTDM-PANEL-THRESH-ORDER-HOOKS]
    auto attachThreshHooks = [this] (juce::Slider& s, int idx)
    {
        s.onDragStart = [this, idx]
        {
            thresholdSliderDragging = true;
            thresholdSliderDraggingIndex = idx;
            pushedGestureActive = { { false, false, false, false } };
        };

        s.onValueChange = [this, idx]
        {
            enforceThresholdOrderingFromSlider (idx);
        };

        s.onDragEnd = [this]
        {
            endPushedThresholdGestures();
            thresholdSliderDragging = false;
            thresholdSliderDraggingIndex = -1;
        };
    };

    attachThreshHooks (t0Slider, 0);
    attachThreshHooks (t1Slider, 1);
    attachThreshHooks (t2Slider, 2);
    attachThreshHooks (t3Slider, 3);
    // [END MTDM-PANEL-THRESH-ORDER-HOOKS]

    // UI poll rate for meters
    startTimerHz (30);
}

MtdmControlPanel::~MtdmControlPanel()
{
    stopTimer();
}

void MtdmControlPanel::configureToggle (juce::ToggleButton& b, const juce::String& text)
{
    b.setButtonText (text);
    b.setColour (juce::ToggleButton::textColourId, juce::Colours::white.withMultipliedAlpha (0.9f));
}

void MtdmControlPanel::configureSliderForParam (juce::Slider& s,
                                                const juce::String& paramID,
                                                juce::Slider::SliderStyle style,
                                                const juce::String& suffix)
{
    s.setSliderStyle (style);
    s.setTextBoxStyle (juce::Slider::TextBoxBelow, false, 64, 18);
    s.setColour (juce::Slider::textBoxTextColourId, juce::Colours::white.withMultipliedAlpha (0.9f));
    s.setColour (juce::Slider::textBoxOutlineColourId, juce::Colours::white.withMultipliedAlpha (0.15f));
    s.setTextValueSuffix (suffix);

    // [BEGIN MTDM-PANEL-SLIDER-RANGE-DOUBLE-FIX]
    if (auto* p = apvts.getParameter (paramID))
    {
        const auto rf = p->getNormalisableRange(); // NormalisableRange<float>

        juce::NormalisableRange<double> rd ((double) rf.start,
                                            (double) rf.end,
                                            (double) rf.interval);

        rd.skew = (double) rf.skew;
        rd.symmetricSkew = rf.symmetricSkew;

        s.setNormalisableRange (rd);
    }
    // [END MTDM-PANEL-SLIDER-RANGE-DOUBLE-FIX]
}

// [BEGIN MTDM-PANEL-THRESH-ORDER-IMPL]
void MtdmControlPanel::endPushedThresholdGestures()
{
    using namespace levelscope::mtdm::ParamIDs;

    const juce::String ids[4] = { t0Lufs, t1Lufs, t2Lufs, t3Lufs };

    for (int i = 0; i < 4; ++i)
    {
        if (! pushedGestureActive[(size_t) i])
            continue;

        if (auto* p = apvts.getParameter (ids[i]))
            p->endChangeGesture();

        pushedGestureActive[(size_t) i] = false;
    }
}

void MtdmControlPanel::enforceThresholdOrderingFromSlider (int changedIndex)
{
    if (thresholdCallbacksSuppressed)
        return;

    if (changedIndex < 0 || changedIndex > 3)
        return;

    using namespace levelscope::mtdm::ParamIDs;
    const juce::String ids[4] = { t0Lufs, t1Lufs, t2Lufs, t3Lufs };

    juce::RangedAudioParameter* params[4] =
    {
        apvts.getParameter (ids[0]),
        apvts.getParameter (ids[1]),
        apvts.getParameter (ids[2]),
        apvts.getParameter (ids[3])
    };

    if (params[0] == nullptr || params[1] == nullptr || params[2] == nullptr || params[3] == nullptr)
        return;

    auto clampSnap = [&] (int idx, float v) -> float
    {
        const auto r = params[idx]->getNormalisableRange();
        v = juce::jlimit ((float) r.start, (float) r.end, v);
        v = r.snapToLegalValue (v);
        return v;
    };

    // Read current parameter values (authoritative)
    float v[4];
    for (int i = 0; i < 4; ++i)
    {
        // convertFrom0to1(getValue()) gives the parameter's current "real" value.
        v[i] = (float) params[i]->convertFrom0to1 (params[i]->getValue());
        v[i] = clampSnap (i, v[i]);
    }

    // Anchor = the value the user is changing (from the slider, so it matches what they see)
    const juce::Slider* sliders[4] = { &t0Slider, &t1Slider, &t2Slider, &t3Slider };
    v[changedIndex] = clampSnap (changedIndex, (float) sliders[changedIndex]->getValue());

    // Push OUTWARDS from the changed index (same UX as handle dragging)
    for (int i = changedIndex + 1; i < 4; ++i)
    {
        const float minAllowed = v[i - 1] + minGapLu;
        if (v[i] < minAllowed)
            v[i] = clampSnap (i, minAllowed);
    }

    for (int i = changedIndex - 1; i >= 0; --i)
    {
        const float maxAllowed = v[i + 1] - minGapLu;
        if (v[i] > maxAllowed)
            v[i] = clampSnap (i, maxAllowed);
    }

    // Apply pushed neighbors only (do not fight the actively changed slider)
    thresholdCallbacksSuppressed = true;

    for (int i = 0; i < 4; ++i)
    {
        if (i == changedIndex)
            continue;

        const float cur = (float) params[i]->convertFrom0to1 (params[i]->getValue());
        const float dst = v[i];

        if (std::abs (dst - cur) < 1.0e-6f)
            continue;

        // Gesture strategy:
        // - while user drags a threshold slider, keep neighbor gestures open (lazy-start)
        // - otherwise (text entry / programmatic), do an immediate begin-set-end
        if (thresholdSliderDragging)
        {
            if (! pushedGestureActive[(size_t) i])
            {
                params[i]->beginChangeGesture();
                pushedGestureActive[(size_t) i] = true;
            }

            params[i]->setValueNotifyingHost (juce::jlimit (0.0f, 1.0f, params[i]->convertTo0to1 (dst)));
        }
        else
        {
            params[i]->beginChangeGesture();
            params[i]->setValueNotifyingHost (juce::jlimit (0.0f, 1.0f, params[i]->convertTo0to1 (dst)));
            params[i]->endChangeGesture();
        }
    }

    thresholdCallbacksSuppressed = false;
}
// [END MTDM-PANEL-THRESH-ORDER-IMPL]

void MtdmControlPanel::timerCallback()
{
    // Poll DSP snapshots (message thread)
    const auto lim  = processor.getLimiterMeteringSnapshot();
    const auto down = processor.getDownwardMeteringSnapshot();

    limGrMeter.setValuesDb (lim.grDbCurrent, lim.grDbHold);
    downGrMeter.setValuesDb (down.grDbCurrent, down.grDbHold);

    limGrMeter.repaint();
    downGrMeter.repaint();
}

void MtdmControlPanel::paint (juce::Graphics& g)
{
    g.fillAll (juce::Colour::fromRGB (12, 22, 36));
    g.setColour (juce::Colours::white.withMultipliedAlpha (0.08f));
    g.drawRect (getLocalBounds());
}

void MtdmControlPanel::resized()
{
    auto r = getLocalBounds().reduced (8);

    // Top row: toggles
    {
        auto top = r.removeFromTop (26);
        mtdmEnabledButton.setBounds (top.removeFromLeft (90));
        downEnabledButton.setBounds (top.removeFromLeft (90));
        limEnabledButton.setBounds  (top.removeFromLeft (90));
    }

    r.removeFromTop (6);

    // Bottom: meters
    auto meters = r.removeFromBottom (52);
    limGrMeter.setBounds  (meters.removeFromTop (24));
    meters.removeFromTop (4);
    downGrMeter.setBounds (meters.removeFromTop (24));

    r.removeFromBottom (6);

    // Middle: 4 threshold sliders (vertical)
    auto slidersArea = r;
    const int colW = juce::jmax (60, slidersArea.getWidth() / 4);

    auto col0 = slidersArea.removeFromLeft (colW);
    auto col1 = slidersArea.removeFromLeft (colW);
    auto col2 = slidersArea.removeFromLeft (colW);
    auto col3 = slidersArea;

    auto place = [] (juce::Rectangle<int> col, juce::Label& lab, juce::Slider& s)
    {
        lab.setBounds (col.removeFromTop (16));
        col.removeFromTop (4);
        s.setBounds (col);
    };

    place (col0, t0Label, t0Slider);
    place (col1, t1Label, t1Slider);
    place (col2, t2Label, t2Slider);
    place (col3, t3Label, t3Slider);
}
// [END MTDM-PANEL-IMPL]

LevelScopeAudioProcessorEditor::LevelScopeAudioProcessorEditor (LevelScopeAudioProcessor& p)
    : AudioProcessorEditor (&p),
            historyComponent (p),
      // [BEGIN UI3A-EDITOR-CTOR-MISSIONCONTROL]
      missionControl (p, historyComponent),
      // [END UI3A-EDITOR-CTOR-MISSIONCONTROL]
      // [BEGIN MTDM-PANEL-EDITOR-CTOR-INIT]
      mtdmPanel (p),
      layoutResizerBar (&layoutHistoryAndPanel, 1, false) // false = horizontal bar (drag up/down)
      // [END MTDM-PANEL-EDITOR-CTOR-INIT]
{
    // [BEGIN UI3A-EDITOR-ADD-MISSIONCONTROL]
    addAndMakeVisible (missionControl);
    addAndMakeVisible (historyComponent);
    // [END UI3A-EDITOR-ADD-MISSIONCONTROL]

    // [BEGIN MTDM-PANEL-EDITOR-ADD]
    addAndMakeVisible (layoutResizerBar);
    addAndMakeVisible (mtdmPanel);

    // Stretchable layout: [0]=history, [1]=resizer, [2]=panel
    layoutHistoryAndPanel.setItemLayout (0, 140, -1.0, -0.75); // history: min 140px, preferred 75%
    layoutHistoryAndPanel.setItemLayout (1, 6,   10,   8);     // resizer bar: fixed-ish
    layoutHistoryAndPanel.setItemLayout (2, 140, -1.0, -0.25); // panel:  min 140px, preferred 25%
    // [END MTDM-PANEL-EDITOR-ADD]

    setResizable (true, true);
    // [BEGIN MTDM-PANEL-EDITOR-RESIZE-LIMITS]
    setResizeLimits (520, 320, 4096, 2048);
    // [END MTDM-PANEL-EDITOR-RESIZE-LIMITS]
    // [BEGIN MTDM-PANEL-EDITOR-DEFAULT-SIZE]
    setSize (900, 520);
    // [END MTDM-PANEL-EDITOR-DEFAULT-SIZE]
}

LevelScopeAudioProcessorEditor::~LevelScopeAudioProcessorEditor() = default;

//==============================================================================

void LevelScopeAudioProcessorEditor::paint (juce::Graphics& g)
{
    g.fillAll (juce::Colours::black);
}

void LevelScopeAudioProcessorEditor::resized()
{
    auto r = getLocalBounds();
    // [BEGIN UI3A-EDITOR-RESIZED-TOPSTRIP]
    const int statsH = missionControl.getPreferredHeight();
    missionControl.setBounds (r.removeFromTop (statsH));
    // [END UI3A-EDITOR-RESIZED-TOPSTRIP]

    // [BEGIN MTDM-PANEL-EDITOR-RESIZED-SPLIT]
    juce::Component* comps[] = { &historyComponent, &layoutResizerBar, &mtdmPanel };

    layoutHistoryAndPanel.layOutComponents (comps,
                                           3,
                                           r.getX(), r.getY(), r.getWidth(), r.getHeight(),
                                           true,   // vertically stacked
                                           true);  // resize other dimension
    // [END MTDM-PANEL-EDITOR-RESIZED-SPLIT]
}
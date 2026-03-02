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

    // [BEGIN UI3A-ROUTINGGRAPHIC-ROWLABELS-FIX]
    const float rowH = r.getHeight() * 0.5f;

    const float labelW = 30.0f; // reserved inside the component for "Det"/"App"
    const float cellsW = juce::jmax (1.0f, r.getWidth() - labelW);
    const float cellW  = (numCh > 0 ? cellsW / (float) numCh : cellsW);

    auto drawRow = [&] (float y0, const juce::String& rowName, auto isActiveFn)
    {
        // Row label inside reserved left strip
        g.setColour (juce::Colours::white.withMultipliedAlpha (0.70f));
        g.setFont (12.0f);

        juce::Rectangle<int> lbl ((int) r.getX(), (int) y0, (int) labelW, (int) rowH);
        g.drawText (rowName, lbl, juce::Justification::centred);

        // Cells start after label strip
        const float x0 = r.getX() + labelW;

        for (int ch = 0; ch < numCh; ++ch)
        {
            const auto t = layout.getTypeOfChannel (ch);
            const bool active = isActiveFn (t);

            juce::Rectangle<float> cell (x0 + cellW * (float) ch + 1.0f,
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
    // [END UI3A-ROUTINGGRAPHIC-ROWLABELS-FIX]

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

    // [BEGIN UI3C-MISSIONCONTROL-SETUP-FOLLOW-TOGGLE]
    setupToggle (toggleMomentary);
    setupToggle (toggleShortTerm);
    setupToggle (toggleGate);
    setupToggle (toggleRolling);
    setupToggle (toggleFollow);
    // [END UI3C-MISSIONCONTROL-SETUP-FOLLOW-TOGGLE]

    // [BEGIN UI3C-MISSIONCONTROL-ADD-FOLLOW-TOGGLE]
    addAndMakeVisible (toggleMomentary);
    addAndMakeVisible (toggleShortTerm);
    addAndMakeVisible (toggleGate);
    addAndMakeVisible (toggleRolling);
    addAndMakeVisible (toggleFollow);
    // [END UI3C-MISSIONCONTROL-ADD-FOLLOW-TOGGLE]

    // [BEGIN UI3C3-MISSIONCONTROL-RLRA-WINDOW-CTOR]
    rollingWindowBox.addItem ("30s",  30);
    rollingWindowBox.addItem ("60s",  60);
    rollingWindowBox.addItem ("120s", 120);

    rollingWindowBox.setJustificationType (juce::Justification::centred);
    rollingWindowBox.setColour (juce::ComboBox::textColourId, juce::Colours::white.withMultipliedAlpha (0.90f));
    rollingWindowBox.setColour (juce::ComboBox::outlineColourId, juce::Colours::white.withMultipliedAlpha (0.18f));
    rollingWindowBox.setColour (juce::ComboBox::backgroundColourId, juce::Colours::black.withMultipliedAlpha (0.15f));

    addAndMakeVisible (rollingWindowBox);

    rollingWindowBox.onChange = [this]
    {
        const int s = rollingWindowBox.getSelectedId();
        if (s == 30 || s == 60 || s == 120)
            processor.setRollingLraWindowSeconds (s);
    };
    // [END UI3C3-MISSIONCONTROL-RLRA-WINDOW-CTOR]

    // [BEGIN UI3C-MISSIONCONTROL-WIRE-FOLLOW]
    toggleMomentary.onClick = [this] { history.setShowMomentaryCurve (toggleMomentary.getToggleState()); };
    toggleShortTerm.onClick = [this] { history.setShowShortTermCurve (toggleShortTerm.getToggleState()); };
    toggleGate.onClick      = [this] { history.setShowGateCurve      (toggleGate.getToggleState()); };
    toggleRolling.onClick   = [this] { history.setShowRollingLraLane (toggleRolling.getToggleState()); };
    toggleFollow.onClick    = [this] { history.setFollowRightEdge    (toggleFollow.getToggleState()); };
    // [END UI3C-MISSIONCONTROL-WIRE-FOLLOW]

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

// [BEGIN UI3C-MISSIONCONTROL-REFRESH-WITH-FOLLOW]
// [BEGIN UI3C3-MISSIONCONTROL-REFRESH-RLRA-WINDOW]
void MissionControlComponent::refreshCurveToggleStatesFromHistory()
{
    toggleMomentary.setToggleState (history.getShowMomentaryCurve(), juce::dontSendNotification);
    toggleShortTerm.setToggleState (history.getShowShortTermCurve(), juce::dontSendNotification);
    toggleGate.setToggleState      (history.getShowGateCurve(), juce::dontSendNotification);
    toggleRolling.setToggleState   (history.getShowRollingLraLane(), juce::dontSendNotification);
    toggleFollow.setToggleState    (history.getFollowRightEdge(), juce::dontSendNotification);

    // Sync rolling window selector from processor state
    const int win = processor.getRollingLraWindowSeconds();
    if (win == 30 || win == 60 || win == 120)
        rollingWindowBox.setSelectedId (win, juce::dontSendNotification);

    // Disable window selector when lane is hidden (optional UX)
    rollingWindowBox.setEnabled (toggleRolling.getToggleState());
}
// [END UI3C3-MISSIONCONTROL-REFRESH-RLRA-WINDOW]
// [END UI3C-MISSIONCONTROL-REFRESH-WITH-FOLLOW]

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

// [BEGIN UI3A-MISSIONCONTROL-RESIZED-RELAYOUT]
void MissionControlComponent::resized()
{
    auto r = getLocalBounds().reduced (8);

    // Right: policy row + routing graphic + curve toggles under it
    const int rightW = juce::jlimit (320, 560, r.getWidth() / 2);
    auto right = r.removeFromRight (rightW);

    // Put curve toggles under the routing graphic (right column bottom)
    // [BEGIN UI3C3-MISSIONCONTROL-TOGGLES-LAYOUT-6]
    auto rightToggles = right.removeFromBottom (22);
    const int gap = 4;

    auto take = [&] (int w)
    {
        auto a = rightToggles.removeFromLeft (w);
        rightToggles.removeFromLeft (gap);
        return a;
    };

    // Compact fixed widths so everything fits nicely:
    toggleMomentary.setBounds (take (48));
    toggleShortTerm.setBounds (take (48));
    toggleGate.setBounds      (take (60));
    toggleRolling.setBounds   (take (60));

    rollingWindowBox.setBounds (take (66)); // 30s/60s/120s

    // Follow takes remaining space
    toggleFollow.setBounds (rightToggles);
    // [END UI3C3-MISSIONCONTROL-TOGGLES-LAYOUT-6]

    right.removeFromBottom (6);

    // Policy row at top of right column
    auto policyRow = right.removeFromTop (22);
    mcPolicyBox.setBounds    (policyRow.removeFromLeft (130));
    dialogDetBox.setBounds   (policyRow.removeFromLeft (60));
    dialogApplyBox.setBounds (policyRow.removeFromLeft (60));
    lfeDetToggle.setBounds   (policyRow.removeFromLeft (80));
    lfeApplyToggle.setBounds (policyRow.removeFromLeft (90));

    right.removeFromTop (6);

    // Remaining right area is routing graphic
    routingGraphic.setBounds (right);

    // Left: presets + targets/current aligned (now has room because toggles moved right)
    auto left = r;

    auto topRow = left.removeFromTop (22);
    savePresetButton.setBounds (topRow.removeFromLeft (120));
    topRow.removeFromLeft (6);
    loadPresetButton.setBounds (topRow.removeFromLeft (120));

    left.removeFromTop (6);

    // Targets/current table with headers
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

    // Current row (this was missing for you; it will now be visible)
    auto currentRow = left.removeFromTop (26);
    hdrCurrentRow.setBounds (currentRow.removeFromLeft (rowHdrW));
    currentILabel.setBounds    (currentRow.removeFromLeft (colW));
    currentPeakLabel.setBounds (currentRow.removeFromLeft (colW));
    currentLraLabel.setBounds  (currentRow.removeFromLeft (colW));
}
// [END UI3A-MISSIONCONTROL-RESIZED-RELAYOUT]

//==============================================================================
// [BEGIN UI4A-MTDM-PANEL-CARDS-IMPL]
static void setSliderRangeFromParam (juce::AudioProcessorValueTreeState& apvts,
                                     const juce::String& paramID,
                                     juce::Slider& s)
{
    if (auto* p = apvts.getParameter (paramID))
    {
        const auto rf = p->getNormalisableRange(); // float
        juce::NormalisableRange<double> rd ((double) rf.start, (double) rf.end, (double) rf.interval);
        rd.skew = (double) rf.skew;
        rd.symmetricSkew = rf.symmetricSkew;
        s.setNormalisableRange (rd);
    }
}

static void styleSlider (juce::Slider& s, const juce::String& suffix)
{
    s.setSliderStyle (juce::Slider::LinearHorizontal);
    s.setTextBoxStyle (juce::Slider::TextBoxRight, false, 80, 18);
    s.setTextValueSuffix (suffix);

    s.setColour (juce::Slider::textBoxTextColourId, juce::Colours::white.withMultipliedAlpha (0.90f));
    s.setColour (juce::Slider::textBoxOutlineColourId, juce::Colours::white.withMultipliedAlpha (0.15f));
}

static void styleLabel (juce::Label& l, const juce::String& text)
{
    l.setText (text, juce::dontSendNotification);
    l.setColour (juce::Label::textColourId, juce::Colours::white.withMultipliedAlpha (0.75f));
    l.setFont (juce::Font (12.0f));
    l.setJustificationType (juce::Justification::centredLeft);
}

//==============================================================================
// Card base
//==============================================================================

MtdmCardComponent::MtdmCardComponent (juce::String titleText)
{
    title.setText (std::move (titleText), juce::dontSendNotification);
    title.setColour (juce::Label::textColourId, juce::Colours::white.withMultipliedAlpha (0.92f));
    title.setFont (juce::Font (14.0f, juce::Font::bold));
    addAndMakeVisible (title);

    setOpaque (false);
}

void MtdmCardComponent::paint (juce::Graphics& g)
{
    auto r = getLocalBounds().toFloat();
    g.setColour (juce::Colours::black.withMultipliedAlpha (0.18f));
    g.fillRoundedRectangle (r, 8.0f);

    g.setColour (juce::Colours::white.withMultipliedAlpha (0.10f));
    g.drawRoundedRectangle (r, 8.0f, 1.0f);
}

// [BEGIN UI4A3-CARD-CONTENTAREA-TITLE-OPTIONAL]
juce::Rectangle<int> MtdmCardComponent::getContentArea() const
{
    auto r = getLocalBounds().reduced (10);

    const int headerH = (title.isVisible() ? 24 : 0);
    if (headerH > 0)
        r.removeFromTop (headerH);

    return r;
}
// [END UI4A3-CARD-CONTENTAREA-TITLE-OPTIONAL]

// [BEGIN UI4A3-CARD-RESIZED-TITLE-OPTIONAL]
void MtdmCardComponent::resized()
{
    if (title.isVisible())
        title.setBounds (getLocalBounds().reduced (10, 6).removeFromTop (18));
}
// [END UI4A3-CARD-RESIZED-TITLE-OPTIONAL]

//==============================================================================
// Levelling placeholder
//==============================================================================

LevellingCard::LevellingCard()
    : MtdmCardComponent ("Levelling (Coming soon)")
{
    info.setText ("Placeholder for gain-riding module.", juce::dontSendNotification);
    info.setColour (juce::Label::textColourId, juce::Colours::white.withMultipliedAlpha (0.60f));
    info.setFont (juce::Font (12.0f));
    addAndMakeVisible (info);
}

void LevellingCard::resized()
{
    MtdmCardComponent::resized();
    info.setBounds (getContentArea().removeFromTop (18));
}

//==============================================================================
// Zones / Thresholds
//==============================================================================

MtdmZonesCard::MtdmZonesCard (LevelScopeAudioProcessor& p)
    : MtdmCardComponent ("Zones / Thresholds"),
      processor (p),
      apvts (p.getAPVTS())
{
    using namespace levelscope::mtdm::ParamIDs;

    addAndMakeVisible (mtdmEnabledButton);
    mtdmEnabledButton.setColour (juce::ToggleButton::textColourId, juce::Colours::white.withMultipliedAlpha (0.90f));
    mtdmEnabledAtt = std::make_unique<ButtonAttachment> (apvts, enabled, mtdmEnabledButton);

    styleLabel (t0Label, "T0 (LUFS)"); styleLabel (t1Label, "T1 (LUFS)");
    styleLabel (t2Label, "T2 (LUFS)"); styleLabel (t3Label, "T3 (LUFS)");
    addAndMakeVisible (t0Label); addAndMakeVisible (t1Label);
    addAndMakeVisible (t2Label); addAndMakeVisible (t3Label);

    styleSlider (t0Slider, " LUFS"); styleSlider (t1Slider, " LUFS");
    styleSlider (t2Slider, " LUFS"); styleSlider (t3Slider, " LUFS");

    setSliderRangeFromParam (apvts, t0Lufs, t0Slider);
    setSliderRangeFromParam (apvts, t1Lufs, t1Slider);
    setSliderRangeFromParam (apvts, t2Lufs, t2Slider);
    setSliderRangeFromParam (apvts, t3Lufs, t3Slider);

    addAndMakeVisible (t0Slider); addAndMakeVisible (t1Slider);
    addAndMakeVisible (t2Slider); addAndMakeVisible (t3Slider);

    t0Att = std::make_unique<SliderAttachment> (apvts, t0Lufs, t0Slider);
    t1Att = std::make_unique<SliderAttachment> (apvts, t1Lufs, t1Slider);
    t2Att = std::make_unique<SliderAttachment> (apvts, t2Lufs, t2Slider);
    t3Att = std::make_unique<SliderAttachment> (apvts, t3Lufs, t3Slider);

    auto attachHooks = [this] (juce::Slider& s, int idx)
    {
        s.onDragStart = [this]
        {
            thresholdSliderDragging = true;
            pushedGestureActive = { { false, false, false, false } };
        };

        s.onValueChange = [this, idx] { enforceThresholdOrderingFromSlider (idx); };

        s.onDragEnd = [this]
        {
            endPushedThresholdGestures();
            thresholdSliderDragging = false;
        };
    };

    attachHooks (t0Slider, 0);
    attachHooks (t1Slider, 1);
    attachHooks (t2Slider, 2);
    attachHooks (t3Slider, 3);
}

MtdmZonesCard::~MtdmZonesCard()
{
    endPushedThresholdGestures();
}

void MtdmZonesCard::endPushedThresholdGestures()
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

void MtdmZonesCard::enforceThresholdOrderingFromSlider (int changedIndex)
{
    if (callbacksSuppressed || changedIndex < 0 || changedIndex > 3)
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

    if (! params[0] || ! params[1] || ! params[2] || ! params[3])
        return;

    auto clampSnap = [&] (int idx, float v) -> float
    {
        const auto r = params[idx]->getNormalisableRange();
        v = juce::jlimit ((float) r.start, (float) r.end, v);
        v = r.snapToLegalValue (v);
        return v;
    };

    const juce::Slider* sliders[4] = { &t0Slider, &t1Slider, &t2Slider, &t3Slider };

    float v[4];
    for (int i = 0; i < 4; ++i)
    {
        v[i] = (float) params[i]->convertFrom0to1 (params[i]->getValue());
        v[i] = clampSnap (i, v[i]);
    }

    v[changedIndex] = clampSnap (changedIndex, (float) sliders[changedIndex]->getValue());

    // push outward from anchor
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

    callbacksSuppressed = true;

    for (int i = 0; i < 4; ++i)
    {
        if (i == changedIndex)
            continue;

        const float cur = (float) params[i]->convertFrom0to1 (params[i]->getValue());
        const float dst = v[i];

        if (std::abs (dst - cur) < 1.0e-6f)
            continue;

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

    callbacksSuppressed = false;
}

void MtdmZonesCard::resized()
{
    MtdmCardComponent::resized();
    auto r = getContentArea();

    mtdmEnabledButton.setBounds (r.removeFromTop (22));

    // [BEGIN UI4A3-ZONES-COMPACT-HIDE]
    const bool compact = (r.getHeight() < 70);

    t0Label.setVisible (! compact);  t0Slider.setVisible (! compact);
    t1Label.setVisible (! compact);  t1Slider.setVisible (! compact);
    t2Label.setVisible (! compact);  t2Slider.setVisible (! compact);
    t3Label.setVisible (! compact);  t3Slider.setVisible (! compact);

    if (compact)
        return;
    // [END UI4A3-ZONES-COMPACT-HIDE]

    // [BEGIN UI4A3-ZONES-SECTION-GAP]
    r.removeFromTop (2);
    // [END UI4A3-ZONES-SECTION-GAP]

    auto row = [&] (juce::Label& lab, juce::Slider& s)
    {
        auto rr = r.removeFromTop (24);
        lab.setBounds (rr.removeFromLeft (90));
        s.setBounds (rr);
        // [BEGIN UI4A3-ROW-GAP-ZERO]
        r.removeFromTop (0);
        // [END UI4A3-ROW-GAP-ZERO]
    };

    row (t0Label, t0Slider);
    row (t1Label, t1Slider);
    row (t2Label, t2Slider);
    row (t3Label, t3Slider);
}

//==============================================================================
// Upward
//==============================================================================

MtdmUpwardCard::MtdmUpwardCard (LevelScopeAudioProcessor& p)
    : MtdmCardComponent ("Upward (Essentials)"),
      apvts (p.getAPVTS())
{
    using namespace levelscope::mtdm::ParamIDs;

    styleLabel (modeLabel, "Mode");
    addAndMakeVisible (modeLabel);
    addAndMakeVisible (modeBox);

    modeBox.addItemList (juce::StringArray { "Spectral", "Broadband" }, 1);
    modeAtt = std::make_unique<ComboAttachment> (apvts, upwardModeChoice, modeBox);

    styleLabel (amountLabel, "Amount");
    styleLabel (maxBoostLabel, "MaxBoost");
    styleLabel (attackLabel, "Attack");
    styleLabel (releaseLabel, "Release");

    addAndMakeVisible (amountLabel); addAndMakeVisible (maxBoostLabel);
    addAndMakeVisible (attackLabel); addAndMakeVisible (releaseLabel);

    styleSlider (amountSlider, "");
    styleSlider (maxBoostSlider, " dB");
    styleSlider (attackSlider, " ms");
    styleSlider (releaseSlider, " ms");

    setSliderRangeFromParam (apvts, sucAmount01, amountSlider);
    setSliderRangeFromParam (apvts, sucMaxBoostDb, maxBoostSlider);
    setSliderRangeFromParam (apvts, sucAttackMs, attackSlider);
    setSliderRangeFromParam (apvts, sucReleaseMs, releaseSlider);

    addAndMakeVisible (amountSlider); addAndMakeVisible (maxBoostSlider);
    addAndMakeVisible (attackSlider); addAndMakeVisible (releaseSlider);

    amountAtt   = std::make_unique<SliderAttachment> (apvts, sucAmount01, amountSlider);
    maxBoostAtt = std::make_unique<SliderAttachment> (apvts, sucMaxBoostDb, maxBoostSlider);
    attackAtt   = std::make_unique<SliderAttachment> (apvts, sucAttackMs, attackSlider);
    releaseAtt  = std::make_unique<SliderAttachment> (apvts, sucReleaseMs, releaseSlider);
}

void MtdmUpwardCard::resized()
{
    MtdmCardComponent::resized();
    auto r = getContentArea();

    // [BEGIN UI4A3-UPWARD-COMPACT-HIDE]
    const bool compact = (r.getHeight() < 70);

    modeLabel.setVisible (! compact);
    modeBox.setVisible (! compact);

    amountLabel.setVisible (! compact);  amountSlider.setVisible (! compact);
    maxBoostLabel.setVisible (! compact);maxBoostSlider.setVisible (! compact);
    attackLabel.setVisible (! compact);  attackSlider.setVisible (! compact);
    releaseLabel.setVisible (! compact); releaseSlider.setVisible (! compact);

    if (compact)
        return;
    // [END UI4A3-UPWARD-COMPACT-HIDE]

    auto rr = r.removeFromTop (24);
    modeLabel.setBounds (rr.removeFromLeft (90));
    modeBox.setBounds (rr);

    // [BEGIN UI4A3-UPWARD-SECTION-GAP]
    r.removeFromTop (2);
    // [END UI4A3-UPWARD-SECTION-GAP]

    auto row = [&] (juce::Label& lab, juce::Slider& s)
    {
        auto r2 = r.removeFromTop (24);
        lab.setBounds (r2.removeFromLeft (90));
        s.setBounds (r2);
        // [BEGIN UI4A3-ROW-GAP-ZERO-UP]
        r.removeFromTop (0);
        // [END UI4A3-ROW-GAP-ZERO-UP]
    };

    row (amountLabel, amountSlider);
    row (maxBoostLabel, maxBoostSlider);
    row (attackLabel, attackSlider);
    row (releaseLabel, releaseSlider);
}

//==============================================================================
// Downward
//==============================================================================

// [BEGIN UI4A3-DOWNWARD-HIDE-TITLE]
MtdmDownwardCard::MtdmDownwardCard (LevelScopeAudioProcessor& p)
    : MtdmCardComponent (""), // title hidden; enabled toggle acts as header
      apvts (p.getAPVTS())
{
    title.setVisible (false);
// [END UI4A3-DOWNWARD-HIDE-TITLE]
    using namespace levelscope::mtdm::ParamIDs;

    enabledButton.setColour (juce::ToggleButton::textColourId, juce::Colours::white.withMultipliedAlpha (0.90f));
    addAndMakeVisible (enabledButton);
    enabledAtt = std::make_unique<ButtonAttachment> (apvts, downEnabled01, enabledButton);

    styleLabel (ratioLabel, "Ratio");
    styleLabel (attackLabel, "Attack");
    styleLabel (releaseLabel, "Release");
    styleLabel (makeupLabel, "Makeup");

    addAndMakeVisible (ratioLabel); addAndMakeVisible (attackLabel);
    addAndMakeVisible (releaseLabel); addAndMakeVisible (makeupLabel);

    styleSlider (ratioSlider, "");
    styleSlider (attackSlider, " ms");
    styleSlider (releaseSlider, " ms");
    styleSlider (makeupSlider, " dB");

    setSliderRangeFromParam (apvts, downRatio, ratioSlider);
    setSliderRangeFromParam (apvts, downAttackMs, attackSlider);
    setSliderRangeFromParam (apvts, downReleaseMs, releaseSlider);
    setSliderRangeFromParam (apvts, downMakeupDb, makeupSlider);

    addAndMakeVisible (ratioSlider); addAndMakeVisible (attackSlider);
    addAndMakeVisible (releaseSlider); addAndMakeVisible (makeupSlider);

    ratioAtt   = std::make_unique<SliderAttachment> (apvts, downRatio, ratioSlider);
    attackAtt  = std::make_unique<SliderAttachment> (apvts, downAttackMs, attackSlider);
    releaseAtt = std::make_unique<SliderAttachment> (apvts, downReleaseMs, releaseSlider);
    makeupAtt  = std::make_unique<SliderAttachment> (apvts, downMakeupDb, makeupSlider);
}

void MtdmDownwardCard::resized()
{
    MtdmCardComponent::resized();
    auto r = getContentArea();

    enabledButton.setBounds (r.removeFromTop (22));

    // [BEGIN UI4A3-DOWNWARD-COMPACT-HIDE]
    const bool compact = (r.getHeight() < 60);

    ratioLabel.setVisible (! compact);   ratioSlider.setVisible (! compact);
    attackLabel.setVisible (! compact);  attackSlider.setVisible (! compact);
    releaseLabel.setVisible (! compact); releaseSlider.setVisible (! compact);
    makeupLabel.setVisible (! compact);  makeupSlider.setVisible (! compact);

    if (compact)
        return;
    // [END UI4A3-DOWNWARD-COMPACT-HIDE]

    // [BEGIN UI4A3-DOWNWARD-SECTION-GAP]
    r.removeFromTop (2);
    // [END UI4A3-DOWNWARD-SECTION-GAP]

    auto row = [&] (juce::Label& lab, juce::Slider& s)
    {
        auto r2 = r.removeFromTop (24);
        lab.setBounds (r2.removeFromLeft (90));
        s.setBounds (r2);
        // [BEGIN UI4A3-ROW-GAP-ZERO-DOWN]
        r.removeFromTop (0);
        // [END UI4A3-ROW-GAP-ZERO-DOWN]
    };

    row (ratioLabel, ratioSlider);
    row (attackLabel, attackSlider);
    row (releaseLabel, releaseSlider);
    row (makeupLabel, makeupSlider);
}

//==============================================================================
// Limiter
//==============================================================================

// [BEGIN UI4A3-LIMITER-HIDE-TITLE]
MtdmLimiterCard::MtdmLimiterCard (LevelScopeAudioProcessor& p)
    : MtdmCardComponent (""), // title hidden; enabled toggle acts as header
      apvts (p.getAPVTS())
{
    title.setVisible (false);
// [END UI4A3-LIMITER-HIDE-TITLE]
    using namespace levelscope::mtdm::ParamIDs;

    enabledButton.setColour (juce::ToggleButton::textColourId, juce::Colours::white.withMultipliedAlpha (0.90f));
    addAndMakeVisible (enabledButton);
    enabledAtt = std::make_unique<ButtonAttachment> (apvts, limEnabled01, enabledButton);

    styleLabel (ceilingLabel, "Ceiling");
    styleLabel (lookLabel, "Lookahead");
    styleLabel (osLabel, "OS");
    styleLabel (attackLabel, "Attack");
    styleLabel (releaseLabel, "Release");
    styleLabel (driveLabel, "Drive");

    addAndMakeVisible (ceilingLabel); addAndMakeVisible (lookLabel);
    addAndMakeVisible (osLabel); addAndMakeVisible (attackLabel);
    addAndMakeVisible (releaseLabel); addAndMakeVisible (driveLabel);

    styleSlider (ceilingSlider, " dBFS");
    styleSlider (lookSlider, " ms");
    styleSlider (attackSlider, " ms");
    styleSlider (releaseSlider, " ms");
    styleSlider (driveSlider, " dB");

    setSliderRangeFromParam (apvts, limCeilingDb, ceilingSlider);
    setSliderRangeFromParam (apvts, limLookaheadMs, lookSlider);
    setSliderRangeFromParam (apvts, limAttackMs, attackSlider);
    setSliderRangeFromParam (apvts, limReleaseMs, releaseSlider);
    setSliderRangeFromParam (apvts, limDriveDb, driveSlider);

    addAndMakeVisible (ceilingSlider);
    addAndMakeVisible (lookSlider);
    addAndMakeVisible (attackSlider);
    addAndMakeVisible (releaseSlider);
    addAndMakeVisible (driveSlider);

    ceilingAtt = std::make_unique<SliderAttachment> (apvts, limCeilingDb, ceilingSlider);
    lookAtt    = std::make_unique<SliderAttachment> (apvts, limLookaheadMs, lookSlider);
    attackAtt  = std::make_unique<SliderAttachment> (apvts, limAttackMs, attackSlider);
    releaseAtt = std::make_unique<SliderAttachment> (apvts, limReleaseMs, releaseSlider);
    driveAtt   = std::make_unique<SliderAttachment> (apvts, limDriveDb, driveSlider);

    addAndMakeVisible (osBox);
    osBox.addItemList (juce::StringArray { "Off", "2x", "4x" }, 1);
    osAtt = std::make_unique<ComboAttachment> (apvts, limOversamplingChoice, osBox);
}

void MtdmLimiterCard::resized()
{
    MtdmCardComponent::resized();
    auto r = getContentArea();

    enabledButton.setBounds (r.removeFromTop (22));

    // [BEGIN UI4A3-LIMITER-COMPACT-HIDE]
    const bool compact = (r.getHeight() < 60);

    ceilingLabel.setVisible (! compact); ceilingSlider.setVisible (! compact);
    lookLabel.setVisible (! compact);    lookSlider.setVisible (! compact);
    osLabel.setVisible (! compact);      osBox.setVisible (! compact);
    attackLabel.setVisible (! compact);  attackSlider.setVisible (! compact);
    releaseLabel.setVisible (! compact); releaseSlider.setVisible (! compact);
    driveLabel.setVisible (! compact);   driveSlider.setVisible (! compact);

    if (compact)
        return;
    // [END UI4A3-LIMITER-COMPACT-HIDE]

    r.removeFromTop (8);

    auto rowS = [&] (juce::Label& lab, juce::Slider& s)
    {
        auto r2 = r.removeFromTop (24);
        lab.setBounds (r2.removeFromLeft (90));
        s.setBounds (r2);
        // [BEGIN UI4A3-ROW-GAP-ZERO-LIM]
        r.removeFromTop (0);
        // [END UI4A3-ROW-GAP-ZERO-LIM]
    };

    rowS (ceilingLabel, ceilingSlider);
    rowS (lookLabel, lookSlider);

    // Oversampling row
    auto rOS = r.removeFromTop (24);
    osLabel.setBounds (rOS.removeFromLeft (90));
    osBox.setBounds (rOS);
    // [BEGIN UI4A3-ROW-GAP-ZERO-LIM-OS]
    r.removeFromTop (0);
    // [END UI4A3-ROW-GAP-ZERO-LIM-OS]

    rowS (attackLabel, attackSlider);
    rowS (releaseLabel, releaseSlider);
    rowS (driveLabel, driveSlider);
}

//==============================================================================
// Content stack
//==============================================================================

// [BEGIN UI4A1-CARDS-RESIZABLE-CTOR]
MtdmCardsContent::MtdmCardsContent (LevelScopeAudioProcessor& p)
    : levelling(),
      zones (p),
      upward (p),
      downward (p),
      limiter (p)
{
    addAndMakeVisible (levelling);
    addAndMakeVisible (bar01);
    addAndMakeVisible (zones);
    addAndMakeVisible (bar12);
    addAndMakeVisible (upward);
    addAndMakeVisible (bar23);
    addAndMakeVisible (downward);
    addAndMakeVisible (bar34);
    addAndMakeVisible (limiter);

    // Layout items: card, bar, card, bar, ...
    // setItemLayout (index, minSize, maxSize, preferredSize)
    // Bars have near-fixed size; cards have flexible size.
    cardsLayout.setItemLayout (0,  60, -1.0,  80); // levelling
    cardsLayout.setItemLayout (1,   6,  10.0,  8); // bar
    cardsLayout.setItemLayout (2,  34, -1.0, 220); // zones: allow collapse (title-only + maybe enabled)
    cardsLayout.setItemLayout (3,   6,  10.0,  8); // bar
    cardsLayout.setItemLayout (4,  34, -1.0, 220); // upward: allow collapse (title-only)
    cardsLayout.setItemLayout (5,   6,  10.0,  8); // bar
    cardsLayout.setItemLayout (6,  26, -1.0, 220); // downward: enabled toggle header
    cardsLayout.setItemLayout (7,   6,  10.0,  8); // bar
    cardsLayout.setItemLayout (8,  26, -1.0, 260); // limiter: enabled toggle header
}
// [END UI4A1-CARDS-RESIZABLE-CTOR]

// [BEGIN UI4A1-CARDS-PREFERRED-HEIGHT]
int MtdmCardsContent::getPreferredHeight() const noexcept
{
    // Used by the Viewport host to size the content.
    // Resizer bars redistribute height; total height comes from the viewport host.
    if (contentPreferredHeightPx > 0)
        return contentPreferredHeightPx;

    // Fallback initial value (before first resized)
    return 900;
}
// [END UI4A1-CARDS-PREFERRED-HEIGHT]

// [BEGIN UI4A1-CARDS-RESIZED-LAYOUT]
void MtdmCardsContent::resized()
{
    juce::Component* comps[] =
    {
        &levelling, &bar01,
        &zones,     &bar12,
        &upward,    &bar23,
        &downward,  &bar34,
        &limiter
    };

    cardsLayout.layOutComponents (comps,
                                 9,
                                 0, 0, getWidth(), getHeight(),
                                 true,  // vertically stacked
                                 true); // resize other dimension too

    contentPreferredHeightPx = juce::jmax (getHeight(), 1);
}
// [END UI4A1-CARDS-RESIZED-LAYOUT]

//==============================================================================
// Public panel component (viewport)
//==============================================================================

MtdmControlPanel::MtdmControlPanel (LevelScopeAudioProcessor& p)
    : content (p)
{
    setOpaque (true);

    viewport.setScrollBarsShown (true, false); // vertical only
    viewport.setViewedComponent (&content, false);

    addAndMakeVisible (viewport);
}

void MtdmControlPanel::paint (juce::Graphics& g)
{
    g.fillAll (juce::Colour::fromRGB (12, 22, 36));
}

void MtdmControlPanel::resized()
{
    viewport.setBounds (getLocalBounds());

    // Ensure content is tall enough to scroll; width matches viewport width.
    const int prefH = content.getPreferredHeight();
    content.setSize (juce::jmax (1, viewport.getWidth() - viewport.getScrollBarThickness()),
                     juce::jmax (prefH, viewport.getHeight()));
}
// [END UI4A-MTDM-PANEL-CARDS-IMPL]

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
    // [BEGIN UI4A-EDITOR-RESIZE-LIMITS-BIGGER]
    setResizeLimits (820, 520, 4096, 2048);
    // [END UI4A-EDITOR-RESIZE-LIMITS-BIGGER]
    // [BEGIN UI4A-EDITOR-DEFAULT-SIZE-BIGGER]
    setSize (1200, 760);
    // [END UI4A-EDITOR-DEFAULT-SIZE-BIGGER]
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
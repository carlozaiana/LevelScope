#include "PluginEditor.h"
// [BEGIN MTDM-PANEL-INCLUDE-PARAMIDS]
#include "Core/Processing/Modules/MultiThresholdDynamicsParamIDs.h"
// [END MTDM-PANEL-INCLUDE-PARAMIDS]

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
      statsComponent (p),
      historyComponent (p),
      // [BEGIN MTDM-PANEL-EDITOR-CTOR-INIT]
      mtdmPanel (p),
      layoutResizerBar (&layoutHistoryAndPanel, 1, false) // false = horizontal bar (drag up/down)
      // [END MTDM-PANEL-EDITOR-CTOR-INIT]
{
    addAndMakeVisible (statsComponent);
    addAndMakeVisible (historyComponent);

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
    const int statsH = statsComponent.getPreferredHeight();

    statsComponent.setBounds (r.removeFromTop (statsH));

    // [BEGIN MTDM-PANEL-EDITOR-RESIZED-SPLIT]
    juce::Component* comps[] = { &historyComponent, &layoutResizerBar, &mtdmPanel };

    layoutHistoryAndPanel.layOutComponents (comps,
                                           3,
                                           r.getX(), r.getY(), r.getWidth(), r.getHeight(),
                                           true,   // vertically stacked
                                           true);  // resize other dimension
    // [END MTDM-PANEL-EDITOR-RESIZED-SPLIT]
}
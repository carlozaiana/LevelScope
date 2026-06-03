#include "StandaloneMainComponent.h"

namespace levelscope::standalone
{
StandaloneMainComponent::StandaloneMainComponent()
{
    titleLabel.setText ("LevelScope Standalone", juce::dontSendNotification);
    titleLabel.setJustificationType (juce::Justification::centredLeft);
    titleLabel.setFont (juce::Font (28.0f, juce::Font::bold));
    addAndMakeVisible (titleLabel);

    subtitleLabel.setText ("Scaffold shell: proposal engine, file ingest, render/export not wired yet.",
                           juce::dontSendNotification);
    subtitleLabel.setJustificationType (juce::Justification::centredLeft);
    addAndMakeVisible (subtitleLabel);

    statusBox.setMultiLine (true);
    statusBox.setReadOnly (true);
    statusBox.setCaretVisible (false);
    statusBox.setScrollbarsShown (true);
    statusBox.setText (
        "Standalone scaffold status\n\n"
        "Implemented in this patch:\n"
        "- Opens a native JUCE app window.\n"
        "- Keeps the VST3 plugin target unchanged.\n"
        "- Links against LevelScopeCore without moving plugin code.\n\n"
        "Intentionally not implemented yet:\n"
        "- Audio file ingest\n"
        "- Offline analysis/render\n"
        "- Proposal engine\n"
        "- Export/report workflow\n",
        false);
    addAndMakeVisible (statusBox);

    setSize (980, 640);
}

void StandaloneMainComponent::paint (juce::Graphics& g)
{
    g.fillAll (juce::Colour (0xff15171c));

    g.setColour (juce::Colour (0xff2a2d35));
    g.drawRoundedRectangle (getLocalBounds().toFloat().reduced (12.0f), 10.0f, 1.0f);
}

void StandaloneMainComponent::resized()
{
    auto bounds = getLocalBounds().reduced (24);

    titleLabel.setBounds (bounds.removeFromTop (42));
    subtitleLabel.setBounds (bounds.removeFromTop (30));

    bounds.removeFromTop (18);
    statusBox.setBounds (bounds);
}
}
#include "LoudnessStatsComponent.h"
#include "PluginProcessor.h"

#include <cmath>

static void initLabel (juce::Label& l)
{
    l.setJustificationType (juce::Justification::centredLeft);
    l.setColour (juce::Label::textColourId, juce::Colours::white.withMultipliedAlpha (0.9f));
    l.setFont (juce::Font (14.0f));
}

juce::String LoudnessStatsComponent::formatLufs (float v)
{
    if (v <= -199.0f)
        return "-- LUFS";

    return juce::String (v, 1) + " LUFS";
}

juce::String LoudnessStatsComponent::formatLu (float v)
{
    if (! std::isfinite (v))
        return "-- LU";

    return juce::String (v, 1) + " LU";
}

LoudnessStatsComponent::LoudnessStatsComponent (LevelScopeAudioProcessor& p)
    : processor (p)
{
    initLabel (integratedLabel);
    initLabel (lraLabel);
    initLabel (rollingLabel);

    addAndMakeVisible (integratedLabel);
    addAndMakeVisible (lraLabel);
    addAndMakeVisible (rollingLabel);

    // Rolling window selector
    rollingWindowBox.addItem ("30s", 30);
    rollingWindowBox.addItem ("60s", 60);
    rollingWindowBox.addItem ("120s", 120);

    rollingWindowBox.setSelectedId (60, juce::dontSendNotification); // default
    processor.setRollingLraWindowSeconds (60);

    rollingWindowBox.onChange = [this]
    {
        const int s = rollingWindowBox.getSelectedId();
        if (s > 0)
            processor.setRollingLraWindowSeconds (s);
    };

    rollingWindowBox.setColour (juce::ComboBox::textColourId, juce::Colours::white.withMultipliedAlpha (0.9f));
    rollingWindowBox.setColour (juce::ComboBox::backgroundColourId, juce::Colours::black.withMultipliedAlpha (0.25f));
    rollingWindowBox.setColour (juce::ComboBox::outlineColourId, juce::Colours::white.withMultipliedAlpha (0.25f));

    addAndMakeVisible (rollingWindowBox);

    startTimerHz (5); // UI update: ~5 Hz
}

void LoudnessStatsComponent::resized()
{
    auto r = getLocalBounds().reduced (6, 2);

    auto row1 = r.removeFromTop (22);
    auto row2 = r.removeFromTop (22);

    {
        auto left = row1.removeFromLeft (row1.getWidth() / 2);
        integratedLabel.setBounds (left);
        lraLabel.setBounds (row1);
    }

    {
        auto boxArea = row2.removeFromRight (72);
        rollingWindowBox.setBounds (boxArea.reduced (2, 2));
        rollingLabel.setBounds (row2);
    }
}

void LoudnessStatsComponent::timerCallback()
{
    const float I = processor.getRunningIntegratedLufs();
    const float LRA = processor.getRunningLraLu();

    const int winS = processor.getRollingLraWindowSeconds();
    const float rolling = processor.getRollingLraLu();

    integratedLabel.setText ("I (running): " + formatLufs (I), juce::dontSendNotification);
    lraLabel.setText ("LRA (running): " + formatLu (LRA), juce::dontSendNotification);

    rollingLabel.setText ("LRA (" + juce::String (winS) + "s): " + formatLu (rolling), juce::dontSendNotification);

    // Keep UI selector and processor in sync if needed
    if (rollingWindowBox.getSelectedId() != winS)
        rollingWindowBox.setSelectedId (winS, juce::dontSendNotification);
}
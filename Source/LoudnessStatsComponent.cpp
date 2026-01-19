#include "LoudnessStatsComponent.h"
#include "PluginProcessor.h"

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

    addAndMakeVisible (integratedLabel);
    addAndMakeVisible (lraLabel);

    startTimerHz (5); // UI update rate: ~5 Hz
}

void LoudnessStatsComponent::resized()
{
    auto r = getLocalBounds().reduced (6, 2);

    auto left = r.removeFromLeft (r.getWidth() / 2);
    integratedLabel.setBounds (left);
    lraLabel.setBounds (r);
}

void LoudnessStatsComponent::timerCallback()
{
    const float I = processor.getRunningIntegratedLufs();
    const float LRA = processor.getRunningLraLu();

    integratedLabel.setText ("I (running): " + formatLufs (I), juce::dontSendNotification);
    lraLabel.setText ("LRA (running): " + formatLu (LRA), juce::dontSendNotification);
}
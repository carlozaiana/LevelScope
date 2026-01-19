#include "PluginEditor.h"

//==============================================================================

LevelScopeAudioProcessorEditor::LevelScopeAudioProcessorEditor (LevelScopeAudioProcessor& p)
    : AudioProcessorEditor (&p),
      statsComponent (p),
      historyComponent (p)
{
    addAndMakeVisible (statsComponent);
    addAndMakeVisible (historyComponent);

    setResizable (true, true);
    setResizeLimits (400, 200, 4096, 2048);
    setSize (800, 400);
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
    const int statsH = 28;

    statsComponent.setBounds (r.removeFromTop (statsH));
    historyComponent.setBounds (r);
}
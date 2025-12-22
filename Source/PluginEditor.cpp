#include "PluginEditor.h"

//==============================================================================

LevelScopeAudioProcessorEditor::LevelScopeAudioProcessorEditor (LevelScopeAudioProcessor& p)
    : AudioProcessorEditor (&p),
      historyComponent (p)
{
    addAndMakeVisible (historyComponent);

    // Initial size; user can resize freely
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
    historyComponent.setBounds (getLocalBounds());
}
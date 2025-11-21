#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"
#include "VolumeHistoryComponent.h"

//==============================================================================

class LevelScopeAudioProcessorEditor : public juce::AudioProcessorEditor
{
public:
    explicit LevelScopeAudioProcessorEditor (LevelScopeAudioProcessor&);
    ~LevelScopeAudioProcessorEditor() override;

    void paint (juce::Graphics&) override;
    void resized() override;

private:
    LevelScopeAudioProcessor& processor;
    VolumeHistoryComponent    historyComponent;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (LevelScopeAudioProcessorEditor)
};
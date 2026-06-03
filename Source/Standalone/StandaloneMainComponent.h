#pragma once

#include <juce_gui_basics/juce_gui_basics.h>

namespace levelscope::standalone
{
class StandaloneMainComponent final : public juce::Component
{
public:
    StandaloneMainComponent();

    void paint (juce::Graphics& g) override;
    void resized() override;

private:
    juce::Label titleLabel;
    juce::Label subtitleLabel;
    juce::TextEditor statusBox;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (StandaloneMainComponent)
};
}
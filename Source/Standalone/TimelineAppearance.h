#pragma once
#include <juce_gui_basics/juce_gui_basics.h>

namespace levelscope::standalone
{
// The existing plugin palette: history (16,30,50), controls (12,22,36),
// cards (10,18,28). Display preferences never enter programme identities.
namespace palette
{
inline juce::Colour graph(){return juce::Colour::fromRGB(16,30,50);}
inline juce::Colour panel(){return juce::Colour::fromRGB(12,22,36);}
inline juce::Colour card(){return juce::Colour::fromRGB(10,18,28);}
}
class TimelineBlend final : public juce::Slider
{
public:
    TimelineBlend()
    {
        setName("Waveform / loudness visibility");setSliderStyle(LinearVertical);
        setTextBoxStyle(NoTextBox,false,0,0);setRange(-1,1,.01);setValue(0,juce::dontSendNotification);
        setDoubleClickReturnValue(true,0);setScrollWheelEnabled(false);
        setTooltip("Up: loudness only. Centre: both. Down: waveform only. Double-click resets both. Display only.");
    }
    float waveformAlpha() const {return float(1-std::max(0.,getValue()));}
    float loudnessAlpha() const {return float(1+std::min(0.,getValue()));}
    void paint(juce::Graphics& g) override
    {
        const float x=getWidth()/2.f,top=8,bottom=std::max(top,float(getHeight()-8));
        g.setColour(juce::Colours::lightgrey.withAlpha(.5f));g.drawLine(x,top,x,bottom,1);
        g.drawLine(x-3,(top+bottom)/2,x+3,(top+bottom)/2,1);
        const float y=bottom-float((getValue()+1)/2)*(bottom-top);
        g.setColour(juce::Colours::white);g.drawLine(x-6,y,x+6,y,2);
    }
};
}

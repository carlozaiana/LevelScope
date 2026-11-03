#pragma once
#include <juce_gui_basics/juce_gui_basics.h>
namespace levelscope::standalone
{
struct TimelineLayout
{
    // The fixed 22 px information line now sits outside this height: 168+22
    // retains the old 190 px default. Travel is 120 px (previously 48).
    static constexpr int minimumComparison=90,defaultComparison=168,maximumComparison=210;
    static constexpr int informationHeight=22,dividerHeight=8,minimumTimeline=220,minimumInspector=290;
    int comparison=defaultComparison,inspector=minimumInspector;
    void moveComparison(int delta){comparison=juce::jlimit(minimumComparison,maximumComparison,comparison+delta);}
    void moveTimeline(int delta,int available,int channels=0,int laneHeight=0){inspector=juce::jlimit(minimumInspector,std::max(minimumInspector,available-std::max(minimumTimeline,137+16*channels)-laneHeight-dividerHeight),inspector-delta);}
    void reset(){comparison=defaultComparison;inspector=minimumInspector;}
};
class TimelineDivider final : public juce::Component
{
public:
    std::function<void(int)> move;
    explicit TimelineDivider(const juce::String& title){setMouseCursor(juce::MouseCursor::UpDownResizeCursor);setWantsKeyboardFocus(true);setTitle(title);}
    void paint(juce::Graphics& g) override {g.setColour(juce::Colour(0xff647487));g.fillRoundedRectangle(float(getWidth()/2-20),3.f,40.f,2.f,1.f);}
    void mouseDown(const juce::MouseEvent& e) override {lastY=e.getScreenY();}
    void mouseDrag(const juce::MouseEvent& e) override {const int y=e.getScreenY();if(move)move(y-lastY);lastY=y;}
    bool keyPressed(const juce::KeyPress& k) override {if(k==juce::KeyPress::upKey||k==juce::KeyPress::downKey){if(move)move(k==juce::KeyPress::upKey?-4:4);return true;}return false;}
private:int lastY=0;
};
}

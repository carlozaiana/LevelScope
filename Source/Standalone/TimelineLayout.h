#pragma once
#include <juce_gui_basics/juce_gui_basics.h>
#include <algorithm>
#include <cmath>
namespace levelscope::standalone
{
inline int64_t timelineFrameAt(double x,int left,int width,double begin,double end,int64_t frames)
{return width<=0?0:std::clamp<int64_t>(int64_t(std::floor(begin+std::clamp((x-left)/std::max(1,width-1),0.,1.)*(end-begin))),0,frames);}
struct TimelineLayout
{
    // The fixed 22 px information line now sits outside this height: 168+22
    // retains the old 190 px default. Travel is 120 px (previously 48).
    static constexpr int minimumComparison=90,defaultComparison=168,maximumComparison=210;
    static constexpr int informationHeight=22,dividerHeight=8,minimumTimeline=220,minimumInspector=150,defaultInspector=290;
    int comparison=defaultComparison,inspector=defaultInspector;
    void moveComparison(int delta){comparison=juce::jlimit(minimumComparison,maximumComparison,comparison+delta);}
    void moveTimeline(int delta,int available,int channels=0,int laneHeight=0){inspector=juce::jlimit(minimumInspector,std::max(minimumInspector,available-std::max(minimumTimeline,137+16*channels)-laneHeight-dividerHeight),inspector-delta);}
    void reset(){comparison=defaultComparison;inspector=defaultInspector;}
};
class TimelineDivider final : public juce::Component,public juce::TooltipClient
{
public:
    std::function<void(int)> move;
    std::function<void()> finished,reset;
    explicit TimelineDivider(const juce::String& title){setMouseCursor(juce::MouseCursor::UpDownResizeCursor);setWantsKeyboardFocus(true);setTitle(title);setName(title);}
    juce::String getTooltip() override {return getTitle()+". Drag vertically or use arrow keys."+(reset?" Double-click resets this height.":juce::String());}
    void paint(juce::Graphics& g) override {g.setColour(juce::Colour(0xff647487));g.fillRoundedRectangle(float(getWidth()/2-20),3.f,40.f,2.f,1.f);}
    void mouseDown(const juce::MouseEvent& e) override {lastY=e.getScreenY();grabKeyboardFocus();}
    void mouseDrag(const juce::MouseEvent& e) override {const int y=e.getScreenY();if(move)move(y-lastY);lastY=y;}
    void mouseUp(const juce::MouseEvent&) override {if(finished)finished();}
    void mouseDoubleClick(const juce::MouseEvent&) override {if(reset)reset();if(finished)finished();}
    bool keyPressed(const juce::KeyPress& k) override {if(k==juce::KeyPress::upKey||k==juce::KeyPress::downKey){if(move)move(k==juce::KeyPress::upKey?-4:4);if(finished)finished();return true;}return false;}
private:int lastY=0;
};
}

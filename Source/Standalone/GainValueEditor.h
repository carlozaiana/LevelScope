#pragma once
#include <juce_gui_basics/juce_gui_basics.h>
namespace levelscope::standalone
{
// Gain text is an unconfirmed entry until Enter. In particular a transport
// Space must not erase a selected number when focus returns from a dialog.
class GainValueEditor final : public juce::TextEditor
{
public:
    GainValueEditor(){juce::Desktop::getInstance().addGlobalMouseListener(&outside);}
    ~GainValueEditor() override {juce::Desktop::getInstance().removeGlobalMouseListener(&outside);}
    std::function<void()> onOutsideClick;
    bool keyPressed(const juce::KeyPress& key) override
    {
        if(key==juce::KeyPress::spaceKey)return true;
        return juce::TextEditor::keyPressed(key);
    }
private:
    struct OutsideClicks final : juce::MouseListener
    {
        GainValueEditor& editor;
        explicit OutsideClicks(GainValueEditor& e):editor(e){}
        void mouseDown(const juce::MouseEvent& e) override
        {if(editor.isShowing()&&e.originalComponent!=&editor&&!editor.isParentOf(e.originalComponent)&&editor.onOutsideClick)editor.onOutsideClick();}
    } outside{*this};
};
}

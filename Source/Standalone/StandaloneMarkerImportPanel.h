#pragma once
#include "StandaloneMarkerImport.h"
#include <juce_gui_basics/juce_gui_basics.h>
#include <future>
namespace levelscope::standalone
{
class MarkerImportWindow final : public juce::DocumentWindow
{
public:
    MarkerImportWindow(StandaloneSessionModel&,bool paste,std::function<bool()> canCommit,std::function<void()> changed);
    void closeButtonPressed() override;
private:
    class Panel;
    Panel* panel=nullptr;
};
}

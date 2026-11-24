#include <memory>

#include <juce_gui_basics/juce_gui_basics.h>

#include "StandaloneMainComponent.h"

namespace
{
class StandaloneMainWindow final : public juce::DocumentWindow
{
public:
    explicit StandaloneMainWindow (const juce::String& name)
        : juce::DocumentWindow (name,
                                levelscope::standalone::palette::panel(),
                                juce::DocumentWindow::allButtons)
    {
        setUsingNativeTitleBar (true);
        setResizable (true, true);
        setResizeLimits (1120, 950, 2400, 1600);
        auto* content=new levelscope::standalone::StandaloneMainComponent({},levelscope::standalone::OutputDeviceSettings(),{},levelscope::standalone::TimelineViewPreferences::defaultFile());setContentOwned(content,true);
        juce::PopupMenu appMenu;appMenu.addItem(levelscope::standalone::StandaloneMainComponent::preferences,juce::String::fromUTF8("Preferences… (⌘,)"));
        juce::MenuBarModel::setMacMainMenu(content,&appMenu);
        centreWithSize (getWidth(), getHeight());
        setVisible (true);
    }

    ~StandaloneMainWindow() override {juce::MenuBarModel::setMacMainMenu(nullptr);}

    void closeButtonPressed() override
    {
        juce::JUCEApplication::getInstance()->systemRequestedQuit();
    }

private:
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (StandaloneMainWindow)
};

class LevelScopeStandaloneApplication final : public juce::JUCEApplication
{
public:
    const juce::String getApplicationName() override
    {
        return JUCE_APPLICATION_NAME_STRING;
    }

    const juce::String getApplicationVersion() override
    {
        return JUCE_APPLICATION_VERSION_STRING;
    }

    bool moreThanOneInstanceAllowed() override
    {
        return true;
    }

    void initialise (const juce::String&) override
    {
        mainWindow = std::make_unique<StandaloneMainWindow> ("LevelScope-Timeline");
    }

    void shutdown() override
    {
        mainWindow = nullptr;
    }

    void systemRequestedQuit() override
    {
        if (auto* content = dynamic_cast<levelscope::standalone::StandaloneMainComponent*>(mainWindow->getContentComponent()))
            content->requestClose([] { quit(); });
        else quit();
    }

    void anotherInstanceStarted (const juce::String&) override
    {
    }

private:
    std::unique_ptr<StandaloneMainWindow> mainWindow;
};
}

START_JUCE_APPLICATION (LevelScopeStandaloneApplication)

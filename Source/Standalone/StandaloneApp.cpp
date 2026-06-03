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
                                juce::Colour (0xff15171c),
                                juce::DocumentWindow::allButtons)
    {
        setUsingNativeTitleBar (true);
        setResizable (true, true);
        setContentOwned (new levelscope::standalone::StandaloneMainComponent(), true);
        centreWithSize (getWidth(), getHeight());
        setVisible (true);
    }

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
        mainWindow = std::make_unique<StandaloneMainWindow> (getApplicationName());
    }

    void shutdown() override
    {
        mainWindow = nullptr;
    }

    void systemRequestedQuit() override
    {
        quit();
    }

    void anotherInstanceStarted (const juce::String&) override
    {
    }

private:
    std::unique_ptr<StandaloneMainWindow> mainWindow;
};
}

START_JUCE_APPLICATION (LevelScopeStandaloneApplication)

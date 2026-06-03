#pragma once

#include <memory>

#include <juce_gui_basics/juce_gui_basics.h>

#include "StandaloneSessionModel.h"

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-WORKFLOW-SHELL]

class StandaloneMainComponent final : public juce::Component
{
public:
    StandaloneMainComponent();

    void paint (juce::Graphics& g) override;
    void resized() override;

private:
    StandaloneSessionModel session;

    juce::Label titleLabel;
    juce::Label subtitleLabel;
    juce::Label workflowStateLabel;
    juce::Label targetProfileLabel;

    juce::ComboBox targetProfileBox;

    juce::TextButton importButton { "Import" };
    juce::TextButton analyzeButton { "Analyze" };
    juce::TextButton currentStateButton { "Current State" };
    juce::TextButton exportButton { "Export" };

    juce::TextButton chooseSourceButton { "Choose Source File..." };
    juce::TextButton clearSourceButton { "Clear Source" };
    juce::TextButton measureSourceButton { "Measure Source (not implemented)" };

    juce::TextEditor sourceStateBox;
    juce::TextEditor targetStateBox;
    juce::TextEditor currentStateBox;
    juce::TextEditor pageDetailBox;

    std::unique_ptr<juce::FileChooser> sourceFileChooser;

    void configureLabels();
    void configureTargetProfileBox();
    void configureNavigationButtons();
    void configureSourceControls();
    void configureReadOnlyBox (juce::TextEditor& box);

    void chooseSourceFile();
    void clearSourceFile();

    void setPage (WorkflowPage page);
    void refreshFromSession();

    static juce::String getSourceFileWildcard();

    juce::String getPageName (WorkflowPage page) const;
    juce::String buildSourceStateText() const;
    juce::String buildTargetStateText() const;
    juce::String buildCurrentStateText() const;
    juce::String buildPageDetailText() const;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (StandaloneMainComponent)
};

// [END LS-STANDALONE-WORKFLOW-SHELL]
}
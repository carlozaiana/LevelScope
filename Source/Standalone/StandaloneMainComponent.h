#pragma once

#include <memory>
#include <future>

#include <juce_gui_basics/juce_gui_basics.h>

#include "StandaloneSessionModel.h"
#include "StandaloneAudioPanel.h"

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-WORKFLOW-SHELL]

class StandaloneMainComponent final : public juce::Component, private juce::Timer
{
public:
    StandaloneMainComponent();
    ~StandaloneMainComponent() override;

    void paint (juce::Graphics& g) override;
    void resized() override;

private:
    StandaloneSessionModel session;
    StandaloneAudition audition;
    StandaloneAudioPanel audioPanel { session, audition, [this] { refreshFromSession(); } };
    std::shared_ptr<JobControl> importControl;
    std::future<MediaValidationResult> importResult;
    JobTicket importTicket;
    juce::File importingFile;
    juce::Label jobLabel;
    juce::TextButton cancelJobButton { "Cancel import" };

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

    juce::TextButton initializeCurrentStateButton { "Initialize Current From Source" };
    juce::TextButton clearCurrentStateButton { "Clear Current State" };
    juce::TextButton remeasureCurrentStateButton { "Re-measure Current State (not implemented)" };

    juce::TextEditor sourceStateBox;
    juce::TextEditor targetStateBox;
    juce::TextEditor currentStateBox;
    juce::TextEditor pageDetailBox;

    std::unique_ptr<juce::FileChooser> sourceFileChooser;

    void configureLabels();
    void configureTargetProfileBox();
    void configureNavigationButtons();
    void configureSourceControls();
    void configureCurrentStateControls();
    void configureReadOnlyBox (juce::TextEditor& box);

    void chooseSourceFile();
    void clearSourceFile();
    void initializeCurrentStateFromSource();
    void clearCurrentState();
    void beginImport (const juce::File&);
    void timerCallback() override;

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

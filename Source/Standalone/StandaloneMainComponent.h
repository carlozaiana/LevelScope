#pragma once
#include "SessionDocument.h"
#include "StandaloneProgrammeImport.h"

#include <memory>
#include <future>

#include <juce_gui_basics/juce_gui_basics.h>

#include "StandaloneSessionModel.h"
#include "StandaloneAudioPanel.h"
#include "StandaloneVideoPanel.h"
#include "StandaloneDiagnosticsPanel.h"
#include "StandaloneSectionsPanel.h"

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-WORKFLOW-SHELL]

class StandaloneMainComponent final : public juce::Component, private juce::Timer
{
public:
    StandaloneMainComponent();
    ~StandaloneMainComponent() override;
    void requestClose (std::function<void()>);

    void paint (juce::Graphics& g) override;
    void resized() override;

private:
    StandaloneSessionModel session;
    juce::TextButton newSessionButton { "New session" }, openSessionButton { "Open session..." },
        saveSessionButton { "Save" }, saveAsSessionButton { "Save as..." },
        relinkAudioButton { "Relink audio..." }, relinkVideoButton { "Relink video..." }, cancelSessionButton { "Cancel file job" };
    juce::Label sessionStatus;
    std::unique_ptr<juce::FileChooser> sessionChooser;
    std::future<sessionfile::Outcome> documentResult;
    std::shared_ptr<JobControl> documentControl;
    std::optional<sessionfile::Document> capturedSave;
    juce::File pendingDocumentFile;
    std::function<void()> afterSave;
    juce::String documentNotice;
    bool relinkingVideo = false;
    juce::String relinkingMember;
    bool migrationSaveApproved = false;
    enum class FileOperation { none, save, open, relink } fileOperation = FileOperation::none;
    void configureSessionControls();
    void pollSessionOperation();
    void refreshSessionControls();
    void confirmDiscard (std::function<void()>);
    void beginSaveSession (bool saveAs, std::function<void()> = {});
    void beginOpenSession (const juce::File&);
    void chooseRelink (bool video, juce::String slot = {});
    void installDocument (sessionfile::Document, bool relinking, std::shared_ptr<sessionfile::Prepared>);
    bool canReplaceSession() const;

    StandaloneAudition audition;
    StandaloneAudioPanel audioPanel { session, audition, [this] { refreshFromSession(); } };
    StandaloneVideoPanel videoPanel { session, audition };
    std::unique_ptr<juce::DocumentWindow> videoWindow;
    juce::TextButton videoButton { "Reference picture..." };
    juce::TextButton diagnosticsButton { "Measured diagnostics..." };
    std::unique_ptr<StandaloneDiagnosticsPanel> diagnosticsPanel;
    std::unique_ptr<StandaloneSectionsPanel> sectionsPanel;
    std::unique_ptr<juce::DocumentWindow> sectionsWindow;
    std::unique_ptr<juce::DocumentWindow> diagnosticsWindow;
    std::shared_ptr<JobControl> importControl;
    std::future<MediaValidationResult> importResult;
    std::future<ServiceResult<AnalysisSnapshot>> analysisResult;
    std::future<proposal::TrimProposal> proposalResult;
    JobTicket proposalTicket;
    MeasurementDomain analysisDomain = MeasurementDomain::source;
    JobTicket analysisTicket;
    JobTicket importTicket;
    juce::File importingFile;
    juce::Label jobLabel;
    juce::TextButton cancelJobButton { "Cancel job" };

    juce::Label titleLabel;
    juce::Label subtitleLabel;
    juce::Label workflowStateLabel;
    juce::Label targetProfileLabel;

    juce::ComboBox targetProfileBox;
    juce::ToggleButton applicabilityConfirmation { "Confirm pre-produced whole-programme linear-audio QC" };

    juce::TextButton importButton { "Import" };
    juce::TextButton analyzeButton { "Analyze" };
    juce::TextButton currentStateButton { "Current State" };
    juce::TextButton exportButton { "Export" };

    std::unique_ptr<ProgrammeImportPanel> programmePanel;
    std::unique_ptr<juce::DocumentWindow> programmeWindow;
    juce::TextButton chooseProgrammeButton { "Import channel files..." };
    void chooseProgramme();
    juce::TextButton chooseSourceButton { "Choose Source File..." };
    juce::TextButton clearSourceButton { "Clear Source" };
    juce::TextButton measureSourceButton { "Measure Source" };

    juce::TextButton initializeCurrentStateButton { "Initialize Current From Source" };
    juce::TextButton clearCurrentStateButton { "Clear Current State" };
    juce::TextButton remeasureCurrentStateButton { "Re-measure Current State" };
    juce::Slider currentGain;
    juce::TextButton draftTrimButton { "Draft global trim" };
    juce::TextButton applyTrimButton { "Apply trim" };
    juce::TextButton discardTrimButton { "Discard draft" };

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
    void beginMeasurement (MeasurementDomain);
    void beginTrimProposal();
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

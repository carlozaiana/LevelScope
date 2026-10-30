#pragma once
#include "SessionDocument.h"
#include "StandaloneProgrammeImport.h"
#include "StandaloneCompositionImport.h"
#include "StandaloneExport.h"
#include "StandaloneProposalPreferences.h"
#include "StandaloneJobRequirements.h"

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
    JobTicket proposalTicket, sectionProposalTicket;
    std::future<SectionProposalOutcome> sectionProposalResult;
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

    JobRequirementsPanel jobRequirementsPanel { session, [this] { refreshFromSession(); } };
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
    std::unique_ptr<CompositionImportPanel> compositionPanel;
    std::unique_ptr<juce::DocumentWindow> compositionWindow;
    juce::TextButton chooseCompositionButton { "Compose M&E + dialogue..." };
    void chooseComposition();
    juce::TextButton chooseSourceButton { "Choose Source File..." };
    juce::TextButton clearSourceButton { "Clear Source" };
    juce::TextButton measureSourceButton { "Measure Source" };

    juce::TextButton initializeCurrentStateButton { "Initialize Current From Source" };
    juce::TextButton clearCurrentStateButton { "Clear Current State" };
    juce::TextButton remeasureCurrentStateButton { "Re-measure Current State" };
    juce::Slider currentGain;
    juce::TextButton sectionProposalButton { "Generate section-led" }, applySectionProposalButton { "Apply proposal" }, discardSectionProposalButton { "Discard proposal" };
    std::unique_ptr<juce::DocumentWindow> preferencesWindow,proposalInspectionWindow;
    juce::TextButton preferencesButton { "Proposal limits..." },inspectProposalButton { "Inspect complete proposal..." };
    juce::TextButton draftTrimButton { "Draft global trim" };
    juce::TextButton applyTrimButton { "Apply trim" };
    juce::TextButton discardTrimButton { "Discard draft" };
    juce::TextButton renderPcmButton { "Render and verify..." }, publishPcmButton { "Publish export" },
        measurementReportButton { "Save measurement report..." }, artifactReportButton { "Save artifact report..." }, discardExportButton { "Discard rendered artifact" };
    std::unique_ptr<juce::FileChooser> exportChooser;
    std::future<ServiceResult<std::shared_ptr<exporting::Candidate>>> renderResult;
    std::future<ServiceResult<juce::File>> exportWriteResult;
    std::shared_ptr<exporting::Candidate> exportCandidate;
    std::optional<exporting::Request> exportRequest;
    JobTicket exportTicket;
    juce::File exportDestination;
    juce::String exportNotice = "No rendered artifact. Measurement-only reporting remains available without an export pass.";
    std::shared_ptr<const VerifiedSource> exportSourceCopy;
    bool publishingExport = false, savingArtifactReport = false;
    void configureExportControls();
    void chooseExportDestination (bool reportOnly, bool artifactOnly = false);
    void publishExport();
    void pollExport();
    void refreshExportControls();
    juce::String exportDetails() const;

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
    void beginSectionProposal();
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

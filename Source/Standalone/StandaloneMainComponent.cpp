#include "StandaloneMainComponent.h"

namespace levelscope::standalone
{
namespace
{
class DiagnosticsWindow final : public juce::DocumentWindow
{
public:
    explicit DiagnosticsWindow (StandaloneDiagnosticsPanel& p)
        : DocumentWindow ("LevelScope - Measured diagnostics", juce::Colour (0xff151b24), DocumentWindow::allButtons), panel (p)
    {
        setUsingNativeTitleBar (true); setResizable (true, false); setResizeLimits (900, 800, 2400, 1600);
        setContentNonOwned (&panel, true); centreWithSize (1040, 850); setVisible (true);
    }
    void closeButtonPressed() override { panel.suspend(); setVisible (false); }
private:
    StandaloneDiagnosticsPanel& panel;
};
class ReferenceWindow final : public juce::DocumentWindow
{
public:
    explicit ReferenceWindow (juce::Component& panel)
        : DocumentWindow ("LevelScope - Reference picture", juce::Colour (0xff181b20), DocumentWindow::allButtons)
    {
        setUsingNativeTitleBar (true); setResizable (true, false); setResizeLimits (620, 480, 2400, 1600);
        setContentNonOwned (&panel, true); centreWithSize (700, 600); setVisible (true);
    }
    void closeButtonPressed() override { setVisible (false); }
};
class SectionsWindow final : public juce::DocumentWindow
{
public:
    explicit SectionsWindow (StandaloneSectionsPanel& p)
        : DocumentWindow ("LevelScope - Programme sections", juce::Colour (0xff151b24), DocumentWindow::allButtons), panel (p)
    {
        setUsingNativeTitleBar (true); setResizable (true, false); setResizeLimits (1000, 820, 2400, 1800);
        setContentNonOwned (&panel, true); centreWithSize (1080, 850); setVisible (true);
    }
    void closeButtonPressed() override { panel.suspend(); setVisible (false); }
private:
    StandaloneSectionsPanel& panel;
};
juce::String metricText (const measurement::Metric& metric, const juce::String& unit)
{
    if (metric.value) return juce::String (*metric.value, 1) + " " + unit;
    switch (metric.status)
    {
        case measurement::Status::silence: return "-- (digital silence)";
        case measurement::Status::belowGate: return "-- (no eligible loudness blocks)";
        default: return "-- (insufficient duration)";
    }
}
juce::String measurementsText (const MeasurementSummary& summary)
{
    if (! summary.analysis) return "Integrated: --\nLRA: --\nTrue peak: --\n";
    const auto& r = summary.analysis->detail;
    juce::String text;
    text << "Integrated: " << metricText (r.integrated, "LUFS") << "\n"
         << "LRA: " << metricText (r.lra, "LU") << (r.lra.value && ! r.lraStable ? " (under 60 s; unstable)" : "") << "\n"
         << "True peak: " << metricText (r.truePeak, "dBTP") << "\n"
         << "Sample peak: " << metricText (r.samplePeak, "dBFS") << "\n"
         << "Max M: " << metricText (r.momentaryMax, "LUFS") << "\n"
         << "Max S: " << metricText (r.shortTermMax, "LUFS") << "\n";
    if (r.integratedRelativeGate) text << "I relative gate: " << juce::String (*r.integratedRelativeGate, 1) << " LUFS\n";
    if (r.lraRelativeGate) text << "LRA relative gate: " << juce::String (*r.lraRelativeGate, 1) << " LUFS\n";
    text << "Absolute gate: -70.0 LUFS\n" << "Engine: " << summary.analysis->measurementVersion << "\n";
    return text;
}
}
// [BEGIN LS-STANDALONE-WORKFLOW-SHELL]

StandaloneMainComponent::StandaloneMainComponent()
{
    configureLabels();
    configureTargetProfileBox();
    configureNavigationButtons();
    configureSourceControls();
    configureCurrentStateControls();

    configureReadOnlyBox (sourceStateBox);
    configureReadOnlyBox (targetStateBox);
    configureReadOnlyBox (currentStateBox);
    configureReadOnlyBox (pageDetailBox);

    addAndMakeVisible (sourceStateBox);
    addAndMakeVisible (targetStateBox);
    addAndMakeVisible (currentStateBox);
    addAndMakeVisible (pageDetailBox);
    addAndMakeVisible (audioPanel);
    addAndMakeVisible (videoButton);
    addAndMakeVisible (diagnosticsButton);
    diagnosticsButton.onClick = [this]
    {
        if (! diagnosticsPanel) diagnosticsPanel = std::make_unique<StandaloneDiagnosticsPanel> (session);
        diagnosticsPanel->onOpenSections = [this]
        {
            if (! sectionsPanel) sectionsPanel = std::make_unique<StandaloneSectionsPanel> (session);
            if (! sectionsWindow) sectionsWindow = std::make_unique<SectionsWindow> (*sectionsPanel);
            sectionsPanel->resume(); sectionsWindow->setVisible (true); sectionsWindow->toFront (true);
        };
        if (! diagnosticsWindow) diagnosticsWindow = std::make_unique<DiagnosticsWindow> (*diagnosticsPanel);
        diagnosticsPanel->resume(); diagnosticsWindow->setVisible (true); diagnosticsWindow->toFront (true);
    };
    videoButton.onClick = [this]
    {
        if (!videoWindow) videoWindow = std::make_unique<ReferenceWindow> (videoPanel);
        videoWindow->setVisible (true); videoWindow->toFront (true);
    };
    addAndMakeVisible (jobLabel);
    addAndMakeVisible (cancelJobButton);
    cancelJobButton.onClick = [this]
    {
        session.job.requestCancel();
        if (importControl) importControl->cancel.store (true);
        refreshFromSession();
    };

    refreshFromSession();

    setSize (1200, 900);
    startTimerHz (10);
}

// [BEGIN LS-STANDALONE-IMPORT-LIFECYCLE]
StandaloneMainComponent::~StandaloneMainComponent()
{
    stopTimer();
    videoWindow.reset();
    diagnosticsWindow.reset(); diagnosticsPanel.reset();
    sectionsWindow.reset(); sectionsPanel.reset();
    if (importControl) importControl->cancel.store (true);
    if (importResult.valid()) importResult.wait();
    if (analysisResult.valid()) analysisResult.wait();
    if (proposalResult.valid()) proposalResult.wait();
    audition.disconnect();
}

void StandaloneMainComponent::beginImport (const juce::File& file)
{
    if (session.job.isActive()) return;
    audition.stop();
    // Replacement keeps the previous validated source until success.
    if (! session.job.begin (JobStage::ingest, session.revision)) return;
    importTicket = session.job.ticket;
    importingFile = file;
    importControl = std::make_shared<JobControl>();
    importResult = std::async (std::launch::async, [file, control = importControl]
    {
        try { return validateMedia (file, *control); }
        catch (const std::exception& error)
        {
            return MediaValidationResult { {}, ErrorCode::decodeFailed, "Import failed: " + juce::String (error.what()) };
        }
    });
    refreshFromSession();
}

void StandaloneMainComponent::timerCallback()
{
    videoPanel.synchronize(); // observes new seek/play generation before audio preroll can finish
    audition.poll (videoPanel.prerollReady());
    audioPanel.refresh();
    if (proposalResult.valid() && proposalResult.wait_for (std::chrono::seconds (0)) == std::future_status::ready)
    {
        const auto result = proposalResult.get();
        // Job history may outlive its inputs. Keep verdicts and actionable advice
        // in trimNotice, which is invalidated with the proposal binding.
        if (session.job.finish (proposalTicket, session.revision, ErrorCode::none,
                                "Proposal calculation complete; see Current State"))
            if (! session.acceptTrimProposal (result))
            {
                session.job.state = JobState::failed;
                session.job.error = ErrorCode::staleRevision;
                session.job.message = "Proposal inputs changed; draft discarded";
            }
        importControl.reset(); refreshFromSession();
    }
    if (analysisResult.valid() && analysisResult.wait_for (std::chrono::seconds (0)) == std::future_status::ready)
    {
        auto result = analysisResult.get();
        if (session.job.finish (analysisTicket, session.revision, result.error, result.message))
        {
            if (! result.value || ! session.acceptMeasurement (std::move (*result.value)))
            {
                session.job.state = JobState::failed;
                session.job.error = ErrorCode::staleRevision;
                session.job.message = "Measurement identity no longer matches; result discarded";
                session.measurementUnsuccessful (analysisDomain, ErrorCode::staleRevision);
            }
        }
        else session.measurementUnsuccessful (analysisDomain, session.job.error);
        importControl.reset();
        refreshFromSession();
    }
    if (importResult.valid() && importResult.wait_for (std::chrono::seconds (0)) == std::future_status::ready)
    {
        const auto result = importResult.get();
        if (session.job.finish (importTicket, session.revision, result.error, result.message))
        {
            if (session.setValidatedSource (importingFile, result.media))
            {
                session.selectedPage = WorkflowPage::analyze;
                audioPanel.sourceChanged();
            }
            else
            {
                session.job.state = JobState::failed;
                session.job.error = ErrorCode::sourceChanged;
                session.job.message = "Source changed before import completed; previous source retained";
            }
        }
        importControl.reset();
        refreshFromSession();
    }
    auto text = session.job.message;
    if (session.job.isActive() && importControl)
        text << "  " << juce::String (juce::roundToInt (100 * importControl->progress.load())) << "%";
    jobLabel.setText (text, juce::dontSendNotification);
}

void StandaloneMainComponent::beginMeasurement (MeasurementDomain domain)
{
    if (session.job.isActive() || session.sourceContentInvalid) return;
    if (domain == MeasurementDomain::currentState && ! session.hasCurrentStateInitialized()) return;
    auto context = session.measurementContext (domain);
    if (OfflineAnalysisService::validateRequest (context, domain).isNotEmpty()) return;
    audition.stop();
    if (! session.job.begin (JobStage::analysis, session.revision)) return;
    analysisTicket = session.job.ticket;
    analysisDomain = domain;
    context.ticket = analysisTicket;
    session.job.message = "Verifying source snapshot and measuring";
    importControl = std::make_shared<JobControl>();
    analysisResult = std::async (std::launch::async, [context, domain, control = importControl]
    {
        return OfflineAnalysisService().analyze (context, domain, *control);
    });
    refreshFromSession();
}

void StandaloneMainComponent::beginTrimProposal()
{
    if (! session.getWorkflowReadiness().canDraftProposal()) return;
    const auto request = session.currentTrimRequest();
    session.discardTrim();
    if (! session.job.begin (JobStage::proposal, session.revision)) return;
    proposalTicket = session.job.ticket;
    session.job.message = "Drafting global trim; audio unchanged";
    importControl = std::make_shared<JobControl>();
    proposalResult = std::async (std::launch::async, [request]
    {
        try { return proposal::draftGlobalTrim (request); }
        catch (...) { return proposal::TrimProposal(); }
    });
    session.selectedPage = WorkflowPage::currentState;
    refreshFromSession();
}
// [END LS-STANDALONE-IMPORT-LIFECYCLE]

void StandaloneMainComponent::configureLabels()
{
    titleLabel.setText ("LevelScope Standalone", juce::dontSendNotification);
    titleLabel.setJustificationType (juce::Justification::centredLeft);
    titleLabel.setFont (juce::Font (28.0f, juce::Font::bold));
    titleLabel.setColour (juce::Label::textColourId, juce::Colours::white);
    addAndMakeVisible (titleLabel);

    subtitleLabel.setText ("Measured Source / Current State rules. Source audition. Complete delivery compliance is not evaluated.",
                           juce::dontSendNotification);
    subtitleLabel.setJustificationType (juce::Justification::centredLeft);
    subtitleLabel.setColour (juce::Label::textColourId, juce::Colour (0xffc6ccd7));
    addAndMakeVisible (subtitleLabel);

    workflowStateLabel.setJustificationType (juce::Justification::centredLeft);
    workflowStateLabel.setColour (juce::Label::textColourId, juce::Colour (0xffc6ccd7));
    addAndMakeVisible (workflowStateLabel);

    targetProfileLabel.setText ("Target profile:", juce::dontSendNotification);
    targetProfileLabel.setJustificationType (juce::Justification::centredRight);
    targetProfileLabel.setColour (juce::Label::textColourId, juce::Colour (0xffc6ccd7));
    addAndMakeVisible (targetProfileLabel);
}

void StandaloneMainComponent::configureTargetProfileBox()
{
    for (int i = 0; i < StandaloneSessionModel::numTargetProfiles; ++i)
        targetProfileBox.addItem (StandaloneSessionModel::getTargetProfile (i).displayName, i + 1);

    targetProfileBox.setSelectedId (session.selectedTargetProfileIndex + 1, juce::dontSendNotification);

    targetProfileBox.onChange = [this]
    {
        const auto selectedId = targetProfileBox.getSelectedId();

        if (selectedId > 0)
            session.setSelectedTargetProfileIndex (selectedId - 1);

        refreshFromSession();
    };

    addAndMakeVisible (targetProfileBox);
    applicabilityConfirmation.setTooltip ("Confirm programme context. Excludes live exceptions, deliberate lower targets, short-form supplements and distribution-specific requirements.");
    applicabilityConfirmation.onClick = [this]
    {
        session.setProgrammeApplicability (applicabilityConfirmation.getToggleState()
            ? proposal::Applicability::preProducedLinearProgramme : proposal::Applicability::unconfirmed);
        refreshFromSession();
    };
    addAndMakeVisible (applicabilityConfirmation);
}

void StandaloneMainComponent::configureNavigationButtons()
{
    importButton.onClick = [this] { setPage (WorkflowPage::importSource); };
    analyzeButton.onClick = [this] { setPage (WorkflowPage::analyze); };
    currentStateButton.onClick = [this] { setPage (WorkflowPage::currentState); };
    exportButton.onClick = [this] { setPage (WorkflowPage::exportResults); };

    addAndMakeVisible (importButton);
    addAndMakeVisible (analyzeButton);
    addAndMakeVisible (currentStateButton);
    addAndMakeVisible (exportButton);
}

void StandaloneMainComponent::configureSourceControls()
{
    chooseSourceButton.onClick = [this] { chooseSourceFile(); };
    clearSourceButton.onClick = [this] { clearSourceFile(); };

    measureSourceButton.onClick = [this] { beginMeasurement (MeasurementDomain::source); };

    addAndMakeVisible (chooseSourceButton);
    addAndMakeVisible (clearSourceButton);
    addAndMakeVisible (measureSourceButton);
}

void StandaloneMainComponent::configureCurrentStateControls()
{
    initializeCurrentStateButton.onClick = [this] { initializeCurrentStateFromSource(); };
    clearCurrentStateButton.onClick = [this] { clearCurrentState(); };

    remeasureCurrentStateButton.onClick = [this] { beginMeasurement (MeasurementDomain::currentState); };
    currentGain.setRange (-60, 24, 0); // Continuous double value; Apply must not quantize to tenths.
    currentGain.setNumDecimalPlacesToDisplay (12);
    currentGain.setSliderStyle (juce::Slider::LinearHorizontal);
    currentGain.setTextBoxStyle (juce::Slider::TextBoxRight, false, 100, 24);
    currentGain.setTextValueSuffix (" dB trim");
    currentGain.setTooltip ("Current State static gain, including LFE. Re-measure after editing. Source audition uses original audio.");
    currentGain.onValueChange = [this] { session.setCurrentGainDb (currentGain.getValue()); refreshFromSession(); };
    addAndMakeVisible (currentGain);

    addAndMakeVisible (initializeCurrentStateButton);
    addAndMakeVisible (clearCurrentStateButton);
    addAndMakeVisible (remeasureCurrentStateButton);
    draftTrimButton.onClick = [this] { beginTrimProposal(); };
    applyTrimButton.onClick = [this] { session.applyTrim(); refreshFromSession(); };
    discardTrimButton.onClick = [this] { session.discardTrim(); refreshFromSession(); };
    addAndMakeVisible (draftTrimButton);
    addAndMakeVisible (applyTrimButton);
    addAndMakeVisible (discardTrimButton);
}

void StandaloneMainComponent::configureReadOnlyBox (juce::TextEditor& box)
{
    box.setMultiLine (true);
    box.setReadOnly (true);
    box.setCaretVisible (false);
    box.setScrollbarsShown (true);

    box.setColour (juce::TextEditor::backgroundColourId, juce::Colour (0xff20242b));
    box.setColour (juce::TextEditor::textColourId, juce::Colour (0xffeef1f7));
    box.setColour (juce::TextEditor::outlineColourId, juce::Colour (0xff3a3f4a));
    box.setColour (juce::TextEditor::focusedOutlineColourId, juce::Colour (0xff6f7788));
}

void StandaloneMainComponent::chooseSourceFile()
{
    sourceFileChooser = std::make_unique<juce::FileChooser> (
        "Choose an audio source file",
        juce::File(),
        getSourceFileWildcard());

    juce::Component::SafePointer<StandaloneMainComponent> safeThis (this);

    sourceFileChooser->launchAsync (
        juce::FileBrowserComponent::openMode | juce::FileBrowserComponent::canSelectFiles,
        [safeThis] (const juce::FileChooser& chooser)
        {
            if (safeThis == nullptr)
                return;

            const auto chosenFile = chooser.getResult();

            if (! chosenFile.existsAsFile())
                return;

            safeThis->beginImport (chosenFile);
        });
}

void StandaloneMainComponent::clearSourceFile()
{
    audition.stop();
    session.job.requestCancel();
    if (importControl) importControl->cancel.store (true);
    session.clearSourceFile();
    audioPanel.sourceChanged();
    refreshFromSession();
}

void StandaloneMainComponent::initializeCurrentStateFromSource()
{
    session.initializeCurrentStateFromSource();
    refreshFromSession();
}

void StandaloneMainComponent::clearCurrentState()
{
    session.clearCurrentState();
    refreshFromSession();
}

void StandaloneMainComponent::setPage (WorkflowPage page)
{
    session.selectedPage = page;
    refreshFromSession();
}

void StandaloneMainComponent::refreshFromSession()
{
    workflowStateLabel.setText ("Workflow page: " + getPageName (session.selectedPage)
                                    + "   |   Global-trim proposals",
                                juce::dontSendNotification);

    targetProfileBox.setSelectedId (session.selectedTargetProfileIndex + 1, juce::dontSendNotification);
    applicabilityConfirmation.setToggleState (session.programmeApplicability == proposal::Applicability::preProducedLinearProgramme, juce::dontSendNotification);
    applicabilityConfirmation.setEnabled (session.sourceDocument.hasSource && ! session.job.isActive());

    clearSourceButton.setEnabled (session.sourceDocument.hasSource || session.job.isActive());
    chooseSourceButton.setEnabled (! session.job.isActive());
    targetProfileBox.setEnabled (! session.job.isActive());
    cancelJobButton.setEnabled (session.job.state == JobState::running);

    initializeCurrentStateButton.setEnabled (! session.job.isActive() && session.canInitializeCurrentStateFromSource());
    clearCurrentStateButton.setEnabled (! session.job.isActive() && session.hasCurrentStateInitialized());
    const auto measurementError = OfflineAnalysisService::validateRequest (session.measurementContext (MeasurementDomain::source), MeasurementDomain::source);
    measureSourceButton.setEnabled (! session.job.isActive() && ! session.sourceContentInvalid && measurementError.isEmpty());
    measureSourceButton.setTooltip (measurementError);
    remeasureCurrentStateButton.setEnabled (! session.job.isActive() && ! session.sourceContentInvalid && session.hasCurrentStateInitialized() && measurementError.isEmpty());
    currentGain.setEnabled (! session.job.isActive() && session.hasCurrentStateInitialized());
    currentGain.setValue (session.currentStateDocument.staticGainDb, juce::dontSendNotification);
    draftTrimButton.setEnabled (session.getWorkflowReadiness().canDraftProposal());
    applyTrimButton.setEnabled (session.canApplyTrim());
    discardTrimButton.setEnabled (! session.job.isActive() && session.trimDraft.has_value());

    sourceStateBox.setText (buildSourceStateText(), false);
    targetStateBox.setText (buildTargetStateText(), false);
    currentStateBox.setText (buildCurrentStateText(), false);
    pageDetailBox.setText (buildPageDetailText(), false);
    if (diagnosticsPanel) diagnosticsPanel->synchronize();
    if (sectionsPanel) sectionsPanel->synchronize();
    audioPanel.refresh();
}

juce::String StandaloneMainComponent::getSourceFileWildcard()
{
    return "*.wav;*.wave;*.aif;*.aiff;*.flac";
}

juce::String StandaloneMainComponent::getPageName (WorkflowPage page) const
{
    switch (page)
    {
        case WorkflowPage::importSource:  return "Import";
        case WorkflowPage::analyze:       return "Analyze";
        case WorkflowPage::currentState:  return "Current State";
        case WorkflowPage::exportResults: return "Export";
        default:                          return "Unknown";
    }
}

juce::String StandaloneMainComponent::buildSourceStateText() const
{
    juce::String text;

    text << "SOURCE\n\n";
    text << evaluationText (session.evaluation (MeasurementDomain::source)) << "\n";
    text << "File: " << session.sourceDocument.displayName << "\n";
    text << "Import status: " << session.sourceDocument.importStatus << "\n";
    text << "Measurement status: " << session.source.statusText << "\n\n";

    if (session.sourceDocument.hasSource)
    {
        text << "Extension: " << session.sourceDocument.extension << "\n";
        text << "Size: " << session.sourceDocument.fileSizeBytes << " bytes\n";
        text << "Path: " << session.sourceDocument.fullPath << "\n\n";
        const auto& media = session.sourceDocument.media;
        text << media.format << " | " << media.sampleRate << " Hz | " << media.channels << " ch | " << media.bitsPerSample << " bit\n";
        text << media.frames << " frames | " << juce::String (media.durationSeconds(), 3) << " seconds\n";
        text << "SHA-256: " << media.sha256 << "\n";
        text << "Channel roles: " << (session.sourceDocument.layout.confirmed ? "confirmed" : "confirmation required") << "\n\n";
    }
    else
    {
        text << "Extension: --\n";
        text << "Size: --\n";
        text << "Path: --\n\n";
    }

    text << measurementsText (session.source);

    return text;
}

juce::String StandaloneMainComponent::buildTargetStateText() const
{
    const auto profile = session.getSelectedTargetProfile();

    juce::String text;

    text << "TARGET\n\n";
    text << "Profile: " << profile.displayName << "\n";
    text << "Family: " << profile.familyName << "\n";
    text << "ID: " << profile.id << "\n\n";

    text << profile.description << "\n\n";

    text << StandaloneTargetProfileCatalog::buildLimitsText (profile) << "\n\n";

    text << "Workflow selection: ";
    text << (profile.isSelectableForWorkflow ? "target family selected" : "no target family selected");
    text << "\n";

    text << "Evaluation: selected measured rules only\n";
    text << "Proposal use: measured global trim only; explicit Apply and remeasurement required\n";

    return text;
}

juce::String StandaloneMainComponent::buildCurrentStateText() const
{
    juce::String text;

    text << "CURRENT STATE\n\n";
    text << evaluationText (session.evaluation (MeasurementDomain::currentState)) << "\n";
    text << "Initialized: " << (session.currentStateDocument.isInitialized ? "yes" : "no") << "\n";
    text << "Status: " << session.currentStateDocument.statusText << "\n";
    text << "Based on Source: " << session.currentStateDocument.basedOnSourceName << "\n";
    text << "Revision: " << session.currentStateDocument.revision << "\n";
    text << "Needs re-measurement: " << (session.currentStateNeedsRemeasurement() ? "yes" : "no") << "\n\n";

    text << "Measurement status: " << session.currentState.statusText << "\n\n";

    text << "Static gain (all channels): " << juce::String (session.currentStateDocument.staticGainDb, 12) << " dB\n\n";
    text << measurementsText (session.currentState);

    return text;
}

juce::String StandaloneMainComponent::buildPageDetailText() const
{
    juce::String text;

    text << getPageName (session.selectedPage).toUpperCase() << " VIEW\n\n";

    switch (session.selectedPage)
    {
        case WorkflowPage::importSource:
            text << "Validated source import.\n\n";
            text << "Implemented now:\n";
            text << "- Choose an audio source file.\n";
            text << "- Decode the entire file on a cancellable worker; record format, rate, frames, channels and SHA-256.\n";
            text << "- Confirm speaker roles separately from file order and physical outputs.\n";
            text << "- Clear the selected source.\n";
            text << "- Successful replacement clears Current State; failure/cancellation keeps the previous source.\n\n";
            text << "Accepted picker patterns:\n";
            text << getSourceFileWildcard() << "\n\n";

            if (session.sourceDocument.hasSource)
            {
                text << "Selected source:\n";
                text << session.sourceDocument.displayName << "\n\n";
                text << "Confirm roles, then use Measure Source.\n";
            }
            else
            {
                text << "No source selected yet.\n";
            }

            break;

        case WorkflowPage::analyze:
        {
            text << "Full-program offline measurement.\n\n";
            text << evaluationText (session.evaluation (MeasurementDomain::source)) << "\n";
            const auto error = OfflineAnalysisService::validateRequest (session.measurementContext (MeasurementDomain::source), MeasurementDomain::source);
            text << (error.isEmpty() ? "Ready to measure Source." : error) << "\n\n";
            text << "Supported: mono, stereo and 5.0/5.1 at 44.1/48/96/192 kHz.\n";
            text << "LFE is excluded from loudness and included in peak measurement.\n";
            text << "Measure Source reads a verified copy of the full original file.\n\n";
            text << measurementsText (session.source);

            break;
        }

        case WorkflowPage::currentState:
            text << session.trimNotice << "\n\n";
            text << evaluationText (session.evaluation (MeasurementDomain::currentState)) << "\n";
            text << "Initialize Current State, adjust static gain, then re-measure.\n\n";
            text << "Static gain applies to every file channel, including LFE.\n";
            text << "Editing gain invalidates the previous Current State measurement.\n";
            text << "Source audition plays the original audio with monitor gain only.\n\n";
            text << measurementsText (session.currentState);
            break;

        case WorkflowPage::exportResults:
        {
            const auto readiness = session.getWorkflowReadiness();

            text << "Workflow readiness scaffold.\n\n";
            text << readiness.buildChecklistText() << "\n\n";
            text << readiness.buildNextActionText() << "\n\n";

            text << "Proposal gate:\n";
            text << "- Inputs satisfied: " << (readiness.areProposalInputsSatisfied() ? "yes" : "no") << "\n";
            text << "- Global-trim proposal engine available: yes\n";
            text << "- Draft proposal: " << (readiness.canDraftProposal() ? "available" : "blocked; inputs required") << "\n\n";

            text << "Export gate:\n";
            text << "- Render/export available: no\n";
            text << "- Export: blocked\n\n";

            text << "Selected Source/Current State rules are evaluated. A rule pass does not verify an encoded export. Global-trim drafts require explicit Apply and measurement; public export remains unavailable.\n";
            break;
        }

        default:
            text << "Unknown workflow page.\n";
            break;
    }

    text << "\n\nSource audition uses original audio with monitor gain only. No processing or compliance claim.";

    return text;
}

void StandaloneMainComponent::paint (juce::Graphics& g)
{
    g.fillAll (juce::Colour (0xff15171c));

    g.setColour (juce::Colour (0xff2a2d35));
    g.drawRoundedRectangle (getLocalBounds().toFloat().reduced (12.0f), 10.0f, 1.0f);
}

void StandaloneMainComponent::resized()
{
    auto bounds = getLocalBounds().reduced (24);

    titleLabel.setBounds (bounds.removeFromTop (40));
    subtitleLabel.setBounds (bounds.removeFromTop (28));

    bounds.removeFromTop (10);

    auto targetRow = bounds.removeFromTop (32);
    targetProfileLabel.setBounds (targetRow.removeFromLeft (120));
    targetProfileBox.setBounds (targetRow.removeFromLeft (360));
    workflowStateLabel.setBounds (targetRow.reduced (16, 0));
    applicabilityConfirmation.setBounds (bounds.removeFromTop (28));

    bounds.removeFromTop (14);

    auto navRow = bounds.removeFromTop (34);
    importButton.setBounds (navRow.removeFromLeft (120));
    navRow.removeFromLeft (8);
    analyzeButton.setBounds (navRow.removeFromLeft (120));
    navRow.removeFromLeft (8);
    currentStateButton.setBounds (navRow.removeFromLeft (140));
    navRow.removeFromLeft (8);
    exportButton.setBounds (navRow.removeFromLeft (120));
    navRow.removeFromLeft (12);
    videoButton.setBounds (navRow.removeFromLeft (180));
    navRow.removeFromLeft (8);
    diagnosticsButton.setBounds (navRow.removeFromLeft (190));

    bounds.removeFromTop (10);

    auto sourceActionRow = bounds.removeFromTop (34);
    chooseSourceButton.setBounds (sourceActionRow.removeFromLeft (190));
    sourceActionRow.removeFromLeft (8);
    clearSourceButton.setBounds (sourceActionRow.removeFromLeft (120));
    sourceActionRow.removeFromLeft (8);
    measureSourceButton.setBounds (sourceActionRow.removeFromLeft (240));

    bounds.removeFromTop (8);

    auto currentStateActionRow = bounds.removeFromTop (34);
    initializeCurrentStateButton.setBounds (currentStateActionRow.removeFromLeft (230));
    currentStateActionRow.removeFromLeft (8);
    clearCurrentStateButton.setBounds (currentStateActionRow.removeFromLeft (150));
    currentStateActionRow.removeFromLeft (8);
    remeasureCurrentStateButton.setBounds (currentStateActionRow.removeFromLeft (300));
    currentStateActionRow.removeFromLeft (8);
    currentGain.setBounds (currentStateActionRow);

    bounds.removeFromTop (8);
    auto proposalRow = bounds.removeFromTop (30);
    draftTrimButton.setBounds (proposalRow.removeFromLeft (180));
    proposalRow.removeFromLeft (8);
    applyTrimButton.setBounds (proposalRow.removeFromLeft (120));
    proposalRow.removeFromLeft (8);
    discardTrimButton.setBounds (proposalRow.removeFromLeft (140));
    bounds.removeFromTop (8);

    auto jobRow = bounds.removeFromTop (28);
    cancelJobButton.setBounds (jobRow.removeFromRight (120));
    jobLabel.setBounds (jobRow);
    bounds.removeFromTop (8);
    auto summaryArea = bounds.removeFromTop (180);
    const auto gap = 12;
    const auto columnWidth = (summaryArea.getWidth() - (2 * gap)) / 3;

    sourceStateBox.setBounds (summaryArea.removeFromLeft (columnWidth));
    summaryArea.removeFromLeft (gap);
    targetStateBox.setBounds (summaryArea.removeFromLeft (columnWidth));
    summaryArea.removeFromLeft (gap);
    currentStateBox.setBounds (summaryArea);

    bounds.removeFromTop (14);

    pageDetailBox.setBounds (bounds.removeFromLeft (bounds.getWidth() * 38 / 100));
    bounds.removeFromLeft (12);
    audioPanel.setBounds (bounds);
}

// [END LS-STANDALONE-WORKFLOW-SHELL]
}

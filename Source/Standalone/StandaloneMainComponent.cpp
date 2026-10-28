#include "StandaloneMainComponent.h"

namespace levelscope::standalone
{
namespace
{
class ProgrammeWindow final : public juce::DocumentWindow
{
public:
    explicit ProgrammeWindow(ProgrammeImportPanel& p):DocumentWindow("Import channel files",juce::Colour(0xff181b20),DocumentWindow::closeButton),panel(p)
    {setUsingNativeTitleBar(true);setContentNonOwned(&p,true);centreWithSize(820,520);setVisible(true);}
    void closeButtonPressed() override {panel.cancel();}
private:ProgrammeImportPanel& panel;
};
class CompositionWindow final : public juce::DocumentWindow
{
public:
    explicit CompositionWindow(CompositionImportPanel& p):DocumentWindow("Compose M&E + dialogue",juce::Colour(0xff181b20),closeButton),panel(p)
    {setUsingNativeTitleBar(true);setContentNonOwned(&p,true);centreWithSize(p.getWidth(),p.getHeight());setVisible(true);}
    void closeButtonPressed() override{panel.cancel();}
private:CompositionImportPanel& panel;
};
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
        setUsingNativeTitleBar (true); setResizable (true, false); setResizeLimits (1000, 700, 2400, 1800);
        viewport.setViewedComponent(&panel,false);setContentNonOwned (&viewport, false);
        centreWithSize (1120, 900); setVisible (true);
    }
    void resized() override
    { DocumentWindow::resized();panel.setSize(juce::jmax(1000,viewport.getWidth()-16),1000); }
    void closeButtonPressed() override { panel.suspend(); setVisible (false); }
private:
    StandaloneSectionsPanel& panel;
    juce::Viewport viewport;
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
    if (summary.analysis->restoredFromSession) text << "Restored measurement; original provenance retained.\n";
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
    configureSessionControls();
    configureLabels();
    configureTargetProfileBox();
    configureNavigationButtons();
    configureSourceControls();
    configureCurrentStateControls();
    configureExportControls();

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
            sectionsPanel->onChanged=[this]{refreshFromSession();};
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

    setSize (1200, 950);
    startTimerHz (10);
}

// [BEGIN LS-STANDALONE-IMPORT-LIFECYCLE]
StandaloneMainComponent::~StandaloneMainComponent()
{
    stopTimer();
    if (documentControl) documentControl->cancel.store(true);
    if (documentResult.valid()) documentResult.wait();
    if (importControl) importControl->cancel.store (true);
    if (renderResult.valid()) renderResult.wait();
    if (exportWriteResult.valid()) exportWriteResult.wait();
    exportCandidate.reset();
    programmeWindow.reset(); programmePanel.reset();
    compositionWindow.reset(); compositionPanel.reset();
    videoWindow.reset();
    diagnosticsWindow.reset(); diagnosticsPanel.reset();
    sectionsWindow.reset(); sectionsPanel.reset();
    if (importControl) importControl->cancel.store (true);
    if (importResult.valid()) importResult.wait();
    if (analysisResult.valid()) analysisResult.wait();
    if (proposalResult.valid()) proposalResult.wait();
    if (sectionProposalResult.valid()) sectionProposalResult.wait();
    audition.disconnect();
}

void StandaloneMainComponent::beginImport (const juce::File& file)
{
    if (session.job.isActive() || session.persistenceLoading) return;
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
    pollExport();
    pollSessionOperation();
    audioPanel.synchronize(); // reject unexpectedly invalidated audition before any completion publishes
    videoPanel.synchronize(); // observes new seek/play generation before audio preroll can finish
    const auto audioStateBeforePoll = audition.getState();
    audition.poll (videoPanel.prerollReady());
    audioPanel.refresh();
    if (audition.getState() != audioStateBeforePoll) refreshFromSession();
    if(sectionProposalResult.valid()&&sectionProposalResult.wait_for(std::chrono::seconds(0))==std::future_status::ready)
    {
        const auto result=sectionProposalResult.get();
        if(session.job.finish(sectionProposalTicket,session.revision,result.error,result.message))
        {
            if(!result.value||!session.acceptSectionProposal(*result.value))
            {session.job.state=JobState::failed;session.job.message="Proposal inputs changed; private candidate discarded";}
        }
        else session.sectionProposalNotice=result.message;
        importControl.reset();refreshFromSession();
    }
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
    refreshSessionControls();
}

void StandaloneMainComponent::beginMeasurement (MeasurementDomain domain)
{
    if (session.job.isActive() || session.sourceContentInvalid || session.persistenceLoading) return;
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
    if (session.persistenceLoading || ! session.getWorkflowReadiness().canDraftProposal()) return;
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
void StandaloneMainComponent::beginSectionProposal()
{
    if(session.persistenceLoading||session.audioEditsBlocked()||session.job.isActive())return;
    session.synchronizeSections();const auto request=session.sectionProposalRequest();
    const auto first=proposal::draftGlobalTrim(request.scalar());
    session.discardSectionProposal();session.selectedPage=WorkflowPage::currentState;
    // No attribution, copy or measurement for existing global-only/no-change outcomes.
    if(first.status!=proposal::TrimStatus::noSolution)
    {
        JobControl control;const auto result=generateSectionProposal(request,control);
        if(result.value)session.acceptSectionProposal(*result.value);else session.sectionProposalNotice=result.message;
        refreshFromSession();return;
    }
    if(!session.job.begin(JobStage::proposal,session.revision))return;
    sectionProposalTicket=session.job.ticket;session.job.message="Assessing one section envelope; Current State unchanged";
    importControl=std::make_shared<JobControl>();
    sectionProposalResult=std::async(std::launch::async,[request,control=importControl]{return generateSectionProposal(request,*control);});
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

    subtitleLabel.setText ("Measured Source / Current State rules. Source / Current State audition. Complete delivery compliance is not evaluated.",
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
    chooseProgrammeButton.onClick = [this] { chooseProgramme(); };
    addAndMakeVisible(chooseProgrammeButton);
    chooseCompositionButton.onClick=[this]{chooseComposition();};addAndMakeVisible(chooseCompositionButton);
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
    currentGain.setTooltip ("Current State static gain, including LFE. Re-measure after editing. Audition selection is separate; Stop before editing processing gain.");
    currentGain.onValueChange = [this] { session.setCurrentGainDb (currentGain.getValue()); refreshFromSession(); };
    addAndMakeVisible (currentGain);

    addAndMakeVisible (initializeCurrentStateButton);
    addAndMakeVisible (clearCurrentStateButton);
    addAndMakeVisible (remeasureCurrentStateButton);
    sectionProposalButton.onClick=[this]{beginSectionProposal();};
    applySectionProposalButton.onClick=[this]{session.applySectionProposal();refreshFromSession();};
    discardSectionProposalButton.onClick=[this]{session.discardSectionProposal();refreshFromSession();};
    for(auto* button:{&sectionProposalButton,&applySectionProposalButton,&discardSectionProposalButton})addAndMakeVisible(*button);
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
    if (session.persistenceLoading) return;
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
    refreshSessionControls();
    refreshExportControls();
    workflowStateLabel.setText ("Workflow page: " + getPageName (session.selectedPage)
                                    + "   |   Global-trim proposals",
                                juce::dontSendNotification);

    targetProfileBox.setSelectedId (session.selectedTargetProfileIndex + 1, juce::dontSendNotification);
    applicabilityConfirmation.setToggleState (session.programmeApplicability == proposal::Applicability::preProducedLinearProgramme, juce::dontSendNotification);
    applicabilityConfirmation.setEnabled (session.sourceDocument.hasSource && ! session.job.isActive() && !session.persistenceLoading);

    clearSourceButton.setEnabled (!session.persistenceLoading && (session.sourceDocument.hasSource || session.job.isActive()));
    chooseProgrammeButton.setEnabled(canReplaceSession());
    chooseCompositionButton.setEnabled(canReplaceSession());
    chooseSourceButton.setEnabled (! session.job.isActive() && !session.persistenceLoading);
    targetProfileBox.setEnabled (! session.job.isActive() && !session.persistenceLoading);
    cancelJobButton.setEnabled (session.job.state == JobState::running);

    initializeCurrentStateButton.setEnabled (! session.audioEditsBlocked() && ! session.job.isActive() && session.canInitializeCurrentStateFromSource());
    clearCurrentStateButton.setEnabled (! session.audioEditsBlocked() && ! session.job.isActive() && session.hasCurrentStateInitialized());
    const auto measurementError = OfflineAnalysisService::validateRequest (session.measurementContext (MeasurementDomain::source), MeasurementDomain::source);
    measureSourceButton.setEnabled (!session.persistenceLoading && ! session.job.isActive() && ! session.sourceContentInvalid && measurementError.isEmpty());
    measureSourceButton.setTooltip (measurementError);
    remeasureCurrentStateButton.setEnabled (!session.persistenceLoading && ! session.job.isActive() && ! session.sourceContentInvalid && session.hasCurrentStateInitialized() && measurementError.isEmpty());
    currentGain.setEnabled (! session.audioEditsBlocked() && ! session.job.isActive() && session.hasCurrentStateInitialized());
    currentGain.setValue (session.currentStateDocument.staticGainDb, juce::dontSendNotification);
    draftTrimButton.setEnabled (!session.persistenceLoading && session.getWorkflowReadiness().canDraftProposal());
    applyTrimButton.setEnabled (session.canApplyTrim());
    sectionProposalButton.setEnabled(!session.persistenceLoading&&!session.audioEditsBlocked()&&!session.job.isActive()&&session.currentState.hasMeasurement);
    applySectionProposalButton.setEnabled(session.canApplySectionProposal());
    discardSectionProposalButton.setEnabled(bool(session.sectionProposal)&&!session.job.isActive());
    discardTrimButton.setEnabled (!session.persistenceLoading && ! session.job.isActive() && session.trimDraft.has_value());

    sourceStateBox.setText (buildSourceStateText(), false);
    targetStateBox.setText (buildTargetStateText(), false);
    currentStateBox.setText (buildCurrentStateText(), false);
    pageDetailBox.setText (buildPageDetailText(), false);
    if (diagnosticsPanel) diagnosticsPanel->synchronize();
    if (sectionsPanel) sectionsPanel->synchronize();
    audioPanel.refresh();
}

void StandaloneMainComponent::chooseProgramme()
{
    if(!canReplaceSession())return;
    confirmDiscard([this]
    {
        if(!canReplaceSession())return;
        audition.stop();session.persistenceLoading=true;
        juce::Component::SafePointer<StandaloneMainComponent> safe(this);
        programmePanel=std::make_unique<ProgrammeImportPanel>([safe](std::optional<ProgrammeAssembly> group)
        {
            juce::MessageManager::callAsync([safe,group=std::move(group)]
            {
                if(!safe)return;
                safe->programmeWindow.reset();safe->programmePanel.reset();safe->session.persistenceLoading=false;
                if(group)
                {
                    if(safe->session.setValidatedProgramme(*group)) {safe->audioPanel.sourceChanged();safe->documentNotice="Aligned programme imported; all routes muted";}
                    else safe->documentNotice="Programme changed before adoption; previous session retained";
                }
                else safe->documentNotice="Programme import cancelled; previous session retained";
                safe->refreshFromSession();
            });
        });
        programmeWindow=std::make_unique<ProgrammeWindow>(*programmePanel);refreshFromSession();
    });
}

void StandaloneMainComponent::chooseComposition()
{
    if(!canReplaceSession())return;
    confirmDiscard([this]
    {
        if(!canReplaceSession())return;
        audition.stop();audition.disconnect();session.persistenceLoading=true;
        juce::Component::SafePointer<StandaloneMainComponent> safe(this);
        compositionPanel=std::make_unique<CompositionImportPanel>(session.sourceDocument.composition,[safe](std::optional<ProgrammeComposition> recipe)
        {
            juce::MessageManager::callAsync([safe,recipe=std::move(recipe)]
            {
                if(!safe)return;safe->compositionWindow.reset();safe->compositionPanel.reset();safe->session.persistenceLoading=false;
                if(recipe)
                {
                    const auto before=safe->session.revision;
                    if(safe->session.setValidatedComposition(*recipe))
                    {
                        safe->audioPanel.sourceChanged();
                        safe->documentNotice=before==safe->session.revision
                            ?"Composition labels or locations updated; valid measurements retained; device closed, routes muted"
                            :"Confirmed composition committed; device closed, routes muted, fresh measurements required";
                    }
                    else safe->documentNotice="Composition changed before adoption; previous session retained";
                }
                else safe->documentNotice="Composition cancelled; previous session retained stopped";
                safe->refreshFromSession();
            });
        });
        compositionWindow=std::make_unique<CompositionWindow>(*compositionPanel);refreshFromSession();
    });
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

    text << (session.sourceDocument.composition?"SOURCE - ORIGINAL CONFIRMED COMPOSITION\n\n":"SOURCE\n\n");
    text << evaluationText (session.evaluation (MeasurementDomain::source)) << "\n";
    text << "Programme: " << session.sourceDocument.displayName << "\n";
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
        if(session.sourceDocument.assembly)
            for(const auto& member:session.sourceDocument.assembly->members)
                text << juce::AudioChannelSet::getAbbreviatedChannelTypeName(member.role) << ": " << member.file.getFileName() << "\n";
        if(session.sourceDocument.composition)
            for(const auto* component:{&session.sourceDocument.composition->me,&session.sourceDocument.composition->dialogue})
            {
                text << component->purpose << " | " << component->language << " | " << component->versionLabel << " | unity gain\n";
                for(const auto& member:component->files())text << "  " << member.file.getFullPathName() << "\n";
            }
        text << ((session.sourceDocument.composition||session.sourceDocument.assembly)?"Programme definition digest: ":"SHA-256: ") << media.sha256 << "\n";
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

    text << (session.sourceDocument.composition?"CURRENT STATE - COMPOSED PROGRAMME PLUS GLOBAL TRIM\n\n":"CURRENT STATE\n\n");
    text << evaluationText (session.evaluation (MeasurementDomain::currentState)) << "\n";
    text << "Initialized: " << (session.currentStateDocument.isInitialized ? "yes" : "no") << "\n";
    text << "Status: " << session.currentStateDocument.statusText << "\n";
    text << "Based on Source: " << session.currentStateDocument.basedOnSourceName << "\n";
    text << "Revision: " << session.currentStateDocument.revision << "\n";
    text << "Needs re-measurement: " << (session.currentStateNeedsRemeasurement() ? "yes" : "no") << "\n\n";

    text << "Measurement status: " << session.currentState.statusText << "\n\n";

    if(session.currentStateDocument.sectionGains) text << "Applied section gains: " << juce::String(int(session.currentStateDocument.sectionGains->intervals.size())) << " intervals; inspect Programme sections. All channels including LFE.\n";
    text << "Global trim (all channels): " << juce::String (session.currentStateDocument.staticGainDb, 12) << " dB\n\n";
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
            text << "- Assemble aligned mono channel files, or explicitly compose complementary M&E and dialogue.\n";
            text << "- Decode the entire file on a cancellable worker; record format, rate, frames, channels and SHA-256.\n";
            text << "- Confirm speaker roles separately from file order and physical outputs.\n";
            text << "- Clear the selected source.\n";
            text << "- New source import clears Current State. Editing existing composition audio retains its gain but invalidates measurements.\n";
            text << "- Failure/cancellation keeps the previous session stopped.\n\n";
            text << "Accepted picker patterns:\n";
            text << getSourceFileWildcard() << "\n\n";

            if (session.sourceDocument.hasSource)
            {
                text << "Selected source:\n";
                text << session.sourceDocument.displayName << "\n\n";
                text << session.sourceDocument.media.layoutProvenance << "\n";
                if (session.sourceDocument.layout.manuallyInterpreted)
                    text << "Explicit manual interpretation selected; original metadata remains recorded.\n";
                if (session.sourceDocument.media.measurementRestriction.isNotEmpty())
                    text << session.sourceDocument.media.measurementRestriction << "\n";
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
            text << "Supported: mono, stereo, 5.0/5.1 and explicit BS.2051 System I 7.1 at 44.1/48/96/192 kHz.\n";
            text << "7.1 requires distinct side/rear roles and confirmation; monitor muting never excludes a channel from measurement.\n";
            text << "LFE is excluded from loudness and included in peak measurement.\n";
            text << "Measure Source reads the full original programme from verified member copies.\n\n";
            text << measurementsText (session.source);

            break;
        }

        case WorkflowPage::currentState:
            text << session.trimNotice << "\n\n" << session.sectionProposalDisplayText() << "\n\n";
            text << evaluationText (session.evaluation (MeasurementDomain::currentState)) << "\n";
            text << "Initialize Current State, adjust static gain, then re-measure.\n\n";
            text << "Static gain applies to every programme channel, including LFE, after any composition.\n";
            text << "Editing gain invalidates the previous Current State measurement.\n";
            text << "Audition Source or the explicitly applied Current State; monitor gain is separate.\n\n";
            text << measurementsText (session.currentState);
            break;

        case WorkflowPage::exportResults:
        {
            text << exportDetails();
            break;
        }

        default:
            text << "Unknown workflow page.\n";
            break;
    }

    text << "\n\nSelected audition uses Source or Current State processing, followed by separate monitor gain. Playback is not measurement or compliance.";

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

    auto sessionRow = bounds.removeFromTop(30);
    for (auto* button : { &newSessionButton, &openSessionButton, &saveSessionButton, &saveAsSessionButton, &relinkAudioButton, &relinkVideoButton, &cancelSessionButton })
        button->setBounds(sessionRow.removeFromLeft(button == &saveSessionButton ? 70 : 135).reduced(2));
    sessionStatus.setBounds(bounds.removeFromTop(24));
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
    chooseProgrammeButton.setBounds(sourceActionRow.removeFromLeft(190));
    sourceActionRow.removeFromLeft(8);
    chooseCompositionButton.setBounds(sourceActionRow.removeFromLeft(210));
    sourceActionRow.removeFromLeft(8);
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
    proposalRow.removeFromLeft(8);sectionProposalButton.setBounds(proposalRow.removeFromLeft(190));
    proposalRow.removeFromLeft(8);applySectionProposalButton.setBounds(proposalRow.removeFromLeft(150));
    proposalRow.removeFromLeft(8);discardSectionProposalButton.setBounds(proposalRow.removeFromLeft(160));
    bounds.removeFromTop (8);

    auto exportRow = bounds.removeFromTop (30);
    for (auto* button : { &renderPcmButton, &publishPcmButton, &measurementReportButton, &artifactReportButton, &discardExportButton })
        button->setBounds (exportRow.removeFromLeft (button == &measurementReportButton ? 215 : 185).reduced (2,0));
    bounds.removeFromTop (8);

    auto jobRow = bounds.removeFromTop (28);
    cancelJobButton.setBounds (jobRow.removeFromRight (120));
    jobLabel.setBounds (jobRow);
    bounds.removeFromTop (8);
    // Keep the audition controls and a usable channel viewport at smaller window
    // heights. The workflow detail box also contains each state's full summary.
    const int availableSummary = juce::jlimit (0, 180, bounds.getHeight() - 14 - 390);
    auto summaryArea = bounds.removeFromTop (availableSummary >= 80 ? availableSummary : 0);
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

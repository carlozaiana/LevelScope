#include "TimelineTransport.h"
#include "StandaloneMainComponent.h"

namespace levelscope::standalone
{
namespace
{
class ProgrammeWindow final : public juce::DocumentWindow
{
public:
    explicit ProgrammeWindow(ProgrammeImportPanel& p):DocumentWindow("Import channel files",palette::panel(),DocumentWindow::closeButton),panel(p)
    {setUsingNativeTitleBar(true);setContentNonOwned(&p,true);centreWithSize(820,520);setVisible(true);}
    void closeButtonPressed() override {panel.cancel();}
private:ProgrammeImportPanel& panel;
};
class CompositionWindow final : public juce::DocumentWindow
{
public:
    explicit CompositionWindow(CompositionImportPanel& p):DocumentWindow("Compose M&E + dialogue",palette::panel(),closeButton),panel(p)
    {setUsingNativeTitleBar(true);setContentNonOwned(&p,true);centreWithSize(p.getWidth(),p.getHeight());setVisible(true);}
    void closeButtonPressed() override{panel.cancel();}
private:CompositionImportPanel& panel;
};
class DiagnosticsWindow final : public juce::DocumentWindow
{
public:
    explicit DiagnosticsWindow (StandaloneDiagnosticsPanel& p)
        : DocumentWindow ("LevelScope - Measured diagnostics", palette::panel(), DocumentWindow::allButtons), panel (p)
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
        : DocumentWindow ("LevelScope - Reference picture", palette::panel(), DocumentWindow::allButtons)
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
        : DocumentWindow ("LevelScope - Programme sections", palette::panel(), DocumentWindow::allButtons), panel (p)
    {
        setUsingNativeTitleBar (true); setResizable (true, false); setResizeLimits (1000, 700, 2400, 1800);
        viewport.setViewedComponent(&panel,false);setContentNonOwned (&viewport, false);
        centreWithSize (1120, 900); setVisible (true);
    }
    void resized() override
    { DocumentWindow::resized();panel.setSize(juce::jmax(1000,viewport.getWidth()-16),1120); }
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
    if (summary.analysis->identicalSourceMeasurementIdentity.isNotEmpty())
        text << "Initial Current measurement reused from the identical Source at 0 dB.\n";
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

StandaloneMainComponent::StandaloneMainComponent(std::shared_ptr<OutputDeviceCatalog> catalog, OutputDeviceSettings settings,std::unique_ptr<juce::AudioIODevice> output)
    : audition(std::move(output)),audioPanel(session, audition, [this] { refreshFromSession(); }, {}, std::move(catalog), std::move(settings))
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
        diagnosticsPanel->onOpenSections=[this]{showSections();};
        if (! diagnosticsWindow) diagnosticsWindow = std::make_unique<DiagnosticsWindow> (*diagnosticsPanel);
        diagnosticsPanel->resume(); diagnosticsWindow->setVisible (true); diagnosticsWindow->toFront (true);
    };
    videoButton.onClick = [this]
    {
        if (!videoWindow) videoWindow = std::make_unique<ReferenceWindow> (videoPanel);
        videoWindow->setVisible (true); videoWindow->toFront (true);
    };
    addAndMakeVisible (jobLabel);
    addAndMakeVisible (jobDetail);
    addAndMakeVisible (jobProgress);
    jobLabel.setComponentID("job-stage");jobDetail.setComponentID("job-detail");jobProgress.setComponentID("job-progress");
    jobLabel.setFont(juce::Font(17.0f, juce::Font::bold));
    jobProgress.setPercentageDisplay(false);
    jobProgress.setColour(juce::ProgressBar::foregroundColourId, juce::Colour(0xff42b8c4));
    addAndMakeVisible (cancelJobButton);
    cancelJobButton.onClick = [this]
    {
        if(programmePanel && programmePanel->preparationControl()) programmePanel->cancel();
        else if(compositionPanel && compositionPanel->preparationControl()) compositionPanel->cancel();
        else if(!session.job.isActive() && !proposalOperation.active()
            && (audition.preparing() || audition.getState()==AuditionState::buffering)) audition.stop();
        else
        {
            if(proposalOperation.active())proposalOperation.cancelRequested=true;
            session.job.requestCancel();
            if (importControl) importControl->cancel.store (true);
        }
        refreshFromSession();
    };

    configureTimelineShell();
    refreshFromSession();

    setSize (1240, 1000);
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
    audioWindow.reset();
    timelineDiagnostics.cancel();
    videoWindow.reset();
    diagnosticsWindow.reset(); diagnosticsPanel.reset();
    sectionsWindow.reset(); sectionsPanel.reset();
    if (importControl) importControl->cancel.store (true);
    if (importResult.valid()) importResult.wait();
    if (analysisResult.valid()) analysisResult.wait();
    if (waveformResult.valid()) waveformResult.wait();
    if (proposalResult.valid()) proposalResult.wait();
    if (sectionProposalResult.valid()) sectionProposalResult.wait();
    audition.disconnect();
}

void StandaloneMainComponent::beginImport (const juce::File& file)
{
    if (session.proposalOperationBusy || session.job.isActive() || session.persistenceLoading) return;
    audition.stop();
    // Replacement keeps the previous validated source until success.
    if (! session.job.begin (JobStage::ingest, session.revision)) return;
    importTicket = session.job.ticket;
    importingFile = file;
    importControl = std::make_shared<JobControl>();
    importControl->feedback = std::make_shared<JobFeedback>();
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
    proposalOperation.observe(session);
    advanceProposalOperation();
    if(sectionProposalResult.valid()&&sectionProposalResult.wait_for(std::chrono::seconds(0))==std::future_status::ready)
    {
        const auto result=sectionProposalResult.get();
        const bool workflow=sectionOperationGeneration!=0;
        const bool live=sectionProposalRuntime==session.sessionId&&(!workflow||
            (proposalOperation.generation==sectionOperationGeneration&&proposalOperation.phase==ProposalOperation::Phase::search&&proposalOperation.ownedBy(session)));
        bool finished=false;
        if(live)finished=session.job.finish(sectionProposalTicket,session.revision,result.error,result.message);
        else if(sectionProposalRuntime==session.sessionId&&session.job.ticket.id==sectionProposalTicket.id&&session.job.stage==JobStage::proposal)
            session.job.finish(sectionProposalTicket,session.revision,ErrorCode::staleRevision,"Proposal operation changed; old completion discarded");
        if(importControl==sectionControl)importControl.reset();sectionControl.reset();
        if(live&&workflow)
        {
            session.discardSectionProposal();
            if(finished)proposalOperation.candidateComplete(session,result);
            else proposalOperation.finish(session,proposalOperation.cancelRequested?ProposalOperation::Phase::cancelled:ProposalOperation::Phase::failed,
                proposalOperation.cancelRequested?"Propose cancelled - processing unchanged":"Proposal did not complete",result.message+"\nCurrent processing was not changed.");
        }
        else if(live)
        {
            if(finished&&(!result.value||!session.acceptSectionProposal(*result.value)))
            {session.job.state=JobState::failed;session.job.error=ErrorCode::staleRevision;session.sectionProposalError=ErrorCode::staleRevision;session.sectionProposalNotice=session.job.message="Proposal inputs changed; private candidate discarded";}
            else if(!finished){session.sectionProposalError=session.job.error;session.sectionProposalNotice=result.message;}
        }
        refreshFromSession();
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
        bool accepted=false;
        const bool workflow=analysisOperationGeneration!=0;
        const bool expectedPhase=proposalOperation.phase==ProposalOperation::Phase::source||proposalOperation.phase==ProposalOperation::Phase::current||proposalOperation.phase==ProposalOperation::Phase::verify;
        const bool live=analysisRuntime==session.sessionId&&(!workflow||
            (proposalOperation.generation==analysisOperationGeneration&&expectedPhase&&proposalOperation.ownedBy(session)));
        if (live&&session.job.finish (analysisTicket, session.revision, result.error, result.message))
        {
            if (!result.value || !session.acceptMeasurement(std::move(*result.value)))
            {
                session.job.state=JobState::failed;session.job.error=ErrorCode::staleRevision;
                session.job.message="Measurement identity no longer matches; result discarded";
                session.measurementUnsuccessful(analysisDomain,ErrorCode::staleRevision);
            }
            else accepted=true;
        }
        else if(live)session.measurementUnsuccessful(analysisDomain,session.job.error);
        else if(analysisRuntime==session.sessionId&&session.job.ticket.id==analysisTicket.id&&session.job.stage==JobStage::analysis)
            session.job.finish(analysisTicket,session.revision,ErrorCode::staleRevision,"Measurement operation changed; old completion discarded");
        if(accepted && waveformCache.accept(analysisWaveformBinding,WaveformBinding::from(session.measurementContext(analysisDomain),analysisDomain)))
        {
            if(analysisDomain==MeasurementDomain::source && session.source.analysis && session.currentState.analysis)
                waveformCache.shareInitialCurrent(*session.source.analysis,*session.currentState.analysis);
        }
        else waveformCache.reject(analysisWaveformBinding);
        if(loudnessProgress){loudnessProgress->abandon();loudnessProgress.reset();}
        if(importControl==analysisControl)importControl.reset();analysisControl.reset();
        if(live&&workflow)proposalOperation.measurementComplete(session,accepted,session.job.message);
        else if(accepted&&analysisDomain==MeasurementDomain::currentState&&proposalOperation.applied
            &&!proposalOperation.active()&&proposalOperation.phase!=ProposalOperation::Phase::stale&&proposalOperation.ownedBy(session))
        {
            // An explicit retry after cancelled/failed verification is still an
            // independent Current measurement of this exact applied choice.
            proposalOperation.phase=ProposalOperation::Phase::verify;proposalOperation.cancelRequested=false;
            proposalOperation.measurementComplete(session,true,session.job.message);
        }
        refreshFromSession();
    }
    if(waveformResult.valid() && waveformResult.wait_for(std::chrono::seconds(0))==std::future_status::ready)
    {
        auto result=waveformResult.get();
        const auto expected=WaveformBinding::from(session.measurementContext(waveformBinding.domain),waveformBinding.domain);
        if(session.job.finish(waveformBinding.ticket,session.revision,result.error,result.message)
            && result.value && waveformCache.accept(waveformBinding,expected))
        {
            if(waveformBinding.domain==MeasurementDomain::source && session.source.analysis && session.currentState.analysis)
                waveformCache.shareInitialCurrent(*session.source.analysis,*session.currentState.analysis);
        }
        else waveformCache.reject(waveformBinding);
        importControl.reset();refreshFromSession();
    }
    if (importResult.valid() && importResult.wait_for (std::chrono::seconds (0)) == std::future_status::ready)
    {
        const auto result = importResult.get();
        if (session.job.finish (importTicket, session.revision, result.error, result.message))
        {
            if (session.setValidatedSource (importingFile, result.media))
            {
                session.selectedPage = WorkflowPage::analyze;
                audioPanel.sourceChanged(true);
                showAudioSetup();
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
    refreshJobFeedback();
    refreshSessionControls();
    refreshTimeline();
    if(!session.proposalOperationBusy && !audition.active() && automaticSource.eligible(session))beginMeasurement(MeasurementDomain::source);
    preparePlaybackWhenIdle();
}

void StandaloneMainComponent::preparePlaybackWhenIdle()
{
    // One attempt for this runtime/request/position. Cancellation or a failure
    // does not start an endless retry loop; explicit Play remains the retry.
    if(session.job.isActive() || session.persistenceLoading || session.audioEditsBlocked()
        || audition.active() || audition.preparing() || audition.playbackReady()
        || !session.source.hasMeasurement || session.sourceContentInvalid
        || programmePanel || compositionPanel)return;
    const auto request=session.auditionRequest();
    if(request.validationError().isNotEmpty() || !request.context.sourceFilesMatch()
        || audition.positionFrames()>=request.context.media.frames)return;
    const auto source=waveformCache.view(WaveformBinding::from(session.measurementContext(MeasurementDomain::source),MeasurementDomain::source));
    if(source.state==WaveformProgress::State::collecting)return;
    const auto key=request.identity()+":"+juce::String(audition.positionFrames());
    if(key==automaticPlaybackAttempt && sameVerifiedAssetOwner(request.context.verifiedAssets,automaticPlaybackAssets))return;
    automaticPlaybackAttempt=key;automaticPlaybackAssets=request.context.verifiedAssets;
    audition.prepare(request);
    refreshJobFeedback();refreshTimeline();
}

void StandaloneMainComponent::beginMeasurement (MeasurementDomain domain,bool workflow)
{
    if(session.pendingBoundaryEdits()){jobLabel.setText(session.pendingBoundaryMessage(),juce::dontSendNotification);return;}
    if ((!workflow&&session.proposalOperationBusy) || session.job.isActive() || session.sourceContentInvalid || session.persistenceLoading) return;
    if (domain == MeasurementDomain::currentState && ! session.hasCurrentStateInitialized()) return;
    auto context = session.measurementContext (domain);
    if (OfflineAnalysisService::validateRequest (context, domain).isNotEmpty()) return;
    if(domain==MeasurementDomain::source)
        automaticSource.attemptedSource(session);
    audition.stop();automaticPlaybackAttempt.clear();
    // Keep command focus on a stable surface before the initiating measurement
    // button is disabled. JUCE must not transfer it to a neighbouring text field.
    auto* focusView=&timeline;
    for(auto* c=juce::Component::getCurrentlyFocusedComponent();c;c=c->getParentComponent())
        if(auto* view=dynamic_cast<MeasurementTimeline*>(c)){focusView=view;break;}
    focusView->focusForTransport();
    if (! session.job.begin (JobStage::analysis, session.revision)) return;
    analysisTicket = session.job.ticket;analysisRuntime=session.sessionId;analysisOperationGeneration=workflow?proposalOperation.generation:0;
    analysisDomain = domain;
    context.ticket = analysisTicket;
    session.job.message = "Verifying source snapshot and measuring";
    importControl = std::make_shared<JobControl>();
    importControl->feedback = std::make_shared<JobFeedback>();
    analysisControl=importControl;
    timeline.setWaveform({},false,{});
    analysisWaveformBinding=WaveformBinding::from(context,domain);
    auto preview=waveformCache.begin(analysisWaveformBinding);
    loudnessProgress=std::make_shared<LoudnessProgress>(analysisWaveformBinding);
    analysisResult = std::async (std::launch::async, [context, domain, control = importControl,preview,loudness=loudnessProgress]
    {
        return OfflineAnalysisService(nullptr,preview.get(),loudness.get()).analyze (context, domain, *control);
    });
    refreshFromSession();
}

void StandaloneMainComponent::beginWaveform()
{
    const auto domain=inspectDomain.getSelectedId()==2?MeasurementDomain::currentState:MeasurementDomain::source;
    const auto& measured=domain==MeasurementDomain::source?session.source:session.currentState;
    if(!measured.hasMeasurement || session.job.isActive() || session.audioEditsBlocked() || session.sourceContentInvalid)return;
    auto context=session.measurementContext(domain);
    if(OfflineAnalysisService::validateRequest(context,domain).isNotEmpty() || !session.job.begin(JobStage::waveform,session.revision))return;
    context.ticket=session.job.ticket;waveformBinding=WaveformBinding::from(context,domain);
    timeline.setWaveform({},false,{});auto preview=waveformCache.begin(waveformBinding);
    importControl=std::make_shared<JobControl>();session.job.message="Verifying source and building waveform; measurements retained";
    waveformResult=std::async(std::launch::async,[context,domain,preview,control=importControl]
    {return OfflineAnalysisService(nullptr,preview.get()).buildWaveform(context,domain,*control);});
    refreshFromSession();
}

void StandaloneMainComponent::beginTrimProposal()
{
    if (session.audioEditsBlocked() || ! session.getWorkflowReadiness().canDraftProposal()) return;
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
void StandaloneMainComponent::beginSectionProposal(bool workflow)
{
    if(session.persistenceLoading||(!workflow&&session.audioEditsBlocked())||session.job.isActive()||audition.active())return;
    session.synchronizeSections();const auto request=workflow?proposalOperation.baseline:session.sectionProposalRequest();
    const auto first=proposal::draftGlobalTrim(request.scalar());
    session.discardSectionProposal();session.selectedPage=WorkflowPage::currentState;
    // No attribution, copy or measurement for existing global-only/no-change outcomes.
    if(!workflow&&!request.regionSettings&&(first.status==proposal::TrimStatus::noChange||(!request.currentContext.profile.limits.hasMaxLoudnessRange&&first.status!=proposal::TrimStatus::noSolution)))
    {
        JobControl control;const auto result=generateSectionProposal(request,control);
        if(result.value)session.acceptSectionProposal(*result.value);else {session.sectionProposalError=result.error;session.sectionProposalNotice=result.message;}
        refreshFromSession();return;
    }
    if(!session.job.begin(JobStage::proposal,session.revision))return;
    sectionProposalTicket=session.job.ticket;sectionProposalRuntime=session.sessionId;sectionOperationGeneration=workflow?proposalOperation.generation:0;session.job.message=request.regionSettings?"Assessing experimental regions; forecasts are not measurements; Current unchanged":request.currentContext.profile.limits.hasMaxLoudnessRange?"Assessing bounded LRA candidates; Current unchanged":"Assessing one section envelope; Current State unchanged";
    importControl=std::make_shared<JobControl>();
    sectionControl=importControl;
    sectionProposalResult=std::async(std::launch::async,[request,control=importControl,workflow]
    {
        auto result=generateSectionProposal(request,*control);
        if(workflow&&result.value&&(*result.value)->status()==SectionProposalStatus::globalOnly)
            result=measureGlobalProposal(request,*result.value,*control);
        return result;
    });
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
    addAndMakeVisible(jobRequirementsPanel);
    qcScope.setText("PASS / FAIL applies to this file and the selected measured limits.",juce::dontSendNotification);
    qcScope.setTooltip("The complete imported file is evaluated automatically. For excerpts, results describe the excerpt; they do not certify the parent programme or unselected delivery requirements.");
    addAndMakeVisible(qcScope);
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
    currentGain.onValueChange = [this] { session.synchronizeSections();auto next=session.manualState();next.globalDb=currentGain.getValue();session.commitManualEdit(std::move(next),"Change global gain");refreshFromSession(); };
    addAndMakeVisible (currentGain);

    addAndMakeVisible (initializeCurrentStateButton);
    addAndMakeVisible (clearCurrentStateButton);
    addAndMakeVisible (remeasureCurrentStateButton);
    preferencesButton.onClick=[this]{ProposalPreferencesPanel::showWindow(preferencesWindow, session, this, [this]{refreshFromSession();});};
    inspectProposalButton.onClick=[this]{reviewProposal(false);};
    addAndMakeVisible(preferencesButton);addAndMakeVisible(inspectProposalButton);
    sectionProposalButton.onClick=[this]{proposeWithActiveMethod();};
    applySectionProposalButton.onClick=[this]{if(gainLane.diagnosticRegionInspection())return;session.applySectionProposal();refreshFromSession();};
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

    box.setColour (juce::TextEditor::backgroundColourId, palette::card());
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
    proposalOperation.observe(session);
    refreshSessionControls();
    refreshExportControls();
    proposalPanel.refresh(importControl?importControl->progress.load():-1);
    workflowStateLabel.setText ("Workflow page: " + getPageName (session.selectedPage)
                                    + "   |   Global-trim proposals",
                                juce::dontSendNotification);

    targetProfileBox.setSelectedId (session.selectedTargetProfileIndex + 1, juce::dontSendNotification);

    clearSourceButton.setEnabled (!session.proposalOperationBusy && !session.pendingBoundaryEdits() && !session.persistenceLoading && (session.sourceDocument.hasSource || session.job.isActive()));
    chooseProgrammeButton.setEnabled(canReplaceSession());
    chooseCompositionButton.setEnabled(canReplaceSession());
    chooseSourceButton.setEnabled (!session.proposalOperationBusy && !session.pendingBoundaryEdits() && ! session.job.isActive() && !session.persistenceLoading);
    targetProfileBox.setEnabled (!session.audioEditsBlocked() && !session.job.isActive());
    jobRequirementsPanel.refresh();
    refreshJobFeedback();

    initializeCurrentStateButton.setEnabled (! session.audioEditsBlocked() && ! session.job.isActive() && session.canInitializeCurrentStateFromSource());
    clearCurrentStateButton.setEnabled (! session.audioEditsBlocked() && ! session.job.isActive() && session.hasCurrentStateInitialized());
    const auto measurementError = OfflineAnalysisService::validateRequest (session.measurementContext (MeasurementDomain::source), MeasurementDomain::source);
    measureSourceButton.setEnabled (!session.proposalOperationBusy && !session.pendingBoundaryEdits() && !session.persistenceLoading && ! session.job.isActive() && ! session.sourceContentInvalid && measurementError.isEmpty());
    measureSourceButton.setTooltip (measurementError);
    remeasureCurrentStateButton.setEnabled (!session.proposalOperationBusy && !session.pendingBoundaryEdits() && !session.persistenceLoading && ! session.job.isActive() && ! session.sourceContentInvalid && session.hasCurrentStateInitialized() && measurementError.isEmpty());
    currentGain.setEnabled (! session.audioEditsBlocked() && ! session.job.isActive() && session.hasCurrentStateInitialized() && !(session.sectionTrimDraft&&session.sectionTrimDraft->rollbackAnnotations));
    currentGain.setValue (session.currentStateDocument.staticGainDb, juce::dontSendNotification);
    draftTrimButton.setEnabled (!session.audioEditsBlocked() && session.getWorkflowReadiness().canDraftProposal());
    applyTrimButton.setEnabled (session.canApplyTrim());
    preferencesButton.setEnabled(!session.persistenceLoading&&!session.audioEditsBlocked()&&!session.job.isActive());
    inspectProposalButton.setEnabled(true);
    sectionProposalButton.setEnabled(!session.persistenceLoading&&!session.audioEditsBlocked()&&!session.job.isActive()&&session.currentState.hasMeasurement);
    applySectionProposalButton.setEnabled(session.canApplySectionProposal()&&!gainLane.diagnosticRegionInspection());
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
            if(!safe || !safe->programmePanel)return;
            const auto confirmedOrder=safe->programmePanel->confirmedFileOrder();
            juce::MessageManager::callAsync([safe,group=std::move(group),confirmedOrder]
            {
                if(!safe)return;
                safe->programmeWindow.reset();safe->programmePanel.reset();safe->session.persistenceLoading=false;
                if(group)
                {
                    if(safe->session.setValidatedProgramme(*group))
                    {
                        safe->audioPanel.sourceChanged();
                        const auto saved=RememberedChannelRoles().remember(confirmedOrder);
                        safe->documentNotice=saved.wasOk()?"Aligned programme imported":saved.getErrorMessage();
                        safe->beginMeasurement(MeasurementDomain::source);
                    }
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
                            ?"Composition labels or locations updated; valid measurements retained; stopped, local monitoring restored when matched"
                            :"Confirmed composition committed; stopped, local monitoring restored when matched; fresh measurements required";
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
        if(media.layoutProvenance.contains("RF64 input:"))text << media.layoutProvenance << "\n";
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
    text << session.jobRequirements.text() << "\n";
    text << "Base profile: " << profile.displayName << "\n";
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
            if(!session.sectionProposal)text << "Global-trim draft: " << session.trimNotice << "\n\n";
            text << session.proposalPolicy.text() << "\n" << session.sectionProposalDisplayText().substring(0,1800) << "\nOpen Inspect complete proposal for all rows and final outcome; this is a summary only.\n" << "\n\n";
            text << evaluationText (session.evaluation (MeasurementDomain::currentState)) << "\n";
            text << "The first Source measurement also measures the untouched Current State. After a gain change, re-measure Current State.\n\n";
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
    g.fillAll (palette::panel());

    g.setColour (juce::Colour (0xff2a2d35));
    g.drawRoundedRectangle (getLocalBounds().toFloat().reduced (12.0f), 10.0f, 1.0f);
}


void StandaloneMainComponent::configureTimelineShell()
{
    addChildComponent(monitorWarning);addChildComponent(monitorSetup);
    monitorWarning.setName("Monitor overload warning");monitorWarning.setFont(juce::Font(15.f,juce::Font::bold));
    monitorWarning.setColour(juce::Label::backgroundColourId,juce::Colour(0xff4a231d));
    monitorWarning.setColour(juce::Label::textColourId,juce::Colour(0xffffdfbf));
    monitorSetup.onClick=[this]{showAudioSetup();};
    audioPanel.onAnalyse=[this]
    {
        beginMeasurement(MeasurementDomain::source);
        if(audioWindow && session.job.isActive())audioWindow->setVisible(false);
    };
    timeline.attachGainLane(gainLane);
    timeline.bindSectionEditing(session,[this](uint64_t id)
    {
        refreshFromSession();
        refreshTimeline();
        for(const auto& row:session.sectionDocument.map().sections)
            if(row.id==id){gainLane.selectFrame(row.begin,true);break;}
    });
    timeline.onTogglePlayback=[this]{if(audition.active()){audition.stop();refreshTimeline();}else playTimeline();};
    timeline.onToggleInspectionAudition=[this]{toggleInspectionAudition();};
    gainLane.onChanged=[this]{refreshFromSession();refreshTimeline();};
    proposalPanel.setOperation(&proposalOperation);
    proposalPanel.onRecover=[this]{proposalOperation.recover(session);timeline.focusForTransport();refreshFromSession();refreshTimeline();};
    proposalPanel.onAlternative=[this](size_t index)
    {if(index!=proposalOperation.selectedChoice&&proposalOperation.select(session,index,false)){timeline.focusForTransport();refreshFromSession();}};
    proposalPanel.onBoundaries=[this]
    {
        if(!proposalOperation.ownedBy(session))return;
        inspectDomain.setSelectedId(2,juce::dontSendNotification);refreshTimeline();
        int64_t frame=0;
        if(const auto& plan=session.currentStateDocument.sectionGains;plan&&plan->intervals.size()>1)frame=plan->intervals[1].begin;
        gainLane.inspectApplied(frame);timeline.highlight(double(std::max<int64_t>(0,frame-48000)),double(std::min<int64_t>(session.sourceDocument.media.frames,frame+48000)));
        pageDetailBox.setText(proposalOperation.boundaryAdvice(),false);
    };
    proposalPanel.onPropose=[this]{proposeWithActiveMethod();refreshTimeline();};
    proposalPanel.onRegions=[this]{showRegionPlanning();};
    lraPassagesButton.onClick=[this]{showLraPassages();};
    proposalPanel.onReview=[this]{reviewProposal(false);};proposalPanel.onHistory=[this]{reviewProposal(true);};
    proposalPanel.onPreferences=[this]{preferencesButton.onClick();};
    proposalPanel.onChanged=[this]{refreshFromSession();refreshTimeline();};
    for(auto* c:std::initializer_list<juce::Component*>{&titleLabel,&subtitleLabel,&workflowStateLabel,&targetProfileLabel,&targetProfileBox,
        &newSessionButton,&openSessionButton,&saveAsSessionButton,&relinkAudioButton,&relinkVideoButton,&cancelSessionButton,
        &importButton,&analyzeButton,&currentStateButton,&exportButton,&chooseSourceButton,&chooseProgrammeButton,&chooseCompositionButton,
        &clearSourceButton,&initializeCurrentStateButton,&clearCurrentStateButton,&sourceStateBox,&targetStateBox,&currentStateBox,
        &renderPcmButton,&publishPcmButton,&measurementReportButton,&artifactReportButton,&discardExportButton,&audioPanel,
        &preferencesButton,&inspectProposalButton,&sectionProposalButton,&applySectionProposalButton,&discardSectionProposalButton,
        &draftTrimButton,&applyTrimButton,&discardTrimButton})c->setVisible(false);
    for(auto* c:std::initializer_list<juce::Component*>{&openMixButton,&toolsButton,&audioSetupButton,&timelinePlay,&timelineStop,&sectionsButton,
        &inspectDomain,&auditionDomainBox,&positionLabel,&timeline,&comparisonDivider,&waveformDivider,&buildWaveformButton,&proposalPanel})addAndMakeVisible(c);
    openMixButton.onClick=[this]{chooseMix();};toolsButton.onClick=[this]{showTools();};audioSetupButton.onClick=[this]{showAudioSetup();};
    timelinePlay.onClick=[this]{playTimeline();};timelineStop.onClick=[this]{timeline.focusForTransport();audition.stop();refreshTimeline();};
    inspectDomain.addItem("Inspect Original",1);inspectDomain.addItem("Inspect Current",2);inspectDomain.setSelectedId(1,juce::dontSendNotification);
    inspectDomain.onChange=[this]{timeline.clearHighlight();refreshTimeline();pageDetailBox.setText(inspectDomain.getText()+": select a passage or inspect the entered time.",false);};
    auditionDomainBox.addItem("Audition Original",1);auditionDomainBox.addItem("Audition Current",2);auditionDomainBox.setSelectedId(1,juce::dontSendNotification);
    auditionDomainBox.onChange=[this]{session.selectAuditionDomain(auditionDomainBox.getSelectedId()==2?MeasurementDomain::currentState:MeasurementDomain::source);audioPanel.synchronize();refreshTimeline();};
    addAndMakeVisible(jobRequirementsPanel.informationLine());
    comparisonDivider.move=[this](int delta){timelineLayout.moveComparison(delta);resized();};
    waveformDivider.move=[this](int delta){timelineLayout.moveTimeline(delta,timelineAvailableHeight,session.sourceDocument.layout.count,gainLane.preferredHeight());resized();};
    timeline.onRequestSeek=[this](int64_t frame,const juce::String& identity)
    {const auto error=TimelineTransport::seek(frame,identity,session,audition);refreshTimeline();return error;};
    sectionsButton.onClick=[this]{showSections();};
    buildWaveformButton.onClick=[this]{beginWaveform();};
    timeline.onInspect=[this](const juce::String& text){pageDetailBox.setText(inspectDomain.getText()+" | "+text,false);};
    setWantsKeyboardFocus(true);
    addAndMakeVisible(sourceFilename);sourceFilename.setName("Programme filename");sourceFilename.setFont(14.f);
    sourceFilename.setMinimumHorizontalScale(1.f);sessionStatus.setFont(11.f);
    addAndMakeVisible(inspectorViewport);inspectorViewport.setViewedComponent(&inspectorContent,false);inspectorViewport.setScrollBarsShown(true,false);
    for(auto* c:std::initializer_list<juce::Component*>{&proposalPanel,&measureSourceButton,&remeasureCurrentStateButton,&currentGain,&qcScope,&buildWaveformButton,&pageDetailBox})inspectorContent.addAndMakeVisible(c);
    inspectDomain.setTooltip("Inspect independently, or A to switch inspection and audition together while stopped.");
    auditionDomainBox.setTooltip("A switches inspection and audition together while stopped. Playback always requires Play.");
}
void StandaloneMainComponent::chooseMix()
{
    juce::PopupMenu menu;menu.addItem(1,"Interleaved",chooseSourceButton.isEnabled());menu.addItem(2,"Separate channel files",chooseProgrammeButton.isEnabled());
    menu.addItem(3,"Compose M&E + dialogue",chooseCompositionButton.isEnabled());
    menu.showMenuAsync(juce::PopupMenu::Options().withTargetComponent(&openMixButton),[safe=juce::Component::SafePointer<StandaloneMainComponent>(this)](int n)
    {if(!safe)return;if(n==1)safe->chooseSourceFile();else if(n==2)safe->chooseProgramme();else if(n==3)safe->chooseComposition();});
}
void StandaloneMainComponent::showTools()
{
    std::vector<juce::TextButton*> actions={&newSessionButton,&openSessionButton,&saveAsSessionButton,&relinkAudioButton,&relinkVideoButton,&cancelSessionButton,
        &clearSourceButton,&initializeCurrentStateButton,&clearCurrentStateButton,&lraPassagesButton,&preferencesButton,&inspectProposalButton,&sectionProposalButton,
        &applySectionProposalButton,&discardSectionProposalButton,&renderPcmButton,&publishPcmButton,&measurementReportButton,&artifactReportButton,&discardExportButton};
    juce::PopupMenu menu;int id=1;for(auto* button:actions)menu.addItem(id++,button->getButtonText(),button->isEnabled());
    menu.addSeparator();menu.addItem(100,"Reset layout");
    menu.showMenuAsync(juce::PopupMenu::Options().withTargetComponent(&toolsButton),[safe=juce::Component::SafePointer<StandaloneMainComponent>(this),actions](int n)
    {if(!safe)return;if(n==100){safe->timelineLayout.reset();safe->resized();}else if(n>0&&n<=int(actions.size())&&actions[size_t(n-1)]->isEnabled())actions[size_t(n-1)]->onClick();});
}
void StandaloneMainComponent::showAudioSetup()
{
    if(!audioWindow)
    {
        audioWindow=std::make_unique<ReferenceWindow>(audioPanel);
        audioWindow->setName("LevelScope - audio device and channel roles");
    }
    audioPanel.setVisible(true);audioWindow->setVisible(true);audioWindow->toFront(true);
}
void StandaloneMainComponent::showLraPassages(bool current)
{
    auto panel=std::make_unique<LraPassagesPanel>(session,current||inspectDomain.getSelectedId()==2?MeasurementDomain::currentState:MeasurementDomain::source);
    panel->onAssess=[this]{if(lraPassagesWindow)lraPassagesWindow->setVisible(false);proposeWithActiveMethod();refreshTimeline();};
    panel->onInspect=[this](MeasurementDomain d,int64_t begin,int64_t end)
    {inspectDomain.setSelectedId(d==MeasurementDomain::source?1:2,juce::dontSendNotification);refreshTimeline();timeline.highlight(double(begin),double(end));};
    lraPassagesWindow=std::make_unique<ProposalWindow>("LRA passages — inspection and candidate estimates",std::move(panel),this);
    lraPassagesWindow->present();
}
void StandaloneMainComponent::proposeWithActiveMethod()
{
    gainLane.cancelPending();if(sectionsPanel)sectionsPanel->cancelPendingGain();
    if(!jobRequirementsPanel.commitFields())return;
    if(proposalOperation.begin(session)){timeline.focusForTransport();if(proposalOperation.phase!=ProposalOperation::Phase::initialize)advanceProposalOperation();refreshFromSession();}
    else {refreshFromSession();refreshTimeline();}
}
void StandaloneMainComponent::advanceProposalOperation()
{
    proposalOperation.observe(session);
    if(!proposalOperation.active()||session.job.isActive()||analysisResult.valid()||sectionProposalResult.valid())return;
    using Phase=ProposalOperation::Phase;
    if(proposalOperation.cancelRequested)
    {
        proposalOperation.finish(session,Phase::cancelled,proposalOperation.applied?"Alternative switch cancelled - previous Current retained":"Propose cancelled - processing unchanged",
            proposalOperation.applied?"The previous applied choice and its independent measurement are unchanged. Choose another alternative or Restore before Propose.":"No proposal was applied.");return;
    }
    if(proposalOperation.phase==Phase::initialize){proposalOperation.prepare(session);return;}
    if(proposalOperation.phase==Phase::apply)
    {if(!proposalOperation.apply(session))return;}
    if(proposalOperation.phase==Phase::source)beginMeasurement(MeasurementDomain::source,true);
    else if(proposalOperation.phase==Phase::current||proposalOperation.phase==Phase::verify)beginMeasurement(MeasurementDomain::currentState,true);
    else if(proposalOperation.phase==Phase::search)beginSectionProposal(true);
    if(!session.job.isActive())proposalOperation.finish(session,Phase::failed,proposalOperation.applied?"Proposal applied - Current unverified":"Propose could not start",
        "The next stage could not start. Review programme roles and requirements, then retry or restore before Propose.");
}
void StandaloneMainComponent::showRegionPlanning()
{
    session.synchronizeSections();
    RegionPlanningPanel::showWindow(regionPlanningWindow,session,this,[this]{refreshFromSession();refreshTimeline();},
        [this]{proposeWithActiveMethod();refreshTimeline();});
}
void StandaloneMainComponent::reviewProposal(bool historical)
{
    if(!historical&&!session.sectionProposal&&proposalOperation.result)
    {
        const auto result=proposalOperation.selected?proposalOperation.selected:proposalOperation.result;const auto runtime=session.sessionId;const auto generation=proposalOperation.generation;
        auto status=[safe=juce::Component::SafePointer<StandaloneMainComponent>(this),runtime,generation]() -> juce::String
        {
            if(!safe||safe->session.sessionId!=runtime||safe->proposalOperation.generation!=generation)return "Historical inspection from an earlier operation; no Apply authority.";
            return safe->proposalOperation.title+". Browsing does not select or apply a plan. Choose an eligible alternative in the main window.";
        };
        std::unique_ptr<juce::Component> view;
        if(result->regionReview())
        {
            auto panel=std::make_unique<RegionProposalReviewPanel>(result,status);
            panel->onInspect=[this,runtime,generation](size_t,int64_t begin,int64_t end)
            {if(session.sessionId==runtime&&proposalOperation.generation==generation)timeline.highlight(double(begin),double(end));};
            view=std::move(panel);
        }
        else view=std::make_unique<ProposalInspectionPanel>(result->text(),status);
        proposalInspectionWindow=std::make_unique<ProposalWindow>("Proposal search - inspect all attempts",std::move(view),this);proposalInspectionWindow->present();return;
    }

    if(historical?!session.lastAppliedExplanation:!session.sectionProposal)return;
    const auto runtime=session.sessionId;const auto proposal=session.sectionProposal;const auto history=session.lastAppliedExplanation;
    auto binding=std::make_shared<std::shared_ptr<const SectionProposalResult>>(proposal);
    auto status=[safe=juce::Component::SafePointer<StandaloneMainComponent>(this),runtime,binding,history,historical]() -> juce::String
    {
        if(!safe||safe->session.sessionId!=runtime)return "Historical inspection snapshot from a previous session opening. Cannot Apply.";
        auto& s=safe->session;
        if(historical)return s.lastAppliedExplanation==history?s.appliedExplanationStatus():"An older applied explanation snapshot. A newer proposal has been applied.";
        if(s.sectionProposal!=*binding)return "Previous proposal snapshot — discarded, replaced or applied. Cannot Apply this snapshot.";
        return s.sectionProposalStale()?"STALE proposal snapshot — inputs changed. Apply unavailable. Generate again.":
            (*binding)->outcomeTitle()+". "+(s.canApplySectionProposal()&&!safe->gainLane.diagnosticRegionInspection()?"Explicit Apply is in the main window.":"Apply unavailable. Inspect the selected attempt and resolve any pending edit/job; Current verification stays separate.");
    };
    std::unique_ptr<juce::Component> view;
    if(!historical&&proposal->regionReview())
    {
        auto panel=std::make_unique<RegionProposalReviewPanel>(proposal,std::move(status));
        panel->onInspect=[this,runtime,binding](size_t attempt,int64_t begin,int64_t end)
        {
            if(session.sessionId!=runtime||session.sectionProposal!=*binding)return;
            gainLane.inspectRegionAttempt(attempt);timeline.highlight(double(begin),double(end));
            gainLane.selectFrame(begin,false);pageDetailBox.setText(gainLane.inspectionText(),false);refreshTimeline();
        };
        auto blocked=[safe=juce::Component::SafePointer<StandaloneMainComponent>(this),runtime,binding]() -> juce::String
        {
            if(!safe||safe->session.sessionId!=runtime||safe->session.sectionProposal!=*binding)return "Previous proposal/session snapshot. Reopen Review for the current result.";
            const auto& s=safe->session;
            if(s.pendingBoundaryEdits()||s.sectionTrimDraft)return "Finish or cancel pending manual edits before selecting an alternative.";
            if(s.audioEditsBlocked()||s.job.isActive()||s.persistenceLoading)return "Stop playback and finish or cancel the current job before selecting.";
            if(s.sectionProposalStale()||s.sourceContentInvalid)return "Stale inputs, requirements or preferences. Generate a fresh proposal.";
            return {};
        };
        panel->configureMeasuredSelection(eventAlternativeSelectionQualified,blocked,
            [safe=juce::Component::SafePointer<StandaloneMainComponent>(this),binding,blocked](const juce::String& resultId,const juce::String& attemptId) -> std::shared_ptr<const SectionProposalResult>
            {
                if(!safe||!safe->eventAlternativeSelectionQualified||blocked().isNotEmpty())return {};
                if(!safe->session.selectMeasuredRegionAlternative(resultId,attemptId))return {};
                *binding=safe->session.sectionProposal;safe->refreshFromSession();safe->refreshTimeline();return *binding;
            });
        view=std::move(panel);
    }
    else view=std::make_unique<ProposalInspectionPanel>(historical?session.appliedExplanationDisplayText():session.sectionProposalDisplayText(),std::move(status),historical&&session.lastAppliedExplanation->regionSettings?RegionProposalReviewPanel::pageBytes:4000);
    proposalInspectionWindow=std::make_unique<ProposalWindow>(historical?"Last applied proposal explanation":"Review complete proposal",std::move(view),this);
    proposalInspectionWindow->present();
}
void StandaloneMainComponent::showSections()
{
    if(!sectionsPanel)
    {
        sectionsPanel=std::make_unique<StandaloneSectionsPanel>(session);
        sectionsPanel->onChanged=[this]{refreshFromSession();refreshTimeline();};
        sectionsPanel->canImportMarkers=[this]{return canReplaceSession();};
        sectionsPanel->onPlay=[this]{playTimeline(&sectionsPanel->timelineView());};
        sectionsPanel->onStop=[this]{audition.stop();refreshTimeline();};
        auto& view=sectionsPanel->timelineView();view.setTransitionShading(timeline.transitionShading());
        view.onRequestSeek=timeline.onRequestSeek;
        view.onTogglePlayback=[this]{if(audition.active())sectionsPanel->onStop();else sectionsPanel->onPlay();};
        view.onToggleInspectionAudition=timeline.onToggleInspectionAudition;
        view.onUndo=[this](bool redo){return performCommand(redo?redoEdit:undoEdit);};
        view.onRemeasure=[this]{performCommand(measureCurrent);};
        view.peakVisibility().setToggleState(timeline.peakVisibility().getToggleState(),juce::dontSendNotification);
        view.peakVisibility().onClick=[this]{timeline.peakVisibility().setToggleState(sectionsPanel->timelineView().peakVisibility().getToggleState(),juce::dontSendNotification);timeline.repaint();sectionsPanel->timelineView().repaint();};
    }
    if(!sectionsWindow)sectionsWindow=std::make_unique<SectionsWindow>(*sectionsPanel);
    sectionsPanel->resume();
    sectionsWindow->setVisible(true);
    sectionsWindow->toFront(true);
    refreshTimeline();
}
void StandaloneMainComponent::playTimeline(MeasurementTimeline* origin)
{
    if(session.proposalOperationBusy)return;
    if(session.pendingBoundaryEdits()){jobLabel.setText(session.pendingBoundaryMessage(),juce::dontSendNotification);return;}
    if(session.job.isActive()||session.persistenceLoading)return;
    (origin?*origin:timeline).focusForTransport();
    if(!audition.connected()){showAudioSetup();return;}
    audioPanel.synchronize();const auto error=audition.play(session.auditionRequest(),session.sourceDocument.routing);
    if(error.isNotEmpty())jobLabel.setText(error,juce::dontSendNotification);refreshTimeline();
}
juce::StringArray StandaloneMainComponent::getMenuBarNames(){return {"File","Edit","View","Transport","Analysis"};}
juce::PopupMenu StandaloneMainComponent::getMenuForIndex(int index,const juce::String&)
{
    juce::PopupMenu menu;const bool idle=!session.audioEditsBlocked()&&!session.job.isActive()&&!session.persistenceLoading;
    const auto utf8=[](const char* text){return juce::String::fromUTF8(text);};
    if(index==0){menu.addItem(newDocument,"New session",newSessionButton.isEnabled());menu.addItem(openDocument,utf8("Open session…"),openSessionButton.isEnabled());menu.addItem(openAudio,utf8("Open mix… (⌘O)"),openMixButton.isEnabled());menu.addSeparator();menu.addItem(saveDocument,utf8("Save session (⌘S)"),saveSessionButton.isEnabled());menu.addItem(saveDocumentAs,utf8("Save session as…"),saveAsSessionButton.isEnabled());}
    if(index==1){menu.addItem(removeBoundaries,"Remove section boundaries only",idle&&session.sectionDocument.map().sections.size()>1);menu.addSeparator();auto* text=dynamic_cast<juce::TextEditor*>(juce::Component::getCurrentlyFocusedComponent());menu.addItem(undoEdit,utf8("Undo (⌘Z)"),text||session.canUndoManual());menu.addItem(redoEdit,utf8("Redo (⇧⌘Z)"),text||session.canUndoManual(true));}
    if(index==2){menu.addItem(transitionBands,"Gain transition shading (Current)",true,timeline.transitionShading());menu.addItem(peakWarnings,"Peak warnings",true,timeline.peakVisibility().getToggleState());menu.addItem(sectionsView,utf8("Sections and gain…"),sectionsButton.isEnabled());}
    if(index==3)menu.addItem(playStop,audition.active()?"Stop (Space)":"Play (Space)",timelineStop.isEnabled()||timelinePlay.isEnabled());
    if(index==4){menu.addItem(measureOriginal,"Measure Original",idle&&measureSourceButton.isEnabled());menu.addItem(measureCurrent,utf8("Re-measure Current State (⌘R)"),idle&&remeasureCurrentStateButton.isEnabled());menu.addSeparator();menu.addItem(proposalPreferences,utf8("Proposal limits…"),preferencesButton.isEnabled());}
    return menu;
}
bool StandaloneMainComponent::performCommand(int id)
{
    if(juce::Component::getCurrentlyModalComponent())return false;
    const auto cancelGainEntries=[this]{gainLane.cancelPending();if(sectionsPanel)sectionsPanel->cancelPendingGain();return true;};
    const auto click=[](juce::Button& button){if(!button.isEnabled()||!button.onClick)return false;button.onClick();return true;};
    switch(id)
    {
        case newDocument:return cancelGainEntries()&&click(newSessionButton);
        case openDocument:return cancelGainEntries()&&click(openSessionButton);
        case openAudio:return cancelGainEntries()&&click(openMixButton);
        case saveDocument:return cancelGainEntries()&&click(saveSessionButton);
        case saveDocumentAs:return cancelGainEntries()&&click(saveAsSessionButton);
        case undoEdit:case redoEdit:
            if(auto* text=dynamic_cast<juce::TextEditor*>(juce::Component::getCurrentlyFocusedComponent()))return id==undoEdit?text->undo():text->redo();
            if(!cancelGainEntries()||!session.undoManual(id==redoEdit))return false;
            refreshFromSession();refreshTimeline();if(sectionsPanel)sectionsPanel->synchronize();return true;
        case transitionBands:
            timeline.setTransitionShading(!timeline.transitionShading());
            if(sectionsPanel)sectionsPanel->timelineView().setTransitionShading(timeline.transitionShading());return true;
        case peakWarnings:
            timeline.peakVisibility().setToggleState(!timeline.peakVisibility().getToggleState(),juce::sendNotificationSync);
            if(sectionsPanel)sectionsPanel->timelineView().peakVisibility().setToggleState(timeline.peakVisibility().getToggleState(),juce::sendNotificationSync);return true;
        case removeBoundaries:
        {
            TimelineGainLane* selected=&gainLane;
            for(auto* c=juce::Component::getCurrentlyFocusedComponent();c;c=c->getParentComponent())
                if(auto* view=dynamic_cast<MeasurementTimeline*>(c);view&&view->attachedGainLane()){selected=view->attachedGainLane();break;}
            const auto changed=selected->removeSelectedBoundaries();if(changed){refreshFromSession();refreshTimeline();}return changed;
        }
        case sectionsView:return cancelGainEntries()&&click(sectionsButton);
        case playStop:if(audition.active()){audition.stop();refreshTimeline();return true;}return cancelGainEntries()&&click(timelinePlay);
        case measureOriginal:return cancelGainEntries()&&click(measureSourceButton);
        case measureCurrent:if(!cancelGainEntries())return false;return click(remeasureCurrentStateButton);
        case preferences:showAudioSetup();return true;
        case proposalPreferences:return click(preferencesButton);
        default:return false;
    }
}
bool StandaloneMainComponent::keyPressed(const juce::KeyPress& key)
{
    auto* focused=juce::Component::getCurrentlyFocusedComponent();
    for(auto* c=focused;c;c=c->getParentComponent())if(dynamic_cast<juce::TextEditor*>(c))return false;
    if(juce::Component::getCurrentlyModalComponent()!=nullptr)return false;
    if(key.getModifiers().isCommandDown())
    {
        if(key.getKeyCode()=='s'||key.getKeyCode()=='S')return performCommand(saveDocument);
        if(key.getKeyCode()=='z'||key.getKeyCode()=='Z')return performCommand(key.getModifiers().isShiftDown()?redoEdit:undoEdit);
        if(key.getKeyCode()=='r'||key.getKeyCode()=='R')return performCommand(measureCurrent);
        if(key.getKeyCode()==',')return performCommand(preferences);
        if(key.getKeyCode()=='o'||key.getKeyCode()=='O')return performCommand(openAudio);
    }
    return timeline.keyPressed(key);
}
void StandaloneMainComponent::toggleInspectionAudition()
{
    if(session.audioEditsBlocked()||session.job.isActive()||session.persistenceLoading)
    {pageDetailBox.setText("Stop playback/preparation or finish the job before switching inspection and audition (A).",false);return;}
    if(!session.sourceDocument.hasSource||!session.hasCurrentStateInitialized())return;
    const auto next=session.auditionDomain==MeasurementDomain::source?MeasurementDomain::currentState:MeasurementDomain::source;
    session.selectAuditionDomain(next);inspectDomain.setSelectedId(next==MeasurementDomain::source?1:2,juce::dontSendNotification);
    if(sectionsPanel)sectionsPanel->selectInspectionDomain(next);
    audioPanel.synchronize();timeline.clearHighlight();refreshTimeline();
    pageDetailBox.setText(inspectDomain.getText()+" + "+auditionDomainBox.getText()+". Stopped; press Play to listen.",false);
}
void StandaloneMainComponent::refreshTimeline()
{
    const auto domain=inspectDomain.getSelectedId()==2?MeasurementDomain::currentState:MeasurementDomain::source;
    const auto sourceBinding=WaveformBinding::from(session.measurementContext(MeasurementDomain::source),MeasurementDomain::source);
    const auto currentBinding=WaveformBinding::from(session.measurementContext(MeasurementDomain::currentState),MeasurementDomain::currentState);
    const bool valid=session.sourceDocument.hasSource && !session.sourceContentInvalid && session.sourceDocument.matchesFileStamps();
    gainLane.synchronize(valid&&session.sourceDocument.layout.confirmed?sourceBinding.runtime+sourceBinding.identity:juce::String{});
    if(gainLane.getHeight()!=gainLane.preferredHeight()){resized();timeline.resized();}
    proposalPanel.setDiagnosticInspection(gainLane.diagnosticRegionInspection());
    timeline.peakVisibility().onClick=[this]{if(sectionsPanel){sectionsPanel->timelineView().peakVisibility().setToggleState(timeline.peakVisibility().getToggleState(),juce::dontSendNotification);sectionsPanel->timelineView().repaint();}timeline.repaint();};
    timeline.onUndo=[this](bool redo){return performCommand(redo?redoEdit:undoEdit);};
    timeline.onRemeasure=[this]{performCommand(measureCurrent);};
    timeline.setProgramme(valid?session.sourceDocument.media.frames:0,valid?int(session.sourceDocument.media.sampleRate):1,valid?sourceBinding.runtime+sourceBinding.identity:juce::String());
    const auto currentPlan=session.currentStateDocument.processingPlan();
    waveformCache.synchronize(sourceBinding,currentBinding,valid,session.hasCurrentStateInitialized(),&currentPlan);
    auto view=waveformCache.view(domain==MeasurementDomain::source?sourceBinding:currentBinding);
    if(!view.data&&valid&&session.hasCurrentStateInitialized()&&domain==MeasurementDomain::currentState)view=waveformCache.gainPreview(sourceBinding);
    if(view.data&&domain==MeasurementDomain::currentState&&session.currentStateNeedsRemeasurement()&&view.accepted){view.previewBase=currentPlan;view.accepted=false;}
    timeline.setGainPreview(view.previewBase,currentPlan);
    juce::StringArray channels;for(int c=0;c<session.sourceDocument.layout.count;++c)
        channels.add(juce::AudioChannelSet::getAbbreviatedChannelTypeName(session.sourceDocument.layout.roles[size_t(c)]));
    timeline.setWaveform(view.data,view.accepted,channels);
    if(loudnessProgress&&(!valid||!session.job.isActive()||session.job.state==JobState::cancelling
        || !loudnessProgress->binding.sameJob(analysisDomain==MeasurementDomain::source?sourceBinding:currentBinding)))
    {loudnessProgress->abandon();loudnessProgress.reset();}
    timeline.setPreview(loudnessProgress?loudnessProgress->view(domain==MeasurementDomain::source?sourceBinding:currentBinding):nullptr);
    const auto& measured=domain==MeasurementDomain::source?session.source:session.currentState;
    const auto limits=session.getSelectedTargetProfile().limits;
    const auto peakLimit=limits.hasTruePeakCeiling?std::optional<double>(limits.truePeakCeilingDbtp):std::nullopt;
    timeline.setPeakTarget(peakLimit,valid&&measured.hasMeasurement&&(domain==MeasurementDomain::source||!session.currentStateNeedsRemeasurement()));
    buildWaveformButton.setEnabled(valid && measured.hasMeasurement && !view.data && !session.job.isActive() && !session.audioEditsBlocked());
    timelineDiagnostics.select(measured.analysis,session.measurementContext(domain),domain);timelineDiagnostics.poll();
    const auto result=timelineDiagnostics.result();auto analysis=result?result->analysis:nullptr;
    timeline.setAnalysis(analysis);
    const auto original=session.source.hasMeasurement&&!session.sourceContentInvalid?session.source.analysis:nullptr;
    const auto reference=original?std::shared_ptr<const measurement::Result>(original,&original->detail):nullptr;
    timeline.setReferenceMeasurement(reference,domain==MeasurementDomain::currentState);
    auditionDomainBox.setSelectedId(session.auditionDomain==MeasurementDomain::source?1:2,juce::dontSendNotification);
    auditionDomainBox.setEnabled(!session.audioEditsBlocked()&&!session.job.isActive());
    timelinePlay.setEnabled(!session.proposalOperationBusy&&!session.pendingBoundaryEdits()&&!session.job.isActive()&&!session.persistenceLoading&&!audition.active()&&session.sourceDocument.hasSource);
    timelineStop.setEnabled(audition.active());
    jobRequirementsPanel.refreshAvailability();
    proposalPanel.refresh(importControl?importControl->progress.load():-1);
    auto seekReason=TimelineTransport::blockedReason(session,audition);
    if(seekReason.isEmpty()&&audition.boundIdentity().isEmpty())seekReason="Waiting for the current audition binding.";
    timeline.setTransport(audition.positionFrames(),seekReason.isEmpty(),seekReason,audition.boundIdentity());
    timeline.refreshSectionEditing();
    if(sectionsWindow && sectionsWindow->isVisible())
    {
        auto& editor=sectionsPanel->timelineView();
        const auto& binding=sectionsPanel->inspectionDomain()==MeasurementDomain::source?sourceBinding:currentBinding;
        auto sectionWave=waveformCache.view(binding);
        if(!sectionWave.data&&valid&&session.hasCurrentStateInitialized()&&sectionsPanel->inspectionDomain()==MeasurementDomain::currentState)sectionWave=waveformCache.gainPreview(sourceBinding);
        if(sectionWave.data&&sectionsPanel->inspectionDomain()==MeasurementDomain::currentState&&session.currentStateNeedsRemeasurement()&&sectionWave.accepted){sectionWave.previewBase=currentPlan;sectionWave.accepted=false;}
        editor.setProgramme(valid?session.sourceDocument.media.frames:0,valid?int(session.sourceDocument.media.sampleRate):1,
                            valid?sourceBinding.runtime+sourceBinding.identity:juce::String());
        editor.setGainPreview(sectionWave.previewBase,currentPlan);
        sectionsPanel->refreshMeasuredView();
        if(sectionsPanel->inspectionDomain()==MeasurementDomain::currentState&&session.currentStateNeedsRemeasurement())editor.setAnalysis(nullptr);
        editor.setWaveform(sectionWave.data,sectionWave.accepted,channels);
        editor.setPreview(loudnessProgress?loudnessProgress->view(binding):nullptr);
        editor.setReferenceMeasurement(reference,sectionsPanel->inspectionDomain()==MeasurementDomain::currentState);
        editor.setPeakTarget(peakLimit,valid&&(sectionsPanel->inspectionDomain()==MeasurementDomain::source?session.source.hasMeasurement:session.currentState.hasMeasurement&&!session.currentStateNeedsRemeasurement()));
        editor.refreshSectionEditing();
        editor.setTransport(audition.positionFrames(),seekReason.isEmpty(),seekReason,audition.boundIdentity());
        sectionsPanel->setPlaybackState(timelinePlay.isEnabled(),timelineStop.isEnabled(),session.auditionDomain==MeasurementDomain::source?"Original":"Current");
    }
    positionLabel.setText(juce::String(audition.positionSeconds(),3)+" s | "+(audition.playbackReady()?"Playback ready":audition.preparing()?"Preparing playback":audition.connected()?"Device open":"Device closed"),juce::dontSendNotification);
    positionLabel.setTooltip(audition.status()+"\n"+(audition.connected()?"Output connected; Play is explicit.":"Output disconnected; offline work remains available."));
    const auto warning=audition.monitorOverloadWarning();
    const bool warningVisible=warning.isNotEmpty(),visibilityChanged=monitorWarning.isVisible()!=warningVisible;
    monitorWarning.setText(warning,juce::dontSendNotification);monitorWarning.setTooltip(warning);
    monitorWarning.setVisible(warningVisible);monitorSetup.setVisible(warningVisible);
    if(visibilityChanged)resized();
}
void StandaloneMainComponent::refreshJobFeedback()
{
    auto control = session.job.isActive() ? importControl : nullptr;
    bool nestedImport = false;
    if(programmePanel && programmePanel->preparationControl()) { control = programmePanel->preparationControl(); nestedImport = true; }
    else if(compositionPanel && compositionPanel->preparationControl()) { control = compositionPanel->preparationControl(); nestedImport = true; }
    const bool cancelling = (proposalOperation.cancelRequested&&proposalOperation.active()) || (control && control->isCancelled()) || (!nestedImport && session.job.state == JobState::cancelling);
    const bool active = nestedImport || session.job.isActive() || proposalOperation.active();
    const bool phased = control && control->feedback;
    auto text = phased ? describeJobFeedback(control->feedback->snapshot(), !nestedImport && session.job.stage == JobStage::analysis,
                                            analysisDomain == MeasurementDomain::currentState, cancelling)
                       : JobFeedbackText {session.job.message, {}, active && control ? double(control->progress.load()) : -1};
    if(!active && session.job.state == JobState::cancelled && (session.job.stage == JobStage::analysis || session.job.stage == JobStage::ingest))
        text = {session.job.message, session.job.stage == JobStage::ingest ? "Previous source retained. Use Open Mix to retry import." : "No incomplete measurement accepted. Use Measure Source or Re-measure Current State when ready.", -1};
    else if(!active && session.job.state == JobState::failed && (session.job.stage == JobStage::analysis || session.job.stage == JobStage::ingest))
        text.detail = "No new completed measurement accepted. Review the reason above and retry when ready.";
    else if(!active && (programmePanel || compositionPanel))
        text = {"Import: review assignments and alignment", "Continue in the import window to confirm the programme, or cancel to keep the previous source.", -1};
    else if(!active && session.sourceDocument.hasSource && !session.sourceDocument.layout.confirmed)
        text = {"Source validated: confirm channel roles", "Open Audio setup to review the proposed channel order, then Confirm roles and analyse.", -1};
    else if(!active && !session.sourceDocument.hasSource && session.job.state == JobState::idle)
        text = {"Ready to import", "Open Mix to begin. Output connection is optional for offline analysis.", -1};
    if(!active&&session.boundaryBatch)text={juce::String(session.pendingBoundaryCount())+(session.pendingBoundaryCount()==1?" boundary change pending":" boundary changes pending"), "Enter: review Move gain with sections / Cancel. Esc: cancel. Zoom and pan remain available; playback and measurement wait for this edit.", -1};
    if(!active&&session.appliedBoundaryDrag)text={"Moving applied region boundary",session.pendingBoundaryMessage()+" The previous measurement describes the unchanged applied audio until release.",-1};
    if(proposalOperation.phase!=ProposalOperation::Phase::idle&&proposalOperation.phase!=ProposalOperation::Phase::stale)
    {
        if(proposalOperation.active())text.detail=text.title+". "+text.detail;
        else text.detail=proposalOperation.detail.replace("\n","  ");
        text.title=proposalOperation.title;
    }
    const bool playbackPreparation=!active && (audition.preparing() || audition.getState()==AuditionState::buffering);
    const bool playbackCancelling=playbackPreparation && audition.getState()==AuditionState::cancelling;
    if(playbackPreparation)
        text={audition.getState()==AuditionState::buffering?"Preparing playback output":audition.preparationStage(),
            audition.backgroundPreparing()?"Playback remains stopped. Play will use this preparation; Cancel leaves offline analysis available.":audition.status(),
            audition.preparing()?double(audition.preparationProgress()):-1};
    jobLabel.setText(text.title, juce::dontSendNotification);jobLabel.setTooltip(text.title);
    jobDetail.setText(text.detail, juce::dontSendNotification);jobDetail.setTooltip(text.detail);
    jobFraction = text.fraction;jobProgress.setVisible(active || playbackPreparation);
    jobProgress.setTextToDisplay(jobFraction < 0 ? "Working..." : juce::String(juce::roundToInt(jobFraction * 100)) + (phased ? "% this pass" : "%"));
    cancelJobButton.setEnabled((active && !cancelling) || (playbackPreparation && !playbackCancelling));
    cancelJobButton.setButtonText(cancelling || playbackCancelling ? "Cancelling..." : playbackPreparation ? "Cancel playback" : "Cancel job");
}
void StandaloneMainComponent::resized()
{
    auto bounds=getLocalBounds().reduced(16);
    auto top=bounds.removeFromTop(32);
    openMixButton.setBounds(top.removeFromLeft(125).reduced(2));videoButton.setBounds(top.removeFromLeft(155).reduced(2));
    saveSessionButton.setBounds(top.removeFromLeft(115).reduced(2));toolsButton.setBounds(top.removeFromLeft(160).reduced(2));
    audioSetupButton.setBounds(top.removeFromRight(150).reduced(2));top.reduce(5,0);sourceFilename.setBounds(top.removeFromTop(17));sessionStatus.setBounds(top);bounds.removeFromTop(7);
    auto feedback = bounds.removeFromTop(56);
    cancelJobButton.setBounds(feedback.removeFromRight(110).withSizeKeepingCentre(106,28));
    auto progressArea = feedback.removeFromRight(165);jobProgress.setBounds(progressArea.reduced(8,18));
    jobLabel.setBounds(feedback.removeFromTop(26));jobDetail.setBounds(feedback);
    jobRequirementsPanel.setBounds(bounds.removeFromTop(timelineLayout.comparison));comparisonDivider.setBounds(bounds.removeFromTop(TimelineLayout::dividerHeight));
    jobRequirementsPanel.informationLine().setBounds(bounds.removeFromTop(TimelineLayout::informationHeight));
    auto transport=bounds.removeFromTop(32);
    inspectDomain.setBounds(transport.removeFromLeft(155).reduced(2));auditionDomainBox.setBounds(transport.removeFromLeft(165).reduced(2));
    timelinePlay.setBounds(transport.removeFromLeft(65).reduced(2));timelineStop.setBounds(transport.removeFromLeft(65).reduced(2));
    audioSetupButton.setTooltip("Choose this Mac's output and optional automatic reconnection; review channel roles and routing. Opening always leaves playback stopped.");
    sectionsButton.setBounds(transport.removeFromRight(160).reduced(2));diagnosticsButton.setBounds(transport.removeFromRight(190).reduced(2));positionLabel.setBounds(transport);
    if(monitorWarning.isVisible())
    {
        auto warning=bounds.removeFromTop(64).reduced(0,2);
        monitorSetup.setBounds(warning.removeFromRight(145).reduced(4,16));monitorWarning.setBounds(warning);
    }
    timelineAvailableHeight=bounds.getHeight();timelineLayout.moveTimeline(0,timelineAvailableHeight,session.sourceDocument.layout.count,gainLane.preferredHeight());
    auto bottom=bounds.removeFromBottom(timelineLayout.inspector);
    waveformDivider.setBounds(bounds.removeFromBottom(TimelineLayout::dividerHeight));
    inspectorViewport.setBounds(bottom);inspectorContent.setSize(std::max(1,bottom.getWidth()-inspectorViewport.getScrollBarThickness()),std::max(TimelineLayout::defaultInspector,bottom.getHeight()));
    bottom=inspectorContent.getLocalBounds();
    proposalPanel.setBounds(bottom.removeFromTop(ProposalWorkflowPanel::height));
    auto processing=bottom.removeFromTop(32);
    measureSourceButton.setBounds(processing.removeFromLeft(150).reduced(2));remeasureCurrentStateButton.setBounds(processing.removeFromLeft(190).reduced(2));
    currentGain.setBounds(processing.removeFromLeft(300).reduced(2));
    auto applicability=bottom.removeFromTop(26);buildWaveformButton.setBounds(applicability.removeFromRight(150));qcScope.setBounds(applicability);
    pageDetailBox.setBounds(bottom.reduced(0,3));timeline.setBounds(bounds.reduced(0,4));
}

// [END LS-STANDALONE-WORKFLOW-SHELL]
}

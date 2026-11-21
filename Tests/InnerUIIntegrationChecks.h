// Native callbacks, private drafts and copied programme data. Never opens output hardware.
template<class T> T* childOfType(juce::Component& c)
{if(auto* result=dynamic_cast<T*>(&c))return result;for(auto* child:c.getChildren())if(auto* result=childOfType<T>(*child))return result;return nullptr;}
void waitFile(StandaloneMainComponent& ui)
{
    const auto deadline=juce::Time::getMillisecondCounterHiRes()+90000;
    do{Access::tick(ui);juce::MessageManager::getInstance()->runDispatchLoopUntil(5);}while(Access::filePending(ui)&&juce::Time::getMillisecondCounterHiRes()<deadline);
    expect(!Access::filePending(ui),"native file operation completes within focused timeout");if(Access::model(ui).job.isActive())waitJob(ui);Access::refresh(ui);
}
juce::Component& modal()
{
    const auto deadline=juce::Time::getMillisecondCounterHiRes()+3000;
    while(!juce::Component::getCurrentlyModalComponent()&&juce::Time::getMillisecondCounterHiRes()<deadline)juce::MessageManager::getInstance()->runDispatchLoopUntil(5);
    auto* result=juce::Component::getCurrentlyModalComponent();expect(result!=nullptr,"expected native modal is visible");return *result;
}
void pressModal(const juce::String& button)
{control<juce::TextButton>(modal(),button).triggerClick();juce::MessageManager::getInstance()->runDispatchLoopUntil(40);}
void finishSaveAs(StandaloneMainComponent& ui,const juce::File& destination)
{
    auto& picker=modal();auto* browser=childOfType<juce::FileBrowserComponent>(picker);expect(browser!=nullptr,"Save As uses real chooser callback");
    browser->setFileName(destination.getFullPathName());juce::MessageManager::getInstance()->runDispatchLoopUntil(50);
    pressModal("Save");waitFile(ui);expect(destination.existsAsFile()&&Access::model(ui).documentFile==destination,"native Save As commits selected file only after successful save");
}
ProposalInspectionPanel& visibleReport(RegionProposalReviewPanel& panel)
{for(auto* c:panel.getChildren())if(auto* p=dynamic_cast<ProposalInspectionPanel*>(c);p&&p->isVisible())return *p;throw std::runtime_error("No visible report");}

void innerEditorRoundTrips(StandaloneMainComponent& ui,const juce::File& root)
{
    auto& s=Access::model(ui);const auto input=sectionInputIdentity(s.measurementContext(MeasurementDomain::source));juce::String error;
    for(int schema=1;schema<=3;++schema)
    {
        auto settings=schema==3?RegionPlanningSettings::innerDefaults(RegionMode::selection,input):schema==2?RegionPlanningSettings::eventDefaults(RegionMode::selection,input):RegionPlanningSettings::defaults(RegionMode::selection,input);
        settings.selection={480001,2400011};settings.allowGlobal=true;settings.leftMarginMs=500;settings.rightMarginMs=1500;
        settings.start=RegionEdgeRange{470000,480000};settings.internalCuts={{1200001,1200001}};settings.protectedLocal={{1400001,1500007}};
        if(schema>=2)settings.eventFeatures.hopMs=5;
        expect(s.setRegionPlanning(settings,error),"install valid legacy/current full definition");Access::refresh(ui);const auto digest=settings.digest();
        Access::settings(ui);auto* panel=Access::settingsPanel(ui);expect(panel!=nullptr,"main opens proposal settings");std::optional<RegionPlanningSettings> draft;
        expect(control<juce::ComboBox>(*panel,"Proposal method").getSelectedId()==schema+1&&panel->readDraft(draft).isEmpty()&&draft->digest()==digest,"opening each schema preserves complete definition, pins and protections");
        text(*panel,"Left margin seconds (0-3, step 0.1)","0.9");click(*panel,"Cancel");expect(s.regionPlanning->digest()==digest,"Cancel does not install private margin edits");
        Access::settings(ui);panel=Access::settingsPanel(ui);expect(panel->readDraft(draft).isEmpty()&&draft->digest()==digest,"reopening cancelled editor reloads full active definition");
        click(*panel,"Confirm settings");expect(s.regionPlanning->digest()==digest,"confirming untouched schema preserves exact digest");Access::closeSettings(ui);
        const auto file=root.getChildFile("schema-"+juce::String(schema)+".levelscope");Access::save(ui,true);finishSaveAs(ui,file);
        JobControl job;auto read=sf::read(file,job);expect(read.document&&read.document->loadedSchema==(schema==3?11:10)&&read.document->regionPlanning->digest()==digest,"saved definition uses prospective format and exact settings");
        Access::open(ui,file);waitFile(ui);expect(s.regionPlanning&&s.regionPlanning->digest()==digest&&Access::stopped(ui),"native reopen restores exact settings and stays stopped");
    }
    Access::open(ui,root.getChildFile("schema-2.levelscope"));waitFile(ui);expect(s.openedSessionSchema==10,"native older-format source opened before real upgrade");
    Access::settings(ui);auto* panel=Access::settingsPanel(ui);choose(*panel,"Proposal method",1);click(*panel,"Confirm settings");expect(!s.regionPlanning,"standard workflow explicitly available");
    choose(*panel,"Proposal method",4);choose(*panel,"Region mode",1);std::optional<RegionPlanningSettings> draft;
    expect(panel->readDraft(draft).isEmpty()&&draft->schema==3&&draft->eventFeatures.innerEntrances&&!draft->allowAutomaticMargins&&draft->leftMarginMs==0&&draft->rightMarginMs==0,"explicit competing activation retains conservative unpermitted automatic margins");
    toggle(*panel,"Allow automatic outward margins",true);text(*panel,"Left margin seconds (0-3, step 0.1)","0.9");text(*panel,"Right margin seconds (0-3, step 0.1)","1.1");
    expect(panel->readDraft(draft).isEmpty()&&draft->leftMarginMs==900&&draft->rightMarginMs==1100&&draft->allowAutomaticMargins,"independent margins retain independent permission");click(*panel,"Confirm settings");Access::closeSettings(ui);
    auto presentation=SessionSavePresentation::forDocument(sf::capture(s),10,s.documentFile);expect(presentation.writer==11&&presentation.upgradesExisting&&presentation.message.contains("format 11")&&presentation.suggestedCopy.getFileName().contains("format-11"),"native compatibility predicts schema11 from unsaved settings");
    const auto original=s.documentFile;const auto hash=juce::SHA256(original).toHexString();
    Access::save(ui);auto* alert=dynamic_cast<juce::AlertWindow*>(&modal());expect(alert&&alert->getName()=="Session format upgrade","native prospective compatibility uses format upgrade dialog");imageOf(*alert,root.getChildFile("format-11-upgrade.png"));pressModal("Cancel");expect(!Access::filePending(ui)&&juce::SHA256(original).toHexString()==hash,"cancel upgrade does no file IO");
    Access::save(ui);pressModal("Save As");finishSaveAs(ui,root.getChildFile("upgraded-copy.levelscope"));expect(juce::SHA256(original).toHexString()==hash,"upgrade Save As preserves original bytes");
    Access::settings(ui);panel=Access::settingsPanel(ui);imageOf(*panel,root.getChildFile("competing-settings.png"));Access::closeSettings(ui);
    const auto retained=s;const auto destination=s.documentFile;const auto savedHash=juce::SHA256(destination).toHexString();s.openedSessionSchema=10;Access::save(ui);modal();s.newDocument();pressModal("Save");
    expect(!Access::filePending(ui)&&!s.sourceDocument.hasSource&&juce::SHA256(destination).toHexString()==savedHash,"old format prompt cannot save into New runtime");s=retained;Access::refresh(ui);
}

void standardCapabilityChecks(StandaloneMainComponent& ui,const juce::File& root,const StandaloneSessionModel& stereo)
{
    auto& s=Access::model(ui);const auto busy=s.auditionBusy;juce::String error;
    for(int channels:{2,1})
    {
        s=channels==2?stereo:programme(root.getChildFile("mono"),1);s.auditionBusy=busy;
        auto policy=s.proposalPolicy;policy.durationThousandths=1000;expect(s.setProposalPolicy(policy,error),"explicit global allowance for compatibility case");
        for(bool withLra:{false,true})
        {
            auto requirements=s.jobRequirements;requirements.maximumLraTenths=withLra?std::optional<int>(200):std::nullopt;expect(s.setJobRequirements(requirements,error),"standard no-LRA / global-only target");Access::refresh(ui);click(Access::proposal(ui),"Propose");
            expect(s.job.isActive()||bool(s.sectionProposal),"single standard Propose uses existing immediate global or asynchronous path without interstitial");waitJob(ui);
            expect(!s.regionPlanning&&s.sectionProposal,"standard mono/stereo returns a proposal without experimental opt-in");
            if(channels==1&&withLra)
            {JobControl control;const auto expected=generateSectionProposal(s.sectionProposalRequest(),control);expect(expected.value&&s.sectionProposal->status()==(*expected.value)->status()&&s.sectionProposal->text().contains("stereo and conventional 5.1 only"),"unsupported mono LRA retains existing service outcome and visible capability explanation");}
            else expect(s.sectionProposal->status()==SectionProposalStatus::globalOnly||s.sectionProposal->status()==SectionProposalStatus::ready||s.sectionProposal->status()==SectionProposalStatus::noChange,"supported standard no-LRA/global-only stereo and mono jobs remain available");
            expect(s.sectionProposal->methodVersion()!=innerSearchVersion,"standard dispatch does not silently opt into competing method");
        }
    }
    s=stereo;s.auditionBusy=busy;Access::refresh(ui);
}

void innerUIChecks(const juce::File& root,const juce::File& supplied)
{
    auto& look=juce::LookAndFeel::getDefaultLookAndFeel();const bool nativeAlerts=look.isUsingNativeAlertWindows();look.setUsingNativeAlertWindows(false);
    struct Restore{juce::LookAndFeel& l;bool native;~Restore(){l.setUsingNativeAlertWindows(native);}} restoreLook{look,nativeAlerts};
    StandaloneMainComponent ui(emptyOutputCatalog(),noOutputSettings());Access::stopTimer(ui);Access::deterministicChooser(ui);ui.setSize(1320,1100);ui.addToDesktop(juce::ComponentPeer::windowHasTitleBar);ui.setVisible(true);
    juce::String errorForInner;auto& s=Access::model(ui);const auto busy=s.auditionBusy;auto fixture=programme(root.getChildFile("fixture"),2);s=fixture;s.auditionBusy=busy;Access::refresh(ui);
    innerEditorRoundTrips(ui,root);standardCapabilityChecks(ui,root,fixture);
    expect(supplied.existsAsFile(),"explicit supplied passing session exists");const auto originalHash=juce::SHA256(supplied).toHexString();JobControl job;auto read=sf::read(supplied,job);
    expect(read.document&&read.document->regionPlanning&&read.document->regionPlanning->schema==3,"supplied checkpoint has active competing settings");
    const juce::File originalAudio(read.document->source.fullPath);expect(originalAudio.existsAsFile(),"supplied audio locator resolves");const auto audioHash=juce::SHA256(originalAudio).toHexString();const auto copyAudio=root.getChildFile("supplied-audio.wav");expect(originalAudio.copyFileTo(copyAudio),"copy audio before independent native checks");
    expect(sf::relinkAudio(*read.document,copyAudio),"relink owned session copy to owned identical audio");const auto copied=root.getChildFile("supplied-copy.levelscope");expect(sf::save(*read.document,copied,job).succeeded(),"write owned source session copy");
    Access::open(ui,copied);waitFile(ui);expect(s.regionPlanning&&s.regionPlanning->schema==3&&s.regionPlanning->allowAutomaticMargins&&s.regionPlanning->leftMarginMs==1000&&s.regionPlanning->rightMarginMs==1000,"native supplied session restores actual automatic 1s permissions");
    expect(s.lastAppliedExplanation&&!s.sectionProposal&&Access::stopped(ui),"reopened history has no live selection or playback authority");const auto historical=encodeAppliedExplanation(*s.lastAppliedExplanation);const auto processing=proposalPlanIdentity(s.currentStateDocument.processingPlan());
    Access::settings(ui);auto* editor=Access::settingsPanel(ui);std::optional<RegionPlanningSettings> draft;expect(editor->readDraft(draft).isEmpty()&&draft->digest()==s.regionPlanning->digest(),"supplied v3 draft no longer downgrades method");
    text(*editor,"Left margin seconds (0-3, step 0.1)","0.9");text(*editor,"Right margin seconds (0-3, step 0.1)","1.1");imageOf(*editor,root.getChildFile("supplied-editor.png"));click(*editor,"Confirm settings");expect(s.regionPlanning->leftMarginMs==900&&s.regionPlanning->rightMarginMs==1100,"actual supplied settings support independent margin confirmation");
    click(*editor,"Discard edits / reload");text(*editor,"Left margin seconds (0-3, step 0.1)","1.0");text(*editor,"Right margin seconds (0-3, step 0.1)","1.0");click(*editor,"Confirm settings");Access::closeSettings(ui);Access::refresh(ui);
    if(!s.source.hasMeasurement){click(ui,"Measure Source");waitJob(ui);}expect(s.source.hasMeasurement,"supplied Source independently verified for integrated build");
    click(ui,"Re-measure Current State");waitJob(ui);expect(s.currentState.hasMeasurement&&s.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass&&proposalPlanIdentity(s.currentStateDocument.processingPlan())==processing,"independent supplied Current verifies unchanged applied plan");
    Access::save(ui);waitFile(ui);
    auto historyOnly=sf::capture(s);historyOnly.regionPlanning.reset();auto presentation=SessionSavePresentation::forDocument(historyOnly,10,copied);expect(presentation.writer==11&&presentation.upgradesExisting,"schema3 retained only in applied history still requires document11");
    expect(s.setRegionPlanning({},errorForInner),"temporarily use standard settings with schema3 history");Access::refresh(ui);Access::save(ui,true);const auto historyCopy=root.getChildFile("history-only.levelscope");finishSaveAs(ui,historyCopy);Access::open(ui,historyCopy);waitFile(ui);
    expect(!s.regionPlanning&&s.openedSessionSchema==11&&s.lastAppliedExplanation&&encodeAppliedExplanation(*s.lastAppliedExplanation)==historical,"native history-only Save As/reopen preserves full explanation and writer11");
    Access::open(ui,copied);waitFile(ui);auto neutral=s.manualState();neutral.gains.reset();neutral.globalDb=0;expect(s.commitManualEdit(neutral,"Neutral test baseline"),"owned copied programme reset for fresh reassessment");Access::refresh(ui);click(ui,"Re-measure Current State");waitJob(ui);
    const auto baseline=s.currentState.analysis;const auto source=s.source.analysis;const auto revision=s.revision;Heartbeat heartbeat;heartbeat.begin();click(Access::proposal(ui),"Propose");waitJob(ui);
    expect(s.sectionProposal&&s.sectionProposal->methodVersion()==innerSearchVersion&&s.sectionProposal->status()==SectionProposalStatus::reviewRequired,"native competing reassessment of supplied baseline requires review");
    const auto result=s.sectionProposal;const auto review=result->regionReview();expect(review&&review->eventEvidence&&review->attempts.size()==3&&!review->selectedAttempt&&!review->recommendedAttempt,"bounded three exact attempts preserve unresolved/partial interpretation");
    expect(review->constraints.settings.digest()==s.regionPlanning->digest()&&review->constraints.settings.allowAutomaticMargins&&review->constraints.settings.leftMarginMs==1000&&review->constraints.settings.rightMarginMs==1000,"actual reassessment binds the full confirmed definition and both explicit margin permissions");
    for(const auto& r:review->constraints.roots)expect(r.allowed.begin>=std::max<int64_t>(0,r.hint.begin-48000)&&r.allowed.end<=std::min<int64_t>(review->constraints.frames,r.hint.end+48000),"effective root area remains clipped inside explicit per-side permission");
    const auto pass=std::find_if(review->attempts.begin(),review->attempts.end(),[](const auto& a){return a.evaluation.verdict==Verdict::pass;});expect(pass!=review->attempts.end(),"supplied controlled comparison produces an owned measured pass");const auto passing=size_t(pass-review->attempts.begin());
    Access::review(ui);auto* panel=Access::reviewPanel(ui);expect(panel!=nullptr,"main opens competing review");
    for(size_t n=0;n<review->attempts.size();++n)
    {choose(*panel,"Measured region attempt",int(n+2));expect(visibleReport(*panel).allText().contains("ACTUAL measured rule result")&&visibleReport(*panel).allText().contains("Boundary choices by root"),"every actual attempt includes typed boundary family and actual results");expect(control<juce::TextButton>(*panel,"Use this measured alternative").isEnabled()==(review->attempts[n].evaluation.verdict==Verdict::pass&&bool(review->attempts[n].receipt)),"selection available only for owned passing exact result");}
    const auto inspected=Access::lane(ui).diagnosticRegionInspection();
    choose(*panel,"Measured region attempt",RegionProposalReviewPanel::sharedEvidenceId);expect(visibleReport(*panel).allText().contains(review->eventEvidence->innerDescription())&&!control<juce::TextButton>(*panel,"Use this measured alternative").isEnabled(),"all immutable rival witnesses accessible without selection authority");imageOf(*panel,root.getChildFile("shared-boundary-evidence.png"));
    choose(*panel,"Measured region attempt",RegionProposalReviewPanel::proxyComparisonsId);for(const auto& restriction:review->restrictions)expect(visibleReport(*panel).allText().contains(restriction),"every omitted/proxy comparison remains accessible");
    expect(!control<juce::TextButton>(*panel,"Use this measured alternative").isEnabled()&&s.sectionProposal==result&&Access::lane(ui).diagnosticRegionInspection()==inspected,"documentary browsing never changes selection or timeline plan");imageOf(*panel,root.getChildFile("proxy-comparisons.png"));
    std::vector<double> switches;for(int n=0;n<30;++n){const auto start=juce::Time::getMillisecondCounterHiRes();choose(*panel,"Measured region attempt",n%3==0?RegionProposalReviewPanel::sharedEvidenceId:n%3==1?RegionProposalReviewPanel::proxyComparisonsId:int(passing+2));switches.push_back(juce::Time::getMillisecondCounterHiRes()-start);juce::MessageManager::getInstance()->runDispatchLoopUntil(20);}std::sort(switches.begin(),switches.end());expect(switches[28]<=100,"competing evidence/review switching p95 <=100ms");auto timing=juce::var(new juce::DynamicObject());timing.getDynamicObject()->setProperty("inner_review_switch_p95_ms",switches[28]);cases.add(timing);heartbeat.finish();
    choose(*panel,"Measured region attempt",int(passing+2));expect(!s.selectMeasuredRegionAlternative("foreign-result",pass->attemptIdentity),"foreign result cannot select measured pass");
    auto policy=s.proposalPolicy;policy.durationThousandths=600;expect(s.setProposalPolicy(policy,errorForInner),"temporary real preference change");panel->refreshSelection();expect(!control<juce::TextButton>(*panel,"Use this measured alternative").isEnabled(),"changed competing preferences make results stale");policy.durationThousandths=200;expect(s.setProposalPolicy(policy,errorForInner),"restore exact assessed policy");panel->refreshSelection();
    click(*panel,"Use this measured alternative");expect(s.sectionProposal->explicitlySelectedAlternative()&&s.currentState.analysis==baseline&&s.source.analysis==source&&s.revision==revision,"explicit actual selection does not change Current or revision");
    const auto selected=s.sectionProposal;
    const auto unchanged=s;for(int guard=0;guard<3;++guard){if(guard==0)s.sectionTrimDraft.emplace();if(guard==1)s.boundaryBatch=BoundaryEditBatch{s.manualState(),s.manualState(),s.sessionId,s.revision};if(guard==2)s.sessionId="replacement-runtime";Access::refresh(ui);expect(!control<juce::TextButton>(Access::proposal(ui),"Apply proposal").isEnabled()&&!s.applySectionProposal(),"pending edit and replacement runtime guards reject competing Apply");s=unchanged;Access::refresh(ui);}
    Access::lane(ui).inspectRegionAttempt((passing+1)%review->attempts.size());Access::refresh(ui);expect(!control<juce::TextButton>(Access::proposal(ui),"Apply proposal").isEnabled(),"Apply requires currently inspected matching selected plan");click(*panel,"Inspect in timeline");Access::refresh(ui);
    s.auditionBusy=[]{return true;};Access::refresh(ui);expect(!control<juce::TextButton>(Access::proposal(ui),"Apply proposal").isEnabled()&&!s.applySectionProposal(),"competing Apply requires stopped transport");s.auditionBusy=busy;Access::refresh(ui);imageOf(*panel,root.getChildFile("competing-selected.png"));imageOf(ui,root.getChildFile("competing-main.png"));
    click(Access::proposal(ui),"Apply proposal");expect(s.currentStateNeedsRemeasurement()&&!s.currentState.hasMeasurement,"separate Apply leaves Current unverified");const auto applied=proposalPlanIdentity(s.currentStateDocument.processingPlan());
    click(ui,"Re-measure Current State");const auto cancelAt=juce::Time::getMillisecondCounterHiRes();Access::cancel(ui);waitJob(ui);expect(juce::Time::getMillisecondCounterHiRes()-cancelAt<=1000&&!s.currentState.hasMeasurement&&proposalPlanIdentity(s.currentStateDocument.processingPlan())==applied,"cancel measurement promptly without losing applied processing");
    click(ui,"Re-measure Current State");waitJob(ui);const auto& actual=s.currentState.analysis->detail;const auto& measured=selected->eventAssessment()->measured();
    expect(actual.integrated.value==measured.integrated.value&&actual.lra.value==measured.lra.value&&actual.channelTruePeaks==measured.channelTruePeaks&&actual.channelSamplePeaks==measured.channelSamplePeaks&&actual.truePeakEighthFrames==measured.truePeakEighthFrames&&actual.samplePeakFrames==measured.samplePeakFrames,"independent competing Current matches exact candidate metrics and all peak positions");
    const auto savedHistory=encodeAppliedExplanation(*s.lastAppliedExplanation);s.discardSectionProposal();Access::closeReview(ui);Access::refresh(ui);Access::save(ui,true);const auto finalFile=root.getChildFile("native-applied.levelscope");finishSaveAs(ui,finalFile);Access::open(ui,finalFile);waitFile(ui);expect(s.lastAppliedExplanation&&encodeAppliedExplanation(*s.lastAppliedExplanation)==savedHistory&&!s.sectionProposal&&s.currentState.hasMeasurement,"native applied history/current save-reopen preserves verification without transient authority");
    s=fixture;s.auditionBusy=busy;auto requirements=s.jobRequirements;requirements.maximumLraTenths=49;auto limits=s.proposalPolicy;limits.durationThousandths=600;
    expect(s.setJobRequirements(requirements,errorForInner)&&s.setProposalPolicy(limits,errorForInner)&&s.setRegionPlanning(RegionPlanningSettings::innerDefaults(RegionMode::automatic,sectionInputIdentity(s.measurementContext(MeasurementDomain::source))),errorForInner),"owned short competing fixture for cancellation and late completion");Access::refresh(ui);
    Access::propose(ui);const auto proposalCancelAt=juce::Time::getMillisecondCounterHiRes();Access::cancel(ui);waitJob(ui);expect(!s.sectionProposal&&s.job.state==JobState::cancelled&&juce::Time::getMillisecondCounterHiRes()-proposalCancelAt<=1000,"competing proposal cancellation stays prompt and grants no authority");
    Access::propose(ui);s.newDocument();const auto replacement=s.sessionId;const auto deadline=juce::Time::getMillisecondCounterHiRes()+90000;
    while(Access::proposalPending(ui)&&juce::Time::getMillisecondCounterHiRes()<deadline){Access::tick(ui);juce::MessageManager::getInstance()->runDispatchLoopUntil(5);}
    expect(!Access::proposalPending(ui)&&s.sessionId==replacement&&!s.sectionProposal&&!s.sourceDocument.hasSource&&!s.currentState.hasMeasurement,"late competing completion after New cannot install old result or measurement");
    expect(juce::SHA256(supplied).toHexString()==originalHash&&juce::SHA256(originalAudio).toHexString()==audioHash,"supplied original session and audio remain byte-identical");expect(Access::stopped(ui),"all competing native checks remain stopped and disconnected");
}

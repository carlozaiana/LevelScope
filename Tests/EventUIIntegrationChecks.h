// Included in RegionUIIntegrationTests' fixture namespace. These tests exercise
// the actual main-window selection callbacks with synthetic PCM and no device.
void eventEditorChecks(StandaloneSessionModel& s,const juce::File& root)
{
    const auto baseline=s.currentState.analysis;const auto revision=s.revision;
    const auto policy=s.proposalPolicy.digest();RegionPlanningPanel panel(s);
    expect(!s.regionPlanning&&control<juce::ComboBox>(panel,"Proposal method").getSelectedId()==1,"event treatment is not a new-session default");
    expect(s.proposalPolicy.durationThousandths==200&&s.proposalPolicy.additionalTenths==30&&s.proposalPolicy.cumulativeTenths==120&&s.proposalPolicy.reserveTenths==5,"conservative duration/gain/reserve defaults unchanged");
    std::optional<RegionPlanningSettings> draft;choose(panel,"Proposal method",2);
    expect(panel.readDraft(draft).isEmpty()&&draft->schema==1&&!draft->allowAutomaticMargins,"choosing M20 alone preserves legacy schema1");
    choose(panel,"Proposal method",3);
    expect(panel.readDraft(draft).isEmpty()&&draft->schema==2&&draft->preferCompleteEvent&&draft->compareAlternatives&&!draft->allowAutomaticMargins&&draft->leftMarginMs==0&&draft->rightMarginMs==0,"explicit event draft requires comparison without granting outward margins");
    expect(!control<juce::ToggleButton>(panel,"Assess alternatives after a passing candidate").isEnabled()&&!control<juce::TextEditor>(panel,"Left margin seconds (0-3, step 0.1)").isEnabled(),"comparison cannot be disabled and unpermitted automatic margins cannot be edited");
    expect(!s.regionPlanning&&s.currentState.analysis==baseline&&s.revision==revision,"private event opt-in changes no active settings or audio");
    toggle(panel,"Allow automatic outward margins",true);text(panel,"Left margin seconds (0-3, step 0.1)","1.0");text(panel,"Right margin seconds (0-3, step 0.1)","1.0");
    expect(panel.readDraft(draft).isEmpty()&&draft->allowAutomaticMargins&&draft->leftMarginMs==1000&&draft->rightMarginMs==1000,"one-second experimental margins require an independent explicit permission");
    toggle(panel,"Allow automatic outward margins",false);expect(panel.readDraft(draft).isEmpty()&&draft->leftMarginMs==0&&draft->rightMarginMs==0,"revoking permission clears outward margin draft");
    click(panel,"Confirm settings");expect(s.regionPlanning&&s.regionPlanning->schema==2&&s.proposalPolicy.digest()==policy&&s.currentState.analysis==baseline,"event confirmation keeps 20 percent policy and baseline");
    auto loaded=*s.regionPlanning;loaded.eventFeatures.hopMs=5;juce::String error;expect(s.setRegionPlanning(loaded,error),"load supported nondefault full feature definition");
    click(panel,"Discard edits / reload");expect(panel.readDraft(draft).isEmpty()&&draft->digest()==loaded.digest(),"opening/reconfirming settings preserves exact loaded feature definition");
    imageOf(panel,root.getChildFile("event-settings.png"));
    choose(panel,"Region mode",3);text(panel,"Temporary begin frame","1449600");text(panel,"Temporary end frame","1728000");
    expect(panel.readDraft(draft).isEmpty()&&!draft->allowGlobal&&!draft->allowAutomaticMargins&&draft->selection.begin==1449600,"selected passage uses exact frames with global permission reset");
    choose(panel,"Region mode",1);choose(panel,"Proposal method",2);click(panel,"Confirm settings");
    expect(s.regionPlanning->schema==1&&!s.regionPlanning->allowAutomaticMargins&&s.regionPlanning->leftMarginMs==0,"explicit preference removal returns to legacy semantics");
    choose(panel,"Proposal method",1);click(panel,"Confirm settings");expect(!s.regionPlanning&&s.currentState.analysis==baseline,"M19 remains available without changing measurement");
}

void eventUIChecks(const juce::File& root)
{
    auto fixture=programme(root.getChildFile("fixture"),2);eventEditorChecks(fixture,root);
    StandaloneMainComponent ui(emptyOutputCatalog(),noOutputSettings());Access::stopTimer(ui);ui.setSize(1320,1100);ui.addToDesktop(juce::ComponentPeer::windowHasTitleBar);ui.setVisible(true);
    auto& s=Access::model(ui);const auto busy=s.auditionBusy;s=fixture;s.auditionBusy=busy;
    auto requirements=s.jobRequirements;requirements.maximumLraTenths=49;juce::String error;expect(s.setJobRequirements(requirements,error),"explicit controlled 4.9 LU requirement");
    auto settings=RegionPlanningSettings::eventDefaults(RegionMode::selection,sectionInputIdentity(s.measurementContext(MeasurementDomain::source)));
    settings.selection={48000*302/10,48000*36};settings.allowGlobal=true;
    expect(s.setRegionPlanning(settings,error),"explicit late-start selection for partial coverage");auto policy=s.proposalPolicy;policy.durationThousandths=600;expect(s.setProposalPolicy(policy,error),"explicit test-only 60 percent allowance");
    Access::refresh(ui);const auto baseline=s.currentState.analysis,source=s.source.analysis;const auto revision=s.revision;const auto domain=s.auditionDomain;
    const auto processing=proposalPlanIdentity(s.currentStateDocument.processingPlan());Heartbeat heartbeat;heartbeat.begin();Access::propose(ui);
    expect(s.job.isActive(),"native Propose starts async event assessment after request capture");waitJob(ui);
    expect(s.sectionProposal&&s.sectionProposal->status()==SectionProposalStatus::reviewRequired,"partial-only native proposal prominently requires review");
    const auto unselected=s.sectionProposal;const auto review=unselected->regionReview();expect(review&&review->eventMode&&!review->selectedAttempt&&!review->recommendedAttempt,"review required has neither recommendation nor selected authority");
    expect(!s.canApplySectionProposal()&&!control<juce::TextButton>(Access::proposal(ui),"Apply proposal").isEnabled(),"passing alternatives do not implicitly authorize Apply");
    const auto pass=std::find_if(review->attempts.begin(),review->attempts.end(),[](const auto& a){return a.evaluation.verdict==Verdict::pass;});expect(pass!=review->attempts.end(),"owned actual passing alternative exists");const auto passing=size_t(pass-review->attempts.begin());
    expect(Access::eventAlternativesEnabled(ui),"final candidate exposes the action after the separately preserved authority gate passed");
    Access::enableEventAlternatives(ui,false);Access::review(ui);auto* panel=Access::reviewPanel(ui);expect(panel!=nullptr,"main Review opens native event panel");choose(*panel,"Measured region attempt",int(passing+2));
    expect(!control<juce::TextButton>(*panel,"Use this measured alternative").isEnabled(),"product qualification gate blocks fallback before authority checks");
    Access::closeReview(ui);Access::enableEventAlternatives(ui,true);Access::review(ui);panel=Access::reviewPanel(ui);expect(panel!=nullptr,"explicit integration test gate exercises actual product callbacks");
    const auto view=[&]() -> ProposalInspectionPanel& {for(auto* c:panel->getChildren())if(auto* p=dynamic_cast<ProposalInspectionPanel*>(c);p&&p->isVisible())return *p;throw std::runtime_error("Missing current attempt report");};
    for(size_t n=0;n<review->attempts.size();++n)
    {
        choose(*panel,"Measured region attempt",int(n+2));Access::refresh(ui);const auto& a=review->attempts[n];
        expect(s.sectionProposal==unselected&&!review->selectedAttempt&&Access::lane(ui).diagnosticRegionInspection(),"browsing each real attempt remains inspection only");
        expect(control<juce::TextButton>(*panel,"Use this measured alternative").isEnabled()==(a.evaluation.verdict==Verdict::pass&&bool(a.receipt)),"only owned measured passes expose selection");
        expect(view().allText().contains(a.coverage.text())&&view().allText().contains("ACTUAL measured rule result")&&view().allText().contains("Global gain"),"complete attempt report includes coverage, actual rules and exact gains");
        expect(!control<juce::TextButton>(Access::proposal(ui),"Apply proposal").isEnabled(),"inspected passing or failed attempts cannot bypass explicit selection");
    }
    choose(*panel,"Measured region attempt",int(passing+2));
    expect(!s.selectMeasuredRegionAlternative("foreign-result",pass->attemptIdentity)&&!s.selectMeasuredRegionAlternative(review->resultIdentity,"invented-attempt"),"foreign result and fabricated attempt cannot mint authority");
    // Mutate only this disposable model, preserving the immutable result. Every
    // blocked control is also invoked directly to exercise its click-time check.
    const auto restore=s;
    const auto rejected=[&](const char* reason)
    {
        panel->refreshSelection();auto& button=control<juce::TextButton>(*panel,"Use this measured alternative");expect(!button.isEnabled(),reason);button.onClick();
        expect(s.sectionProposal==unselected&&s.revision==revision,"blocked selection preserves result, revision and audio");s=restore;Access::refresh(ui);panel->refreshSelection();
    };
    s.auditionBusy=[]{return true;};rejected("selection checks stopped transport");
    s.persistenceLoading=true;rejected("selection blocked during session loading");
    s.sectionTrimDraft.emplace();rejected("pending manual gain draft blocks alternative selection");
    s.boundaryBatch=BoundaryEditBatch{s.manualState(),s.manualState(),s.sessionId,s.revision};rejected("pending boundary batch blocks alternative selection");
    expect(s.job.begin(JobStage::analysis,s.revision),"simulate active independent job");rejected("active job blocks selection");
    auto changed=requirements;changed.maximumLraTenths=48;expect(s.setJobRequirements(changed,error),"change real target requirement");rejected("changed requirements make selection stale");
    auto changedSettings=settings;changedSettings.eventFeatures.hopMs=5;expect(s.setRegionPlanning(changedSettings,error),"change real feature settings");rejected("changed feature definition makes selection stale");
    auto changedPolicy=policy;changedPolicy.durationThousandths=200;expect(s.setProposalPolicy(changedPolicy,error),"change real proposal allowance");rejected("changed policy makes selection stale");
    s.sourceContentInvalid=true;rejected("invalidated Source blocks selection");
    s.sessionId="replacement-session";rejected("old review window cannot select into New/Open runtime");
    s.sectionProposal.reset();panel->refreshSelection();expect(!control<juce::TextButton>(*panel,"Use this measured alternative").isEnabled(),"discarded result disables old review action");s=restore;Access::refresh(ui);panel->refreshSelection();
    // A selection may be made from this review while a different timeline
    // attempt is displayed. That display must survive the immutable replacement.
    const size_t other=review->attempts.size()>1?(passing+1)%review->attempts.size():passing;
    Access::lane(ui).inspectRegionAttempt(other);Access::refresh(ui);click(*panel,"Use this measured alternative");
    const auto selected=s.sectionProposal;expect(selected!=unselected&&selected->explicitlySelectedAlternative()&&selected->status()==SectionProposalStatus::ready,"native action creates immutable explicitly selected result");
    expect(selected->regionReview()->resultIdentity==review->resultIdentity&&selected->regionReview()->selectedAttempt==passing&&!selected->regionReview()->recommendedAttempt,"owned UUID preserved; chosen partial pass does not become recommendation");
    expect(selected->outcomeTitle().contains(eventCoverageName(pass->coverage.category))&&selected->outcomeTitle().contains("Selected passing alternative"),"selected fallback retains exact coverage class in outcome");
    expect(s.revision==revision&&s.currentState.analysis==baseline&&s.source.analysis==source&&proposalPlanIdentity(s.currentStateDocument.processingPlan())==processing&&s.auditionDomain==domain&&Access::stopped(ui),"selection leaves samples, baseline, audition domain and stopped/disconnected output unchanged");
    if(other!=passing)expect(Access::lane(ui).diagnosticRegionInspection()&&!control<juce::TextButton>(Access::proposal(ui),"Apply proposal").isEnabled(),"immutable selection replacement retains different diagnostic inspection and blocks Apply");
    click(*panel,"Inspect in timeline");Access::refresh(ui);expect(!Access::lane(ui).diagnosticRegionInspection()&&control<juce::TextButton>(Access::proposal(ui),"Apply proposal").isEnabled(),"inspecting selected measured attempt restores explicit Apply");
    bool freshStatus=false;for(auto* child:view().getChildren())if(auto* label=dynamic_cast<juce::Label*>(child))freshStatus|=label->getText().contains("Selected passing alternative")&&label->getText().contains("Explicit Apply is in the main window");
    expect(freshStatus,"selection and inspection refresh the report's status immediately without stale Review required text");
    imageOf(*panel,root.getChildFile("selected-partial-review.png"));imageOf(ui,root.getChildFile("selected-partial-main.png"));
    const auto summaryView=[&]() -> juce::TextEditor* {for(auto* child:panel->getChildren())if(auto* e=dynamic_cast<juce::TextEditor*>(child);e&&e->isVisible()&&e->getName()=="Measured alternative coverage and gains")return e;return nullptr;};
    const auto* selectedSummary=summaryView();
    std::vector<double> switches;for(int k=0;k<40;++k){const auto start=juce::Time::getMillisecondCounterHiRes();choose(*panel,"Measured region attempt",k%2?int(passing+2):1);switches.push_back(juce::Time::getMillisecondCounterHiRes()-start);juce::MessageManager::getInstance()->runDispatchLoopUntil(20);}
    std::sort(switches.begin(),switches.end());expect(switches[37]<=100,"native complete/actual review switching p95 <=100 ms");auto timing=juce::var(new juce::DynamicObject());timing.getDynamicObject()->setProperty("event_review_switch_p95_ms",switches[37]);cases.add(timing);heartbeat.finish();
    expect(selectedSummary&&summaryView()==selectedSummary,"switching reuses the selected attempt's shaped coverage summary");
    choose(*panel,"Measured region attempt",1);const auto* completeSummary=summaryView();const auto completeText=completeSummary?completeSummary->getText():juce::String();
    choose(*panel,"Measured region attempt",int(passing+2));choose(*panel,"Measured region attempt",1);
    expect(completeSummary&&completeSummary!=selectedSummary&&summaryView()==completeSummary&&summaryView()->getText()==completeText,"complete explanation retains its separate summary without dropping status text");
    choose(*panel,"Measured region attempt",int(passing+2));click(*panel,"Inspect in timeline");Access::refresh(ui);
    s.auditionBusy=[]{return true;};Access::refresh(ui);expect(!control<juce::TextButton>(Access::proposal(ui),"Apply proposal").isEnabled()&&!s.applySectionProposal(),"selected actual pass still cannot Apply during transport activity");s.auditionBusy=busy;Access::refresh(ui);
    click(Access::proposal(ui),"Apply proposal");expect(s.currentStateNeedsRemeasurement()&&!s.currentState.hasMeasurement&&s.lastAppliedExplanation,"explicit Apply makes Current unverified and captures history");
    const auto applied=proposalPlanIdentity(s.currentStateDocument.processingPlan());const auto history=encodeAppliedExplanation(*s.lastAppliedExplanation);
    click(ui,"Re-measure Current State");Access::cancel(ui);const auto cancelledAt=juce::Time::getMillisecondCounterHiRes();waitJob(ui);
    expect(!s.currentState.hasMeasurement&&s.currentStateNeedsRemeasurement()&&proposalPlanIdentity(s.currentStateDocument.processingPlan())==applied,"cancelled independent measurement retains applied gains and honest unverified Current");
    const auto cancelMs=juce::Time::getMillisecondCounterHiRes()-cancelledAt;expect(cancelMs<=1000,"native remeasurement cancellation <=1 second");auto cancel=juce::var(new juce::DynamicObject());cancel.getDynamicObject()->setProperty("event_remeasurement_cancel_ms",cancelMs);cases.add(cancel);
    click(ui,"Re-measure Current State");waitJob(ui);expect(s.currentState.hasMeasurement&&!s.currentStateNeedsRemeasurement(),"independent native remeasurement verifies Current");
    const auto& actual=s.currentState.analysis->detail;const auto& measured=selected->eventAssessment()->measured();
    expect(actual.integrated.value==measured.integrated.value&&actual.lra.value==measured.lra.value&&actual.channelTruePeaks==measured.channelTruePeaks&&actual.channelSamplePeaks==measured.channelSamplePeaks&&actual.truePeakEighthFrames==measured.truePeakEighthFrames&&actual.samplePeakFrames==measured.samplePeakFrames,"native independent Current matches actual candidate scalar and every channel peak/location exactly");
    const auto verified=s.currentState.analysis;Access::propose(ui);waitJob(ui);
    expect(s.sectionProposal&&s.sectionProposal->status()==SectionProposalStatus::noChange&&!control<juce::TextButton>(Access::proposal(ui),"Apply proposal").isEnabled()&&s.currentState.analysis==verified,"already verified pass shows no-change without Apply or another measurement");
    JobControl job;const auto file=root.getChildFile("event-selected.levelscope");expect(sf::save(sf::capture(s),file,job).succeeded(),"new event UI history saved");auto read=sf::read(file,job);expect(read.document&&read.document->loadedSchema==10,"event UI writes schema10");auto opened=sf::verify(*read.document,job);expect(opened.prepared&&opened.prepared->model.lastAppliedExplanation&&!opened.prepared->model.sectionProposal,"reopen restores complete history without runtime selection authority");
    Access::closeReview(ui);s=opened.prepared->model;s.auditionBusy=busy;Access::refresh(ui);expect(encodeAppliedExplanation(*s.lastAppliedExplanation)==history&&s.regionPlanning->digest()==settings.digest(),"reopened explanation and full settings roundtrip exactly");
    Access::review(ui,true);expect(!Access::reviewPanel(ui)&&!control<juce::TextButton>(Access::proposal(ui),"Apply proposal").isEnabled(),"main historical inspection has no live alternative-selection or Apply capability");
    s.setCurrentGainDb(s.currentStateDocument.staticGainDb+.1);Access::refresh(ui);expect(Access::historyStatus(ui).containsIgnoreCase("historical")&&s.currentStateNeedsRemeasurement()&&encodeAppliedExplanation(*s.lastAppliedExplanation)==history,"later processing marks history historical without losing evidence or reviving verification");Access::closeReview(ui);
    // New invocation identity and late completion after New are actual worker
    // paths. No receipt is synthesized to stand in for measured authority.
    s=fixture;s.auditionBusy=busy;expect(s.setJobRequirements(requirements,error)&&s.setRegionPlanning(settings,error)&&s.setProposalPolicy(policy,error),"restore owned fixture for cancellation and runtime guards");Access::refresh(ui);
    Access::propose(ui);const auto proposalCancelAt=juce::Time::getMillisecondCounterHiRes();Access::cancel(ui);waitJob(ui);expect(!s.sectionProposal&&s.job.state==JobState::cancelled,"cancelled native event request creates no Apply authority");
    const auto proposalCancelMs=juce::Time::getMillisecondCounterHiRes()-proposalCancelAt;expect(proposalCancelMs<=1000,"native event proposal cancellation <=1 second");auto cancellation=juce::var(new juce::DynamicObject());cancellation.getDynamicObject()->setProperty("event_proposal_cancel_ms",proposalCancelMs);cases.add(cancellation);
    Access::propose(ui);s.newDocument();const auto newRuntime=s.sessionId;const auto until=juce::Time::getMillisecondCounterHiRes()+90000;
    while(Access::proposalPending(ui)&&juce::Time::getMillisecondCounterHiRes()<until){Access::tick(ui);juce::MessageManager::getInstance()->runDispatchLoopUntil(5);}
    expect(!Access::proposalPending(ui)&&s.sessionId==newRuntime&&!s.sectionProposal&&!s.currentState.hasMeasurement&&!s.sourceDocument.hasSource,"late native completion after New cannot populate replacement session or restore selection");
    s=fixture;s.auditionBusy=busy;auto automatic=RegionPlanningSettings::eventDefaults(RegionMode::automatic,settings.inputIdentity);automatic.allowAutomaticMargins=true;automatic.leftMarginMs=automatic.rightMarginMs=1000;
    expect(s.setRegionPlanning(automatic,error)&&s.setProposalPolicy(policy,error)&&s.setJobRequirements(requirements,error),"explicit complete-event automatic experiment");Access::refresh(ui);Access::propose(ui);waitJob(ui);
    expect(s.sectionProposal&&s.sectionProposal->status()==SectionProposalStatus::ready&&s.sectionProposal->regionReview()->recommendedAttempt&&s.sectionProposal->regionReview()->selectedAttempt==s.sectionProposal->regionReview()->recommendedAttempt,"supported-complete actual candidate is both recommended and initially selected");
    expect(!s.sectionProposal->explicitlySelectedAlternative()&&!Access::lane(ui).diagnosticRegionInspection()&&control<juce::TextButton>(Access::proposal(ui),"Apply proposal").isEnabled(),"supported-complete recommendation retains ordinary explicit Apply");
    expect(s.sectionProposal->regionReview()->resultIdentity!=review->resultIdentity,"new proposal invocation has a distinct immutable result UUID");
    Access::review(ui);panel=Access::reviewPanel(ui);choose(*panel,"Measured region attempt",int(*s.sectionProposal->regionReview()->selectedAttempt+2));expect(!control<juce::TextButton>(*panel,"Use this measured alternative").isEnabled(),"already selected complete recommendation does not offer redundant reselection");imageOf(*panel,root.getChildFile("complete-event-review.png"));Access::closeReview(ui);
    auto globalRequirements=requirements;globalRequirements.maximumLraTenths=200;expect(s.setJobRequirements(globalRequirements,error),"global-only controlled requirement");Access::refresh(ui);Access::propose(ui);waitJob(ui);
    expect(s.sectionProposal&&s.sectionProposal->status()==SectionProposalStatus::ready&&s.sectionProposal->regionReview()->selectedAttempt,"actual global-only event result ready");
    const auto globalIndex=*s.sectionProposal->regionReview()->selectedAttempt;expect(s.sectionProposal->regionReview()->attempts[globalIndex].family=="measured global-only"&&regionAttemptAuthority(*s.sectionProposal,globalIndex).contains("no local event treatment"),"actual global-only result is explicitly labelled without invented local coverage");
    auto impossible=requirements;impossible.maximumLraTenths=1;expect(s.setJobRequirements(impossible,error),"deliberately unreachable requirement for native failure presentation");Access::refresh(ui);Access::propose(ui);waitJob(ui);
    expect(s.sectionProposal&&s.sectionProposal->status()==SectionProposalStatus::noSolution&&!control<juce::TextButton>(Access::proposal(ui),"Apply proposal").isEnabled(),"bounded measured failure remains visible with no Apply");
    Access::review(ui);panel=Access::reviewPanel(ui);expect(panel!=nullptr,"failed actual attempts remain reviewable");
    for(size_t n=0;n<s.sectionProposal->regionReview()->attempts.size();++n){choose(*panel,"Measured region attempt",int(n+2));expect(!control<juce::TextButton>(*panel,"Use this measured alternative").isEnabled(),"failed native actual attempt cannot be selected");}
    imageOf(*panel,root.getChildFile("failed-event-review.png"));Access::closeReview(ui);
    // Exact class labels are also checked for documentary-only states. They
    // cannot manufacture a private receipt or a passing selection action.
    for(auto category:{EventCoverageClass::supportedComplete,EventCoverageClass::knownPartial,EventCoverageClass::uncertain,EventCoverageClass::partialAndUncertain})
    {RegionAttempt a;a.coverage.category=category;expect(eventAttemptCoverage(a)==eventCoverageName(category),"coverage class labels preserve supported, partial, uncertain and mixed distinctions");}
    RegionAttempt global;global.family="measured global-only";expect(eventAttemptCoverage(global).contains("no local event treatment"),"global-only label never invents local event coverage");
    expect(Access::stopped(ui),"entire integrated event UI pass leaves output disconnected and transport stopped");
}

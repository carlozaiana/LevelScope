// Main-window orchestration with real asynchronous measurements, no device IO.
template<class Predicate> void waitOperation(StandaloneMainComponent& ui,Predicate done)
{
    const auto deadline=juce::Time::getMillisecondCounterHiRes()+120000;
    while(!done()&&juce::Time::getMillisecondCounterHiRes()<deadline)
    {Access::tick(ui);juce::MessageManager::getInstance()->runDispatchLoopUntil(5);}
    expect(done(),"workflow reaches the expected asynchronous phase within timeout");Access::refresh(ui);
}
void completeOperation(StandaloneMainComponent& ui)
{waitOperation(ui,[&]{return !Access::operation(ui).active()&&!Access::workers(ui);});}
void installWorkflowFixture(StandaloneMainComponent& ui,const StandaloneSessionModel& fixture)
{auto busy=Access::model(ui).auditionBusy;Access::model(ui)=fixture;Access::model(ui).auditionBusy=busy;Access::operation(ui)=ProposalOperation();Access::refresh(ui);}
void proposalOperationChecks(const juce::File& root,const juce::File& supplied)
{
    using Phase=ProposalOperation::Phase;
    StandaloneMainComponent ui(emptyOutputCatalog(),noOutputSettings());Access::stopTimer(ui);ui.setSize(1360,1100);
    ui.addToDesktop(juce::ComponentPeer::windowHasTitleBar);ui.setVisible(true);
    auto& s=Access::model(ui);auto& op=Access::operation(ui);juce::String error;JobControl job;
    std::vector<RegionReviewScore> shortlist;
    for(double excess:{.020,.014,.007})
    {RegionReviewScore score;score.complete=true;score.lraExcess=excess;score.failedRules=1;score.activity.affectedFrames=int64_t(shortlist.size()+1);score.id=juce::String(shortlist.size());shortlist.push_back(score);
        expect(selectRegionReview(shortlist)==std::optional<size_t>(shortlist.size()==3?1:0),"anchored failure shortlist can return an earlier non-winning attempt after the third failure");}
    const auto stereo=programme(root.getChildFile("stereo"),2);
    const auto setScalar=[&]
    {
        auto req=s.jobRequirements;req.maximumLraTenths.reset();expect(s.setJobRequirements(req,error),"scalar requirements accepted");
        auto policy=s.proposalPolicy;policy.durationThousandths=1000;expect(s.setProposalPolicy(policy,error),"global duration allowance explicit");
    };
    for(int channels:{2,1,6})
    {
        installWorkflowFixture(ui,channels==2?stereo:programme(root.getChildFile(channels==1?"mono":"surround"),channels));setScalar();
        auto manual=s.manualState();manual.globalDb=-1;
        processing::SectionGainPlan local;local.inputIdentity=manual.map.inputIdentity;local.frames=manual.map.frames;local.sampleRate=manual.map.sampleRate;
        local.intervals={{0,480000,0},{480000,960000,-1},{960000,local.frames,0}};manual.gains=local;
        expect(s.commitManualEdit(manual,"manual baseline"),"manual gains precede one-click Propose");const auto provenance=s.currentStateDocument.sectionGainAnnotations.digest();
        Access::refresh(ui);click(Access::proposal(ui),"Propose");
        expect(op.phase==Phase::current&&s.manualState()==manual,"Propose visibly remeasures actual manual Current without resetting gains");
        expect(!ui.performCommand(StandaloneMainComponent::newDocument)&&!ui.performCommand(StandaloneMainComponent::measureCurrent),"New and competing measurement blocked throughout operation");
        waitOperation(ui,[&]{return op.phase==Phase::apply||!op.active();});
        expect(op.phase==Phase::apply&&op.pending&&op.pending->assessment()&&op.pending->assessment()->domain==MeasurementDomain::privateCandidate,"forecast-only global plan gets a complete private measurement before Apply");
        expect(s.manualState()==manual&&s.currentState.hasMeasurement&&!op.applied,"private candidate does not alter live processing or Current measurement");
        const auto candidate=op.pending->assessment()->identity;
        completeOperation(ui);
        expect(op.phase==Phase::verified&&s.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass,"main Propose independently verifies an applied measured pass");
        expect(s.currentState.analysis->domain==MeasurementDomain::currentState&&s.currentState.analysis->identity!=candidate,"Current PASS comes from distinct accepted Current-domain measurement");
        expect(s.currentStateDocument.sectionGains->audioKey()==manual.gains->audioKey()&&s.currentStateDocument.sectionGainAnnotations.digest()==provenance,"scalar application preserves local envelope and provenance");
        expect(s.lastAppliedExplanation&&s.lastAppliedExplanation->before.staticGainDb==-1&&s.lastAppliedExplanation->candidate&&s.lastAppliedExplanation->narrative.contains("selected by Propose"),"complete history identifies actual baseline, private measurement and automatic selection");
        expect(op.canRecover(s)&&Access::stopped(ui)&&!s.sectionProposal,"recovery available, transient Apply authority cleared, playback never started");
        const auto applied=s.manualState();const auto history=encodeAppliedExplanation(*s.lastAppliedExplanation);
        expect(s.undoManual()&&s.manualState()==manual&&!s.currentState.hasMeasurement&&!s.lastAppliedExplanation,"Undo restores pre-Propose processing and prior history without restoring PASS");
        expect(s.undoManual(true)&&s.manualState()==applied&&encodeAppliedExplanation(*s.lastAppliedExplanation)==history&&!s.currentState.hasMeasurement,"Redo restores processing and complete explanation only");
        op.observe(s);expect(op.phase==Phase::stale&&!op.canChoose(s)&&!op.canRecover(s),"manual Undo/Redo expires replacement authority");
        const auto file=root.getChildFile("scalar-"+juce::String(channels)+".levelscope");expect(sf::save(sf::capture(s),file,job).succeeded(),"save scalar applied history");
        auto read=sf::read(file,job);auto verified=sf::verify(*read.document,job);expect(verified.prepared&&encodeAppliedExplanation(*verified.prepared->model.lastAppliedExplanation)==history&&!verified.prepared->model.currentState.hasMeasurement,"save/reopen retains history separately from unverified Current");
        Access::workflow(ui);completeOperation(ui);expect(op.phase==Phase::noChange&&s.currentState.hasMeasurement,"already passing applied audio is prepared and finishes without redundant correction");
        auto lra=s.jobRequirements;lra.maximumLraTenths=200;expect(s.setJobRequirements(lra,error),"Standard legacy LRA requirement retained for capability check");
        expect(s.setCurrentGainDb(s.currentStateDocument.staticGainDb-1),"prepare unchanged-method scalar correction under LRA ceiling");Access::refresh(ui);Access::workflow(ui);completeOperation(ui);
        expect(op.phase==Phase::verified&&op.selected->methodVersion()==lraTrimMethodVersion&&s.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass,"Standard LRA/global supported path uses one-click workflow without dispatch changes");
    }
    installWorkflowFixture(ui,stereo);setScalar();const auto original=s.manualState();
    Access::workflow(ui);Access::cancel(ui);completeOperation(ui);
    expect(op.phase==Phase::cancelled&&s.manualState()==original&&!op.applied,"cancellation before Apply preserves processing");
    Access::workflow(ui);waitOperation(ui,[&]{return op.phase==Phase::verify||!op.active();});
    expect(op.phase==Phase::verify&&op.applied&&!s.currentState.hasMeasurement,"Apply immediately invalidates Current before independent verification");
    Access::cancel(ui);completeOperation(ui);
    expect(op.phase==Phase::cancelled&&op.applied&&!s.currentState.hasMeasurement&&op.canRecover(s),"cancelled post-Apply verification retains applied-unverified state and recovery");
    const auto cancelledChoice=op.selected;const auto steps=s.manualHistory->undoCount();
    expect(ui.performCommand(StandaloneMainComponent::measureCurrent),"ordinary remeasurement is available after cancelled verification");waitJob(ui);
    expect(op.phase==Phase::verified&&op.selected==cancelledChoice&&s.currentState.hasMeasurement&&s.manualHistory->undoCount()==steps&&Access::proposal(ui).outcomeText().contains("independently verified"),"explicit Current retry reconciles workflow verdict without changing choice or Undo history");
    click(Access::proposal(ui),"Restore before Propose");expect(s.manualState()==original&&!s.currentState.hasMeasurement,"main recovery restores original processing without old measurement authority");

    // First import still adopts identical Source into untouched Current.
    installWorkflowFixture(ui,StandaloneSessionModel());expect(s.setValidatedSource(juce::File(stereo.sourceDocument.fullPath),stereo.sourceDocument.media),"fresh import model");
    s.sourceDocument.layout=stereo.sourceDocument.layout;s.setSelectedTargetProfileIndex(1);setScalar();Access::refresh(ui);Access::workflow(ui);
    expect(op.phase==Phase::initialize&&s.hasCurrentStateInitialized(),"first Propose exposes Current initialization before unmeasured Source preparation");
    Access::tick(ui);expect(op.phase==Phase::source,"initial Propose continues visible Source preparation");completeOperation(ui);
    expect(s.source.hasMeasurement&&op.baseline.current&&op.baseline.current->identicalSourceMeasurementIdentity==s.source.analysis->identity,"first Source measurement initialises identical Current inside operation");

    // Both success and cancelled/error callbacks must reject an old runtime.
    installWorkflowFixture(ui,stereo);setScalar();s.setCurrentGainDb(-3);Access::refresh(ui);Access::workflow(ui);
    const auto oldRuntime=s.sessionId;s.newDocument();s.job.begin(JobStage::analysis,s.revision);
    waitOperation(ui,[&]{return !Access::workers(ui);});
    expect(s.sessionId!=oldRuntime&&!s.source.hasMeasurement&&!s.currentState.hasMeasurement&&s.job.isActive(),"late old-runtime measurement cannot finish or mutate a new document job");
    installWorkflowFixture(ui,stereo);setScalar();Access::workflow(ui);const auto sameRuntime=s.sessionId;
    ++s.revision;completeOperation(ui);expect(s.sessionId==sameRuntime&&op.phase==Phase::stale&&s.manualState()==original&&!s.lastAppliedExplanation,"revision-changed completion is discarded before application");

    // Schema-1 immutable receipts and common-baseline alternative replacement.
    installWorkflowFixture(ui,stereo);auto req=s.jobRequirements;req.maximumLraTenths=49;expect(s.setJobRequirements(req,error),"schema1 LRA requirement");
    auto settings=RegionPlanningSettings::defaults(RegionMode::automatic,sectionInputIdentity(s.measurementContext(MeasurementDomain::source)));settings.compareAlternatives=true;
    expect(s.setRegionPlanning(settings,error),"schema1 explicit comparisons enabled");Access::refresh(ui);Access::workflow(ui);completeOperation(ui);
    expect(op.phase==Phase::verified&&op.result&&op.result->regionReview()&&op.choices.size()>1,"schema1 exposes every eligible measured alternative in main workflow");
    const auto initialResult=op.result;const auto request=op.baseline;const auto before=op.original;
    auto recomputed=generateSectionProposal(request,job);expect(recomputed.value&&(*recomputed.value)->regionReview()->resultIdentity!=initialResult->regionReview()->resultIdentity,"same-input recomputation has a distinct generation identity");
    const auto first=op.choices.front();expect(!selectWorkflowProposal(request,*recomputed.value,first.resultId,first.attemptId,false).value,"old result/attempt pair cannot select same-plan recomputation");
    const auto count=op.choices.size();
    const size_t other=(op.selectedChoice+1)%count;const auto selectedBeforeCancel=op.selectedChoice;const auto measuredBeforeCancel=s.currentState.analysis;
    choose(Access::proposal(ui),"Measured passing alternatives",int(other+1));Access::cancel(ui);completeOperation(ui);
    expect(op.selectedChoice==selectedBeforeCancel&&s.currentState.analysis==measuredBeforeCancel&&op.canChoose(s),"cancel before alternative commit preserves choice and verified Current, allowing retry");
    const auto foreign=selectWorkflowProposal(request,*recomputed.value,(*recomputed.value)->regionReview()->resultIdentity,(*recomputed.value)->regionReview()->attempts[(*recomputed.value)->regionReview()->retainedAttempt].attemptIdentity,true);
    expect(foreign.value&&!s.applyWorkflowProposal(request,initialResult,*foreign.value,op.owned,op.ownedRevision,op.requirements),"model commit rejects a measured selection owned by another same-input result generation");
    for(size_t index=0;index<count;++index)
    {
        if(index==op.selectedChoice)continue;
        choose(Access::proposal(ui),"Measured passing alternatives",int(index+1));completeOperation(ui);
        expect(op.phase==Phase::verified&&op.selected->eventAssessment()&&!op.selected->assessment(),"schema1 selected alternative exposes its owned compact receipt, never original winner snapshot");
        const auto& receipt=*op.selected->eventAssessment();
        expect(s.currentState.analysis->detail.truePeak.value==receipt.measured().truePeak.value&&s.currentState.analysis->detail.truePeakEighthFrames==receipt.measured().truePeakEighthFrames,"independent replacement measurement matches selected alternative peaks and positions exactly");
        expect(s.lastAppliedExplanation->beforeIdentity==proposalPlanIdentity(request.currentContext.processing)&&s.lastAppliedExplanation->candidate->identity==receipt.snapshotIdentity(),"replacement history retains common baseline and selected candidate evidence");
        expect(s.manualState().globalDb==op.selected->plan().staticGainDb&&s.currentStateDocument.sectionGains->audioKey()==op.selected->plan().sectionGains->audioKey(),"alternative replaces absolute plan rather than adding corrections");
    }
    const auto altHistory=encodeAppliedExplanation(*s.lastAppliedExplanation);const auto altFile=root.getChildFile("alternative.levelscope");expect(sf::save(sf::capture(s),altFile,job).succeeded(),"save alternative history");
    auto altRead=sf::read(altFile,job);auto altOpen=sf::verify(*altRead.document,job);expect(altOpen.prepared&&encodeAppliedExplanation(*altOpen.prepared->model.lastAppliedExplanation)==altHistory&&altOpen.prepared->model.currentState.hasMeasurement&&!altOpen.prepared->model.sectionProposal,"reopen restores matching local Current and historical alternative only");
    expect(op.recover(s)&&s.manualState()==before&&!s.currentState.hasMeasurement,"recovery after several alternatives restores common original baseline");
    Access::workflow(ui);completeOperation(ui);const auto oldSettings=*s.regionPlanning;
    auto changed=oldSettings;changed.needThresholdTenths=15;expect(s.setRegionPlanning(changed,error),"change settings after verification");op.observe(s);
    expect(s.setRegionPlanning(oldSettings,error),"restore original settings definition");op.observe(s);
    expect(op.phase==Phase::stale&&!op.canChoose(s)&&!op.canRecover(s),"restoring equal settings cannot resurrect expired operation authority");

    installWorkflowFixture(ui,stereo);req=s.jobRequirements;req.maximumLraTenths=0;expect(s.setJobRequirements(req,error),"deliberately restrictive requirement");
    const auto unchanged=s.manualState();Access::refresh(ui);Access::workflow(ui);completeOperation(ui);
    expect(op.phase==Phase::noPass&&!op.applied&&s.manualState()==unchanged&&op.detail.contains("do not prove"),"no eligible measured pass preserves Current and describes finite search limitation");

    // Copied development programme: advisory boundary evidence must not veto a measured pass.
    const auto hash=juce::SHA256(supplied).toHexString();auto suppliedRead=sf::read(supplied,job);expect(bool(suppliedRead.document),"supplied accepted development session read");
    auto suppliedAudio=juce::File(suppliedRead.document->source.fullPath);if(!suppliedAudio.existsAsFile())suppliedAudio=supplied.getSiblingFile("supplied-audio.wav");
    auto audioCopy=root.getChildFile("supplied-audio.wav");expect(suppliedAudio.copyFileTo(audioCopy)&&sf::relinkAudio(*suppliedRead.document,audioCopy),"operate on copied supplied audio only");
    auto loaded=sf::verify(*suppliedRead.document,job);expect(bool(loaded.prepared),"copied development session verifies");installWorkflowFixture(ui,loaded.prepared->model);
    auto neutral=s.manualState();neutral.gains.reset();neutral.globalDb=0;expect(s.commitManualEdit(neutral,"owned neutral test baseline"),"neutral processing for supplied reassessment");
    Access::refresh(ui);Heartbeat heartbeat;heartbeat.begin();Access::workflow(ui);completeOperation(ui);heartbeat.finish();
    expect(op.phase==Phase::verified&&op.advisory&&s.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass,"measured pass with uncertain boundaries automatically applies and independently verifies");
    expect(Access::proposal(ui).outcomeText().contains("check region boundaries")&&op.detail.contains("s:"),"persistent prominent boundary advisory identifies affected intervals");
    const auto selected=op.selectedChoice;Access::review(ui);auto* panel=Access::reviewPanel(ui);expect(panel!=nullptr,"optional complete search inspection available after Apply");
    for(size_t n=0;n<op.result->regionReview()->attempts.size();++n)choose(*panel,"Measured region attempt",int(n+2));
    expect(op.selectedChoice==selected&&!s.sectionProposal&&Access::stopped(ui),"browsing attempts never switches applied alternative or starts playback");Access::closeReview(ui);
    click(Access::proposal(ui),"Inspect region boundaries");expect(Access::inspection(ui)==2&&!Access::lane(ui).selectsSections(),"main advisory opens applied Current region handles");
    for(int n=0;n<20;++n){Access::tick(ui);juce::MessageManager::getInstance()->runDispatchLoopUntil(10);}
    imageOf(ui,root.getChildFile("verified-boundary-advisory.png"));
    const auto appliedHistory=encodeAppliedExplanation(*s.lastAppliedExplanation);const auto cut=s.currentStateDocument.sectionGains->intervals[1].begin;
    expect(s.beginAppliedBoundary(1)&&s.previewAppliedBoundary(s.appliedBoundaryDrag->token,cut+20)&&s.commitAppliedBoundary(s.appliedBoundaryDrag->token),"user can refine applied boundary after successful workflow");Access::tick(ui);
    expect(op.phase==Phase::stale&&!op.canChoose(s)&&!s.currentState.hasMeasurement&&encodeAppliedExplanation(*s.lastAppliedExplanation)==appliedHistory,"manual boundary adjustment expires alternatives and PASS while preserving historical explanation");
    expect(juce::SHA256(supplied).toHexString()==hash&&Access::stopped(ui),"original supplied session preserved and no physical device opened");ui.removeFromDesktop();
}

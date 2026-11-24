// Focused native step-one workflow; no physical outputs or audible testing.
void unchangedSectionAudio(const StandaloneSessionModel& s,const ProcessingPlan& before)
{
    const auto after=s.currentStateDocument.processingPlan();
    juce::WavAudioFormat wav;std::unique_ptr<juce::AudioFormatReader> reader(wav.createReaderFor(s.sourceDocument.file.createInputStream().release(),true));
    expect(bool(reader),"open owned PCM for annotation-only sample comparison");
    ProgrammeGainProcessor a(before),b(after);const int channels=int(reader->numChannels);
    juce::AudioBuffer<float> dry(channels,4096),left(channels,4096),right(channels,4096);bool exact=true;
    for(int64_t first=0;first<reader->lengthInSamples;first+=4096)
    {
        const int count=int(std::min<int64_t>(4096,reader->lengthInSamples-first));
        exact&=reader->read(&dry,0,count,first,true,true);left.makeCopyOf(dry);right.makeCopyOf(dry);
        exact&=a.process(left.getArrayOfWritePointers(),channels,count,first)&&b.process(right.getArrayOfWritePointers(),channels,count,first);
        for(int ch=0;ch<channels;++ch)exact&=std::memcmp(left.getReadPointer(ch),right.getReadPointer(ch),size_t(count)*sizeof(float))==0;
    }
    expect(exact,"creating sections preserves every processed sample including transition ramps");
}
void currentAndSectionsUIChecks(const juce::File& root)
{
    using Phase=ProposalOperation::Phase;
    StandaloneMainComponent ui(emptyOutputCatalog(),noOutputSettings());Access::stopTimer(ui);ui.setSize(1360,1100);
    ui.addToDesktop(juce::ComponentPeer::windowHasTitleBar);ui.setVisible(true);
    auto& s=Access::model(ui);auto& op=Access::operation(ui);auto& lane=Access::lane(ui);auto& timeline=Access::timeline(ui);juce::String error;
    const auto ready=[&]{return control<juce::TextButton>(Access::proposal(ui),"Propose").isEnabled();};
    Access::refresh(ui);expect(!ready()&&control<juce::Label>(Access::proposal(ui),"Propose readiness").getText().contains("Open a mix"),"empty document shows a short Propose blocker beside the controls");
    const auto scalar=[&]{auto req=s.jobRequirements;req.maximumLraTenths.reset();req.maximumMomentaryTenths.reset();req.maximumShortTermTenths.reset();expect(s.setJobRequirements(req,error),"scalar requirements retain IL and true peak");auto policy=s.proposalPolicy;policy.durationThousandths=1000;expect(s.setProposalPolicy(policy,error),"explicit global duration allowance");};
    for(int channels:{1,2,6})
    {
        const auto dir=root.getChildFile("channels-"+juce::String(channels));dir.createDirectory();
        const auto base=programme(dir.getChildFile("fixture"),channels);installWorkflowFixture(ui,base);scalar();
        const auto source=s.source.analysis;const auto copies=s.verifiedAssets->counts();
        expect(!control<juce::TextButton>(lane,"Create sections from applied regions").isEnabled(),"unity Current has no local region boundaries to copy");
        Access::clearCurrent(ui);Access::refresh(ui);
        expect(!s.hasCurrentStateInitialized()&&s.source.analysis==source&&ready(),"Clear Current retains Original and leaves Propose actionable");
        expect(control<juce::Label>(Access::proposal(ui),"Propose readiness").getText().contains("create Current"),"main controls explain that Propose initializes only absent Current");
        click(Access::proposal(ui),"Propose");
        expect(op.phase==Phase::initialize&&s.hasCurrentStateInitialized()&&s.currentStateDocument.staticGainDb==0&&!s.currentStateDocument.sectionGains,"explicit Propose initializes unity Current and exposes preparation before search");
        expect(Access::stage(ui).contains("initialized from Original")&&Access::cancellationAvailable(ui)&&Access::stopped(ui),"initialization preparation is prominent, cancellable and stopped");
        Access::cancel(ui);completeOperation(ui);
        expect(op.phase==Phase::cancelled&&!op.applied&&op.detail.contains("retained at unity")&&s.hasCurrentStateInitialized()&&s.currentStateDocument.staticGainDb==0&&!s.currentStateDocument.sectionGains,"cancellation before search retains initialized unity without applying a proposal");
        Access::clearCurrent(ui);Access::refresh(ui);
        if(channels==1)
        {
            auto req=s.jobRequirements;req.maximumLraTenths=49;expect(s.setJobRequirements(req,error),"confirmed mono uses established LRA requirements");
            auto settings=RegionPlanningSettings::innerDefaults(RegionMode::automatic,sectionInputIdentity(s.measurementContext(MeasurementDomain::source)));
            expect(s.setRegionPlanning(settings,error),"remembered competing-boundary settings accepted");
        }
        const auto method=s.regionPlanning?s.regionPlanning->digest():juce::String();
        click(Access::proposal(ui),"Propose");completeOperation(ui);
        expect((op.phase==Phase::verified||op.phase==Phase::noChange)&&s.currentState.hasMeasurement&&s.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass,"cleared Current continues through existing proposal and independent passing verification");
        expect((s.regionPlanning?s.regionPlanning->digest():juce::String())==method&&Access::stopped(ui),"initialization retains method preferences and never starts playback");
        expect(s.verifiedAssets->counts().copies==copies.copies&&s.verifiedAssets->counts().copiedBytes==copies.copiedBytes,"cleared-Current Propose reuses the session verified asset");
        imageOf(ui,dir.getChildFile("propose-from-empty.png"));

        installWorkflowFixture(ui,base);scalar();
        auto manual=s.manualState();manual.globalDb=-1;processing::SectionGainPlan local;local.inputIdentity=manual.map.inputIdentity;local.frames=manual.map.frames;local.sampleRate=manual.map.sampleRate;
        local.intervals={{0,720000,0},{720000,1440000,-2},{1440000,local.frames,0}};manual.gains=local;
        expect(s.commitManualEdit(manual,"Deliberate native gain edits"),"existing Current has deliberately unmeasured local and global edits");Access::refresh(ui);
        click(Access::proposal(ui),"Propose");expect(op.phase==Phase::current&&s.manualState()==manual&&!op.initializedCurrent,"existing edited Current is remeasured and never initialized from Source");
        Access::cancel(ui);completeOperation(ui);expect(s.manualState()==manual&&!op.applied,"cancelled preparation preserves all existing manual gains");
        click(Access::proposal(ui),"Propose");completeOperation(ui);
        expect(op.phase==Phase::verified&&s.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass&&op.original==manual,"normal Propose verifies the actual edited starting point");
        expect(s.currentStateDocument.sectionGains->audioKey()==manual.gains->audioKey(),"scalar correction keeps the existing local envelope");
        const auto verified=s.currentState.analysis;const auto history=encodeAppliedExplanation(*s.lastAppliedExplanation);
        // Exercise an annotation round trip before the UI gets an observe tick.
        expect(s.createSectionsFromAppliedRegions(error)&&s.undoManual(),"annotation creation and immediate Undo are available as one action");op.observe(s);
        expect(op.phase==Phase::stale&&!op.canRecover(s)&&!op.canChoose(s)&&s.currentState.analysis==verified,"quick annotation Undo cannot revive a previous operation or invalidate its measured audio");

        Access::inspectCurrent(ui);Access::refresh(ui);choose(lane,"Gain lane selection",2);
        const auto cut=s.currentStateDocument.sectionGains->intervals[1].begin;timeline.setViewRange(cut-48000,cut+48000);timeline.focusForTransport();
        const juce::Point<float> begin{float(timeline.graphBounds().getX()+(timeline.graphBounds().getWidth()-1)*.5),float(timeline.graphBounds().getCentreY())},end=begin.translated(15,0);
        const auto event=[&](juce::Point<float> p,bool down){const auto t=juce::Time::getCurrentTime();return juce::MouseEvent{juce::Desktop::getInstance().getMainMouseSource(),p,juce::ModifierKeys(down?juce::ModifierKeys::leftButtonModifier:0),0,0,0,0,0,&timeline,&timeline,t,p,t,1,false};};
        timeline.mouseDown(event(begin,true));timeline.mouseDrag(event(end,true));timeline.mouseUp(event(end,false));
        const auto moved=s.currentStateDocument.sectionGains->intervals[1].begin;
        expect(moved!=cut&&s.currentStateNeedsRemeasurement(),"actual graph drag changes the applied cut before section creation");
        expect(ui.performCommand(StandaloneMainComponent::measureCurrent),"remeasure the manually adjusted Current");waitJob(ui);
        expect(s.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass,"adjusted Current independently passes before annotation creation");
        auto editorial=s.manualState();auto& map=editorial.map;const auto firstId=map.sections.front().id;
        map.sections={{firstId,0,500000,"Opening"},{map.nextId++,500000,map.frames,"Body"}};
        expect(s.commitManualEdit(editorial,"Existing named annotations"),"existing named section boundaries coexist with processing regions");Access::refresh(ui);
        const auto before=s.manualState();const auto pass=s.currentState.analysis;const auto revision=s.revision;const auto count=s.manualHistory->undoCount();const auto processing=s.currentStateDocument.processingPlan();const auto asset=s.verifiedAssets->counts();
        click(lane,"Create sections from applied regions");
        expect(lane.selectsSections()&&lane.selectedCount()==1,"one-time copy exposes the new sections in the existing gain editor");
        expect(s.currentState.analysis==pass&&s.revision==revision&&s.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass,"section creation retains the exact measured Current PASS");
        expect(proposalPlanIdentity(s.currentStateDocument.processingPlan())==proposalPlanIdentity(processing)&&encodeAppliedExplanation(*s.lastAppliedExplanation)==history,"section copy preserves applied gains, ramps and complete historical explanation");
        for(const auto& old:before.map.sections)
        {auto it=std::find_if(s.sectionDocument.map().sections.begin(),s.sectionDocument.map().sections.end(),[&](const auto& n){return n.id==old.id;});expect(it!=s.sectionDocument.map().sections.end()&&it->name==old.name&&it->begin==old.begin,"existing section name and starting cut retained");}
        expect(std::any_of(s.sectionDocument.map().sections.begin(),s.sectionDocument.map().sections.end(),[&](const auto& n){return n.begin==moved;})&&s.manualHistory->undoCount()==count+1,"copied cut is the adjusted nominal position, in one Undo step");
        unchangedSectionAudio(s,processing);
        const auto created=s.manualState();
        lane.selectFrame(0,true);lane.toggleSection(moved+5000);const auto selectionCount=lane.selectedCount();
        expect(selectionCount==2,"non-adjacent sections can remain selected during repeat copy");
        click(lane,"Create sections from applied regions");
        expect(lane.selectedCount()==selectionCount&&lane.selectsSections(),"matching copy leaves the existing multi-selection intact");
        expect(s.manualState()==created&&s.manualHistory->undoCount()==count+1&&s.currentState.analysis==pass,"repeat invocation is a no-op without another Undo or measurement");
        expect(s.verifiedAssets->counts().copies==asset.copies&&s.verifiedAssets->counts().copiedBytes==asset.copiedBytes,"copying annotations cannot start source preparation or another copy");
        timeline.focusForTransport();expect(ui.performCommand(StandaloneMainComponent::undoEdit)&&s.manualState()==before&&s.currentState.analysis==pass,"native Undo removes copied annotations while preserving PASS and audio");
        expect(ui.performCommand(StandaloneMainComponent::redoEdit)&&s.manualState()==created&&s.currentState.analysis==pass,"native Redo restores annotations while preserving PASS and audio");
        JobControl save;const auto file=dir.getChildFile("sections.levelscope");expect(sf::save(sf::capture(s),file,save).succeeded(),"save created sections and measured Current");
        Access::open(ui,file);waitFile(ui);expect(s.manualState()==created&&s.currentState.hasMeasurement&&s.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass,"native reopening preserves copied sections, processing and measured PASS");
        choose(lane,"Gain lane selection",1);lane.selectFrame(moved+5000,true);
        auto& field=control<juce::TextEditor>(lane,"Timeline section gain dB");expect(field.getText().getDoubleValue()==-2,"section field displays the existing absolute local gain");
        const auto global=s.currentStateDocument.staticGainDb;const auto beforeEntry=s.manualState();
        field.setText("-",false);field.onTextChange();
        auto& copyButton=control<juce::TextButton>(lane,"Create sections from applied regions");
        expect(!copyButton.isEnabled(),"unfinished gain entry disables annotation copying");copyButton.onClick();
        expect(s.manualState()==beforeEntry,"direct section-copy callback cannot apply an unfinished gain entry");field.onEscapeKey();
        field.setText("-1",false);field.onTextChange();field.onReturnKey();
        expect(section_edit::uniform(&*s.currentStateDocument.sectionGains,moved+5000,moved+6000)==-1&&s.currentStateDocument.staticGainDb==global,"deliberate section edit replaces -2 with -1 dB rather than adding twice");
        expect(!s.currentState.hasMeasurement&&s.currentStateNeedsRemeasurement(),"later gain edit invalidates PASS immediately");
        expect(ui.performCommand(StandaloneMainComponent::measureCurrent),"explicit remeasurement verifies the deliberate section edit");waitJob(ui);
        expect(s.currentState.hasMeasurement&&!s.currentStateNeedsRemeasurement()&&Access::stopped(ui),"section edit receives its own completed measurement without playback");
        imageOf(ui,dir.getChildFile("sections-and-gain.png"));
    }
    const auto good=s;
    Access::clearCurrent(ui);JobControl cacheSave;const auto cleared=root.getChildFile("cleared.levelscope");
    expect(sf::save(sf::capture(s),cleared,cacheSave).succeeded(),"save an absent Current with verified Source cache");
    Access::open(ui,cleared);waitFile(ui);
    expect(!s.hasCurrentStateInitialized()&&s.source.analysis&&s.source.analysis->restoredFromSession,"reopen retains the explicitly absent Current and cache provenance");
    click(Access::proposal(ui),"Propose");
    expect(op.phase==Phase::initialize&&!s.currentState.hasMeasurement,"restored Source cache cannot become fresh Current measurement authority");
    Access::tick(ui);expect(op.phase==Phase::current,"reopened absent Current requires independent preparation measurement");
    Access::cancel(ui);completeOperation(ui);
    expect(op.phase==Phase::cancelled&&!s.currentState.hasMeasurement&&s.currentStateDocument.staticGainDb==0,"cancelled reopened preparation retains unverified unity Current");
    Access::clearCurrent(ui);Access::refresh(ui);click(Access::proposal(ui),"Propose");
    const auto previousRuntime=s.sessionId;s.newDocument();Access::tick(ui);
    expect(s.sessionId!=previousRuntime&&op.phase==Phase::idle&&!s.hasCurrentStateInitialized()&&!Access::workers(ui),"new runtime expires initialization before any old work can start");

    // Genuine prerequisites stay guarded; an absent Current is not itself one.
    installWorkflowFixture(ui,good);Access::clearCurrent(ui);s.sourceDocument.layout.confirmed=false;Access::refresh(ui);
    expect(!ready()&&control<juce::Label>(Access::proposal(ui),"Propose readiness").getText().contains("Confirm channel roles"),"unconfirmed roles expose a short action beside Propose");Access::workflow(ui);expect(!s.hasCurrentStateInitialized(),"blocked Propose cannot initialize Current");
    installWorkflowFixture(ui,good);Access::clearCurrent(ui);s.auditionBusy=[]{return true;};Access::refresh(ui);expect(!ready(),"playback blocks absent-Current initialization");Access::workflow(ui);expect(!s.hasCurrentStateInitialized(),"direct callback cannot bypass stopped transport");s.auditionBusy={};
    s.persistenceLoading=true;Access::refresh(ui);expect(!ready(),"session loading blocks Propose");s.persistenceLoading=false;
    s.sourceContentInvalid=true;Access::refresh(ui);expect(!ready()&&Access::proposal(ui).actionText().contains("Resolve"),"changed source remains blocked with visible remedy");s.sourceContentInvalid=false;
    ui.removeFromDesktop();
}

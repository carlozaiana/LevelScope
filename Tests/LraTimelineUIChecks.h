void lraTimelineUIChecks(const juce::File& root)
{
    // Explicit view-only fixture: simultaneously gate-rejected supports. This
    // does not mint a session measurement or test estimator mathematics.
    {
        auto a=std::make_shared<proposal::ProgramAnalysis>();a->frames=48000*61;a->sampleRate=48000;a->identity="renderer-fixture";a->snapshotIdentity=a->identity;
        a->lra={std::nullopt,measurement::Status::belowGate};a->lraGate=-45;
        a->lraWindows={{0,144000,1e-9,0,measurement::OfflineMeasurement::loudness(1e-9),false,false},{4800,148800,1e-6,0,measurement::OfflineMeasurement::loudness(1e-6),false,false},{9600,153600,0,0,std::nullopt,false,false}};
        TimelineHeights h;h.lra={76,76,76};LraTimelineLanes renderer;renderer.setHeights(&h);renderer.setSize(900,renderer.preferredHeight());
        auto index=proposal::lra_timeline::Index::build(a);renderer.measured(index,ebuR128ProgrammeQc(),Applicability::preProducedLinearProgramme,"Original");renderer.viewport(0,144000,12,16);
        expect(renderer.status(1).contains("No eligible percentile population")&&!renderer.status(1).contains("Unstable"),"no percentile population is distinct from unstable available LRA");
        expect((renderer.displayedColumns().front().flags&(proposal::lra_timeline::absoluteRejected|proposal::lra_timeline::relativeRejected))==(proposal::lra_timeline::absoluteRejected|proposal::lra_timeline::relativeRejected),"overlapping gate categories reach one display column");
        juce::Image image(juce::Image::ARGB,900,renderer.getHeight(),true);juce::Graphics graphics(image);renderer.paint(graphics);const auto y=renderer.laneBounds(2).getY()+35;
        expect(image.getPixelAt(13,y+15)==juce::Colour(0xffb4bcc6)&&image.getPixelAt(13,y+25)==juce::Colour(0xff758392),"absolute and relative gate evidence have separate visible strips and colours in a mixed column");
        auto stream=root.getChildFile("gate-strip-renderer.png").createOutputStream();juce::PNGImageFormat png;expect(stream&&png.writeImageToStream(image,*stream),"mixed-gate rendering preserved");
        auto shortData=std::make_shared<proposal::ProgramAnalysis>(*a);shortData->frames=48000*2;shortData->lra.status=measurement::Status::insufficientDuration;
        renderer.measured(proposal::lra_timeline::Index::build(shortData),ebuR128ProgrammeQc(),Applicability::preProducedLinearProgramme,"Original");expect(renderer.status(2).contains("Insufficient programme duration"),"too-short material has its own diagnostic state");
    }
    const auto preferences=root.getChildFile("owned-view-settings.json");
    StandaloneMainComponent ui(emptyOutputCatalog(),noOutputSettings(),{},preferences);Access::stopTimer(ui);ui.setSize(1360,1100);
    ui.addToDesktop(juce::ComponentPeer::windowHasTitleBar);ui.setTopLeftPosition(30,30);ui.setAlwaysOnTop(true);ui.setVisible(true);ui.toFront(true);
    auto& s=Access::model(ui);auto& graph=Access::timeline(ui);auto& gain=Access::lane(ui);auto& lanes=graph.lraDiagnostics();
    const auto pump=[]{juce::MessageManager::getInstance()->runDispatchLoopUntil(10);};
    const auto waitMeasured=[&](MeasurementTimeline& t,const std::function<void()>& refresh)
    {
        const auto until=juce::Time::getMillisecondCounterHiRes()+10000;
        do{refresh();pump();}while(!t.lraDiagnostics().measuredIndex()&&juce::Time::getMillisecondCounterHiRes()<until);
        expect(bool(t.lraDiagnostics().measuredIndex()),"matching completed measured LRA index reaches the native timeline");
    };
    Access::refresh(ui);expect(!lanes.measuredIndex()&&lanes.status(0).contains("Confirm roles"),"empty document displays concise unavailable state");
    installWorkflowFixture(ui,programme(root.getChildFile("fixture"),2,true));juce::String error;
    auto settings=RegionPlanningSettings::innerDefaults(RegionMode::automatic,sectionInputIdentity(s.measurementContext(MeasurementDomain::source)));
    expect(s.setRegionPlanning(settings,error),"existing five context lanes available for resizing");
    auto noCeiling=s.jobRequirements;noCeiling.maximumLraTenths.reset();expect(s.setJobRequirements(noCeiling,error),"explicit no-LRA-target fixture");
    Access::refresh(ui);waitMeasured(graph,[&]{Access::refreshTimeline(ui);});
    const auto sourceIndex=lanes.measuredIndex();
    expect(lanes.status(0).contains("Set an LRA target")&&lanes.status(1).contains("P95")&&lanes.status(2).contains("P10"),"no target suppresses potential-passage guidance but retains measured percentile lanes");
    auto req=s.jobRequirements;req.maximumLraTenths=std::max(1,int(std::floor(*s.source.analysis->detail.lra.value*10))-10);expect(s.setJobRequirements(req,error),"controlled LRA target below this fixture's measured baseline");Access::refresh(ui);
    expect(lanes.measuredIndex()==sourceIndex&&lanes.measuredTarget().state==proposal::lra_timeline::TargetState::exceeds,"target edits query the accepted measurement without remeasuring");
    expect(std::any_of(lanes.displayedColumns().begin(),lanes.displayedColumns().end(),[](const auto& c){return (c.flags&proposal::lra_timeline::targetPassage)!=0;}),"measured potential-passage flags are drawn for the selected ceiling");
    graph.setViewRange(480000,2400000);const auto view=graph.viewRange();
    const auto before=s.manualState();const auto source=s.source.analysis,current=s.currentState.analysis;const auto revision=s.revision;const auto undo=(s.manualHistory?s.manualHistory->undoCount():0);const auto audioPosition=Access::audio(ui).positionFrames();
    const auto unchanged=[&]{expect(s.manualState()==before&&s.source.analysis==source&&s.currentState.analysis==current&&s.revision==revision&&(s.manualHistory?s.manualHistory->undoCount():0)==undo&&Access::audio(ui).positionFrames()==audioPosition,"view resize/inspection preserves exact audio, measurements, history and transport");expect(graph.viewRange()==view,"view resize preserves exact horizontal zoom and scroll");};
    const auto values=[](const TimelineHeights& h){std::vector<int> v{h.graph,h.gain};v.insert(v.end(),h.context.begin(),h.context.end());v.insert(v.end(),h.lra.begin(),h.lra.end());return v;};
    const auto drag=[&](TimelineDivider& d,int delta)
    {
        const auto now=juce::Time::getCurrentTime();const juce::Point<float> p{float(d.getWidth()/2),4},q=p.translated(0,float(delta));
        const auto event=[&](juce::Point<float> point,int mods){return juce::MouseEvent{juce::Desktop::getInstance().getMainMouseSource(),point,juce::ModifierKeys(mods),0,0,0,0,0,&d,&d,now,p,now,1,false};};
        d.mouseDown(event(p,juce::ModifierKeys::leftButtonModifier));d.mouseDrag(event(q,juce::ModifierKeys::leftButtonModifier));d.mouseUp(event(q,0));pump();
    };
    std::vector<TimelineDivider*> dividers{&graph.graphHeightDivider(),&graph.gainHeightDivider()};
    for(size_t n=0;n<5;++n)dividers.push_back(gain.contextDivider(n));for(size_t n=0;n<3;++n)dividers.push_back(&lanes.divider(n));
    for(size_t n=0;n<dividers.size();++n)
    {
        const auto old=values(graph.laneHeights());expect(dividers[n]&&dividers[n]->isVisible(),"every displayed graph/lane has its own direct divider");drag(*dividers[n],24);const auto changed=values(graph.laneHeights());
        for(size_t k=0;k<old.size();++k)expect(changed[k]==old[k]+(k==n?24:0),"one divider changes only its own chosen height");
        expect(graph.graphBounds().getHeight()==graph.laneHeights().graph&&gain.plotBounds().getHeight()==graph.laneHeights().gain,"graph and gain heights follow their own independent settings");unchanged();
    }
    const auto sizes=graph.laneHeights();drag(graph.graphHeightDivider(),-100000);expect(graph.laneHeights().graph==80&&graph.laneHeights().gain==sizes.gain,"graph minimum protects a usable plot without shrinking another lane");drag(graph.graphHeightDivider(),sizes.graph-80);
    drag(lanes.divider(0),-100000);expect(graph.laneHeights().lra[0]==44&&graph.laneHeights().lra[1]==sizes.lra[1],"LRA minimum retains title and plotted evidence while other lanes stay fixed");drag(lanes.divider(0),sizes.lra[0]-44);
    const auto stored=graph.laneHeights();expect(TimelineViewPreferences(preferences).main==stored,"completed gestures save local main-window heights");
    const auto graphHeight=graph.graphBounds().getHeight();Access::moveDivider(ui,70);expect(graph.graphBounds().getHeight()==graphHeight&&graph.laneHeights()==stored,"outer divider changes timeline viewport capacity, not graph or lane preferences");unchanged();
    ui.setSize(1120,950);Access::refresh(ui);auto& viewport=Access::timelineViewport(ui);
    expect(ui.getLocalBounds().contains(viewport.getBounds())&&viewport.getHeight()>100&&graph.getHeight()==graph.preferredHeight(),"small main window keeps timeline viewport in bounds and full content scrollable");
    expect(viewport.getVerticalScrollBar().isVisible(),"overflow has a visible vertical scrollbar");
    viewport.setViewPosition(0,graph.getHeight());expect(viewport.getViewArea().getBottom()>=lanes.getBottom(),"last diagnostic lane is reachable by vertical scrolling");
    viewport.setViewPosition(0,0);expect(viewport.getViewArea().getY()==0,"graph/navigation controls remain reachable at the top");unchanged();
    const auto q=lanes.queryCount();juce::Image scratch(juce::Image::ARGB,graph.getWidth(),graph.getHeight(),true);juce::Graphics paint(scratch);lanes.paint(paint);lanes.paint(paint);expect(lanes.queryCount()==q,"painting consumes cached columns without rebuilding diagnostics");
    const auto columns=sourceIndex->columns(int64_t(view.first),int64_t(view.second),graph.graphBounds().getWidth(),lanes.measuredTarget());
    expect(columns.size()==lanes.displayedColumns().size()&&std::equal(columns.begin(),columns.end(),lanes.displayedColumns().begin(),[](const auto& a,const auto& b){return a.flags==b.flags&&a.minLufs==b.minLufs&&a.maxLufs==b.maxLufs;}),"drawn columns use the exact shared time range, width and Logic categories");
    const auto text=lanes.inspectionAt(48000*20,1);expect(text.contains("Measured 3 s LRA windows")&&text.contains("LUFS")&&text.contains("not per-section compliance")&&text.contains("inaudibility"),"exact window inspection identifies measured evidence and its interpretation limits");
    const juce::Point<float> p{float(graph.graphBounds().getCentreX()),float(lanes.laneBounds(1).getY()+40)};const auto now=juce::Time::getCurrentTime();
    lanes.mouseDown(juce::MouseEvent{juce::Desktop::getInstance().getMainMouseSource(),p,{},0,0,0,0,0,&lanes,&lanes,now,p,now,1,false});expect(graph.inspectionText().contains("Measured 3 s LRA windows")&&graph.inspectionText().contains("samples ["),"diagnostic click opens exact read-only inspection");unchanged();
    for(const auto px:{float(graph.graphBounds().getCentreX()),float(graph.graphBounds().getRight()-1)})
    {const juce::Point<float> point{px,p.y};lanes.mouseDown(juce::MouseEvent{juce::Desktop::getInstance().getMainMouseSource(),point,{},0,0,0,0,0,&lanes,&lanes,now,point,now,1,false});expect(lanes.inspectedFrame()==std::min(s.sourceDocument.media.frames-1,timelineFrameAt(px,graph.graphBounds().getX(),graph.graphBounds().getWidth(),view.first,view.second,s.sourceDocument.media.frames)),"graph and LRA inspection use the same centre/right-edge sample mapping");}
    graph.setViewRange(0,double(s.sourceDocument.media.frames));
    {const juce::Point<float> eof{float(graph.graphBounds().getRight()-1),p.y};lanes.mouseDown(juce::MouseEvent{juce::Desktop::getInstance().getMainMouseSource(),eof,{},0,0,0,0,0,&lanes,&lanes,now,eof,now,1,false});expect(lanes.inspectedFrame()==s.sourceDocument.media.frames-1,"rightmost full-programme pixel inspects the last actual sample, not exclusive EOF");expect(graph.inspectionText().contains("terminal padding")&&graph.inspectionText().contains("display clipped to programme"),"EOF inspection identifies measured padding and actual source support");}
    graph.setViewRange(view.first,view.second);
    Access::refresh(ui);expect(Access::inspectionDetail(ui).getText().contains("Measured 3 s LRA windows"),"routine main refresh retains a clicked LRA report for reading");
    auto targetAgain=s.jobRequirements;targetAgain.maximumLraTenths=*targetAgain.maximumLraTenths+1;expect(s.setJobRequirements(targetAgain,error),"change target after LRA inspection");Access::refresh(ui);expect(Access::inspectionDetail(ui).getText().contains("LRA inspection expired")&&!lanes.inspectedFrame(),"target change expires the LRA-owned main report and inspection marker");
    auto& e=control<juce::TextEditor>(gain,"Timeline section gain dB");choose(gain,"Gain lane selection",1);gain.selectFrame(600000,true);e.grabKeyboardFocus();e.setText("-7",false);e.onTextChange();drag(graph.graphHeightDivider(),12);expect(s.manualState()==before&&e.getText()!="-7","divider interaction cancels an unconfirmed gain without applying it");
    graph.setViewRange(view.first,view.second);imageOf(ui,root.getChildFile("main-overflow.png"));imageOf(graph,root.getChildFile("all-lanes.png"));

    auto& panel=Access::sections(ui);auto& st=panel.timelineView();waitMeasured(st,[&]{panel.synchronize();Access::refreshTimeline(ui);});
    st.setViewRange(view.first,view.second);const auto sectionView=st.viewRange();const auto sectionInitial=st.laneHeights();
    drag(st.graphHeightDivider(),32);drag(st.gainHeightDivider(),16);for(size_t n=0;n<5;++n)drag(*st.attachedGainLane()->contextDivider(n),8);for(size_t n=0;n<3;++n)drag(st.lraDiagnostics().divider(n),12);
    expect(st.laneHeights().graph==sectionInitial.graph+32&&st.laneHeights().gain==sectionInitial.gain+16&&graph.laneHeights()!=st.laneHeights(),"Sections graph and lane sizes are independent from main-window sizes");
    expect(st.viewRange()==sectionView&&s.manualState()==before&&s.currentState.analysis==current,"Sections resizing preserves time view and measured audio");
    expect(TimelineViewPreferences(preferences).sections==st.laneHeights(),"Sections heights saved separately");
    imageOf(panel,root.getChildFile("sections-lanes.png"));
    panel.selectInspectionDomain(MeasurementDomain::currentState);waitMeasured(st,[&]{panel.synchronize();});
    auto& sl=st.lraDiagnostics();const juce::Point<float> sp{float(st.graphBounds().getCentreX()),float(sl.laneBounds(1).getY()+40)};
    const auto inspectSection=[&]{sl.mouseDown(juce::MouseEvent{juce::Desktop::getInstance().getMainMouseSource(),sp,{},0,0,0,0,0,&sl,&sl,now,sp,now,1,false});};inspectSection();
    panel.synchronize();expect(st.inspectionText().contains("Measured 3 s")&&control<juce::TextEditor>(panel,"Section window inspection").getText().contains("Measured 3 s"),"Sections click opens its own measured report and routine refresh keeps it readable");
    auto changedTarget=s.jobRequirements;changedTarget.maximumLraTenths=*changedTarget.maximumLraTenths+1;s.setJobRequirements(changedTarget,error);panel.synchronize();expect(st.inspectionText().contains("expired")&&control<juce::TextEditor>(panel,"Section window inspection").getText().contains("expired")&&!sl.inspectedFrame(),"Sections target change expires the copied LRA report");
    inspectSection();s.setCurrentGainDb(s.currentStateDocument.staticGainDb+.1);panel.synchronize();expect(!sl.measuredIndex()&&st.inspectionText().contains("expired")&&control<juce::TextEditor>(panel,"Section window inspection").getText().contains("expired"),"Sections Current gain edit clears both coloured evidence and copied report");
    measuredComposition(s,MeasurementDomain::currentState);panel.synchronize();waitMeasured(st,[&]{panel.synchronize();});inspectSection();panel.selectInspectionDomain(MeasurementDomain::source);panel.synchronize();expect(st.inspectionText().contains("expired"),"Sections domain switch expires its previous LRA report");
    Access::hideSections(ui);
    {StandaloneMainComponent reopened(emptyOutputCatalog(),noOutputSettings(),{},preferences);Access::stopTimer(reopened);expect(Access::timeline(reopened).laneHeights()==graph.laneHeights(),"a new app shell restores local view heights without portable session data");}

    Access::inspectCurrent(ui);waitMeasured(graph,[&]{Access::refreshTimeline(ui);});expect(lanes.status(1).contains("Current"),"Current domain clearly labelled");
    expect(s.setCurrentGainDb(s.currentStateDocument.staticGainDb+.5),"deliberate processing edit for stale display check");Access::refresh(ui);
    expect(!lanes.measuredIndex()&&lanes.displayedColumns().empty()&&lanes.status(0).contains("Needs remeasurement"),"edited Current clears old measured diagnostics immediately");
    Access::inspectOriginal(ui);waitMeasured(graph,[&]{Access::refreshTimeline(ui);});expect(lanes.status(1).contains("Original")&&lanes.measuredIndex()->measured().snapshotIdentity==sourceIndex->measured().snapshotIdentity,"Original retains its own valid measured diagnostics after Current edits");
    s.sourceContentInvalid=true;Access::refresh(ui);expect(!lanes.measuredIndex()&&lanes.status(0).contains("resolve"),"unresolved source clears all measured lanes");s.sourceContentInvalid=false;Access::refresh(ui);waitMeasured(graph,[&]{Access::refreshTimeline(ui);});
    const auto live=s.manualState();const auto heightsBeforeOpen=graph.laneHeights();JobControl save;const auto sessionFile=root.getChildFile("owned.levelscope");expect(sf::save(sf::capture(s),sessionFile,save).succeeded(),"save owned edited session with local display settings kept separate");
    Access::open(ui,sessionFile);waitFile(ui);expect(graph.laneHeights()==heightsBeforeOpen&&s.manualState()==live,"Open retains local view heights and exact edited processing");
    Access::inspectCurrent(ui);Access::refresh(ui);expect(!lanes.measuredIndex(),"reopened unverified Current does not acquire measured diagnostic authority");
    expect(ui.performCommand(StandaloneMainComponent::measureCurrent),"explicit Current remeasurement starts normally");Access::refresh(ui);expect(!lanes.measuredIndex()&&lanes.status(0).contains("progress"),"in-progress edited Current has no completed diagnostic highlights");waitJob(ui);waitMeasured(graph,[&]{Access::refreshTimeline(ui);});expect(lanes.status(1).contains("Current"),"independent remeasurement restores matching Current diagnostics");
    Access::resetLayout(ui);expect(graph.laneHeights()==TimelineHeights{}&&TimelineViewPreferences(preferences).main==TimelineHeights{}&&TimelineViewPreferences(preferences).sections==TimelineHeights{},"Reset layout restores and saves main and Sections height defaults");
    auto invalid=preferences.getSiblingFile("invalid-view.json");invalid.replaceWithText("{\"schema\":1,\"main\":{\"graph\":\"bad\"}}");expect(TimelineViewPreferences(invalid).main==TimelineHeights{},"invalid view preferences use safe defaults");
    expect(Access::stopped(ui),"diagnostic and resize checks never connect output or start playback");
    ui.removeFromDesktop();
}

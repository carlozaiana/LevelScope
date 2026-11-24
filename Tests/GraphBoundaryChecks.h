// Focused graph entry-point checks; processing transactions remain the existing ones.
void graphBoundaryChecks(const juce::File& root)
{
    auto req=fixture(root.getChildFile("graph-boundary.wav"),48000,2,240000);StandaloneSessionModel s;
    expect(s.setValidatedSource(req.context.sourceFile,req.context.media),"graph fixture imported");s.sourceDocument.layout=req.context.fileLayout;
    std::shared_ptr<const waveform::Overview> overview;measure(s,Domain::source,&overview);s.synchronizeSections();
    juce::String error;expect(s.sectionDocument.split(s.sectionDocument.map().sections.front().id,48003,error),"annotation deliberately coincides with a processing cut");
    auto before=s.manualState();before.globalDb=2;
    processing::SectionGainPlan plan;plan.inputIdentity=before.map.inputIdentity;plan.sampleRate=48000;plan.frames=240000;
    plan.intervals={{0,48003,0},{48003,96000,-4},{96000,144000,-2},{144000,240000,0}};before.gains=plan;
    expect(s.commitManualEdit(before,"Graph fixture processing"),"graph fixture has actual applied gains");measure(s,Domain::currentState);before=s.manualState();
    const auto verified=s.currentState.analysis,original=s.source.analysis;const auto revision=s.revision,undoCount=s.manualHistory->undoCount();
    TimelineGainLane lane(s);MeasurementTimeline ui;ui.setSize(1150,680);ui.attachGainLane(lane);ui.setProgramme(240000,48000,"graph");
    ui.setWaveform(overview,true,{"L","R"});lane.synchronize("graph");ui.bindSectionEditing(s,[&](uint64_t){lane.synchronize("graph");});ui.setReferenceMeasurement({},true);
    juce::ComboBox* mode=nullptr;for(auto* c:lane.getChildren())if(auto* b=dynamic_cast<juce::ComboBox*>(c))mode=b;
    expect(mode!=nullptr,"graph mode selector exists");mode->setSelectedId(2,juce::sendNotificationSync);
    ui.addToDesktop(juce::ComponentPeer::windowHasTitleBar);ui.setVisible(true);ui.toFront(true);pump(20);ui.focusForTransport();
    const auto x=[&](int64_t frame){const auto v=ui.viewRange();return float(ui.graphBounds().getX()+(frame-v.first)*(ui.graphBounds().getWidth()-1)/(v.second-v.first));};
    const auto point=[&](int64_t f,bool graph=true){return juce::Point<float>{x(f)+2,float(graph?ui.graphBounds().getCentreY():lane.getY()+lane.plotBounds().getY()+4)};};
    const auto event=[&](juce::Point<float> p,int mods){auto t=juce::Time::getCurrentTime();return juce::MouseEvent{juce::Desktop::getInstance().getMainMouseSource(),p,juce::ModifierKeys(mods),0,0,0,0,0,&ui,&ui,t,p,t,1,false};};
    const auto down=[&](juce::Point<float> p){ui.mouseDown(mouse(ui,p));};
    const auto move=[&](juce::Point<float> p){ui.mouseDrag(mouse(ui,p));};
    const auto up=[&](juce::Point<float> p){ui.mouseUp(mouse(ui,p,false));};
    const auto refresh=[&]{lane.synchronize("graph");ui.refreshSectionEditing();};
    expect(ui.appliedGraphBoundaryHandle(point(48003))==1&&lane.appliedBoundaryHandle(point(48003,false).toInt()-lane.getPosition())==1,"graph and gain lane target the same shared internal cut");
    expect(!ui.appliedGraphBoundaryHandle({x(0),float(ui.graphBounds().getCentreY())})&&!ui.appliedGraphBoundaryHandle({x(240000),float(ui.graphBounds().getCentreY())}),"graph has no movable programme endpoints");
    ui.mouseMove(mouse(ui,point(48003),false));expect(ui.getMouseCursor()==juce::MouseCursor::LeftRightResizeCursor,"graph boundary has a horizontal resize affordance");
    auto p=point(48003),q=p.translated(27,0);down(p);up(p);
    expect(s.manualState()==before&&s.currentState.analysis==verified&&s.revision==revision&&s.manualHistory->undoCount()==undoCount,"off-centre no-movement graph click retains exact sample, measurement and Undo");
    down(p);move(q);move(p);up(p);
    expect(s.manualState()==before&&s.currentState.analysis==verified&&s.manualHistory->undoCount()==undoCount,"out-and-back graph gesture is an exact no-op");
    down(p);move(q);const auto preview=s.displayedSectionGains();const auto cut=preview->intervals[1].begin;
    expect(s.appliedBoundaryDrag&&s.appliedBoundaryDrag->right==1&&s.manualState()==before&&s.currentState.analysis==verified,"graph preview uses existing transaction without changing applied audio or measurement");
    expect(cut==48003+int64_t(std::llround(27.*240000/(ui.graphBounds().getWidth()-1)))&&preview->intervals[0].end==cut,"graph drag uses exact original cut plus sample-position delta");
    expect(lane.appliedBoundaryAt(x(cut))==1&&ui.appliedGraphBoundaryHandle({x(cut),p.y})==1&&ui.transitionSpans().front()==std::make_pair(cut-preview->radius(1),cut+preview->radius(1)),"both views and transition band follow the same full edited plan");
    ui.setTransitionShading(false);expect(!ui.transitionSpans().empty(),"active boundary preview shows its ramp even when normal transition shading is hidden");ui.setTransitionShading(true);
    imageOf(ui,root.getChildFile("graph-boundary-preview.png"));
    expect(ui.keyPressed(juce::KeyPress(juce::KeyPress::escapeKey)),"Escape cancels graph drag before release");up(q);
    expect(!s.appliedBoundaryDrag&&s.manualState()==before&&s.currentState.analysis==verified&&s.revision==revision,"Escape and later release preserve processing and verification");
    down(p);move(q);up(q);const auto fromGraph=s.manualState();
    expect(s.manualHistory->undoCount()==undoCount+1&&s.currentStateNeedsRemeasurement()&&!s.currentState.hasMeasurement&&s.source.analysis==original,"graph release is one Undo and invalidates Current only");
    expect(fromGraph.globalDb==before.globalDb&&proposal::sectionGeometryKey(fromGraph.map)==proposal::sectionGeometryKey(before.map)&&fromGraph.gains->intervals[1].db==-4&&fromGraph.gains->intervals[2].begin==96000,"coincident annotation, both gains, global trim and other cuts remain fixed");
    expect(s.undoManual()&&s.manualState()==before,"Undo restores graph-edited envelope");refresh();
    auto lp=point(48003,false),lq=lp.translated(27,0);down(lp);move(lq);up(lq);
    expect(s.manualState()==fromGraph,"gain-lane and graph gestures produce byte-equivalent manual processing state");
    expect(s.undoManual()&&s.undoManual(true)&&s.manualState()==fromGraph,"Redo restores the same graph-equivalent cut");s.undoManual();refresh();
    ui.setViewRange(47980,48080);p=point(48003);q=p.translated(float(3.*(ui.graphBounds().getWidth()-1)/100.),0);
    down(p);move(q);up(q);expect(s.currentStateDocument.sectionGains->intervals[1].begin==48006,"deep overview zoom moves an applied cut by exactly three audio samples");s.undoManual();refresh();
    ui.setViewRange(0,240000);p=point(48003);q=p.translated(20,0);
    for(int reason=0;reason<5;++reason)
    {
        down(p);move(q);expect(s.appliedBoundaryDrag.has_value(),"graph cancellation fixture starts a real preview");
        if(reason==0)ui.setViewRange(0,120000);
        if(reason==1)ui.focusLost(juce::Component::focusChangedDirectly);
        if(reason==2)ui.setVisible(false);
        if(reason==3)ui.setSize(1160,680);
        if(reason==4)ui.setReferenceMeasurement({},false);
        up(q);expect(!s.appliedBoundaryDrag&&s.manualState()==before,"viewport, focus, hide, resize and Original inspection cancel graph preview");
        ui.setVisible(true);ui.setSize(1150,680);ui.setReferenceMeasurement({},true);ui.setViewRange(0,240000);ui.focusForTransport();refresh();p=point(48003);q=p.translated(20,0);
    }
    down(p);move(q);mode->setSelectedId(1,juce::sendNotificationSync);up(q);
    expect(!s.appliedBoundaryDrag&&s.manualState()==before,"changing selection mode cancels unfinished graph edit");
    expect(!ui.appliedGraphBoundaryHandle(point(48003)),"section mode cannot grab applied processing boundaries");
    ui.refreshSectionEditing();ui.mouseDown(event(p,juce::ModifierKeys::leftButtonModifier|juce::ModifierKeys::commandModifier));ui.mouseDrag(event(q,juce::ModifierKeys::leftButtonModifier|juce::ModifierKeys::commandModifier));ui.mouseUp(event(q,juce::ModifierKeys::commandModifier));
    expect(s.boundaryBatch&&!s.appliedBoundaryDrag&&s.manualState()==before,"coincident line in section mode retains the existing Command boundary-batch transaction");ui.keyPressed(juce::KeyPress(juce::KeyPress::escapeKey));refresh();
    mode->setSelectedId(2,juce::sendNotificationSync);ui.setReferenceMeasurement({},false);p=point(48003);down(p);move(p.translated(10,0));up(p.translated(10,0));
    expect(!ui.appliedGraphBoundaryHandle(p)&&s.manualState()==before&&!s.appliedBoundaryDrag&&!ui.editPrompt()->active(),"Original graph cannot edit Current or silently move coincident section markers in region mode");
    ui.setReferenceMeasurement({},true);ui.setViewRange(24000,144000);p=point(48003);const auto viewport=ui.viewRange();ui.mouseDown(event(p,juce::ModifierKeys::leftButtonModifier|juce::ModifierKeys::altModifier));ui.mouseDrag(event(p.translated(15,0),juce::ModifierKeys::leftButtonModifier|juce::ModifierKeys::altModifier));ui.mouseUp(event(p.translated(15,0),juce::ModifierKeys::altModifier));
    expect(!s.appliedBoundaryDrag&&s.manualState()==before&&ui.viewRange()!=viewport,"Alt-drag pans even over an applied boundary");
    int seeks=0;ui.onRequestSeek=[&](int64_t,const juce::String&){++seeks;return juce::String();};ui.setTransport(0,true,{},"graph-transport");click(ui,ui.rulerBounds().getCentre().toFloat());
    expect(seeks==1&&s.manualState()==before,"timecode ruler seeking does not grab graph handles");
    p=point(48003);const auto range=ui.viewRange();const auto frac=(p.x-ui.graphBounds().getX())/(ui.graphBounds().getWidth()-1);const auto anchor=range.first+frac*(range.second-range.first);
    juce::MouseWheelDetails wheel;wheel.deltaY=.1f;ui.mouseWheelMove(event(p,juce::ModifierKeys::commandModifier),wheel);const auto zoomed=ui.viewRange();
    expect(std::abs(zoomed.first+frac*(zoomed.second-zoomed.first)-anchor)<.01&&s.manualState()==before,"Command-wheel keeps its pointer anchor without editing the region");
    ui.setViewRange(0,240000);p=point(48003);q=p.translated(20,0);
    for(int busy=0;busy<3;++busy)
    {
        if(busy==0)s.auditionBusy=[]{return true;};if(busy==1)s.persistenceLoading=true;if(busy==2)s.job.begin(JobStage::analysis,s.revision);
        down(p);move(q);up(q);expect(!s.appliedBoundaryDrag&&s.manualState()==before,"graph cannot bypass playback, load or measurement guards");
        s.auditionBusy={};s.persistenceLoading=false;if(busy==2)s.job.finish(s.job.ticket,s.revision,ErrorCode::none,"done");refresh();
    }
    down(p);move(q);++s.revision;up(q);expect(!s.appliedBoundaryDrag&&s.manualState()==before,"stale processing revision cannot commit graph gesture");refresh();
    imageOf(ui,root.getChildFile("graph-boundary-handles.png"));ui.removeFromDesktop();
}

// Native reproduction uses owned fixtures; never the user's running session.
void gainInputNavigationChecks(const juce::File& root)
{
    StandaloneMainComponent ui(emptyOutputCatalog(),noOutputSettings());Access::stopTimer(ui);ui.setSize(1360,1100);
    ui.addToDesktop(juce::ComponentPeer::windowHasTitleBar);ui.setTopLeftPosition(40,40);ui.setAlwaysOnTop(true);ui.setVisible(true);ui.toFront(true);
    auto& s=Access::model(ui);auto& lane=Access::lane(ui);auto& graph=Access::timeline(ui);
    installWorkflowFixture(ui,programme(root.getChildFile("fixture"),2));
    auto state=s.manualState();juce::String error;
    expect(s.sectionDocument.split(1,480000,error)&&s.sectionDocument.split(2,960000,error)&&s.sectionDocument.split(3,1440000,error),"four editable sections");
    state=s.manualState();state.globalDb=1.742;processing::SectionGainPlan gains;gains.inputIdentity=state.map.inputIdentity;gains.frames=state.map.frames;gains.sampleRate=48000;
    gains.intervals={{0,480000,-1},{480000,960000,-3},{960000,1440000,-2},{1440000,gains.frames,0}};state.gains=gains;
    expect(s.commitManualEdit(state,"Owned gain fixture"),"distinct adjacent applied gains");measuredComposition(s,MeasurementDomain::currentState);Access::inspectCurrent(ui);Access::refresh(ui);
    const auto pump=[] {juce::MessageManager::getInstance()->runDispatchLoopUntil(35);};
    const auto mouse=[](MeasurementTimeline& t,juce::Point<float> p,int mods){const auto now=juce::Time::getCurrentTime();return juce::MouseEvent{juce::Desktop::getInstance().getMainMouseSource(),p,juce::ModifierKeys(mods),0,0,0,0,0,&t,&t,now,p,now,1,false};};
    const auto x=[](MeasurementTimeline& t,int64_t f){const auto v=t.viewRange();return float(t.graphBounds().getX()+(f-v.first)*(t.graphBounds().getWidth()-1)/(v.second-v.first));};
    const auto entry=[&](TimelineGainLane& l,const juce::String& value){auto& e=control<juce::TextEditor>(l,"Timeline section gain dB");e.grabKeyboardFocus();e.setText(value,false);e.onTextChange();};
    choose(lane,"Gain lane selection",1);lane.selectFrame(600000,true);graph.setViewRange(400000,1100000);
    const auto before=s.manualState();const auto verified=s.currentState.analysis;const auto count=s.manualHistory->undoCount();const auto revision=s.revision;
    for(const auto& value:juce::StringArray{" ","", "wrong", "99", "-7.5"})
    {
        entry(lane,value);const juce::Point<float> p{x(graph,700000),float(graph.rulerBounds().getCentreY())};graph.mouseDown(mouse(graph,p,juce::ModifierKeys::leftButtonModifier));graph.mouseUp(mouse(graph,p,0));pump();
        expect(Access::audio(ui).positionFrames()>699990&&Access::audio(ui).positionFrames()<700010,"outside ruler click seeks on the first click after unconfirmed gain");
        expect(s.manualState()==before&&s.currentState.analysis==verified&&s.revision==revision&&s.manualHistory->undoCount()==count,"outside click cancels valid or invalid text without audio, verification or history mutation");
        expect(control<juce::TextEditor>(lane,"Timeline section gain dB").getText().getDoubleValue()==-3,"applied value restored before navigation");
    }
    auto& field=control<juce::TextEditor>(lane,"Timeline section gain dB");
    // Dispatch a real JUCE peer mouse click, including global listeners, to a
    // button which deliberately does not steal keyboard focus.
    auto& zoomButton=control<juce::TextButton>(graph,"Zoom +");const auto takesFocus=zoomButton.getMouseClickGrabsKeyboardFocus();zoomButton.setMouseClickGrabsKeyboardFocus(false);
    entry(lane,"-7");const auto priorView=graph.viewRange();const auto z=ui.getLocalPoint(&zoomButton,zoomButton.getLocalBounds().getCentre()).toFloat();
    ui.toFront(true);pump();
    const auto now=juce::Time::currentTimeMillis();
    ui.getPeer()->handleMouseEvent(juce::MouseInputSource::InputSourceType::mouse,z,juce::ModifierKeys(),0,0,now);
    ui.getPeer()->handleMouseEvent(juce::MouseInputSource::InputSourceType::mouse,z,juce::ModifierKeys(juce::ModifierKeys::leftButtonModifier),0,0,now+1);
    ui.getPeer()->handleMouseEvent(juce::MouseInputSource::InputSourceType::mouse,z,juce::ModifierKeys(),0,0,now+2);pump();
    expect(graph.viewRange()!=priorView&&field.getText().getDoubleValue()==-3&&s.manualState()==before,"outside click cancels even when the clicked control retains editor focus, then performs zoom");zoomButton.setMouseClickGrabsKeyboardFocus(takesFocus);
    graph.setViewRange(400000,1500000);entry(lane," ");const juce::Point<float> inspectPoint{x(graph,1300000),float(graph.graphBounds().getCentreY())};graph.mouseDown(mouse(graph,inspectPoint,juce::ModifierKeys::leftButtonModifier));graph.mouseUp(mouse(graph,inspectPoint,0));
    expect(lane.selectedSectionId()==state.map.sections[2].id&&s.manualState()==before,"first outside graph click cancels whitespace and selects the clicked section");
    lane.selectFrame(600000,true);graph.setViewRange(400000,1100000);entry(lane,"-9");const auto p=graph.graphBounds().getCentre().toFloat(),q=p.translated(20,0);const auto panView=graph.viewRange();
    graph.mouseDown(mouse(graph,p,juce::ModifierKeys::leftButtonModifier|juce::ModifierKeys::altModifier));graph.mouseDrag(mouse(graph,q,juce::ModifierKeys::leftButtonModifier|juce::ModifierKeys::altModifier));graph.mouseUp(mouse(graph,q,juce::ModifierKeys::altModifier));
    expect(graph.viewRange()!=panView&&s.manualState()==before&&field.getText().getDoubleValue()==-3,"Alt-pan cancels valid entry and pans on the same gesture");
    for(const auto& value:juce::StringArray{" ","wrong","99"})
    {
        entry(lane,value);field.onReturnKey();expect(s.manualState()==before&&lane.inspectionText().isNotEmpty(),"invalid Enter never applies");
        bool message=false;for(auto* c:lane.getChildren())if(auto* l=dynamic_cast<juce::Label*>(c))message|=l->getText().contains("Gain not applied");
        expect(message,"invalid Enter has concise visible validation");
        expect(ui.getPeer()->handleKeyPress(juce::KeyPress(juce::KeyPress::escapeKey)),"native Escape handles invalid focused gain");pump();
        expect(field.getText().getDoubleValue()==-3&&s.currentState.analysis==verified&&graph.hasKeyboardFocus(false),"Escape restores applied gain and reliable timeline focus");
    }
    entry(lane,"-7");graph.focusForTransport();pump();expect(field.getText().getDoubleValue()==-3&&s.manualState()==before,"focus loss cancels valid text");
    lane.selectFrame(100000,true);lane.toggleSection(600000);expect(lane.selectedCount()==2&&field.getText().isEmpty()&&field.getTextToShowWhenEmpty()=="Mixed","heterogeneous multi-selection displays Mixed");
    entry(lane,"-8");field.onEscapeKey();pump();expect(lane.selectedCount()==2&&field.getText().isEmpty()&&field.getTextToShowWhenEmpty()=="Mixed"&&s.manualState()==before,"Escape restores Mixed and both individual gains");
    entry(lane,"-9");choose(lane,"Gain lane selection",2);expect(lane.selectsAppliedRegions()&&s.manualState()==before&&s.currentState.analysis==verified,"selection dropdown cancels unconfirmed text without applying");choose(lane,"Gain lane selection",1);lane.selectFrame(600000,true);
    entry(lane,"-9");click(lane,"Create sections from applied regions");expect(s.manualState()==before&&s.currentState.analysis==verified,"matching section copy cancels entry without changing any applied gain");
    entry(lane,"-9");Access::workflow(ui);expect(Access::operation(ui).original==before,"Propose baseline excludes unconfirmed valid gain");Access::cancel(ui);completeOperation(ui);expect(s.manualState()==before,"cancelled Propose preserves applied gains");
    lane.selectFrame(600000,true);entry(lane,"-2.5");field.onReturnKey();pump();expect(section_edit::uniform(&*s.currentStateDocument.sectionGains,500000,900000)==-2.5&&!s.currentState.hasMeasurement&&s.manualHistory->undoCount()==count+1,"Enter confirms one ordinary undoable gain edit and invalidates Current");
    expect(ui.performCommand(StandaloneMainComponent::undoEdit)&&s.manualState()==before,"Undo restores confirmed gain");measuredComposition(s,MeasurementDomain::currentState);Access::refresh(ui);

    // Main and Sections window: both edges, every popup choice, exact viewport.
    for(bool separate:{false,true})
    {
        std::unique_ptr<StandaloneSectionsPanel> panel;
        if(separate){panel=std::make_unique<StandaloneSectionsPanel>(s);panel->suspend();panel->setSize(1120,1000);panel->addToDesktop(juce::ComponentPeer::windowHasTitleBar);panel->setVisible(true);panel->selectInspectionDomain(MeasurementDomain::currentState);}
        auto& t=separate?panel->timelineView():graph;auto& l=*t.attachedGainLane();if(separate){t.setProgramme(s.sourceDocument.media.frames,48000,"owned-sections");t.setReferenceMeasurement({},true);t.refreshSectionEditing();}
        const auto refresh=[&]{if(panel)panel->synchronize();else Access::refresh(ui);};
        refresh();choose(l,"Gain lane selection",1);
        for(bool left:{true,false})for(int answer:{0,1,2})
        {
            refresh();l.selectFrame(600000,true);if(panel)control<juce::ListBox>(*panel,"Programme sections").selectRow(1);
            const auto old=s.manualState();const auto measurement=s.currentState.analysis;const auto undo=s.manualHistory->undoCount();
            const auto edge=left?480000:960000;t.setViewRange(edge-10000,edge+15000);const auto view=t.viewRange();
            auto& e=control<juce::TextEditor>(l,"Timeline section gain dB");entry(l,answer==0?" ":"-8");
            const juce::Point<float> a{x(t,edge),float(t.graphBounds().getCentreY())},b{x(t,edge+2000),a.y};
            t.mouseDown(mouse(t,a,juce::ModifierKeys::leftButtonModifier));t.mouseDrag(mouse(t,b,juce::ModifierKeys::leftButtonModifier));t.mouseUp(mouse(t,b,0));
            expect(t.editPrompt()&&t.editPrompt()->active()&&s.manualState()==old,"first boundary drag cancels unconfirmed gain then opens normal confirmation");
            t.editPrompt()->window()->exitModalState(answer);for(int n=0;n<60&&t.editPrompt()->active();++n)pump();pump();refresh();pump();
            expect(!t.editPrompt()->active()&&t.viewRange()==view,"Move, Keep and Cancel preserve exact zoom and scroll after asynchronous refresh");
            expect(l.selectedSectionId()==old.map.sections[1].id,"left and right edges retain the intended selected section");
            expect(t.hasKeyboardFocus(false),"boundary popup returns keyboard focus to graph, not numeric field");
            const auto gainAfterDialog=e.getText();int toggles=0;auto oldToggle=t.onTogglePlayback;t.onTogglePlayback=[&]{++toggles;};
            expect(t.getPeer()->handleKeyPress(juce::KeyPress(juce::KeyPress::spaceKey)),"Space handled after popup");pump();expect(toggles==1&&e.getText()==gainAfterDialog,"Space after popup targets transport without replacing gain text");t.onTogglePlayback=oldToggle;
            if(answer==0)expect(s.manualState()==old&&s.currentState.analysis==measurement&&s.manualHistory->undoCount()==undo,"Cancel preserves processing, measurement and Undo");
            else
            {
                const auto edited=s.manualState();expect(edited.map.sections[0].end==(left?edited.map.sections[1].begin:old.map.sections[0].end)&&s.manualHistory->undoCount()==undo+1,"confirmed boundary retains existing one-edit transaction");
                if(answer==1)expect(edited.gains->audioKey()==old.gains->audioKey()&&s.currentState.analysis==measurement,"Keep gain preserves exact processing and measurement");
                else expect(edited.gains->audioKey()!=old.gains->audioKey()&&!s.currentState.hasMeasurement,"Move gain changes real envelope and invalidates Current");
                expect(t.keyPressed(juce::KeyPress('Z',juce::ModifierKeys::commandModifier,0))&&s.manualState()==old&&t.viewRange()==view,"Undo preserves viewport and restores boundary operation");
                expect(t.keyPressed(juce::KeyPress('Z',juce::ModifierKeys::commandModifier|juce::ModifierKeys::shiftModifier,0))&&s.manualState()==edited&&t.viewRange()==view,"Redo preserves viewport and reapplies only boundary operation");
                expect(s.undoManual()&&s.manualState()==old,"restore fixture for next popup choice");measuredComposition(s,MeasurementDomain::currentState);refresh();
            }
            e.grabKeyboardFocus();e.selectAll();const auto textBefore=e.getText();t.getPeer()->handleKeyPress(juce::KeyPress(juce::KeyPress::spaceKey));pump();
            expect(e.getText()==textBefore,"bare Space in deliberately focused numeric editor cannot erase the selected value");t.focusForTransport();
        }
        t.setViewRange(460000,500000);click(t,"Whole programme");expect(t.viewRange().first==0&&t.viewRange().second==s.sourceDocument.media.frames,"explicit Fit retains normal navigation");click(t,"Zoom +");expect(t.viewRange().second-t.viewRange().first<s.sourceDocument.media.frames,"explicit zoom still changes range");
        imageOf(separate?static_cast<juce::Component&>(*panel):static_cast<juce::Component&>(ui),root.getChildFile(separate?"sections-stable.png":"main-stable.png"));if(panel)panel->removeFromDesktop();
    }
    expect(Access::stopped(ui),"all checks leave transport stopped and open no physical output");ui.removeFromDesktop();
}

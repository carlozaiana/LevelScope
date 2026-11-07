// Annotation deletion and applied-region removal exercise the real shared lane/model.
void deletionChecks(const juce::File& root)
{
    auto req=fixture(root.getChildFile("delete.wav"),48000,2,240013);StandaloneSessionModel s;
    expect(s.setValidatedSource(req.context.sourceFile,req.context.media),"deletion programme loaded");s.sourceDocument.layout=req.context.fileLayout;
    expect(s.initializeCurrentStateFromSource(),"deletion Current initialized");s.synchronizeSections();juce::String error;auto& doc=s.sectionDocument;
    expect(doc.split(doc.map().sections.front().id,96000,error)&&doc.split(doc.map().sections.back().id,192000,error),"three deletion annotations");
    s.setCurrentGainDb(2);expect(s.beginSectionTrimDraft()&&s.setSectionDraftGain(0,-1)&&s.setSectionDraftGain(1,-3)&&s.setSectionDraftGain(2,2)&&s.applySectionTrimDraft(),"known three-region applied plan");
    const auto measure=[&](Domain domain){expect(s.job.begin(JobStage::analysis,s.revision),"deletion measurement starts");auto context=s.measurementContext(domain);JobControl control;auto result=OfflineAnalysisService{}.analyze(context,domain,control);expect(result.value&&s.job.finish(context.ticket,s.revision,result.error,result.message)&&s.acceptMeasurement(*result.value),"deletion measurement independently accepted");};
    measure(Domain::source);measure(Domain::currentState);
    const auto verified=s.currentState.analysis,original=s.source.analysis;const auto plan=s.currentStateDocument.sectionGains->audioKey();const auto revision=s.revision;
    TimelineGainLane lane(s);MeasurementTimeline ui;ui.setSize(1150,640);ui.attachGainLane(lane);ui.setProgramme(req.context.media.frames,48000,"delete");ui.setWaveform({},false,{"L","R"});lane.synchronize("delete");
    ui.addToDesktop(juce::ComponentPeer::windowHasTitleBar);ui.setVisible(true);ui.focusForTransport();
    const auto select=[&](int64_t frame,int kind){for(auto* c:lane.getChildren())if(auto* box=dynamic_cast<juce::ComboBox*>(c))box->setSelectedId(kind,juce::sendNotificationSync);lane.selectAt(frame,lane.plotBounds().getCentreY());ui.focusForTransport();};
    select(100000,1);const auto deletedId=doc.map().sections[1].id;const auto nextId=doc.map().sections[2].id;
    expect(ui.keyPressed(juce::KeyPress(juce::KeyPress::backspaceKey))&&doc.map().sections.size()==2&&doc.map().sections[0].end==192000&&doc.map().sections[1].id==nextId,"Backspace removes selected interior annotation into predecessor, retaining complete coverage");
    expect(s.currentState.analysis==verified&&s.source.analysis==original&&s.revision==revision&&s.currentStateDocument.sectionGains->audioKey()==plan,"annotation deletion preserves exact applied plan and both measurements");
    expect(!ui.keyPressed(juce::KeyPress(juce::KeyPress::backspaceKey))&&doc.map().sections.size()==2,"held/repeated Backspace cannot cascade without another gain-lane selection");
    button(lane,"Undo delete")->onClick();expect(doc.map().sections.size()==3&&doc.map().sections[1].id==deletedId,"Undo delete restores original annotation IDs/bounds");
    select(0,1);ui.keyPressed(juce::KeyPress(juce::KeyPress::backspaceKey));expect(doc.map().sections.size()==2&&doc.map().sections.front().id==deletedId&&doc.map().sections.front().begin==0,"first annotation deletes into successor, retaining successor identity");button(lane,"Undo delete")->onClick();
    select(220000,1);ui.keyPressed(juce::KeyPress(juce::KeyPress::backspaceKey));expect(doc.map().sections.size()==2&&doc.map().sections.back().end==doc.map().frames,"last annotation deletes into predecessor without truncating audio");button(lane,"Undo delete")->onClick();
    select(100000,1);auto* field=laneEditor(lane);field->grabKeyboardFocus();expect(!ui.keyPressed(juce::KeyPress(juce::KeyPress::backspaceKey))&&doc.map().sections.size()==3,"Backspace in gain text entry cannot delete a section");
    select(100000,2);ui.keyPressed(juce::KeyPress(juce::KeyPress::backspaceKey));
    expect(s.sectionTrimDraft&&s.sectionTrimDraft->processingRegions&&s.sectionTrimDraft->plan.intervals[1].db==0&&s.sectionTrimDraft->plan.intervals[0].db==-1&&s.sectionTrimDraft->plan.intervals[2].db==2,"applied-region deletion drafts only selected local gain removal");
    expect(s.currentState.analysis==verified&&s.currentStateDocument.sectionGains->audioKey()==plan&&s.currentStateDocument.staticGainDb==2,"removal draft preserves measurement, applied plan and global gain");
    expect(!laneEditor(lane)->isEnabled()&&button(lane,"Apply gains")->isEnabled(),"region removal has explicit Apply/Discard and cannot enter the section-index editor");
    imageOf(ui,root.getChildFile("region-removal-draft.png"));button(lane,"Discard gains")->onClick();expect(!s.sectionTrimDraft&&s.currentState.analysis==verified,"discard region removal retains verified Current");
    // Different annotation geometry must never make a region edit reset unrelated applied gains.
    expect(doc.moveEnd(doc.map().sections[0].id,90000,error),"move annotation independently of processing");lane.synchronize("delete");select(100000,2);ui.keyPressed(juce::KeyPress(juce::KeyPress::backspaceKey));
    expect(s.sectionTrimDraft&&s.sectionTrimDraft->plan.intervals[1].begin==96000&&s.sectionTrimDraft->plan.intervals[0].db==-1,"region removal uses applied interval geometry despite changed annotations");
    button(lane,"Apply gains")->onClick();expect(!s.sectionTrimDraft&&!s.currentState.hasMeasurement&&s.currentStateNeedsRemeasurement()&&s.source.analysis==original&&s.currentStateDocument.staticGainDb==2,"Apply removal invalidates only Current and retains global gain/Source");
    const auto& applied=*s.currentStateDocument.sectionGains;expect(applied.intervals.size()==3&&applied.intervals[0].db==-1&&applied.intervals[1].db==0&&applied.intervals[2].db==2,"Apply preserves neighbouring plateaus and leaves selected local gain at unity");
    select(100000,2);ui.keyPressed(juce::KeyPress(juce::KeyPress::backspaceKey));expect(!s.sectionTrimDraft,"already-unity region removal is a no-op");
    select(200000,2);s.job.begin(JobStage::analysis,s.revision);ui.keyPressed(juce::KeyPress(juce::KeyPress::backspaceKey));expect(!s.sectionTrimDraft,"active job blocks gain deletion");s.job.finish(s.job.ticket,s.revision,ErrorCode::none,"done");
    select(200000,2);ui.keyPressed(juce::KeyPress(juce::KeyPress::backspaceKey));expect(s.sectionTrimDraft.has_value(),"second removal draft");doc.moveEnd(doc.map().sections[0].id,91000,error);lane.synchronize("delete");expect(!button(lane,"Apply gains")->isEnabled(),"changed annotation geometry makes removal draft stale");button(lane,"Discard gains")->onClick();
    ui.removeFromDesktop();
}

// Focused regression of the user's draw -> gain -> resize/delete -> redraw workflow.
void manualWorkflowChecks(const juce::File& root)
{
    auto req=fixture(root.getChildFile("manual-workflow.wav"),48000,2,240013);StandaloneSessionModel s;
    expect(s.setValidatedSource(req.context.sourceFile,req.context.media),"manual workflow source");s.sourceDocument.layout=req.context.fileLayout;s.setSelectedTargetProfileIndex(1);
    expect(s.programmeApplicability==proposal::Applicability::preProducedLinearProgramme,"fresh imports evaluate supported file limits without confirmation");
    measure(s,Domain::source);if(!s.hasCurrentStateInitialized())s.initializeCurrentStateFromSource();
    expect(s.evaluation(Domain::source).verdict!=proposal::Verdict::notEvaluated,"completed Source is evaluated without checking an applicability control");
    s.setCurrentGainDb(4);measure(s,Domain::currentState);s.synchronizeSections();juce::String error;
    auto& document=s.sectionDocument;expect(document.split(1,96000,error)&&document.split(document.map().sections.back().id,192000,error),"draw manual section inside programme");
    const auto id=document.map().sections[1].id;
    TimelineGainLane lane(s);MeasurementTimeline ui;ui.setSize(1150,650);ui.attachGainLane(lane);ui.setProgramme(req.context.media.frames,48000,"manual");ui.setWaveform({},false,{"L","R"});
    ui.bindSectionEditing(s,[&](uint64_t){lane.synchronize("manual");});lane.synchronize("manual");
    ui.addToDesktop(juce::ComponentPeer::windowHasTitleBar);ui.setVisible(true);ui.focusForTransport();
    auto* editor=laneEditor(lane);expect(button(lane,"Set draft")==nullptr,"redundant Set draft button removed");
    const auto select=[&](int64_t frame){lane.selectAt(frame,lane.plotBounds().getCentreY());ui.focusForTransport();};
    const auto answer=[&](SectionEditPrompt& prompt,int choice,const char* png=nullptr)
    {
        expect(prompt.window()!=nullptr,"section edit has explicit gain choices");if(png)imageOf(*prompt.window(),root.getChildFile(png));
        prompt.window()->exitModalState(choice);juce::MessageManager::getInstance()->runDispatchLoopUntil(20);
        expect(!prompt.active()&&juce::Component::getCurrentlyModalComponent()==nullptr,"section choice closes completely");ui.focusForTransport();
    };
    select(120000);const auto verified=s.currentState.analysis;
    editor->setText("-2",false);editor->onTextChange();editor->onFocusLost();
    expect(s.sectionDraftFresh()&&section_edit::uniform(&s.sectionTrimDraft->plan,96000,192000)==-2&&s.currentState.analysis==verified&&!s.currentStateDocument.sectionGains,"leaving numeric field drafts -2 without changing Current");
    button(lane,"Apply gains")->onClick();
    expect(s.currentStateDocument.staticGainDb==4&&s.currentStateDocument.sectionGains->gainAt(120000)==-2&&s.currentStateNeedsRemeasurement(),"Apply yields global +4/local -2 and unverified Current");
    measure(s,Domain::currentState);lane.synchronize("manual");const auto applied=s.currentStateDocument.sectionGains->audioKey();const auto measured=s.currentState.analysis;
    lane.selectAt(120000,lane.plotBounds().getBottom()+5);
    expect(editor->isEnabled()&&lane.inspectionText().contains("Section"),"clicking applied readout preserves manual section mode");
    // Cancellation and annotation-only edits never rewrite the gain envelope.
    ui.changeSectionEnd(document.map().sections[0].id,100000);answer(*ui.editPrompt(),0,"boundary-gain-choice.png");
    expect(document.map().sections[0].end==96000&&!s.sectionTrimDraft&&s.currentState.analysis==measured,"Cancel leaves annotation, audio and verification intact");
    ui.changeSectionEnd(document.map().sections[0].id,100000);answer(*ui.editPrompt(),1);
    expect(document.map().sections[0].end==100000&&s.currentStateDocument.sectionGains->audioKey()==applied&&s.currentState.analysis==measured,"Keep gain in place changes only annotation boundary");
    // Split within the adjusted annotation: the applied envelope remains exact.
    uint64_t nested=0;expect(document.defineInterval(120000,140000,nested,error),"draw nested section");s.documentChanged();lane.synchronize("manual");select(125000);
    expect(editor->getText()=="-2"&&s.currentStateDocument.sectionGains->audioKey()==applied,"nested section inherits existing local plateau without moving ramps");
    stageGain(lane,"-3");
    expect(section_edit::uniform(&s.sectionTrimDraft->plan,96000,120000)==-2&&section_edit::uniform(&s.sectionTrimDraft->plan,140000,192000)==-2&&section_edit::uniform(&s.sectionTrimDraft->plan,0,96000)==0,"editing nested section preserves every other existing gain");
    button(lane,"Discard gains")->onClick();
    // This section crosses the old processing cut at 96000.
    select(92000);expect(editor->getText().isEmpty()&&editor->getTextToShowWhenEmpty()=="Mixed","section crossing multiple applied gains displays Mixed");
    // Move a real gain boundary and confirm that both parts remain a single draft.
    ui.changeSectionEnd(document.map().sections[0].id,96000);answer(*ui.editPrompt(),1);
    const auto before=document.map();ui.changeSectionEnd(before.sections[0].id,104000);answer(*ui.editPrompt(),2);
    expect(s.sectionTrimDraft&&s.sectionTrimDraft->rollbackAnnotations&&document.map().sections[0].end==104000&&s.currentStateDocument.sectionGains->audioKey()==applied&&s.currentState.analysis==measured,"Move gain stages paired geometry with applied Current unchanged");
    expect(section_edit::uniform(&s.sectionTrimDraft->plan,96000,104000)==0&&section_edit::uniform(&s.sectionTrimDraft->plan,104000,192000)==-2,"expanding side fills only the swept nominal span");
    expect(!s.setCurrentGainDb(5),"paired edit prevents an unrelated global change before Apply/Discard");
    bool saveBlocked=false;try{sessionfile::capture(s);}catch(const std::exception&){saveBlocked=true;}expect(saveBlocked,"cannot save provisional geometry as applied processing");
    proposal::markers::Preview markerPreview;juce::String markerError;
    expect(!preflightMarkerImport(s,markerPreview,markerError)&&markerError.contains("paired section/gain"),"marker import cannot replace provisional paired annotations");
    imageOf(ui,root.getChildFile("paired-boundary-draft.png"));s.discardSectionTrimDraft();lane.synchronize("manual");ui.refreshSectionEditing();
    expect(document.map().sections[0].end==96000&&s.currentState.analysis==measured&&!s.sectionTrimDraft,"Discard restores section boundary and exact measured Current");
    ui.changeSectionEnd(document.map().sections[0].id,104000);answer(*ui.editPrompt(),2);expect(s.applySectionTrimDraft(),"Apply paired boundary edit");lane.synchronize("manual");
    expect(document.map().sections[0].end==104000&&s.currentStateDocument.sectionGains->gainAt(99000)==0&&s.currentStateDocument.sectionGains->gainAt(130000)==-2&&s.currentStateDocument.staticGainDb==4&&s.currentStateNeedsRemeasurement(),"paired Apply commits matching geometry/gain and invalidates only Current");
    measure(s,Domain::currentState);lane.synchronize("manual");const auto newMeasured=s.currentState.analysis;
    // Delete the nested section and local gain together, then Discard both.
    select(125000);expect(lane.deleteSelected(),"delete selected nested annotation");answer(lane.editPrompt(),2,"delete-gain-choice.png");
    expect(s.sectionTrimDraft&&s.sectionTrimDraft->rollbackAnnotations&&std::none_of(document.map().sections.begin(),document.map().sections.end(),[&](const auto& r){return r.id==nested;}),"delete-with-gain shows provisional annotation removal");
    expect(section_edit::uniform(&s.sectionTrimDraft->plan,120000,140000)==0&&s.currentState.analysis==newMeasured,"delete-with-gain drafts local unity without changing verified Current");
    s.discardSectionTrimDraft();lane.synchronize("manual");expect(std::any_of(document.map().sections.begin(),document.map().sections.end(),[&](const auto& r){return r.id==nested;})&&s.currentState.analysis==newMeasured,"Discard brings deleted annotation and applied gain back together");
    select(125000);lane.deleteSelected();answer(lane.editPrompt(),2);expect(s.applySectionTrimDraft(),"apply section plus gain removal");lane.synchronize("manual");
    JobControl job;auto saved=sessionfile::capture(s);const auto file=root.getChildFile("manual-paired.levelscope");
    expect(sessionfile::save(saved,file,job).succeeded(),"save final paired manual edit with matching processing provenance");auto read=sessionfile::read(file,job);auto checked=sessionfile::verify(*read.document,job);
    expect(checked.prepared&&checked.prepared->model.currentStateNeedsRemeasurement()&&checked.prepared->model.currentStateDocument.sectionGains->audioKey()==s.currentStateDocument.sectionGains->audioKey(),"save/reopen retains exact edited processing and unverified Current");
    saved.applicability=proposal::Applicability::unconfirmed;StandaloneSessionModel legacy;
    expect(sessionfile::install(saved,legacy,error)&&legacy.programmeApplicability==proposal::Applicability::preProducedLinearProgramme,"legacy unchecked sessions adopt automatic file-limit evaluation");
    saved.applicability=proposal::Applicability::liveProgramme;
    expect(sessionfile::install(saved,legacy,error)&&legacy.programmeApplicability==proposal::Applicability::liveProgramme&&legacy.evaluation(Domain::currentState).verdict!=proposal::Verdict::pass,"explicit unsupported legacy applicability is not silently promoted to a pass");
    // A modal answer cannot alter a newer programme/map.
    ui.changeSectionEnd(document.map().sections.front().id,105000);const auto pending=document.map().digest();
    expect(document.rename(document.map().sections.front().id,"Changed while question open",error),"external edit while question open");
    const auto newer=document.map().digest();answer(*ui.editPrompt(),2);
    expect(newer!=pending&&document.map().digest()==newer&&!s.sectionTrimDraft,"stale modal answer cannot edit newer annotations or audio");
    StandaloneSectionsPanel other(s);other.suspend();juce::TextEditor* field=nullptr;
    for(auto* child:other.getChildren())if(child->getName()=="Section draft gain dB")field=dynamic_cast<juce::TextEditor*>(child);
    expect(field&&field->isEnabled()&&!button(other,"Set draft gain")&&!button(other,"Draft section gains"),"Sections window offers direct gain entry without draft buttons");
    const auto oldPlan=s.currentStateDocument.sectionGains->audioKey();field->setText("-1.25",false);field->onTextChange();field->onFocusLost();
    expect(s.sectionDraftFresh()&&section_edit::uniform(&s.sectionTrimDraft->plan,0,104000)==-1.25&&s.currentStateDocument.sectionGains->audioKey()==oldPlan,"Sections focus loss stages gain with applied processing unchanged");
    button(other,"Discard draft")->onClick();expect(!s.sectionTrimDraft&&s.currentStateDocument.sectionGains->audioKey()==oldPlan,"Sections Discard preserves applied plan");
    ui.removeFromDesktop();
}

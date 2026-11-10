void manualWorkflowChecks(const juce::File& root)
{
    auto req=fixture(root.getChildFile("manual.wav"),48000,2,240013);StandaloneSessionModel s;
    expect(s.programmeApplicability==proposal::Applicability::preProducedLinearProgramme,"automatic file-limit scope needs no confirmation");
    expect(s.setValidatedSource(req.context.sourceFile,req.context.media),"manual fixture imported");s.sourceDocument.layout=req.context.fileLayout;measure(s,Domain::source);s.synchronizeSections();
    const auto original=s.source.analysis;auto& document=s.sectionDocument;juce::String error;
    TimelineGainLane lane(s);MeasurementTimeline ui;ui.setSize(1150,640);ui.attachGainLane(lane);ui.bindSectionEditing(s,[&](uint64_t){lane.synchronize("manual");});ui.setProgramme(240013,48000,"manual");lane.synchronize("manual");
    ui.addToDesktop(juce::ComponentPeer::windowHasTitleBar);ui.setVisible(true);ui.toFront(true);pump(20);ui.focusForTransport();
    auto state=s.manualState();state.globalDb=4;expect(s.commitManualEdit(state,"Global trim"),"global trim participates in Undo");lane.synchronize("manual");
    const auto graph=ui.graphBounds();draw(ui,float(graph.getX()+graph.getWidth()*.2),float(graph.getX()+graph.getWidth()*.4));
    expect(document.map().sections.size()==3,"draw section commits immediately");const auto middle=document.map().sections[1];lane.selectFrame(middle.begin,false);stageGain(lane,"-2");
    const auto gainKey=s.currentStateDocument.sectionGains->audioKey();
    auto* field=laneEditor(lane);field->grabKeyboardFocus();field->setText("-3",false);field->onTextChange();
    const bool returnHandled=ui.getPeer()->handleKeyPress(juce::KeyPress(juce::KeyPress::returnKey));
    const auto returnStart=juce::Time::getMillisecondCounterHiRes();while(!ui.hasKeyboardFocus(false)&&juce::Time::getMillisecondCounterHiRes()-returnStart<1500)pump(10);
    expect(returnHandled&&ui.hasKeyboardFocus(false),"native Return commits gain and returns focus to the timeline");
    expect(ui.getPeer()->handleKeyPress(juce::KeyPress('Z',juce::ModifierKeys::commandModifier,0))&&s.currentStateDocument.sectionGains->audioKey()==gainKey,"native Command-Z immediately after Return undoes processing, not just editor text");

    draw(ui,float(graph.getX()+graph.getWidth()*.6),float(graph.getX()+graph.getWidth()*.8));
    expect(document.map().sections.size()==5&&!s.sectionTrimDraft&&s.currentStateDocument.sectionGains->audioKey()==gainKey&&s.currentStateNeedsRemeasurement(),"draw a second section immediately after gain entry without Apply or measurement");
    const auto undoCount=s.manualHistory->undoCount();expect(s.undoManual()&&document.map().sections.size()==3&&s.currentStateDocument.sectionGains->audioKey()==gainKey,"Undo drawing preserves audio");expect(s.undoManual(true)&&document.map().sections.size()==5&&s.manualHistory->undoCount()==undoCount,"Redo drawing restores exact sections");
    const auto before=s.manualState();
    const auto answer=[&](int value){auto* prompt=ui.editPrompt();expect(prompt->window()!=nullptr,"boundary movement offers keep/move gain");prompt->window()->exitModalState(value);const auto start=juce::Time::getMillisecondCounterHiRes();while(prompt->active()&&juce::Time::getMillisecondCounterHiRes()-start<1500)pump(10);expect(!prompt->active(),"boundary choice dismisses");};
    const auto boundary=document.map().sections.front().end;ui.changeSectionEnd(document.map().sections.front().id,boundary+500);imageOf(*ui.editPrompt()->window(),root.getChildFile("move-boundary-choices.png"));answer(1);
    expect(s.currentStateDocument.sectionGains->audioKey()==gainKey&&document.map().sections.front().end==boundary+500,"Keep gain in place moves only annotations");s.undoManual();ui.refreshSectionEditing();
    ui.changeSectionEnd(document.map().sections.front().id,boundary+500);answer(2);
    expect(section_edit::uniform(&*s.currentStateDocument.sectionGains,boundary,boundary+500)==0&&!s.sectionTrimDraft,"Move gain immediately extends the expanding side");expect(s.undoManual()&&s.manualState()==before,"boundary Undo restores both geometry and exact processing");
    // Hundreds of journal entries remain available without retaining programme data in RAM.
    const auto id=document.map().sections[1].id;const auto count=s.manualHistory->undoCount();
    for(int n=0;n<260;++n)expect(s.setManualSectionGains({id},n%2?-3:-4),"journal accepts another direct edit");
    expect(s.manualHistory->undoCount()==count+260,"Undo has no arbitrary 100-step cap");
    for(int n=0;n<260;++n)expect(s.undoManual(),"journal reads and undoes earlier edit");
    expect(s.manualState()==before&&s.source.analysis==original&&!s.currentState.hasMeasurement,"long Undo restores exact processing without falsely restoring measurement authority");
    expect(s.setManualSectionGains({id},-1)&&!s.canUndoManual(true),"new edit after Undo truncates the redo branch");
    auto saved=sessionfile::capture(s);JobControl job;const auto path=root.getChildFile("direct.levelscope");
    expect(sessionfile::save(saved,path,job).succeeded(),"direct edits can be saved without a draft resolution step");auto read=sessionfile::read(path,job);auto checked=sessionfile::verify(*read.document,job);
    expect(checked.prepared&&checked.prepared->model.currentStateNeedsRemeasurement()&&checked.prepared->model.currentStateDocument.sectionGains->audioKey()==s.currentStateDocument.sectionGains->audioKey()&&!checked.prepared->model.canUndoManual(),"reopen preserves exact gains, independent verification and a fresh session-local history");
    saved.applicability=proposal::Applicability::unconfirmed;StandaloneSessionModel legacy;expect(sessionfile::install(saved,legacy,error)&&legacy.programmeApplicability==proposal::Applicability::preProducedLinearProgramme,"legacy unchecked session adopts automatic file-limit scope");
    saved.applicability=proposal::Applicability::liveProgramme;expect(sessionfile::install(saved,legacy,error)&&legacy.programmeApplicability==proposal::Applicability::liveProgramme,"unsupported saved scope is not silently promoted");
    // Gain-lane and right-waveform marquee share fully enclosed section semantics.
    lane.synchronize("manual");const auto& sections=document.map().sections;const auto begin=sections[1].begin,end=sections[3].end;
    const auto x=[&](int64_t f){return graph.getX()+float(double(f)*(graph.getWidth()-1)/240013.);};
    const auto rightMouse=[&](float px){const auto t=juce::Time::getCurrentTime();juce::Point<float> p{px,float(graph.getCentreY())};return juce::MouseEvent{juce::Desktop::getInstance().getMainMouseSource(),p,juce::ModifierKeys(juce::ModifierKeys::rightButtonModifier),0,0,0,0,0,&ui,&ui,t,p,t,1,false};};
    ui.refreshSectionEditing();ui.mouseDown(rightMouse(x(begin)-.1f));ui.mouseDrag(rightMouse(x(end)+.1f));ui.mouseUp(rightMouse(x(end)+.1f));
    expect(lane.selectedCount()==3,"right-drag selects fully enclosed waveform sections without drawing or panning");
    imageOf(ui,root.getChildFile("multi-selection.png"));const auto selectedState=s.manualState();
    expect(lane.deleteSelected()&&!lane.editPrompt().active()&&document.map().sections.size()==1&&!s.currentStateDocument.sectionGains&&s.currentStateDocument.staticGainDb==4,"marquee deletion resets all enclosed sections and merges zero neighbours in one action");
    expect(s.undoManual()&&s.manualState()==selectedState,"one Undo restores the entire multi-section deletion");ui.removeFromDesktop();
}
void gainPreviewChecks(const juce::File& root)
{
    WaveformCache cache;WaveformBinding source{"preview-runtime","original",{1,1},Domain::source},current{"preview-runtime","unity",{1,1},Domain::currentState};
    auto progress=cache.begin(source);auto data=progress->begin({512,48000,2});juce::AudioBuffer<float> samples(2,512);samples.clear();samples.setSample(0,4,-.8f);samples.setSample(0,200,1.2f);
    expect(data&&data->append(samples.getArrayOfReadPointers(),2,512,0),"preview uses actual measured-sample overview");progress->finish(data);expect(cache.accept(source,source),"complete source cache accepted");
    ProcessingPlan unity;cache.synchronize(source,current,true,true,&unity);const auto memory=cache.retainedBytes();const auto fallback=cache.gainPreview(source);
    expect(fallback.data==data&&fallback.previewBase&&cache.retainedBytes()==memory,"preview shares source overview without another slab");
    ProcessingPlan changed;processing::SectionGainPlan gain;gain.inputIdentity=juce::String::repeatedString("a",64);gain.sampleRate=48000;gain.frames=512;gain.intervals={{0,128,0},{128,512,-6}};changed.sectionGains=gain;
    MeasurementTimeline ui;ui.setSize(1150,470);ui.setProgramme(512,48000,"preview");ui.setWaveform(data,true,{"L","R"});const auto before=ui.displayedWaveform();ui.setReferenceMeasurement({},true);ui.setGainPreview(fallback.previewBase,changed);ui.setPeakTarget({},false);
    expect(ui.isGainPreview()&&ui.peakVisibility().getToggleState()&&ui.peakVisibility().isEnabled(),"unverified gain preview retains default-on warnings without a target");
    const auto& after=ui.displayedWaveform();bool kept=false,scaled=false;
    for(size_t n=0;n<after.size();++n)if(after[n].available())
    {const auto& a=after[n].channels[0];const auto& b=before[n].channels[0];if(b.minimumFrame==4)kept|=a.minimum==b.minimum&&a.minimumFrame==4;if(b.maximumFrame==200)scaled|=std::abs(a.maximum-float(1.2*std::pow(10.,gain.gainAt(200)/20.)))<1e-6&&a.maximumFrame==200;}
    expect(kept&&scaled&&cache.retainedBytes()==memory,"preview preserves unchanged extrema and positions while applying the actual transition envelope at retained samples");
    expect(!ui.sampleExceedsTarget(.99,100)&&!ui.sampleExceedsTarget(1.,100)&&ui.sampleExceedsTarget(1.01,100),"No Limit uses strict greater-than-zero dBFS for provisional sample warnings");
    expect(ui.rangePeakText(0,512).contains("preview"),"preview cannot be inspected as an exact current regional measurement");imageOf(ui,root.getChildFile("gain-preview.png"));
    auto stale=source;stale.runtime="reopened";cache.synchronize(stale,current,true,true,&changed);expect(!cache.gainPreview(stale).data,"preview cannot survive a different session/source identity");
}

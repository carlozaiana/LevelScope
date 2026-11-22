void appliedBoundaryUIChecks(const juce::File& root,const juce::File& supplied)
{
    const auto originalHash=juce::SHA256(supplied).toHexString();JobControl job;
    auto read=sf::read(supplied,job);expect(read.document&&read.document->regionPlanning&&read.document->regionPlanning->schema==3,"applied-boundary native check uses the supplied schema-3 session");
    auto audio=juce::File(read.document->source.fullPath);
    if(!audio.existsAsFile())audio=supplied.getSiblingFile("supplied-audio.wav");
    expect(audio.existsAsFile(),"explicit original or staged supplied audio exists");
    const auto audioHash=juce::SHA256(audio).toHexString();const auto audioCopy=root.getChildFile("audio-copy.wav");
    expect(audio.copyFileTo(audioCopy)&&juce::SHA256(audioCopy).toHexString()==audioHash&&sf::relinkAudio(*read.document,audioCopy),"byte-identical audio copy supplies the native session check");
    const auto copy=root.getChildFile("editable-copy.levelscope");expect(sf::save(*read.document,copy,job).succeeded(),"native checks use an owned session copy");
    StandaloneMainComponent ui(emptyOutputCatalog(),noOutputSettings());Access::stopTimer(ui);Access::deterministicChooser(ui);
    ui.setSize(1320,1100);ui.addToDesktop(juce::ComponentPeer::windowHasTitleBar);ui.setVisible(true);
    Access::open(ui,copy);waitFile(ui);auto& s=Access::model(ui);auto& lane=Access::lane(ui);auto& timeline=Access::timeline(ui);
    expect(s.lastAppliedExplanation&&s.currentStateDocument.sectionGains&&s.currentStateDocument.sectionGains->intervals.size()>1,"reopened applied regions and complete measured explanation available");
    const auto history=encodeAppliedExplanation(*s.lastAppliedExplanation);const auto settings=s.regionPlanning->digest();const auto before=s.manualState();
    const auto cut=before.gains->intervals[1].begin;
    choose(lane,"Gain lane selection",2);timeline.setViewRange(std::max<int64_t>(0,cut-48000),std::min<int64_t>(before.gains->frames,cut+48000));timeline.focusForTransport();
    const auto range=timeline.viewRange();const auto px=float(timeline.graphBounds().getX()+(cut-range.first)*(timeline.graphBounds().getWidth()-1)/(range.second-range.first));
    const juce::Point<float> start{px,float(lane.getY()+lane.plotBounds().getY()+4)},end{px+32,start.y};
    const auto event=[&](juce::Point<float> p,bool down){const auto t=juce::Time::getCurrentTime();return juce::MouseEvent{juce::Desktop::getInstance().getMainMouseSource(),p,juce::ModifierKeys(down?juce::ModifierKeys::leftButtonModifier:0),0,0,0,0,0,&timeline,&timeline,t,p,t,1,false};};
    timeline.mouseDown(event(start,true));timeline.mouseDrag(event(end,true));Access::refresh(ui);
    expect(s.appliedBoundaryDrag&&s.manualState()==before,"main-window refresh retains a transient applied-envelope preview");
    expect(!ui.performCommand(StandaloneMainComponent::measureCurrent)&&!ui.performCommand(StandaloneMainComponent::newDocument)&&!ui.performCommand(StandaloneMainComponent::saveDocument)&&Access::stopped(ui),"drag guards main measurement, New, Save and stopped transport");
    auto& setup=Access::setup(ui);setup.refresh();
    for(const auto* name:{"Prepare audition","Play selected","Restart","Seek"})
    {auto& button=control<juce::TextButton>(setup,name);expect(!button.isEnabled(),"audio-setup transport also waits for the boundary gesture");button.onClick();}
    expect(!Access::audio(ui).preparing()&&Access::stopped(ui)&&s.appliedBoundaryDrag,"direct audio-setup callbacks cannot bypass pending-drag guards");
    imageOf(ui,root.getChildFile("main-applied-boundary-preview.png"));
    timeline.mouseUp(event(end,false));Access::refresh(ui);const auto edited=s.manualState();
    expect(!s.appliedBoundaryDrag&&edited.gains->audioKey()!=before.gains->audioKey()&&edited.globalDb==before.globalDb&&s.currentStateNeedsRemeasurement()&&!s.currentState.hasMeasurement,"main release edits actual processing and immediately removes Current verification");
    expect(s.appliedExplanationStatus().contains("processing edited since application")&&encodeAppliedExplanation(*s.lastAppliedExplanation)==history&&s.regionPlanning->digest()==settings,"measured explanation becomes historical while complete history and planner permissions remain byte-identical");
    expect(!s.canApplySectionProposal()&&!control<juce::TextButton>(Access::proposal(ui),"Apply proposal").isEnabled(),"historical measured candidate cannot apply to edited Current");
    imageOf(ui,root.getChildFile("main-applied-boundary-committed.png"));
    Access::save(ui);waitFile(ui);expect(!s.documentDirty(),"native save commits the edited envelope");
    Access::open(ui,copy);waitFile(ui);
    expect(s.manualState()==edited&&s.lastAppliedExplanation&&encodeAppliedExplanation(*s.lastAppliedExplanation)==history&&s.regionPlanning->digest()==settings&&s.appliedExplanationStatus().contains("processing edited since application")&&s.currentStateNeedsRemeasurement(),"native save/reopen preserves edited schema-3 processing and complete historical explanation");
    expect(juce::SHA256(supplied).toHexString()==originalHash&&juce::SHA256(audio).toHexString()==audioHash&&Access::stopped(ui),"supplied session unchanged; no physical device opened or playback started");
    ui.removeFromDesktop();
}

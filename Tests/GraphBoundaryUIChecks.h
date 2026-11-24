// Main-window graph gestures with the accepted session source and playback path.
// Small synthetic inputs, injected callbacks only; no device or latency benchmark.
void graphBoundaryPlayback(const juce::File& root,int channels)
{
    constexpr int frames=240013;
    auto r=fixture(root.getChildFile("graph-playback-"+juce::String(channels)+".wav"),48000,channels,frames);
    auto device=std::make_unique<Device>();auto* d=device.get();
    StandaloneMainComponent ui(emptyOutputCatalog(),noOutputSettings(),std::move(device));
    Access::stopTimer(ui);ui.setSize(1320,1100);ui.addToDesktop(juce::ComponentPeer::windowHasTitleBar);ui.setVisible(true);
    auto& s=Access::model(ui);auto& audio=Access::audio(ui);auto& graph=Access::timeline(ui);auto& lane=Access::lane(ui);
    expect(s.setValidatedSource(r.context.sourceFile,r.context.media),"graph playback fixture imported");s.sourceDocument.layout=r.context.fileLayout;
    Access::refresh(ui);Access::measureSource(ui);waitJob(ui);waitPreparation(audio);Access::tick(ui);
    expect(s.source.hasMeasurement&&Access::overviewComplete(ui)&&d->opens==0&&!audio.active(),"initial analysis and preparation remain device-free");
    s.synchronizeSections();auto state=s.manualState();state.globalDb=2;
    processing::SectionGainPlan plan;plan.inputIdentity=state.map.inputIdentity;plan.sampleRate=48000;plan.frames=frames;
    plan.intervals={{0,48003,0},{48003,96000,-4},{96000,144000,-2},{144000,frames,0}};state.gains=plan;
    expect(s.commitManualEdit(state,"Graph playback fixture"),"native graph has an applied envelope");measuredComposition(s,MeasurementDomain::currentState);
    expect(s.selectAuditionDomain(MeasurementDomain::currentState),"audition explicitly selects Current");
    Access::inspectCurrent(ui);Access::refresh(ui);Access::tick(ui);waitPreparation(audio);Access::tick(ui);
    const auto counts=s.verifiedAssets->counts();const auto original=s.source.analysis;const auto prior=s.currentState.analysis;
    expect(counts.copies==1&&audio.playbackReady()&&!audio.active(),"Current shares the single verified source asset before editing");
    choose(lane,"Gain lane selection",2);graph.setViewRange(0,frames);graph.focusForTransport();
    const auto cut=state.gains->intervals[1].begin;
    const juce::Point<float> start{float(graph.graphBounds().getX()+double(cut)*(graph.graphBounds().getWidth()-1)/frames),float(graph.graphBounds().getCentreY())},end=start.translated(35,0);
    const auto event=[&](juce::Point<float> p,bool down){const auto t=juce::Time::getCurrentTime();return juce::MouseEvent{juce::Desktop::getInstance().getMainMouseSource(),p,juce::ModifierKeys(down?juce::ModifierKeys::leftButtonModifier:0),0,0,0,0,0,&graph,&graph,t,p,t,1,false};};
    graph.mouseDown(event(start,true));graph.mouseDrag(event(end,true));
    expect(s.appliedBoundaryDrag&&s.currentState.analysis==prior,"native graph previews without measuring or mutating Current");
    graph.mouseUp(event(end,false));Access::refresh(ui);Access::tick(ui);
    expect(s.currentStateNeedsRemeasurement()&&!s.currentState.hasMeasurement&&s.source.analysis==original&&s.currentStateDocument.sectionGains->intervals[1].begin!=cut,"native graph commit changes the actual cut and immediately invalidates Current");
    expect(!audio.active()&&d->opens==0,"editing and background preparation cannot authorize sound");

    // Decode the fixture and evaluate the committed processor directly. Do not
    // measure Current before Play: the interaction must work while unverified.
    Capture expected(channels,frames),played(channels,frames);
    juce::WavAudioFormat wav;std::unique_ptr<juce::AudioFormatReader> reader(wav.createReaderFor(r.context.sourceFile.createInputStream().release(),true));
    expect(reader&&reader->read(&expected.data,0,frames,0,true,true),"read exact input samples for edited-envelope reference");
    ProgrammeGainProcessor processor(s.currentStateDocument.processingPlan());
    expect(processor.process(expected.data.getArrayOfWritePointers(),channels,frames,0),"evaluate the actual committed global gain and transition envelope");expected.next=frames;
    expect(audio.connect(48000).isEmpty(),"explicitly connect injected output for graph listening check");
    for(int c=0;c<channels;++c)expect(s.setOutputRoute(c,c),"explicit injected route");
    Access::refresh(ui);audio.setCaptureObserver(&played);Access::play(ui);awaitPlaying(audio,*d);drain(audio,*d,997);audio.setCaptureObserver(nullptr);
    expect(identical(expected,played),"next explicit Play is sample-exact edited Current across all transitions without prior remeasurement");
    expect(s.currentStateNeedsRemeasurement()&&!s.currentState.hasMeasurement,"listening does not confer measurement validity");
    const auto after=s.verifiedAssets->counts();
    expect(after.copies==counts.copies&&after.copiedBytes==counts.copiedBytes,"graph editing and next Play require no additional full source copy or verification");
    expect(ui.performCommand(StandaloneMainComponent::measureCurrent),"normal Re-measure Current remains available after graph listening");waitJob(ui);
    expect(s.currentState.hasMeasurement&&!s.currentStateNeedsRemeasurement(),"independent normal Current measurement verifies the edited processing");
    imageOf(ui,root.getChildFile("graph-current-"+juce::String(channels)+".png"));audio.disconnect();ui.removeFromDesktop();
}

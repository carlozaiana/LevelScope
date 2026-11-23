// Included by the existing native integration harness. All outputs are injected.
void playbackEngineChecks(const juce::File& root)
{
    auto r=fixture(root.getChildFile("warmup.wav"),48000,2,48017);
    auto device=std::make_unique<Device>();auto* d=device.get();StandaloneAudition audio(std::move(device));auto map=routes(2);
    audio.synchronize(r,map);
    expect(audio.prepare(r).isEmpty()&&audio.backgroundPreparing()&&!audio.active(),"background reader preparation is separate from transport intent");
    expect(d->opens==0&&!d->isPlaying(),"background preparation never opens or starts an output");
    expect(audio.prepare(r).isEmpty(),"same in-flight preparation is joined without a duplicate");
    expect(audio.connect(48000).isEmpty()&&audio.preparing(),"explicit Connect during background preparation leaves that preparation running");
    expect(audio.play(r,map).isEmpty()&&audio.active(),"explicit Play joins existing preparation");
    expect(!d->isPlaying(),"Play cannot start before the prepared request is adopted");
    audio.stop();expect(!audio.active(),"Stop clears transport intent immediately while preparation cleanup is outstanding");waitPreparation(audio);
    expect(!audio.active()&&!d->isPlaying()&&!audio.playbackReady(),"Stop revokes delayed Play and discards cancelled prefix");
    expect(audio.prepare(r).isEmpty(),"explicit retry after cancellation starts preparation");waitPreparation(audio);
    const auto* bytes=audio.preparedCopy();expect(bytes&&audio.playbackReady()&&!audio.active()&&!d->isPlaying(),"warm reader and bounded prefix ready while stopped");
    expect(audio.seek(12345,r,map).isEmpty()&&!audio.playbackReady(),"stopped seek invalidates old prefix without starting output");
    expect(audio.prepare(r).isEmpty(),"new position prepares independently");
    expect(audio.seek(23456,r,map).isEmpty(),"stopped seek can cancel background preparation");waitPreparation(audio);
    expect(audio.positionFrames()==23456&&!d->isPlaying(),"late preparation completion cannot restore old seek position");
    expect(audio.prepare(r).isEmpty(),"new seek position can be warmed");waitPreparation(audio);
    expect(audio.preparedCopy()==bytes,"warmup at another position retains verified source bytes");
    expect(audio.play(r,map).isEmpty(),"prepared playback begins only from explicit Play");
    // A picture preroll can delay audio start; a prefilled reader cannot bypass it.
    audio.poll(false);expect(!d->isPlaying()&&audio.getState()==AuditionState::buffering,"warm prefix retains picture preroll gate");
    awaitPlaying(audio,*d);audio.stop();

    r.domain=MeasurementDomain::currentState;r.context.processing.staticGainDb=-4;++r.context.ticket.sessionRevision;
    processing::SectionGainPlan local;local.sampleRate=48000;local.frames=r.context.media.frames;local.inputIdentity=sectionInputIdentity(r.context);
    local.intervals={{0,1800,0},{1800,4100,-6},{4100,48017,-2}};
    r.context.processing.schemaVersion=2;r.context.processing.engineVersion=processing::SectionGainPlan::version;r.context.processing.sectionGains=local;
    audio.synchronize(r,map);expect(audio.prepare(r).isEmpty(),"edited Current gets separate processing preparation");waitPreparation(audio);
    Capture expected(2,48017),played(2,48017);JobControl measured;
    expect(OfflineAnalysisService(&expected).analyze(r.context,r.domain,measured).value.has_value(),"independent Current reference measurement available");
    audio.seek(0,r,map);audio.setCaptureObserver(&played);audio.play(r,map);awaitPlaying(audio,*d);drain(audio,*d,997);
    expect(identical(expected,played),"prepared prefix and continuation exactly match Current global/local gains and transitions across the prefix boundary");
    expect(audio.preparedCopy()==bytes,"processing changes never require new bytes in the same runtime");
    audio.setCaptureObserver(nullptr);audio.seek(0,r,map);audio.prepare(r);waitPreparation(audio);
    auto stale=r;++r.context.ticket.sessionRevision;r.context.processing.staticGainDb=-2;audio.synchronize(r,map);
    expect(!audio.playbackReady()&&audio.play(stale,map).isNotEmpty(),"manual processing revision rejects stale prepared request");
    audio.prepare(r);audio.play(r,map);auto remap=map;remap.fileToDevice[0]=3;audio.synchronize(r,remap);waitPreparation(audio);
    expect(!audio.active()&&!d->isPlaying(),"routing change revokes pending Play before completion");
    audio.prepare(r);audio.play(r,remap);audio.disconnect();waitPreparation(audio);
    expect(!audio.connected()&&!audio.active()&&!d->isPlaying(),"Disconnect cancels explicit Play without automatic reconnection");
    audio.prepare(r);waitPreparation(audio);auto oldCopy=audio.sharedPreparedCopy();r.context.sessionId="different-document-same-audio";audio.synchronize(r,remap);
    expect(!audio.playbackReady()&&!audio.preparedCopy(),"new runtime cannot inherit identical audio trust or a warm reader");
    audio.prepare(r);waitPreparation(audio);expect(audio.preparedCopy()&&audio.preparedCopy()!=oldCopy.get(),"standalone new runtime establishes its own verified source");
    audio.prepare(r);audio.stop();waitPreparation(audio);
    expect(!d->isPlaying(),"preparation cancellation never authorizes sound");
    auto owner=std::make_shared<SessionVerifiedAssets>();r.context.verifiedAssets=owner;
    audio.synchronize(r,remap);audio.prepare(r);waitPreparation(audio);
    expect(audio.playbackReady()&&owner->counts().copies==1,"session-bound preparation establishes one shared verified asset");
    const auto held=audio.sharedPreparedCopy();owner=std::make_shared<SessionVerifiedAssets>();r.context.verifiedAssets=owner;
    expect(audio.prepare(r).isNotEmpty(),"same audio identity from another owner cannot bypass the synchronized binding");
    audio.synchronize(r,remap);expect(!audio.playbackReady()&&!audio.preparedCopy(),"replaced store inside unchanged runtime revokes prior prepared trust");
    audio.prepare(r);waitPreparation(audio);expect(owner->counts().copies==1&&audio.preparedCopy()!=held.get(),"replacement session store establishes fresh trust for identical bytes");
    owner.reset();audio.synchronize(r,remap);
    expect(r.validationError().isNotEmpty()&&!audio.playbackReady()&&audio.prepare(r).isNotEmpty(),"expired session store cannot authorize a retained warm reader");
}

void playbackNativeChecks(const juce::File& root,int channels)
{
    auto r=fixture(root.getChildFile("native-"+juce::String(channels)+".wav"),48000,channels,240013);
    auto device=std::make_unique<Device>();auto* d=device.get();
    auto ui=std::make_unique<StandaloneMainComponent>(emptyOutputCatalog(),noOutputSettings(),std::move(device));
    Access::stopTimer(*ui);ui->setSize(1320,1100);ui->addToDesktop(juce::ComponentPeer::windowHasTitleBar);ui->setVisible(true);
    auto& s=Access::model(*ui);auto& audio=Access::audio(*ui);
    expect(s.setValidatedSource(r.context.sourceFile,r.context.media),"native programme imported");s.sourceDocument.layout=r.context.fileLayout;
    Access::refresh(*ui);Access::measureSource(*ui);
    const auto deadline=juce::Time::getMillisecondCounterHiRes()+15000;
    while(s.job.isActive()&&juce::Time::getMillisecondCounterHiRes()<deadline){Access::tick(*ui);juce::Thread::sleep(1);}
    expect(s.source.hasMeasurement&&Access::overviewComplete(*ui),"automatic preparation waits for completed initial Source measurement and waveform");
    expect(audio.backgroundPreparing()&&d->opens==0&&!audio.active(),"main schedules background preparation without setup or an output");
    expect(Access::stage(*ui).containsIgnoreCase("playback")&&Access::cancellationAvailable(*ui)&&Access::playAvailable(*ui),"main prominently exposes preparation, cancellation and explicit Play");
    {auto stream=root.getChildFile("playback-preparing-"+juce::String(channels)+".png").createOutputStream();expect(stream&&juce::PNGImageFormat().writeImageToStream(ui->createComponentSnapshot(ui->getLocalBounds()),*stream),"save native background preparation capture");}
    const auto original=s.source.analysis;const auto current=s.currentState.analysis;
    Access::cancel(*ui);waitPreparation(audio);Access::tick(*ui);
    expect(!audio.preparing()&&!audio.playbackReady()&&d->opens==0,"cancelled automatic preparation does not silently restart");
    expect(s.source.analysis==original&&s.currentState.analysis==current,"playback cancellation cannot change measurement authority");
    s.verifiedAssets=std::make_shared<SessionVerifiedAssets>();Access::tick(*ui);
    expect(audio.backgroundPreparing()&&d->opens==0,"relink owner change re-arms background preparation even when audio identity and position are unchanged");
    Access::cancel(*ui);waitPreparation(audio);Access::tick(*ui);
    expect(!audio.preparing(),"cancellation remains suppressed for the same replacement owner");
    expect(audio.connect(48000).isEmpty(),"native test explicitly opens injected output");
    for(int c=0;c<channels;++c)expect(s.setOutputRoute(c,c),"native test explicitly assigns output");
    Access::refresh(*ui);Access::play(*ui);
    expect(audio.preparing()&&audio.active()&&!d->isPlaying(),"main Play retries cancelled preparation with explicit pending intent");
    Access::cancel(*ui);waitPreparation(audio);Access::tick(*ui);
    expect(!d->isPlaying()&&!audio.active(),"main cancellation also revokes explicit pending Play");
    expect(s.selectAuditionDomain(MeasurementDomain::currentState)&&s.setCurrentGainDb(-3),"manual edit is available while background preparation is stopped");
    Access::refresh(*ui);
    Access::tick(*ui);expect(audio.backgroundPreparing()&&s.currentStateNeedsRemeasurement(),"edited Current prepares current processing without measuring or clearing unverified status");
    Access::play(*ui);awaitPlaying(audio,*d);
    juce::AudioBuffer<float> block(maxChannels,256);d->tick(block);expect(block.getMagnitude(0,256)>0,"explicit main Play delivers edited Current to injected output");audio.stop();
    const auto identity=s.auditionRequest().identity();expect(TimelineTransport::seek(12000,identity,s,audio).isEmpty(),"timeline seeking remains stopped and changes position");
    Access::tick(*ui);expect(audio.backgroundPreparing()&&!audio.active(),"new stopped position warms without sound");
    const auto before=audio.positionFrames();click(Access::timeline(*ui),"Zoom +");Access::refreshTimeline(*ui);
    expect(audio.positionFrames()==before&&!d->isPlaying(),"inspection zoom cannot start or move audition");
    waitPreparation(audio);Access::tick(*ui);
    expect(audio.playbackReady()&&!audio.active()&&s.currentStateNeedsRemeasurement(),"ready does not imply Current verified or playing");
    {auto stream=root.getChildFile("playback-ready-"+juce::String(channels)+".png").createOutputStream();expect(stream&&juce::PNGImageFormat().writeImageToStream(ui->createComponentSnapshot(ui->getLocalBounds()),*stream),"save native playback readiness capture");}
    const auto saved=root.getChildFile("playback-"+juce::String(channels)+".levelscope");JobControl save;
    expect(sessionfile::save(sessionfile::capture(s),saved,save).succeeded(),"save programme while prepared and stopped");
    auto previousStore=s.verifiedAssets;const auto oldRuntime=s.sessionId;
    audio.seek(15000,s.auditionRequest(),s.sourceDocument.routing);Access::tick(*ui);
    expect(audio.backgroundPreparing(),"native replacement test has an outstanding prepared worker");
    Access::open(*ui,saved);
    const auto openedUntil=juce::Time::getMillisecondCounterHiRes()+15000;
    do{Access::tick(*ui);juce::Thread::sleep(1);}while((Access::filePending(*ui)||s.persistenceLoading)&&juce::Time::getMillisecondCounterHiRes()<openedUntil);
    expect(s.sessionId!=oldRuntime&&s.verifiedAssets!=previousStore&&s.currentStateDocument.staticGainDb==-3,"native Open installs a new runtime/store while preserving saved processing");
    expect(!d->isPlaying()&&!audio.connected(),"Open during background preparation stays disconnected with reconnection off");
    expect(!audio.active()&&s.currentStateNeedsRemeasurement(),"reopening cannot turn playback readiness into Current verification");
    waitPreparation(audio);Access::tick(*ui);
    // Reopened documents may require overview reconstruction, but can prepare
    // source playback without a device or a setup-window action.
    expect(audio.playbackReady()||audio.backgroundPreparing(),"reopened measured Source automatically prepares in its own runtime");
    s.savedDocumentRevision=s.documentRevision;Access::refresh(*ui);
    expect(ui->performCommand(StandaloneMainComponent::newDocument),"native New action accepts stopped background/readied playback");
    waitPreparation(audio);Access::tick(*ui);
    expect(!s.sourceDocument.hasSource&&!audio.active()&&!audio.playbackReady()&&!d->isPlaying(),"New discards prior preparation and delayed intent without sound");
    expect(s.setValidatedSource(r.context.sourceFile,r.context.media),"native Clear fixture reimports programme");s.sourceDocument.layout=r.context.fileLayout;
    Access::refresh(*ui);expect(audio.prepare(s.auditionRequest()).isEmpty()&&audio.backgroundPreparing(),"native Clear has a background worker to cancel");
    Access::clearSource(*ui);
    expect(!s.sourceDocument.hasSource&&!audio.active()&&!d->isPlaying(),"Clear Source cancels preparation and clears the programme immediately");
    waitPreparation(audio);Access::tick(*ui);
    expect(!s.sourceDocument.hasSource&&!audio.playbackReady()&&!d->isPlaying(),"cancelled worker cannot restore the cleared programme or authorize output");
}

void playbackPreparationChecks(const juce::File& root)
{
    playbackEngineChecks(root);
    for(int channels:{2,6})playbackNativeChecks(root,channels);
}

#pragma once
namespace
{
void compositionSafety(const juce::File& root)
{
    auto f=fixtureComposition(root,48000,CompositionPreset::surroundMonoCentre,true,48000*4);
    auto s=compositionModel(f.recipe);measuredComposition(s,MeasurementDomain::source);measuredComposition(s,MeasurementDomain::currentState);
    auto original=f.recipe;
    for(int fault=0;fault<8;++fault)
    {
        auto invalid=original;
        if(fault==0)invalid.confirmation.clear();if(fault==1)invalid.dialogue.offsetFrames=1;
        if(fault==2)invalid.dialogue.media.frames--;if(fault==3)invalid.dialogue.media.sampleRate=44100;
        if(fault==4)invalid.dialogue.purpose="analysis-reference";if(fault==5)invalid.dialogue.layout.roles[0]=juce::AudioChannelSet::left;
        if(fault==6)invalid.preset=static_cast<CompositionPreset>(99);if(fault==7)invalid.me.layout.roles[0]=juce::AudioChannelSet::right;
        expect(invalid.validate().isNotEmpty()&&!s.setValidatedComposition(invalid),"invalid mapping purpose alignment or identity cannot replace valid session");
    }
    auto context=s.measurementContext(MeasurementDomain::source);auto mixed=context;mixed.assembly=original.me.assembly;
    expect(!mixed.isValid(),"cannot combine legacy assembly and composition representations");
    auto labelled=original;labelled.dialogue.language="Italian";labelled.dialogue.versionLabel="review <v2> & mix";
    const auto revision=s.revision;const auto source=s.source.analysis,current=s.currentState.analysis;
    expect(s.setValidatedComposition(labelled)&&s.revision==revision&&s.source.analysis==source&&s.currentState.analysis==current,"labels preserve audio revision and valid measurements");
    expect(s.sourceDocument.composition->confirmation==original.confirmation,"labels preserve exact alignment confirmation");
    const auto id=OfflineAnalysisService::identityFor(s.measurementContext(MeasurementDomain::source),MeasurementDomain::source,s.sourceDocument.media.sha256);
    s.referenceVideo.startOnAudio=video::Time(-1001,30000);s.sourceDocument.routing.fileToDevice[0]=0;s.documentChanged();
    expect(id==OfflineAnalysisService::identityFor(s.measurementContext(MeasurementDomain::source),MeasurementDomain::source,s.sourceDocument.media.sha256),"video offsets and monitor routes preserve composition measurement identity");
    auto draft=proposal::draftGlobalTrim(s.currentTrimRequest());expect(draft.status==proposal::TrimStatus::draft&&s.acceptTrimProposal(draft),"M3 drafts measured composed programme correction");
    const auto beforeGain=s.currentStateDocument.staticGainDb;expect(beforeGain==0&&s.currentState.hasMeasurement,"draft alone never changes audio");
    expect(s.applyTrim()&&s.currentStateDocument.staticGainDb!=beforeGain&&!s.currentState.hasMeasurement,"explicit Apply changes only Current plan and invalidates measured Current");
    const double applied=s.currentStateDocument.staticGainDb;s.measurementUnsuccessful(MeasurementDomain::currentState,ErrorCode::cancelled);
    expect(s.currentStateDocument.staticGainDb==applied&&!s.currentState.hasMeasurement&&s.currentStateDocument.needsRemeasurement,"cancelled verification leaves applied gain unverified without rollback");
    measuredComposition(s,MeasurementDomain::currentState);expect(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"actual composed Current remeasurement passes");
    auto noChange=proposal::draftGlobalTrim(s.currentTrimRequest());const auto passing=s.currentState.analysis;const auto audioRevision=s.revision;
    expect(noChange.status==proposal::TrimStatus::noChange&&s.acceptTrimProposal(noChange)&&s.applyTrim()
        &&s.revision==audioRevision&&s.currentState.analysis==passing&&s.currentStateDocument.staticGainDb==applied,"already-passes preserves measured Current and exact gain");
    JobControl job;auto doc=sf::capture(s);auto alternate=root.getChildFile("dialogue-identical.wav");
    expect(original.dialogue.file.copyFileTo(alternate)&&sf::relinkAudioMember(doc,"dialogue/file",alternate),"composition member relink changes locator only");
    auto reopened=sf::verify(doc,job);expect(reopened.document&&reopened.document->audioResolved&&reopened.document->source.composition->confirmation==original.confirmation,"identical verified relink retains composition confirmation");
    auto member=original.me.assembly->members.front();auto hidden=member.file.getSiblingFile("temporarily-missing.wav");expect(member.file.moveFileTo(hidden),"move disposable required member for missing test");
    auto missing=sf::verify(sf::capture(s),job);expect(missing.document&&!missing.document->audioResolved&&missing.prepared->model.sourceContentInvalid
        &&!missing.prepared->model.source.hasMeasurement&&missing.prepared->model.currentStateDocument.staticGainDb==applied,"missing one component member recovers settings but cannot retain a current pass");
    expect(!prepareVerifiedSource(s.measurementContext(MeasurementDomain::source),job).value,"missing required member cannot become substitute silence");
    expect(hidden.moveFileTo(member.file),"restore only disposable test member");
    auto bad=doc;expect(sf::relinkAudioMember(bad,"dialogue/file",f.oracle),"different-content replacement staged for validation");
    auto rejected=sf::verify(bad,job);expect(rejected.document&&!rejected.document->audioResolved,"matching path selection cannot authorize different content");
    auto changed=original;
    // Change one exact decoded sample, with unchanged rate and frame count.
    juce::AudioBuffer<float> altered(1,48000*4);altered.clear();altered.setSample(0,100,0.125f);const auto changedFile=root.getChildFile("changed-dialogue.wav");floatFile(changedFile,48000,altered);
    auto changedMedia=validateMedia(changedFile,job);expect(changedMedia.succeeded(),"changed same-length dialogue inspected");changed.dialogue.file=changedFile;changed.dialogue.media=changedMedia.media;
    expect(changed.validate().isNotEmpty()&&changed.definitionIdentity()!=original.definitionIdentity(),"same rate/count different content requires new confirmation");
    expect(changed.confirmComplementaryBalancedAlignment()&&s.setValidatedComposition(changed)&&s.currentStateDocument.staticGainDb==applied
        &&!s.source.hasMeasurement&&!s.currentState.hasMeasurement,"explicit replacement invalidates both measured domains and retains applied gain");
    auto stale=*source;expect(!s.acceptMeasurement(stale),"old composition completion cannot become current");
    auto mismatch=trimRequest(source.get(),passing.get(),s.measurementContext(MeasurementDomain::source),s.measurementContext(MeasurementDomain::currentState),s.getSelectedTargetProfile(),s.programmeApplicability);
    expect(proposal::draftGlobalTrim(mismatch).status==proposal::TrimStatus::unavailable,"individually valid old measurements cannot satisfy a new composition");
    SnapshotStorage noSpace;noSpace.directory=root;noSpace.availableBytes=[](const auto&){return juce::int64(0);};
    expect(prepareVerifiedSource(context,job,noSpace).error==ErrorCode::insufficientDiskSpace,"composition temporary space failure explicit");
    struct BrokenWriter final:SnapshotWriter{bool write(const void*,size_t)override{return false;}bool flush()override{return false;}};
    SnapshotStorage broken;broken.directory=root;broken.createWriter=[](const auto&){return std::make_unique<BrokenWriter>();};
    expect(prepareVerifiedSource(context,job,broken).error==ErrorCode::writeFailed,"component copy write failure fails entire composition");
    JobControl cancelled;cancelled.cancel=true;expect(prepareVerifiedSource(context,cancelled).error==ErrorCode::cancelled,"cancelled composition preparation publishes no partial reader");
    expect(root.findChildFiles(juce::File::findFiles,false,"levelscope-*.wav").isEmpty(),"owned failed source copies removed");
    // A saved snapshot is coherent even after an edit or a runtime replacement.
    auto saved=sf::capture(s);const auto target=root.getChildFile("coherent.levelscope");expect(s.setCurrentGainDb(applied-1),"edit after save capture");
    expect(sf::save(saved,target,job).succeeded()&&sf::saved(s,saved,target)&&s.documentDirty(),"saving captured revision leaves later edit unsaved");
    StandaloneSessionModel other;expect(!sf::saved(other,saved,target),"late save completion cannot mark a newly opened runtime saved");
    auto recovered=sf::read(target,job);expect(recovered.document&&recovered.document->current.staticGainDb==applied,"saved core keeps its captured processing precision");
    auto html=ex::reportHtml(ex::report(ex::capture(s)));expect(html.contains("Source: original confirmed")&&!html.contains("<v2>"),"report distinguishes composition and escapes labels");
}
void compositionOverload(const juce::File& root)
{
    expect(root.createDirectory().wasOk(),"overload fixture directory");juce::AudioBuffer<float> a(2,48000);for(int c=0;c<2;++c)for(int n=0;n<48000;++n)a.setSample(c,n,0.75f);
    ProgrammeComposition recipe;recipe.me=component(root,"me",48000,a,false);recipe.dialogue=component(root,"dialogue",48000,a,false);expect(recipe.confirmComplementaryBalancedAlignment(),"known test sum confirmed");
    auto s=compositionModel(recipe);Capture capture(2,48000);const auto measured=measuredComposition(s,MeasurementDomain::source,&capture);
    expect(capture.data.getSample(0,0)==1.5f&&measured.detail.samplePeak.value&&*measured.detail.samplePeak.value>0,"finite over-full-scale sum remains measurable without hidden attenuation");
    auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));auto request=s.auditionRequest();auto routing=routes(2);
    audition.synchronize(request,routing);expect(audition.connect(48000).isEmpty(),"overload uses injected device");audition.setMonitorGainDb(0);
    // 1024 frames lets the existing smooth startup monitor ramp reach overload.
    expect(audition.play(request,routing,true).isEmpty(),"explicit composed Source playback");awaitPlaying(audition,*injected);juce::AudioBuffer<float> output(maxChannels,1024);injected->tick(output);audition.poll();
    expect(audition.getState()==AuditionState::failed&&output.getMagnitude(0,1024)==0&&audition.positionFrames()==0,"sum overload silences entire callback and preserves presented audio position");
    const auto* copy=audition.preparedCopy();audition.setMonitorGainDb(-15);audition.poll();expect(!audition.active(),"monitor reduction never automatically resumes");
    Capture recovered(2,48000);audition.setCaptureObserver(&recovered);expect(audition.play(request,routing,true).isEmpty(),"explicit lower-monitor restart");awaitPlaying(audition,*injected);drain(audition,*injected,256);
    expect(audition.preparedCopy()==copy&&recovered.data.getSample(0,0)==1.5f,"overload recovery reuses copies and preserves original processing samples");audition.disconnect();
    measuredComposition(s,MeasurementDomain::currentState);JobControl job;expect(!ex::render(ex::capture(s),root,job).value,"overloaded failing Current cannot bypass M10 publication gates");
}
}

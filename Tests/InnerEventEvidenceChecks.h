// Included inside the event planning test namespace.
void inspectInnerMusic(const juce::File& root,const juce::File& file)
{
    JobControl job;
    auto media=validateMedia(file,job);
    expect(media.succeeded()&&media.media.sha256=="7ba6404041c94c909df290a6829a0fe92ce594329a1d575490c2bb4a5b6cc8f0","unchanged music development PCM");
    StandaloneSessionModel session;
    expect(session.setValidatedSource(file,media.media),"set verified source");
    session.sourceDocument.layout=layoutFor(2);session.setSelectedTargetProfileIndex(1);
    session.setProgrammeApplicability(Applicability::preProducedLinearProgramme);
    auto measured=measuredComposition(session,MeasurementDomain::source);
    session.synchronizeSections();
    auto analysis=analyzeProgramme(session.currentState.analysis->detail,session.currentState.analysis->identity);
    auto passages=findLraPassages(analysis,18);
    const auto context=session.measurementContext(MeasurementDomain::currentState);
    Plan before;before.frames=media.media.frames;before.sampleRate=int(media.media.sampleRate);
    before.inputIdentity=sectionInputIdentity(context);before.intervals={{0,before.frames,0}};
    auto settings=RegionPlanningSettings::innerDefaults(RegionMode::automatic,before.inputIdentity);
    settings.allowAutomaticMargins=true;settings.leftMarginMs=settings.rightMarginMs=1000;
    auto constraints=resolveRegionConstraints(settings,session.sectionDocument.map(),*passages,before);
    auto coarse=formRegions(constraints,*analysis,18);
    auto roles=OfflineAnalysisService::measurementRoles(context.fileLayout);
    auto planId=regionPlanIdentity(before,0);
    EventFeatureCollector features(constraints,coarse,*analysis,roles,18,"evidence-stage",measured.identity,planId,&measured.detail);
    PredictionCollector bins(before.sampleRate,roles,before.frames,"evidence-stage",analysis->snapshotIdentity,planId);
    juce::AudioFormatManager formats;formats.registerBasicFormats();
    std::unique_ptr<juce::AudioFormatReader> reader(formats.createReaderFor(file));
    expect(bool(reader)&&root.createDirectory().wasOk(),"owned evidence directory and read-only source");
    juce::AudioBuffer<float> audio(2,4096);
    for(int64_t pos=0;pos<before.frames;)
    {
        const int count=int(std::min<int64_t>(4096,before.frames-pos));
        if(!reader->read(audio.getArrayOfWritePointers(),2,pos,count))throw std::runtime_error("Incomplete development input");
        features.push(audio.getArrayOfReadPointers(),audio.getArrayOfReadPointers(),count);
        bins.push(audio.getArrayOfReadPointers(),count);pos+=count;
    }
    auto sidecar=bins.finish();auto evidence=features.finish(sidecar.get());
    auto ledger=evidence->innerDescription();
    auto formed=formEventRegions(constraints,*analysis,18,*evidence);
    for(const auto& geoms:formed.alternatives)for(const auto& g:geoms)
        ledger+="Geometry root "+juce::String(g.root+1)+" ["+juce::String(g.support.begin)+", "+juce::String(g.support.end)+") "+g.reason+"\n";
    expect(root.getChildFile("evidence.txt").replaceWithText(ledger),"retain frozen witness and formation ledger before assertions");
    expect(evidence->conflictingEntrances>0&&evidence->omittedWitnesses==0,"record occupied contenders independently from summary-cap losses");
    int found=0;
    for(const auto& g:evidence->groups)if(g.innerNominee&&(g.root==2||g.root==0))
    {
        ++found;const bool first=g.root==2;
        expect(g.onset.begin<(first?12.283:52.615)*48000&&g.onset.end>(first?12.28:52.4)*48000,"real pre-hit inner onset");
        const auto& options=formed.alternatives[g.root];
        auto alternative=std::find_if(options.begin(),options.end(),[](const auto& x){return x.family==BoundaryFamily::innerAfter;});
        expect(alternative!=options.end(),"legal inner after-return geometry reaches formation");
        auto decision=assignment(alternative->support);decision.geometry=*alternative;ProposalPolicy policy;
        const auto plan=constructRegionPlan(before,0,constraints,policy,{decision},0);
        expect(plan.error.isEmpty()&&plan.envelope.gainAt(g.onset.begin)==-2&&plan.envelope.gainAt(g.lastAttack.end-1)==-2
            &&plan.envelope.gainAt(g.decay.end-1)==-2,"canonical fixed ramps reach the onset and hold through the last hit and full decay");
        expect(plan.envelope.gainAt(alternative->support.begin)==0&&plan.envelope.gainAt(alternative->support.end-1)==0,"actual two exterior ramp endpoints");
        auto automated=before;
        const auto mid=g.onset.end+(g.lastAttack.end-g.onset.end)/2;
        automated.intervals={{0,mid,-.75},{mid,before.frames,-1.25}};
        const auto overExisting=constructRegionPlan(automated,0,constraints,policy,{decision},0);
        expect(overExisting.error.isEmpty(),"existing automation accepts the local overlay");
        for(auto frame:{g.onset.begin,mid-1,mid,g.lastAttack.end-1,g.decay.end-1})
            expect(automated.gainAt(frame)-overExisting.envelope.gainAt(frame)==2,"canonical additional gain maintained through existing automation and both event edges");
        auto earlier=std::find_if(options.begin(),options.end(),[](const auto& x){return x.family==BoundaryFamily::innerDecay;});
        expect(earlier!=options.end()&&earlier->support.end<alternative->support.end,"distinct earlier return survives formation");
        decision.geometry=*earlier;const auto shortPlan=constructRegionPlan(before,0,constraints,policy,{decision},0);
        expect(shortPlan.error.isEmpty()&&shortPlan.envelope.gainAt(g.lastAttack.end-1)==-2&&shortPlan.envelope.gainAt(g.inDecay.end-1)==-2
            &&shortPlan.envelope.gainAt(g.decay.end-1)>-2&&!assessEventCoverage(before,shortPlan,constraints,*evidence,*analysis,18).complete(),"earlier fixed ramp covers the final hit and exposes the independently recorded later ending");
        const auto coverage=assessEventCoverage(before,plan,constraints,*evidence,*analysis,18);
        expect(!coverage.complete(),"parent uncertainty or later ending survives inner nomination");
        expect(g.lastAttack.end>=(first?30.25:57.3)*48000&&g.endingSupported,"last comparable hit and independent sealed ending");
    }
    expect(found==2,"both need-bearing inner passages retained");
}

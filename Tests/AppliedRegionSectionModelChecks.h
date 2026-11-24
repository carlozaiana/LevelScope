// Included by the section processing runner; uses its real measurement service and PCM oracle.
namespace {
void checkProposeInitialization(StandaloneSessionModel s)
{
    juce::String error;
    s.clearCurrentState();
    const auto source = s.source.analysis;
    expect(source && !source->restoredFromSession && !s.currentStateDocument.isInitialized, "cleared Current with fresh Source");
    auto blocked = s;
    blocked.auditionBusy = []{ return true; };
    expect(!blocked.ensureCurrentForPropose(error) && !blocked.currentStateDocument.isInitialized, "playing prevents initialization");
    blocked = s;
    expect(blocked.job.begin(JobStage::analysis,blocked.revision) && !blocked.ensureCurrentForPropose(error), "active job prevents initialization");
    blocked = s; blocked.sourceContentInvalid = true;
    expect(!blocked.ensureCurrentForPropose(error) && !blocked.currentStateDocument.isInitialized, "changed source cannot initialize");
    blocked = s; blocked.sourceDocument.layout.confirmed = false;
    expect(!blocked.ensureCurrentForPropose(error) && !blocked.currentStateDocument.isInitialized, "unconfirmed roles cannot initialize");
    StandaloneSessionModel empty;
    expect(!empty.ensureCurrentForPropose(error) && !empty.currentStateDocument.isInitialized, "missing source cannot initialize");
    expect(s.ensureCurrentForPropose(error) && error.isEmpty() && s.source.analysis == source && s.currentState.hasMeasurement
           && s.currentState.analysis->identicalSourceMeasurementIdentity == source->identity
           && s.currentState.analysis->domain == MeasurementDomain::currentState, "absent Current reuses proven identical live Source at unity");
    const auto unity = s.currentState.analysis;
    const auto revision = s.revision;
    expect(s.ensureCurrentForPropose(error) && s.revision == revision && s.currentState.analysis == unity, "existing measured Current stays exact");
    expect(s.setCurrentGainDb(-1.25) && !s.currentState.hasMeasurement, "manual gain leaves Current unmeasured");
    const auto before = s.manualState();
    const auto edited = s.revision;
    expect(s.ensureCurrentForPropose(error) && s.manualState()==before && s.revision==edited && !s.currentState.hasMeasurement,
           "existing edited and unverified Current cannot be reset or assigned Source measurements");

    s.clearCurrentState();
    s.source = MeasurementSummary();
    expect(s.ensureCurrentForPropose(error) && s.awaitsInitialCurrent() && !s.currentState.hasMeasurement,
           "valid absent Current initializes while awaiting a fresh Source measurement");
    measuredComposition(s,MeasurementDomain::source);
    expect(!s.awaitsInitialCurrent() && s.currentState.hasMeasurement &&
           s.currentState.analysis->identicalSourceMeasurementIdentity==s.source.analysis->identity, "newly accepted Source verifies identical unity Current");
}

void appliedRegionSectionChecks(const juce::File& root)
{
    constexpr int rate=48000;
    juce::String error;
    juce::AudioBuffer<float> raw;
    auto s=inputModel(root,rate,2,0,raw);
    measuredComposition(s,MeasurementDomain::source);
    const auto source=s.source.analysis;
    checkProposeInitialization(s);
    s.synchronizeSections();
    const auto noLocal=s.manualState();
    expect(!s.createSectionsFromAppliedRegions(error) && s.manualState()==noLocal && !error.isEmpty(), "no local gain yields a truthful, inert result");
    auto requirements=proposal::JobRequirements::editable(s.getSelectedTargetProfile());
    requirements.integratedTenths=juce::roundToInt(*source->detail.integrated.value*10)+10;
    expect(s.setJobRequirements(requirements,error), "set an achievable one-LU level change");
    auto policy=s.proposalPolicy;policy.durationThousandths=1000;
    expect(s.setProposalPolicy(policy,error), "allow a global correction to record real applied history");
    if (!s.currentState.hasMeasurement) measuredComposition(s,MeasurementDomain::currentState);
    expect(s.acceptTrimProposal(proposal::draftGlobalTrim(s.currentTrimRequest())) && s.applyTrim() && s.lastAppliedExplanation,
           "install a real scalar proposal and its history");
    const auto history=encodeAppliedExplanation(*s.lastAppliedExplanation);
    auto& doc=s.sectionDocument;
    expect(doc.rename(1,"Opening",error) && doc.split(1,48006,error) && doc.rename(2,"Named phrase",error)
           && doc.split(2,96000,error) && doc.rename(3,"Ending",error), "retain named original sections, including one cut inside a ramp");
    processing::SectionGainPlan plan;
    plan.sampleRate=rate;plan.frames=raw.getNumSamples();plan.inputIdentity=doc.map().inputIdentity;
    plan.intervals={{0,47000,0},{47000,96000,-3},{96000,144000,0},{144000,plan.frames,-1}};
    auto next=s.manualState();next.gains=plan;
    expect(s.commitManualEdit(next,"Local treatment"), "apply local treatment");
    const auto token=s.beginAppliedBoundary(1);
    expect(token && s.previewAppliedBoundary(token,48000) && !s.createSectionsFromAppliedRegions(error)
           && s.commitAppliedBoundary(token), "snapshot rejects a pending drag and preserves its subsequently committed cut");
    const auto before=s.manualState();
    auto blocked=s;
    expect(blocked.job.begin(JobStage::analysis,blocked.revision) && !blocked.createSectionsFromAppliedRegions(error)
           && blocked.manualState()==before, "active job blocks snapshot without mutation");
    blocked=s;blocked.auditionBusy=[]{return true;};
    expect(!blocked.createSectionsFromAppliedRegions(error) && blocked.manualState()==before, "playback blocks snapshot");
    blocked=s;blocked.currentStateDocument.sectionGains->inputIdentity=juce::String::repeatedString("0",64);
    expect(!blocked.createSectionsFromAppliedRegions(error) && blocked.sectionDocument.map().digest()==before.map.digest(),
           "cut list with stale input identity rejected before changing annotations");
    blocked=s;blocked.currentStateDocument.sectionGains->intervals.back().end--;
    expect(!blocked.createSectionsFromAppliedRegions(error) && blocked.sectionDocument.map().digest()==before.map.digest(),
           "invalid final cut rejects the complete transaction");
    auto unverified=s;unverified.manualHistory.reset();
    expect(unverified.createSectionsFromAppliedRegions(error) && !unverified.currentState.hasMeasurement &&
           unverified.currentStateNeedsRemeasurement(), "creation cannot verify unmeasured processing");
    Capture capture(2,raw.getNumSamples());
    measuredComposition(s,MeasurementDomain::currentState,&capture);
    const auto measured=s.currentState.analysis;
    requirements.integratedTenths=juce::roundToInt(*measured->detail.integrated.value*10);
    expect(s.setJobRequirements(requirements,error) && s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,
           "independently measured nonunity processing meets selected limits");
    const auto context=s.measurementContext(MeasurementDomain::currentState);
    const auto verifyAudio=[&] {
        expect(s.source.analysis==source && s.currentState.analysis==measured && !s.currentStateNeedsRemeasurement()
               && s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass, "verified snapshots and PASS survive annotation transaction");
        const auto current=s.measurementContext(MeasurementDomain::currentState);
        expect(before.globalDb==current.processing.staticGainDb && before.gains->audioKey()==current.processing.sectionGains->audioKey()
               && current.processing.sectionGains->radius(1)==2400 && current.ticket.sessionRevision==context.ticket.sessionRevision,
               "exact global, canonical envelope, 50 ms transition and audio revision preserved");
        juce::AudioBuffer<float> processed(raw);
        ProgrammeGainProcessor processor(current.processing);
        expect(processor.process(processed.getArrayOfWritePointers(),2,raw.getNumSamples(),0)
               && sameSamples(processed,capture.data) && sameSamples(processed,oracleAudio(raw,context.processing)),
               "all processed float samples match prior real Current and independent oracle, including the close annotation cut");
        expect(encodeAppliedExplanation(*s.lastAppliedExplanation)==history, "full applied explanation remains historical evidence");
    };
    verifyAudio();
    const auto baseline=s.sectionProposalRequest();
    JobControl job;
    const auto candidate=generateSectionProposal(baseline,job);
    expect(bool(candidate.value) && s.acceptSectionProposal(*candidate.value), "obtain a genuine bound pre-snapshot result");
    expect(!s.createSectionsFromAppliedRegions(error) && s.manualState()==before, "outstanding unapplied proposals require resolution");
    s.discardSectionProposal();
    const auto count=s.manualHistory->undoCount();
    expect(s.createSectionsFromAppliedRegions(error) && s.manualHistory->undoCount()==count+1, "single recoverable transaction inserts both missing cuts");
    const auto after=s.manualState();
    const auto& rows=after.map.sections;
    expect(rows.size()==5 && rows[0].begin==0 && rows[0].end==48000 && rows[1].end==48006 && rows[2].end==96000 &&
           rows[3].end==144000 && rows[4].end==plan.frames, "exact sorted union, manually moved cut, existing cut deduplicated");
    expect(rows[0].id==1 && rows[0].name=="Opening" && rows[2].id==2 && rows[2].name=="Named phrase" &&
           rows[3].id==3 && rows[3].name=="Ending" && rows[1].id==4 && rows[4].id==5 && after.map.nextId==6,
           "existing names and identities survive and only new subdivisions consume IDs");
    verifyAudio();
    const auto serial=s.documentRevision, generation=s.sectionPlanningGeneration;
    expect(s.createSectionsFromAppliedRegions(error) && error.contains("already") && s.manualState()==after &&
           s.documentRevision==serial && s.sectionPlanningGeneration==generation && s.manualHistory->undoCount()==count+1,
           "repeated action makes no edit, history entry or invalidation");
    expect(!s.acceptSectionProposal(*candidate.value) && !s.proposalOperationMatches(baseline,before,context.ticket.sessionRevision,s.jobRequirements.digest()),
           "previous proposal authority unavailable after snapshot");
    expect(s.undoManual() && s.manualState()==before && !s.acceptSectionProposal(*candidate.value) &&
           !s.proposalOperationMatches(baseline,before,context.ticket.sessionRevision,s.jobRequirements.digest()),
           "one Undo restores entire original map and cannot revive the old result after geometry matches");
    verifyAudio();
    expect(s.undoManual(true) && s.manualState()==after && !s.acceptSectionProposal(*candidate.value), "one Redo restores snapshot only");
    verifyAudio();

    const auto saved=root.getChildFile("sections-from-regions.levelscope");
    expect(sf::save(sf::capture(s),saved,job).succeeded(),"save snapshot");
    auto disk=sf::read(saved,job);expect(bool(disk.document),"read snapshot");
    auto opened=sf::verify(*disk.document,job);
    expect(opened.prepared && opened.prepared->model.manualState()==after &&
           opened.prepared->model.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass &&
           encodeAppliedExplanation(*opened.prepared->model.lastAppliedExplanation)==history,"reopened sections, original processing, historical explanation and verified PASS");
    auto& restored=opened.prepared->model;
    restored.clearCurrentState();
    expect(restored.source.analysis->restoredFromSession && restored.ensureCurrentForPropose(error) && !restored.currentState.hasMeasurement,
           "restored Source alone cannot manufacture a fresh identical-Source receipt");
    expect(s.setManualSectionGains({2},-2.0) && !s.currentState.hasMeasurement && s.currentStateNeedsRemeasurement() &&
           s.currentStateDocument.staticGainDb==before.globalDb && s.currentStateDocument.sectionGains->gainAt(70000)==-2.0 &&
           s.currentStateDocument.sectionGains->gainAt(160000)==-1.0, "deliberate section level edit sets absolute local gain and keeps unrelated processing");

    for (const auto capacity : {true,false})
    {
        auto map=before.map;
        if (capacity)
        {
            map.sections.clear();
            for(int i=0;i<10000;++i)map.sections.push_back({uint64_t(i+1),i,i==9999?plan.frames:i+1,"Existing"});
            map.nextId=10001;
        }
        else map.nextId=std::numeric_limits<uint64_t>::max()-1;
        expect(s.sectionDocument.restore(map,map.inputIdentity),"install valid near-limit annotation map");
        const auto initial=s.manualState();const auto journal=s.manualHistory->undoCount();
        expect(s.appliedRegionSectionsBlocked().isNotEmpty() && !s.createSectionsFromAppliedRegions(error) &&
               s.manualState()==initial && s.manualHistory->undoCount()==journal,"map/ID limit rejects entire command without partial insertion");
    }
    s.sourceContentInvalid=true;
    expect(!s.createSectionsFromAppliedRegions(error) && error.contains("source"), "invalidated source blocks edits");
}
}

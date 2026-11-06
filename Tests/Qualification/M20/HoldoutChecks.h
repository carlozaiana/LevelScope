// Supplemental engineering harness. The user-approved protocol was committed
// before any measurement of these independent inputs; never a tuning fixture.
juce::var holdoutMetrics(const measurement::Result& d)
{
    auto value=juce::var(new juce::DynamicObject());auto* o=value.getDynamicObject();
    const char* names[]={"integrated_lufs","true_peak_dbtp","lra_lu","maximum_m_lufs","maximum_s_lufs","sample_peak_dbfs"};
    const measurement::Metric* values[]={&d.integrated,&d.truePeak,&d.lra,&d.momentaryMax,&d.shortTermMax,&d.samplePeak};
    for(int n=0;n<6;++n)o->setProperty(names[n],values[n]->value?juce::var(*values[n]->value):juce::var());return value;
}
juce::var holdoutPlan(const ProcessingPlan& p)
{
    auto value=juce::var(new juce::DynamicObject());auto* o=value.getDynamicObject();o->setProperty("identity",proposalPlanIdentity(p));o->setProperty("global_db",p.staticGainDb);
    juce::Array<juce::var> rows;if(p.sectionGains)for(const auto& part:p.sectionGains->intervals)
    {auto row=juce::var(new juce::DynamicObject());auto* r=row.getDynamicObject();r->setProperty("begin_frame",juce::String(part.begin));r->setProperty("end_frame_exclusive",juce::String(part.end));r->setProperty("plateau_db",part.db);rows.add(row);}
    o->setProperty("nominal_intervals",rows);return value;
}
void controlledHoldouts(const juce::File& inputs,const juce::File& root)
{
    struct Input{const char* name;const char* hash;int channels;};
    const Input files[]={
        {"DialogAmbiance.wav","602485eef95774826a575b0069d358e6723a31f6928dbd6e04dd0df44c98eb90",2},
        {"MusicInBackAndForegroundSomeDialog.wav","7ba6404041c94c909df290a6829a0fe92ce594329a1d575490c2bb4a5b6cc8f0",2},
        {"5.1_Busy.wav","4fa30bfd9a8fbb64e68fc00eaf11f8adf81ead8502d05ce639b22cc07645d7b6",6}};
    for(const auto& input:files)
    {
        const auto file=inputs.getChildFile(input.name),directory=root.getChildFile(file.getFileNameWithoutExtension());expect(directory.createDirectory().wasOk(),"owned independent-review output directory");
        JobControl job;auto media=validateMedia(file,job);expect(media.succeeded()&&media.media.sha256==input.hash&&media.media.channels==input.channels&&media.media.sampleRate==48000,"independent input matches the precommitted content and explicit iXML role order");
        StandaloneSessionModel baseline;expect(baseline.setValidatedSource(file,media.media),"adopt independent source read-only");baseline.sourceDocument.layout=layoutFor(input.channels);baseline.setSelectedTargetProfileIndex(1);baseline.setProgrammeApplicability(Applicability::preProducedLinearProgramme);
        auto source=measuredComposition(baseline,MeasurementDomain::source);expect(baseline.currentState.hasMeasurement,"initial identical Current uses completed Source provenance");baseline.synchronizeSections();
        const int ceiling=source.detail.lra.value?int(std::floor((*source.detail.lra.value-1.)*10.)):0;
        if(!source.detail.lraStable||ceiling<=0)
        {
            auto unavailable=juce::var(new juce::DynamicObject());auto* u=unavailable.getDynamicObject();u->setProperty("file",file.getFullPathName());u->setProperty("content_sha256",input.hash);u->setProperty("baseline",holdoutMetrics(source.detail));u->setProperty("outcome","No stable positive one-LU-reduced LRA ceiling; neither method assessed");cases.add(unavailable);directory.getChildFile("not-assessed.json").replaceWithText(juce::JSON::toString(unavailable,true));continue;
        }
        auto requirements=baseline.jobRequirements;requirements.integratedTenths=-230;requirements.integratedToleranceTenths=2;requirements.truePeakTenths=-10;requirements.maximumLraTenths=ceiling;requirements.maximumMomentaryTenths.reset();requirements.maximumShortTermTenths.reset();requirements.provenance="Independent engineering comparison; no delivery-compliance claim";
        juce::String error;expect(baseline.setJobRequirements(requirements,error),"identical predeclared numerical requirements for both methods");
        ProposalPolicy policy;expect(policy.additionalTenths==30&&policy.cumulativeTenths==120&&policy.durationThousandths==200&&policy.reserveTenths==5&&!policy.regionLimit&&baseline.setProposalPolicy(policy,error),"frozen preservation preferences unchanged");
        for(bool experimental:{false,true})
        {
            auto s=baseline;const juce::String method=experimental?"m20":"m19";const auto destination=directory.getChildFile(method);expect(destination.createDirectory().wasOk(),"separate method evidence directory");
            auto settings=RegionPlanningSettings::defaults(RegionMode::automatic,sectionInputIdentity(s.measurementContext(MeasurementDomain::source)));
            expect(s.setRegionPlanning(experimental?std::optional<RegionPlanningSettings>(settings):std::nullopt,error),"explicit method selection, independent unchanged baseline");
            const auto identity=s.currentState.analysis->identity,applied=proposalPlanIdentity(s.currentStateDocument.processingPlan());
            const auto start=juce::Time::getMillisecondCounterHiRes();JobControl proposalJob;auto outcome=generateSectionProposal(s.sectionProposalRequest(),proposalJob);
            expect(outcome.value&&(*outcome.value)->passes()<=3,"complete bounded engineering outcome without extending the private-pass budget");const auto result=*outcome.value;
            expect(s.currentState.analysis->identity==identity&&proposalPlanIdentity(s.currentStateDocument.processingPlan())==applied,"assessment leaves original Current and processing unchanged");expect(s.acceptSectionProposal(result),"actual outcome matches immutable request binding");
            auto record=juce::var(new juce::DynamicObject());auto* row=record.getDynamicObject();row->setProperty("file",file.getFullPathName());row->setProperty("content_sha256",input.hash);row->setProperty("channels",input.channels);row->setProperty("frames",juce::String(media.media.frames));row->setProperty("method",method);row->setProperty("requirements",requirements.encode());row->setProperty("preferences",policy.encode());row->setProperty("settings",experimental?settings.encode():juce::var());row->setProperty("baseline",holdoutMetrics(source.detail));row->setProperty("baseline_identity",identity);row->setProperty("outcome",result->outcomeTitle());row->setProperty("reason",result->outcomeReason());row->setProperty("private_passes",result->passes());row->setProperty("assessment_seconds",(juce::Time::getMillisecondCounterHiRes()-start)/1000.);row->setProperty("apply_available",s.canApplySectionProposal());
            if(result->assessment()){row->setProperty("retained_actual_candidate",holdoutMetrics(result->assessment()->detail));row->setProperty("retained_actual_evaluation",proposalMeasurementRecord(*result->assessment(),s.getSelectedTargetProfile(),s.programmeApplicability).evaluation);}
            if(auto preview=result->previewPlan())row->setProperty("inspection_only_plan",holdoutPlan(*preview));
            if(const auto review=result->regionReview())
            {
                expect(review->proxyEvaluations<=64&&review->planningBytes<=64*1024*1024,"held-out M20 result stays within frozen proxy and additional-memory budgets");row->setProperty("proxy_count",int(review->proxyEvaluations));row->setProperty("planning_seconds",review->planningSeconds);row->setProperty("planning_bytes",int64_t(review->planningBytes));row->setProperty("retained_attempt",review->retainedAttempt<review->attempts.size()?int(review->retainedAttempt):-1);
                juce::Array<juce::var> attempts;
                for(const auto& a:review->attempts)
                {
                    auto value=juce::var(new juce::DynamicObject());auto* v=value.getDynamicObject();v->setProperty("identity",a.plan.identity);v->setProperty("actual_evaluation",evaluationText(a.evaluation));v->setProperty("actual_pass",a.evaluation.verdict==Verdict::pass);juce::Array<juce::var> actual;for(const auto& metric:a.measured)actual.add(metric.value?juce::var(*metric.value):juce::var());v->setProperty("actual_il_tp_lra_m_s_sample_peak",actual);v->setProperty("predicted_integrated",a.prediction.integrated?juce::var(*a.prediction.integrated):juce::var());v->setProperty("predicted_lra",a.prediction.lra?juce::var(*a.prediction.lra):juce::var());v->setProperty("forecast_limitations",a.prediction.limitations);v->setProperty("affected_frames",juce::String(a.plan.activity.affectedFrames));v->setProperty("connected_regions",int(a.plan.activity.support.size()));v->setProperty("activity",regionActivityDescription(a.plan,48000,media.media.frames));attempts.add(value);
                }
                row->setProperty("complete_measured_attempts",attempts);
            }
            destination.getChildFile("proposal.txt").replaceWithText(result->text());
            if(s.canApplySectionProposal())
            {
                expect(bool(result->assessment()),"applicable held-out correction has a complete actual candidate");
                expect(s.applySectionProposal()&&s.currentStateNeedsRemeasurement(),"owned model Apply leaves Current unverified");auto verified=measuredComposition(s,MeasurementDomain::currentState);
                expect(s.evaluation(MeasurementDomain::currentState).verdict==Verdict::pass,"independent Current passes every declared rule after actual Apply");
                expect(juce::JSON::toString(holdoutMetrics(verified.detail),true)==juce::JSON::toString(holdoutMetrics(result->assessment()->detail),true),"independent Current exactly reproduces all six privately measured scalar metrics");row->setProperty("independent_current",holdoutMetrics(verified.detail));
            }
            else expect(!s.applySectionProposal()&&proposalPlanIdentity(s.currentStateDocument.processingPlan())==applied,"measured failure or unavailable outcome cannot Apply into the owned model");
            expect(sf::save(sf::capture(s),destination.getChildFile("inspection.levelscope"),job).succeeded(),"save owned inspection state with verification authority kept accurate");
            destination.getChildFile("result.json").replaceWithText(juce::JSON::toString(record,true));cases.add(record);std::cout<<input.name<<" "<<method<<": "<<result->outcomeTitle()<<std::endl;
        }
        auto after=validateMedia(file,job);expect(after.succeeded()&&after.media.sha256==input.hash,"original independent audio remains byte-exact after both methods");
    }
}

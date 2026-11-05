// Private qualification only. User audio and diagnostic captures stay outside Git.
namespace
{
void realProgramme(const juce::File& root,const juce::File& file)
{
    JobControl control;const auto media=validateMedia(file,control);expect(media.succeeded(),"bounded real programme validated");
    expect(media.media.sampleRate==48000&&media.media.channels==6&&media.media.frames==48000*480,"fixed eight-minute cinema excerpt, not complete film");
    StandaloneSessionModel s;s.setValidatedSource(file,media.media);s.sourceDocument.layout=layoutFor(6);s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
    const auto source=measuredComposition(s,MeasurementDomain::source);s.initializeCurrentStateFromSource();if(!s.currentState.hasMeasurement)measuredComposition(s,MeasurementDomain::currentState);s.synchronizeSections();
    auto requirement=s.jobRequirements;requirement.maximumLraTenths=150;requirement.provenance="M19 bounded cinema-excerpt qualification example: explicit 15.0 LU ceiling, not a claimed delivery specification";juce::String error;expect(s.setJobRequirements(requirement,error),"fixed real-programme example ceiling selected with provenance");
    const auto request=s.sectionProposalRequest();const auto revision=s.revision;const auto baseline=s.currentState.analysis;
    auto description=proposal::findLraPassages(proposal::analyzeProgramme(baseline->detail,baseline->identity),15.);
    const auto started=juce::Time::getMillisecondCounterHiRes();auto result=generateSectionProposal(request,control);expect(bool(result.value),"bounded real proposal completes without source mutation");auto candidate=*result.value;
    expect(candidate->passes()<=3&&s.currentState.analysis==baseline&&s.revision==revision,"real programme assessment preserves Current and pass budget");
    expect(root.getChildFile("proposal.txt").replaceWithText(candidate->text()),"record real outcome and exact constraints");
    boundaryCaptures(root.getChildFile("listening"),request,*candidate);
    auto record=juce::var(new juce::DynamicObject());auto* r=record.getDynamicObject();r->setProperty("excerpt_sha256",media.media.sha256);r->setProperty("source_integrated_lufs",*source.detail.integrated.value);r->setProperty("source_lra_lu",*source.detail.lra.value);r->setProperty("source_true_peak_dbtp",*source.detail.truePeak.value);r->setProperty("p10",*source.detail.p10);r->setProperty("p95",*source.detail.p95);r->setProperty("target_lra_lu",15.);r->setProperty("private_passes",candidate->passes());r->setProperty("assessment_seconds",(juce::Time::getMillisecondCounterHiRes()-started)/1000);r->setProperty("outcome",candidate->outcomeTitle());r->setProperty("reason",candidate->outcomeReason());r->setProperty("full_film_repeated",false);
    for(bool high:{false,true})
    {
        juce::Array<juce::var> passages;const auto& rows=high?description->high:description->low;
        for(size_t n=0;n<std::min<size_t>(20,rows.size());++n)
        {
            const auto& p=rows[n];auto row=juce::var(new juce::DynamicObject());auto* o=row.getDynamicObject();o->setProperty("begin_frame",juce::int64(p.begin));o->setProperty("end_frame",juce::int64(p.end));o->setProperty("begin_original_frame",juce::int64(77760000+p.begin));o->setProperty("end_original_frame",juce::int64(77760000+p.end));o->setProperty("score_lu_windows",p.score);o->setProperty("included_windows",juce::int64(p.windows));o->setProperty("padding_windows",juce::int64(p.paddedWindows));passages.add(row);
        }
        r->setProperty(high?"high_passages":"quiet_passages",passages);
    }
    if(candidate->status()==SectionProposalStatus::ready)
    {
        expect(s.acceptSectionProposal(candidate)&&s.applySectionProposal(),"explicit disposable real-programme Apply only for exact passing candidate");
        auto verified=measuredComposition(s,MeasurementDomain::currentState);expect(s.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass&&verified.detail.lra.value==candidate->assessment()->detail.lra.value,"independent real Current verifies actual complete plan");
        r->setProperty("independent_current_lra",*verified.detail.lra.value);
    }
    else
    {
        expect(s.acceptSectionProposal(candidate)&&!s.canApplySectionProposal()&&s.currentState.analysis==baseline,"policy/method-limited real result cannot Apply; no tuning to force a pass");
        r->setProperty("current_unchanged",true);
    }
    expect(root.getChildFile("real-programme.json").replaceWithText(juce::JSON::toString(record,true)),"record bounded real evidence without source media");
    expect(root.getChildFile("measurement-only.json").replaceWithText(juce::JSON::toString(ex::report(ex::capture(s)),true)),"honest Source/Current report without fabricated artifact");
    cases.add(record);
}
}

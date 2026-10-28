#pragma once
namespace
{
void cinemaEvaluation(const juce::File& root,const juce::File& file)
{
    const auto started=compositionNow();JobControl job;auto validation=validateMedia(file,job);
    expect(validation.succeeded()&&validation.media.channels==6&&validation.media.sampleRate==48000&&validation.media.bitsPerSample==24,"cinema test payload is supported six-channel 48 kHz PCM24");
    StandaloneSessionModel s;expect(s.setValidatedSource(file,validation.media),"adopt separate cinema evaluation copy");
    s.sourceDocument.layout=layoutFor(6);expect(s.confirmLayout()&&s.initializeCurrentStateFromSource(),"confirm user-stated L R C LFE Ls Rs and whole-programme zero-based range");
    s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);s.synchronizeSections();
    expect(s.sectionDocument.map().sections.size()==1,"user chose one whole-programme section; no inferred scenes");
    auto source=std::async(std::launch::async,[&]{return measuredComposition(s,MeasurementDomain::source);});compositionAwait(source);
    auto current=std::async(std::launch::async,[&]{return measuredComposition(s,MeasurementDomain::currentState);});compositionAwait(current);
    const auto request=s.sectionProposalRequest();auto pending=std::async(std::launch::async,[&]{return generateSectionProposal(request,job);});const auto outcome=compositionAwait(pending);
    expect(outcome.value.has_value(),"single cinema proposal evaluation completes with an honest bounded outcome");const auto result=*outcome.value;
    expect(result->passes()<=3&&!s.currentStateDocument.sectionGains&&s.currentStateDocument.staticGainDb==0,"single cinema attempt preserves applied processing and pass budget");
    expect(root.getChildFile("proposal.txt").replaceWithText(result->text()),"retain exact cinema outcome and policy explanation");
    const auto report=ex::report(ex::capture(s),nullptr);expect(root.getChildFile("measurement-only.json").replaceWithText(juce::JSON::toString(report,true)),"cinema report contains actual Source/Current observations without an artifact claim");
    auto record=juce::var(new juce::DynamicObject());auto* r=record.getDynamicObject();r->setProperty("evaluated_copy_sha256",validation.media.sha256);
    r->setProperty("frames",juce::String(validation.media.frames));r->setProperty("sample_rate",48000);r->setProperty("channels",6);r->setProperty("section_count",1);
    r->setProperty("private_passes",result->passes());r->setProperty("outcome",result->status()==SectionProposalStatus::noSolution?"No solution found within this policy":result->status()==SectionProposalStatus::ready?"Private candidate ready, not applied":result->status()==SectionProposalStatus::globalOnly?"Global-only prediction, not applied":result->status()==SectionProposalStatus::noChange?"No change":"Unavailable");
    r->setProperty("elapsed_ms",compositionNow()-started);r->setProperty("original_modified",false);r->setProperty("audio_played",false);
    r->setProperty("external_reference","User-reported Nuendo, 28 October 2026; exact version/settings and range equivalence unconfirmed. Not exact vectors or an engine-error verdict. Cubase 10 historical screenshot is supplementary.");
    expect(root.getChildFile("evaluation.json").replaceWithText(juce::JSON::toString(record,true)),"retain cinema provenance limits and attempt count");
}
}

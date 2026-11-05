// New method over existing programme readers; no historical suite recursion.
namespace
{
void inputCases(const juce::File& root)
{
    for(int channels:{2,6})
    {
        const auto dir=root.getChildFile(juce::String(channels));auto original=lraFixture(dir,channels);requireLraReduction(original);
        JobControl job;auto reference=generateSectionProposal(original.sectionProposalRequest(),job);
        expect(reference.value&&(*reference.value)->status()==SectionProposalStatus::ready,"interleaved reference candidate available");
        auto verified=prepareVerifiedSource(original.measurementContext(MeasurementDomain::source),job);expect(bool(verified.value),"reader for exact programme variants");
        auto reader=(*verified.value)->createReader();const int frames=int(original.sourceDocument.media.frames);
        juce::AudioBuffer<float> pcm(channels,frames);expect(reader->read(pcm.getArrayOfWritePointers(),channels,0,frames),"complete short fixture decoded without requantization");
        const auto members=component(dir,"members",48000,pcm,true);StandaloneSessionModel grouped;
        expect(grouped.setValidatedProgramme(*members.assembly),"confirmed separate-file stereo/5.1 adopted");
        juce::AudioBuffer<float> half(pcm);half.applyGain(.5f);ProgrammeComposition composition;
        composition.preset=channels==2?CompositionPreset::stereoStereo:CompositionPreset::surroundSurround;
        composition.me=component(dir,"me",48000,half,false);composition.dialogue=component(dir,"dialogue",48000,half,false);
        expect(composition.confirmComplementaryBalancedAlignment(),"independently generated complementary halves confirmed");
        StandaloneSessionModel composed;expect(composed.setValidatedComposition(composition),"explicit complementary composition adopted");
        for(auto* session:{&grouped,&composed})
        {
            session->setSelectedTargetProfileIndex(1);session->setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
            measuredComposition(*session,MeasurementDomain::source);session->initializeCurrentStateFromSource();if(!session->currentState.hasMeasurement)measuredComposition(*session,MeasurementDomain::currentState);session->synchronizeSections();requireLraReduction(*session);
            const auto baseline=session->currentState.analysis,oracle=original.currentState.analysis;
            expect(baseline->detail.integrated.value==oracle->detail.integrated.value&&baseline->detail.lra.value==oracle->detail.lra.value&&baseline->detail.channelTruePeaks==oracle->detail.channelTruePeaks,"separate and composed readers preserve exact baseline measurements");
            auto result=generateSectionProposal(session->sectionProposalRequest(),job);expect(result.value&&(*result.value)->status()==SectionProposalStatus::ready,"new method works on confirmed group and composition");
            const auto p=*result.value,r=*reference.value;const auto& a=p->plan().sectionGains->intervals;const auto& b=r->plan().sectionGains->intervals;bool same=a.size()==b.size();
            for(size_t n=0;same&&n<a.size();++n)same=a[n].begin==b[n].begin&&a[n].end==b[n].end&&a[n].db==b[n].db;
            expect(same&&p->plan().staticGainDb==r->plan().staticGainDb&&p->assessment()->detail.integrated.value==r->assessment()->detail.integrated.value&&p->assessment()->detail.lra.value==r->assessment()->detail.lra.value&&p->assessment()->detail.truePeakEighthFrames==r->assessment()->detail.truePeakEighthFrames,"equivalent complete samples yield identical gains, regions, measurements and peak positions");
            expect(p->identity()!=r->identity()&&!session->acceptSectionProposal(r),"distinct input definitions never share proposal authority despite sample equality");
            expect(session->acceptSectionProposal(p)&&session->applySectionProposal(),"own exact candidate explicitly applies");
            measuredComposition(*session,MeasurementDomain::currentState);expect(session->evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"independent grouped/composed Current passes selected requirements");
            const auto path=dir.getChildFile(session==&grouped?"group.levelscope":"composition.levelscope");expect(sf::save(sf::capture(*session),path,job).succeeded(),"persist input definition and historical LRA processing");
            auto saved=sf::read(path,job);expect(bool(saved.document),"read grouped/composed session");auto reopened=sf::verify(*saved.document,job);
            expect(reopened.prepared&&reopened.prepared->model.lastAppliedExplanation&&!reopened.prepared->model.sectionProposal,"reopen preserves history and exact components without candidate authority");
        }
    }
}
}

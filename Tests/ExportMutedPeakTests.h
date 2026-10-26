#pragma once
namespace
{
void exportMutedPeaks()
{
    Temp temp;constexpr int frames=192000;
    for(int owner:{3,5})
    {
        juce::AudioBuffer<float> data(6,frames);for(int ch=0;ch<6;++ch)for(int n=0;n<frames;++n)data.setSample(ch,n,.008f*float(std::sin(n*(.071+ch*.013))));
        const int peak=frames/2+19;data.setSample(owner,peak,.7f);
        for(int i=0;i<4;++i)
        {
            measurement::OfflineMeasurement meter(48000,OfflineAnalysisService::measurementRoles(layoutFor(6)));meter.push(data.getArrayOfReadPointers(),frames);const auto result=meter.finish();
            const double desired=*result.integrated.value+22.1;data.setSample(owner,peak,float(data.getSample(owner,peak)*std::pow(10.,(desired-*result.truePeak.value)/20.)));
        }
        const auto file=temp.directory.getChildFile("peak-"+juce::String(owner)+".wav");auto stream=file.createOutputStream();juce::WavAudioFormat wav;
        std::unique_ptr<juce::AudioFormatWriter> writer(wav.createWriterFor(stream.get(),48000,6,24,{},0));expect(writer!=nullptr,"highest peak fixture writer");stream.release();expect(writer->writeFromAudioSampleBuffer(data,0,frames),"highest peak fixture payload");writer.reset();
        JobControl control;auto validation=validateMedia(file,control);expect(validation.succeeded(),"highest peak fixture verified");StandaloneSessionModel s;
        expect(s.setValidatedSource(file,validation.media)&&s.confirmLayout()&&s.initializeCurrentStateFromSource(),"conventional 5.1 peak fixture confirmed");s.setSelectedTargetProfileIndex(1);s.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
        measure(s,MeasurementDomain::source);measure(s,MeasurementDomain::currentState);const auto initial=s.currentState.analysis;const auto& peaks=initial->detail.channelTruePeaks;
        expect(std::max_element(peaks.begin(),peaks.end())-peaks.begin()==owner,"unique highest true peak belongs to LFE or surround");
        const auto draft=proposal::draftGlobalTrim(s.currentTrimRequest());const auto evaluation=s.evaluation(MeasurementDomain::currentState).identity;
        expect(draft.status==proposal::TrimStatus::draft&&draft.predictedTruePeak&&std::abs(*draft.predictedTruePeak+1)<1e-12,"all-channel maximum constrains proposed gain");
        for(int ch=0;ch<6;++ch)s.setOutputRoute(ch,ch);s.setOutputRoute(owner,-1);
        expect(s.currentState.analysis==initial&&s.evaluation(MeasurementDomain::currentState).identity==evaluation&&proposal::draftGlobalTrim(s.currentTrimRequest()).identity==draft.identity,"muting highest LFE/surround changes no measurements, evaluation or peak-constrained proposal");
        // Deliberate test trim provides a small strict-ceiling margin for dither; no automatic product correction.
        s.setCurrentGainDb(*draft.totalGainDb-.001);Capture offline(6,frames);measure(s,MeasurementDomain::currentState,&offline);const auto request=ex::capture(s);Capture captured(6,frames);
        auto a=ex::render(request,temp.directory,control,{},&captured);expect(a.value&&identical(offline,captured)&&captured.data.getMagnitude(owner,0,frames)>0,"muted peak owner remains present exactly in render/offline samples");
        s.setOutputRoute(owner,owner);auto b=ex::render(ex::capture(s),temp.directory,control,{},nullptr,(*a.value)->verifiedSource());
        expect(b.value&&(*b.value)->hash()==(*a.value)->hash()&&ex::audioIdentity(request)==ex::audioIdentity(ex::capture(s)),"monitor mute/unmute cannot alter PCM bytes or export identity");
    }
}
}

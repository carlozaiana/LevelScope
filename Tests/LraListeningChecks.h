// Diagnostic float captures only: not a published PCM24 export or delivery approval.
namespace
{
void boundaryCaptures(const juce::File& root,const SectionProposalRequest& request,const SectionProposalResult& result)
{
    if(!result.assessment()||!result.lraSelection()||result.lraSelection()->selected.empty())return;
    expect(root.createDirectory().wasOk(),"owned boundary-capture directory");
    JobControl job;auto verified=prepareVerifiedSource(request.currentContext,job);expect(bool(verified.value),"verified reader for bounded diagnostic captures");auto reader=(*verified.value)->createReader();
    const auto& region=result.lraSelection()->selected.front().passage;const int rate=int(request.currentContext.media.sampleRate),channels=request.currentContext.media.channels;
    juce::Array<juce::var> records;
    for(bool ending:{false,true})
    {
        const auto boundary=ending?region.end:region.begin;
        const int64_t first=std::clamp<int64_t>(boundary-int64_t(rate)*4,0,request.currentContext.media.frames-int64_t(rate)*8);
        juce::AudioBuffer<float> raw(channels,rate*8);expect(reader->read(raw.getArrayOfWritePointers(),channels,first,raw.getNumSamples()),"bounded eight-second onset/return read");
        for(bool after:{false,true})
        {
            const auto& plan=after?result.assessment()->context.processing:request.currentContext.processing;juce::AudioBuffer<float> data(raw);
            ProgrammeGainProcessor processor(plan);expect(processor.process(data.getArrayOfWritePointers(),channels,data.getNumSamples(),first),"capture uses exact complete measured plan, without monitor gain or dither");
            const auto name=juce::String(ending?"return":"onset")+(after?"-candidate":"-baseline");const auto file=root.getChildFile(name+".wav");floatFile(file,rate,data);
            auto row=juce::var(new juce::DynamicObject());auto* r=row.getDynamicObject();r->setProperty("file",file.getFileName());r->setProperty("sha256",juce::SHA256(file).toHexString());r->setProperty("begin_source_frame",juce::int64(first));r->setProperty("end_source_frame",juce::int64(first+data.getNumSamples()));r->setProperty("boundary_source_frame",juce::int64(boundary));r->setProperty("processing_plan",proposalPlanIdentity(plan));r->setProperty("channels",channels);r->setProperty("maximum_absolute_sample",data.getMagnitude(0,data.getNumSamples()));
            if(channels==6)
            {
                juce::AudioBuffer<float> pair(2,data.getNumSamples());pair.copyFrom(0,0,data,0,0,data.getNumSamples());pair.copyFrom(1,0,data,1,0,data.getNumSamples());
                const auto partial=root.getChildFile(name+"-LEFT-RIGHT-ONLY.wav");floatFile(partial,rate,pair);r->setProperty("optional_partial_capture",partial.getFileName());r->setProperty("partial_sha256",juce::SHA256(partial).toHexString());r->setProperty("partial_note","Original Left/Right channels only; no downmix. Does not qualify Centre/LFE/surround listening or complete 5.1 transparency.");
            }
            records.add(row);
        }
    }
    auto receipt=juce::var(new juce::DynamicObject());auto* r=receipt.getDynamicObject();r->setProperty("purpose","Diagnostic boundary listening preparation only; no audio device opened");r->setProperty("source_sha256",request.currentContext.media.sha256);r->setProperty("candidate_meets_requirements",result.status()==SectionProposalStatus::ready);r->setProperty("candidate_verification",result.assessment()->identity);r->setProperty("sample_rate",rate);r->setProperty("dither",false);r->setProperty("published_export",false);r->setProperty("captures",records);
    expect(root.getChildFile("manifest.json").replaceWithText(juce::JSON::toString(receipt,true)),"diagnostic listening identity and limitations recorded");
}
}

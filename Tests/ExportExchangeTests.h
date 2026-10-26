#pragma once
namespace
{
void exchangeExports(const juce::File& root)
{
    auto bundles=root.findChildFiles(juce::File::findDirectories,false,"verified-*");bundles.sort();
    expect(bundles.size()==16,"complete opposite-platform bundle matrix");
    for(const auto& folder:bundles)
    {
        const auto saved=juce::JSON::parse(folder.getChildFile("report.json"));const auto old=saved["artifact"];
        const auto audio=folder.getChildFile("programme.wav");JobControl control;
        expect(pcm::hashFile(audio,control)==old["sha256"].toString(),"exchanged exact PCM bytes match original artifact hash");
        auto media=validateMedia(audio,control);expect(media.succeeded(),"exchanged PCM independently imported");
        StandaloneSessionModel session;expect(session.setValidatedSource(audio,media.media)&&session.confirmLayout(),"exchanged explicit metadata roles confirmed");
        session.setSelectedTargetProfileIndex(1);session.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
        expect(!session.source.hasMeasurement,"saved report verdict never restores a measurement");
        measure(session,MeasurementDomain::source);
        const auto fresh=ex::report(ex::capture(session))["source"];
        for(const char* name:{"integrated","true_peak","sample_peak","lra","momentary_max","short_term_max"})
        {
            const auto a=old["metrics"][name],b=fresh["metrics"][name];
            expect(a["status"].toString()==b["status"].toString(),"exchanged metric status preserved");
            if(a["status"].toString()=="available")expect(std::abs(double(a["value"])-double(b["value"]))<=(juce::String(name)=="lra"?.1:.01),"exchanged measured metric within documented cross-platform allowance");
        }
        expect(fresh["frames"].toString()==old["frames"].toString()&&fresh["sample_rate"]==old["sample_rate"]
            &&fresh["measured_content_sha256"].toString()==old["sha256"].toString(),"exchanged frame/rate/content identity exact");
        const auto* a=old["per_channel_peaks"].getArray();const auto* b=fresh["per_channel_peaks"].getArray();
        expect(a&&b&&a->size()==b->size(),"exchanged all-channel peak records");
        for(int i=0;i<a->size();++i)for(const char* key:{"channel","role","sample_peak_frame","true_peak_eighth_frame"})
            expect((*a)[i][key].toString()==(*b)[i][key].toString(),"exchanged ordered channel identity and peak positions exact");
        expect(session.evaluation(MeasurementDomain::source).verdict==proposal::Verdict::pass,"exchanged PCM passes freshly measured selected rules");
        auto row=fresh;row.getDynamicObject()->setProperty("bundle",folder.getFileName());cases.add(row);
    }
}
}

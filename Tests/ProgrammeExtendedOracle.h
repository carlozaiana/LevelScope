#pragma once

namespace
{
juce::Array<juce::var> extendedOracleEvidence;

Fixture makeExtendedOracle(const juce::File& root,int channels,int rate,bool silence)
{
    expect(root.createDirectory().wasOk(),"extended oracle directory");
    Fixture f;f.rate=rate;f.frames=rate*(silence?4:64)+17;
    f.interleaved=root.getChildFile("independent-interleaved.wav");
    juce::WavAudioFormat wav;juce::AiffAudioFormat aiff;juce::FlacAudioFormat flac;
    auto stream=f.interleaved.createOutputStream();
    std::unique_ptr<juce::AudioFormatWriter> interleaved(wav.createWriterFor(stream.get(),rate,channels,24,{},0));
    expect(interleaved!=nullptr,"extended independent interleaved writer");stream.release();
    std::vector<juce::File> files;std::vector<std::unique_ptr<juce::AudioFormatWriter>> writers;
    for(int ch=0;ch<channels;++ch)
    {
        auto file=root.getChildFile("mono-"+juce::String(ch)+(ch%3==0?".wav":ch%3==1?".aiff":".flac"));
        files.push_back(file);auto out=file.createOutputStream();
        juce::AudioFormat* format=ch%3==0?static_cast<juce::AudioFormat*>(&wav):ch%3==1?static_cast<juce::AudioFormat*>(&aiff):static_cast<juce::AudioFormat*>(&flac);
        std::unique_ptr<juce::AudioFormatWriter> writer(format->createWriterFor(out.get(),rate,1,ch%2?16:24,{},0));
        expect(writer!=nullptr,"extended mono writer");out.release();writers.push_back(std::move(writer));
    }
    // Bounded fixture generation; both representations are written independently
    // from this known PCM16-grid signal, never decoded from one another.
    juce::AudioBuffer<float> block(channels,4096),mono(1,4096);
    for(int begin=0;begin<f.frames;begin+=4096)
    {
        const int count=std::min(4096,f.frames-begin);
        for(int ch=0;ch<channels;++ch)for(int n=0;n<count;++n)
        {
            const int frame=begin+n;
            const int amplitude=silence||frame<8*rate||frame>=56*rate?0:frame<24*rate?150:frame<40*rate?5000:1500;
            const int phase=(frame*(ch+1)+ch*7)%128;
            const int triangle=phase<64?phase-32:96-phase;
            block.setSample(ch,n,static_cast<float>(amplitude*(ch+1)*triangle/(channels*32))/32768.f);
        }
        if(!interleaved->writeFromAudioSampleBuffer(block,0,count))throw std::runtime_error("extended independent oracle block write");
        for(int ch=0;ch<channels;++ch)
        {mono.copyFrom(0,0,block,ch,0,count);if(!writers[ch]->writeFromAudioSampleBuffer(mono,0,count))throw std::runtime_error("extended mono block write");}
    }
    interleaved.reset();writers.clear();JobControl c;auto inspected=inspectProgrammeFiles(files,c);
    expect(inspected.error==ErrorCode::none,"extended group validates");f.group=std::move(inspected.group);
    const auto roles=assemblyRoles(channels);for(int ch=0;ch<channels;++ch)f.group.members[ch].role=roles[ch];
    expect(f.group.confirmAlignment(),"known generated extended alignment confirmed");return f;
}

juce::var extendedMetric(const measurement::Metric& metric)
{
    auto v=juce::var(new juce::DynamicObject());v.getDynamicObject()->setProperty("status",static_cast<int>(metric.status));
    v.getDynamicObject()->setProperty("value",metric.value?juce::var(*metric.value):juce::var());return v;
}

void extendedOracles(const juce::File& root)
{
    for(int channels:{2,6})for(int rate:{44100,48000,96000,192000})for(bool silence:{false,true})
    {
        auto dir=root.getChildFile(juce::String(channels)+"-"+juce::String(rate)+(silence?"-silence":"-dynamic"));
        auto f=makeExtendedOracle(dir,channels,rate,silence);auto s=model(f.group);s.setCurrentGainDb(-3.305074118944);
        for(auto domain:{MeasurementDomain::source,MeasurementDomain::currentState})
        {
            measure(s,domain);const auto& actual=(domain==MeasurementDomain::source?s.source:s.currentState).analysis->detail;
            JobControl c;auto oracle=OfflineAnalysisService().analyze(interleavedContext(f,s.measurementContext(domain)),domain,c);
            expect(oracle.value&&equalMetrics(actual,oracle.value->detail),"extended full-resolution M/S/LRA gates percentiles and all-channel peaks equal interleaved oracle");
            if(silence)
                expect(!actual.integrated.value&&!actual.lra.value&&actual.integratedIncluded==0&&actual.lraIncluded==0,"silent assembly retains unavailable loudness and no included gate windows");
            else
                expect(actual.lraStable&&actual.lra.value&&*actual.lra.value>1&&actual.integratedIncluded>0&&actual.lraIncluded>0
                    &&actual.integratedIncluded<actual.momentary.size()&&actual.lraIncluded<actual.lraWindows.size(),"dynamic oracle exercises stable nonzero LRA and excluded windows");
            auto row=juce::var(new juce::DynamicObject());auto* o=row.getDynamicObject();
            o->setProperty("channels",channels);o->setProperty("rate",rate);o->setProperty("frames",f.frames);o->setProperty("silence",silence);
            o->setProperty("domain",domain==MeasurementDomain::source?"Source":"Current State");o->setProperty("gain_db",domain==MeasurementDomain::source?0:s.currentStateDocument.staticGainDb);
            o->setProperty("lra_stable",actual.lraStable);o->setProperty("integrated_included",static_cast<juce::int64>(actual.integratedIncluded));o->setProperty("lra_included",static_cast<juce::int64>(actual.lraIncluded));
            auto metrics=juce::var(new juce::DynamicObject());
            metrics.getDynamicObject()->setProperty("integrated",extendedMetric(actual.integrated));metrics.getDynamicObject()->setProperty("lra",extendedMetric(actual.lra));
            metrics.getDynamicObject()->setProperty("momentary",extendedMetric(actual.momentaryMax));metrics.getDynamicObject()->setProperty("short_term",extendedMetric(actual.shortTermMax));
            metrics.getDynamicObject()->setProperty("sample_peak",extendedMetric(actual.samplePeak));metrics.getDynamicObject()->setProperty("true_peak",extendedMetric(actual.truePeak));o->setProperty("metrics",metrics);
            juce::Array<juce::var> windows,sampleFrames,trueFrames;
            for(int family=0;family<3;++family)
            {
                const auto& list=family==0?actual.momentary:family==1?actual.shortTerm:actual.lraWindows;
                for(const auto& w:list){juce::Array<juce::var> item;item.add(family);item.add(static_cast<juce::int64>(w.endFrame));item.add(w.energy);windows.add(item);}
            }
            for(auto frame:actual.samplePeakFrames)sampleFrames.add(static_cast<juce::int64>(frame));
            for(auto frame:actual.truePeakEighthFrames)trueFrames.add(static_cast<juce::int64>(frame));
            o->setProperty("windows",windows);o->setProperty("sample_peak_frames",sampleFrames);o->setProperty("true_peak_eighth_frames",trueFrames);
            juce::Array<juce::var> samplePeaks,truePeaks;for(auto value:actual.channelSamplePeaks)samplePeaks.add(value);for(auto value:actual.channelTruePeaks)truePeaks.add(value);
            o->setProperty("sample_peaks",samplePeaks);o->setProperty("true_peaks",truePeaks);
            o->setProperty("integrated_gate",actual.integratedRelativeGate?juce::var(*actual.integratedRelativeGate):juce::var());
            o->setProperty("lra_gate",actual.lraRelativeGate?juce::var(*actual.lraRelativeGate):juce::var());
            o->setProperty("p10",actual.p10?juce::var(*actual.p10):juce::var());o->setProperty("p95",actual.p95?juce::var(*actual.p95):juce::var());
            o->setProperty("same_platform_full_result_exact",true);extendedOracleEvidence.add(row);
            std::cout<<"Extended oracle "<<channels<<" channels "<<rate<<" Hz "<<(silence?"silence":"dynamic")<<" "<<(domain==MeasurementDomain::source?"Source":"Current")<<" passed"<<std::endl;
        }
        expect(dir.deleteRecursively(),"extended oracle owned fixtures cleaned");
    }
    expect(extendedOracleEvidence.size()==32,"complete extended oracle rate/layout/domain/silence matrix");
}
}

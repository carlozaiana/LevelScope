#pragma once
// Included inside RegionPlanningTests' private namespace. No real device access.
double independentEnvelope(const Plan& p,int64_t frame)
{
    double db=0;for(const auto& row:p.intervals)if(frame>=row.begin&&frame<row.end){db=row.db;break;}
    for(size_t n=1;n<p.intervals.size();++n)
    {
        const auto& left=p.intervals[n-1];const auto& right=p.intervals[n];
        const auto radius=std::min({int64_t(std::llround(p.sampleRate*.05)),(left.end-left.begin)/2,(right.end-right.begin)/2});
        const auto begin=right.begin-radius,end=right.begin+radius;
        if(frame==begin)db=left.db;else if(frame==end-1)db=right.db;
        else if(frame>begin&&frame<end-1)db=left.db+(right.db-left.db)*(double(frame-begin)/double(end-begin-1));
    }
    return db;
}
void audioCase(const juce::File& root,int rate,int channels,int mode)
{
    expect(root.createDirectory().wasOk(),"owned short audio parity root");const int frames=8*rate;
    juce::AudioBuffer<float> raw(channels,frames);StandaloneSessionModel s;
    if(mode==2)
    {
        auto fixture=fixtureComposition(root,rate,CompositionPreset::surroundMonoCentre,true,frames);raw=fixture.expected;s=compositionModel(fixture.recipe);
    }
    else
    {
        for(int ch=0;ch<channels;++ch)for(int n=0;n<frames;++n)raw.setSample(ch,n,float(.01*(ch+1)*std::sin(2*juce::MathConstants<double>::pi*(701+ch*97)*n/rate)));
        auto input=component(root,"programme",rate,raw,mode==1);
        if(input.assembly)expect(s.setValidatedProgramme(*input.assembly),"exact separate-file stereo group adopted");
        else {expect(s.setValidatedSource(input.file,input.media),"exact interleaved parity source adopted");s.sourceDocument.layout=input.layout;}
        expect(s.initializeCurrentStateFromSource(),"initialize parity Current");
    }
    auto context=s.measurementContext(MeasurementDomain::currentState);Plan before;before.sampleRate=rate;before.frames=frames;before.inputIdentity=sectionInputIdentity(context);before.intervals={{0,frames,0}};
    auto c=selected(before,{2*rate,6*rate});c.settings.minimumIntervalMs=500;ProposalPolicy policy;policy.durationThousandths=1000;
    auto plan=constructRegionPlan(before,0,c,policy,{assignment(c.roots[0].hint,{3*rate,5*rate},{1.123456789,0,2.234567891})},0);
    expect(plan.error.isEmpty(),"three-part full-precision envelope is valid at each rate");
    s.currentStateDocument.sectionGains=plan.envelope;expect(s.setCurrentGainDb(-.123456789),"explicit global stage invalidates Current");
    Capture offline(channels,frames);auto measured=measuredComposition(s,MeasurementDomain::currentState,&offline);expect(offline.valid&&offline.next==frames,"complete all-channel offline capture");
    juce::AudioBuffer<float> expected(channels,frames);bool exterior=true;
    for(int n=0;n<frames;++n)
    {
        const double db=independentEnvelope(plan.envelope,n);if(n<2*rate||n>=6*rate)exterior&=db==0;
        for(int ch=0;ch<channels;++ch){const float local=db==0?raw.getSample(ch,n):float(double(raw.getSample(ch,n))*std::pow(10.,db/20.));expected.setSample(ch,n,float(double(local)*std::pow(10.,-.123456789/20.)));}
    }
    expect(exterior&&sameSamples(expected,offline.data),"independent ramp oracle matches every processed sample including LFE and unchanged exterior");
    auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));OutputRouting muted;s.selectAuditionDomain(MeasurementDomain::currentState);auto request=s.auditionRequest();audition.synchronize(request,muted);expect(audition.connect(rate).isEmpty(),"injected device only for new geometry");
    for(int64_t first:{int64_t(0),int64_t(2*rate-2048),int64_t(3*rate-2048),int64_t(5*rate-2048),int64_t(frames-4096)})
    {
        expect(audition.seek(first,request,muted).isEmpty(),"audio-only exact-frame seek");Capture capture(channels,4096,first);audition.setCaptureObserver(&capture);expect(audition.play(request,muted,false,true).isEmpty(),"muted injected playback");awaitPlaying(audition,*injected);juce::AudioBuffer<float> output(maxChannels,256);
        for(int n=0;n<16;++n){injected->tick(output);audition.poll();std::this_thread::sleep_for(std::chrono::milliseconds(1));}audition.stop();audition.setCaptureObserver(nullptr);
        bool same=capture.valid&&capture.next==first+4096;for(int ch=0;ch<channels&&same;++ch)same=std::memcmp(capture.data.getReadPointer(ch),expected.getReadPointer(ch,int(first)),4096*sizeof(float))==0;
        expect(same&&audition.underruns()==0,"all-channel pre-monitor capture equals offline and independent oracle across seek/ramp/EOF");
    }
    audition.disconnect();expect(measured.detail.frames==frames,"no processing duration change");
}
void audioCases(const juce::File& root)
{
    for(int rate:{44100,48000,96000,192000})for(int channels:{2,6})audioCase(root.getChildFile(juce::String(rate)+"-"+juce::String(channels)),rate,channels,0);
    audioCase(root.getChildFile("group-stereo"),48000,2,1);audioCase(root.getChildFile("composition"),48000,6,2);
}

// Display-only regressions: no device, analysis job, proposal or processing edit.
void appearanceChecks(const juce::File& root)
{
    MeasurementTimeline ui;ui.setSize(1150,540);ui.setProgramme(240000,48000,"appearance");ui.setWaveform({},false,{"L","R"});
    auto a=std::make_shared<proposal::ProgramAnalysis>();a->sampleRate=48000;a->frames=240000;a->integratedGate=-40;a->lraGate=-50;
    auto original=std::make_shared<measurement::Result>();original->sampleRate=48000;original->frames=a->frames;
    for(int n=0;n<4;++n)
    {proposal::DiagnosticWindow w;w.end=19200+n*4800;w.begin=w.end-19200;w.energy=n==2?.02:.01;w.lufs=measurement::OfflineMeasurement::loudness(w.energy);w.included=true;a->momentary.push_back(w);original->momentary.push_back({w.end,.01});}
    for(int n=0;n<3;++n)
    {proposal::DiagnosticWindow w;w.begin=n*4800;w.end=w.begin+144000;w.energy=.01;w.lufs=n==1?-60:-30;w.included=n!=1;a->lraWindows.push_back(w);w.included=false;a->shortTerm.push_back(w);}
    ui.setAnalysis(a);ui.setReferenceMeasurement(original,true);
    expect(std::any_of(ui.changedMomentaryColumns().begin(),ui.changedMomentaryColumns().end(),[](bool b){return b;}),"Current-only comparison marks differing measured M samples");
    size_t included=0,excluded=0;for(const auto& bin:ui.shortTermGatePlot()){included+=bin.included.count;excluded+=bin.excluded.count;}
    expect(included==2&&excluded==1&&a->shortTerm[1].lufs==-60&&!a->shortTerm[1].included,"short-term colour retains exact measured LRA membership without mutating ungated measurements");
    ui.setReferenceMeasurement(original,false);expect(ui.changedMomentaryColumns().empty()&&a->integratedGate==-40,"Original suppresses red comparison only; measured relative gate is unchanged");
    auto preview=std::make_shared<LoudnessPreview>(240000,48000);preview->window(false,{19200,.02});preview->publish(24000,true);ui.setPreview(preview);
    expect(ui.shortTermGatePlot().empty()&&ui.changedMomentaryColumns().empty(),"incomplete previews never acquire completed gate membership or difference highlights");ui.setPreview({});
    auto& blend=ui.visibilityBlend();expect(blend.waveformAlpha()==1&&blend.loudnessAlpha()==1,"centre preserves original full visibility of both layers");
    blend.setValue(1,juce::sendNotificationSync);expect(blend.waveformAlpha()==0&&blend.loudnessAlpha()==1,"up hides only waveform");imageOf(ui,root.getChildFile("loudness-only.png"));
    blend.setValue(-1,juce::sendNotificationSync);expect(blend.waveformAlpha()==1&&blend.loudnessAlpha()==0,"down hides only loudness");
    blend.setValue(.5,juce::sendNotificationSync);expect(blend.waveformAlpha()==.5&&blend.loudnessAlpha()==1,"upper half fades waveform continuously without dimming loudness");blend.setValue(0,juce::dontSendNotification);
    ui.setViewRange(24000,200000);ui.highlight(60000,80000,false);const auto before=ui.viewRange();const auto box=ui.graphBounds();const float px=box.getX()+.75f*(box.getWidth()-1);
    const auto t=juce::Time::getCurrentTime();const juce::Point<float> p{px,float(box.getCentreY())};
    juce::MouseEvent e{juce::Desktop::getInstance().getMainMouseSource(),p,juce::ModifierKeys(juce::ModifierKeys::commandModifier),0,0,0,0,0,&ui,&ui,t,p,t,1,false};
    juce::MouseWheelDetails wheel{};wheel.deltaY=1;ui.mouseWheelMove(e,wheel);const auto after=ui.viewRange();
    expect(std::abs((before.first+.75*(before.second-before.first))-(after.first+.75*(after.second-after.first)))<.1,"Command-wheel anchors the frame at the pointer, independent of selected section");
    ui.setProgramme(512,48000,"peaks");auto overview=waveform::Overview::create({512,48000,2},std::make_shared<waveform::Budget>());
    juce::AudioBuffer<float> samples(2,512);samples.clear();samples.setSample(0,4,-.8f);samples.setSample(0,75,.5f);samples.setSample(1,75,.25f);
    expect(overview->append(samples.getArrayOfReadPointers(),2,512,0)&&overview->finish(),"known PCM extrema fixture");ui.setWaveform(overview,true,{"L","R"});
    expect(ui.rangePeakText(0,128).contains("selected sample peak -1.938 dBFS at frame 4"),"sample peak readout retains absolute magnitude and actual minimum sample position");
    const auto boundary=ui.rangePeakText(50,80);expect(boundary.contains("covering overview peak")&&boundary.contains("boundary bins/coverage")&&boundary.contains("Regional true peak is not available"),"boundary overlap is explicit and PCM is never labelled true peak");
    ui.setWaveform(overview,false,{"L","R"});expect(ui.rangePeakText(0,128).contains("INCOMPLETE"),"unaccepted waveform cannot assert completed selected peak");ui.setWaveform(overview,true,{"L","R"});
    blend.setValue(-1,juce::sendNotificationSync);imageOf(ui,root.getChildFile("waveform-only.png"));blend.setValue(0,juce::dontSendNotification);imageOf(ui,root.getChildFile("blend-centre.png"));
    auto surround=waveform::Overview::create({512,48000,6},std::make_shared<waveform::Budget>());juce::AudioBuffer<float> six(6,512);six.clear();
    for(int ch=0;ch<6;++ch)six.setSample(ch,80+ch,float(.1*(ch+1)));
    expect(surround->append(six.getArrayOfReadPointers(),6,512,0)&&surround->finish(),"six-channel sample overview fixture");ui.setWaveform(surround,true,{"L","R","C","LFE","Ls","Rs"});
    const auto sixPeaks=ui.rangePeakText(64,128);expect(sixPeaks.contains("LFE: selected sample peak")&&sixPeaks.contains("Rs: selected sample peak")&&sixPeaks.contains("frame 85"),"5.1 selection retains each role and its own sample peak position");
    imageOf(ui,root.getChildFile("surround-blend.png"));
}

// Display-only regressions: no device, analysis job, proposal or processing edit.
void appearanceChecks(const juce::File& root)
{
    MeasurementTimeline ui;ui.setSize(1150,540);ui.setProgramme(240000,48000,"appearance");ui.setWaveform({},false,{"L","R"});
    auto a=std::make_shared<proposal::ProgramAnalysis>();a->sampleRate=48000;a->frames=240000;a->integratedGate=-40;a->lraGate=-50;
    auto original=std::make_shared<measurement::Result>();original->sampleRate=48000;original->frames=a->frames;
    const double sourceEnergy[]{1.,.01,.1,.03,.002,.001},currentEnergy[]{2.,.1,.01,.03,.003,.001};
    measurement::Result current;
    for(int n=0;n<6;++n){original->momentary.push_back({19200+n*4800,sourceEnergy[n]});current.momentary.push_back({19200+n*4800,currentEnergy[n]});}
    const auto sourceGate=measurement::relativeGateEnergy(original->momentary,false),currentGate=measurement::relativeGateEnergy(current.momentary,false);
    original->integratedRelativeGate=measurement::OfflineMeasurement::loudness(*sourceGate);a->integratedGate=measurement::OfflineMeasurement::loudness(*currentGate);
    for(const auto& m:current.momentary)
    {proposal::DiagnosticWindow w;w.end=m.endFrame;w.begin=w.end-19200;w.energy=m.energy;w.lufs=measurement::OfflineMeasurement::loudness(w.energy);w.included=measurement::includedInPopulation(w.energy,*currentGate,false);a->momentary.push_back(w);}
    for(int n=0;n<3;++n)
    {proposal::DiagnosticWindow w;w.begin=n*4800;w.end=w.begin+144000;w.energy=.01;w.lufs=n==1?-60:-30;w.included=n!=1;a->lraWindows.push_back(w);w.included=false;a->shortTerm.push_back(w);}
    ui.setAnalysis(a);ui.setReferenceMeasurement(original,true);
    const auto redAt=[&](double seconds){const auto range=ui.viewRange();const auto& bins=ui.changedMomentaryColumns();const auto n=size_t((seconds*48000-range.first)*bins.size()/(range.second-range.first));return n<bins.size()&&bool(bins[n]);};
    expect(redAt(.45)&&redAt(.55),"red marks both excluded-to-included and included-to-excluded measured IL membership changes");
    expect(redAt(.65),"unchanged momentary energy is red when its measured gate membership changes because the programme gate moves");
    expect(!redAt(.35)&&!redAt(.75)&&!redAt(.85),"changed energy with unchanged membership and unchanged energy are not red");
    ui.setViewRange(.42*48000,.48*48000);expect(redAt(.45),"red cadence cell survives a viewport wholly inside the cell before its window end");ui.setViewRange(0,240000);
    auto uniform=std::make_shared<proposal::ProgramAnalysis>(*a);uniform->integratedGate=*original->integratedRelativeGate+measurement::OfflineMeasurement::loudness(4.)-measurement::OfflineMeasurement::loudness(1.);
    for(size_t n=0;n<uniform->momentary.size();++n){auto& w=uniform->momentary[n];w.energy=sourceEnergy[n]*4;w.included=measurement::includedInPopulation(sourceEnergy[n],*sourceGate,false);}
    ui.setAnalysis(uniform);expect(std::none_of(ui.changedMomentaryColumns().begin(),ui.changedMomentaryColumns().end(),[](bool b){return b;}),"uniform gain with unchanged population never turns the relative line red");ui.setAnalysis(a);
    ui.setViewRange(0,48000);imageOf(ui,root.getChildFile("relative-membership.png"));ui.setViewRange(0,240000);
    ui.setReferenceMeasurement({},true);expect(ui.changedMomentaryColumns().empty(),"missing Original measurement cannot create membership highlights");ui.setReferenceMeasurement(original,true);
    size_t included=0,excluded=0;for(const auto& bin:ui.shortTermGatePlot()){included+=bin.included.count;excluded+=bin.excluded.count;}
    expect(included==2&&excluded==1&&a->shortTerm[1].lufs==-60&&!a->shortTerm[1].included,"short-term colour retains exact measured LRA membership without mutating ungated measurements");
    ui.setReferenceMeasurement(original,false);expect(ui.changedMomentaryColumns().empty()&&a->integratedGate==measurement::OfflineMeasurement::loudness(*currentGate),"Original suppresses red comparison only; measured relative gate is unchanged");
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
    ui.setPeakTarget(-6,true);expect(!ui.peakVisibility().getToggleState()&&ui.peakVisibility().isEnabled(),"over-target overlay defaults off with a measured finite target");
    ui.peakVisibility().setToggleState(true,juce::sendNotificationSync);
    expect(ui.sampleExceedsTarget(-.8,4)&&!ui.sampleExceedsTarget(.5,75)&&!ui.sampleExceedsTarget(std::pow(10.,-6./20),75),"peak overlay uses magnitude, exact retained position and strict numeric ceiling without rounded comparisons");
    const auto hoverBox=ui.graphBounds();
    const juce::Point<float> peakPoint{float(hoverBox.getX()+4.*(hoverBox.getWidth()-1)/512.),float(hoverBox.getY()+hoverBox.getHeight()*.25+.8*.44*hoverBox.getHeight()/2.)};
    const auto hover=ui.peakHoverText(peakPoint);
    expect(hover.contains("Original | L")&&hover.contains("Sample peak: -1.938 dBFS")&&hover.contains("-4.062 dB")&&hover.contains("4 frames")&&!hover.contains("Above numeric ceiling"),"hover identifies actual negative sample, role, exact frame and additional reduction without redundant delta");
    ui.setReferenceMeasurement({},true);expect(ui.peakHoverText(peakPoint).startsWith("Current"),"hover explicitly identifies Current without adding monitor/global gain again");
    blend.setValue(1,juce::sendNotificationSync);expect(ui.peakHoverText(peakPoint).isEmpty(),"hidden waveform markers have no hover target");blend.setValue(0,juce::dontSendNotification);
    ui.setPeakTarget(-6,false);expect(ui.peakHoverText(peakPoint).isEmpty(),"unverified measurement cannot supply hover guidance");ui.setPeakTarget(-6,true);
    juce::TooltipWindow tip;tip.displayTip({120,120},ui.peakHoverText(peakPoint));imageOf(tip,root.getChildFile("peak-hover.png"));tip.hideTip();
    imageOf(ui,root.getChildFile("over-target-peaks.png"));
    ui.setViewRange(50,80);expect(!ui.sampleExceedsTarget(-.8,4),"covering overview extremum outside visible passage cannot create a false in-view marker");ui.setViewRange(0,512);
    ui.setWaveform(overview,false,{"L","R"});expect(!ui.sampleExceedsTarget(-.8,4),"incomplete waveform cannot claim a target exceedance");ui.setWaveform(overview,true,{"L","R"});
    ui.setPeakTarget(-6,false);expect(!ui.peakVisibility().isEnabled()&&!ui.sampleExceedsTarget(-.8,4),"unverified Current suppresses peak verdict overlay");
    ui.setPeakTarget({},true);expect(!ui.peakVisibility().isEnabled()&&!ui.sampleExceedsTarget(-.8,4)&&ui.peakVisibility().getTooltip().contains("No Limit"),"No Limit disables highlight without inventing a ceiling");
    ui.setPeakTarget(-1,true);expect(!ui.sampleExceedsTarget(-.8,4),"target changes immediately update the numeric threshold");ui.peakVisibility().setToggleState(false,juce::dontSendNotification);
    blend.setValue(-1,juce::sendNotificationSync);imageOf(ui,root.getChildFile("waveform-only.png"));blend.setValue(0,juce::dontSendNotification);imageOf(ui,root.getChildFile("blend-centre.png"));
    auto surround=waveform::Overview::create({512,48000,6},std::make_shared<waveform::Budget>());juce::AudioBuffer<float> six(6,512);six.clear();
    for(int ch=0;ch<6;++ch)six.setSample(ch,80+ch,float(.1*(ch+1)));
    expect(surround->append(six.getArrayOfReadPointers(),6,512,0)&&surround->finish(),"six-channel sample overview fixture");ui.setWaveform(surround,true,{"L","R","C","LFE","Ls","Rs"});
    const auto sixPeaks=ui.rangePeakText(64,128);expect(sixPeaks.contains("LFE: selected sample peak")&&sixPeaks.contains("Rs: selected sample peak")&&sixPeaks.contains("frame 85"),"5.1 selection retains each role and its own sample peak position");
    imageOf(ui,root.getChildFile("surround-blend.png"));
    juce::AudioBuffer<float> intersample(2,5000);
    for(int c=0;c<2;++c)for(int n=0;n<5000;++n)intersample.setSample(c,n,float(.9*std::sin(juce::MathConstants<double>::twoPi*.25*n+juce::MathConstants<double>::pi*.25)));
    measurement::OfflineMeasurement meter(48000,{measurement::Role::left,measurement::Role::right});meter.push(intersample.getArrayOfReadPointers(),5000);const auto& result=meter.finish();
    expect(result.samplePeak.value&&result.truePeak.value&&*result.samplePeak.value< -3&&*result.truePeak.value> -3,"actual measured inter-sample fixture exceeds TP target while PCM samples do not");
    auto intersamplePcm=waveform::Overview::create({5000,48000,2},std::make_shared<waveform::Budget>());expect(intersamplePcm->append(intersample.getArrayOfReadPointers(),2,5000,0)&&intersamplePcm->finish(),"matching true-peak-only overview");
    ui.setProgramme(5000,48000,"true-peak-only");ui.setWaveform(intersamplePcm,true,{"L","R"});ui.setAnalysis(proposal::analyzeProgramme(result,"measured-true-peak-only"));ui.setPeakTarget(-3,true);ui.peakVisibility().setToggleState(true,juce::sendNotificationSync);
    const auto measuredAnalysis=proposal::analyzeProgramme(result,"measured-true-peak-only");const auto& tp=measuredAnalysis->peaks.front();const auto tpBox=ui.graphBounds();
    const juce::Point<float> tpPoint{float(tpBox.getX()+(tp.truePeakEighthFrame/8.)*(tpBox.getWidth()-1)/5000.)+4,float(tpBox.getY()+4)};
    expect(ui.peakHoverText(tpPoint).contains("Measured programme true peak:")&&ui.peakHoverText(tpPoint).contains("/8 frames"),"TP marker hover uses measured true peak and fractional sample position");
    expect(!ui.sampleExceedsTarget(intersample.getSample(0,100),100),"sample markers never pretend to be inter-sample true peaks");imageOf(ui,root.getChildFile("true-peak-only-warning.png"));

}

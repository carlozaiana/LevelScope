#pragma once
namespace
{
void compositionMutedPeaks(const juce::File& root)
{
    for(int owner:{3,5})
    {
        const auto folder=root.getChildFile(juce::String(owner));expect(folder.createDirectory().wasOk(),"muted composition peak fixture");const int frames=192000,peakFrame=frames/2+19;juce::AudioBuffer<float> sum(6,frames);
        for(int ch=0;ch<6;++ch)for(int n=0;n<frames;++n)sum.setSample(ch,n,.008f*static_cast<float>(std::sin(n*(.071+ch*.013))));sum.setSample(owner,peakFrame,.7f);
        for(int iteration=0;iteration<4;++iteration){measurement::OfflineMeasurement m(48000,OfflineAnalysisService::measurementRoles(layoutFor(6)));m.push(sum.getArrayOfReadPointers(),frames);auto r=m.finish();sum.setSample(owner,peakFrame,static_cast<float>(sum.getSample(owner,peakFrame)*std::pow(10.,(*r.integrated.value+22.1-*r.truePeak.value)/20.)));}
        juce::AudioBuffer<float> half(sum);half.applyGain(.5f);ProgrammeComposition recipe;recipe.preset=CompositionPreset::surroundSurround;recipe.me=component(folder,"me",48000,half,false);recipe.dialogue=component(folder,"dialogue",48000,half,false);expect(recipe.confirmComplementaryBalancedAlignment(),"known two-part peak composition confirmed");
        auto s=compositionModel(recipe);measuredComposition(s,MeasurementDomain::source);Capture offline(6,frames);auto current=measuredComposition(s,MeasurementDomain::currentState,&offline);
        expect(std::max_element(current.detail.channelTruePeaks.begin(),current.detail.channelTruePeaks.end())-current.detail.channelTruePeaks.begin()==owner,"composed unique highest true peak belongs to LFE or surround");
        const auto draft=proposal::draftGlobalTrim(s.currentTrimRequest());expect(draft.status==proposal::TrimStatus::draft&&draft.predictedTruePeak&&std::abs(*draft.predictedTruePeak+1)<1e-12&&*draft.deltaDb<-23-*current.detail.integrated.value,"composed LFE or surround peak constrains predicted global trim");
        const auto evaluation=s.evaluation(MeasurementDomain::currentState).identity;const auto source=s.source.analysis,measured=s.currentState.analysis;
        for(int ch=0;ch<6;++ch)s.setOutputRoute(ch,ch==owner?-1:ch);
        expect(s.source.analysis==source&&s.currentState.analysis==measured&&s.evaluation(MeasurementDomain::currentState).identity==evaluation&&proposal::draftGlobalTrim(s.currentTrimRequest()).identity==draft.identity,"muting composed peak owner leaves measured rules and peak-constrained proposal unchanged");
        s.selectAuditionDomain(MeasurementDomain::currentState);auto request=s.auditionRequest();auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));audition.synchronize(request,s.sourceDocument.routing);expect(audition.connect(48000).isEmpty(),"muted peak injected device only");Capture playback(6,frames);audition.setCaptureObserver(&playback);audition.setMonitorGainDb(-30);
        expect(audition.play(request,s.sourceDocument.routing,true).isEmpty(),"explicit partial-monitoring composition play");awaitPlaying(audition,*injected);drain(audition,*injected,1024);expect(identical(offline,playback)&&sameSamples(sum,playback.data),"muted composed LFE/surround remains exact in every pre-monitor sample");audition.disconnect();
    }
}
void compositionFormats(const juce::File& root)
{
    expect(root.createDirectory().wasOk(),"mixed-format composition fixtures");
    for(int channels:{2,6})for(int format=0;format<8;++format)
    {
        const auto folder=root.getChildFile(juce::String(channels)+"-"+juce::String(format));expect(folder.createDirectory().wasOk(),"mixed-format case directory");
        const int bits=format<3?16+8*format:format<6?16+8*(format-3):16+8*(format-6);
        const auto suffix=format<3?".wav":format<6?".aiff":".flac";
        AuditionRequest a;
        if(format==5)
        {
            const auto file=folder.getChildFile("me.aiff");auto out=file.createOutputStream();expect(out!=nullptr,"independent AIFF PCM32 fixture stream");
            const int bytes=48000*channels*4;out->write("FORM",4);out->writeIntBigEndian(bytes+46);out->write("AIFFCOMM",8);out->writeIntBigEndian(18);out->writeShortBigEndian(static_cast<short>(channels));out->writeIntBigEndian(48000);out->writeShortBigEndian(32);out->writeShortBigEndian(16398);out->writeInt64BigEndian(static_cast<juce::int64>(uint64_t(48000)<<48));out->write("SSND",4);out->writeIntBigEndian(bytes+8);out->writeIntBigEndian(0);out->writeIntBigEndian(0);
            for(int n=0;n<48000;++n)for(int ch=0;ch<channels;++ch)out->writeIntBigEndian(static_cast<int32_t>((((n+ch*7)%128)-64)*100000));
            out->flush();expect(out->getStatus().wasOk(),"AIFF PCM32 fixture complete");out.reset();JobControl job;const auto v=validateMedia(file,job);expect(v.succeeded(),"AIFF PCM32 component validated");a=AuditionRequest::source(file,v.media,layoutFor(channels));
        }
        else a=fixture(folder.getChildFile("me"+juce::String(suffix)),48000,channels,48000,0,bits);
        auto b=fixture(folder.getChildFile("dialogue.wav"),48000,channels,48000,0,24);
        ProgrammeComposition recipe;recipe.preset=channels==2?CompositionPreset::stereoStereo:CompositionPreset::surroundSurround;
        for(int index=0;index<2;++index){const auto& request=index==0?a:b;auto& p=index==0?recipe.me:recipe.dialogue;p.id=index==0?"me":"dialogue";p.purpose=index==0?"music-effects":"dialogue";p.file=request.context.sourceFile;p.media=request.context.media;p.layout=request.context.fileLayout;}
        expect(recipe.confirmComplementaryBalancedAlignment(),"qualified component encodings accepted with explicit roles");auto s=compositionModel(recipe);Capture actual(channels,48000);measuredComposition(s,MeasurementDomain::source,&actual);
        juce::AudioFormatManager formats;registerSourceFormats(formats);std::unique_ptr<juce::AudioFormatReader> ra(formats.createReaderFor(a.context.sourceFile)),rb(formats.createReaderFor(b.context.sourceFile));
        juce::AudioBuffer<float> x(channels,48000),y(channels,48000);expect(ra&&rb&&ra->read(x.getArrayOfWritePointers(),channels,0,48000)&&rb->read(y.getArrayOfWritePointers(),channels,0,48000),"independent native decoder operands");
        for(int ch=0;ch<channels;++ch)for(int n=0;n<48000;++n)x.setSample(ch,n,static_cast<float>(double(x.getSample(ch,n))+double(y.getSample(ch,n))));
        expect(sameSamples(x,actual.data),"mixed WAV AIFF FLAC PCM16/24/32 composed samples equal independent decoded sum");
    }
}
void editCompositionCore(const juce::File& from,const juce::File& to,const std::function<void(juce::var&)>& edit)
{
    juce::MemoryBlock bytes;expect(from.loadFileAsData(bytes),"read composition session envelope");
    const auto* data=static_cast<const uint8_t*>(bytes.getData());uint64_t length=0;for(int i=0;i<8;++i)length|=uint64_t(data[8+i])<<(8*i);
    auto core=juce::JSON::parse(juce::String::fromUTF8(reinterpret_cast<const char*>(data+80),static_cast<int>(length)));edit(core);
    const auto json=juce::JSON::toString(core,true).toStdString();const auto hash=juce::SHA256(json.data(),json.size()).toHexString();
    auto out=to.createOutputStream();expect(out!=nullptr,"edited core fixture stream");out->setPosition(0);out->truncate();out->write("LSSESS01",8);out->writeInt64(static_cast<juce::int64>(json.size()));out->write(hash.toRawUTF8(),64);out->write(json.data(),json.size());out->write(data+80+length,bytes.getSize()-80-length);out->flush();expect(out->getStatus().wasOk(),"edited core fixture flushed");
}
void compositionEdges(const juce::File& root)
{
    auto f=fixtureComposition(root,48000,CompositionPreset::surroundSurround,false,48000*4);auto s=compositionModel(f.recipe);JobControl job;
    measuredComposition(s,MeasurementDomain::source);measuredComposition(s,MeasurementDomain::currentState);
    const auto file=root.getChildFile("original.levelscope");expect(sf::save(sf::capture(s),file,job).succeeded(),"composition serialization fixture");
    const auto bad=root.getChildFile("edited.levelscope");
    std::vector<std::function<void(juce::var&)>> corrupt{
        [](auto& v){v["programme"]["composition"].getDynamicObject()->setProperty("engine","unknown");},
        [](auto& v){v["programme"]["composition"].getDynamicObject()->setProperty("preset",3);},
        [](auto& v){v["programme"]["composition"]["dialogue"].getDynamicObject()->setProperty("purpose","analysis-reference");},
        [](auto& v){v["programme"]["composition"]["dialogue"].getDynamicObject()->setProperty("offset","1");},
        [](auto& v){v["programme"]["composition"].getDynamicObject()->setProperty("confirmation",juce::String::repeatedString("a",64));},
        [](auto& v){v["current"]["componentPlans"].getArray()->add("future-DSP");},
        [](auto& v){v["current"].getDynamicObject()->setProperty("wholeProgrammeGraph","future-graph");}
    };
    for(const auto& edit:corrupt){editCompositionCore(file,bad,edit);expect(!sf::read(bad,job).document,"unsupported or mismatched composition core rejected safely");}
    editCompositionCore(file,bad,[](auto& v){v["caches"].getArray()->set(1,"independently damaged");});auto recovered=sf::read(bad,job);
    expect(recovered.document&&recovered.document->source.composition&&recovered.document->sourceCache&&!recovered.document->currentCache,"damaged Current cache preserves composition settings and valid Source cache");
    auto read=sf::read(file,job);auto d=*read.document;auto incompatible=std::make_shared<AnalysisSnapshot>(*d.currentCache);incompatible->context.platformAndFloatingPointMode="other-platform";
    incompatible->identity=OfflineAnalysisService::identityFor(incompatible->context,incompatible->domain,incompatible->measuredContentHash);d.currentCache=incompatible;
    auto verified=sf::verify(d,job);expect(verified.prepared&&!verified.prepared->model.currentState.hasMeasurement&&verified.prepared->model.source.hasMeasurement,"foreign-platform composition cache discarded independently");
    auto duplicate=f.recipe;duplicate.dialogue.file=duplicate.me.file;duplicate.dialogue.media=duplicate.me.media;
    expect(duplicate.confirmComplementaryBalancedAlignment()&&!duplicate.matchesFileStamps()&&!s.setValidatedComposition(duplicate),"duplicate physical component cannot be implicitly double-counted");
    auto permuted=f.recipe;std::swap(permuted.dialogue.layout.roles[0],permuted.dialogue.layout.roles[5]);
    expect(permuted.validate().isNotEmpty()&&permuted.confirmComplementaryBalancedAlignment(),"role change requires fresh exact confirmation");
    auto p=compositionModel(permuted);Capture capture(6,48000*4);measuredComposition(p,MeasurementDomain::source,&capture);
    auto me=ProgrammeReader::open({f.recipe.me.file},f.recipe.me.media),dialogue=ProgrammeReader::open({f.recipe.dialogue.file},f.recipe.dialogue.media);
    juce::AudioBuffer<float> a(6,48000*4),b(6,48000*4);expect(me&&dialogue&&me->read(a.getArrayOfWritePointers(),6,0,a.getNumSamples())&&dialogue->read(b.getArrayOfWritePointers(),6,0,b.getNumSamples()),"independent role permutation operands read");
    for(int c=0;c<6;++c)for(int n=0;n<a.getNumSamples();++n)a.setSample(c,n,static_cast<float>(double(a.getSample(c,n))+double(b.getSample(c==0?5:c==5?0:c,n))));
    expect(sameSamples(a,capture.data),"explicit roles map actual operands without reordering original component identity");
    for(int kind=0;kind<4;++kind)
    {
        auto path=root.getChildFile("numeric-"+juce::String(kind));expect(path.createDirectory().wasOk(),"numeric composition fixture folder");juce::AudioBuffer<float> x(2,48000),y(2,48000);
        for(int c=0;c<2;++c)for(int n=0;n<48000;++n){x.setSample(c,n,kind==0?.25f:std::numeric_limits<float>::max());y.setSample(c,n,kind==0?-.25f:std::numeric_limits<float>::max());}
        if(kind==2)x.setSample(0,99,std::numeric_limits<float>::quiet_NaN());if(kind==3)y.setSample(1,499,std::numeric_limits<float>::infinity());
        if(kind>=2)
        {
            const auto invalidFile=path.getChildFile("non-finite.wav");floatFile(invalidFile,48000,kind==2?x:y);
            expect(validateMedia(invalidFile,job).error==ErrorCode::nonFiniteAudio,"nonfinite component rejected before it can replace confirmed composition");continue;
        }
        ProgrammeComposition recipe;recipe.me=component(path,"me",48000,x,false);recipe.dialogue=component(path,"dialogue",48000,y,false);expect(recipe.confirmComplementaryBalancedAlignment(),"numeric fixture roles explicitly confirmed");
        auto model=compositionModel(recipe);auto result=OfflineAnalysisService().analyze(model.measurementContext(MeasurementDomain::source),MeasurementDomain::source,job);
        if(kind==0)expect(result.value&&!result.value->detail.integrated.value&&!result.value->detail.truePeak.value,"exact phase cancellation remains unavailable silence, never a measured zero loudness");
        else expect(!result.value,"nonfinite operand or float32 sum overflow fails measurement without partial success");
    }
    const auto context=s.measurementContext(MeasurementDomain::source);SnapshotStorage storage;storage.directory=root;
    struct FlushFailure:SnapshotWriter{bool write(const void*,size_t)override{return true;}bool flush()override{return false;}};
    storage.createWriter=[](const auto&){return std::make_unique<FlushFailure>();};
    expect(prepareVerifiedSource(context,job,storage).error==ErrorCode::writeFailed,"component flush failure rejects entire verified composition");
    JobControl cancel;auto worker=std::async(std::launch::async,[&]{return prepareVerifiedSource(context,cancel);});
    while(cancel.progress.load()==0&&worker.wait_for(std::chrono::milliseconds(1))!=std::future_status::ready){}
    cancel.cancel=true;const auto stopped=worker.get();expect(stopped.error==ErrorCode::cancelled,"active composition copy cancellation returns no partial session");
}
}

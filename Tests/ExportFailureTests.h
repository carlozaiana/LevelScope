#pragma once
namespace
{
void replaceByte(const juce::File& file,int64_t pos,uint8_t value)
{const int fd=::open(file.getFullPathName().toRawUTF8(),O_WRONLY);expect(fd>=0&&::pwrite(fd,&value,1,pos)==1,"controlled owned-byte tamper");::close(fd);}
void failureTests()
{
    Temp temp;auto s=measured(temp.directory,2,48000);const auto request=ex::capture(s);JobControl hashControl;
    const auto original=pcm::hashFile(request.current.sourceFile,hashControl);const auto source=s.source.analysis,current=s.currentState.analysis;
    for(const char* stage:{"prepare","audio-create","audio-header","render-block","audio-write","audio-flush","audio-close","verify","verified"})
    for(bool cancel:{false,true})
    {
        JobControl control;pcm::Hooks hook;
        if(cancel)hook.reached=[&](const char* name){if(juce::String(name)==stage)control.cancel=true;};
        else hook.fail=[&](const char* name){return juce::String(name)==stage;};
        auto result=ex::render(request,temp.directory,control,hook);
        expect(!result.value&&(!cancel||result.error==ErrorCode::cancelled),"render stage failure/cancellation publishes no candidate");
        expect(temp.directory.findChildFiles(juce::File::findDirectories,false,".levelscope-export-*").isEmpty(),"failed render cleans owned staging");
    }
    JobControl c;pcm::Hooks noSpace;noSpace.availableBytes=[](const auto&){return int64_t(0);};
    expect(ex::render(request,temp.directory,c,noSpace).error==ErrorCode::insufficientDiskSpace,"destination space fail before render");
    ex::Specification needed;expect(ex::specification(request,needed).isEmpty(),"space-test export specification");
    pcm::Hooks combined;combined.availableBytes=[&](const juce::File&){return static_cast<int64_t>(needed.format.fileBytes+80ULL*1024*1024+request.current.media.bytes/2);};
    expect(ex::render(request,temp.directory,c,combined).error==ErrorCode::insufficientDiskSpace,"same-volume source-copy plus output disk budget checked together");
    auto rendered=ex::render(request,temp.directory,c);expect(rendered.value.has_value(),"publication failure fixture verified");auto candidate=*rendered.value;
    const auto approval=ex::approvalIdentity(*candidate,request);auto destination=temp.directory.getChildFile("published");
    for(const char* stage:{"publish-recheck","report-json","report-html","report-manifest","report-write","report-flush","report-close","publish-rename"})for(bool cancel:{false,true})
    {
        JobControl control;pcm::Hooks hook;
        if(cancel)hook.reached=[&](const char* name){if(juce::String(name)==stage)control.cancel=true;};else hook.fail=[&](const char* name){return juce::String(name)==stage;};
        auto result=ex::publish(*candidate,request,destination,approval,control,hook);
        expect(!result.value&&!destination.exists()&&candidate->file().existsAsFile(),"publication failure/cancel retains verified audio without final set");
        expect(!candidate->folder().getChildFile("report.json").exists(),"failed publication leaves no misleading report");
    }
    for(const char* name:{"report.json","report.html","manifest.json","programme.wav"})
    {
        JobControl control;pcm::Hooks hook;hook.reached=[&](const char* stage)
        {if(juce::String(stage)=="publish-rename")replaceByte(candidate->folder().getChildFile(name),0,0x7f);};
        auto result=ex::publish(*candidate,request,destination,approval,control,hook);
        expect(!result.value&&!destination.exists(),"last-stage audio/report/manifest tamper blocks publication");
    }
    candidate.reset();rendered.value.reset();rendered=ex::render(request,temp.directory,c);expect(rendered.value.has_value(),"fresh measured artifact after tamper");candidate=*rendered.value;
    pcm::Hooks race;race.reached=[&](const char* stage){if(juce::String(stage)=="publish-rename")destination.createDirectory();};
    expect(!ex::publish(*candidate,request,destination,ex::approvalIdentity(*candidate,request),c,race).value&&destination.isDirectory(),"exclusive destination race never overwrites");destination.deleteFile();
    const auto symlink=temp.directory.getChildFile("source-alias"),hardlink=temp.directory.getChildFile("source-hardlink");
    expect(::symlink(request.current.sourceFile.getFullPathName().toRawUTF8(),symlink.getFullPathName().toRawUTF8())==0
        &&::link(request.current.sourceFile.getFullPathName().toRawUTF8(),hardlink.getFullPathName().toRawUTF8())==0,"owned source alias fixtures");
    for(const auto& alias:{request.current.sourceFile,symlink,hardlink})
        expect(!ex::publish(*candidate,request,alias,ex::approvalIdentity(*candidate,request),c).value
            &&!ex::measurementReport(request,alias,c).value,"audio/report destinations never overwrite a source or alias");
    expect(pcm::hashFile(request.current.sourceFile,c)==original,"rejected aliases preserve source bytes");
    auto stale=request;stale.current.processing.staticGainDb+=.01;
    expect(!candidate->matchesAudio(stale)&&ex::approvalIdentity(*candidate,stale).isEmpty(),"processing change invalidates render and approval");
    stale=request;stale.current.sessionId="another-opening";expect(!candidate->matchesAudio(stale),"runtime reopening cannot inherit approval");
    stale=request;stale.applicability=proposal::Applicability::unconfirmed;
    expect(candidate->matchesAudio(stale)&&ex::approvalIdentity(*candidate,stale).isEmpty(),"applicability reevaluates unchanged artifact without rerender");
    stale=request;stale.current.profile.definitionVersion="test-compatible-revision";
    expect(candidate->matchesAudio(stale)&&ex::artifactEvaluation(*candidate,stale).verdict==proposal::Verdict::pass&&ex::approvalIdentity(*candidate,stale)!=ex::approvalIdentity(*candidate,request),"compatible profile version reevaluates the same measured artifact and refreshes approval");
    stale=request;stale.documentRevision++;stale.video.startOnAudio=video::Time(-1001,30000);
    expect(candidate->matchesAudio(stale)&&ex::approvalIdentity(*candidate,stale)!=ex::approvalIdentity(*candidate,request),"video/report revision changes approval but preserves measured audio");
    expect(!ex::publish(*candidate,stale,destination,ex::approvalIdentity(*candidate,request),c).value,"obsolete report approval rejected");
    expect(pcm::hashFile(request.current.sourceFile,hashControl)==original,"render/publication failures preserve original source bytes");
    const auto member=request.current.sourceFile;replaceByte(member,member.getSize()-1,0x7f);
    expect(!ex::publish(*candidate,request,destination,ex::approvalIdentity(*candidate,request),c).value&&!destination.exists(),"changed original content blocks publication");
    expect(s.source.analysis==source&&s.currentState.analysis==current,"export failures never modify session measurements or processing");
    // A complete, deliberately altered encoded signal must be evaluated on its actual bytes.
    Temp failTemp;auto failingModel=measured(failTemp.directory,2,48000);const auto failRequest=ex::capture(failingModel);pcm::Hooks alter;
    alter.reached=[&](const char* stage){if(juce::String(stage)=="verify")
    {const auto dirs=failTemp.directory.findChildFiles(juce::File::findDirectories,false,".levelscope-export-*");expect(dirs.size()==1,"identify owned fault injection");const auto f=dirs[0].getChildFile("programme.wav");auto stream=f.createOutputStream();stream->setPosition(68);std::vector<char> zeros(static_cast<size_t>(f.getSize()-68),0);expect(stream->write(zeros.data(),zeros.size()),"zero encoded payload fault");}};
    auto failure=ex::render(failRequest,failTemp.directory,c,alter);expect(failure.value.has_value(),"complete failing output still receives actual measurement");
    expect(ex::artifactEvaluation(**failure.value,failRequest).verdict!=proposal::Verdict::pass&&ex::approvalIdentity(**failure.value,failRequest).isEmpty(),"passing Current never overrides actual artifact failure");
    expect(ex::report(failRequest,failure.value->get())["artifact"]["metrics"]["integrated"]["value"].isVoid(),"silent artifact unavailable loudness not invented zero");
    expect(failingModel.currentStateDocument.staticGainDb==failRequest.current.processing.staticGainDb,"failed verification retains exact gain, no retry/correction");
    const auto failedReport=ex::artifactReport(**failure.value,failRequest,failTemp.directory.getChildFile("actual-failure-report"),c);
    expect(failedReport.value.has_value()&&!failedReport.value->getChildFile("programme.wav").exists(),"actual failing artifact report can be saved explicitly without publishing audio");
    failure.value.reset();pcm::Hooks half;
    half.reached=[&](const char* stage){if(juce::String(stage)=="verify")
    {
        const auto dirs=failTemp.directory.findChildFiles(juce::File::findDirectories,false,".levelscope-export-*");expect(dirs.size()==1,"identify measured-failure test staging");const auto file=dirs[0].getChildFile("programme.wav");juce::MemoryBlock bytes;expect(file.loadFileAsData(bytes),"read complete injected artifact");auto* data=static_cast<uint8_t*>(bytes.getData());
        for(size_t pos=68;pos+2<bytes.getSize();pos+=3){const auto value=static_cast<uint32_t>(pcmCode(data+pos)/2);data[pos]=uint8_t(value);data[pos+1]=uint8_t(value>>8);data[pos+2]=uint8_t(value>>16);}
        expect(file.replaceWithData(bytes.getData(),bytes.getSize()),"inject changed completed PCM signal");
    }};
    auto changed=ex::render(failRequest,failTemp.directory,c,half);
    expect(changed.value&&ex::artifactEvaluation(**changed.value,failRequest).verdict==proposal::Verdict::fail&&ex::approvalIdentity(**changed.value,failRequest).isEmpty(),"actual measured artifact FAIL overrides passing Current; publication unavailable");
    const auto noArtifact=ex::report(failRequest);expect(!noArtifact.hasProperty("artifact"),"measurement report independent of failing artifact");
    // Reports must remain useful before rendering, for invalid/unavailable/failing Current.
    for(int kind=0;kind<4;++kind)
    {
        auto r=failRequest;if(kind==0)r.currentAnalysis.reset();if(kind==1){r.current.processing.staticGainDb-=1;}
        if(kind==2)r.applicability=proposal::Applicability::unconfirmed;if(kind==3){r.sourceAnalysis.reset();r.currentAnalysis.reset();r.source={};r.current={};}
        const auto d=failTemp.directory.getChildFile("measurement-"+juce::String(kind));
        expect(ex::readiness(r).isNotEmpty()&&ex::measurementReport(r,d,c).value.has_value(),"blocked audio does not block bounded measurement-only report");
        expect(!juce::JSON::parse(d.getChildFile("report.json")).hasProperty("artifact")&&!d.getChildFile("programme.wav").exists(),"report-only has no invented artifact/hash/verdict");
    }
    Temp unsupported;auto wide=measured(unsupported.directory,8,48000);const auto wideRequest=ex::capture(wide);
    expect(ex::readiness(wideRequest).isNotEmpty()&&wide.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass,"qualified 7.1 Current can pass measurement while export is explicitly unsupported");
    auto wideReport=ex::measurementReport(wideRequest,unsupported.directory.getChildFile("report-only"),c);
    expect(wideReport.value&&!juce::JSON::parse(wideReport.value->getChildFile("report.json")).hasProperty("artifact"),"unsupported export layout retains measurement-only reporting");
    const auto reportDest=failTemp.directory.getChildFile("failed-report");pcm::Hooks reportFailure;reportFailure.fail=[](const char* name){return juce::String(name)=="report-html";};
    expect(!ex::measurementReport(failRequest,reportDest,c,reportFailure).value&&!reportDest.exists(),"report-only write failure exposes no partial final set");
    auto escaped=failRequest;escaped.documentId="<script>&\"';";const auto json=ex::reportJson(ex::report(escaped));const auto html=ex::reportHtml(ex::report(escaped));
    expect(!html.contains("<script>")&&html.contains("&lt;script&gt;")&&juce::JSON::parse(json)["document_id"].toString()==escaped.documentId,"report HTML escapes user text; JSON roundtrips");
    expect(static_cast<double>(juce::JSON::parse(json)["processing_gain_db"])==escaped.current.processing.staticGainDb,"report gain retains full double precision");
}
void mappingAndParity(const juce::File& root)
{
    for(int channels:{2,6})for(int rate:{44100,48000,96000,192000})
    {
        Temp temp;auto s=measured(temp.directory,channels,rate);
        // A legal non-canonical input order must be permuted explicitly, never by count.
        auto roles=s.sourceDocument.layout;std::reverse(roles.roles.begin(),roles.roles.begin()+channels);
        for(int ch=0;ch<channels;++ch)s.setChannelRole(ch,roles.roles[ch]);
        expect(s.confirmLayout()&&s.initializeCurrentStateFromSource(),"noncanonical roles explicitly confirmed");
        measure(s,MeasurementDomain::source);s.setCurrentGainDb(-23-*s.source.integratedLufs);
        Capture measuredAudio(channels,rate*4);measure(s,MeasurementDomain::currentState,&measuredAudio);
        const auto request=ex::capture(s);Capture renderedAudio(channels,rate*4);JobControl control;
        auto rendered=ex::render(request,root,control,{},&renderedAudio);
        if(!rendered.value)std::cerr<<rendered.message<<'\n';expect(rendered.value.has_value()&&identical(measuredAudio,renderedAudio),"export pre-quantization exact offline sample parity at every channel/rate");
        auto candidate=*rendered.value;bool reverse=true;for(int ch=0;ch<channels;++ch)reverse&=candidate->spec().outputToInput[ch]==channels-1-ch;expect(reverse,"explicit canonical output permutation preserves input role order");
        auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));s.selectAuditionDomain(MeasurementDomain::currentState);auto req=s.auditionRequest();OutputRouting muted;Capture heard(channels,rate*4);audition.setCaptureObserver(&heard);audition.synchronize(req,muted);expect(audition.connect(rate).isEmpty(),"injected device only for export parity");expect(audition.play(req,muted,true,true).isEmpty(),"injected muted Current play");awaitPlaying(audition,*injected);drain(audition,*injected,997);audition.disconnect();
        expect(identical(heard,renderedAudio)&&ex::audioIdentity(ex::capture(s))==ex::audioIdentity(request),"muted all-channel audition equals render; monitor never alters export identity");
        const auto dest=root.getChildFile("verified-"+juce::String(channels)+"-"+juce::String(rate));
        expect(ex::publish(*candidate,request,dest,ex::approvalIdentity(*candidate,request),control).value.has_value(),"mapped real PCM independently measured and published");
        const auto parsed=juce::JSON::parse(dest.getChildFile("report.json"));expect(parsed["artifact"]["sha256"].toString()==pcm::hashFile(dest.getChildFile("programme.wav"),control),"published report binds exact actual PCM bytes");
        if(channels==2&&rate==48000)
        {
            namespace sf=levelscope::standalone::sessionfile;
            const auto imported=validateMedia(dest.getChildFile("programme.wav"),control);StandaloneSessionModel ui;
            expect(imported.succeeded()&&ui.setValidatedSource(dest.getChildFile("programme.wav"),imported.media)&&ui.confirmLayout()&&ui.initializeCurrentStateFromSource(),"reimport exact stereo output as a new explicit programme");
            ui.setSelectedTargetProfileIndex(1);ui.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
            measure(ui,MeasurementDomain::source);measure(ui,MeasurementDomain::currentState);
            expect(sf::save(sf::capture(ui),root.getChildFile("ui-ready.levelscope"),control).succeeded(),"coherent same-build UI fixture has measured Source/Current at zero gain");
        }
        auto artifactRow=parsed["artifact"];artifactRow.getDynamicObject()->setProperty("input_kind","interleaved");cases.add(artifactRow);
        // Assemble identical decoded channels as a separate-file programme.
        auto reader=ProgrammeReader::open({request.current.sourceFile},request.current.media);juce::AudioBuffer<float> originals(channels,rate*4);
        expect(reader&&reader->read(originals.getArrayOfWritePointers(),channels,0,rate*4),"read original group equivalence fixture");std::vector<juce::File> files;
        for(int ch=0;ch<channels;++ch)
        {
            auto file=temp.directory.getChildFile("mono-"+juce::String(ch)+".wav");auto stream=file.createOutputStream();juce::WavAudioFormat wav;
            std::unique_ptr<juce::AudioFormatWriter> writer(wav.createWriterFor(stream.get(),rate,1,24,{},0));expect(writer!=nullptr,"group equivalent mono writer");stream.release();
            juce::AudioBuffer<float> mono(1,rate*4);mono.copyFrom(0,0,originals,ch,0,rate*4);expect(writer->writeFromAudioSampleBuffer(mono,0,rate*4),"group equivalent mono payload");writer.reset();files.push_back(file);
        }
        auto group=inspectProgrammeFiles(files,control);expect(group.error==ErrorCode::none,"equivalent group verified completely");
        for(int ch=0;ch<channels;++ch)group.group.members[ch].role=roles.roles[ch];expect(group.group.orderByRoles()&&group.group.confirmAlignment(),"group exact mapping and alignment explicitly confirmed");
        StandaloneSessionModel grouped;expect(grouped.setValidatedProgramme(group.group)&&grouped.initializeCurrentStateFromSource(),"grouped export session");grouped.setSelectedTargetProfileIndex(1);grouped.setProgrammeApplicability(proposal::Applicability::preProducedLinearProgramme);
        measure(grouped,MeasurementDomain::source);grouped.setCurrentGainDb(request.current.processing.staticGainDb);Capture groupOffline(channels,rate*4);measure(grouped,MeasurementDomain::currentState,&groupOffline);
        Capture groupRender(channels,rate*4);const auto gr=ex::capture(grouped);auto groupCandidate=ex::render(gr,root,control,{},&groupRender);expect(groupCandidate.value&&identical(groupOffline,groupRender),"group export exact offline pre-quantization parity");
        bool same=true;for(int ch=0;ch<channels;++ch)same&=std::memcmp(groupRender.data.getReadPointer(ch),renderedAudio.data.getReadPointer(channels-1-ch),size_t(rate*4)*sizeof(float))==0;
        expect(same,"equivalent group/interleaved processing samples match in canonical speaker order");
        const auto groupDest=root.getChildFile("verified-group-"+juce::String(channels)+"-"+juce::String(rate));expect(ex::publish(**groupCandidate.value,gr,groupDest,ex::approvalIdentity(**groupCandidate.value,gr),control).value.has_value(),"group stereo/5.1 export published after actual verification");
        auto groupRow=juce::JSON::parse(groupDest.getChildFile("report.json"))["artifact"];groupRow.getDynamicObject()->setProperty("input_kind","mono-group");cases.add(groupRow);

    }
}
}

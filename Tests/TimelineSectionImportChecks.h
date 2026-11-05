#include "Standalone/StandaloneAudioPanel.h"
#include "Standalone/StandaloneProgrammeImport.h"
#include "Standalone/SessionDocument.h"
#include "Standalone/Pcm24Wave.h"

namespace
{
void pump(int ms=50){juce::MessageManager::getInstance()->runDispatchLoopUntil(ms);}
void presetDismissal(const juce::File& root)
{
    StandaloneSessionModel s;
    const auto folder=root.getChildFile("Target presets");
    JobRequirementsPanel panel(s,[]{},folder);
    for(bool save:{true,false})
    {
        button(panel,"Save as preset...")->onClick();
        auto* alert=dynamic_cast<juce::AlertWindow*>(juce::Component::getCurrentlyModalComponent());
        expect(alert!=nullptr,"preset naming dialog opened");
        juce::Component::SafePointer<juce::AlertWindow> lifetime(alert);
        alert->getTextEditor("name")->setText(save?"Delivery A":"Cancelled B");
        button(*alert,save?"Save":"Cancel")->onClick();
        pump();
        expect(!lifetime && juce::Component::getCurrentlyModalComponent()==nullptr,"Save and Cancel both destroy preset dialog and release modality");
        expect(folder.getNumberOfChildFiles(juce::File::findFiles,"*.levelscope-target")==1 && s.jobRequirements.name=="Delivery A","Save persisted one named preset and Cancel left it intact");
    }
}
void draw(MeasurementTimeline& view,float start,float end)
{
    const auto y=float(view.graphBounds().getCentreY());
    view.mouseDown(mouse(view,{start,y}));
    view.mouseDrag(mouse(view,{end,y}));
    view.mouseUp(mouse(view,{end,y},false));
}
void drawSections(const juce::File& root)
{
    auto mix=fixture(root.getChildFile("draw.wav"),48000,2,48000);
    StandaloneSessionModel s;
    expect(s.setValidatedSource(mix.context.sourceFile,mix.context.media),"drawing fixture imported");
    s.sourceDocument.layout=mix.context.fileLayout;s.synchronizeSections();
    measure(s,Domain::source);
    const auto original=s.source.analysis,current=s.currentState.analysis;
    const auto audioRevision=s.revision,docRevision=s.documentRevision;
    MeasurementTimeline main;
    uint64_t lastSelected=0;
    main.bindSectionEditing(s,[&](uint64_t id){lastSelected=id;});
    StandaloneSectionsPanel sections(s);
    auto& secondary=sections.timelineView();
    WaveformProgress pcm(std::make_shared<waveform::Budget>());JobControl control;
    expect(OfflineAnalysisService(nullptr,&pcm).buildWaveform(s.measurementContext(Domain::source),Domain::source,control).value.has_value(),"real waveform for section gestures");
    for(auto* view:{&main,&secondary})
    {
        view->setSize(1146,400);view->setProgramme(48000,48000,"draw");view->setWaveform(pcm.overview(),true,{"L","R"});view->refreshSectionEditing();
    }
    const auto x=float(main.graphBounds().getX());
    draw(main,x+500,x+200);
    const auto& map=s.sectionDocument.map();
    expect(map.sections.size()==3&&map.sections[1].id==lastSelected&&s.documentRevision==docRevision+1,"reverse waveform drag commits one document edit at release");
    const auto width=main.graphBounds().getWidth()-1;
    expect(map.sections[1].begin==200*48000/width && map.sections[1].end==500*48000/width,"selected range snaps to exact source frames");
    secondary.refreshSectionEditing();
    imageOf(main,root.getChildFile("drawn-section.png"));
    draw(secondary,x+500,x+560);
    expect(map.sections[1].end==560*48000/width&&map.sections[2].begin==map.sections[1].end,"dragging the same boundary in Sections changes both neighbours");
    draw(secondary,x+620,x+780);
    expect(map.sections.size()==5,"Sections waveform can define another section directly");
    main.refreshSectionEditing();const auto digest=map.digest();
    draw(main,x+100,x+650);
    expect(map.digest()==digest&&main.navigationStatus().contains("one existing section"),"cross-boundary draw preserves named sections");
    const auto y=float(main.graphBounds().getCentreY());
    main.mouseDown(mouse(main,{x+20,y}));
    main.mouseDrag(mouse(main,{x+140,y}));
    juce::String error;expect(s.sectionDocument.rename(map.sections.front().id,"Edited elsewhere",error),"external annotation update");
    main.mouseUp(mouse(main,{x+140,y},false));
    expect(map.sections.size()==5,"in-flight drag with stale annotation revision rejected");
    expect(s.revision==audioRevision&&s.source.analysis==original&&s.currentState.analysis==current,"drawing and moving boundaries leaves accepted measurements and audio revision intact");

    auto io=std::make_unique<Device>();auto* injected=io.get();StandaloneAudition transport(std::move(io));
    s.selectAuditionDomain(Domain::source);auto req=s.auditionRequest();auto routing=routes(2);transport.synchronize(req,routing);
    secondary.onRequestSeek=[&](int64_t f,const juce::String& id){return TimelineTransport::seek(f,id,s,transport);};
    auto refreshTransport=[&]{auto why=TimelineTransport::blockedReason(s,transport);secondary.setTransport(transport.positionFrames(),why.isEmpty(),why,transport.boundIdentity());};
    refreshTransport();auto point=secondary.rulerBounds().getCentre().toFloat();
    click(secondary,point);
    expect(transport.positionFrames()>0 && !transport.active() && injected->opens==0,"Sections ruler seeks with device closed and transport stopped");
    expect(transport.connect(48000).isEmpty(),"connect injected Sections output");transport.synchronize(req,routing);
    sections.onPlay=[&]{expect(transport.play(req,routing).isEmpty(),"Sections Play starts shared audition");refreshTransport();};
    sections.onStop=[&]{transport.stop();refreshTransport();};
    secondary.onTogglePlayback=[&]{if(transport.active())sections.onStop();else sections.onPlay();};
    sections.setPlaybackState(true,false,"Original");
    button(sections,"Play")->onClick();awaitPlaying(transport,*injected);
    click(secondary,point.translated(20,0));
    expect(transport.getState()==AuditionState::playing,"ruler while playing preserves shared audition");
    sections.setPlaybackState(false,true,"Original");button(sections,"Stop")->onClick();
    expect(!transport.active(),"Sections Stop stops shared transport");
    sections.keyPressed(juce::KeyPress(juce::KeyPress::spaceKey));awaitPlaying(transport,*injected);
    sections.keyPressed(juce::KeyPress(juce::KeyPress::spaceKey));
    expect(!transport.active(),"Sections Space toggles playback and stop");
    transport.disconnect();sections.suspend();
}
MediaInfo labelledWave(const juce::File& file,int channels)
{
    pcm24::Format format;juce::String error;
    expect(pcm24::describe(48000,channels,4800,format,error),"describe labelled PCM");
    auto data=pcm24::header(format);
    data.resize(size_t(format.fileBytes),0);
    expect(file.replaceWithData(data.data(),data.size()),"write WAVEFORMATEXTENSIBLE fixture");
    JobControl verify;auto media=validateMedia(file,verify);
    expect(media.succeeded() && media.media.hasFileRoleLabels,"actual speaker mask recognized as reliable file labels");
    return media.media;
}
void rememberedOrders(const juce::File& root)
{
    using C=juce::AudioChannelSet;
    const auto folder=root.getChildFile("roles");
    RememberedChannelRoles memory(folder);
    auto stereo=fixture(root.getChildFile("unlabelled.wav"),48000,2,4800);
    auto six=fixture(root.getChildFile("unlabelled.aiff"),48000,6,4800);
    expect(!stereo.context.media.hasFileRoleLabels&&!six.context.media.hasFileRoleLabels,"unlabelled formats have no asserted role order");
    const auto initial=memory.suggest(stereo.context.media);
    expect(!initial.order.confirmed && initial.explanation.contains("standard"),"count proposes a reviewable order only");
    StandaloneSessionModel s;
    auto device=std::make_unique<Device>();auto* injected=device.get();StandaloneAudition audition(std::move(device));
    StandaloneAudioPanel audio(s,audition,[]{},folder,emptyOutputCatalog(),noOutputSettings());int analyses=0;
    audio.onAnalyse=[&]{expect(s.sourceDocument.layout.confirmed,"analysis starts only after confirmation");++analyses;};
    auto pick=[&](const juce::String& name,ChannelRole role)
    {
        auto find=[&](auto&& self,juce::Component& parent)->juce::ComboBox*
        {
            for(auto* c:parent.getChildren())
            {
                if(c->getName()==name)return dynamic_cast<juce::ComboBox*>(c);
                if(auto* match=self(self,*c))return match;
            }
            return nullptr;
        };
        auto* box=find(find,audio);expect(box&&box->isEnabled(),"Edit order enables the selected file-channel editor");
        const auto index=std::find(supportedRoles.begin(),supportedRoles.end(),role)-supportedRoles.begin();
        box->setSelectedId(int(index)+2,juce::sendNotificationSync);
    };
    auto import=[&](const AuditionRequest& request)
    {
        expect(s.setValidatedSource(request.context.sourceFile,request.context.media),"fresh clip accepted");
        audio.sourceChanged(true);
        expect(!s.sourceDocument.layout.confirmed&&!AutomaticSourceAnalysis().eligible(s),"proposal alone cannot start analysis");
    };
    import(stereo);
    button(audio,"Edit order")->onClick();
    pick("Role for file channel 1",C::right);pick("Role for file channel 2",C::left);audio.refresh();
    button(audio,"Confirm roles and analyse")->onClick();
    expect(analyses==1&&memory.suggest(stereo.context.media).order.roles[0]==C::right,"one click confirms custom stereo order, remembers it and requests analysis");
    const auto project=root.getChildFile("channel-roles.levelscope");JobControl control;
    expect(sessionfile::save(sessionfile::capture(s),project,control).succeeded(),"save reversed stereo project");
    import(six);button(audio,"Edit order")->onClick();
    pick("Role for file channel 1",C::centre);pick("Role for file channel 3",C::left);audio.refresh();button(audio,"Confirm roles and analyse")->onClick();
    expect(analyses==2&&memory.suggest(six.context.media).order.roles[0]==C::centre,"5.1 has a separate explicitly confirmed order");
    import(stereo);expect(s.sourceDocument.layout.roles[0]==C::right,"alternate back to stereo recalls its own order");
    import(six);expect(s.sourceDocument.layout.roles[0]==C::centre,"alternate back to 5.1 recalls its own order");
    const auto labelledFile=root.getChildFile("declared-stereo.wav");
    const auto labelled=labelledWave(labelledFile,2);const auto conflict=memory.suggest(labelled);
    expect(conflict.order.roles[0]==C::left&&!conflict.order.confirmed&&conflict.explanation.contains("Differs from remembered")&&conflict.explanation.contains("R / L"),"file mask wins conflict and both orders are visible before confirmation");
    const auto labelledSix=labelledWave(root.getChildFile("declared-surround.wav"),6);
    expect(memory.suggest(labelledSix).order.roles[0]==C::left && memory.suggest(labelledSix).explanation.contains("Differs"),"5.1 metadata also takes priority");
    const auto native=fixture(root.getChildFile("native.flac"),48000,6,4800);
    expect(native.context.media.hasFileRoleLabels&&memory.suggest(native.context.media).order.roles[0]==C::left,"native FLAC 5.1 convention precedes remembered order");
    #if JUCE_MAC
    const auto aiff=root.getChildFile("labels.aiff");
    fixture(aiff,48000,2,4800);
    juce::MemoryBlock original;expect(aiff.loadFileAsData(original),"read owned AIFF fixture");
    juce::MemoryOutputStream output;
    output.write(original.getData(),4);
    output.writeIntBigEndian(int(original.getSize()+60-8));
    output.write(static_cast<const char*>(original.getData())+8,original.getSize()-8);
    output.write("CHAN",4);output.writeIntBigEndian(52);output.writeIntBigEndian(0);output.writeIntBigEndian(0);output.writeIntBigEndian(2);
    for(int label:{2,1}){output.writeIntBigEndian(label);for(int j=0;j<4;++j)output.writeIntBigEndian(0);}
    expect(aiff.replaceWithData(output.getData(),output.getDataSize()),"write explicit ordered AIFF labels");
    auto aiffResult=validateMedia(aiff,control);
    expect(aiffResult.succeeded()&&aiffResult.media.hasFileRoleLabels&&aiffResult.media.suggestedLayout.roles[0]==C::right,"AIFF CHAN description order survives decoding");
    #endif
    expect(s.setValidatedSource(labelledFile,labelled),"adopt labelled file");audio.sourceChanged(true);
    audio.setSize(700,600);imageOf(audio,root.getChildFile("role-conflict.png"));
    button(audio,"Confirm roles and analyse")->onClick();
    auto reopened=sessionfile::read(project,control);expect(bool(reopened.document),"read saved assignment");
    auto verified=sessionfile::verify(*reopened.document,control);juce::String error;
    expect(verified.document&&sessionfile::install(*verified.document,s,error),"open saved project against verified media");
    audio.sourceChanged();
    expect(s.sourceDocument.layout.confirmed&&s.sourceDocument.layout.roles[0]==C::right&&memory.suggest(stereo.context.media).order.roles[0]==C::left,
           "project retains its saved roles instead of importing or overwriting current Mac preferences");
    expect(injected->opens==0&&s.sourceDocument.routing.fileToDevice[0]==-1,"roles never reconnect or assign physical outputs");

    ProgrammeImportPanel* panelPtr=nullptr;bool accepted=false;ChannelLayout displayed;
    ProgrammeImportPanel panel([&](std::optional<ProgrammeAssembly> group)
    {
        expect(group&&group->validate().isEmpty()&&!group->alignmentConfirmation.isEmpty(),"complete canonical assembly has explicit alignment");
        displayed=panelPtr->confirmedFileOrder();accepted=true;
    },folder);
    panelPtr=&panel;
    juce::Array<juce::File> mono;
    for(int n=0;n<2;++n){auto file=root.getChildFile("mono-"+juce::String(n)+".wav");fixture(file,48000,1,4800);mono.add(file);}
    panel.inspectFiles(mono);auto* edit=button(panel,"Edit order");
    for(int i=0;i<100 && !edit->isEnabled();++i)pump(25);
    auto* confirm=button(panel,"Confirm roles and analyse");
    juce::String why;
    for(auto* child:panel.getChildren())if(auto* label=dynamic_cast<juce::Label*>(child))why+=" | "+label->getText();
    if(!edit->isEnabled())std::cerr<<why<<"\n";
    expect(edit->isEnabled()&&!confirm->isEnabled(),"separate-file suggestion still needs alignment");
    juce::ToggleButton* aligned=nullptr;
    for(auto* component:panel.getChildren())if(auto* toggle=dynamic_cast<juce::ToggleButton*>(component))aligned=toggle;
    expect(aligned && aligned->isEnabled(),"alignment control available after group inspection");
    imageOf(panel,root.getChildFile("separate-role-review.png"));
    aligned->setToggleState(true,juce::dontSendNotification);aligned->onClick();
    expect(confirm->isEnabled(),"explicit alignment permits confirmation");confirm->onClick();
    expect(accepted && displayed.confirmed && displayed.roles[0]==C::left,"separate files confirm the displayed file-slot order");
}
}

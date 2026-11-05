namespace
{
juce::Component* namedChild(juce::Component& c,const juce::String& name)
{if(c.getName()==name)return &c;for(auto* child:c.getChildren())if(auto* found=namedChild(*child,name))return found;return nullptr;}
juce::Slider& monitorSlider(StandaloneAudioPanel& ui)
{for(auto* c:ui.getChildren())if(auto* s=dynamic_cast<juce::Slider*>(c))return *s;throw std::runtime_error("monitor slider");}
juce::ComboBox& routeBox(StandaloneAudioPanel& ui,int channel)
{auto* box=dynamic_cast<juce::ComboBox*>(namedChild(ui,"Output for file channel "+juce::String(channel+1)));expect(box!=nullptr,"find named output route control");return *box;}
void setRoute(StandaloneAudioPanel& ui,int channel,int output)
{auto& box=routeBox(ui,channel);expect(box.isEnabled(),"routing edit available for confirmed roles");box.setSelectedId(output+2,juce::sendNotificationSync);}
juce::String monitoringText(StandaloneAudioPanel& ui)
{return dynamic_cast<juce::Label*>(ui.findChildWithID("monitoring-memory-status"))->getText();}
void chooseOutput(StandaloneAudioPanel& ui,int number)
{auto* box=dynamic_cast<juce::ComboBox*>(namedChild(ui,"Output interface"));expect(box!=nullptr,"output chooser exists");box->setSelectedId(number,juce::sendNotificationSync);}
bool muted(const StandaloneSessionModel& s)
{return std::all_of(s.sourceDocument.routing.fileToDevice.begin(),s.sourceDocument.routing.fileToDevice.end(),[](int n){return n==-1;});}
void adoptInput(StandaloneSessionModel& s,const AuditionRequest& input,bool confirmed=true)
{expect(s.setValidatedSource(input.context.sourceFile,input.context.media),"adopt new verified programme");s.sourceDocument.layout=input.context.fileLayout;s.sourceDocument.layout.confirmed=confirmed;}
void compareRoles(const StandaloneSessionModel& s,const std::map<ChannelRole,int>& intended)
{for(int n=0;n<s.sourceDocument.layout.count;++n)expect(s.sourceDocument.routing.fileToDevice[n]==intended.at(MonitoringSettings::speakerRole(s.sourceDocument.layout.roles[n])),"each file channel routes by its confirmed speaker role, including deliberate mute");}
struct MonitoringPlaybackCatalog final : OutputDeviceCatalog
{
    OutputDeviceIdentity identity{"Injected","Software monitoring","software-output"};Device* device=nullptr;
    std::vector<OutputDeviceChoice> scan() override {return {{identity,{}}};}
    juce::StringArray outputConfiguration(const OutputDeviceIdentity&) override {return {"1","2","3","4","5","6"};}
    std::unique_ptr<juce::AudioIODevice> create(const OutputDeviceIdentity&) override
    {auto result=std::make_unique<Device>();device=result.get();device->outputCount=6;return result;}
};
void restoredSoftwareOutput(const juce::File& root)
{
    using C=juce::AudioChannelSet;root.createDirectory();const auto input=fixture(root.getChildFile("software.wav"),48000,6,16384);
    auto catalog=std::make_shared<MonitoringPlaybackCatalog>();auto pref=settings(root);pref.lastOutput=catalog->identity;pref.automaticReconnect=true;pref.save();
    MonitoringSettings memory(pref.file());OutputRouting original;original.fileToDevice={0,1,2,-1,4,5,-1,-1,-1,-1,-1,-1};
    expect(memory.rememberRouting(catalog->identity,catalog->outputConfiguration(catalog->identity),input.context.fileLayout,original)&&memory.rememberVolume(catalog->identity,-24),"persist role routing and volume for software callback inspection");
    StandaloneSessionModel s;adoptInput(s,input);s.sourceDocument.layout.roles={C::centre,C::LFE,C::rightSurround,C::left,C::leftSurround,C::right};
    auto audioOwner=std::make_unique<StandaloneAudition>();auto& audio=*audioOwner;StandaloneAudioPanel ui(s,audio,[]{},root.getChildFile("roles"),catalog,settings(root));
    expect(audio.connected()&&!audio.active()&&catalog->device&&!catalog->device->isPlaying(),"restoration opens only the injected device and never starts its callback");
    Capture capture(6,16384);expect(audio.setCaptureObserver(&capture),"capture actual consumed programme before monitor gain");
    button(ui,"Play selected")->onClick();awaitPlaying(audio,*catalog->device);
    juce::AudioBuffer<float> output(6,512);bool mapped=true,seen=false;const float level=juce::Decibels::decibelsToGain(-24.0f);
    const int sourceForOutput[]{3,5,0,-1,4,2};const auto deadline=juce::Time::getMillisecondCounterHiRes()+10000;
    while(audio.active()&&juce::Time::getMillisecondCounterHiRes()<deadline)
    {
        const auto first=capture.next;catalog->device->tick(output);
        if(first>=8192&&capture.next>first)
        {
            seen=true;
            for(int channel=0;channel<6;++channel)for(int n=0;n<int(capture.next-first);++n)
            {const int source=sourceForOutput[channel];const float expected=source<0?0:capture.data.getSample(source,int(first)+n)*level;mapped&=std::abs(output.getSample(channel,n)-expected)<0.000001f;}
        }
        audio.poll();std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    expect(seen&&mapped&&capture.valid&&capture.next==16384&&audio.getState()==AuditionState::finished,"actual software output follows reordered speaker roles, muted LFE and restored -24 dB volume after the existing ramp");
    audio.disconnect();
}
void routingMemory(const juce::File& root)
{
    using C=juce::AudioChannelSet;
    root.createDirectory();const auto stereo=fixture(root.getChildFile("stereo.wav"),48000,2,48013),surround=fixture(root.getChildFile("surround.wav"),48000,6,48013);
    auto catalog=std::make_shared<SettingsCatalog>();catalog->calls->outputs={"Main L","Main R","Centre","LFE","Surround L","Surround R","Phones L","Phones R"};
    const auto originalOutputs=catalog->calls->outputs;const auto studio=catalog->devices[0].identity,other=catalog->devices[1].identity;
    auto local=settings(root);local.lastOutput=studio;local.save();
    StandaloneSessionModel s;adoptInput(s,stereo);s.initializeCurrentStateFromSource();s.setCurrentGainDb(-2);
    auto audio=std::make_unique<StandaloneAudition>();
    auto ui=std::make_unique<StandaloneAudioPanel>(s,*audio,[]{},root.getChildFile("roles"),catalog,settings(root));ui->setSize(700,760);
    expect(!audio->connected()&&catalog->calls->opens==0&&muted(s)&&monitorSlider(*ui).getValue()==-30,"remembered selection without saved monitoring stays muted at -30 dB and does not connect");
    const auto revision=s.revision,document=s.documentRevision;const auto identity=OfflineAnalysisService::identityFor(s.measurementContext(MeasurementDomain::currentState),MeasurementDomain::currentState,s.sourceDocument.media.sha256);
    setRoute(*ui,0,6);setRoute(*ui,1,7);monitorSlider(*ui).setValue(-18.5,juce::sendNotificationSync);
    // JUCE may defer listener callbacks; user changes must survive a routine refresh.
    auto& pendingRoute=routeBox(*ui,0);pendingRoute.setSelectedId(5+2,juce::sendNotificationAsync);ui->refresh();
    expect(pendingRoute.getSelectedId()==7,"routine refresh preserves a pending user route choice");pendingRoute.onChange();setRoute(*ui,0,6);
    monitorSlider(*ui).setValue(-18,juce::sendNotificationAsync);
    expect(MonitoringSettings(settings(root).file()).volume(studio)==-18,"explicit monitor change is captured before deferred JUCE listeners");
    monitorSlider(*ui).setValue(-18.5,juce::sendNotificationAsync);
    MonitoringSettings memory(settings(root).file());const auto stableMemory=memory.file().loadFileAsString();
    expect(memory.volume(studio)==-18.5&&memory.restore(studio,originalOutputs,s.sourceDocument.layout).matched,"explicit route and volume changes persist immediately while disconnected");
    for(int n=0;n<100;++n)ui->refresh();
    expect(s.sourceDocument.routing.fileToDevice[0]==6&&monitorSlider(*ui).getValue()==-18.5&&memory.file().loadFileAsString()==stableMemory,"routine refresh neither overwrites current adjustments nor rewrites remembered settings");
    expect(s.revision==revision&&s.documentRevision==document&&s.currentStateDocument.staticGainDb==-2&&OfflineAnalysisService::identityFor(s.measurementContext(MeasurementDomain::currentState),MeasurementDomain::currentState,s.sourceDocument.media.sha256)==identity,"monitoring edits do not dirty portable programme, processing or measurement identity");
    auto snap=sessionfile::capture(s);JobControl job;const auto project=root.getChildFile("saved.levelscope");expect(sessionfile::save(snap,project,job).succeeded(),"save actual portable session with local monitoring configured");
    expect(!project.loadFileAsString().contains("monitorDb")&&!project.loadFileAsString().contains("Main L"),"portable session excludes monitor volume and device output metadata");
    button(*ui,"Connect")->onClick();expect(audio->connected()&&monitorSlider(*ui).getValue()==-18.5&&s.sourceDocument.routing.fileToDevice[0]==6,"explicit Connect preserves restored routing and volume");stopped(*audio,*catalog);
    button(*ui,"Disconnect")->onClick();const auto opened=catalog->calls->opens;
    dynamic_cast<juce::ComboBox*>(namedChild(*ui,"Output interface"))->onChange();routeBox(*ui,0).onChange();
    expect(!audio->connected()&&catalog->calls->opens==opened&&memory.file().loadFileAsString()==stableMemory,"queued internal control notifications cannot reconnect or replace remembered user settings");
    s.newDocument();ui->sourceChanged();expect(monitorSlider(*ui).getValue()==-18.5&&memory.file().loadFileAsString()==stableMemory,"New and internal disconnect/reset never overwrite user monitor volume");
    adoptInput(s,surround,false);ui->sourceChanged();expect(muted(s)&&!routeBox(*ui,0).isEnabled(),"unconfirmed role suggestions cannot restore or store routing");
    button(*ui,"Confirm roles")->onClick();expect(muted(s)&&monitoringText(*ui).contains("No saved 5.1"),"5.1 without memory stays muted independently of existing stereo mapping");
    const std::map<ChannelRole,int> surroundMap{{C::left,0},{C::right,1},{C::centre,2},{C::LFE,-1},{C::leftSurround,4},{C::rightSurround,5}};
    for(int n=0;n<s.sourceDocument.layout.count;++n)setRoute(*ui,n,surroundMap.at(s.sourceDocument.layout.roles[n]));
    // Reordered 5.1, with the supported side-surround aliases, keeps speaker destinations.
    adoptInput(s,surround);s.sourceDocument.layout.roles={C::centre,C::LFE,C::rightSurroundSide,C::left,C::leftSurroundSide,C::right};ui->sourceChanged();compareRoles(s,surroundMap);
    expect(!audio->connected()&&catalog->calls->opens==opened&&monitorSlider(*ui).getValue()==-18.5,"routing restoration never defeats explicit Disconnect or changes volume for another layout");
    ui->setSize(700,760);auto png=root.getChildFile("reordered-routing.png").createOutputStream();expect(png&&juce::PNGImageFormat().writeImageToStream(ui->createComponentSnapshot(ui->getLocalBounds()),*png),"render reordered role routing and monitor-memory review");png.reset();
    adoptInput(s,stereo);ui->sourceChanged();compareRoles(s,{{C::left,6},{C::right,7}});
    for(int n=0;n<2;++n)setRoute(*ui,n,-1);expect(MonitoringSettings(settings(root).file()).restore(studio,originalOutputs,s.sourceDocument.layout).matched,"all-muted stereo configuration is saved, not treated as absence");
    adoptInput(s,surround);ui->sourceChanged();compareRoles(s,surroundMap);adoptInput(s,stereo);ui->sourceChanged();expect(muted(s),"alternating jobs restores deliberate all-muted stereo without assigning any output");
    // Switching interfaces carries neither routes nor level from the previous interface.
    chooseOutput(*ui,2);expect(audio->connected()&&muted(s)&&monitorSlider(*ui).getValue()==-30,"unsaved interface starts muted with independent -30 dB default");
    expect(!MonitoringSettings(settings(root).file()).volume(other),"initial default on another interface is never persisted as a user volume");
    setRoute(*ui,0,0);setRoute(*ui,1,1);monitorSlider(*ui).setValue(-9,juce::sendNotificationSync);
    chooseOutput(*ui,1);expect(muted(s)&&monitorSlider(*ui).getValue()==-18.5,"return to studio restores all-muted routes and its last user volume");chooseOutput(*ui,2);compareRoles(s,{{C::left,0},{C::right,1}});expect(monitorSlider(*ui).getValue()==-9,"return to other interface restores its own independent volume");
    chooseOutput(*ui,1);button(*ui,"Disconnect")->onClick();
    // Open the actual saved session through its existing verifier/installer.
    const auto read=sessionfile::read(project,job);expect(bool(read.document),"read saved portable session");const auto verified=sessionfile::verify(*read.document,job);expect(bool(verified.document),"verify reopened programme");
    juce::String reason;expect(sessionfile::install(*verified.document,s,reason),"install reopened session");ui->sourceChanged();
    expect(muted(s)&&monitorSlider(*ui).getValue()==-18.5&&s.currentStateDocument.staticGainDb==-2&&!audio->connected(),"reopened session uses latest local all-muted routing/volume and retains portable processing");
    ui.reset();audio.reset();
    // Real application-object restart, with connection preference off and then on.
    for(bool reconnect:{false,true})
    {
        auto preference=settings(root);preference.automaticReconnect=reconnect;preference.save();
        StandaloneSessionModel reopened;expect(sessionfile::install(*verified.document,reopened,reason),"fresh application installs saved programme");
        auto outputOwner=std::make_unique<StandaloneAudition>();auto& output=*outputOwner;StandaloneAudioPanel panel(reopened,output,[]{},root.getChildFile("roles"),catalog,settings(root));
        expect(output.connected()==reconnect&&muted(reopened)&&monitorSlider(panel).getValue()==-18.5,"restart restores monitoring with reconnection preference on or off, without unmuting");stopped(output,*catalog);
        button(panel,"Disconnect")->onClick();const auto before=catalog->calls->opens;
        reopened.newDocument();panel.sourceChanged();adoptInput(reopened,surround);panel.sourceChanged();compareRoles(reopened,surroundMap);
        expect(!output.connected()&&catalog->calls->opens==before,"New/import after explicit Disconnect restores monitoring without automatic connection");
    }
    // Same interface identity but changed ordered output list is never guessed.
    auto pref=settings(root);pref.automaticReconnect=false;pref.save();StandaloneSessionModel changed;adoptInput(changed,surround);
    catalog->calls->outputs.set(2,"Renamed centre output");auto outputOwner=std::make_unique<StandaloneAudition>();auto& output=*outputOwner;StandaloneAudioPanel panel(changed,output,[]{},root.getChildFile("roles"),catalog,settings(root));
    expect(muted(changed)&&monitoringText(panel).contains("configuration changed")&&monitoringText(panel).contains("Saved:")&&monitoringText(panel).contains("Available:"),"changed output names mute saved routing and expose old/new lists for review");
    const auto storedBefore=MonitoringSettings(settings(root).file()).file().loadFileAsString();for(int n=0;n<10;++n)panel.refresh();
    expect(MonitoringSettings(settings(root).file()).file().loadFileAsString()==storedBefore,"missing/changed-output fallback never saves an accidental all-muted replacement");
    catalog->calls->outputs=originalOutputs;button(panel,"Scan devices")->onClick();compareRoles(changed,surroundMap);
    catalog->calls->outputs.remove(7);button(panel,"Scan devices")->onClick();expect(muted(changed),"output-count change invalidates even unused physical outputs");
    button(panel,"Remember routing")->onClick();expect(monitoringText(panel).contains("Remembered"),"explicit review can save a deliberately all-muted replacement configuration");
    panel.sourceChanged();expect(muted(changed)&&monitoringText(panel).contains("Restored"),"explicit all-muted replacement restores on next programme event");
    catalog->devices.erase(catalog->devices.begin());button(panel,"Scan devices")->onClick();expect(muted(changed)&&!output.connected()&&monitoringText(panel).contains("No uniquely identified"),"missing remembered interface never substitutes another output");
    catalog->devices.insert(catalog->devices.begin(),{studio,{}});catalog->calls->outputs=originalOutputs;button(panel,"Scan devices")->onClick();
    expect(muted(changed)&&monitoringText(panel).contains("configuration changed"),"returning different output configuration still requires review");
    // Restore original saved mapping, then observe a live configuration change while open.
    for(int n=0;n<changed.sourceDocument.layout.count;++n)setRoute(panel,n,surroundMap.at(changed.sourceDocument.layout.roles[n]));button(panel,"Connect")->onClick();
    catalog->calls->outputs.set(0,"Changed live output");panel.refresh();expect(muted(changed)&&!output.active()&&monitoringText(panel).contains("configuration changed"),"live output metadata change stops transport and removes stale routes");
    // Independent machine namespace, and finite supported volume range / corrupt data fallback.
    MonitoringSettings anotherMac(settings(root,"Mac-B").file());expect(!anotherMac.volume(studio)&&!anotherMac.restore(studio,originalOutputs,surround.context.fileLayout).matched,"another Mac never inherits routing or monitor volume");
    MonitoringSettings store(settings(root).file());expect(!store.rememberVolume(studio,1)&&!store.rememberVolume(studio,-61)&&!store.rememberVolume(studio,std::numeric_limits<double>::quiet_NaN()),"invalid monitor levels cannot replace a remembered user value");
    expect(store.rememberVolume(studio,-60)&&MonitoringSettings(settings(root).file()).volume(studio)==-60&&store.rememberVolume(studio,0)&&MonitoringSettings(settings(root).file()).volume(studio)==0,"both existing monitor range endpoints round-trip exactly");
    auto parsed=juce::JSON::parse(store.file());auto* entries=parsed.getDynamicObject()->getProperty("interfaces").getArray();expect(entries&&!entries->isEmpty(),"monitor settings have bounded versioned records");
    for(auto& e:*entries)if(e["uid"].toString()==studio.uid)e.getDynamicObject()->setProperty("monitorDb",200);
    store.file().replaceWithText(juce::JSON::toString(parsed));MonitoringSettings corrupt(settings(root).file());expect(!corrupt.volume(studio)&&corrupt.volume(other)==-9,"invalid saved volume falls back independently without erasing another interface");
    StandaloneSessionModel invalidVolume;adoptInput(invalidVolume,stereo);auto invalidOutputOwner=std::make_unique<StandaloneAudition>();auto& invalidOutput=*invalidOutputOwner;
    StandaloneAudioPanel fallback(invalidVolume,invalidOutput,[]{},root.getChildFile("roles"),catalog,settings(root));
    expect(monitorSlider(fallback).getValue()==-30&&!MonitoringSettings(settings(root).file()).volume(studio),"invalid saved volume uses -30 dB without persisting the internal fallback");
    stopped(invalidOutput,*catalog);stopped(output,*catalog);
    restoredSoftwareOutput(root.getChildFile("software-playback"));
}
}

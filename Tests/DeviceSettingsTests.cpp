#include "AuditionTestSupport.h"
#include "Standalone/SessionDocument.h"

namespace
{
struct ConnectionCounts { int creates=0, opens=0, starts=0; bool busy=false; juce::Array<double> rates{44100,48000,96000}; };
struct SettingsDevice final : juce::AudioIODevice
{
    SettingsDevice(const OutputDeviceIdentity& id,std::shared_ptr<ConnectionCounts> counts)
        :AudioIODevice(id.name,id.type),calls(std::move(counts)){}
    std::shared_ptr<ConnectionCounts> calls;bool opened=false;double rate=0;juce::BigInteger outputs;
    juce::StringArray getOutputChannelNames()override{return {"Left output","Right output"};}
    juce::StringArray getInputChannelNames()override{return {};}
    juce::Array<double> getAvailableSampleRates()override{return calls->rates;}
    juce::Array<int> getAvailableBufferSizes()override{return {256};}
    int getDefaultBufferSize()override{return 256;}
    juce::String open(const juce::BigInteger& inputs,const juce::BigInteger& o,double r,int)override
    {expect(inputs.isZero(),"connection opens output only");++calls->opens;opened=!calls->busy;outputs=o;rate=r;return calls->busy?"Device is busy (injected)":"";}
    void close()override{opened=false;}
    bool isOpen()override{return opened;}
    void start(juce::AudioIODeviceCallback*)override{++calls->starts;}
    void stop()override{}
    bool isPlaying()override{return false;}
    juce::String getLastError()override{return {};}
    int getCurrentBufferSizeSamples()override{return 256;}
    double getCurrentSampleRate()override{return rate;}
    int getCurrentBitDepth()override{return 32;}
    juce::BigInteger getActiveOutputChannels()const override{return outputs;}
    juce::BigInteger getActiveInputChannels()const override{return {};}
    int getOutputLatencyInSamples()override{return 0;}
    int getInputLatencyInSamples()override{return 0;}
};
struct SettingsCatalog final : OutputDeviceCatalog
{
    std::shared_ptr<ConnectionCounts> calls=std::make_shared<ConnectionCounts>();
    std::vector<OutputDeviceChoice> devices{{{"Injected","Studio output","uid-studio"},{}},{{"Injected","Other output","uid-other"},{}}};
    bool disappearOnCreate=false;
    std::vector<OutputDeviceChoice> scan()override{return devices;}
    std::unique_ptr<juce::AudioIODevice> create(const OutputDeviceIdentity& id)override
    {++calls->creates;if(disappearOnCreate)return {};return std::make_unique<SettingsDevice>(id,calls);}
};
OutputDeviceSettings settings(const juce::File& root,const juce::String& machine="Mac-A") {return OutputDeviceSettings(root,machine);}
void stopped(const StandaloneAudition& audio,const SettingsCatalog& catalog)
{expect(!audio.active()&&catalog.calls->starts==0,"every open/reconnect leaves transport stopped without starting callbacks");}
void persistence(const juce::File& root)
{
    auto catalog=std::make_shared<SettingsCatalog>();const auto first=catalog->devices[0].identity;
    StandaloneAudition audio;OutputDeviceConnection c(audio,catalog,settings(root));
    expect(!c.selection()&&!c.automaticReconnect()&&!audio.connected(),"fresh Mac has no guessed/default output; auto reconnect off");
    c.select(first);expect(c.awaitingRate()&&catalog->calls->opens==0&&!settings(root).lastOutput,"explicit selection waits for programme rate; not remembered before success");
    c.programmeChanged(48000);expect(audio.connected()&&audio.connectedRate()==48000&&catalog->calls->opens==1,"explicit pending selection opens at validated programme rate");stopped(audio,*catalog);
    expect(settings(root).lastOutput==std::optional<OutputDeviceIdentity>(first),"only successfully opened exact device identity saved");
    expect(!settings(root,"Mac-B").lastOutput&&!settings(root,"Mac-B").automaticReconnect,"other Mac has separate output memory and preference");
    StandaloneAudition restored;OutputDeviceConnection off(restored,catalog,settings(root));off.programmeChanged(48000);
    expect(off.selection()==std::optional<OutputDeviceIdentity>(first)&&!restored.connected()&&catalog->calls->opens==1,"restored selection remains disconnected with preference off");
    catalog->devices[0].identity.name="Renamed studio";off.rescan();
    expect(off.selection()->name=="Renamed studio"&&!restored.connected(),"stable UID restores renamed output visibly without opening");
    off.select(*off.selection());expect(restored.connected()&&catalog->calls->opens==2,"explicit selection opens even with automatic preference off");stopped(restored,*catalog);
    catalog->calls->busy=true;off.select(catalog->devices[1].identity);
    expect(!restored.connected()&&off.status().containsIgnoreCase("busy")&&off.status().contains("offline analysis"),"failed open reports native busy reason without starting output");
    expect(settings(root).lastOutput->uid==first.uid,"failed output does not replace last successfully opened interface");
    catalog->calls->busy=false;off.connect();expect(restored.connected()&&settings(root).lastOutput->uid=="uid-other","successful retry becomes remembered device");stopped(restored,*catalog);
    off.setAutomaticReconnect(true);expect(settings(root).automaticReconnect,"opt-in preference persisted for this Mac");
    expect(!settings(root,"Mac-B").automaticReconnect,"opt-in does not spread to another Mac");
}
void connectionIntent(const juce::File& root)
{
    auto catalog=std::make_shared<SettingsCatalog>();auto saved=settings(root);saved.lastOutput=catalog->devices[0].identity;saved.automaticReconnect=true;saved.save();
    StandaloneAudition audio;OutputDeviceConnection c(audio,catalog,settings(root));
    expect(!audio.connected()&&catalog->calls->opens==0,"opt-in restoration waits for programme rate");
    c.programmeChanged(48000);expect(audio.connected()&&catalog->calls->opens==1,"opt-in connects when programme becomes known");stopped(audio,*catalog);
    c.programmeChanged(48000);expect(catalog->calls->opens==1,"same-rate source lifecycle does not reopen an already connected device");
    c.programmeChanged(96000);expect(audio.connectedRate()==96000&&catalog->calls->opens==2,"opt-in closes and reconnects a changed compatible programme rate");stopped(audio,*catalog);
    c.setAutomaticReconnect(false);expect(audio.connected(),"turning preference off leaves an explicitly usable current connection");
    c.programmeChanged(44100);expect(!audio.connected()&&catalog->calls->opens==2,"rate change with preference off waits for explicit Connect");
    c.connect();expect(audio.connectedRate()==44100&&catalog->calls->opens==3,"explicit Connect uses changed programme rate");
    c.setAutomaticReconnect(true);c.disconnect();const int opened=catalog->calls->opens;
    c.programmeChanged(0);c.programmeChanged(48000);c.rescan();c.setAutomaticReconnect(false);c.setAutomaticReconnect(true);
    expect(!audio.connected()&&catalog->calls->opens==opened&&c.automaticSuppressed(),"Disconnect suppresses New/Open/rate/scan/preference automatic attempts for application session");
    expect(c.status().contains("restarts"),"Disconnect lifetime explained beside controls");
    c.connect();expect(audio.connected()&&c.automaticSuppressed(),"manual Connect remains possible without removing Disconnect suppression");
    c.programmeChanged(96000);expect(!audio.connected()&&catalog->calls->opens==opened+1,"later rate change still respects suppression after manual Connect");
    c.programmeChanged(0);c.select(catalog->devices[0].identity);expect(c.awaitingRate(),"manual selection may await a new source even after Disconnect");
    c.programmeChanged(48000);expect(audio.connected()&&catalog->calls->opens==opened+2,"pending explicit selection is honoured despite automatic suppression");
    StandaloneAudition restarted;OutputDeviceConnection next(restarted,catalog,settings(root));next.programmeChanged(48000);
    expect(restarted.connected()&&!next.automaticSuppressed(),"new application session restores opted-in reconnect after previous explicit Disconnect");stopped(restarted,*catalog);
}
void unavailable(const juce::File& root)
{
    auto catalog=std::make_shared<SettingsCatalog>();auto saved=settings(root);saved.lastOutput=catalog->devices[0].identity;saved.automaticReconnect=true;saved.save();
    const auto known=*saved.lastOutput;catalog->devices[0].identity.uid="replacement-same-name";
    StandaloneAudition audio;OutputDeviceConnection c(audio,catalog,settings(root));c.programmeChanged(48000);
    expect(!audio.connected()&&catalog->calls->opens==0&&c.status().contains("Missing"),"missing remembered UID never falls back to same name or another/default output");
    catalog->devices={{known,{}},{known,{}}};c.rescan();expect(c.status().contains("Ambiguous")&&catalog->calls->opens==0,"duplicate matching identities cannot connect");
    catalog->devices={{known,"Duplicate CoreAudio name"}};c.rescan();expect(c.status().contains("Ambiguous")&&catalog->calls->opens==0,"backend identification ambiguity blocks opening");
    catalog->devices={{known,{}}};catalog->calls->rates={44100};c.rescan();
    expect(!audio.connected()&&catalog->calls->opens==0&&c.status().contains("Incompatible")&&c.status().contains("48000"),"unsupported programme rate reports exact rate without opening/resampling");
    catalog->calls->rates={44100,48000};catalog->calls->busy=true;c.rescan();
    expect(!audio.connected()&&catalog->calls->opens==1&&c.status().containsIgnoreCase("busy")&&c.status().contains("offline analysis"),"busy output leaves offline analysis available and reason visible");
    expect(catalog->calls->opens==1,"failed automatic attempt has no timer retry loop");
    catalog->calls->busy=false;catalog->disappearOnCreate=true;c.connect();
    expect(!audio.connected()&&catalog->calls->opens==1&&c.status().contains("became unavailable"),"device disappearing between discovery and construction fails without substitution");
    catalog->disappearOnCreate=false;c.connect();expect(audio.connected(),"explicit retry recovers once the chosen device returns");stopped(audio,*catalog);
    audio.disconnect();expect(c.status().contains("disconnected")&&!c.status().startsWith("Connected"),"external device closure cannot leave a stale connected status");
    c.programmeChanged(std::numeric_limits<double>::quiet_NaN());expect(!audio.connected(),"unknown/nonfinite rate closes output and never opens");
    auto nameOnly=settings(root,"NameOnly");nameOnly.lastOutput=OutputDeviceIdentity{"Injected","Studio output",{}};nameOnly.save();
    StandaloneAudition named;OutputDeviceConnection byName(named,catalog,settings(root,"NameOnly"));
    expect(byName.selection()->uid==known.uid&&!named.connected(),"unique legacy name may acquire available UID but stays disconnected");
    catalog->devices.push_back({{"Injected","Studio output","different"},{}});
    StandaloneAudition ambiguous;OutputDeviceConnection duplicate(ambiguous,catalog,settings(root,"NameOnly"));duplicate.programmeChanged(48000);duplicate.connect();
    expect(!ambiguous.connected()&&duplicate.status().contains("Ambiguous"),"name fallback requires unique match");
}
void malformed(const juce::File& root)
{
    auto store=settings(root);expect(store.file().getParentDirectory().createDirectory().wasOk(),"settings directory");
    expect(store.file().replaceWithText("{\"schema\":1,\"automaticReconnect\":true,\"machine\":\"another-Mac\"}"),"write foreign settings fixture");
    auto foreign=settings(root);expect(!foreign.lastOutput&&!foreign.automaticReconnect&&foreign.warning.isNotEmpty(),"copied or malformed machine settings fail safely with automatic preference off");
    expect(store.file().replaceWithText("not json"),"write malformed fixture");auto bad=settings(root);
    expect(!bad.automaticReconnect&&!bad.lastOutput&&!bad.warning.isEmpty(),"corrupt preference never opens/defaults a device");
    auto unavailableId=OutputDeviceSettings(root,{});expect(unavailableId.file()==juce::File{}&&!unavailableId.automaticReconnect,"missing machine ID does not use another Mac's settings");
}
juce::TextButton* button(juce::Component& component,const juce::String& name)
{for(auto* child:component.getChildren())if(auto* b=dynamic_cast<juce::TextButton*>(child);b&&b->getButtonText()==name)return b;return nullptr;}
void panel(const juce::File& root)
{
    auto input=fixture(root.getChildFile("panel.wav"),48000,2,4800);StandaloneSessionModel s;
    expect(s.setValidatedSource(input.context.sourceFile,input.context.media),"UI source validated");s.sourceDocument.layout=input.context.fileLayout;
    s.initializeCurrentStateFromSource();s.setCurrentGainDb(-3);const auto revision=s.revision,document=s.documentRevision;const auto sourceRoles=s.sourceDocument.layout.roles;
    auto catalog=std::make_shared<SettingsCatalog>();auto saved=settings(root);saved.lastOutput=catalog->devices[0].identity;saved.automaticReconnect=true;saved.save();
    StandaloneAudition audio;StandaloneAudioPanel ui(s,audio,[]{},root.getChildFile("roles"),catalog,settings(root));ui.setSize(700,660);ui.resized();
    juce::ComboBox* outputs=nullptr;juce::Slider* gain=nullptr;juce::ToggleButton* automatic=nullptr;
    for(auto* child:ui.getChildren())
    {
        if(auto* box=dynamic_cast<juce::ComboBox*>(child);box&&box->getName()=="Output interface")outputs=box;
        if(auto* slider=dynamic_cast<juce::Slider*>(child))gain=slider;
        if(auto* toggle=dynamic_cast<juce::ToggleButton*>(child))automatic=toggle;
    }
    expect(outputs&&outputs->getText().contains("Studio output")&&automatic&&automatic->getToggleState(),"remembered selection and opt-in visibly restored in panel");
    expect(audio.connected()&&gain&&gain->getValue()==-30,"automatic UI open resets safe monitor level");
    expect(std::all_of(s.sourceDocument.routing.fileToDevice.begin(),s.sourceDocument.routing.fileToDevice.end(),[](int route){return route<0;}),"automatic UI open leaves every physical route muted");
    expect(s.revision==revision&&s.documentRevision==document&&s.sourceDocument.layout.roles==sourceRoles&&s.currentStateDocument.staticGainDb==-3,"device preferences preserve roles, processing, document and audio identities");
    gain->setValue(0,juce::sendNotificationSync);button(ui,"Disconnect")->onClick();button(ui,"Connect")->onClick();
    expect(gain->getValue()==-30&&audio.connected(),"explicit reconnect restores safe monitor level");stopped(audio,*catalog);
    s.job.begin(JobStage::analysis,s.revision);ui.refresh();const int opened=catalog->calls->opens;
    expect(!outputs->isEnabled()&&!automatic->isEnabled()&&!button(ui,"Connect")->isEnabled(),"active analysis locks connection edits while remaining offline-capable");
    button(ui,"Connect")->onClick();expect(catalog->calls->opens==opened,"disabled Connect handler cannot bypass active job guard");s.job.finish(s.job.ticket,s.revision,ErrorCode::none,"done");ui.refresh();
    auto stream=root.getChildFile("audio-setup.png").createOutputStream();expect(stream&&juce::PNGImageFormat().writeImageToStream(ui.createComponentSnapshot(ui.getLocalBounds()),*stream),"write device panel review image");
    button(ui,"Disconnect")->onClick();s.newDocument();ui.sourceChanged();expect(!audio.connected()&&catalog->calls->opens==opened,"New project does not defeat explicit Disconnect");
    s.setValidatedSource(input.context.sourceFile,input.context.media);s.sourceDocument.layout=input.context.fileLayout;ui.sourceChanged();
    expect(!audio.connected()&&catalog->calls->opens==opened,"fresh import in same application session still respects Disconnect");
    catalog->calls->busy=true;button(ui,"Connect")->onClick();const auto failedOpens=catalog->calls->opens;
    for(int i=0;i<100;++i)ui.refresh();
    expect(!audio.connected()&&catalog->calls->opens==failedOpens,"repeated UI refresh does not retry a busy output");
}
void projectReopen(const juce::File& root)
{
    root.createDirectory();auto input=fixture(root.getChildFile("project.wav"),48000,2,4800);StandaloneSessionModel s;
    s.setValidatedSource(input.context.sourceFile,input.context.media);s.sourceDocument.layout=input.context.fileLayout;s.initializeCurrentStateFromSource();s.setCurrentGainDb(-2);
    auto document=sessionfile::capture(s);JobControl job;const auto file=root.getChildFile("project.levelscope");
    expect(sessionfile::save(document,file,job).succeeded(),"save actual session before remembered-device reopen");
    const auto read=sessionfile::read(file,job);expect(read.document.has_value(),"read saved session independently of device settings");
    const auto verified=sessionfile::verify(*read.document,job);expect(verified.document.has_value(),"verify saved project media");
    auto catalog=std::make_shared<SettingsCatalog>();auto saved=settings(root,"Mac-A");saved.lastOutput=catalog->devices[0].identity;saved.automaticReconnect=true;saved.save();
    for(const auto& mac:{juce::String("Mac-A"),juce::String("Mac-B")})
    {
        StandaloneSessionModel reopened;juce::String reason;
        expect(sessionfile::install(*verified.document,reopened,reason),"reopen keeps project processing and role assignments");
        StandaloneAudition audio;StandaloneAudioPanel ui(reopened,audio,[]{},root.getChildFile("roles"),catalog,settings(root,mac));
        expect(audio.connected()==(mac=="Mac-A"),"reopened project uses this Mac's opt-in connection preference only");stopped(audio,*catalog);
        expect(reopened.currentStateDocument.staticGainDb==-2&&reopened.sourceDocument.layout.roles==s.sourceDocument.layout.roles
            &&reopened.sourceDocument.layout.confirmed,"reconnection preserves reopened processing and confirmed roles");
        expect(std::all_of(reopened.sourceDocument.routing.fileToDevice.begin(),reopened.sourceDocument.routing.fileToDevice.end(),[](int route){return route<0;}),"reopened project retains muted physical routes with preference on or off");
    }
}
void discovery()
{
    // Read-only metadata lookup on the actual OS; no create/open/start call.
    const auto outputs=OutputDeviceCatalog::system()->scan();int core=0,stable=0,ambiguous=0;
    for(const auto& entry:outputs)
    {
        expect(entry.identity.type.isNotEmpty()&&entry.identity.name.isNotEmpty(),"native catalog gives visible backend and output name");
        if(entry.identity.type=="CoreAudio")++core;
        if(entry.identity.uid.isNotEmpty()&&entry.problem.isEmpty())++stable;
        if(entry.problem.isNotEmpty())++ambiguous;
    }
    expect(core==0||stable>0,"actual CoreAudio enumeration resolves a stable output UID without opening hardware");
    auto data=juce::var(new juce::DynamicObject());data.getDynamicObject()->setProperty("read_only_outputs",int(outputs.size()));
    data.getDynamicObject()->setProperty("stable_identities",stable);data.getDynamicObject()->setProperty("ambiguous_entries",ambiguous);cases.add(data);
}
void preparationGuard(const juce::File& root)
{
    root.createDirectory();auto request=fixture(root.getChildFile("prepare.wav"),48000,2,48000);
    auto catalog=std::make_shared<SettingsCatalog>();StandaloneAudition audio;OutputDeviceConnection c(audio,catalog,settings(root));
    c.select(catalog->devices[0].identity);c.programmeChanged(48000);c.disconnect();
    audio.sourceChanged(request.context.media);audio.synchronize(request,routes(2));expect(audio.prepare(request).isEmpty()&&audio.active(),"injected preparation active for connection guard");
    const auto selected=c.selection();const int opened=catalog->calls->opens;
    c.connect();c.select(catalog->devices[1].identity);c.rescan();
    expect(c.selection()==selected&&catalog->calls->opens==opened&&audio.active(),"queued connection/selection/scan cannot interrupt active preparation");
    audio.stop();waitPreparation(audio);c.programmeChanged(48000);
    expect(!audio.connected()&&catalog->calls->opens==opened,"rejected preparation-time Connect creates no later pending open");
    stopped(audio,*catalog);
}
}
int main(int argc,char** argv)
{
    juce::ScopedJuceInitialiser_GUI gui;
    if(argc!=2)return 2;const juce::File root(argv[1]);if(root.exists())return 2;root.createDirectory();bool passed=false;juce::String error;
    try{persistence(root.getChildFile("memory"));connectionIntent(root.getChildFile("intent"));unavailable(root.getChildFile("unavailable"));malformed(root.getChildFile("malformed"));auto ui=root.getChildFile("ui");ui.createDirectory();panel(ui);projectReopen(root.getChildFile("project"));preparationGuard(root.getChildFile("preparation"));discovery();passed=true;}
    catch(const std::exception& e){error=e.what();}
    auto report=juce::var(new juce::DynamicObject());auto* r=report.getDynamicObject();r->setProperty("passed",passed);r->setProperty("error",error);r->setProperty("checks",checks);
    r->setProperty("build",LEVELSCOPE_MEASUREMENT_BUILD);r->setProperty("os",juce::SystemStats::getOperatingSystemName());r->setProperty("physical_devices_opened",false);
    r->setProperty("catalog_observations",cases);
    root.getChildFile("result.json").replaceWithText(juce::JSON::toString(report,true));std::cout<<(passed?"PASS":"FAIL")<<": "<<checks.size()<<" checks "<<error<<std::endl;return passed?0:1;
}

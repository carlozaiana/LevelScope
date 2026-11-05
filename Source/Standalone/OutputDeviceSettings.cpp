#include "OutputDeviceSettings.h"
#include <juce_cryptography/juce_cryptography.h>
#if JUCE_MAC
#include <CoreAudio/CoreAudio.h>
#endif

namespace levelscope::standalone
{
namespace
{
#if JUCE_MAC
juce::String stringProperty(AudioObjectID object, AudioObjectPropertySelector selector)
{
    AudioObjectPropertyAddress address{selector, kAudioObjectPropertyScopeGlobal, 0};
    CFStringRef value = nullptr; UInt32 bytes = sizeof value;
    if (AudioObjectGetPropertyData(object, &address, 0, nullptr, &bytes, &value) != noErr || !value) return {};
    const auto text = juce::String::fromCFString(value); CFRelease(value); return text;
}
std::vector<OutputDeviceIdentity> coreAudioIdentities()
{
    AudioObjectPropertyAddress address{kAudioHardwarePropertyDevices, kAudioObjectPropertyScopeGlobal, 0};
    UInt32 bytes = 0;
    if (AudioObjectGetPropertyDataSize(kAudioObjectSystemObject, &address, 0, nullptr, &bytes) != noErr
        || bytes > 4096 * sizeof(AudioDeviceID)) return {};
    std::vector<AudioDeviceID> devices(bytes / sizeof(AudioDeviceID));
    if (AudioObjectGetPropertyData(kAudioObjectSystemObject, &address, 0, nullptr, &bytes, devices.data()) != noErr) return {};
    std::vector<OutputDeviceIdentity> result;
    for (auto device : devices)
    {
        address = {kAudioDevicePropertyStreamConfiguration, kAudioDevicePropertyScopeOutput, 0}; bytes = 0;
        if (AudioObjectGetPropertyDataSize(device, &address, 0, nullptr, &bytes) != noErr
            || bytes < sizeof(AudioBufferList) || bytes > 1024 * 1024) continue;
        juce::HeapBlock<char> storage(bytes);
        if (AudioObjectGetPropertyData(device, &address, 0, nullptr, &bytes, storage.getData()) != noErr) continue;
        const auto* buffers = reinterpret_cast<const AudioBufferList*>(storage.getData());
        if (buffers->mNumberBuffers > (bytes - offsetof(AudioBufferList, mBuffers)) / sizeof(AudioBuffer)) continue;
        bool outputs = false;
        for (UInt32 i = 0; i < buffers->mNumberBuffers; ++i) outputs |= buffers->mBuffers[i].mNumberChannels > 0;
        if (outputs) result.push_back({"CoreAudio", stringProperty(device, kAudioObjectPropertyName),
                                                   stringProperty(device, kAudioDevicePropertyDeviceUID)});
    }
    return result;
}
#endif

class SystemOutputCatalog final : public OutputDeviceCatalog
{
public:
    std::vector<OutputDeviceChoice> scan() override
    {
        if (!manager) manager = std::make_unique<juce::AudioDeviceManager>();
        std::vector<OutputDeviceChoice> result;
#if JUCE_MAC
        const auto identities = coreAudioIdentities();
#endif
        for (auto* type : manager->getAvailableDeviceTypes())
        {
            type->scanForDevices();
            for (const auto& name : type->getDeviceNames(false))
            {
                OutputDeviceChoice choice{{type->getTypeName(), name, {}}, {}};
#if JUCE_MAC
                if (choice.identity.type == "CoreAudio")
                {
                    int matches = 0;
                    for (const auto& identity : identities) if (identity.name == name)
                    { choice.identity.uid = identity.uid; ++matches; }
                    // JUCE disambiguates duplicate names with numeric suffixes, but
                    // exposes no public HAL ID. Never infer a UID from those suffixes.
                    if (matches != 1)
                        choice.problem = "Ambiguous or changed CoreAudio name. Give each output a unique name in Audio MIDI Setup, then scan again.";
                }
#endif
                result.push_back(std::move(choice));
            }
        }
        return result;
    }
    std::unique_ptr<juce::AudioIODevice> create(const OutputDeviceIdentity& wanted) override
    {
        const auto current = scan();
        int matches = 0;
        for (const auto& entry : current) if (entry.identity == wanted && entry.problem.isEmpty()) ++matches;
        if (matches != 1) return {}; // changed between selection and construction
        for (auto* type : manager->getAvailableDeviceTypes())
            if (type->getTypeName() == wanted.type)
                return std::unique_ptr<juce::AudioIODevice>(type->createDevice(wanted.name, {}));
        return {};
    }
private:
    std::unique_ptr<juce::AudioDeviceManager> manager;
};
bool validText(const juce::var& text) { return text.isString() && text.toString().length() <= 1024; }
}

std::shared_ptr<OutputDeviceCatalog> OutputDeviceCatalog::system() { return std::make_shared<SystemOutputCatalog>(); }

OutputDeviceSettings::OutputDeviceSettings(juce::File directory, juce::String machineId)
{
    if (machineId.isEmpty())
    { warning = "This Mac's identity is unavailable; output settings will not be saved or restored."; return; }
    machineKey = juce::SHA256(machineId.toRawUTF8(), size_t(machineId.getNumBytesAsUTF8())).toHexString();
    if (directory == juce::File{}) directory = juce::File::getSpecialLocation(juce::File::userApplicationDataDirectory)
        .getChildFile("LevelScope/Output Settings");
    path = directory.getChildFile(machineKey + ".json");
    if (!path.existsAsFile()) return;
    const auto data = path.getSize() <= 16384 ? juce::JSON::parse(path) : juce::var();
    const auto last = data["lastOutput"];
    if (!data.isObject() || int(data["schema"]) != 1 || data["machine"].toString() != machineKey
        || !data["automaticReconnect"].isBool()
        || (!last.isVoid() && (!last.isObject() || !validText(last["type"]) || !validText(last["name"])
            || !validText(last["uid"]) || last["type"].toString().isEmpty() || last["name"].toString().isEmpty())))
    { warning = "Saved output settings are invalid; automatic reconnection is off. Select an output explicitly."; return; }
    if (last.isObject()) lastOutput = OutputDeviceIdentity{last["type"], last["name"], last["uid"]};
    automaticReconnect = bool(data["automaticReconnect"]);
}
void OutputDeviceSettings::save()
{
    if (path == juce::File{}) return;
    auto data = juce::var(new juce::DynamicObject()); auto* object = data.getDynamicObject();
    object->setProperty("schema", 1); object->setProperty("machine", machineKey);
    object->setProperty("automaticReconnect", automaticReconnect);
    if (lastOutput)
    {
        auto last = juce::var(new juce::DynamicObject()); auto* entry = last.getDynamicObject();
        entry->setProperty("type", lastOutput->type); entry->setProperty("name", lastOutput->name); entry->setProperty("uid", lastOutput->uid);
        object->setProperty("lastOutput", last);
    }
    juce::TemporaryFile temporary(path);
    if (path.getParentDirectory().createDirectory().failed()
        || !temporary.getFile().replaceWithText(juce::JSON::toString(data)) || !temporary.overwriteTargetFileWithTemporary())
        warning = "Could not save output settings on this Mac. This connection can still be used.";
    else warning.clear();
}

OutputDeviceConnection::OutputDeviceConnection(StandaloneAudition& engine, std::shared_ptr<OutputDeviceCatalog> provider, OutputDeviceSettings stored)
    : audio(engine), catalog(provider ? std::move(provider) : OutputDeviceCatalog::system()), settings(std::move(stored)), selected(settings.lastOutput)
{ available = catalog->scan(); describeDisconnected(); }

juce::String OutputDeviceConnection::status() const
{
    if(message.startsWith("Connected ") && selected)
    {
        if(!audio.connected())return "Selected " + selected->name + "; disconnected. Connect or scan to retry; offline analysis is available.";
        return "Connected " + selected->name + " at " + juce::String(audio.connectedRate(),0) + " Hz. "
            + (audio.active() ? "Audition is separate from measurement." : "Playback stopped; review output routing.");
    }
    return message;
}

std::optional<OutputDeviceIdentity> OutputDeviceConnection::resolve()
{
    if (!selected) { message = "No output selected. Offline analysis is available."; return {}; }
    std::vector<const OutputDeviceChoice*> matches;
    for (const auto& entry : available)
        if (entry.identity.type == selected->type && (selected->uid.isNotEmpty()
            ? entry.identity.uid == selected->uid : entry.identity.name == selected->name)) matches.push_back(&entry);
    if (matches.empty())
    { message = "Missing output: " + selected->name + ". Reconnect it and scan, or choose another output. Offline analysis is available."; return {}; }
    if (matches.size() != 1 || matches.front()->problem.isNotEmpty())
    {
        message = "Ambiguous output: " + selected->name + ". " + (matches.size() == 1 ? matches.front()->problem
            : "More than one device matches. Choose a uniquely identified output.") + " Offline analysis is available.";
        return {};
    }
    return matches.front()->identity;
}
void OutputDeviceConnection::describeDisconnected()
{
    const auto found = resolve();
    if (!found)
    {
        if(suppressed)message += " Automatic reconnection paused until LevelScope restarts.";
        return;
    }
    selected = *found;
    message = "Selected " + selected->name + (selected->uid.isNotEmpty() ? " (device identity). " : " (matched by name; no stable ID). ");
    if (suppressed) message += "Automatic reconnection paused until LevelScope restarts. Connect or select explicitly.";
    else if (explicitPending && programmeRate <= 0) message += "Will connect when the programme sample rate is known; playback stays stopped.";
    else if (settings.automaticReconnect && programmeRate <= 0) message += "Automatic reconnection awaits the programme sample rate.";
    else message += "Disconnected. Connect or select explicitly. Offline analysis is available.";
}
void OutputDeviceConnection::attempt()
{
    if (audio.active()) { message = "Stop playback or finish preparation before connecting an output."; return; }
    if (programmeRate <= 0) { describeDisconnected(); return; }
    explicitPending = false;
    if(audio.connected()) { audio.disconnect(); ++revision; }
    available = catalog->scan(); const auto found = resolve();
    if (!found) return;
    auto device = catalog->create(*found);
    if (!device) { message = "Output became unavailable: " + found->name + ". Scan or select again; offline analysis is available."; return; }
    selected = *found;
    audio.selectOutputDevice(std::move(device)); ++revision;
    const auto error = audio.connect(programmeRate);
    if (error.isNotEmpty())
    { message = "Cannot connect " + selected->name + ": " + error + ". Choose another output or retry; offline analysis is available."; return; }
    settings.lastOutput = selected; settings.save();
    message = "Connected " + selected->name + " at " + juce::String(programmeRate, 0) + " Hz. Playback stopped; review output routing.";
}
void OutputDeviceConnection::select(const OutputDeviceIdentity& choice)
{
    if (audio.active()) { message = "Stop playback or finish preparation before selecting an output."; return; }
    audio.selectOutputDevice({}); ++revision; selected = choice; explicitPending = true; attempt();
}
void OutputDeviceConnection::connect()
{
    if(audio.active()) { message = "Stop playback or finish preparation before connecting an output."; return; }
    explicitPending = true; attempt();
}
void OutputDeviceConnection::disconnect()
{
    suppressed = true; explicitPending = false; audio.disconnect(); ++revision; describeDisconnected();
}
void OutputDeviceConnection::setAutomaticReconnect(bool enabled)
{
    settings.automaticReconnect = enabled; settings.save();
    if (!audio.connected() && enabled && !suppressed) attempt();
    else if (!audio.connected()) describeDisconnected();
}
void OutputDeviceConnection::programmeChanged(double rate)
{
    programmeRate = std::isfinite(rate) && rate > 0 ? rate : 0;
    if (audio.connected() && (programmeRate <= 0 || std::abs(audio.connectedRate() - programmeRate) > .01))
    { audio.disconnect(); ++revision; }
    if (!audio.connected())
    {
        if (explicitPending || (settings.automaticReconnect && !suppressed)) attempt();
        else describeDisconnected();
    }
}
void OutputDeviceConnection::rescan()
{
    if(audio.active()) { message = "Stop playback or finish preparation before scanning outputs."; return; }
    available = catalog->scan();
    const auto found = resolve();
    if (!found) { if (audio.connected()) { audio.disconnect(); ++revision; } return; }
    selected = *found;
    if (!audio.connected())
    {
        if (explicitPending || (settings.automaticReconnect && !suppressed)) attempt();
        else describeDisconnected();
    }
}
}

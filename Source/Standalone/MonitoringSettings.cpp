#include "MonitoringSettings.h"

namespace levelscope::standalone
{
namespace
{
bool text(const juce::var& v) { return v.isString() && v.toString().isNotEmpty() && v.toString().length()<=1024; }
bool outputsValid(const juce::StringArray& names)
{
    return !names.isEmpty() && names.size()<=256 && std::all_of(names.begin(),names.end(),[](const auto& n){return n.isNotEmpty()&&n.length()<=256;});
}
bool volumeValid(double db) { return std::isfinite(db)&&db>=-60&&db<=0; }
juce::StringArray readOutputs(const juce::var& v)
{
    juce::StringArray result;if(auto* a=v.getArray();a&&a->size()<=256)for(const auto& n:*a){if(!text(n))return {};result.add(n.toString());}return result;
}
juce::var writeOutputs(const juce::StringArray& a) { juce::Array<juce::var> values;for(const auto& n:a)values.add(n);return values; }
}
juce::String MonitoringSettings::identityKey(const OutputDeviceIdentity& id)
{
    // Length-delimited strings avoid ambiguous keys and never match a UID to a name.
    const auto identity=id.uid.isNotEmpty()?"uid:"+id.uid:"name:"+id.name;
    return juce::String(id.type.length())+":"+id.type+identity;
}
ChannelRole MonitoringSettings::speakerRole(ChannelRole role)
{
    if(role==juce::AudioChannelSet::leftSurroundSide)return juce::AudioChannelSet::leftSurround;
    if(role==juce::AudioChannelSet::rightSurroundSide)return juce::AudioChannelSet::rightSurround;
    return role;
}
juce::String MonitoringSettings::layoutKey(const ChannelLayout& layout)
{
    if(!layout.isValid() || (layout.count!=2&&layout.count!=6))return {};
    using C=juce::AudioChannelSet;
    std::vector<ChannelRole> expected{C::left,C::right};if(layout.count==6)expected={C::left,C::right,C::centre,C::LFE,C::leftSurround,C::rightSurround};
    for(int n=0;n<layout.count;++n)
    {auto found=std::find(expected.begin(),expected.end(),speakerRole(layout.roles[n]));if(found==expected.end())return {};expected.erase(found);}
    return layout.count==2?"stereo":"5.1";
}
MonitoringSettings::MonitoringSettings(const juce::File& connectionFile)
{
    if(connectionFile==juce::File{}){notice="This Mac's identity is unavailable; monitoring settings cannot be remembered.";return;}
    machine=connectionFile.getFileNameWithoutExtension();path=connectionFile.getSiblingFile(machine+".monitoring.json");
    if(!path.existsAsFile())return;
    const auto data=path.getSize()<=1024*1024?juce::JSON::parse(path):juce::var();const auto* list=data["interfaces"].getArray();
    if(!data.isObject()||int(data["schema"])!=1||data["machine"].toString()!=machine||!list||list->size()>32)
    {notice="Invalid local monitoring settings; routes remain muted and monitor uses its initial default.";return;}
    for(const auto& v:*list)
    {
        if(!text(v["type"])||!text(v["name"])||!v["uid"].isString()||v["uid"].toString().length()>1024)
        {notice="Invalid interface in monitoring settings; review routing and volume.";continue;}
        Interface item;item.identity={v["type"],v["name"],v["uid"]};const auto key=identityKey(item.identity);
        if(interfaces.count(key)){interfaces.clear();notice="Duplicate monitoring interface settings; review routing and volume.";return;}
        const auto gain=v["monitorDb"];
        if(!gain.isVoid())
        {if((gain.isDouble()||gain.isInt())&&volumeValid(double(gain)))item.db=double(gain);else notice="Invalid saved monitor volume; using -30 dB until explicitly adjusted.";}
        for(const auto& name:{juce::String("stereo"),juce::String("5.1")})
        {
            const auto route=v[juce::Identifier(name)];if(route.isVoid())continue;MonitoringRoutes saved;saved.outputs=readOutputs(route["outputs"]);
            const auto* roles=route["roles"].getArray();ChannelLayout layout;OutputRouting mapping;
            bool valid=outputsValid(saved.outputs)&&roles&&roles->size()==(name=="stereo"?2:6);
            if(valid)
            {
                layout.count=roles->size();
                for(int n=0;n<roles->size();++n)
                {
                    const auto& r=(*roles)[n];if(!r["speaker"].isInt()||!r["output"].isInt()){valid=false;break;}
                    layout.roles[n]=static_cast<ChannelRole>(int(r["speaker"]));mapping.fileToDevice[n]=int(r["output"]);
                    saved.roles[speakerRole(layout.roles[n])]=mapping.fileToDevice[n];
                }
                valid=valid&&layoutKey(layout)==name&&mapping.validate(layout.count,std::min(maxChannels,saved.outputs.size()),true).isEmpty()&&saved.roles.size()==size_t(layout.count);
            }
            if(valid)item.layouts[name]=std::move(saved);else notice="Invalid saved routing; review the affected layout and assign outputs explicitly.";
        }
        interfaces.emplace(key,std::move(item));
    }
}
std::optional<double> MonitoringSettings::volume(const OutputDeviceIdentity& id) const
{const auto it=interfaces.find(identityKey(id));return it==interfaces.end()?std::optional<double>{}:it->second.db;}
MonitoringRestore MonitoringSettings::restore(const OutputDeviceIdentity& id,const juce::StringArray& outputs,const ChannelLayout& layout) const
{
    MonitoringRestore result;const auto name=layoutKey(layout);
    if(!layout.confirmed){result.explanation="Confirm speaker roles before restoring routing.";return result;}
    if(name.isEmpty()){result.explanation="Routing memory supports stereo and 5.1; configure this layout manually.";return result;}
    if(!outputsValid(outputs)){result.explanation="Output list unavailable; routes muted. Scan devices or Connect to verify this interface. Offline analysis is available.";return result;}
    const auto device=interfaces.find(identityKey(id));
    if(device==interfaces.end()||!device->second.layouts.count(name))
    {result.explanation="No saved "+name+" routing for this interface; initial routes muted. Edit outputs or Remember routing to save deliberate mutes.";return result;}
    const auto& saved=device->second.layouts.at(name);
    if(saved.outputs!=outputs)
    {
        result.explanation="Saved output configuration changed; routes muted. Review outputs below; edits or Remember routing save the new configuration.";
        result.explanation+="\nSaved: "+saved.outputs.joinIntoString(" / ")+"\nAvailable: "+outputs.joinIntoString(" / ");return result;
    }
    for(int n=0;n<layout.count;++n)result.routing.fileToDevice[n]=saved.roles.at(speakerRole(layout.roles[n]));
    result.matched=true;result.explanation="Restored "+name+" routing for "+id.name+", including all deliberately muted roles.";return result;
}
MonitoringSettings::Interface* MonitoringSettings::entry(const OutputDeviceIdentity& id)
{
    if(id.type.isEmpty()||id.name.isEmpty()||id.type.length()>1024||id.name.length()>1024||id.uid.length()>1024)
    {notice="Interface identity is unavailable; monitoring settings cannot be remembered.";return nullptr;}
    const auto key=identityKey(id);
    if(!interfaces.count(key)&&interfaces.size()>=32){notice="Monitoring memory is full (32 interfaces); current settings remain usable but cannot be saved.";return nullptr;}
    auto& item=interfaces[key];item.identity=id;return &item;
}
bool MonitoringSettings::rememberVolume(const OutputDeviceIdentity& id,double db)
{
    if(!volumeValid(db))return false;
    if(auto* item=entry(id)){item->db=db;return save();}return false;
}
bool MonitoringSettings::rememberRouting(const OutputDeviceIdentity& id,const juce::StringArray& outputs,const ChannelLayout& layout,const OutputRouting& routing)
{
    const auto name=layoutKey(layout);
    if(!layout.confirmed||name.isEmpty()||!outputsValid(outputs)||routing.validate(layout.count,std::min(maxChannels,outputs.size()),true).isNotEmpty())return false;
    if(auto* item=entry(id))
    {
        MonitoringRoutes saved;saved.outputs=outputs;
        for(int n=0;n<layout.count;++n)saved.roles[speakerRole(layout.roles[n])]=routing.fileToDevice[n];
        item->layouts[name]=std::move(saved);return save();
    }
    return false;
}
bool MonitoringSettings::save()
{
    if(path==juce::File{})return false;
    auto data=juce::var(new juce::DynamicObject());data.getDynamicObject()->setProperty("schema",1);data.getDynamicObject()->setProperty("machine",machine);
    juce::Array<juce::var> entries;
    for(const auto& pair:interfaces)
    {
        const auto& item=pair.second;auto v=juce::var(new juce::DynamicObject());auto* object=v.getDynamicObject();
        object->setProperty("type",item.identity.type);object->setProperty("name",item.identity.name);object->setProperty("uid",item.identity.uid);
        if(item.db)object->setProperty("monitorDb",*item.db);
        for(const auto& layout:item.layouts)
        {
            auto r=juce::var(new juce::DynamicObject());r.getDynamicObject()->setProperty("outputs",writeOutputs(layout.second.outputs));juce::Array<juce::var> roles;
            for(const auto& entry:layout.second.roles){auto role=juce::var(new juce::DynamicObject());role.getDynamicObject()->setProperty("speaker",int(entry.first));role.getDynamicObject()->setProperty("output",entry.second);roles.add(role);}
            r.getDynamicObject()->setProperty("roles",roles);object->setProperty(layout.first,r);
        }
        entries.add(v);
    }
    data.getDynamicObject()->setProperty("interfaces",entries);const auto json=juce::JSON::toString(data);juce::TemporaryFile temp(path);
    if(json.getNumBytesAsUTF8()>1024*1024 || path.getParentDirectory().createDirectory().failed() || !temp.getFile().replaceWithText(json) || !temp.overwriteTargetFileWithTemporary())
    {notice="Could not save local monitoring settings; current adjustments remain usable. Review storage and try again.";return false;}
    notice.clear();return true;
}
}

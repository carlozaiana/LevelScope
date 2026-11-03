#include "RememberedChannelRoles.h"

namespace levelscope::standalone
{
namespace
{
ChannelRole assemblyRole(ChannelRole role)
{
    using C=juce::AudioChannelSet;
    if(role==C::leftSurroundSide)return C::leftSurround;
    if(role==C::rightSurroundSide)return C::rightSurround;
    return role;
}
bool eligible(const ChannelLayout& order)
{
    auto expected=assemblyRoles(order.count);
    if(expected.empty() || !order.isValid())return false;
    for(int i=0;i<order.count;++i)
    {
        const auto it=std::find(expected.begin(),expected.end(),assemblyRole(order.roles[i]));
        if(it==expected.end())return false;
        expected.erase(it);
    }
    return true;
}
juce::String names(const ChannelLayout& layout)
{
    juce::StringArray values;
    for(int i=0;i<layout.count;++i)
        values.add(juce::AudioChannelSet::getAbbreviatedChannelTypeName(layout.roles[i]));
    return values.joinIntoString(" / ");
}
ChannelRoleSuggestion suggested(int count,const std::optional<ChannelLayout>& remembered)
{
    ChannelRoleSuggestion result;
    result.order.count=count;
    if(remembered)
    {
        result.order=*remembered;
        result.explanation="Last explicitly confirmed "+juce::String(count==2?"stereo":"5.1")+" order on this Mac";
    }
    else
    {
        const auto roles=assemblyRoles(count);
        std::copy(roles.begin(),roles.end(),result.order.roles.begin());
        result.explanation="Suggested standard order for "+juce::String(count)+" channels";
    }
    return result;
}
void describeLabels(ChannelRoleSuggestion& result,const std::optional<ChannelLayout>& memory,const juce::String& origin)
{
    result.explanation=origin+" — "+names(result.order);
    if(memory && !std::equal(memory->roles.begin(),memory->roles.begin()+memory->count,result.order.roles.begin()))
        result.explanation+=". Differs from remembered order — "+names(*memory);
}
}
RememberedChannelRoles::RememberedChannelRoles(juce::File folder)
    :directory(folder==juce::File{} ? juce::File::getSpecialLocation(juce::File::userApplicationDataDirectory)
        .getChildFile("LevelScope/Confirmed Channel Orders") : folder) {}
juce::File RememberedChannelRoles::file(int count) const
{
    return directory.getChildFile(count==2?"stereo.json":"5.1.json");
}
std::optional<ChannelLayout> RememberedChannelRoles::read(int count) const
{
    if(count!=2 && count!=6)return {};
    const auto path=file(count);
    if(!path.existsAsFile() || path.getSize()>1024)return {};
    const auto data=juce::JSON::parse(path);
    auto* channels=data.getArray();
    if(!channels || channels->size()!=count)return {};
    ChannelLayout result;
    result.count=count;
    for(int i=0;i<count;++i)
    {
        if(!(*channels)[i].isInt())return {};
        result.roles[i]=static_cast<ChannelRole>(int((*channels)[i]));
    }
    if(!eligible(result))return {};
    return result;
}
juce::Result RememberedChannelRoles::remember(const ChannelLayout& layout) const
{
    if(!layout.confirmed || !eligible(layout))return juce::Result::ok();
    if(auto status=directory.createDirectory();status.failed())return status;
    juce::Array<juce::var> channels;
    for(int i=0;i<layout.count;++i)channels.add(int(layout.roles[i]));
    juce::TemporaryFile temp(file(layout.count));
    if(!temp.getFile().replaceWithText(juce::JSON::toString(juce::var(channels))) || !temp.overwriteTargetFileWithTemporary())
        return juce::Result::fail("Could not save this channel order on this Mac.");
    return juce::Result::ok();
}
ChannelRoleSuggestion RememberedChannelRoles::suggest(const MediaInfo& media) const
{
    if(media.channels!=2 && media.channels!=6)return {};
    const auto memory=read(media.channels);
    auto result=suggested(media.channels,memory);
    if(media.hasFileRoleLabels)
    {
        result.order=media.suggestedLayout;
        describeLabels(result,memory,media.layoutProvenance);
    }
    else result.explanation+=" — "+names(result.order);
    result.order.confirmed=false;
    return result;
}
ChannelRoleSuggestion RememberedChannelRoles::suggest(const ProgrammeAssembly& group) const
{
    const int count=int(group.members.size());
    if(count!=2 && count!=6)return {};
    const auto memory=read(count);
    auto result=suggested(count,memory);
    const auto standard=assemblyRoles(count);
    auto available=standard;
    bool any=false;
    for(int i=0;i<count;++i)
    {
        const auto& media=group.members[size_t(i)].media;
        if(!media.hasFileRoleLabels)continue;
        any=true;
        result.order.roles[i]=assemblyRole(media.suggestedLayout.roles[0]);
        const auto it=std::find(available.begin(),available.end(),result.order.roles[i]);
        if(it!=available.end())available.erase(it);
    }
    if(any)
    {
        // Keep explicitly labelled slots fixed. Fill the others from the
        // proposed order, using each still-available programme role once.
        for(int i=0;i<count;++i)
        {
            if(group.members[size_t(i)].media.hasFileRoleLabels)continue;
            auto it=std::find(available.begin(),available.end(),assemblyRole(result.order.roles[i]));
            if(it==available.end())it=available.begin();
            if(it!=available.end()){result.order.roles[i]=*it;available.erase(it);}
        }
        describeLabels(result,memory,"File labels first, remaining slots from "+result.explanation);
        if(!eligible(result.order))result.explanation+=". Conflicting member labels require Edit order.";
    }
    else
    {
        for(int i=0;i<count;++i)result.order.roles[i]=assemblyRole(result.order.roles[i]);
        result.explanation+=" — "+names(result.order);
    }
    return result;
}
}

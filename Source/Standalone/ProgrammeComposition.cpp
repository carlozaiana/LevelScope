#include "ProgrammeComposition.h"
#include <juce_cryptography/juce_cryptography.h>
#include <limits>
#include <sstream>
#include <iomanip>
#include <locale>

namespace levelscope::standalone
{
namespace
{
juce::String hash (const std::string& s) { return juce::SHA256 (s.data(), s.size()).toHexString(); }
ChannelRole conventional (ChannelRole r)
{
    using C = juce::AudioChannelSet;
    return r == C::leftSurroundSide ? C::leftSurround : r == C::rightSurroundSide ? C::rightSurround : r;
}
bool hasRoles (const ChannelLayout& layout, const std::vector<ChannelRole>& roles)
{
    if (!layout.confirmed || !layout.isValid() || layout.count != static_cast<int>(roles.size())) return false;
    for (auto role : roles)
    {
        int count = 0;
        for (int i=0; i<layout.count; ++i) if (conventional(layout.roles[i]) == role) ++count;
        if (count != 1) return false;
    }
    return true;
}
}
std::vector<ProgrammeMember> ProgrammeComponent::files() const
{
    if (assembly) return assembly->members;
    return {{ "file", file, media, layout.roles[0], 0 }};
}
juce::String ProgrammeComponent::validate() const
{
    if (id.isEmpty() || id.containsChar('/') || (purpose != "music-effects" && purpose != "dialogue")) return "Explicit supported component identity and purpose required";
    if (offsetFrames != 0) return "Only shared sample-zero alignment is supported";
    if (media.frames <= 0 || media.bytes <= 0 || media.sha256.length() != 64 || !media.sha256.containsOnly("0123456789abcdef")
        || media.channels != layout.count || !layout.confirmed || !layout.isValid() || media.measurementRestriction.isNotEmpty())
        return "Verified component media and explicitly confirmed roles required";
    if (media.sampleRate != 44100 && media.sampleRate != 48000 && media.sampleRate != 96000 && media.sampleRate != 192000)
        return "Composition supports 44.1/48/96/192 kHz";
    if (assembly)
    {
        auto error = assembly->validate(); if (error.isNotEmpty()) return error;
        const auto m = assembly->summary();
        if (m.sha256 != media.sha256 || m.frames != media.frames || m.sampleRate != media.sampleRate
            || m.bytes != media.bytes || m.channels != media.channels || m.bitsPerSample != media.bitsPerSample || m.format != media.format)
            return "Component assembly does not match its verified summary";
        for (int i=0;i<layout.count;++i) if (layout.roles[i] != assembly->members[static_cast<size_t>(i)].role)
            return "Component assembly role identity differs";
    }
    else
    {
        const bool flac = media.format.containsIgnoreCase("flac");
        if ((!flac && !media.format.containsIgnoreCase("wav") && !media.format.containsIgnoreCase("aiff"))
            || (media.bitsPerSample != 16 && media.bitsPerSample != 24 && !(media.bitsPerSample == 32 && !flac)))
            return "Component requires a qualified PCM/float WAV, AIFF or FLAC source";
    }
    return {};
}
juce::String ProgrammeComponent::audioIdentity() const
{
    std::ostringstream s; s.imbue(std::locale::classic());
    s << "programme-component-v1\n" << purpose << '\n' << media.sha256 << '\n' << media.format << '\n'
      << std::hexfloat << media.sampleRate << ':' << media.frames << ':' << media.channels << ':' << media.bitsPerSample
      << ':' << media.bytes << ':' << offsetFrames << '\n';
    for (int i=0;i<std::min(layout.count,maxChannels);++i) s << static_cast<int>(layout.roles[i]) << ',';
    if (assembly) s << '\n' << assembly->definitionIdentity();
    return hash(s.str());
}
bool ProgrammeComponent::matchesFileStamps() const
{ return assembly ? assembly->matchesFileStamps() : media.matchesFileStamp(file); }
ChannelLayout ProgrammeComposition::layout() const
{
    ChannelLayout result; const auto roles = assemblyRoles(preset == CompositionPreset::stereoStereo ? 2 : 6);
    result.count = static_cast<int>(roles.size());
    std::copy(roles.begin(),roles.end(),result.roles.begin()); result.manuallyInterpreted=true;
    result.confirmed = confirmation.isNotEmpty() && confirmation == definitionIdentity(); return result;
}
std::vector<int> ProgrammeComposition::map (bool d) const
{
    const auto roles = assemblyRoles(preset == CompositionPreset::stereoStereo ? 2 : 6);
    const auto& input = d ? dialogue.layout : me.layout; std::vector<int> result;
    for (auto role:roles)
    {
        int index=-1;
        for(int i=0;i<std::min(input.count,maxChannels);++i) if(conventional(input.roles[i])==role) index=i;
        result.push_back(index);
    }
    return result;
}
juce::String ProgrammeComposition::definitionIdentity() const
{
    return hash(std::string(engineVersion)+"\n"+std::to_string(static_cast<int>(preset))+"\n"
        +me.audioIdentity().toStdString()+"\n"+dialogue.audioIdentity().toStdString()+"\nunity-components;sample-zero;ordered-role-sum\n");
}
juce::String ProgrammeComposition::validate (bool confirmed) const
{
    if (preset != CompositionPreset::stereoStereo && preset != CompositionPreset::surroundMonoCentre
        && preset != CompositionPreset::surroundSurround) return "Unsupported composition preset";
    for(const auto* c:{&me,&dialogue}) {const auto error=c->validate();if(error.isNotEmpty())return error;}
    if(me.id==dialogue.id || me.purpose!="music-effects" || dialogue.purpose!="dialogue") return "Distinct M&E and dialogue component slots required";
    if(me.media.sampleRate!=dialogue.media.sampleRate || me.media.frames!=dialogue.media.frames) return "Components must have identical sample rates and exact frame counts";
    if(me.media.bytes>std::numeric_limits<juce::int64>::max()-dialogue.media.bytes) return "Composition byte count overflow";
    const auto roles=assemblyRoles(preset==CompositionPreset::stereoStereo?2:6);
    const auto dr=preset==CompositionPreset::surroundMonoCentre?std::vector<ChannelRole>{juce::AudioChannelSet::centre}:roles;
    if(!hasRoles(me.layout,roles)||!hasRoles(dialogue.layout,dr)) return "Explicit roles do not match the selected composition mapping";
    if(confirmed&&confirmation!=definitionIdentity()) return "Confirm these exact complementary, correctly balanced stems are intended to be combined and aligned at sample zero";
    return {};
}
bool ProgrammeComposition::confirmComplementaryBalancedAlignment()
{ if(validate(false).isNotEmpty())return false;confirmation=definitionIdentity();return true; }
std::vector<ProgrammeMember> ProgrammeComposition::files() const
{
    std::vector<ProgrammeMember> result;
    for(const auto* c:{&me,&dialogue})for(auto m:c->files()){m.id=c->id+"/"+m.id;result.push_back(std::move(m));}
    return result;
}
bool ProgrammeComposition::matchesFileStamps() const
{
    if(!me.matchesFileStamps()||!dialogue.matchesFileStamps())return false;
    const auto members=files();
    for(size_t i=0;i<members.size();++i)for(size_t j=0;j<i;++j)if(samePhysicalFile(members[i].file,members[j].file))return false;
    return true;
}
MediaInfo ProgrammeComposition::summary() const
{
    MediaInfo m; m.sampleRate=me.media.sampleRate;m.frames=me.media.frames;m.channels=preset==CompositionPreset::stereoStereo?2:6;
    if(me.media.bytes>0&&dialogue.media.bytes>0&&me.media.bytes<=std::numeric_limits<juce::int64>::max()-dialogue.media.bytes)
        m.bytes=me.media.bytes+dialogue.media.bytes;
    m.bitsPerSample=32;m.format="Composed programme";m.sha256=definitionIdentity();m.suggestedLayout=layout();
    m.layoutProvenance="Explicit aligned M&E plus dialogue; unity components; source definition digest, not a file hash";return m;
}
}

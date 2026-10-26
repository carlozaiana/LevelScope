#include "ProgrammeSource.h"
#include <juce_cryptography/juce_cryptography.h>
#include <sys/stat.h>
#include <limits>
#include <sstream>
#include <iomanip>
#include <locale>
#include <set>
#include <cmath>

namespace levelscope::standalone
{
std::vector<ChannelRole> assemblyRoles (int n)
{
    using C = juce::AudioChannelSet;
    if (n == 2) return { C::left, C::right };
    if (n == 6) return { C::left, C::right, C::centre, C::LFE, C::leftSurround, C::rightSurround };
    return {};
}
bool samePhysicalFile (const juce::File& a, const juce::File& b)
{
    if (a == b) return true;
    struct stat x{}, y{};
    return ::stat (a.getFullPathName().toRawUTF8(), &x) == 0
        && ::stat (b.getFullPathName().toRawUTF8(), &y) == 0 && x.st_dev == y.st_dev && x.st_ino == y.st_ino;
}
juce::String ProgrammeAssembly::definitionIdentity() const
{
    std::ostringstream key; key.imbue (std::locale::classic());
    key << "mono-channel-assembly-v" << version << '\n' << members.size() << '\n';
    for (const auto& m : members)
        key << m.media.sha256 << '\n' << m.media.format << '\n' << std::hexfloat << m.media.sampleRate
            << ':' << m.media.frames << ':' << m.media.channels << ':' << m.media.bitsPerSample
            << ':' << m.media.bytes << ':' << static_cast<int>(m.role) << ":0:" << m.offsetFrames << '\n';
    const auto bytes = key.str();
    return juce::SHA256 (bytes.data(), bytes.size()).toHexString();
}
juce::String ProgrammeAssembly::validate (bool confirmation) const
{
    const auto roles = assemblyRoles (static_cast<int>(members.size()));
    if (roles.empty()) return "A programme requires two mono files (stereo) or six mono files (5.1)";
    const auto& first = members.front().media;
    if (first.sampleRate != 44100 && first.sampleRate != 48000 && first.sampleRate != 96000 && first.sampleRate != 192000)
        return "Channel assembly supports 44.1/48/96/192 kHz";
    juce::int64 bytes = 0; std::set<juce::String> ids;
    for (size_t i = 0; i < members.size(); ++i)
    {
        const auto& m = members[i]; const auto& info = m.media;
        const auto prefix = "Channel " + juce::String(static_cast<int>(i)+1) + " (" + m.file.getFileName() + "): ";
        if (m.id.isEmpty() || !ids.insert(m.id).second) return prefix + "unique file slot required";
        if (m.role != roles[i]) return prefix + "explicit canonical stereo/5.1 role required";
        if (m.offsetFrames != 0) return prefix + "only sample-zero alignment is supported";
        if (info.channels != 1 || info.frames <= 0 || info.bytes <= 0 || info.sha256.length() != 64
            || !info.sha256.containsOnly("0123456789abcdef") || info.measurementRestriction.isNotEmpty())
            return prefix + "verified supported mono media required";
        const bool flac=info.format.containsIgnoreCase("flac");
        if ((!flac&&!info.format.containsIgnoreCase("wav")&&!info.format.containsIgnoreCase("aiff"))
            || (info.bitsPerSample!=16&&info.bitsPerSample!=24&&!(info.bitsPerSample==32&&!flac)))
            return prefix + "assembly requires qualified WAV/AIFF PCM16/24/32 or float32 WAV, or FLAC PCM16/24";
        if (info.sampleRate != first.sampleRate || info.frames != first.frames)
            return prefix + "sample rate and exact frame count must match every member";
        if (bytes > std::numeric_limits<juce::int64>::max() - info.bytes) return "Programme byte count overflow";
        bytes += info.bytes;
        // Identity validation is filesystem-independent. Alias checks occur when
        // importing/verifying/copying, not when evaluating an immutable snapshot.
    }
    if (confirmation && alignmentConfirmation != definitionIdentity())
        return "Confirm that these exact files and channel assignments form the same aligned mix; sample zero is shared";
    return {};
}
bool ProgrammeAssembly::matchesFileStamps() const
{
    for (size_t i = 0; i < members.size(); ++i)
    {
        if (!members[i].media.matchesFileStamp(members[i].file)) return false;
        for (size_t j = 0; j < i; ++j) if (samePhysicalFile(members[i].file, members[j].file)) return false;
    }
    return !members.empty();
}
ChannelLayout ProgrammeAssembly::layout() const
{
    ChannelLayout result; result.count = static_cast<int>(members.size());
    for (int i = 0; i < std::min(result.count, maxChannels); ++i) result.roles[i] = members[static_cast<size_t>(i)].role;
    result.manuallyInterpreted = true;
    result.confirmed = validate().isEmpty();
    return result;
}
MediaInfo ProgrammeAssembly::summary() const
{
    MediaInfo result;
    if (!members.empty()) { result.sampleRate = members.front().media.sampleRate; result.frames = members.front().media.frames; }
    result.channels = static_cast<int>(members.size());
    for (const auto& m : members) if (m.media.bytes > 0 && result.bytes <= std::numeric_limits<juce::int64>::max()-m.media.bytes) result.bytes += m.media.bytes;
    result.bitsPerSample = 32; // reader presents float32; each member retains its original encoding
    result.format = "Mono-file programme";
    result.sha256 = definitionIdentity(); result.suggestedLayout = layout();
    result.layoutProvenance = "Explicit channel assembly; sample zero in every member; metadata offsets not applied";
    return result;
}
bool ProgrammeAssembly::confirmAlignment()
{
    if (validate(false).isNotEmpty()) return false;
    alignmentConfirmation = definitionIdentity(); return true;
}
bool ProgrammeAssembly::orderByRoles()
{
    const auto roles=assemblyRoles(static_cast<int>(members.size()));
    if(roles.empty())return false;
    std::vector<ProgrammeMember> ordered;
    for(auto role:roles)
    {
        auto first=std::find_if(members.begin(),members.end(),[&](const auto& m){return m.role==role;});
        if(first==members.end()||std::count_if(members.begin(),members.end(),[&](const auto& m){return m.role==role;})!=1)return false;
        ordered.push_back(*first);
    }
    members=std::move(ordered);return true;
}
ProgrammeInspection inspectProgrammeFiles (const std::vector<juce::File>& files, JobControl& control)
{
    ProgrammeInspection result;
    if(assemblyRoles(static_cast<int>(files.size())).empty())return{{},ErrorCode::invalidRequest,"Select exactly two or six mono files"};
    for(size_t i=0;i<files.size();++i)
    {
        if(control.isCancelled())return{{},ErrorCode::cancelled,"Programme import cancelled; previous session retained"};
        for(size_t j=0;j<i;++j)if(samePhysicalFile(files[i],files[j]))return{{},ErrorCode::invalidRequest,files[i].getFileName()+": same underlying file selected more than once"};
        JobControl child;child.parent=&control;child.parentProgress=&control.progress;child.progressBase=static_cast<float>(i)/files.size();child.progressScale=1.0f/files.size();
        auto valid=validateMedia(files[i],child);
        if(!valid.succeeded())return{{},valid.error,files[i].getFileName()+": "+valid.message};
        if(valid.media.channels!=1)return{{},ErrorCode::unsupportedChannels,files[i].getFileName()+": mono file required"};
        if(!result.group.members.empty()&&(valid.media.sampleRate!=result.group.members[0].media.sampleRate||valid.media.frames!=result.group.members[0].media.frames))
            return{{},ErrorCode::invalidMetadata,files[i].getFileName()+": sample rate and exact frame count must match"};
        result.group.members.push_back({"audio-"+juce::String(static_cast<int>(i)+1),files[i],valid.media,juce::AudioChannelSet::unknown,0});
    }
    // Validate the format contract independently of the as-yet-unassigned roles.
    auto probe=result.group;const auto roles=assemblyRoles(static_cast<int>(files.size()));
    for(size_t i=0;i<roles.size();++i)probe.members[i].role=roles[i];
    const auto error=probe.validate(false);
    if(error.isNotEmpty())return{{},ErrorCode::unsupportedMedia,error};
    if(!probe.matchesFileStamps())return{{},ErrorCode::sourceChanged,"A member changed during group inspection"};
    control.setProgress(1);result.message="All members verified. Assign roles and explicitly confirm alignment.";return result;
}
std::unique_ptr<ProgrammeReader> ProgrammeReader::open (const std::vector<juce::File>& files, const MediaInfo& expected)
{
    if (files.empty() || (files.size() != 1 && static_cast<int>(files.size()) != expected.channels)) return {};
    auto result = std::make_unique<ProgrammeReader>();
    juce::AudioFormatManager formats; registerSourceFormats(formats);
    for (const auto& file : files)
    {
        std::unique_ptr<juce::AudioFormatReader> r (formats.createReaderFor(file.createInputStream()));
        if (!r || r->sampleRate != expected.sampleRate || r->lengthInSamples != expected.frames
            || static_cast<int>(r->numChannels) != (files.size() == 1 ? expected.channels : 1)) return {};
        result->hasFlac |= r->getFormatName().containsIgnoreCase("flac");
        result->readers.push_back(std::move(r));
    }
    result->sampleRate = expected.sampleRate; result->lengthInSamples = expected.frames;
    result->numChannels = static_cast<unsigned>(expected.channels);
    return result;
}
bool ProgrammeReader::read (float* const* output, int channels, juce::int64 start, int count)
{
    if (channels != static_cast<int>(numChannels) || start < 0 || count < 0 || start > lengthInSamples
        || count > lengthInSamples - start) return false;
    if (meReader)
    {
        // Worker-only, bounded scratch; partitioning cannot change sum precision.
        for(int offset=0;offset<count;)
        {
            const int n=std::min(4096,count-offset);
            if(!meReader->read(meBuffer.getArrayOfWritePointers(),meBuffer.getNumChannels(),start+offset,n)
                ||!dialogueReader->read(dialogueBuffer.getArrayOfWritePointers(),dialogueBuffer.getNumChannels(),start+offset,n))return false;
            for(int c=0;c<channels;++c)for(int i=0;i<n;++i)
            {
                const float a=meBuffer.getSample(meMapping[static_cast<size_t>(c)],i);
                const int d=dialogueMapping[static_cast<size_t>(c)];
                const float b=d<0?0.0f:dialogueBuffer.getSample(d,i);
                const float x=d<0?a:static_cast<float>(static_cast<double>(a)+static_cast<double>(b));
                if(!std::isfinite(a)||!std::isfinite(b)||!std::isfinite(x))return false;
                output[c][offset+i]=x;
            }
            offset+=n;
        }
        return true;
    }
    if (readers.size() == 1) return readers.front()->read(output, channels, start, count);
    for (int channel = 0; channel < channels; ++channel)
        if (!readers[static_cast<size_t>(channel)]->read(output + channel, 1, start, count)) return false;
    return true;
}
std::unique_ptr<ProgrammeReader> ProgrammeReader::compose (std::unique_ptr<ProgrammeReader> me,
    std::unique_ptr<ProgrammeReader> dialogue,const std::vector<int>& meMap,const std::vector<int>& dialogueMap)
{
    if(!me||!dialogue||me->sampleRate!=dialogue->sampleRate||me->lengthInSamples!=dialogue->lengthInSamples
        ||meMap.empty()||meMap.size()>maxChannels||meMap.size()!=dialogueMap.size())return {};
    for(size_t c=0;c<meMap.size();++c)
        if(meMap[c]<0||meMap[c]>=static_cast<int>(me->numChannels)||dialogueMap[c]<-1||dialogueMap[c]>=static_cast<int>(dialogue->numChannels))return {};
    auto r=std::make_unique<ProgrammeReader>();r->sampleRate=me->sampleRate;r->lengthInSamples=me->lengthInSamples;
    r->numChannels=static_cast<unsigned>(meMap.size());r->hasFlac=me->hasFlac||dialogue->hasFlac;
    r->meMapping=meMap;r->dialogueMapping=dialogueMap;
    r->meBuffer.setSize(static_cast<int>(me->numChannels),4096);r->dialogueBuffer.setSize(static_cast<int>(dialogue->numChannels),4096);
    r->meReader=std::move(me);r->dialogueReader=std::move(dialogue);return r;
}
}

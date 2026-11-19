#pragma once
#include "SectionGainEditing.h"

namespace levelscope::standalone
{
// Session-local disk journal: no step-count limit and no retained audio, waveform,
// measurement snapshots or growing in-memory index. Only one bounded edit is read.
struct ManualEditState
{
    proposal::SectionMap map;
    std::optional<processing::SectionGainPlan> gains;
    double globalDb=0;
    juce::MemoryBlock bytes() const
    {
        juce::MemoryOutputStream out;
        out.writeString(map.inputIdentity);out.writeInt(map.sampleRate);out.writeInt64(map.frames);out.writeInt64(map.nextId);
        out.writeInt(int(map.sections.size()));
        for(const auto& s:map.sections){out.writeInt64(s.id);out.writeInt64(s.begin);out.writeInt64(s.end);out.writeString(s.name);}
        out.writeDouble(globalDb);out.writeInt(gains?int(gains->intervals.size()):0);
        if(gains)for(const auto& g:gains->intervals){out.writeInt64(g.begin);out.writeInt64(g.end);out.writeDouble(g.db);}
        return out.getMemoryBlock();
    }
    static std::optional<ManualEditState> read(const juce::MemoryBlock& data)
    {
        juce::MemoryInputStream in(data,false);ManualEditState s;
        s.map.inputIdentity=in.readString();s.map.sampleRate=in.readInt();s.map.frames=in.readInt64();s.map.nextId=uint64_t(in.readInt64());
        auto n=in.readInt();if(n<1||n>10000)return {};
        for(int i=0;i<n;++i)s.map.sections.push_back({uint64_t(in.readInt64()),in.readInt64(),in.readInt64(),in.readString()});
        s.globalDb=in.readDouble();n=in.readInt();if(n<0||n>10000||!std::isfinite(s.globalDb))return {};
        if(n){processing::SectionGainPlan p;p.inputIdentity=s.map.inputIdentity;p.frames=s.map.frames;p.sampleRate=s.map.sampleRate;
            for(int i=0;i<n;++i)p.intervals.push_back({in.readInt64(),in.readInt64(),in.readDouble()});
            if(p.validationError().isNotEmpty())return {};s.gains=std::move(p);}
        if(!s.map.valid()||in.getPosition()!=int64_t(data.getSize()))return {};return s;
    }
    bool operator==(const ManualEditState& b) const {return bytes()==b.bytes();}
};
class ManualEditHistory
{
public:
    explicit ManualEditHistory(juce::String key):runtime(std::move(key)),
        directory(juce::File::getSpecialLocation(juce::File::tempDirectory).getChildFile("LevelScope-edits-"+juce::Uuid().toString())) {}
    ~ManualEditHistory(){directory.deleteRecursively();}
    const juce::String runtime;
    struct Edit {ManualEditState before,after;juce::String label;};
    bool push(const ManualEditState& before,const ManualEditState& after,const juce::String& label)
    {
        if(before==after)return true;
        // A proposal/import outside this journal starts a new manual-edit chain.
        if(cursor){auto e=peek(false);if(!e||!(e->after==before)){cursor=0;count=0;}}
        if(directory.createDirectory().failed())return false;
        juce::FileOutputStream data(directory.getChildFile("edits")),index(directory.getChildFile("index"));
        if(!data.openedOk()||!index.openedOk())return false;
        const auto offset=data.getPosition();auto a=before.bytes(),b=after.bytes();
        data.writeString(label);data.writeInt(int(a.getSize()));data.write(a.getData(),a.getSize());data.writeInt(int(b.getSize()));data.write(b.getData(),b.getSize());data.flush();
        if(data.getStatus().failed()||!index.setPosition(int64_t(cursor)*8))return false;
        index.writeInt64(offset);index.flush();if(index.getStatus().failed())return false;
        count=++cursor;return true;
    }
    std::optional<Edit> peek(bool redo) const
    {
        if(redo?cursor>=count:cursor==0)return {};
        juce::FileInputStream index(directory.getChildFile("index")),data(directory.getChildFile("edits"));
        if(!index.openedOk()||!data.openedOk()||!index.setPosition(int64_t(redo?cursor:cursor-1)*8))return {};
        const auto offset=index.readInt64();if(offset<0||!data.setPosition(offset))return {};
        Edit e;e.label=data.readString();
        const auto read=[&]()->std::optional<ManualEditState>{auto n=data.readInt();if(n<=0||n>8*1024*1024)return {};juce::MemoryBlock b;if(data.readIntoMemoryBlock(b,size_t(n))!=size_t(n))return {};return ManualEditState::read(b);};
        auto a=read(),b=read();if(!a||!b)return {};e.before=std::move(*a);e.after=std::move(*b);return e;
    }
    void stepped(bool redo){if(redo)++cursor;else --cursor;}
    uint64_t undoCount() const {return cursor;}
private:
    juce::File directory;
    uint64_t cursor=0,count=0;
};
}

#pragma once
#include <juce_gui_basics/juce_gui_basics.h>
#include <array>

namespace levelscope::standalone
{
// Local display sizes only. Never serialized into a programme or processing identity.
struct TimelineHeights
{
    int graph=180,gain=48;
    std::array<int,5> context{{18,18,18,18,18}};
    std::array<int,3> lra{{58,58,58}};
    void clamp()
    {
        graph=juce::jlimit(80,1600,graph);gain=juce::jlimit(32,800,gain);
        for(auto& h:context)h=juce::jlimit(16,240,h);
        for(auto& h:lra)h=juce::jlimit(44,500,h);
    }
    bool operator==(const TimelineHeights& o)const{return graph==o.graph&&gain==o.gain&&context==o.context&&lra==o.lra;}
    bool operator!=(const TimelineHeights& o)const{return !(*this==o);}
};
class TimelineViewPreferences
{
public:
    explicit TimelineViewPreferences(juce::File file={}):path(std::move(file))
    {
        if(!path.existsAsFile()||path.getSize()>8192)return;
        const auto data=juce::JSON::parse(path);if(!data.isObject()||int(data["schema"])!=1)return;
        read(data["main"],main);read(data["sections"],sections);
    }
    static juce::File defaultFile(){return juce::File::getSpecialLocation(juce::File::userApplicationDataDirectory).getChildFile("LevelScope/View Settings/timeline.json");}
    TimelineHeights main,sections;
    void save()
    {
        if(path==juce::File{})return;
        auto data=juce::var(new juce::DynamicObject());auto* o=data.getDynamicObject();
        o->setProperty("schema",1);o->setProperty("main",write(main));o->setProperty("sections",write(sections));
        juce::TemporaryFile temporary(path);
        if(path.getParentDirectory().createDirectory().wasOk()&&temporary.getFile().replaceWithText(juce::JSON::toString(data)))temporary.overwriteTargetFileWithTemporary();
    }
private:
    juce::File path;
    static void read(const juce::var& v,TimelineHeights& h)
    {
        if(!v.isObject()||!v["graph"].isInt()||!v["gain"].isInt())return;
        const auto* c=v["context"].getArray();const auto* l=v["lra"].getArray();if(!c||c->size()!=5||!l||l->size()!=3)return;
        for(const auto& n:*c)if(!n.isInt())return;for(const auto& n:*l)if(!n.isInt())return;
        h.graph=int(v["graph"]);h.gain=int(v["gain"]);for(int i=0;i<5;++i)h.context[size_t(i)]=int((*c)[i]);for(int i=0;i<3;++i)h.lra[size_t(i)]=int((*l)[i]);h.clamp();
    }
    static juce::var write(TimelineHeights h)
    {
        h.clamp();auto v=juce::var(new juce::DynamicObject());auto* o=v.getDynamicObject();o->setProperty("graph",h.graph);o->setProperty("gain",h.gain);
        juce::Array<juce::var> c,l;for(auto n:h.context)c.add(n);for(auto n:h.lra)l.add(n);o->setProperty("context",c);o->setProperty("lra",l);return v;
    }
};
}

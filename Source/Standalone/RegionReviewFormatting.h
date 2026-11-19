#pragma once
#include "Core/Proposal/RegionPlanning.h"
#include <iomanip>
#include <limits>
#include <sstream>

namespace levelscope::standalone
{
inline juce::String regionExact(double value)
{std::ostringstream s;s.imbue(std::locale::classic());s<<std::setprecision(std::numeric_limits<double>::max_digits10)<<value;return s.str();}
// Documentary presentation only. No ranking, geometry or acceptance changes.
inline juce::String regionActivityDescription(const proposal::RegionPlan& p,int rate,int64_t frames,bool full=true)
{
    const auto& a=p.activity;juce::String text="Independent activity accounting (not an audio-quality score)\n";
    text+="New splits per original root:";size_t total=0;for(size_t n=0;n<a.splits.size();++n){text+=" root "+juce::String(n+1)+" = "+juce::String(a.splits[n])+";";total+=a.splits[n];}
    size_t decisions=0;for(const auto& part:p.assignments)decisions+=part.reductions.size();
    text+="\nTotal splits "+juce::String(total)+"; decision parts including unchanged gaps "+juce::String(decisions)+"; connected changed regions "+juce::String(a.support.size())+".\n";
    text+="Newly affected duration including old/new ramps: "+juce::String(a.affectedFrames)+" / "+juce::String(frames)+" programme frames";
    if(rate>0)text+=" = "+regionExact(double(a.affectedFrames)/rate)+" s";if(frames>0)text+=" = "+regionExact(100.*double(a.affectedFrames)/frames)+"%";text+=".\n";
    text+="Nominal engine intervals "+juce::String(a.intervalCount)+" / 10000; actual transition count "+juce::String(a.transitionCount)+"; longest changed run "+juce::String(a.longestRun)+" frames; margin activity "+juce::String(a.marginFrames)+" frames.\n";
    text+="Maximum added attenuation "+regionExact(a.maximumAdditional)+" dB; maximum cumulative attenuation "+regionExact(a.maximumCumulative)+" dB; absolute local dB*frames "+regionExact(a.absoluteDbFrames)+".\n";
    auto unity=p.envelope;unity.intervals={{0,unity.frames,0}};
    const auto cumulative=proposal::changedEnvelopeSupport(unity,p.envelope);int64_t cumulativeFrames=0;for(const auto& span:cumulative)cumulativeFrames+=span.second-span.first;
    text+="Cumulative non-unity local-envelope support: "+juce::String(cumulativeFrames)+" frames in "+juce::String(cumulative.size())+" components; distinct from newly affected duration.\n";
    text+="Global gain "+regionExact(p.globalDb)+" dB; delta "+regionExact(a.globalDelta)+" dB - separate whole-programme change, excluded from the local duration percentage.\n";
    if(!full)return text;
    text+="Complete changed support (half-open source frames, old/new ramps included):\n";
    for(auto s:a.support)text+="["+juce::String(s.begin)+", "+juce::String(s.end)+")\n";
    text+="Complete nominal plan and fixed transitions (source frames; these are not changed-support edges):\n";
    for(size_t n=0;n<p.envelope.intervals.size();++n)
    {
        const auto& i=p.envelope.intervals[n];text+="["+juce::String(i.begin)+", "+juce::String(i.end)+") plateau "+regionExact(i.db)+" dB";
        if(n){const auto radius=p.envelope.radius(n);text+="; incoming ramp ["+juce::String(i.begin-radius)+", "+juce::String(i.begin+radius)+"), nominal cut "+juce::String(i.begin)+", radius "+juce::String(radius);}
        text+="\n";
    }
    return text;
}
}

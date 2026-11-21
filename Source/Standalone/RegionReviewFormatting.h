#pragma once
#include "Core/Proposal/RegionPlanning.h"
#include "StandaloneRegionProposal.h"
#include <iomanip>
#include <limits>
#include <sstream>

namespace levelscope::standalone
{
inline juce::String regionExact(double value)
{std::ostringstream s;s.imbue(std::locale::classic());s<<std::setprecision(std::numeric_limits<double>::max_digits10)<<value;return s.str();}
inline juce::String eventAttemptCoverage(const RegionAttempt& a)
{return a.family=="measured global-only"?"Measured global-only; no local event treatment":juce::String(proposal::eventCoverageName(a.coverage.category));}
inline juce::String regionBoundaryFamily(proposal::BoundaryFamily family)
{
    return family==proposal::BoundaryFamily::innerAfter?"Inner beginning / later background return":
        family==proposal::BoundaryFamily::innerDecay?"Inner beginning / earlier in-decay return":"Original or measured-need boundaries";
}
inline juce::String eventCoverageSummary(const proposal::EventCoverage& coverage)
{
    const auto full=coverage.text();
    return full.substring(0,1200)+"\nSummary only; the complete attempt report preserves every event and reason.";
}
inline juce::String regionAttemptAuthority(const SectionProposalResult& result,size_t n)
{
    const auto r=result.regionReview();if(!r||n>=r->attempts.size())return "No measured attempt";
    const auto& a=r->attempts[n];if(!r->eventMode)return n==r->retainedAttempt&&a.evaluation.verdict==proposal::Verdict::pass?"Retained measured recommendation":"Diagnostic only";
    const auto label=eventAttemptCoverage(a);
    if(a.evaluation.verdict!=proposal::Verdict::pass)return "Measured failure / unavailable - "+label;
    if(r->selectedAttempt==n&&result.explicitlySelectedAlternative())return "Selected passing alternative - "+label;
    if(r->recommendedAttempt==n)return "Recommended measured pass - "+label;
    return "Review only - measured pass - "+label;
}
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

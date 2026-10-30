#include "SectionTrimMethod.h"
#include <cmath>
#include <stdexcept>

namespace levelscope::proposal
{
namespace
{
double originalReductionTarget(double original,double reduction)
{
    double target=original-reduction;
    // A representable plateau must itself respect the additional cap. Move
    // inward by one ULP if subtraction rounded the actual difference outward.
    if(original-target>reduction)target=std::nextafter(target,original);
    return target;
}
}

std::vector<std::pair<int64_t,int64_t>> changedEnvelopeSupport(const processing::SectionGainPlan& before,
                             const processing::SectionGainPlan& after)
{
    if(before.validationError().isNotEmpty()||after.validationError().isNotEmpty()
        ||before.inputIdentity!=after.inputIdentity||before.frames!=after.frames||before.sampleRate!=after.sampleRate)
        throw std::invalid_argument("Envelope support requires matching valid programme");
    const auto a=before.canonical(),b=after.canonical();
    if(a.audioKey()==b.audioKey())return {};
    // Partition into canonical plateau/ramp segments. Different ramp definitions
    // count their complete support, including equal endpoints: a conservative,
    // versioned support measure, not a count of numerically unequal samples.
    struct Piece {int64_t begin,end,rampBegin=0,rampEnd=0;double left=0,right=0;};
    const auto pieces=[](const processing::SectionGainPlan& p)
    {
        std::vector<Piece> v;int64_t begin=0;
        for(size_t n=1;n<p.intervals.size();++n)
        {
            const auto radius=p.radius(n),start=p.intervals[n].begin-radius,end=p.intervals[n].begin+radius;
            if(begin<start)v.push_back({begin,start,0,0,p.intervals[n-1].db,p.intervals[n-1].db});
            v.push_back({start,end,start,end,p.intervals[n-1].db,p.intervals[n].db});begin=end;
        }
        if(begin<p.frames)v.push_back({begin,p.frames,0,0,p.intervals.back().db,p.intervals.back().db});
        return v;
    };
    const auto x=pieces(a),y=pieces(b);size_t i=0,j=0;int64_t begin=0;
    std::vector<std::pair<int64_t,int64_t>> spans;
    while(i<x.size()&&j<y.size())
    {
        const auto& u=x[i];const auto& v=y[j];const auto end=std::min(u.end,v.end);
        const bool equal=u.rampBegin==v.rampBegin&&u.rampEnd==v.rampEnd&&u.left==v.left&&u.right==v.right;
        if(!equal)
        {
            if(!spans.empty()&&spans.back().second==begin)spans.back().second=end;
            else spans.emplace_back(begin,end);
        }
        begin=end;if(end==u.end)++i;if(end==v.end)++j;
    }
    return spans;
}
int64_t changedEnvelopeExtent(const processing::SectionGainPlan& before,const processing::SectionGainPlan& after)
{
    int64_t total=0;for(auto span:changedEnvelopeSupport(before,after))total+=span.second-span.first;return total;
}
SectionTrimSeed selectSectionTrims(const TrimRequest& request,const SectionMap& map,
    const std::optional<processing::SectionGainPlan>& existing,const SectionPeakMap& peaks,const ProposalPolicy& policy)
{
    SectionTrimSeed out;out.policy=policy;
    if(policy.validationError().isNotEmpty()){out.reason=policy.validationError();return out;}
    const auto global=draftGlobalTrim(request);
    if(global.status!=TrimStatus::noSolution)
    {out.reason="Local processing is unavailable: preserve no-change/global-only decisions or resolve invalid measurements first.";return out;}
    if(!map.valid()||peaks.geometry!=sectionGeometryKey(map)||peaks.snapshotIdentity!=request.current.snapshotIdentity
        ||peaks.sections.size()!=map.sections.size()||peaks.channels<1||peaks.channels>measurement::channelCapacity)
    {out.reason="Section peak attribution is missing or belongs to different measured input/geometry.";return out;}
    const double loudness=*request.current.integrated.value;
    const auto& limits=request.profile.limits;
    const double required=limits.integratedLufs-loudness;
    const double low=std::max(-60-request.currentGainDb,limits.integratedLufs-limits.integratedToleranceLu-loudness);
    const double high=std::min(24-request.currentGainDb,limits.integratedLufs+limits.integratedToleranceLu-loudness);
    if(low>high)
    {out.reason="No solution within this policy: required global correction exceeds the supported -60 to +24 dB total gain range.";return out;}
    out.plan.sampleRate=map.sampleRate;out.plan.frames=map.frames;out.plan.inputIdentity=map.inputIdentity;
    if(existing&&(existing->validationError().isNotEmpty()||existing->inputIdentity!=map.inputIdentity
        ||existing->sampleRate!=map.sampleRate||existing->frames!=map.frames))
    {out.reason="Existing section processing does not match this programme.";return out;}
    for(size_t n=0;n<map.sections.size();++n)
    {
        const auto& row=map.sections[n];double original=0;
        if(existing)
        {
            const auto it=std::lower_bound(existing->intervals.begin(),existing->intervals.end(),row.end,[](const auto& i,int64_t end){return i.end<end;});
            if(it==existing->intervals.end()||it->begin>row.begin){out.reason="Annotation geometry cannot inherit the applied envelope. No manual gains remapped.";return out;}
            original=it->db;
        }
        out.plan.intervals.push_back({row.begin,row.end,original});
    }
    const auto before=out.plan;
    const double additional=policy.additionalTenths/10.,cumulative=policy.cumulativeTenths/10.,reserve=policy.reserveTenths/10.;
    for(const auto& row:before.intervals)if(row.db < -cumulative)
    {out.reason="Local proposal unavailable: existing manual envelope exceeds selected cumulative attenuation. Settings unchanged; choose compatible limits or explicitly change/reverify processing.";return out;}
    if(additional==0||policy.durationThousandths==0)
    {out.reason="Local changes disabled by the selected additional attenuation or duration preference.";return out;}
    for(size_t n=0;n<map.sections.size();++n)
    {
        const auto& observed=peaks.sections[n];
        if(observed.truePeak.size()!=size_t(peaks.channels)||observed.sample.size()!=size_t(peaks.channels))
        {out.reason="Incomplete all-channel section attribution.";return out;}
        size_t ch=0;
        for(size_t k=0;k<observed.truePeak.size();++k)
        {
            if(!std::isfinite(observed.truePeak[k].magnitude)||observed.truePeak[k].magnitude<0
                ||!std::isfinite(observed.sample[k].magnitude)||observed.sample[k].magnitude<0)
            {out.reason="Non-finite/invalid peak attribution cannot select a correction.";return out;}
            if(observed.truePeak[k].magnitude>observed.truePeak[ch].magnitude)ch=k;
        }
        const auto& peak=observed.truePeak[ch];
        if(peak.magnitude==0)continue; // silence is not a finite zero-dB measurement
        const double db=20*std::log10(peak.magnitude),deficit=db+required-limits.truePeakCeilingDbtp;
        if(deficit<=0)continue;
        const auto& section=map.sections[n];auto& row=out.plan.intervals[n];
        SectionCorrection correction;correction.id=section.id;correction.begin=row.begin;correction.end=row.end;
        correction.channel=int(ch);correction.peakEighthFrame=peak.eighthFrame;correction.padding=peak.padding;
        correction.originalDb=row.db;correction.peakDbtp=db;correction.deficitDb=deficit;
        correction.requestedDb=deficit+reserve;correction.reductionDb=std::min({additional,std::max(0.,row.db+cumulative),correction.requestedDb});
        correction.reserveDb=std::min(reserve,std::max(0.,correction.reductionDb-deficit));
        correction.uncoveredDb=std::max(0.,deficit-correction.reductionDb);
        row.db=std::max(-cumulative,originalReductionTarget(row.db,correction.reductionDb));
        correction.reductionDb=correction.originalDb-row.db;
        correction.reserveDb=std::min(reserve,std::max(0.,correction.reductionDb-deficit));
        correction.uncoveredDb=std::max(0.,deficit-correction.reductionDb);
        correction.resultingDb=row.db;out.corrections.push_back(correction);
    }
    if(out.corrections.empty()) {out.reason="No local peak-constrained section identified; no section envelope proposed.";return out;}
    const auto error=out.plan.validationError();
    if(error.isNotEmpty()){out.reason="No solution within this policy: "+error+". Existing manual settings preserved.";return out;}
    const auto support=changedEnvelopeSupport(before,out.plan);out.processingRegions=support.size();out.affectedFrames=0;
    for(auto span:support)out.affectedFrames+=span.second-span.first;
    if(out.affectedFrames==0){out.reason="No changed local envelope is available within these gain limits. Uncovered deficits are shown; no private candidate requested.";return out;}
    juce::String limiting;
    if(policy.regionLimit&&out.processingRegions>size_t(*policy.regionLimit))
        limiting << juce::String(out.processingRegions) << " changed processing regions exceeds " << *policy.regionLimit << ". ";
    // Three-hour/192 kHz frame capacity makes both products safe in int64.
    if(out.affectedFrames*1000>map.frames*policy.durationThousandths)
        limiting << "Local envelope/transition support affects " << juce::String(100.*double(out.affectedFrames)/double(map.frames),6)
            << "%, above " << juce::String(policy.durationThousandths/10.,1) << "%. ";
    if(limiting.isNotEmpty())
    {out.reason="No solution within these preferences: "+limiting+"Complete set retained; no subset selected. These are preservation preferences, not delivery rules. Global trim extent is separate.";return out;}
    out.available=true;out.reason="One bounded envelope selected; exact private whole-programme assessment is still required.";return out;
}
juce::String sectionSeedText(const SectionTrimSeed& s)
{
    juce::MemoryOutputStream text; text << s.reason << "\n" << s.policy.text();
    text << "Changed annotation rows: " << juce::String(std::count_if(s.corrections.begin(),s.corrections.end(),[](const auto& row){return row.originalDb!=row.resultingDb;})) << "; inspected correction rows: " << juce::String(s.corrections.size()) << "; connected processing regions: " << juce::String(s.processingRegions) << "; local support "
         << (s.affectedFrames<0?juce::String("not assessed after the earlier rejection"):juce::String(s.affectedFrames)+" frames");
    if(s.affectedFrames>=0&&s.plan.sampleRate>0&&s.plan.frames>0)
        text << " / " << juce::String(double(s.affectedFrames)/s.plan.sampleRate,6) << " seconds / "
             << juce::String(100.*double(s.affectedFrames)/double(s.plan.frames),6) << "% of programme";
    text << "; global trim acts on the whole programme.\n";
    for(const auto& row:s.corrections)
    {
        text << "Section " << juce::String(row.id) << " [" << row.begin << ", " << row.end << ") source frames; channel "
             << row.channel+1 << "; TP " << juce::String(row.peakDbtp,9) << " dBTP at " << row.peakEighthFrame << "/8 frames"
             << (row.padding?" (filter padding attribution)":"") << "\n"
             << "Gain " << juce::String(row.originalDb,12) << " -> " << juce::String(row.resultingDb,12)
             << " dB; requested reduction " << juce::String(row.requestedDb,12) << ", planned " << juce::String(row.reductionDb,12)
             << "; reserve achieved " << juce::String(row.reserveDb,12) << "/" << juce::String(s.policy.reserveTenths/10.,1) << " dB; uncovered deficit " << juce::String(row.uncoveredDb,12) << " dB.\n";
    }
    if(s.plan.validationError().isEmpty())
    {
        auto unity=s.plan;for(auto& row:unity.intervals)row.db=0;
        text << "Cumulative non-unity local envelope support: " << changedEnvelopeExtent(unity,s.plan)
             << " source frames (includes pre-existing manual processing).\n";
        const auto p=s.plan.canonical();
        double largestReduction=0,cumulativeAttenuation=0;
        for(const auto& row:s.corrections)largestReduction=std::max(largestReduction,row.reductionDb);
        for(const auto& row:p.intervals)cumulativeAttenuation=std::max(cumulativeAttenuation,-row.db);
        text << "Maximum additional local reduction: " << juce::String(largestReduction,12) << " dB; maximum resulting cumulative local attenuation: " << juce::String(cumulativeAttenuation,12) << " dB.\n";
        text << "Canonical plateaus: " << juce::String(p.intervals.size()) << "; transitions: " << juce::String(p.intervals.size()-1) << "\n";
        double maximumStep=0,maximumSlew=0;
        for(size_t n=1;n<p.intervals.size();++n)
        {
            const auto step=std::abs(p.intervals[n].db-p.intervals[n-1].db);
            maximumStep=std::max(maximumStep,step);
            maximumSlew=std::max(maximumSlew,step*double(p.sampleRate)/double(2*p.radius(n)-1));
        }
        text << "Largest neighbouring plateau difference: " << juce::String(maximumStep,12)
             << " dB; largest transition slope: " << juce::String(maximumSlew,12)
             << " dB/s (descriptive, not additional acceptance rules).\n";
    }
    return text.toUTF8();
}
}

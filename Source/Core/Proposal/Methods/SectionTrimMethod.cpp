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

int64_t changedEnvelopeExtent(const processing::SectionGainPlan& before,
                             const processing::SectionGainPlan& after)
{
    if(before.validationError().isNotEmpty()||after.validationError().isNotEmpty()
        ||before.inputIdentity!=after.inputIdentity||before.frames!=after.frames||before.sampleRate!=after.sampleRate
        ||before.intervals.size()!=after.intervals.size())throw std::invalid_argument("Envelope extent requires matching valid geometry");
    std::vector<std::pair<int64_t,int64_t>> changed,extent;
    for(size_t n=0;n<before.intervals.size();++n)
    {
        const auto& a=before.intervals[n];const auto& b=after.intervals[n];
        if(a.begin!=b.begin||a.end!=b.end)throw std::invalid_argument("Envelope extent geometry changed");
        if(a.db!=b.db)changed.emplace_back(a.begin,a.end);
    }
    extent=changed;
    for(const auto* source:{&before,&after})
    {
        const auto p=source->canonical();
        for(size_t n=1;n<p.intervals.size();++n)
        {
            const auto radius=p.radius(n),begin=p.intervals[n].begin-radius,end=p.intervals[n].begin+radius;
            if(std::any_of(changed.begin(),changed.end(),[&](auto span){return span.first<end&&span.second>begin;}))
                extent.emplace_back(begin,end);
        }
    }
    // Count complete transition support conservatively, including ramp endpoints.
    std::sort(extent.begin(),extent.end());int64_t total=0,begin=0,end=0;
    for(auto span:extent)
        if(span.first>end){total+=end-begin;begin=span.first;end=span.second;}
        else end=std::max(end,span.second);
    return total+end-begin;
}
SectionTrimSeed selectSectionTrims(const TrimRequest& request,const SectionMap& map,
    const std::optional<processing::SectionGainPlan>& existing,const SectionPeakMap& peaks)
{
    SectionTrimSeed out;
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
            const auto it=std::find_if(existing->intervals.begin(),existing->intervals.end(),[&](const auto& i){return i.begin<=row.begin&&i.end>=row.end;});
            if(it==existing->intervals.end()){out.reason="Annotation geometry cannot inherit the applied envelope. No manual gains remapped.";return out;}
            original=it->db;
        }
        out.plan.intervals.push_back({row.begin,row.end,original});
    }
    const auto before=out.plan;
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
        correction.requestedDb=deficit+.5;correction.reductionDb=std::min(3.,correction.requestedDb);
        correction.reserveDb=std::min(.5,std::max(0.,correction.reductionDb-deficit));
        correction.uncoveredDb=std::max(0.,deficit-correction.reductionDb);
        row.db=originalReductionTarget(row.db,correction.reductionDb);
        correction.reductionDb=correction.originalDb-row.db;
        correction.reserveDb=std::min(.5,std::max(0.,correction.reductionDb-deficit));
        correction.uncoveredDb=std::max(0.,deficit-correction.reductionDb);
        correction.resultingDb=row.db;out.corrections.push_back(correction);
    }
    if(out.corrections.empty()) {out.reason="No local peak-constrained section identified; no section envelope proposed.";return out;}
    const auto error=out.plan.validationError();
    if(error.isNotEmpty()){out.reason="No solution within this policy: "+error+". Existing manual settings preserved.";return out;}
    if(out.corrections.size()>8)
    {out.reason="No solution within this policy: "+juce::String(out.corrections.size())+" changed sections exceeds eight. No set truncated.";return out;}
    out.affectedFrames=changedEnvelopeExtent(before,out.plan);
    if(out.affectedFrames*5>map.frames)
    {out.reason="No solution within this policy: local envelope/transition support affects "+juce::String(100.*double(out.affectedFrames)/double(map.frames),6)+"%, above 20%. Global trim extent is separate.";return out;}
    out.available=true;out.reason="One bounded envelope selected; exact private whole-programme assessment is still required.";return out;
}
juce::String sectionSeedText(const SectionTrimSeed& s)
{
    juce::String text=s.reason+"\n";
    text << "Changed sections: " << juce::String(s.corrections.size()) << " (not peak-event counts); local support "
         << (s.affectedFrames<0?juce::String("not assessed after the earlier rejection"):juce::String(s.affectedFrames)+" frames")
         << "; global trim acts on the whole programme.\n";
    size_t displayed=0;
    for(const auto& row:s.corrections)
    {
        if(displayed++==32){text << "Further row details omitted from this bounded text view; all " << juce::String(s.corrections.size()) << " corrections remain in full-resolution data. The complete set was rejected, not truncated into a plan.\n";break;}
        text << "Section " << juce::String(row.id) << " [" << row.begin << ", " << row.end << ") source frames; channel "
             << row.channel+1 << "; TP " << juce::String(row.peakDbtp,9) << " dBTP at " << row.peakEighthFrame << "/8 frames"
             << (row.padding?" (filter padding attribution)":"") << "\n"
             << "Gain " << juce::String(row.originalDb,12) << " -> " << juce::String(row.resultingDb,12)
             << " dB; requested reduction " << juce::String(row.requestedDb,12) << ", applied " << juce::String(row.reductionDb,12)
             << "; reserve achieved " << juce::String(row.reserveDb,12) << "/0.5 dB; uncovered deficit " << juce::String(row.uncoveredDb,12) << " dB.\n";
    }
    if(s.plan.validationError().isEmpty())
    {
        auto unity=s.plan;for(auto& row:unity.intervals)row.db=0;
        text << "Cumulative non-unity local envelope support: " << changedEnvelopeExtent(unity,s.plan)
             << " source frames (includes pre-existing manual processing).\n";
        const auto p=s.plan.canonical();
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
    return text;
}
}

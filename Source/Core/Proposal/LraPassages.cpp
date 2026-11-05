#include "LraPassages.h"
#include "Core/Measurement/GatePopulation.h"
#include <iomanip>
#include <sstream>

namespace levelscope::proposal
{
namespace
{
using Plan=processing::SectionGainPlan;
void require(bool b,const char* message){if(!b)throw std::invalid_argument(message);}
std::vector<LraPassage> group(const ProgramAnalysis& a,double level,bool upper,const std::function<bool()>& cancelled)
{
    std::vector<LraPassage> out;std::vector<double> values;
    const auto finish=[&]{if(!values.empty()){std::sort(values.begin(),values.end());out.back().p95=values[measurement::percentileIndex(values.size(),.95)];values.clear();}};
    size_t index=0;
    for(const auto& w:a.lraWindows)
    {
        if((index++%256)==0&&cancelled())return {};
        if(!w.included||!w.lufs||(upper?*w.lufs<=level:*w.lufs>level))continue;
        const auto begin=std::clamp<int64_t>(w.begin,0,a.frames),end=std::clamp<int64_t>(w.end,0,a.frames);
        if(begin>=end)continue;
        if(out.empty()||begin>out.back().end){finish();out.push_back({begin,end});}
        auto& row=out.back();row.end=std::max(row.end,end);++row.windows;
        row.paddedWindows+=w.padding?1:0;row.score+=upper?*w.lufs-level:level-*w.lufs;values.push_back(*w.lufs);
    }
    finish();
    std::sort(out.begin(),out.end(),[](const auto& x,const auto& y)
    {if(x.score!=y.score)return x.score>y.score;if(x.windows!=y.windows)return x.windows>y.windows;if(x.begin!=y.begin)return x.begin<y.begin;return x.end<y.end;});
    return out;
}
Plan overlay(const Plan& before,const std::vector<LraRegionGain>& changes)
{
    std::vector<int64_t> cuts{0,before.frames};for(const auto& i:before.intervals)cuts.push_back(i.end);
    for(const auto& c:changes){cuts.push_back(c.passage.begin);cuts.push_back(c.passage.end);}
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    Plan p=before;p.intervals.clear();size_t old=0;
    auto sorted=changes;std::sort(sorted.begin(),sorted.end(),[](const auto& a,const auto& b){return a.passage.begin<b.passage.begin;});size_t change=0;
    for(size_t n=1;n<cuts.size();++n)
    {
        const auto begin=cuts[n-1],end=cuts[n];while(before.intervals[old].end<=begin)++old;
        while(change<sorted.size()&&sorted[change].passage.end<=begin)++change;
        const double original=before.intervals[old].db;
        const double reduction=change<sorted.size()&&sorted[change].passage.begin<=begin?sorted[change].reduction:0;
        double gain=original-reduction;
        if(original-gain>reduction)gain=std::nextafter(gain,original);
        p.intervals.push_back({begin,end,gain});
    }
    return p.canonical();
}
double capFor(const Plan& before,const LraPassage& r,const ProposalPolicy& policy)
{
    double cap=policy.additionalTenths/10.;
    for(const auto& i:before.intervals)if(i.begin<r.end&&r.begin<i.end)cap=std::min(cap,std::max(0.,i.db+policy.cumulativeTenths/10.));
    return cap;
}
void updateRow(LraRegionGain& r,double reserve)
{
    const double deficit=std::max(0.,r.requested-reserve);
    r.reserve=std::min(reserve,std::max(0.,r.reduction-deficit));r.uncovered=std::max(0.,deficit-r.reduction);
}
}
size_t LraPassages::storageBytes()const
{return sizeof(*this)+(analysis?analysis->storageBytes():0)+(low.capacity()+high.capacity())*sizeof(LraPassage);}
std::shared_ptr<const LraPassages> findLraPassages(std::shared_ptr<const ProgramAnalysis> a,std::optional<double> ceiling,const std::function<bool()>& cancelled)
{
    require(bool(a),"Passages require validated measured diagnostics");
    require(!ceiling||(std::isfinite(*ceiling)&&*ceiling>=0&&*ceiling<=100),"Invalid LRA search ceiling");
    if(cancelled())return {};
    auto r=std::make_shared<LraPassages>();r->analysis=std::move(a);r->ceiling=ceiling;
    std::ostringstream key;key.imbue(std::locale::classic());key<<lraPassagesVersion<<';'<<r->analysis->identity.toStdString()<<';'<<std::hexfloat;
    if(ceiling)key<<*ceiling;else key<<"descriptive";auto bytes=key.str();r->identity=juce::SHA256(bytes.data(),bytes.size()).toHexString();
    for(const auto& w:r->analysis->lraWindows){if(!w.lufs)++r->unavailable;else if(w.included)++r->included;else ++r->excluded;}
    if(!r->analysis->p10||!r->analysis->p95||!r->analysis->lra.value)
    {r->notice="LRA population unavailable. Silence/excluded windows are not zero-valued measurements.";return r;}
    r->upperLevel=ceiling?*r->analysis->p10+*ceiling:*r->analysis->p95;
    r->low=group(*r->analysis,*r->analysis->p10,false,cancelled);
    r->high=group(*r->analysis,*r->upperLevel,true,cancelled);if(cancelled())return {};
    r->notice="Passages to investigate, not additive LRA contributions or defective scenes. Overlapping three-second windows are grouped by actual support; EOF padding is labelled.";
    if(!r->analysis->lraStable)r->notice+=" Short programme: available LRA is not stable enough for a changing LRA proposal.";
    if(!ceiling)r->notice+=" Descriptive only: no job LRA ceiling selected.";
    return r;
}
juce::String validateLraEnvelope(const Plan& original,const Plan& next,const ProposalPolicy& policy,int64_t& affected,size_t& regions,double& maximum)
{
    affected=0;regions=0;maximum=0;
    if(policy.validationError().isNotEmpty())return "Policy: "+policy.validationError();
    if(original.validationError().isNotEmpty()||next.validationError().isNotEmpty())return "Method: unsupported envelope geometry/capacity.";
    const auto a=original.canonical(),b=next.canonical();
    if(a.inputIdentity!=b.inputIdentity||a.frames!=b.frames||a.sampleRate!=b.sampleRate)return "Method: envelope input mismatch.";
    for(const auto& i:b.intervals)if(i.db< -policy.cumulativeTenths/10.)return "Policy: cumulative local attenuation cap.";
    std::vector<int64_t> points{0,a.frames-1};
    for(const auto* p:{&a,&b})for(size_t n=1;n<p->intervals.size();++n)
    {
        auto begin=p->intervals[n].begin-p->radius(n),end=p->intervals[n].begin+p->radius(n);
        for(auto f:{begin-1,begin,end-1,end})if(f>=0&&f<a.frames)points.push_back(f);
    }
    for(auto f:points)
    {
        const double delta=a.gainAt(f)-b.gainAt(f);maximum=std::max(maximum,delta);
        if(delta<0)return "Method: existing transition geometry would create a local boost; no neighbouring adjustment allowed.";
        if(delta>policy.additionalTenths/10.)return "Policy: actual transition attenuation exceeds additional gain cap.";
    }
    const auto spans=changedEnvelopeSupport(a,b);regions=spans.size();for(auto span:spans)affected+=span.second-span.first;
    if(policy.regionLimit&&regions>size_t(*policy.regionLimit))return "Policy: changed processing-region count limit.";
    if(affected*1000>a.frames*policy.durationThousandths)return "Policy: affected-duration limit (including complete transition support).";
    return {};
}
LraTrimSelection selectLraTrims(const LraPassages& passages,const Plan& before,const ProposalPolicy& policy,const std::function<bool()>& cancelled)
{
    LraTrimSelection out;out.before=before.canonical();out.plan=out.before;
    if(policy.validationError().isNotEmpty()||before.validationError().isNotEmpty())
    {out.restrictions.push_back("Invalid policy or baseline processing plan.");return out;}
    if(!passages.ceiling||!passages.upperLevel||!passages.analysis->lraStable)
    {out.restrictions.push_back("Method: a selected LRA ceiling and stable measured population are required.");return out;}
    for(const auto& i:before.intervals)if(i.db< -policy.cumulativeTenths/10.)
    {out.restrictions.push_back("Policy: existing attenuation already exceeds cumulative preference.");return out;}
    for(const auto& passage:passages.high)
    {
        if(cancelled())return {};
        LraRegionGain row;row.passage=passage;row.requested=std::max(0.,passage.p95-*passages.upperLevel)+policy.reserveTenths/10.;
        row.reduction=std::min(row.requested,capFor(out.before,passage,policy));updateRow(row,policy.reserveTenths/10.);
        const auto label="["+juce::String(passage.begin)+", "+juce::String(passage.end)+") ";
        if(row.reduction<=0){out.restrictions.push_back(label+"Policy: additional/cumulative gain cap leaves no reduction.");continue;}
        auto trial=out.selected;trial.push_back(row);auto plan=overlay(out.before,trial);
        int64_t affected=0;size_t count=0;double maximum=0;
        const auto error=validateLraEnvelope(out.before,plan,policy,affected,count,maximum);
        if(error.isNotEmpty()){out.restrictions.push_back(label+error);continue;}
        if(plan.audioKey()==out.plan.audioKey())continue;
        out.plan=std::move(plan);out.selected=std::move(trial);out.affectedFrames=affected;out.processingRegions=count;out.maximumAdditional=maximum;
        if(row.reduction<row.requested)out.restrictions.push_back(label+"Policy: gain cap limits requested attenuation/reserve.");
    }
    out.available=!out.selected.empty();
    if(!out.available&&out.restrictions.empty())out.restrictions.push_back("Method: no high-population region selected by this attenuation-only rule.");
    return out;
}
LraTrimSelection refineLraTrims(const LraTrimSelection& first,double deficit,const ProposalPolicy& policy,const std::function<bool()>& cancelled)
{
    auto out=first;out.available=false;
    if(!std::isfinite(deficit)||deficit<=0)return out;
    for(auto& row:out.selected)
    {
        if(cancelled())return {};
        row.requested+=deficit;row.reduction=std::min(row.reduction+deficit,capFor(first.before,row.passage,policy));updateRow(row,policy.reserveTenths/10.);
    }
    out.plan=overlay(first.before,out.selected);
    const auto error=validateLraEnvelope(first.before,out.plan,policy,out.affectedFrames,out.processingRegions,out.maximumAdditional);
    if(error.isNotEmpty())out.restrictions.push_back("Refinement: "+error);
    else out.available=out.plan.audioKey()!=first.plan.audioKey();
    if(!out.available)out.restrictions.push_back("Policy/method: no deeper valid envelope within the same regions and original limits.");
    return out;
}
juce::String LraTrimSelection::text()const
{
    juce::String s="Planning estimates only. Connected changed regions "+juce::String(processingRegions)+"; actual envelope/transition support "+juce::String(affectedFrames)+" frames";
    if(before.frames>0)s+=" / "+juce::String(100.*double(affectedFrames)/before.frames,6)+"%";
    s+="; maximum additional reduction "+juce::String(maximumAdditional,12)+" dB.\n";
    for(const auto& r:selected)s+="Region ["+juce::String(r.passage.begin)+", "+juce::String(r.passage.end)+") source frames; score "+juce::String(r.passage.score,12)+" LU*windows; "+juce::String(r.passage.windows)+" windows; requested "+juce::String(r.requested,12)+" dB, reduction "+juce::String(r.reduction,12)+", reserve "+juce::String(r.reserve,12)+", uncovered overshoot "+juce::String(r.uncovered,12)+" dB.\n";
    for(const auto& why:restrictions)s+=why+"\n";
    return s;
}
}

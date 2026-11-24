#include "RegionPlanning.h"
#include "Core/Measurement/GatePopulation.h"
#include <cmath>
#include <stdexcept>
#include <tuple>
#include <limits>
#include <numeric>

namespace levelscope::proposal
{
namespace
{
constexpr double pi=3.1415926535897932384626433832795;
void require(bool ok,const char* why){if(!ok)throw std::invalid_argument(why);}
using Plan=processing::SectionGainPlan;
std::vector<int64_t> pieces(const Plan& a,const Plan& b)
{
    std::vector<int64_t> out{0,a.frames};
    for(const auto* p:{&a,&b})for(size_t n=1;n<p->intervals.size();++n)
    {out.push_back(p->intervals[n].begin-p->radius(n));out.push_back(p->intervals[n].begin+p->radius(n));}
    std::sort(out.begin(),out.end());out.erase(std::unique(out.begin(),out.end()),out.end());return out;
}
double ratioSum(int64_t begin,int64_t end,const Plan& before,const Plan& after,double globalDelta,const std::vector<int64_t>& edges)
{
    double sum=0;auto edge=std::upper_bound(edges.begin(),edges.end(),begin);
    while(begin<end)
    {
        const auto stop=edge==edges.end()?end:std::min(end,*edge);const double count=double(stop-begin);
        const double first=(after.gainAt(begin)-before.gainAt(begin)+globalDelta)*std::log(10.)/10;
        const double last=(after.gainAt(stop-1)-before.gainAt(stop-1)+globalDelta)*std::log(10.)/10;
        const double step=count>1?(last-first)/(count-1):0;
        sum+=std::exp(first)*(step==0?count:std::expm1(step*count)/std::expm1(step));begin=stop;if(edge!=edges.end()&&begin==*edge)++edge;
    }
    return sum;
}
auto preservation(const RegionReviewScore& s)
{
    const auto splits=std::accumulate(s.activity.splits.begin(),s.activity.splits.end(),size_t(0));
    return std::tuple<int64_t,double,double,size_t,double>{s.activity.affectedFrames,s.activity.absoluteDbFrames,s.activity.maximumAdditional,splits,std::abs(s.activity.globalDelta)};
}
void deficit(RegionReviewScore& s,std::optional<double> value,double upper,std::optional<double> lower={},bool lra=false)
{
    if(!value||!std::isfinite(*value)){s.complete=false;return;}
    const auto excess=std::max(0.,*value-upper)+(lower?std::max(0.,*lower-*value):0.);
    if(excess>0)++s.failedRules;
    if(lra)s.lraExcess=excess;
    else {s.otherMaximum=std::max(s.otherMaximum,excess);s.otherSum+=excess;}
}
}
size_t PredictionSidecar::storageBytes()const
{return sizeof(*this)+bins.capacity()*sizeof(PredictionBin)+roles.capacity()*sizeof(measurement::Role);}
double PredictionCollector::Filter::process(double x)
{const double y=b0*x+z1;z1=b1*x-a1*y+z2;z2=b2*x-a2*y;return y;}
PredictionCollector::PredictionCollector(int rate,std::vector<measurement::Role> roles,int64_t frames,juce::String binding,juce::String snapshot,juce::String plan)
{
    require(measurement::OfflineMeasurement::supports(rate,roles)&&(roles.size()==1||roles.size()==2||roles.size()==6)&&frames>0&&frames<=int64_t(rate)*10800,"Unsupported prediction format/duration");
    data=std::make_shared<PredictionSidecar>();data->sampleRate=rate;data->channels=int(roles.size());data->roles=roles;data->frames=frames;data->binding=std::move(binding);data->baselineSnapshot=std::move(snapshot);data->baselinePlan=std::move(plan);
    data->bins.reserve(size_t((frames+rate/10-1)/(rate/10)+17));filters.resize(roles.size());
    // Deliberately independent auxiliary copy of the qualified m7.1 coefficient
    // expressions. No authoritative meter refactor or observer/state mutation.
    for(size_t n=0;n<roles.size();++n)
    {
        auto& c=filters[n];using R=measurement::Role;c.weight=roles[n]==R::lfe?0:(roles[n]==R::surroundLeft||roles[n]==R::surroundRight?1.41:1);
        if(rate==48000){c.shelf={1.53512485958697,-2.69169618940638,1.19839281085285,-1.69065929318241,.73248077421585};c.highpass={1,-2,1,-1.99004745483398,.99007225036621};}
        else
        {
            double k=std::tan(pi*1681.974450955533/rate),q=.7071752369554196;
            const double vh=std::pow(10.,3.999843853973347/20),vb=std::pow(vh,.4996667741545416);double a0=1+k/q+k*k;
            c.shelf={(vh+vb*k/q+k*k)/a0,2*(k*k-vh)/a0,(vh-vb*k/q+k*k)/a0,2*(k*k-1)/a0,(1-k/q+k*k)/a0};
            k=std::tan(pi*38.13547087602444/rate);q=.5003270373238773;a0=1+k/q+k*k;c.highpass={1,-2,1,2*(k*k-1)/a0,(1-k/q+k*k)/a0};
        }
    }
}
void PredictionCollector::append(double energy,const float* const* input,int offset,bool drain)
{
    require(std::isfinite(energy)&&energy>=0,"Prediction energy overflow/non-finite");
    const auto step=data->sampleRate/10;
    if(data->bins.empty()||position==data->bins.back().end||data->bins.back().drain!=drain)
    {
        PredictionBin bin;bin.begin=position;bin.end=std::min<int64_t>((position/step+1)*step,drain?data->frames+data->sampleRate*3/2:data->frames);bin.drain=drain;data->bins.push_back(bin);
    }
    auto& bin=data->bins.back();bin.weightedEnergy+=energy;
    require(std::isfinite(bin.weightedEnergy),"Prediction bin overflow");
    if(!drain)for(int ch=0;ch<data->channels;++ch)bin.sampleMax[size_t(ch)]=std::max(bin.sampleMax[size_t(ch)],std::abs(double(input[ch][offset])));
    ++position;
}
void PredictionCollector::push(const float* const* input,int count)
{
    require(!finished&&input&&count>=0&&count<=data->frames-position,"Invalid prediction collector chunk");
    for(int ch=0;ch<data->channels;++ch)require(input[ch]!=nullptr,"Missing prediction channel");
    for(int n=0;n<count;++n)
    {
        double energy=0;
        for(size_t ch=0;ch<filters.size();++ch){const double x=input[ch][n];require(std::isfinite(x),"Non-finite prediction input");auto& c=filters[ch];const double y=c.highpass.process(c.shelf.process(x));energy+=c.weight*y*y;}
        append(energy,input,n,false);
    }
}
std::shared_ptr<const PredictionSidecar> PredictionCollector::finish()
{
    if(finished)return data;require(position==data->frames,"Incomplete auxiliary analysis");
    for(int n=0;n<data->sampleRate*3/2;++n)
    {double energy=0;for(auto& c:filters){const double y=c.highpass.process(c.shelf.process(0));energy+=c.weight*y*y;}append(energy,nullptr,0,true);}
    finished=true;return data;
}

RegionPrediction predictRegionPlan(const ProgramAnalysis& measured,const PredictionSidecar& bins,const Plan& before,double beforeGlobal,const RegionPlan& plan,const std::function<bool()>& cancel)
{
    RegionPrediction out;
    require(plan.error.isEmpty()&&before.validationError().isEmpty()&&plan.envelope.validationError().isEmpty(),"Invalid prediction plan");
    require(measured.snapshotIdentity==bins.baselineSnapshot&&bins.baselinePlan==regionPlanIdentity(before,beforeGlobal)&&bins.frames==measured.frames&&bins.sampleRate==measured.sampleRate,"Prediction sidecar belongs to another exact baseline");
    out.identical=regionPlanIdentity(before,beforeGlobal)==regionPlanIdentity(plan.envelope,plan.globalDb);
    out.limitations="Approximate window-energy proxy: uniform bin energy, filter/gain noncommutation and float-stage rounding. Sampled M/S cannot certify maxima; true peak unavailable. Empirical errors are not proven bounds.";
    if(out.identical)
    {
        out.integrated=measured.integrated.value;out.lra=measured.lra.value;out.p10=measured.p10;out.p95=measured.p95;out.integratedGate=measured.integratedGate;out.lraGate=measured.lraGate;
        out.momentarySampled=measured.momentaryMax.value;out.shortTermSampled=measured.shortTermMax.value;return out;
    }
    const auto edges=pieces(before,plan.envelope);const double global=plan.globalDb-beforeGlobal;
    std::vector<double> deltaPrefix(bins.bins.size()+1,0),delta(bins.bins.size(),0);
    for(size_t n=0;n<bins.bins.size();++n)
    {
        if((n%128)==0&&cancel()){out.cancelled=true;return out;}const auto& bin=bins.bins[n];const double count=double(bin.end-bin.begin);
        double scale=1;
        if(bin.drain)
        {
            // Global scaling has a defined linear-filter tail. Local changes may
            // alter IIR history: retain baseline tail shape and flag uncertainty.
            scale=std::pow(10.,global/10.);
            out.eofUncertain=before.audioKey()!=plan.envelope.audioKey();
        }
        else scale=ratioSum(bin.begin,bin.end,before,plan.envelope,global,edges)/count;
        delta[n]=bin.weightedEnergy*(scale-1);deltaPrefix[n+1]=deltaPrefix[n]+delta[n];
    }
    const auto deltaFor=[&](int64_t begin,int64_t end)
    {
        auto first=std::upper_bound(bins.bins.begin(),bins.bins.end(),begin,[](int64_t f,const auto& b){return f<b.end;});
        auto last=std::lower_bound(bins.bins.begin(),bins.bins.end(),end,[](const auto& b,int64_t f){return b.begin<f;});
        if(first==bins.bins.end()||first>=last)return 0.;size_t i=size_t(first-bins.bins.begin()),j=size_t(last-bins.bins.begin());double sum=deltaPrefix[j]-deltaPrefix[i];
        const auto partial=[&](size_t n,int64_t b,int64_t e)
        {const auto& bin=bins.bins[n];if(b>=e)return 0.;if(bin.drain)return delta[n]*double(e-b)/double(bin.end-bin.begin);return bin.weightedEnergy*(ratioSum(b,e,before,plan.envelope,global,edges)-double(e-b))/double(bin.end-bin.begin);};
        if(begin>first->begin)sum-=partial(i,first->begin,begin);
        if(end<bins.bins[j-1].end)sum-=partial(j-1,end,bins.bins[j-1].end);return sum;
    };
    const auto transform=[&](const std::vector<DiagnosticWindow>& windows)
    {
        std::vector<measurement::Window> result;result.reserve(windows.size());size_t n=0;
        for(const auto& w:windows)
        {
            if((n++%256)==0&&cancel()){out.cancelled=true;return std::vector<measurement::Window>{};}
            const double e=w.energy+deltaFor(w.begin,w.end)/double(w.end-w.begin);require(std::isfinite(e),"Prediction energy overflow");result.push_back({w.end,std::max(0.,e)});
        }
        return result;
    };
    auto m=transform(measured.momentary),s=transform(measured.shortTerm),lra=transform(measured.lraWindows);if(out.cancelled)return out;
    const auto integratedGate=measurement::relativeGateEnergy(m,false,cancel),lraGate=measurement::relativeGateEnergy(lra,true,cancel);
    if(cancel()){out.cancelled=true;return out;}
    if(integratedGate)
    {
        out.integratedGate=measurement::OfflineMeasurement::loudness(*integratedGate);double sum=0;size_t count=0;
        for(size_t n=0;n<m.size();++n){const bool included=measurement::includedInPopulation(m[n].energy,*integratedGate,false);if(included){sum+=m[n].energy;++count;}out.integratedEntered+=included&&!measured.momentary[n].included;out.integratedLeft+=!included&&measured.momentary[n].included;}
        if(count&&sum>0)out.integrated=measurement::OfflineMeasurement::loudness(sum/count);
    }
    if(lraGate)
    {
        out.lraGate=measurement::OfflineMeasurement::loudness(*lraGate);std::vector<double> values;
        for(size_t n=0;n<lra.size();++n){const bool included=measurement::includedInPopulation(lra[n].energy,*lraGate,true);if(included)values.push_back(measurement::OfflineMeasurement::loudness(lra[n].energy));out.lraEntered+=included&&!measured.lraWindows[n].included;out.lraLeft+=!included&&measured.lraWindows[n].included;}
        if(!values.empty()){std::sort(values.begin(),values.end());out.p10=values[measurement::percentileIndex(values.size(),.1)];out.p95=values[measurement::percentileIndex(values.size(),.95)];out.lra=*out.p95-*out.p10;}
    }
    const auto max=[](const auto& windows)->std::optional<double>{double e=0;for(auto w:windows)e=std::max(e,w.energy);if(e>0)return measurement::OfflineMeasurement::loudness(e);return {};};
    out.momentarySampled=max(m);out.shortTermSampled=max(s);out.gateSensitive=out.integratedEntered||out.integratedLeft||out.lraEntered||out.lraLeft;
    if(out.eofUncertain)out.limitations+=" EOF drain retains global-scaled baseline shape; locally transformed filter-tail state is uncertain, not measured silence.";
    if(out.gateSensitive)out.limitations+=" Full-population gate membership changes; constant-gain/local-LRA shortcuts are invalid.";return out;
}

RegionReviewScore regionMeasuredScore(const Evaluation& e,const RegionActivity& activity,const juce::String& id)
{
    RegionReviewScore s;s.activity=activity;s.id=id;s.complete=e.verdict==Verdict::pass||e.verdict==Verdict::fail;bool lra=false;
    for(const auto& r:e.rules)
    {
        if(!r.upper||r.verdict==Verdict::notEvaluated||r.verdict==Verdict::indeterminate){s.complete=false;continue;}
        deficit(s,r.measured,*r.upper,r.lower,r.id=="loudness-range");lra|=r.id=="loudness-range";
    }
    s.complete&=lra;s.allPass=s.complete&&e.verdict==Verdict::pass&&s.failedRules==0;return s;
}
RegionReviewScore regionPredictedScore(const RegionPrediction& p,const TargetProfile& t,const RegionActivity& activity,const juce::String& id)
{
    RegionReviewScore s;s.complete=true;s.activity=activity;s.id=id;const auto& l=t.limits;
    deficit(s,p.integrated,l.integratedLufs+l.integratedToleranceLu,l.integratedLufs-l.integratedToleranceLu);
    deficit(s,p.lra,l.maxLoudnessRangeLu,{},true);
    if(l.hasMaxMomentary&&p.momentarySampled)deficit(s,p.momentarySampled,l.maxMomentaryLufs);
    if(l.hasMaxShortTerm&&p.shortTermSampled)deficit(s,p.shortTermSampled,l.maxShortTermLufs);
    // Unknown TP and sampled M/S do not confer pass authority. A loose upper bound
    // exceeding the ceiling is inconclusive and never a numerical failure score.
    s.allPass=false;return s;
}
std::optional<size_t> selectRegionReview(const std::vector<RegionReviewScore>& scores,double band)
{
    require(std::isfinite(band)&&band>=0,"Invalid review shortlist band");std::optional<size_t> best;bool passing=false;double minimum=std::numeric_limits<double>::infinity();
    for(const auto& s:scores)if(s.complete){passing|=s.allPass;minimum=std::min(minimum,s.lraExcess);}
    for(size_t n=0;n<scores.size();++n)
    {
        const auto& s=scores[n];if(!s.complete||(passing?!s.allPass:s.lraExcess>minimum+band))continue;
        if(!best){best=n;continue;}const auto& b=scores[*best];
        if(passing?std::tuple{preservation(s),s.id}<std::tuple{preservation(b),b.id}:std::tuple{s.failedRules,s.otherMaximum,s.otherSum,preservation(s),s.lraExcess,s.id}<std::tuple{b.failedRules,b.otherMaximum,b.otherSum,preservation(b),b.lraExcess,b.id})best=n;
    }
    return best;
}
}

#include "RegionPlanning.h"
#include "Core/Measurement/GatePopulation.h"
#include <juce_cryptography/juce_cryptography.h>
#include <numeric>
#include <set>
#include <sstream>
#include <iomanip>
#include <limits>

namespace levelscope::proposal
{
namespace
{
using Plan=processing::SectionGainPlan;
using V=juce::var;
void require(bool yes,const char* why){if(!yes)throw std::invalid_argument(why);}
juce::String hash(const juce::String& s){return juce::SHA256(s.toRawUTF8(),s.getNumBytesAsUTF8()).toHexString();}
bool overlap(RegionSpan a,RegionSpan b){return a.begin<b.end&&b.begin<a.end;}
int64_t durationFrames(int ms,int rate){return (int64_t(ms)*rate+999)/1000;}
V object(){return V(new juce::DynamicObject());}
void put(V& o,const char* key,const V& value){o.getDynamicObject()->setProperty(key,value);}
V span(RegionSpan s){auto v=object();put(v,"begin",juce::String(s.begin));put(v,"end",juce::String(s.end));return v;}
V range(RegionEdgeRange r){return span({r.minimum,r.maximum});}
int64_t integer(const V& v)
{
    require(v.isString(),"Region frame/ID must be an exact decimal string");
    auto s=v.toString();require(s.isNotEmpty()&&s.length()<=19&&s.containsOnly("0123456789"),"Invalid exact region coordinate");
    const auto n=s.getLargeIntValue();require(n>=0&&juce::String(n)==s,"Region coordinate overflow or noncanonical integer");return n;
}
RegionSpan span(const V& v){require(v.isObject()&&v.hasProperty("begin")&&v.hasProperty("end")&&v.getDynamicObject()->getProperties().size()==2,"Missing or unknown region span fields");return {integer(v["begin"]),integer(v["end"])};}
RegionEdgeRange range(const V& v){auto s=span(v);return {s.begin,s.end};}
int number(const V& v){require(v.isInt(),"Region control must be an integer");return int(v);}
bool flag(const V& v){require(v.isBool(),"Region permission must be explicit");return bool(v);}
bool hashValid(const juce::String& s){return s.length()==64&&s.containsOnly("0123456789abcdef");}
std::vector<int64_t> envelopePoints(const Plan& a,const Plan& b)
{
    std::vector<int64_t> points{0,a.frames};
    for(const auto* p:{&a,&b})for(size_t n=1;n<p->intervals.size();++n)
    {points.push_back(p->intervals[n].begin-p->radius(n));points.push_back(p->intervals[n].begin+p->radius(n));}
    std::sort(points.begin(),points.end());points.erase(std::unique(points.begin(),points.end()),points.end());return points;
}
std::vector<double> decisions(const RegionAssignment& a)
{
    auto d=a.reductions;
    d.erase(std::unique(d.begin(),d.end()),d.end());
    while(!d.empty()&&d.front()==0)d.erase(d.begin());
    while(!d.empty()&&d.back()==0)d.pop_back();return d;
}
Plan overlay(const Plan& before,const std::vector<RegionAssignment>& assignments)
{
    struct Delta{int64_t begin,end;double reduction;};std::vector<Delta> deltas;
    std::vector<int64_t> cuts{0,before.frames};for(auto i:before.intervals)cuts.push_back(i.end);
    for(const auto& a:assignments)
    {
        std::vector<int64_t> edges{a.nominal.begin};edges.insert(edges.end(),a.geometry.cuts.begin(),a.geometry.cuts.end());edges.push_back(a.nominal.end);
        require(edges.size()==a.reductions.size()+1,"Correction partition count mismatch");
        for(size_t n=0;n<a.reductions.size();++n)
        {
            require(edges[n]>=0&&edges[n]<edges[n+1]&&edges[n+1]<=before.frames,"Invalid nominal correction partition");
            require(std::isfinite(a.reductions[n])&&a.reductions[n]>=0,"Invalid additional reduction");
            if(a.reductions[n]!=0){deltas.push_back({edges[n],edges[n+1],a.reductions[n]});cuts.push_back(edges[n]);cuts.push_back(edges[n+1]);}
        }
    }
    std::sort(deltas.begin(),deltas.end(),[](auto a,auto b){return a.begin<b.begin;});
    for(size_t n=1;n<deltas.size();++n)require(deltas[n-1].end<=deltas[n].begin,"Overlapping correction assignments");
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    Plan out=before;out.intervals.clear();size_t old=0,delta=0;
    for(size_t n=1;n<cuts.size();++n)
    {
        const auto begin=cuts[n-1],end=cuts[n];while(before.intervals[old].end<=begin)++old;
        while(delta<deltas.size()&&deltas[delta].end<=begin)++delta;
        const double original=before.intervals[old].db;
        const double reduction=delta<deltas.size()&&deltas[delta].begin<=begin?deltas[delta].reduction:0;
        double gain=original-reduction;if(original-gain>reduction)gain=std::nextafter(gain,original);
        out.intervals.push_back({begin,end,gain});
    }
    return out.canonical();
}
bool supported(RegionSpan s,const std::vector<PlanningRoot>& roots)
{
    auto position=s.begin;std::vector<RegionSpan> allowed;for(const auto& r:roots)allowed.push_back(r.allowed);
    std::sort(allowed.begin(),allowed.end(),[](auto a,auto b){return a.begin<b.begin;});
    for(auto a:allowed){if(a.begin>position)break;if(a.end>position)position=a.end;if(position>=s.end)return true;}return false;
}
}

RegionPlanningSettings RegionPlanningSettings::defaults(RegionMode mode,const juce::String& input)
{RegionPlanningSettings s;s.mode=mode;s.inputIdentity=input;s.allowGlobal=mode==RegionMode::automatic;return s;}
RegionPlanningSettings RegionPlanningSettings::eventDefaults(RegionMode mode,const juce::String& input)
{auto s=defaults(mode,input);s.schema=2;s.preferCompleteEvent=true;s.compareAlternatives=true;return s;}
RegionPlanningSettings RegionPlanningSettings::innerDefaults(RegionMode mode,const juce::String& input)
{auto settings=eventDefaults(mode,input);settings.schema=3;settings.eventFeatures.innerEntrances=true;return settings;}
juce::String RegionPlanningSettings::validationError()const
{
    if((schema<1||schema>3||eventFeatures.innerEntrances!=(schema==3))||!hashValid(inputIdentity)||int(mode)<0||int(mode)>2)return "Unsupported region settings version, input or mode";
    if(schema==1&&(preferCompleteEvent||allowAutomaticMargins||eventFeatures.hopMs!=10))return "Legacy region settings cannot carry event permissions";
    if(schema>=2&&(!preferCompleteEvent||!compareAlternatives||eventFeatures.validationError().isNotEmpty()))return "Event preference requires explicit valid features and comparative assessment";
    if(maximumSplits<0||maximumSplits>4||leftMarginMs<0||leftMarginMs>3000||leftMarginMs%100||rightMarginMs<0||rightMarginMs>3000||rightMarginMs%100
        ||minimumIntervalMs<500||minimumIntervalMs>10000||minimumIntervalMs%100||separationMs<0||separationMs>3000||separationMs%100
        ||needThresholdTenths<5||needThresholdTenths>30||persistenceMs<300||persistenceMs>3000||persistenceMs%100)return "Region controls outside documented range/steps";
    if(protectionMethod!="local-envelope-v1")return "Unsupported mandatory protection; only no additional local-envelope change is implemented";
    if(mode==RegionMode::automatic&&(start||end||!internalCuts.empty()))return "Automatic roots do not accept selected-root pins";
    if(mode==RegionMode::automatic&&(leftMarginMs||rightMarginMs)&&(schema==1||!allowAutomaticMargins))return "Automatic outward margins need explicit versioned permission";
    if(mode!=RegionMode::automatic&&allowAutomaticMargins)return "Automatic margin permission does not apply to selected roots";
    if(mode==RegionMode::selectedSection&&sectionId==0)return "Select an original section";
    if(mode==RegionMode::selection&&(selection.begin<0||selection.begin>=selection.end))return "Invalid temporary selection";
    if(internalCuts.size()>4||(!allowSplitting&&!internalCuts.empty())||internalCuts.size()>size_t(maximumSplits))return "Internal constraints exceed original-root split allowance";
    for(const auto* r:{&start,&end})if(*r&&((*r)->minimum<0||(*r)->maximum<(*r)->minimum))return "Contradictory support-edge range";
    int64_t last=-1;for(auto r:internalCuts){if(r.minimum<0||r.maximum<r.minimum||r.minimum<=last)return "Internal cut ranges must be ordered and disjoint";last=r.maximum;}
    if(protectedLocal.size()>128)return "Too many protected intervals (maximum 128)";
    for(auto p:protectedLocal)if(p.begin<0||p.end<=p.begin)return "Invalid protected local interval";
    return {};
}
V RegionPlanningSettings::encode()const
{
    require(validationError().isEmpty(),validationError().toRawUTF8());auto v=object();
    put(v,"schema",schema);put(v,"method",methodVersion());put(v,"constraints",schema>=2?eventConstraintsVersion:regionConstraintsVersion);put(v,"formation",formationVersion());put(v,"proxy",regionProxyVersion);put(v,"rank",schema>=2?eventRankVersion:regionRankVersion);
    if(schema>=2){put(v,"boundaryPreference","prefer-complete-event");put(v,"automaticMargins",allowAutomaticMargins);put(v,"eventFeatures",eventFeatures.encode());put(v,"coverage",coverageVersion());put(v,"schedule",scheduleVersion());}
    put(v,"input",inputIdentity);put(v,"mode",int(mode));put(v,"section",juce::String(sectionId));put(v,"selection",span(selection));
    put(v,"splitting",allowSplitting);put(v,"refinement",allowRefinement);put(v,"fillGap",fillShortGap);put(v,"global",allowGlobal);put(v,"compare",compareAlternatives);
    put(v,"splits",maximumSplits);put(v,"leftMs",leftMarginMs);put(v,"rightMs",rightMarginMs);put(v,"minimumMs",minimumIntervalMs);put(v,"separationMs",separationMs);put(v,"needTenths",needThresholdTenths);put(v,"persistenceMs",persistenceMs);
    put(v,"start",start?range(*start):V());put(v,"end",end?range(*end):V());
    juce::Array<V> cuts,protect;for(auto r:internalCuts)cuts.add(range(r));for(auto p:protectedLocal)protect.add(span(p));put(v,"cuts",cuts);put(v,"protection",protectionMethod);put(v,"protected",protect);return v;
}
RegionPlanningSettings RegionPlanningSettings::decode(const V& v)
{
    require(v.isObject(),"Missing region definition");
    auto s=RegionPlanningSettings{};s.schema=number(v["schema"]);require(s.schema>=1&&s.schema<=3,"Unsupported region schema");
    for(auto pair:{std::pair<const char*,const char*>{"method",s.methodVersion()},{"constraints",s.schema>=2?eventConstraintsVersion:regionConstraintsVersion},{"formation",s.formationVersion()},{"proxy",regionProxyVersion},{"rank",s.schema>=2?eventRankVersion:regionRankVersion}})
        require(v[pair.first].toString()==pair.second,"Unknown mandatory region method version");
    if(s.schema>=2){require(v["boundaryPreference"].toString()=="prefer-complete-event"&&v["coverage"].toString()==s.coverageVersion()&&v["schedule"].toString()==s.scheduleVersion(),"Unknown event treatment semantics");s.preferCompleteEvent=true;s.allowAutomaticMargins=flag(v["automaticMargins"]);s.eventFeatures=EventFeatureDefinition::decode(v["eventFeatures"]);}
    s.inputIdentity=v["input"].toString();s.mode=RegionMode(number(v["mode"]));s.sectionId=uint64_t(integer(v["section"]));s.selection=span(v["selection"]);
    s.allowSplitting=flag(v["splitting"]);s.allowRefinement=flag(v["refinement"]);s.fillShortGap=flag(v["fillGap"]);s.allowGlobal=flag(v["global"]);s.compareAlternatives=flag(v["compare"]);
    s.maximumSplits=number(v["splits"]);s.leftMarginMs=number(v["leftMs"]);s.rightMarginMs=number(v["rightMs"]);s.minimumIntervalMs=number(v["minimumMs"]);s.separationMs=number(v["separationMs"]);s.needThresholdTenths=number(v["needTenths"]);s.persistenceMs=number(v["persistenceMs"]);
    if(!v.hasProperty("start")||!v.hasProperty("end"))throw std::invalid_argument("Missing region edge permissions");
    if(!v["start"].isVoid())s.start=range(v["start"]);if(!v["end"].isVoid())s.end=range(v["end"]);
    require(v["cuts"].isArray()&&v["cuts"].size()<=4&&v["protected"].isArray()&&v["protected"].size()<=128,"Invalid region constraint arrays");
    for(const auto& r:*v["cuts"].getArray())s.internalCuts.push_back(range(r));for(const auto& p:*v["protected"].getArray())s.protectedLocal.push_back(span(p));s.protectionMethod=v["protection"].toString();
    require(s.validationError().isEmpty(),s.validationError().toRawUTF8());
    // Unknown keys may represent mandatory semantics. Never silently drop them.
    const auto expected=s.encode();
    require(v.getDynamicObject()->getProperties().size()==expected.getDynamicObject()->getProperties().size(),"Unknown region definition fields");
    for(const auto& property:v.getDynamicObject()->getProperties())require(expected.hasProperty(property.name),"Unknown region definition field");return s;
}
juce::String RegionPlanningSettings::digest()const
{if(validationError().isNotEmpty())return hash("invalid-region-settings;"+validationError());return hash(juce::JSON::toString(encode(),true));}

RegionConstraints resolveRegionConstraints(const RegionPlanningSettings& s,const SectionMap& map,const LraPassages& p,const Plan& before)
{
    RegionConstraints c;c.settings=s;c.sampleRate=before.sampleRate;c.frames=before.frames;c.error=s.validationError();
    if(c.error.isNotEmpty())return c;
    if(before.validationError().isNotEmpty()||s.inputIdentity!=before.inputIdentity||!p.analysis||p.analysis->frames!=before.frames||p.analysis->sampleRate!=before.sampleRate)
    {c.error="Region request does not describe the measured programme/baseline";return c;}
    for(auto span:s.protectedLocal)if(span.end>before.frames){c.error="Protected material lies beyond the programme";return c;}
    if(s.mode==RegionMode::automatic)
    {
        std::vector<size_t> efficiency(p.high.size());std::iota(efficiency.begin(),efficiency.end(),0);
        std::sort(efficiency.begin(),efficiency.end(),[&](size_t x,size_t y){const auto&a=p.high[x];const auto&b=p.high[y];const double u=a.score/double(a.end-a.begin),v=b.score/double(b.end-b.begin);if(u!=v)return u>v;if(a.windows!=b.windows)return a.windows>b.windows;return a.begin<b.begin;});
        std::set<size_t> used;size_t scoreIndex=0,densityIndex=0;
        while(c.roots.size()<8&&used.size()<p.high.size())
        {
            for(int order=0;order<2&&c.roots.size()<8;++order)
            {
                auto& cursor=order==0?scoreIndex:densityIndex;
                while(cursor<p.high.size()&&used.count(order==0?cursor:efficiency[cursor]))++cursor;
                if(cursor==p.high.size())continue;const size_t n=order==0?cursor:efficiency[cursor];++cursor;used.insert(n);
                const auto& seed=p.high[n];PlanningRoot r;r.id="analysis:"+p.identity+":"+juce::String(seed.begin)+":"+juce::String(seed.end);r.hint=r.allowed={seed.begin,seed.end};r.score=seed.score;c.roots.push_back(r);
            }
        }
        c.omittedRoots=p.high.size()-c.roots.size();
        // Ownership includes frozen neighbours omitted by the eight-root search cap.
        // It is resolved once, never enlarged after a root/candidate is skipped.
        for(auto& r:c.roots)
        {
            r.ownership={0,c.frames};
            for(const auto& h:p.high)
            {
                if(h.end<=r.hint.begin)r.ownership.begin=std::max(r.ownership.begin,h.end+(r.hint.begin-h.end)/2);
                if(h.begin>=r.hint.end)r.ownership.end=std::min(r.ownership.end,r.hint.end+(h.begin-r.hint.end)/2);
            }
            r.requestedAllowed=r.allowed;
            if(s.schema>=2&&s.allowAutomaticMargins)
            {
                r.requestedAllowed={std::max<int64_t>(0,r.hint.begin-int64_t(s.leftMarginMs)*c.sampleRate/1000),std::min(c.frames,r.hint.end+int64_t(s.rightMarginMs)*c.sampleRate/1000)};
                r.allowed={std::max(r.requestedAllowed.begin,r.ownership.begin),std::min(r.requestedAllowed.end,r.ownership.end)};
                if(!(r.allowed==r.requestedAllowed))r.permissionNotice="Requested outward margin clipped by frozen midpoint ownership";
            }
        }
    }
    else
    {
        PlanningRoot r;r.userSelected=true;r.hint=s.selection;
        if(s.mode==RegionMode::selectedSection)
        {
            if(!map.valid()||map.inputIdentity!=before.inputIdentity){c.error="Selected section map belongs to another programme";return c;}
            auto it=std::find_if(map.sections.begin(),map.sections.end(),[&](const auto&a){return a.id==s.sectionId;});
            if(it==map.sections.end()){c.error="Original selected section no longer exists";return c;}r.hint={it->begin,it->end};
        }
        if(r.hint.begin<0||r.hint.begin>=r.hint.end||r.hint.end>c.frames){c.error="Selection is outside the programme";return c;}
        r.id=(s.mode==RegionMode::selectedSection?"section:"+juce::String(s.sectionId):juce::String("selection"))+":"+juce::String(r.hint.begin)+":"+juce::String(r.hint.end);
        r.allowed={std::max<int64_t>(0,r.hint.begin-int64_t(s.leftMarginMs)*c.sampleRate/1000),std::min(c.frames,r.hint.end+int64_t(s.rightMarginMs)*c.sampleRate/1000)};r.requestedAllowed=r.allowed;r.ownership={0,c.frames};r.internalCuts=s.internalCuts;c.roots.push_back(r);
    }
    const auto minimum=durationFrames(s.minimumIntervalMs,c.sampleRate);
    for(auto& r:c.roots)
    {
        r.start={r.allowed.begin,r.allowed.end-minimum};r.end={r.allowed.begin+minimum,r.allowed.end};
        if(!s.allowRefinement){r.start={r.hint.begin,r.hint.begin};r.end={r.hint.end,r.hint.end};}
        if(r.userSelected)
        {
            if(s.start){r.start.minimum=std::max(r.start.minimum,s.start->minimum);r.start.maximum=std::min(r.start.maximum,s.start->maximum);}
            if(s.end){r.end.minimum=std::max(r.end.minimum,s.end->minimum);r.end.maximum=std::min(r.end.maximum,s.end->maximum);}
        }
        if(r.start.minimum>r.start.maximum||r.end.minimum>r.end.maximum){c.error="Empty effective support-edge permissions/minimum duration";return c;}
        for(auto pin:r.internalCuts)if(pin.minimum<=r.allowed.begin||pin.maximum>=r.allowed.end){c.error="Internal cut constraint outside original root permission";return c;}
    }
    juce::String key=s.digest()+";"+before.audioKey();for(const auto&r:c.roots)key+=";"+r.id+":"+juce::String(r.allowed.begin)+":"+juce::String(r.allowed.end);
    c.identity=hash(key);return c;
}

juce::String regionPlanIdentity(const Plan& p,double global)
{std::ostringstream k;k.imbue(std::locale::classic());k<<std::hexfloat<<p.audioKey()<<';'<<(global==0?0.:global);return hash(k.str());}
juce::String validateRegionPlan(const Plan& original,double baselineGlobal,const RegionConstraints& c,const ProposalPolicy& policy,RegionPlan& plan)
{
    if(c.error.isNotEmpty())return c.error;
    if(!std::isfinite(plan.globalDb)||plan.globalDb< -60||plan.globalDb>24||(!c.settings.allowGlobal&&plan.globalDb!=baselineGlobal))return "Whole-programme global trim is not permitted or outside parameter limits";
    auto& a=plan.activity;a={};a.splits.resize(c.roots.size());a.globalDelta=plan.globalDb-baselineGlobal;
    auto error=validateLraEnvelope(original,plan.envelope,policy,a.affectedFrames,a.intervalCount,a.maximumAdditional);if(error.isNotEmpty())return error;
    const auto before=original.canonical(),after=plan.envelope.canonical();a.intervalCount=after.intervals.size();a.transitionCount=after.intervals.size()-1;
    for(auto i:after.intervals)a.maximumCumulative=std::max(a.maximumCumulative,std::max(0.,-i.db));
    auto raw=changedEnvelopeSupport(before,after);for(auto [begin,end]:raw)
    {
        RegionSpan d{begin,end};if(!supported(d,c.roots))return "Complete changed support (old/new ramps included) exceeds permitted local area";
        for(auto p:c.settings.protectedLocal)if(overlap(d,p))return "Complete changed support overlaps protected local material";
        a.support.push_back(d);a.longestRun=std::max(a.longestRun,end-begin);
    }
    const auto separation=durationFrames(c.settings.separationMs,c.sampleRate);
    for(size_t n=1;n<a.support.size();++n)if(a.support[n].begin-a.support[n-1].end<separation)return "Unchanged separation after ramps is shorter than requested";
    std::set<size_t> seen;
    for(const auto& assignment:plan.assignments)
    {
        const size_t root=assignment.geometry.root;if(root>=c.roots.size()||!seen.insert(root).second)return "Missing or duplicate original planning root";
        const auto decision=decisions(assignment);if(decision.empty())continue;
        a.splits[root]=decision.size()-1;
        if(a.splits[root]>size_t(c.settings.allowSplitting?c.settings.maximumSplits:0))return "Original-root split allowance exceeded (unchanged interior parts count)";
        const auto& r=c.roots[root];std::vector<RegionSpan> support;
        for(auto d:a.support)if(overlap(d,r.allowed))support.push_back({std::max(d.begin,r.allowed.begin),std::min(d.end,r.allowed.end)});
        if(support.empty())return "Correction assignment has no supported audible geometry";
        if(!r.start.contains(support.front().begin)||!r.end.contains(support.back().end))return "Actual changed-support edges violate exact ranges/pins";
        if(support.front().begin!=assignment.geometry.support.begin||support.back().end!=assignment.geometry.support.end)return "Fixed processor cannot represent requested exact support edges";
        if(r.userSelected&&!std::any_of(support.begin(),support.end(),[&](auto d){return overlap(d,r.hint);}))return "Margin correction must intersect its original selection";
        for(auto d:support){a.marginFrames+=std::max<int64_t>(0,std::min(d.end,r.hint.begin)-d.begin)+std::max<int64_t>(0,d.end-std::max(d.begin,r.hint.end));}
        std::vector<int64_t> edges{assignment.nominal.begin};edges.insert(edges.end(),assignment.geometry.cuts.begin(),assignment.geometry.cuts.end());edges.push_back(assignment.nominal.end);
        if(edges.size()!=assignment.reductions.size()+1)return "Incomplete correction partition";
        for(auto pin:r.internalCuts)if(std::none_of(assignment.geometry.cuts.begin(),assignment.geometry.cuts.end(),[&](auto f){return pin.contains(f);}))return "Mandatory internal divider range/pin missing";
        // Merge equal correction decisions for minimum-duration accounting; baseline
        // cuts remain in the resulting audio plan and are never flattened here.
        for(size_t n=0;n<assignment.reductions.size();)
        {
            size_t last=n+1;while(last<assignment.reductions.size()&&assignment.reductions[last]==assignment.reductions[n])++last;
            if(assignment.reductions[n]!=0&&edges[last]-edges[n]<durationFrames(c.settings.minimumIntervalMs,c.sampleRate))return "Nonzero nominal correction interval shorter than minimum";n=last;
        }
    }
    if(overlay(before,plan.assignments).audioKey()!=after.audioKey())return "Plan does not match recorded root assignments";
    const auto points=envelopePoints(before,after);for(size_t n=1;n<points.size();++n)
    {
        const auto begin=points[n-1],end=points[n];if(begin>=end)continue;
        const double left=std::abs(before.gainAt(begin)-after.gainAt(begin));
        const double right=std::abs(before.gainAt(end-1)-after.gainAt(end-1));a.absoluteDbFrames+=double(end-begin)*(left+right)/2;
    }
    plan.identity=regionPlanIdentity(after,plan.globalDb);return {};
}
RegionPlan constructRegionPlan(const Plan& before,double globalBefore,const RegionConstraints& c,const ProposalPolicy& policy,std::vector<RegionAssignment> assignments,double global)
{
    RegionPlan p;p.globalDb=global;p.envelope=before.canonical();p.assignments=std::move(assignments);
    try
    {
        if(c.error.isNotEmpty()){p.error=c.error;return p;}
        const int64_t radius=std::llround(c.sampleRate*.05);
        for(auto& a:p.assignments)
        {
            if(a.geometry.root>=c.roots.size())throw std::invalid_argument("Invalid original planning root");
            a.nominal={a.geometry.support.begin==0?0:a.geometry.support.begin+radius,a.geometry.support.end==before.frames?before.frames:a.geometry.support.end-radius};
        }
        for(int attempt=0;attempt<8;++attempt)
        {
            p.envelope=overlay(before,p.assignments);require(p.envelope.validationError().isEmpty(),"Unrepresentable fixed-ramp geometry");
            auto support=changedEnvelopeSupport(before,p.envelope);bool changed=false;
            for(auto& a:p.assignments)
            {
                if(decisions(a).empty())continue;const auto& r=c.roots[a.geometry.root];
                auto first=support.end(),last=support.end();
                // Associate by nominal span, including the nearby old/new ramps.
                for(auto it=support.begin();it!=support.end();++it)if(it->first<a.nominal.end+radius&&a.nominal.begin-radius<it->second)
                {if(first==support.end())first=it;last=it;}
                if(first==support.end())throw std::invalid_argument("No changed support for nonzero assignment");
                const auto left=a.geometry.support.begin-first->first,right=a.geometry.support.end-last->second;
                if(left||right){a.nominal.begin+=left;a.nominal.end+=right;changed=true;}
                if(a.nominal.begin<r.allowed.begin||a.nominal.end>r.allowed.end||a.nominal.begin>=a.nominal.end)throw std::invalid_argument("Cannot fit fixed ramps to exact support permissions");
            }
            if(!changed)break;
        }
        p.envelope=overlay(before,p.assignments);p.error=validateRegionPlan(before,globalBefore,c,policy,p);
    }
    catch(const std::exception& e){p.error=e.what();}
    return p;
}

RegionGainNeed regionGainNeed(const Plan& before,const ProgramAnalysis& a,RegionSpan span,const ProposalPolicy& policy,double ceiling,bool low)
{
    RegionGainNeed out;out.choices={0};out.reserveRequested=policy.reserveTenths/10.;
    double headroom=policy.additionalTenths/10.;for(auto i:before.intervals)if(overlap(span,{i.begin,i.end}))headroom=std::min(headroom,std::max(0.,i.db+policy.cumulativeTenths/10.));
    out.headroom=headroom;if(!a.p10||headroom<=0)return out;std::vector<double> need;
    for(const auto& w:a.lraWindows){const auto centre=w.begin+(w.end-w.begin)/2;if(centre>=span.begin&&centre<span.end&&w.included&&w.lufs&&*w.lufs>*a.p10+ceiling)need.push_back(*w.lufs-*a.p10-ceiling);}
    if(need.empty()&&!low)return out;std::sort(need.begin(),need.end());
    out.deficit=need.empty()?0:need[measurement::percentileIndex(need.size(),.95)];out.requested=out.deficit+out.reserveRequested;
    const double desired=out.requested;
    const double r=std::min(headroom,desired);std::vector<double> values{0,r/2,r,std::min(headroom,r+.5),headroom};
    std::sort(values.begin(),values.end());values.erase(std::unique(values.begin(),values.end()),values.end());out.choices=std::move(values);return out;
}
std::vector<double> regionGainChoices(const Plan& before,const ProgramAnalysis& a,RegionSpan span,const ProposalPolicy& policy,double ceiling,bool low)
{return regionGainNeed(before,a,span,policy,ceiling,low).choices;}
}

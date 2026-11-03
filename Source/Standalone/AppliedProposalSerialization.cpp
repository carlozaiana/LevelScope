#include "AppliedProposalExplanation.h"
#include <cstring>
#include <stdexcept>
namespace levelscope::standalone
{
namespace
{
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
struct Writer
{
    std::vector<uint8_t> bytes;
    void u(uint64_t v){for(int n=0;n<8;++n)bytes.push_back(uint8_t(v>>(n*8)));bound();}
    void f(double v){require(std::isfinite(v),"Non-finite explanation value");uint64_t bits;std::memcpy(&bits,&v,8);u(bits);}
    void t(const juce::String& s){const auto raw=s.toStdString();u(raw.size());bytes.insert(bytes.end(),raw.begin(),raw.end());bound();}
    void bound(){require(bytes.size()<=AppliedProposalExplanation::maximumBytes,"Complete applied explanation exceeds its 8 MiB budget. No proposal applied.");}
};
struct Reader
{
    const std::vector<uint8_t>& bytes;JobControl& control;size_t position=0;
    uint64_t u(){require(!control.isCancelled(),"Explanation reading cancelled");require(position<=bytes.size()&&bytes.size()-position>=8,"Truncated applied explanation");uint64_t v=0;for(int n=0;n<8;++n)v|=uint64_t(bytes[position++])<<(n*8);return v;}
    int64_t i(){auto v=u();int64_t s;std::memcpy(&s,&v,8);return s;}
    int count(int maximum){const auto n=u();require(n<=uint64_t(maximum),"Explanation count outside supported bounds");return int(n);}
    bool b(){return count(1)!=0;}
    double f(){const auto bits=u();double v;std::memcpy(&v,&bits,8);require(std::isfinite(v),"Non-finite explanation value");return v;}
    juce::String t(size_t maximum=65536){const auto n=u();require(n<=maximum&&position<=bytes.size()&&n<=bytes.size()-position,"Invalid explanation text length");auto s=juce::String::fromUTF8(reinterpret_cast<const char*>(bytes.data()+position),int(n));require(s.getNumBytesAsUTF8()==n,"Invalid explanation UTF-8");position+=n;return s;}
};
void optional(Writer& w,std::optional<double> v){w.u(bool(v));if(v)w.f(*v);}
std::optional<double> optional(Reader& r){return r.b()?std::optional<double>(r.f()):std::nullopt;}
bool hash(const juce::String& s){return s.length()==64&&s.containsOnly("0123456789abcdef");}
void plan(Writer& w,const ProcessingPlan& p)
{
    w.u(p.schemaVersion);w.t(p.engineVersion);w.f(p.staticGainDb);w.t(p.graphJson);w.u(bool(p.sectionGains));
    if(p.sectionGains){const auto& g=*p.sectionGains;w.u(g.schema);w.t(g.inputIdentity);w.u(g.sampleRate);w.u(g.frames);w.u(g.intervals.size());for(const auto& i:g.intervals){w.u(i.begin);w.u(i.end);w.f(i.db);}}
}
ProcessingPlan plan(Reader& r)
{
    ProcessingPlan p;p.schemaVersion=r.count(2);p.engineVersion=r.t(128);p.staticGainDb=r.f();p.graphJson=r.t(0);
    if(r.b())
    {
        processing::SectionGainPlan g;g.schema=r.count(1);g.inputIdentity=r.t(64);g.sampleRate=r.count(192000);g.frames=r.i();const int n=r.count(10000);
        for(int i=0;i<n;++i)g.intervals.push_back({r.i(),r.i(),r.f()});
        require(g.validationError().isEmpty(),"Invalid explanation section plan");p.sectionGains=std::move(g);
    }
    require(p.staticGainDb>=-60&&p.staticGainDb<=24,"Invalid explanation global gain");
    require((p.schemaVersion==1&&p.engineVersion=="static-gain-v1"&&!p.sectionGains)||(p.schemaVersion==2&&p.engineVersion==processing::SectionGainPlan::version&&p.sectionGains),"Unsupported explanation processing version");return p;
}
void annotations(Writer& w,const proposal::SectionMap& m)
{
    w.u(!m.sections.empty());if(m.sections.empty())return;
    w.t(m.inputIdentity);w.u(m.sampleRate);w.u(m.frames);w.u(m.revision);w.u(m.nextId);w.u(m.sections.size());
    for(const auto& a:m.sections){w.u(a.id);w.u(a.begin);w.u(a.end);w.t(a.name);}
}
proposal::SectionMap annotations(Reader& r)
{
    proposal::SectionMap m;if(!r.b())return m;
    m.inputIdentity=r.t(64);m.sampleRate=r.count(192000);m.frames=r.i();m.revision=r.u();m.nextId=r.u();const int n=r.count(10000);
    for(int i=0;i<n;++i)m.sections.push_back({r.u(),r.i(),r.i(),r.t(512)});
    require(m.valid(),"Invalid explanation annotations");return m;
}
void target(Writer& w,const proposal::TargetProfile& p)
{
    w.u(p.schemaVersion);w.u(p.isSelectableForWorkflow);w.u(p.userSpecified);w.u(int(p.method));w.u(int(p.applicability));
    for(const auto* s:{&p.id,&p.displayName,&p.familyName,&p.description,&p.definitionVersion,&p.standardEdition,&p.referenceUrl,&p.integratedProvenance,&p.truePeakProvenance,&p.lraProvenance,&p.momentaryProvenance,&p.shortTermProvenance,&p.meterToleranceNote,&p.baseProfileDigest,&p.jobRequirementsDigest,&p.jobRequirementsError})w.t(*s);
    w.f(p.truePeakMeterToleranceDb);w.f(p.truePeakMeterBandwidthHz);const auto& l=p.limits;
    for(bool b:{l.hasIntegratedLoudnessTarget,l.hasMaxLoudnessRange,l.hasTruePeakCeiling,l.hasMaxMomentary,l.hasMaxShortTerm,l.valuesAreAuthoritative})w.u(b);
    for(double v:{l.integratedLufs,l.integratedToleranceLu,l.maxLoudnessRangeLu,l.truePeakCeilingDbtp,l.maxMomentaryLufs,l.maxShortTermLufs})w.f(v);
    w.t(l.validationNote);w.t(p.digest());
}
proposal::TargetProfile target(Reader& r)
{
    proposal::TargetProfile p;p.schemaVersion=r.count(2);p.isSelectableForWorkflow=r.b();p.userSpecified=r.b();p.method=proposal::MeasurementMethod(r.count(1));p.applicability=proposal::Applicability(r.count(5));
    for(auto* s:{&p.id,&p.displayName,&p.familyName,&p.description,&p.definitionVersion,&p.standardEdition,&p.referenceUrl,&p.integratedProvenance,&p.truePeakProvenance,&p.lraProvenance,&p.momentaryProvenance,&p.shortTermProvenance,&p.meterToleranceNote,&p.baseProfileDigest,&p.jobRequirementsDigest,&p.jobRequirementsError})*s=r.t();
    p.truePeakMeterToleranceDb=r.f();p.truePeakMeterBandwidthHz=r.f();auto& l=p.limits;
    for(auto* b:{&l.hasIntegratedLoudnessTarget,&l.hasMaxLoudnessRange,&l.hasTruePeakCeiling,&l.hasMaxMomentary,&l.hasMaxShortTerm,&l.valuesAreAuthoritative})*b=r.b();
    for(auto* v:{&l.integratedLufs,&l.integratedToleranceLu,&l.maxLoudnessRangeLu,&l.truePeakCeilingDbtp,&l.maxMomentaryLufs,&l.maxShortTermLufs})*v=r.f();
    l.validationNote=r.t();require(p.digest()==r.t(64)&&p.validationError().isEmpty(),"Invalid explanation target definition");return p;
}
void measured(Writer& w,const ProposalMeasurementRecord& m)
{
    for(const auto* s:{&m.identity,&m.content,&m.method,&m.build,&m.platform,&m.originSession,&m.identicalSource,&m.evaluation})w.t(*s);
    w.u(m.originJob);w.u(m.originRevision);w.u(int(m.domain));for(const auto& metric:m.metrics){w.u(int(metric.status));optional(w,metric.value);}
}
ProposalMeasurementRecord measured(Reader& r)
{
    ProposalMeasurementRecord m;for(auto* s:{&m.identity,&m.content,&m.method,&m.build,&m.platform,&m.originSession,&m.identicalSource,&m.evaluation})*s=r.t();
    m.originJob=r.u();m.originRevision=r.u();m.domain=MeasurementDomain(r.count(3));
    for(auto& metric:m.metrics){metric.status=measurement::Status(r.count(3));metric.value=optional(r);require((metric.status==measurement::Status::available)==bool(metric.value),"Inconsistent historical metric status");}
    require(hash(m.identity)&&hash(m.content),"Invalid historical measurement identity");return m;
}
void validate(const AppliedProposalExplanation& h)
{
    require(h.roles.confirmed&&h.roles.isValid()&&h.sampleRate>0&&h.frames>0,"Invalid historical programme format or roles");
    require(h.programmeId.length()==32&&hash(h.inputIdentity)&&hash(h.binding),"Invalid explanation programme binding");
    require(h.beforeIdentity==proposalPlanIdentity(h.before)&&h.afterIdentity==proposalPlanIdentity(h.after)&&h.beforeIdentity!=h.afterIdentity,"Invalid explanation before/after identities");
    require(h.preferences.validationError().isEmpty()&&h.requirements.validationError().isEmpty(),"Invalid historical preferences or requirements");
    require(h.source.domain==MeasurementDomain::source&&h.baseline.domain==MeasurementDomain::currentState&&(!h.candidate||h.candidate->domain==MeasurementDomain::privateCandidate),"Invalid historical measurement domains");
    require(h.source.content==h.baseline.content&&(!h.candidate||h.candidate->content==h.source.content),"Historical measurements belong to different source content");
    require(h.annotations.sections.empty()||(h.annotations.valid()&&h.annotations.inputIdentity==h.inputIdentity&&h.annotations.sampleRate==h.sampleRate&&h.annotations.frames==h.frames),"Historical annotations belong to different input");
    for(const auto* p:{&h.before,&h.after})if(p->sectionGains)require(p->sectionGains->inputIdentity==h.inputIdentity&&p->sectionGains->sampleRate==h.sampleRate&&p->sectionGains->frames==h.frames,"Historical plan belongs to different input");
    require(h.selection.corrections.size()<=10000&&h.passes>=0&&h.passes<=3,"Historical proposal exceeds supported bounds");
    int64_t last=-1;
    for(const auto& row:h.selection.corrections)
    {
        const auto a=std::lower_bound(h.annotations.sections.begin(),h.annotations.sections.end(),row.begin,[](const auto& a,int64_t b){return a.begin<b;});
        require(row.begin>last&&a!=h.annotations.sections.end()&&a->id==row.id&&a->begin==row.begin&&a->end==row.end,"Historical correction has no original annotation");last=row.begin;
        require(row.channel>=0&&row.channel<maxChannels&&row.originalDb>=-12&&row.originalDb<=12&&row.resultingDb>=-12&&row.resultingDb<=12,"Invalid historical correction gain/channel");
    }
}
}
std::vector<uint8_t> encodeAppliedExplanation(const AppliedProposalExplanation& h)
{
    validate(h);Writer w;w.u(1);
    for(const auto* s:{&h.programmeId,&h.inputIdentity,&h.binding,&h.method,&h.narrative,&h.beforeIdentity,&h.afterIdentity})w.t(*s);
    plan(w,h.before);plan(w,h.after);w.u(h.sampleRate);w.u(h.frames);w.u(h.roles.count);for(int n=0;n<h.roles.count;++n)w.u(int(h.roles.roles[size_t(n)]));target(w,h.target);
    w.t(juce::JSON::toString(h.requirements.encode(),true));w.t(h.requirements.digest());w.t(juce::JSON::toString(h.preferences.encode(),true));w.t(h.preferences.digest());w.u(int(h.applicability));
    annotations(w,h.annotations);
    const auto& s=h.selection;w.u(s.available);w.u(s.processingRegions);w.u(s.affectedFrames);w.t(s.reason);
    w.u(!s.plan.intervals.empty());if(!s.plan.intervals.empty()){ProcessingPlan p;p.schemaVersion=2;p.engineVersion=processing::SectionGainPlan::version;p.sectionGains=s.plan;plan(w,p);}
    w.u(s.corrections.size());for(const auto& c:s.corrections)
    {
        w.u(c.id);w.u(c.begin);w.u(c.end);w.u(c.peakEighthFrame);w.u(c.channel);w.u(c.padding);
        for(double v:{c.originalDb,c.resultingDb,c.peakDbtp,c.deficitDb,c.requestedDb,c.reductionDb,c.reserveDb,c.uncoveredDb})w.f(v);
    }
    const auto& g=h.global;w.u(int(g.status));w.t(g.identity);w.t(g.explanation);w.f(g.originalGainDb);
    for(auto v:{g.totalGainDb,g.deltaDb,g.predictedIntegrated,g.predictedTruePeak})optional(w,v);
    measured(w,h.source);measured(w,h.baseline);w.u(bool(h.candidate));if(h.candidate)measured(w,*h.candidate);w.u(h.passes);
    return std::move(w.bytes);
}
std::shared_ptr<const AppliedProposalExplanation> decodeAppliedExplanation(const std::vector<uint8_t>& bytes,JobControl& control)
{
    require(bytes.size()<=AppliedProposalExplanation::maximumBytes,"Applied explanation exceeds 8 MiB");Reader r{bytes,control};require(r.u()==1,"Unsupported applied explanation version");auto h=std::make_shared<AppliedProposalExplanation>();
    for(auto* s:{&h->programmeId,&h->inputIdentity,&h->binding,&h->method,&h->narrative,&h->beforeIdentity,&h->afterIdentity})*s=r.t();
    h->before=plan(r);h->after=plan(r);h->sampleRate=r.count(192000);h->frames=r.i();h->roles.count=r.count(maxChannels);h->roles.confirmed=true;
    for(int n=0;n<h->roles.count;++n)h->roles.roles[size_t(n)]=ChannelRole(r.count(2048));h->target=target(r);
    h->requirements=proposal::JobRequirements::fromJson(r.t());require(h->requirements.digest()==r.t(64),"Historical requirements digest mismatch");
    h->preferences=proposal::ProposalPolicy::decode(juce::JSON::parse(r.t()));require(h->preferences.digest()==r.t(64),"Historical preferences digest mismatch");h->applicability=proposal::Applicability(r.count(5));
    h->annotations=annotations(r);auto& s=h->selection;s.policy=h->preferences;s.available=r.b();s.processingRegions=r.count(20000);s.affectedFrames=r.i();s.reason=r.t();
    if(r.b()){auto p=plan(r);require(bool(p.sectionGains),"Missing historical selection plan");s.plan=*p.sectionGains;}
    const int n=r.count(10000);for(int i=0;i<n;++i)
    {
        proposal::SectionCorrection c;c.id=r.u();c.begin=r.i();c.end=r.i();c.peakEighthFrame=r.i();c.channel=r.count(maxChannels-1);c.padding=r.b();
        for(auto* v:{&c.originalDb,&c.resultingDb,&c.peakDbtp,&c.deficitDb,&c.requestedDb,&c.reductionDb,&c.reserveDb,&c.uncoveredDb})*v=r.f();s.corrections.push_back(c);
    }
    auto& g=h->global;g.status=proposal::TrimStatus(r.count(3));g.identity=r.t(64);g.explanation=r.t();g.originalGainDb=r.f();
    for(auto* v:{&g.totalGainDb,&g.deltaDb,&g.predictedIntegrated,&g.predictedTruePeak})*v=optional(r);
    h->source=measured(r);h->baseline=measured(r);if(r.b())h->candidate=measured(r);h->passes=r.count(3);
    require(r.position==bytes.size(),"Trailing applied explanation data");validate(*h);return h;
}
}

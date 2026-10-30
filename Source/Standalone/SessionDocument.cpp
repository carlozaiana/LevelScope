#include "SessionDocument.h"
#include <juce_cryptography/juce_cryptography.h>
#include <charconv>
#include <cstring>
#include <set>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cerrno>

namespace levelscope::standalone::sessionfile
{
namespace
{
using Bytes = std::vector<uint8_t>;
using V = juce::var;
constexpr uint64_t coreLimit = 4 * 1024 * 1024;
struct Failure { ErrorCode code; juce::String message; };
void require (bool ok, const juce::String& why) { if (!ok) throw Failure { ErrorCode::invalidMetadata, why }; }
void check (JobControl& c) { if (c.isCancelled()) throw Failure { ErrorCode::cancelled, "Session operation cancelled; previous session preserved" }; }
juce::String digest (const void* p, size_t n) { return juce::SHA256 (p, n).toHexString(); }
juce::String digest (const Bytes& b) { return digest (b.data(), b.size()); }
bool hashValid (const juce::String& h) { return h.length() == 64 && h.containsOnly ("0123456789abcdef"); }
V object() { return V (new juce::DynamicObject()); }
void put (V& v, const char* k, V x) { v.getDynamicObject()->setProperty (k, std::move(x)); }
V get (const V& v, const char* k) { require (v.isObject() && v.hasProperty(k), "Missing session field: " + juce::String(k)); return v[k]; }
juce::String string (const V& v, int max = 32768) { require(v.isString() && v.toString().length() <= max, "Invalid session text"); return v.toString(); }
juce::String str (const V& v, const char* k, int max = 32768) { return string(get(v,k),max); }
bool boolean (const V& v) { require(v.isBool(), "Invalid session boolean"); return static_cast<bool>(v); }
int number (const V& v, int lo, int hi) { require(v.isInt(), "Invalid session integer"); int n=v; require(n>=lo&&n<=hi,"Session integer outside supported range"); return n; }
template<class T> T integer (const V& v)
{
    auto s=string(v,32).toStdString(); T n{}; auto r=std::from_chars(s.data(),s.data()+s.size(),n);
    require(!s.empty()&&r.ec==std::errc()&&r.ptr==s.data()+s.size(),"Invalid exact session integer"); return n;
}
V signedValue (int64_t n) { return juce::String(n); }
V unsignedValue (uint64_t n) { return juce::String(n); }
V realValue (double n)
{
    require(std::isfinite(n),"Non-finite session number"); uint64_t bits; std::memcpy(&bits,&n,8);
    return juce::String::toHexString(static_cast<juce::int64>(bits)).paddedLeft('0',16);
}
double real (const V& v)
{
    auto s=string(v,16); require(s.length()==16&&s.containsOnly("0123456789abcdef"),"Invalid exact binary64 number");
    uint64_t bits=static_cast<uint64_t>(s.getHexValue64()); double n; std::memcpy(&n,&bits,8);
    require(std::isfinite(n),"Non-finite session number"); return n;
}
const juce::Array<V>& array (const V& v, int max)
{ require(v.isArray()&&v.size()<=max,"Invalid or oversized session array"); return *v.getArray(); }
// Bound JSON recursion and require our canonical encoding. Duplicate object keys,
// trailing data and alternate ambiguous encodings cannot survive this equality.
V parse (const std::string& bytes)
{
    require(bytes.size()<=coreLimit,"Session core exceeds 4 MiB");
    int depth=0; bool quoted=false,escaped=false;
    for (char c:bytes)
    {
        if(quoted) { if(escaped) escaped=false; else if(c=='\\') escaped=true; else if(c=='"') quoted=false; }
        else if(c=='"') quoted=true;
        else if(c=='{'||c=='[') { require(++depth<=24,"Session nesting too deep"); }
        else if(c=='}'||c==']') require(--depth>=0,"Invalid session nesting");
    }
    require(!quoted&&depth==0,"Truncated session core");
    V v; auto result=juce::JSON::parse(juce::String::fromUTF8(bytes.data(),static_cast<int>(bytes.size())),v);
    require(result.wasOk()&&v.isObject(),"Invalid session JSON");
    require(juce::JSON::toString(v,true).toStdString()==bytes,"Non-canonical or duplicate-key session core"); return v;
}
void u64 (Bytes& b,uint64_t n) { for(int i=0;i<8;++i)b.push_back(static_cast<uint8_t>(n>>(8*i))); }
void f64 (Bytes& b,double n) { require(std::isfinite(n),"Non-finite cache value"); uint64_t x;std::memcpy(&x,&n,8);u64(b,x); }
void text (Bytes& b,const juce::String& s) { auto x=s.toStdString();u64(b,x.size());b.insert(b.end(),x.begin(),x.end()); }
struct Reader
{
    const Bytes& b; size_t p=0; JobControl& c;
    uint64_t u() { check(c);require(p<=b.size()&&b.size()-p>=8,"Truncated measurement cache");uint64_t n=0;for(int i=0;i<8;++i)n|=uint64_t(b[p++])<<(8*i);return n; }
    int64_t i() { auto n=u();int64_t s;std::memcpy(&s,&n,8);return s; }
    double f() { auto n=u();double d;std::memcpy(&d,&n,8);require(std::isfinite(d),"Non-finite cache metric");return d; }
    juce::String t() { auto n=u();require(n<=coreLimit&&p<=b.size()&&n<=b.size()-p,"Invalid cache text");auto s=juce::String::fromUTF8(reinterpret_cast<const char*>(b.data()+p),static_cast<int>(n));p+=n;return s; }
};
V layoutValue (const ChannelLayout& l)
{
    auto v=object();put(v,"count",l.count);put(v,"confirmed",l.confirmed);put(v,"manual",l.manuallyInterpreted);
    juce::Array<V> roles;for(int n=0;n<l.count;++n)roles.add(static_cast<int>(l.roles[n]));put(v,"roles",roles);return v;
}
ChannelLayout layout (const V& v)
{
    ChannelLayout l;l.count=number(get(v,"count"),1,maxChannels);l.confirmed=boolean(get(v,"confirmed"));l.manuallyInterpreted=boolean(get(v,"manual"));
    auto a=array(get(v,"roles"),maxChannels);require(a.size()==l.count,"Channel count/role mismatch");
    for(int n=0;n<l.count;++n) { int r=number(a[n],0,2048);l.roles[n]=static_cast<ChannelRole>(r);
        require(std::find(supportedRoles.begin(),supportedRoles.end(),l.roles[n])!=supportedRoles.end() || (!l.confirmed&&(r==0||r>=static_cast<int>(juce::AudioChannelSet::discreteChannel0))),"Unsupported role token"); }
    require(!l.confirmed||l.isValid(),"Invalid confirmed channel interpretation");return l;
}
V mediaValue (const MediaInfo& m)
{
    auto v=object();put(v,"rate",realValue(m.sampleRate));put(v,"frames",signedValue(m.frames));put(v,"bytes",signedValue(m.bytes));put(v,"channels",m.channels);put(v,"bits",m.bitsPerSample);
    put(v,"modified",signedValue(m.modified.toMilliseconds()));put(v,"format",m.format);put(v,"sha256",m.sha256);put(v,"provenance",m.layoutProvenance);put(v,"restriction",m.measurementRestriction);return v;
}
MediaInfo media (const V& v)
{
    MediaInfo m;m.sampleRate=real(get(v,"rate"));m.frames=integer<int64_t>(get(v,"frames"));m.bytes=integer<int64_t>(get(v,"bytes"));m.channels=number(get(v,"channels"),1,maxChannels);m.bitsPerSample=number(get(v,"bits"),1,32);
    m.modified=juce::Time(integer<int64_t>(get(v,"modified")));m.format=str(v,"format",128);m.sha256=str(v,"sha256",64);m.layoutProvenance=str(v,"provenance");m.measurementRestriction=str(v,"restriction");
    require(m.sampleRate>0&&m.sampleRate<=768000&&m.sampleRate==std::floor(m.sampleRate)&&m.frames>0&&m.frames<=INT64_MAX/8-64&&m.bytes>0&&hashValid(m.sha256),"Invalid source format/identity");return m;
}
V assemblyValue (const ProgrammeAssembly& group)
{
    auto v=object();put(v,"schema",ProgrammeAssembly::version);put(v,"confirmation",group.alignmentConfirmation);
    juce::Array<V> rows;
    for(const auto& m:group.members){auto row=object();put(row,"id",m.id);put(row,"locator",m.file.getFullPathName());put(row,"media",mediaValue(m.media));put(row,"role",static_cast<int>(m.role));put(row,"offset",signedValue(m.offsetFrames));rows.add(row);}
    put(v,"members",rows);return v;
}
ProgrammeAssembly assembly (const V& v)
{
    ProgrammeAssembly group;number(get(v,"schema"),2,2);group.alignmentConfirmation=str(v,"confirmation",64);
    for(const auto& row:array(get(v,"members"),6))
    {
        ProgrammeMember m;m.id=str(row,"id",128);auto path=str(row,"locator");require(juce::File::isAbsolutePath(path),"Invalid member locator");m.file=juce::File(path);m.media=media(get(row,"media"));m.role=static_cast<ChannelRole>(number(get(row,"role"),0,2048));m.offsetFrames=integer<int64_t>(get(row,"offset"));group.members.push_back(std::move(m));
    }
    require(group.validate().isEmpty(),group.validate());return group;
}
V componentValue (const ProgrammeComponent& c)
{
    auto v=object();put(v,"id",c.id);put(v,"purpose",c.purpose);put(v,"language",c.language);put(v,"versionLabel",c.versionLabel);
    put(v,"locator",c.file.getFullPathName());put(v,"media",mediaValue(c.media));put(v,"layout",layoutValue(c.layout));
    put(v,"offset",signedValue(c.offsetFrames));put(v,"assembly",c.assembly?assemblyValue(*c.assembly):V());return v;
}
ProgrammeComponent component (const V& v)
{
    ProgrammeComponent c;c.id=str(v,"id",128);c.purpose=str(v,"purpose",128);c.language=str(v,"language",128);c.versionLabel=str(v,"versionLabel",256);
    const auto path=str(v,"locator");require(path.isEmpty()||juce::File::isAbsolutePath(path),"Invalid component locator");
    c.file=path.isEmpty()?juce::File():juce::File(path);c.media=media(get(v,"media"));c.layout=layout(get(v,"layout"));c.offsetFrames=integer<int64_t>(get(v,"offset"));
    if(!get(v,"assembly").isVoid())c.assembly=assembly(get(v,"assembly"));
    require(c.assembly||!path.isEmpty(),"Missing component file locator");require(c.validate().isEmpty(),c.validate());return c;
}
V compositionValue (const ProgrammeComposition& c)
{
    auto v=object();put(v,"schema",ProgrammeComposition::schemaVersion);put(v,"engine",ProgrammeComposition::engineVersion);
    put(v,"preset",static_cast<int>(c.preset));put(v,"me",componentValue(c.me));put(v,"dialogue",componentValue(c.dialogue));put(v,"confirmation",c.confirmation);return v;
}
ProgrammeComposition composition (const V& v)
{
    ProgrammeComposition c;number(get(v,"schema"),1,1);require(str(v,"engine")==ProgrammeComposition::engineVersion,"Unsupported composition engine");
    c.preset=static_cast<CompositionPreset>(number(get(v,"preset"),0,2));c.me=component(get(v,"me"));c.dialogue=component(get(v,"dialogue"));c.confirmation=str(v,"confirmation",64);
    require(c.validate().isEmpty(),c.validate());return c;
}
V sectionsValue (const proposal::SectionMap& s)
{
    if(s.sections.empty())return V();
    require(s.valid(),"Invalid section annotations");auto v=object();put(v,"schema",s.schemaVersion);put(v,"input",s.inputIdentity);put(v,"rate",s.sampleRate);put(v,"frames",signedValue(s.frames));put(v,"revision",unsignedValue(s.revision));put(v,"nextId",unsignedValue(s.nextId));
    juce::Array<V> rows;for(auto& r:s.sections){juce::Array<V> a;a.add(unsignedValue(r.id));a.add(signedValue(r.begin));a.add(signedValue(r.end));a.add(r.name);rows.add(a);}put(v,"rows",rows);return v;
}
proposal::SectionMap sections (const V& v)
{
    proposal::SectionMap s;if(v.isVoid())return s;s.schemaVersion=number(get(v,"schema"),1,1);s.inputIdentity=str(v,"input",64);s.sampleRate=number(get(v,"rate"),1,192000);s.frames=integer<int64_t>(get(v,"frames"));s.revision=integer<uint64_t>(get(v,"revision"));s.nextId=integer<uint64_t>(get(v,"nextId"));
    for(auto& row:array(get(v,"rows"),10000)){auto a=array(row,4);require(a.size()==4,"Invalid section row");s.sections.push_back({integer<uint64_t>(a[0]),integer<int64_t>(a[1]),integer<int64_t>(a[2]),string(a[3],128)});}
    require(s.valid(),"Invalid section IDs, names or coverage");return s;
}
V gainPlanValue(const ProcessingPlan& p)
{
    auto v=object();put(v,"gain",realValue(p.staticGainDb));put(v,"planSchema",p.schemaVersion);put(v,"engine",p.engineVersion);
    if(p.sectionGains)
    {
        const auto& g=*p.sectionGains;auto a=object();put(a,"schema",g.schema);put(a,"rate",g.sampleRate);put(a,"frames",signedValue(g.frames));put(a,"input",g.inputIdentity);
        juce::Array<V> rows;for(const auto& i:g.intervals){juce::Array<V> row;row.add(signedValue(i.begin));row.add(signedValue(i.end));row.add(realValue(i.db));rows.add(row);}put(a,"intervals",rows);put(v,"sectionGains",a);
    }
    else put(v,"sectionGains",V());
    return v;
}
ProcessingPlan gainPlan(const V& v,bool modern)
{
    ProcessingPlan p;p.staticGainDb=real(get(v,"gain"));p.schemaVersion=number(get(v,"planSchema"),1,modern?2:1);p.engineVersion=str(v,"engine",128);
    if(modern&&!get(v,"sectionGains").isVoid())
    {
        const auto a=get(v,"sectionGains");processing::SectionGainPlan g;g.schema=number(get(a,"schema"),1,1);g.sampleRate=number(get(a,"rate"),1,192000);g.frames=integer<int64_t>(get(a,"frames"));g.inputIdentity=str(a,"input",64);
        for(const auto& row:array(get(a,"intervals"),10000)){const auto& r=array(row,3);require(r.size()==3,"Invalid section-gain row");g.intervals.push_back({integer<int64_t>(r[0]),integer<int64_t>(r[1]),real(r[2])});}
        require(g.validationError().isEmpty()&&!g.isUnity(),"Invalid applied section-gain plan");p.sectionGains=g.canonical();
    }
    require(StaticGainProcessor::supports(p)||(p.schemaVersion==2&&p.engineVersion==processing::SectionGainPlan::version&&p.sectionGains),"Unknown processing version");return p;
}
ServiceContext context (const Document& d,MeasurementDomain domain)
{
    ServiceContext c;c.sessionId=d.runtimeId;c.sourceFile=d.source.file;c.media=d.source.media;c.assembly=d.source.assembly;c.composition=d.source.composition;c.fileLayout=d.source.layout;c.range={0,c.media.frames};
    c.processing=domain==MeasurementDomain::source?ProcessingPlan{}:d.current.processingPlan();return c;
}
void validateCore (const Document& d)
{
    if(d.markerImport){const auto receipt=proposal::markers::decode(proposal::markers::encode(*d.markerImport));require(d.sections.valid()&&receipt.imported.inputIdentity==d.sections.inputIdentity&&receipt.imported.sampleRate==d.sections.sampleRate&&receipt.imported.frames==d.sections.frames,"Marker provenance belongs to different programme");}

    const auto uuidValid = [](const juce::String& id) { return id.length()==32 && id.containsOnly("0123456789abcdef"); };
    require(uuidValid(d.id)&&uuidValid(d.programmeId),"Invalid persistent identity");
    if(d.source.assembly||d.source.composition)require(context(d,MeasurementDomain::source).isValid(),"Invalid confirmed programme definition");
    require(d.revision<UINT64_MAX,"Exhausted document revision");
    require(std::isfinite(d.current.staticGainDb)&&d.current.staticGainDb>=-60&&d.current.staticGainDb<=24,"Unsupported processing gain");
    require(!d.current.isInitialized||(d.source.hasSource&&d.source.layout.confirmed&&d.source.layout.isValid()),"Current State requires confirmed source roles");
    require(d.current.isInitialized||(d.current.staticGainDb==0&&!d.current.sectionGains),"Uninitialized processing plan");
    if(d.current.sectionGains)
    {
        const auto c=context(d,MeasurementDomain::currentState);require(ProgrammeGainProcessor::validationError(c).isEmpty(),"Applied section plan does not match this programme");
        require(d.current.sectionGainAnnotations.valid()&&d.current.sectionGainAnnotations.inputIdentity==d.current.sectionGains->inputIdentity
            &&d.current.sectionGainAnnotations.frames==d.current.sectionGains->frames&&d.current.sectionGainAnnotations.sampleRate==d.current.sectionGains->sampleRate,"Missing or mismatched applied annotation provenance");
        for(const auto& interval:d.current.sectionGains->intervals)
            require(std::any_of(d.current.sectionGainAnnotations.sections.begin(),d.current.sectionGainAnnotations.sections.end(),
                [&](const auto& annotation){return annotation.end==interval.end;}),"Applied gain boundary is absent from its annotation provenance");
    }
    else require(d.current.sectionGainAnnotations.sections.empty(),"Provenance without applied section gains");
    require(d.audition==MeasurementDomain::source||d.audition==MeasurementDomain::currentState,"Unsupported audition selection");
    if(!d.sections.sections.empty())require(d.sections.valid()&&d.sections.inputIdentity==sectionInputIdentity(context(d,MeasurementDomain::source)),"Sections do not describe this programme");
    require(d.profile.id.isNotEmpty()&&hashValid(d.profile.digest),"Invalid target reference");
    require(d.reference.file==juce::File()||hashValid(d.reference.sha256),"Reference video identity unavailable; wait for loading or clear the reference");
}
void optionalValue (Bytes& b,std::optional<double> d) {u64(b,d.has_value());if(d)f64(b,*d);}
std::optional<double> optionalValue (Reader& r) {auto present=r.u();require(present<=1,"Invalid optional metric");return present?std::optional<double>(r.f()):std::nullopt;}
void metric (Bytes& b,const measurement::Metric& m){u64(b,static_cast<uint64_t>(m.status));optionalValue(b,m.value);}
measurement::Metric metric (Reader& r)
{
    auto s=r.u();require(s<=3,"Invalid metric status");measurement::Metric m{optionalValue(r),static_cast<measurement::Status>(s)};
    require((m.status==measurement::Status::available)==m.value.has_value(),"Metric availability mismatch");return m;
}
Bytes cacheBytes (const AnalysisSnapshot& a,JobControl& control)
{
    Bytes b;const auto& c=a.context;auto v=gainPlanValue(c.processing);
    put(v,"schema",4);put(v,"assembly",c.assembly?assemblyValue(*c.assembly):V());put(v,"composition",c.composition?compositionValue(*c.composition):V());put(v,"domain",static_cast<int>(a.domain));put(v,"media",mediaValue(c.media));put(v,"layout",layoutValue(c.fileLayout));put(v,"begin",signedValue(c.range.begin));put(v,"end",signedValue(c.range.end));
    put(v,"gain",realValue(c.processing.staticGainDb));put(v,"planSchema",c.processing.schemaVersion);put(v,"engine",c.processing.engineVersion);put(v,"build",c.buildVersion);put(v,"platform",c.platformAndFloatingPointMode);put(v,"method",a.measurementVersion);put(v,"identity",a.identity);put(v,"content",a.measuredContentHash);
    put(v,"session",a.originSessionId.isEmpty()?c.sessionId:a.originSessionId);put(v,"job",unsignedValue(a.originSessionId.isEmpty()?c.ticket.id:a.originJobId));put(v,"revision",unsignedValue(a.originSessionId.isEmpty()?c.ticket.sessionRevision:a.originSessionRevision));
    text(b,juce::JSON::toString(v,true));const auto& d=a.detail;u64(b,d.sampleRate);u64(b,d.frames);
    for(auto* m:{&d.integrated,&d.lra,&d.momentaryMax,&d.shortTermMax,&d.samplePeak,&d.truePeak})metric(b,*m);
    for(auto x:{d.integratedRelativeGate,d.lraRelativeGate,d.p10,d.p95})optionalValue(b,x);
    u64(b,d.integratedIncluded);u64(b,d.lraIncluded);u64(b,d.lraStable);
    for(auto* w:{&d.momentary,&d.shortTerm,&d.lraWindows}){u64(b,w->size());for(auto& x:*w){check(control);u64(b,x.endFrame);f64(b,x.energy);}}
    for(auto* a1:{&d.channelSamplePeaks,&d.channelTruePeaks}){u64(b,a1->size());for(auto x:*a1)f64(b,x);}
    for(auto* a1:{&d.samplePeakFrames,&d.truePeakEighthFrames}){u64(b,a1->size());for(auto x:*a1)u64(b,x);}
    require(b.size()<=maximumBytes/2,"Measurement cache too large");return b;
}
std::shared_ptr<const AnalysisSnapshot> cache (const Bytes& b,MeasurementDomain domain,JobControl& control)
{
    Reader r{b,0,control};auto v=parse(r.t().toStdString());const int cacheSchema=number(get(v,"schema"),1,4);
    auto a=std::make_shared<AnalysisSnapshot>();a->domain=static_cast<MeasurementDomain>(number(get(v,"domain"),0,1));require(a->domain==domain,"Cache domain mismatch");auto& c=a->context;
    if(cacheSchema>=2&&!get(v,"assembly").isVoid())c.assembly=assembly(get(v,"assembly"));
    if(cacheSchema>=3&&!get(v,"composition").isVoid())c.composition=composition(get(v,"composition"));
    c.media=media(get(v,"media"));c.fileLayout=layout(get(v,"layout"));c.range={integer<int64_t>(get(v,"begin")),integer<int64_t>(get(v,"end"))};c.processing=gainPlan(v,cacheSchema>=4);c.buildVersion=str(v,"build",128);c.platformAndFloatingPointMode=str(v,"platform",256);a->measurementVersion=str(v,"method",128);a->identity=str(v,"identity",64);a->measuredContentHash=str(v,"content",64);
    c.sessionId=str(v,"session",128);c.ticket={integer<uint64_t>(get(v,"job")),integer<uint64_t>(get(v,"revision"))};a->originSessionId=c.sessionId;a->originJobId=c.ticket.id;a->originSessionRevision=c.ticket.sessionRevision;a->restoredFromSession=true;
    auto& d=a->detail;auto rate=r.u();require(rate<=192000,"Invalid cached rate");d.sampleRate=static_cast<int>(rate);d.frames=r.i();require(d.frames>0&&d.frames==c.media.frames&&d.sampleRate==c.media.sampleRate,"Cached duration mismatch");
    for(auto* m:{&d.integrated,&d.lra,&d.momentaryMax,&d.shortTermMax,&d.samplePeak,&d.truePeak})*m=metric(r);
    for(auto* x:{&d.integratedRelativeGate,&d.lraRelativeGate,&d.p10,&d.p95})*x=optionalValue(r);
    d.integratedIncluded=r.u();d.lraIncluded=r.u();auto stable=r.u();require(stable<=1,"Invalid LRA stability");d.lraStable=stable;
    for(auto* w:{&d.momentary,&d.shortTerm,&d.lraWindows})
    {
        auto n=r.u();require(n<=1000000&&n<=(b.size()-r.p)/16,"Oversized cached window vector");w->reserve(static_cast<size_t>(n));int64_t last=0;
        for(uint64_t i=0;i<n;++i){auto end=r.i();auto energy=r.f();require(end>last&&end<=d.frames+int64_t(d.sampleRate)*3&&energy>=0,"Invalid cached window");w->push_back({end,energy});last=end;}
    }
    for(auto* a1:{&d.channelSamplePeaks,&d.channelTruePeaks}){auto n=r.u();require(n==static_cast<uint64_t>(c.media.channels),"Cached peak channel mismatch");for(uint64_t i=0;i<n;++i){auto x=r.f();require(x>=0,"Negative peak");a1->push_back(x);}}
    for(auto* a1:{&d.samplePeakFrames,&d.truePeakEighthFrames}){auto n=r.u();require(n==static_cast<uint64_t>(c.media.channels),"Cached peak position count mismatch");for(uint64_t i=0;i<n;++i){auto x=r.i();require(x>=-512&&x<=(d.frames+64)*8,"Invalid cached peak position");a1->push_back(x);}}
    require(r.p==b.size()&&d.integratedIncluded<=d.momentary.size()&&d.lraIncluded<=d.lraWindows.size(),"Malformed measurement cache");
    a->measured={d.integrated.value,d.lra.value,d.truePeak.value,d.samplePeak.value};
    require(a->measurementVersion==measurement::version&&hashValid(a->identity)&&a->identity==OfflineAnalysisService::identityFor(c,domain,a->measuredContentHash),"Incompatible measurement identity");
    require(measuredProgrammeForSnapshot(a.get(),c,domain,proposal::Applicability::unconfirmed).currentIdentityMatches,"Invalid measured cache binding");return a;
}
V coreValue (const Document& d,const std::array<Bytes,2>& caches)
{
    auto v=object();put(v,"schema",schemaVersion);put(v,"documentId",d.id);put(v,"revision",unsignedValue(d.revision));
    auto p=object();put(p,"schema",programmeAssemblyVersion);put(p,"id",d.programmeId);put(p,"mode",d.source.composition?"composition":d.source.assembly?"mono_files":"interleaved");put(p,"purpose","programme_audio");put(p,"alignment","0/1");
    juce::Array<V> files,channels;
    if(d.source.hasSource)
    {
        if(d.source.composition)put(p,"composition",compositionValue(*d.source.composition));
        else if(d.source.assembly)
        {
            put(p,"assembly",assemblyValue(*d.source.assembly));
            for(const auto& m:d.source.assembly->members)
            {
                auto f=object();put(f,"id",m.id);put(f,"locator",m.file.getFullPathName());put(f,"media",mediaValue(m.media));files.add(f);
                auto a=object();put(a,"fileId",m.id);put(a,"fileChannel",0);put(a,"role",static_cast<int>(m.role));channels.add(a);
            }
        }
        else {
        auto f=object();put(f,"id","audio-1");put(f,"locator",d.source.file.getFullPathName());put(f,"media",mediaValue(d.source.media));files.add(f);
        for(int n=0;n<d.source.layout.count;++n){auto a=object();put(a,"fileId","audio-1");put(a,"fileChannel",n);put(a,"role",static_cast<int>(d.source.layout.roles[n]));channels.add(a);}
        }
        put(p,"layout",layoutValue(d.source.layout));put(p,"rangeBegin","0");put(p,"rangeEnd",signedValue(d.source.media.frames));
    }
    put(p,"files",files);put(p,"channels",channels);put(v,"programme",p);
    auto cur=gainPlanValue(d.current.processingPlan());put(cur,"initialized",d.current.isInitialized);put(cur,"annotationProvenance",sectionsValue(d.current.sectionGainAnnotations));
    // Explicit extension boundary: future component plans are Current processing,
    // never mutations of original Source. M11 only accepts the empty identity set.
    put(cur,"componentPlans",juce::Array<V>());put(cur,"wholeProgrammeGraph","");put(v,"current",cur);
    auto pr=object();put(pr,"id",d.profile.id);put(pr,"version",d.profile.version);put(pr,"digest",d.profile.digest);put(pr,"applicability",static_cast<int>(d.applicability));put(v,"profile",pr);
    put(v,"jobRequirements",d.jobRequirements.encode());
    put(v,"proposalPolicy",d.proposalPolicy.encode());put(v,"lastAppliedProposalPolicy",d.lastAppliedProposalPolicy?d.lastAppliedProposalPolicy->encode():V());put(v,"lastAppliedProposalBinding",d.lastAppliedProposalBinding);
    put(v,"sections",sectionsValue(d.sections));put(v,"markerImport",d.markerImport?proposal::markers::encode(*d.markerImport):V());put(v,"audition",static_cast<int>(d.audition));
    auto video=object();put(video,"locator",d.reference.file.getFullPathName());put(video,"sha256",d.reference.sha256);put(video,"stream",d.reference.streamId);put(video,"offset",juce::String(d.reference.startOnAudio.text()));put(v,"video",video);
    juce::Array<V> records;uint64_t offset=0;for(int n=0;n<2;++n){auto a=object();put(a,"domain",n);put(a,"offset",unsignedValue(offset));put(a,"size",unsignedValue(caches[n].size()));put(a,"sha256",digest(caches[n]));records.add(a);offset+=caches[n].size();}put(v,"caches",records);return v;
}
Document core (const V& v)
{
    Document d;d.loadedSchema=number(get(v,"schema"),1,7);d.id=str(v,"documentId",36);d.revision=integer<uint64_t>(get(v,"revision"));auto p=get(v,"programme");require(number(get(p,"schema"),1,7)==d.loadedSchema,"Programme/core schema mismatch");d.programmeId=str(p,"id",36);
    const auto mode=str(p,"mode");require((mode=="interleaved"||(mode=="mono_files"&&d.loadedSchema>=2)||(mode=="composition"&&d.loadedSchema>=3))&&str(p,"purpose")=="programme_audio"&&str(p,"alignment")=="0/1","Unsupported programme assembly");
    auto files=array(get(p,"files"),mode=="interleaved"?1:6);auto channels=array(get(p,"channels"),maxChannels);
    if(mode=="composition")
    {
        require(files.isEmpty()&&channels.isEmpty(),"Composition members/mappings belong to its explicit recipe");
        d.source.setValidatedComposition(composition(get(p,"composition")));
        const auto l=layout(get(p,"layout"));require(l.confirmed&&l.count==d.source.layout.count,"Composition layout mismatch");
        for(int n=0;n<l.count;++n)require(l.roles[n]==d.source.layout.roles[n],"Composition output role mismatch");
        require(integer<int64_t>(get(p,"rangeBegin"))==0&&integer<int64_t>(get(p,"rangeEnd"))==d.source.media.frames,"Only complete compositions supported");
    }
    else if(!files.isEmpty())
    {
        if(mode=="mono_files")
        {
            const auto group=assembly(get(p,"assembly"));d.source.setValidatedProgramme(group);
            require(files.size()==static_cast<int>(group.members.size())&&channels.size()==files.size(),"Incomplete programme members/assignments");
            for(int n=0;n<files.size();++n)
            {
                const auto& m=group.members[static_cast<size_t>(n)];
                require(str(files[n],"id")==m.id&&str(files[n],"locator")==m.file.getFullPathName()
                    &&juce::JSON::toString(get(files[n],"media"),true)==juce::JSON::toString(mediaValue(m.media),true)
                    &&str(channels[n],"fileId")==m.id&&number(get(channels[n],"fileChannel"),0,0)==0
                    &&number(get(channels[n],"role"),0,2048)==static_cast<int>(m.role),"Programme member/assignment mismatch");
            }
            const auto l=layout(get(p,"layout"));require(l.confirmed&&l.count==d.source.layout.count,"Invalid group layout");
            for(int n=0;n<l.count;++n)require(l.roles[n]==d.source.layout.roles[n],"Group role mismatch");
        }
        else {
        auto f=files[0];require(str(f,"id")=="audio-1","Invalid source file ID");auto path=str(f,"locator");require(juce::File::isAbsolutePath(path),"Invalid source locator");d.source.file=juce::File(path);d.source.fullPath=path;d.source.displayName=d.source.file.getFileName();d.source.hasSource=true;d.source.media=media(get(f,"media"));d.source.fileSizeBytes=d.source.media.bytes;d.source.extension=d.source.file.getFileExtension().trimCharactersAtStart(".");d.source.layout=layout(get(p,"layout"));
        require(d.source.layout.count==d.source.media.channels&&channels.size()==d.source.media.channels,"Programme channel count mismatch");
        for(int n=0;n<channels.size();++n)require(str(channels[n],"fileId")=="audio-1"&&number(get(channels[n],"fileChannel"),0,maxChannels-1)==n&&number(get(channels[n],"role"),0,2048)==static_cast<int>(d.source.layout.roles[n]),"Invalid ordered file channel assignment");
        }
        require(integer<int64_t>(get(p,"rangeBegin"))==0&&integer<int64_t>(get(p,"rangeEnd"))==d.source.media.frames,"Only whole-programme sessions are supported");
    }
    else require(channels.isEmpty()&&mode=="interleaved","Channels reference a missing file record");
    if(d.loadedSchema<6)d.notices.add("Older session migrated in memory to default proposal preferences; media, annotations and exact applied processing preserved. Saving writes schema 7 and requires an M17 reader; use Save As to retain the original older-readable copy. Old transient proposals are not restored.");
    else
    {
        d.proposalPolicy=proposal::ProposalPolicy::decode(v["proposalPolicy"]);
        if(d.proposalPolicy.validationError().isNotEmpty())d.notices.add(d.proposalPolicy.validationError());
        if(!get(v,"lastAppliedProposalPolicy").isVoid())d.lastAppliedProposalPolicy=proposal::ProposalPolicy::decode(get(v,"lastAppliedProposalPolicy"));
        d.lastAppliedProposalBinding=str(v,"lastAppliedProposalBinding",64);
    }
    if(d.loadedSchema>=7)
    {
        try {d.jobRequirements=proposal::JobRequirements::decode(v["jobRequirements"]);}
        catch(const std::exception& e){d.jobRequirements=proposal::JobRequirements::unavailable(e.what());d.notices.add("Job requirements unavailable; core settings/media recovered. Explicit replacement required; no dependent approval is valid.");}
    }
    else d.notices.add("No added LRA requirement in this older session. Saving writes schema 7 and requires an M17 reader; use Save As to retain the original.");
    auto cur=get(v,"current");const auto plan=gainPlan(cur,d.loadedSchema>=4);d.current.isInitialized=boolean(get(cur,"initialized"));d.current.staticGainDb=plan.staticGainDb;d.current.sectionGains=plan.sectionGains;d.current.needsRemeasurement=d.current.isInitialized;
    if(d.loadedSchema>=4)d.current.sectionGainAnnotations=sections(get(cur,"annotationProvenance"));
    if(d.loadedSchema>=3)require(array(get(cur,"componentPlans"),0).isEmpty()&&str(cur,"wholeProgrammeGraph").isEmpty(),"Component DSP and arbitrary whole-programme graphs are not implemented in M12");
    auto pr=get(v,"profile");d.profile={str(pr,"id",256),str(pr,"version",256),str(pr,"digest",64)};d.applicability=static_cast<proposal::Applicability>(number(get(pr,"applicability"),0,5));d.sections=sections(get(v,"sections"));d.audition=static_cast<MeasurementDomain>(number(get(v,"audition"),0,1));
    if(d.loadedSchema>=5&&!get(v,"markerImport").isVoid())d.markerImport=proposal::markers::decode(get(v,"markerImport"));
    auto vid=get(v,"video");auto path=str(vid,"locator");require(path.isEmpty()||juce::File::isAbsolutePath(path),"Invalid video locator");d.reference.file=path.isEmpty()?juce::File():juce::File(path);d.reference.sha256=str(vid,"sha256",64);d.reference.streamId=number(get(vid,"stream"),0,INT32_MAX);d.reference.startOnAudio=video::Time::parse(str(vid,"offset",80).toStdString());validateCore(d);return d;
}
bool aliases (const juce::File& a,const juce::File& b)
{
    if(b==juce::File())return false;if(a==b||a.getLinkedTarget()==b.getLinkedTarget())return true;
    struct stat x{},y{};return ::stat(a.getFullPathName().toRawUTF8(),&x)==0&&::stat(b.getFullPathName().toRawUTF8(),&y)==0&&x.st_dev==y.st_dev&&x.st_ino==y.st_ino;
}
bool aliasesMedia (const juce::File& file,const Document& d)
{
    if(aliases(file,d.source.file)||aliases(file,d.reference.file))return true;
    if(d.source.assembly)for(const auto& m:d.source.assembly->members)if(aliases(file,m.file))return true;
    if(d.source.composition)for(const auto& m:d.source.composition->files())if(aliases(file,m.file))return true;
    return false;
}
void hook (const FileHooks& h,const char* stage){if(h.reached)h.reached(stage);if(h.fail&&h.fail(stage))throw Failure{juce::String(stage)=="space"?ErrorCode::insufficientDiskSpace:ErrorCode::writeFailed,"Injected session "+juce::String(stage)+" failure"};}
template<class F> Outcome run (F f)
{try{return f();}catch(const Failure& e){return{{},e.code,e.message};}catch(const std::exception& e){return{{},ErrorCode::invalidMetadata,e.what()};}catch(...){return{{},ErrorCode::invalidMetadata,"Session operation failed"};}}
}
Document capture (StandaloneSessionModel& s)
{
    if(s.hasActionableSectionProposal())throw std::runtime_error("Discard the transient proposal before saving; private candidate measurements are not saved as Current authority");
    if(s.sectionTrimDraft)throw std::runtime_error("Apply or discard the section-gain draft before saving; draft edits are not saved as applied audio");
    s.synchronizeSections();Document d;d.id=s.documentId;d.programmeId=s.programmeId;d.runtimeId=s.sessionId;d.revision=s.documentRevision;d.source=s.sourceDocument;d.current=s.currentStateDocument;d.applicability=s.programmeApplicability;d.sections=s.sectionDocument.map();d.markerImport=s.markerImport;d.jobRequirements=s.jobRequirements;d.proposalPolicy=s.proposalPolicy;d.lastAppliedProposalPolicy=s.lastAppliedProposalPolicy;d.lastAppliedProposalBinding=s.lastAppliedProposalBinding;d.reference=s.referenceVideo;d.audition=s.auditionDomain;
    const auto p=s.getTargetProfile(s.selectedTargetProfileIndex);d.profile={p.id,p.definitionVersion,p.digest()};if(s.restoredProfileId.isNotEmpty())d.profile={s.restoredProfileId,s.restoredProfileVersion,s.restoredProfileDigest};
    for(auto domain:{MeasurementDomain::source,MeasurementDomain::currentState})
    {
        const auto& m=domain==MeasurementDomain::source?s.source:s.currentState;
        if(!s.sourceContentInvalid&&m.hasMeasurement&&m.analysis&&(domain!=MeasurementDomain::currentState||!s.currentStateNeedsRemeasurement())
           &&m.analysis->context.sessionId==s.sessionId && m.analysis->measurementVersion==measurement::version
           &&m.analysis->domain==domain && m.analysis->measuredContentHash==s.sourceDocument.media.sha256
           &&m.analysis->identity==OfflineAnalysisService::identityFor(s.measurementContext(domain),domain,s.sourceDocument.media.sha256))
            (domain==MeasurementDomain::source?d.sourceCache:d.currentCache)=m.analysis;
    }
    if(s.sourceContentInvalid){d.sourceCache=s.unresolvedSourceCache;d.currentCache=s.unresolvedCurrentCache;}
    return d;
}
bool coreFits(const Document& d,juce::String& error)
{
    try{validateCore(d);std::array<Bytes,2> empty;require(juce::JSON::toString(coreValue(d,empty),true).getNumBytesAsUTF8()+4096<=coreLimit,"Annotations/provenance exceed 4 MiB session-core budget");return true;}
    catch(const std::exception& e){error=e.what();return false;}catch(const Failure& e){error=e.message;return false;}
}
Outcome save (const Document& d,const juce::File& file,JobControl& c,const FileHooks& h)
{
    return run([&]() -> Outcome
    {
        check(c);validateCore(d);require(file.hasFileExtension("levelscope"),"Use a .levelscope session filename");require(!aliasesMedia(file,d),"A session cannot overwrite source audio or reference picture");
        std::array<Bytes,2> caches;if(d.sourceCache)caches[0]=cacheBytes(*d.sourceCache,c);c.setProgress(.2f);if(d.currentCache)caches[1]=cacheBytes(*d.currentCache,c);
        auto json=juce::JSON::toString(coreValue(d,caches),true).toStdString();require(json.size()<=coreLimit,"Session core exceeds 4 MiB");
        uint64_t size=80+json.size()+caches[0].size()+caches[1].size();require(size<=maximumBytes,"Session exceeds 128 MiB");
        hook(h,"space");if(file.getParentDirectory().getBytesFreeOnVolume()<static_cast<int64_t>(size))throw Failure{ErrorCode::insufficientDiskSpace,"Insufficient disk space for session save"};
        auto temp=file.getSiblingFile("."+file.getFileName()+"."+juce::Uuid().toString()+".partial");
        struct Owned {juce::File file;int fd=-1;bool created=false;~Owned(){if(fd>=0)::close(fd);if(created)file.deleteFile();}} owned{temp};
        if(h.beforeCreate)h.beforeCreate(temp);
        hook(h,"create");owned.fd=::open(temp.getFullPathName().toRawUTF8(),O_WRONLY|O_CREAT|O_EXCL,0600);if(owned.fd<0)throw Failure{ErrorCode::writeFailed,"Cannot create session temporary file"};owned.created=true;
        uint64_t written=0;auto write=[&](const void* data,size_t n){auto* p=static_cast<const uint8_t*>(data);while(n){check(c);hook(h,"write");auto count=std::min<size_t>(n,65536);auto z=::write(owned.fd,p,count);if(z<0&&errno==EINTR)continue;if(z<=0)throw Failure{ErrorCode::writeFailed,"Session write failed; previous file retained"};p+=z;n-=static_cast<size_t>(z);written+=z;c.setProgress(.4f+.55f*static_cast<float>(double(written)/double(size)));}};
        Bytes header{'L','S','S','E','S','S','0','1'};u64(header,json.size());auto hash=digest(json.data(),json.size()).toStdString();header.insert(header.end(),hash.begin(),hash.end());write(header.data(),header.size());write(json.data(),json.size());for(auto& b:caches)write(b.data(),b.size());
        check(c);hook(h,"flush");if(::fsync(owned.fd)!=0)throw Failure{ErrorCode::writeFailed,"Session flush failed; previous file retained"};
        if(::close(owned.fd)!=0){owned.fd=-1;throw Failure{ErrorCode::writeFailed,"Session close failed; previous file retained"};}owned.fd=-1;
        check(c);hook(h,"replace");require(!aliasesMedia(file,d),"Destination now aliases source media");
        if(::rename(temp.getFullPathName().toRawUTF8(),file.getFullPathName().toRawUTF8())!=0)throw Failure{ErrorCode::writeFailed,"Atomic session replacement failed; previous file retained"};
        owned.created=false; // rename transferred ownership; never unlink a reused temporary path.
        c.setProgress(1);return{{},ErrorCode::none,"Session revision saved"};
    });
}
Outcome read (const juce::File& file,JobControl& c)
{
    return run([&]() -> Outcome
    {
        check(c);juce::FileInputStream input(file);require(input.openedOk(),"Cannot read session file");auto size=input.getTotalLength();require(size>=80&&size<=static_cast<int64_t>(maximumBytes),"Invalid or oversized session file");
        Bytes header(80);require(input.read(header.data(),80)==80&&std::memcmp(header.data(),"LSSESS01",8)==0,"Unsupported session container");Reader r{header,8,c};auto length=r.u();require(length<=coreLimit&&length<=static_cast<uint64_t>(size-80),"Truncated or oversized session core");
        std::string json(static_cast<size_t>(length),'\0');require(input.read(json.data(),static_cast<int>(length))==static_cast<int>(length),"Truncated session core");
        auto expected=juce::String::fromUTF8(reinterpret_cast<const char*>(header.data()+16),64);require(expected==digest(json.data(),json.size()),"Session core checksum mismatch");auto v=parse(json);auto d=core(v);c.setProgress(.15f);
        // Cache directory is non-authoritative: independently bad entries are discarded.
        auto entries=v["caches"];
        if(!entries.isArray()||entries.size()!=2)d.notices.add("Measurement caches unavailable: invalid cache directory; remeasure.");
        else for(int i=0;i<2;++i)
        {
            check(c);try
            {
                const auto e=entries[i];require(number(get(e,"domain"),0,1)==i,"Invalid cache domain");auto offset=integer<uint64_t>(get(e,"offset")),n=integer<uint64_t>(get(e,"size"));
                if(n==0)continue;require(n<=maximumBytes/2&&offset<=maximumBytes&&80+length+offset<=static_cast<uint64_t>(size)&&n<=static_cast<uint64_t>(size)-(80+length+offset),"Missing or truncated cache");require(hashValid(str(e,"sha256",64)),"Invalid cache checksum");
                input.setPosition(static_cast<int64_t>(80+length+offset));Bytes bytes(static_cast<size_t>(n));size_t done=0;while(done<bytes.size()){check(c);auto count=std::min<size_t>(65536,bytes.size()-done);require(input.read(bytes.data()+done,static_cast<int>(count))==static_cast<int>(count),"Cache read failed");done+=count;}
                require(digest(bytes)==str(e,"sha256",64),"Measurement cache checksum mismatch");auto a=cache(bytes,static_cast<MeasurementDomain>(i),c);(i==0?d.sourceCache:d.currentCache)=std::move(a);
            }
            catch(const Failure& e){if(e.code==ErrorCode::cancelled)throw;d.notices.add((i==0?"Source":"Current State")+juce::String(" cache discarded: ")+e.message+". Remeasure.");}
            catch(const std::exception& e){d.notices.add("Cache discarded: "+juce::String(e.what())+". Remeasure.");}
        }
        c.setProgress(1);return{std::move(d),ErrorCode::none,"Session core read; media verification required"};
    });
}
Outcome verify (Document d,JobControl& c)
{
    return run([&]() -> Outcome
    {
        check(c);validateCore(d);d.verifiedAudioMembers.clear();c.setProgress(0);d.audioResolved=!d.source.hasSource;d.videoResolved=d.reference.file==juce::File();
        if(d.source.hasSource)
        {
            const auto same=[](const MediaInfo& a,const MediaInfo& b){return a.sha256==b.sha256&&a.sampleRate==b.sampleRate&&a.frames==b.frames&&a.channels==b.channels&&a.bitsPerSample==b.bitsPerSample&&a.bytes==b.bytes&&a.format==b.format;};
            if(d.source.composition)
            {
                auto& recipe=*d.source.composition;d.audioResolved=true;
                const auto total=recipe.files().size();size_t index=0;
                std::vector<juce::File> seen;
                auto verifyMember=[&](const juce::String& id,const juce::File& file,MediaInfo& info)
                {
                    JobControl child;child.parent=&c;child.parentProgress=&c.progress;
                    child.progressBase=.7f*static_cast<float>(index++)/total;child.progressScale=.7f/total;
                    const auto valid=validateMedia(file,child);check(c);bool ok=valid.succeeded()&&same(valid.media,info);
                    for(const auto& prior:seen)if(samePhysicalFile(file,prior))ok=false;
                    seen.push_back(file);
                    if(ok){info=valid.media;d.verifiedAudioMembers.add(id);}
                    else{d.audioResolved=false;d.notices.add("Unresolved "+id+" ("+file.getFileName()+"): missing, changed, unreadable or duplicate component member. Relink identical verified content.");}
                };
                for(auto* component:{&recipe.me,&recipe.dialogue})
                {
                    if(component->assembly)
                    {
                        for(auto& member:component->assembly->members)verifyMember(component->id+"/"+member.id,member.file,member.media);
                        component->media=component->assembly->summary();
                    }
                    else verifyMember(component->id+"/file",component->file,component->media);
                }
                d.source.media=recipe.summary();
                if(!d.audioResolved)d.notices.add("Whole programme incomplete: no partial audition, processing or measurement; settings retained.");
            }
            else if(d.source.assembly)
            {
                auto& group=*d.source.assembly;d.audioResolved=true;
                for(size_t i=0;i<group.members.size();++i)
                {
                    auto& member=group.members[i];JobControl child;child.parent=&c;child.parentProgress=&c.progress;child.progressBase=.7f*static_cast<float>(i)/group.members.size();child.progressScale=.7f/group.members.size();
                    const auto valid=validateMedia(member.file,child);check(c);
                    bool ok=valid.succeeded()&&same(valid.media,member.media);
                    for(size_t j=0;j<i;++j)if(samePhysicalFile(member.file,group.members[j].file))ok=false;
                    if(ok){member.media=valid.media;d.verifiedAudioMembers.add(member.id);}
                    else{d.audioResolved=false;d.notices.add("Unresolved "+member.id+" ("+member.file.getFileName()+"): missing, changed, unreadable or duplicate file. Relink identical verified content.");}
                }
                d.source.media=group.summary();
                if(!d.audioResolved)d.notices.add("Whole programme incomplete: no partial audition, processing or measurement; settings retained.");
            }
            else
            {
                JobControl audioControl;audioControl.parent=&c;audioControl.parentProgress=&c.progress;audioControl.progressScale=.7f;
                auto valid=validateMedia(d.source.file,audioControl);check(c);
                if(valid.succeeded()&&same(valid.media,d.source.media)){d.source.media=valid.media;d.audioResolved=true;}
                else d.notices.add("Audio unresolved: missing, changed or unsupported source. Relink identical content; measurements unavailable.");
            }
            // Group caches remain dormant until every member verifies; preserve the M8 single-file policy.
            if(!d.audioResolved&&!d.source.assembly&&!d.source.composition){d.sourceCache.reset();d.currentCache.reset();}
        }
        if(d.reference.file!=juce::File())
        {
            c.setProgress(.75f);
            // Existing decoder verifies file identity and supported stream; embedded audio is never decoded.
            try {auto decoder=video::makeDecoder();auto asset=decoder->load(d.reference.file,d.reference.streamId,[&]{return c.isCancelled();});check(c);d.videoResolved=asset&&asset->sha256==d.reference.sha256;if(d.videoResolved){d.verifiedVideoBytes=asset->bytes;d.verifiedVideoModified=asset->modifiedMs;}}
            catch(const Failure&){throw;}catch(...){check(c);d.videoResolved=false;}
            if(!d.videoResolved)d.notices.add("Reference picture unresolved: missing, changed or unsupported. Audio identity is unchanged.");
        }
        auto prepared=std::make_shared<Prepared>();
        require(install(d,prepared->model,prepared->notice),prepared->notice);
        // This potentially large digest work belongs to the worker, not UI adoption.
        prepared->model.evaluation(MeasurementDomain::source);
        prepared->model.evaluation(MeasurementDomain::currentState);
        check(c);c.setProgress(1);return{std::move(d),ErrorCode::none,"Session media verification complete",std::move(prepared)};
    });
}
bool install (const Document& d,StandaloneSessionModel& s,juce::String& error)
{
    try
    {
        validateCore(d);require(!s.audioEditsBlocked()&&!s.job.isActive(),"Stop audition and finish/cancel the current job before opening a session");
        require(!d.source.hasSource || !d.audioResolved || d.source.matchesFileStamps(),"Source changed after verification; previous session retained");
        require(d.reference.file==juce::File() || !d.videoResolved || (d.reference.file.existsAsFile() && d.reference.file.getSize()==d.verifiedVideoBytes && d.reference.file.getLastModificationTime().toMilliseconds()==d.verifiedVideoModified),"Reference picture changed after verification; previous session retained");
        StandaloneSessionModel next;next.openedSessionSchema=d.loadedSchema;next.documentId=d.id;next.programmeId=d.programmeId;next.documentRevision=next.savedDocumentRevision=d.revision;next.sourceDocument=d.source;next.sourceDocument.routing=OutputRouting();next.currentStateDocument=d.current;next.currentStateDocument.needsRemeasurement=d.current.isInitialized;next.programmeApplicability=d.applicability;next.auditionDomain=d.audition;next.markerImport=d.markerImport;next.jobRequirements=d.jobRequirements;next.proposalPolicy=d.proposalPolicy;next.lastAppliedProposalPolicy=d.lastAppliedProposalPolicy;next.lastAppliedProposalBinding=d.lastAppliedProposalBinding;next.referenceVideo=d.reference;next.videoUnresolved=!d.videoResolved;
        if(!d.audioResolved){next.unresolvedSourceCache=d.sourceCache;next.unresolvedCurrentCache=d.currentCache;}
        next.sourceContentInvalid=d.source.hasSource&&!d.audioResolved;next.sourceDocument.importStatus=next.sourceContentInvalid?"Unresolved programme; relink identical audio":"Restored source; saved channel interpretation retained";
        next.currentStateDocument.sourcePathSnapshot=d.source.file.getFullPathName();next.currentStateDocument.basedOnSourceName=d.source.displayName;
        if(!d.sections.sections.empty())require(next.sectionDocument.restore(d.sections,sectionInputIdentity(next.measurementContext(MeasurementDomain::source))),"Invalid restored sections");
        bool profileFound=false;for(int i=0;i<next.numTargetProfiles;++i){auto p=next.getTargetProfile(i);if(p.id==d.profile.id&&p.definitionVersion==d.profile.version&&p.digest()==d.profile.digest){next.selectedTargetProfileIndex=i;profileFound=true;break;}}
        error=d.notices.joinIntoString("\n");if(!profileFound){next.restoredProfileId=d.profile.id;next.restoredProfileVersion=d.profile.version;next.restoredProfileDigest=d.profile.digest;error+="\nSaved profile unavailable: "+d.profile.id+" @ "+d.profile.version+"; choose an authoritative target before evaluation.";}
        for(auto domain:{MeasurementDomain::source,MeasurementDomain::currentState})
        {
            auto a=domain==MeasurementDomain::source?d.sourceCache:d.currentCache;if(!a||!d.audioResolved)continue;
            auto expected=next.measurementContext(domain);
            if(a->identity!=OfflineAnalysisService::identityFor(expected,domain,expected.media.sha256)||a->measurementVersion!=measurement::version||(domain==MeasurementDomain::currentState&&!d.current.isInitialized))
            {error+="\n"+juce::String(domain==MeasurementDomain::source?"Source":"Current State")+" cache incompatible with this build/platform/input. Remeasure.";continue;}
            auto restored=*a;restored.context.sessionId=next.sessionId;restored.context.ticket={0,next.revision};restored.context.sourceFile=next.sourceDocument.file;restored.context.media.modified=next.sourceDocument.media.modified;restored.context.assembly=next.sourceDocument.assembly;restored.context.composition=next.sourceDocument.composition;restored.restoredFromSession=true;
            require(measuredProgrammeForSnapshot(&restored,expected,domain,next.programmeApplicability).currentIdentityMatches,"Restored measurement failed identity validation");
            auto& m=domain==MeasurementDomain::source?next.source:next.currentState;m.hasMeasurement=true;m.analysis=std::make_shared<const AnalysisSnapshot>(std::move(restored));m.integratedLufs=a->measured.integratedLufs;m.loudnessRangeLu=a->measured.loudnessRangeLu;m.truePeakDbtp=a->measured.truePeakDbtp;m.statusText="Restored measurement; source and compatibility verified";
            if(domain==MeasurementDomain::currentState)next.currentStateDocument.needsRemeasurement=false;
        }
        if(next.currentStateDocument.isInitialized)next.currentStateDocument.statusText=next.currentState.hasMeasurement?"Restored measured Current State":"UNVERIFIED: restored gain retained; remeasure";
        next.auditionBusy=s.auditionBusy;s=std::move(next);return true;
    }
    catch(const Failure& f){error=f.message;return false;}catch(const std::exception& e){error=e.what();return false;}
}
bool adopt (const Document& d,Prepared& prepared,StandaloneSessionModel& s,juce::String& error)
{
    if(prepared.consumed||prepared.model.documentId!=d.id||prepared.model.documentRevision!=d.revision)
    {error="Prepared session is stale or belongs to a different document";return false;}
    if(s.audioEditsBlocked()||s.job.isActive()) {error="Stop and finish the current job before opening";return false;}
    if((d.source.hasSource&&d.audioResolved&&!d.source.matchesFileStamps())
       ||(d.reference.file!=juce::File()&&d.videoResolved&&(!d.reference.file.existsAsFile()
          ||d.reference.file.getSize()!=d.verifiedVideoBytes||d.reference.file.getLastModificationTime().toMilliseconds()!=d.verifiedVideoModified)))
    {error="Media changed after verification; previous session retained";return false;}
    prepared.model.auditionBusy=s.auditionBusy;s=std::move(prepared.model);prepared.consumed=true;error=prepared.notice;return true;
}
bool saved (StandaloneSessionModel& s,const Document& d,const juce::File& file)
{
    if(s.sessionId!=d.runtimeId||s.documentId!=d.id||d.revision>s.documentRevision)return false;
    s.documentFile=file;s.openedSessionSchema=schemaVersion;s.savedDocumentRevision=d.revision;return true;
}
bool relinkAudio (Document& d,const juce::File& f)
{if(!d.source.hasSource||d.source.assembly||d.source.composition)return false;d.source.file=f;d.source.fullPath=f.getFullPathName();d.source.displayName=f.getFileName();d.audioResolved=false;return true;}
bool relinkAudioMember (Document& d,const juce::String& id,const juce::File& file)
{
    if(d.source.composition)
    {
        for(auto* component:{&d.source.composition->me,&d.source.composition->dialogue})
        {
            if(component->assembly)
            {for(auto& m:component->assembly->members)if(id==component->id+"/"+m.id){m.file=file;d.audioResolved=false;return true;}}
            else if(id==component->id+"/file"){component->file=file;d.audioResolved=false;return true;}
        }
        return false;
    }
    if(!d.source.assembly)return false;
    for(auto& member:d.source.assembly->members)if(member.id==id){member.file=file;d.audioResolved=false;return true;}
    return false;
}
bool relinkVideo (Document& d,const juce::File& f)
{if(d.reference.file==juce::File())return false;d.reference.file=f;d.videoResolved=false;return true;}
}

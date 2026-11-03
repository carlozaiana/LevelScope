#include "JobRequirements.h"
#include <juce_cryptography/juce_cryptography.h>
#include <stdexcept>
#include <limits>

namespace levelscope::proposal
{
namespace
{
void require(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
juce::var definition(const JobRequirements& r)
{
    auto* o=new juce::DynamicObject();juce::var v(o);
    o->setProperty("format","levelscope-job-requirements");o->setProperty("schema",r.schema);
    o->setProperty("id",r.id);o->setProperty("revision",r.revision);o->setProperty("name",r.name);
    o->setProperty("base_id",r.baseId);o->setProperty("base_version",r.baseVersion);o->setProperty("base_digest",r.baseDigest);
    o->setProperty("method",r.method);o->setProperty("applicability",r.applicability);
    o->setProperty("provenance_kind","user-specified-job-requirement");o->setProperty("provenance",r.provenance);
    o->setProperty("maximum_lra_tenths",r.maximumLraTenths?juce::var(*r.maximumLraTenths):juce::var());
    if(r.schema == 2)
    {
        o->setProperty("integrated_tenths",r.integratedTenths?juce::var(*r.integratedTenths):juce::var());
        o->setProperty("integrated_tolerance_tenths",r.integratedToleranceTenths);
        o->setProperty("true_peak_tenths",r.truePeakTenths?juce::var(*r.truePeakTenths):juce::var());
        o->setProperty("maximum_momentary_tenths",r.maximumMomentaryTenths?juce::var(*r.maximumMomentaryTenths):juce::var());
        o->setProperty("maximum_short_term_tenths",r.maximumShortTermTenths?juce::var(*r.maximumShortTermTenths):juce::var());
    }
    o->setProperty("unavailable_reason",r.unavailableReason);return v;
}
}
juce::String JobRequirements::validationError() const
{
    if(unavailableReason.isNotEmpty())return "Job requirements unavailable: "+unavailableReason;
    if((schema!=1&&schema!=2)||revision<1||revision==std::numeric_limits<int>::max()||id.length()!=32||!id.containsOnly("0123456789abcdef"))return "Unsupported job-requirements identity/version.";
    if(name.trim().isEmpty()||name.length()>80||provenance.getNumBytesAsUTF8()>2048)return "Invalid requirements name/provenance size.";
    const auto base=ebuR128ProgrammeQc();
    if(baseId!=base.id||baseVersion!=base.definitionVersion||baseDigest!=base.digest())return "Unsupported base profile identity/version/digest.";
    if(method!="EBU Tech 3342-2023 whole-programme LRA"||applicability!="pre-produced-whole-programme-linear-audio")return "Unsupported LRA method/applicability.";
    if(maximumLraTenths&&(*maximumLraTenths<0||*maximumLraTenths>1000))return "Maximum LRA must be 0.0 to 100.0 LU in exact tenths.";
    if(schema == 2)
    {
        for(auto value : {integratedTenths, truePeakTenths, maximumMomentaryTenths, maximumShortTermTenths})
            if(value && (*value < -1000 || *value > 240)) return "Target values must be -100.0 to +24.0 in exact tenths, or No Limit.";
        if(integratedToleranceTenths < 0 || integratedToleranceTenths > 1000) return "IL tolerance must be 0.0 to 100.0 LU.";
    }
    return {};
}
juce::String JobRequirements::digest() const
{
    const auto bytes=juce::JSON::toString(definition(*this),true).toStdString();
    return juce::SHA256(bytes.data(),bytes.size()).toHexString();
}
juce::var JobRequirements::encode() const
{auto v=definition(*this);v.getDynamicObject()->setProperty("digest",digest());return v;}
JobRequirements JobRequirements::decode(const juce::var& v)
{
    const auto* o=v.getDynamicObject();require(o!=nullptr,"Requirements must be an object");
    const char* keys[]={"format","schema","id","revision","name","base_id","base_version","base_digest","method","applicability","provenance_kind","provenance","maximum_lra_tenths","unavailable_reason","digest"};
    const bool modern = o->getProperty("schema").isInt() && int(o->getProperty("schema")) == 2;
    require(o->getProperties().size()==(modern?20:15),"Missing or unsupported requirements fields");
    for(auto key:keys)require(o->hasProperty(key),"Missing requirements field");
    const auto str=[&](const char* k){const auto x=o->getProperty(k);require(x.isString(),"Requirements string has invalid type");return x.toString();};
    const auto integer=[&](const char* k){const auto x=o->getProperty(k);require(x.isInt(),"Requirements integer has invalid type");return int(x);};
    require(str("format")=="levelscope-job-requirements"&&str("provenance_kind")=="user-specified-job-requirement","Unsupported requirements format/authority");
    JobRequirements r;r.schema=integer("schema");r.revision=integer("revision");r.id=str("id");r.name=str("name");
    r.baseId=str("base_id");r.baseVersion=str("base_version");r.baseDigest=str("base_digest");r.method=str("method");r.applicability=str("applicability");r.provenance=str("provenance");r.unavailableReason=str("unavailable_reason");
    if(!o->getProperty("maximum_lra_tenths").isVoid())r.maximumLraTenths=integer("maximum_lra_tenths");
    if(modern)
    {
        const auto optional=[&](const char* key)->std::optional<int>
        {require(o->hasProperty(key),"Missing target field");return o->getProperty(key).isVoid()?std::nullopt:std::optional<int>(integer(key));};
        r.integratedTenths=optional("integrated_tenths");r.integratedToleranceTenths=integer("integrated_tolerance_tenths");
        r.truePeakTenths=optional("true_peak_tenths");r.maximumMomentaryTenths=optional("maximum_momentary_tenths");
        r.maximumShortTermTenths=optional("maximum_short_term_tenths");
    }
    require(str("digest")==r.digest(),"Requirements digest mismatch");
    const auto error=r.validationError();if(error.isNotEmpty())throw std::runtime_error(error.toStdString());
    require(juce::JSON::toString(v,true).getNumBytesAsUTF8()<=maximumBytes,"Requirements exceed file budget");return r;
}
JobRequirements JobRequirements::fromJson(const juce::String& s)
{
    require(s.getNumBytesAsUTF8()<=maximumBytes,"Requirements exceed 64 KiB budget");
    // The portable format is one flat object. Bound parser depth and reject
    // duplicate-key ambiguity before handing text to the general JSON parser.
    bool quoted=false,escaped=false,started=false,finished=false;int depth=0,colons=0;
    for(auto c:s)
    {
        if(finished){require(juce::CharacterFunctions::isWhitespace(c),"Trailing preset data is unsupported");continue;}
        if(!started)
        {
            if(juce::CharacterFunctions::isWhitespace(c))continue;
            require(c=='{',"Preset must start with one object");started=true;depth=1;continue;
        }
        if(quoted){if(escaped)escaped=false;else if(c=='\\')escaped=true;else if(c=='"')quoted=false;continue;}
        if(c=='"'){quoted=true;continue;}
        if(c=='[')throw std::runtime_error("Arrays are not supported in requirements presets");
        if(c=='{'&&++depth>1)throw std::runtime_error("Nested requirements data is unsupported");
        if(c=='}'){--depth;require(depth==0,"Malformed requirements object");finished=true;}
        if(c==':')++colons;
    }
    require(finished&&depth==0&&!quoted&&(colons==15||colons==20),"Invalid or duplicate requirements fields");
    juce::var v;const auto result=juce::JSON::parse(s,v);require(result.wasOk(),"Invalid requirements JSON");require(v.getDynamicObject() && v.getDynamicObject()->getProperties().size()==colons,"Duplicate requirements fields");return decode(v);
}
JobRequirements JobRequirements::unavailable(const juce::String& reason)
{JobRequirements r;r.schema=0;r.unavailableReason=reason.substring(0,512);return r;}
std::optional<int> JobRequirements::parseCeiling(const juce::String& input)
{
    const auto s=input.trim();if(s.isEmpty()||s.length()>5)return {};
    const auto dot=s.indexOfChar('.');
    int value=0;
    if(dot>=0)
    {
        if(dot==0||dot!=s.length()-2||!s.substring(0,dot).containsOnly("0123456789")||!s.substring(dot+1).containsOnly("0123456789"))return {};
        value=s.substring(0,dot).getIntValue()*10+s.substring(dot+1).getIntValue();
    }
    else {if(!s.containsOnly("0123456789"))return {};value=s.getIntValue()*10;}
    return value<=1000?std::optional<int>(value):std::nullopt;
}
std::optional<int> JobRequirements::parseSignedTenths(const juce::String& input)
{
    auto s=input.trim();const bool negative=s.startsWithChar('-');
    if(negative||s.startsWithChar('+'))s=s.substring(1);
    const auto value=parseCeiling(s);if(!value)return {};
    const int signedValue=negative?-*value:*value;
    return signedValue<=240?std::optional<int>(signedValue):std::nullopt;
}
JobRequirements JobRequirements::editable(const TargetProfile& p)
{
    JobRequirements r;r.schema=2;r.id=juce::Uuid().toString();r.name="Unsaved preset";
    const auto& l=p.limits;
    r.integratedTenths=l.hasIntegratedLoudnessTarget?std::optional<int>(juce::roundToInt(l.integratedLufs*10)):std::nullopt;
    r.integratedToleranceTenths=juce::roundToInt(l.integratedToleranceLu*10);
    r.truePeakTenths=l.hasTruePeakCeiling?std::optional<int>(juce::roundToInt(l.truePeakCeilingDbtp*10)):std::nullopt;
    r.maximumLraTenths=l.hasMaxLoudnessRange?std::optional<int>(juce::roundToInt(l.maxLoudnessRangeLu*10)):std::nullopt;
    r.maximumMomentaryTenths=l.hasMaxMomentary?std::optional<int>(juce::roundToInt(l.maxMomentaryLufs*10)):std::nullopt;
    r.maximumShortTermTenths=l.hasMaxShortTerm?std::optional<int>(juce::roundToInt(l.maxShortTermLufs*10)):std::nullopt;
    return r;
}
TargetProfile JobRequirements::applyTo(const TargetProfile& base) const
{
    auto p=base;p.baseProfileDigest=base.digest();p.jobRequirementsDigest=digest();p.jobRequirementsError=validationError();
    if(schema==2)
    {
        p.userSpecified=true;p.isSelectableForWorkflow=true;p.limits.valuesAreAuthoritative=false;
        p.id="user-target-"+id;p.definitionVersion=juce::String(revision);p.displayName=name;p.familyName="User target";
        p.method=MeasurementMethod::bs1770FullProgramme;p.applicability=Applicability::preProducedLinearProgramme;
        p.standardEdition="User-specified whole-programme constraints";p.referenceUrl=provenance;
        p.description="User-specified measured scalar constraints for a complete programme.";
        p.integratedProvenance=p.truePeakProvenance=p.lraProvenance=p.momentaryProvenance=p.shortTermProvenance=
            "User-specified target: "+name+". "+provenance;
        p.truePeakMeterToleranceDb=p.truePeakMeterBandwidthHz=0;p.meterToleranceNote={};
        p.limits.validationNote="Selected user constraints; broadcaster provenance is not independently verified.";
        auto& l=p.limits;l.hasIntegratedLoudnessTarget=bool(integratedTenths);l.integratedLufs=integratedTenths.value_or(0)/10.;
        l.integratedToleranceLu=integratedToleranceTenths/10.;l.hasTruePeakCeiling=bool(truePeakTenths);l.truePeakCeilingDbtp=truePeakTenths.value_or(0)/10.;
        l.hasMaxLoudnessRange=bool(maximumLraTenths);l.maxLoudnessRangeLu=maximumLraTenths.value_or(0)/10.;
        l.hasMaxMomentary=bool(maximumMomentaryTenths);l.maxMomentaryLufs=maximumMomentaryTenths.value_or(0)/10.;
        l.hasMaxShortTerm=bool(maximumShortTermTenths);l.maxShortTermLufs=maximumShortTermTenths.value_or(0)/10.;
        return p;
    }
    if(p.jobRequirementsError.isEmpty()&&(base.id!=baseId||base.definitionVersion!=baseVersion||base.digest()!=baseDigest))p.jobRequirementsError="Select the compatible EBU programme profile for these job requirements.";
    if(p.jobRequirementsError.isEmpty()&&maximumLraTenths)
    {
        p.limits.hasMaxLoudnessRange=true;p.limits.maxLoudnessRangeLu=*maximumLraTenths/10.;
        p.lraProvenance="User-specified job requirement: "+name+"; "+method+"; "+(provenance.isEmpty()?juce::String("no external document specified"):provenance)+". User provenance is not independently verified; not an EBU-mandated ceiling.";
    }
    return p;
}
juce::String JobRequirements::text() const
{
    const auto error=validationError();
    if(schema==2)
    {
        const auto value=[](std::optional<int> v){return v?juce::String(*v/10.,1):juce::String("No Limit");};
        return "Target: "+name+" | "+error+"\nIL: "+value(integratedTenths)+" LUFS +/- "+juce::String(integratedToleranceTenths/10.,1)
            +" LU; true peak: "+value(truePeakTenths)+" dBTP; LRA: "+value(maximumLraTenths)+" LU\nMax M: "+value(maximumMomentaryTenths)
            +" LUFS; Max S: "+value(maximumShortTermTenths)+" LUFS\nUser-specified. "+provenance+"\nDefinition: "+digest()+"\n";
    }
    return "Job: "+name+" | revision "+juce::String(revision)+"\n"+(error.isNotEmpty()?error:maximumLraTenths?"User-specified maximum LRA: "+juce::String(*maximumLraTenths/10.,1)+" LU (inclusive).":"No LRA acceptance limit specified.")
        +"\n"+method+"; short-programme LRA under 60 seconds cannot satisfy a required ceiling.\n"+"Provenance: "+(provenance.isEmpty()?juce::String("user-specified; no external document specified"):provenance)+"\nRequirements: "+digest()+"\n";
}
}

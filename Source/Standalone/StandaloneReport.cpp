#include "StandaloneExport.h"
#include <sstream>
#include <iomanip>
#include <locale>
#include <cmath>
#include <stdexcept>

namespace levelscope::standalone::exporting
{
namespace
{
juce::var object() { return juce::var(new juce::DynamicObject()); }
void set(juce::var& o,const char* key,const juce::var& v) { o.getDynamicObject()->setProperty(key,v); }
juce::String exact(double value)
{std::ostringstream out;out.imbue(std::locale::classic());out<<std::hexfloat<<value;return out.str();}
juce::String escaped(const juce::String& s)
{
    // Single pass: repeated String::replace copies the full multi-megabyte
    // proposal report once per quote/angle bracket and is quadratic.
    const auto input=s.toStdString();std::string output;output.reserve(input.size());
    for(char c:input)switch(c)
    {case '&':output+="&amp;";break;case '<':output+="&lt;";break;case '>':output+="&gt;";break;
     case '\"':output+="&quot;";break;case '\'':output+="&#39;";break;default:output+=c;break;}
    return juce::String::fromUTF8(output.data(),int(output.size()));
}
const char* status(measurement::Status s)
{
    switch(s){case measurement::Status::available:return "available";case measurement::Status::silence:return "silence";
    case measurement::Status::belowGate:return "below_gate";default:return "insufficient_duration";}
}
juce::var metric(const measurement::Metric& m,const char* unit)
{
    auto o=object();const bool usable=m.status==measurement::Status::available&&m.value&&std::isfinite(*m.value);
    set(o,"status",m.value&&!std::isfinite(*m.value)?"invalid_non_finite":m.status==measurement::Status::available&&!m.value?"invalid_missing_value":status(m.status));set(o,"unit",unit);
    set(o,"value",usable?juce::var(*m.value):juce::var());
    if(usable)set(o,"value_hex",exact(*m.value));return o;
}
juce::var evaluation(const proposal::Evaluation& e)
{
    auto o=object();set(o,"verdict",proposal::verdictName(e.verdict));set(o,"message",e.message);
    set(o,"identity",e.identity);set(o,"profile_digest",e.profileDigest);set(o,"job_requirements_digest",e.jobRequirementsDigest);set(o,"evaluated_rules_digest",e.rulesDigest);set(o,"snapshot_identity",e.snapshotIdentity);
    juce::Array<juce::var> rules;
    for(const auto& r:e.rules)
    {
        auto row=object();set(row,"id",r.id);set(row,"verdict",proposal::verdictName(r.verdict));set(row,"reason",r.reason);
        set(row,"provenance",r.provenance);set(row,"units",r.units);
        for(const auto& v:std::initializer_list<std::pair<const char*,std::optional<double>>>{{"measured",r.measured},{"lower",r.lower},{"upper",r.upper},{"delta",r.delta}})
            set(row,v.first,v.second&&std::isfinite(*v.second)?juce::var(*v.second):juce::var());
        rules.add(row);
    }
    set(o,"rules",rules);return o;
}
juce::var snapshot(const AnalysisSnapshot* a,const ServiceContext& expected,MeasurementDomain domain,
                   const TargetProfile& profile,proposal::Applicability applicability,const juce::String& mediaIssue={})
{
    auto o=object();const auto input=measuredProgrammeForSnapshot(a,expected,domain,applicability);
    const bool valid=a&&input.currentIdentityMatches&&mediaIssue.isEmpty();
    set(o,"domain",domain==MeasurementDomain::source?"Source":domain==MeasurementDomain::currentState?"Current State":"Exported artifact");
    if(expected.composition)set(o,"input_description",domain==MeasurementDomain::source?"Original confirmed composed programme":"Composed programme with explicit Current processing plan");
    set(o,"measurement_available",valid);
    set(o,"status",valid?"valid_measured_snapshot":a?"stale_or_incompatible":"unavailable");
    set(o,"evaluation",evaluation(evaluateSnapshot(valid?a:nullptr,expected,domain,profile,applicability)));
    if(!valid){set(o,"reason",mediaIssue.isNotEmpty()?mediaIssue:a?"Snapshot does not match the exact current input/processing context":"No valid measurement for this state");return o;}
    set(o,"measurement_identity",a->identity);set(o,(expected.composition||expected.assembly)?"measured_programme_definition_digest":"measured_content_sha256",a->measuredContentHash);
    set(o,"measurement_method",a->measurementVersion);set(o,"build",a->context.buildVersion);
    set(o,"platform",a->context.platformAndFloatingPointMode);set(o,"original_runtime",a->originSessionId);
    set(o,"restored_from_session",a->restoredFromSession);set(o,"frames",juce::String(a->detail.frames));
    set(o,"sample_rate",a->detail.sampleRate);set(o,"lra_stable",a->detail.lraStable);
    auto values=object();set(values,"integrated",metric(a->detail.integrated,"LUFS"));set(values,"true_peak",metric(a->detail.truePeak,"dBTP"));
    set(values,"sample_peak",metric(a->detail.samplePeak,"dBFS"));set(values,"lra",metric(a->detail.lra,"LU"));
    set(values,"momentary_max",metric(a->detail.momentaryMax,"LUFS"));set(values,"short_term_max",metric(a->detail.shortTermMax,"LUFS"));set(o,"metrics",values);
    juce::Array<juce::var> peaks;
    for(size_t i=0;i<a->detail.channelTruePeaks.size();++i)
    {
        auto row=object();set(row,"channel",static_cast<int>(i+1));
        if(i<static_cast<size_t>(a->context.fileLayout.count))set(row,"role",juce::AudioChannelSet::getAbbreviatedChannelTypeName(a->context.fileLayout.roles[i]));
        const double sp=i<a->detail.channelSamplePeaks.size()?a->detail.channelSamplePeaks[i]:std::numeric_limits<double>::quiet_NaN(),tp=a->detail.channelTruePeaks[i];
        set(row,"sample_peak",metric({sp>0?std::optional<double>(20*std::log10(sp)):std::nullopt,sp>0?measurement::Status::available:measurement::Status::silence},"dBFS"));
        if(!std::isfinite(sp)){auto absent=object();set(absent,"status",i<a->detail.channelSamplePeaks.size()?"invalid_non_finite":"unavailable");set(absent,"value",{});set(row,"sample_peak",absent);}
        set(row,"true_peak",metric({tp>0?std::optional<double>(20*std::log10(tp)):std::nullopt,tp>0?measurement::Status::available:measurement::Status::silence},"dBTP"));
        if(!std::isfinite(tp)){auto absent=object();set(absent,"status","invalid_non_finite");set(absent,"value",{});set(row,"true_peak",absent);}
        if(i<a->detail.samplePeakFrames.size())set(row,"sample_peak_frame",juce::String(a->detail.samplePeakFrames[i]));
        if(i<a->detail.truePeakEighthFrames.size())set(row,"true_peak_eighth_frame",juce::String(a->detail.truePeakEighthFrames[i]));peaks.add(row);
    }
    set(o,"per_channel_peaks",peaks);return o;
}
void jsonValue(std::ostringstream& out,const juce::var& v,int depth)
{
    if(depth>24)throw std::runtime_error("Report nesting exceeds limit");
    if(v.isVoid()||v.isUndefined()){out<<"null";return;}
    if(v.isBool()){out<<(static_cast<bool>(v)?"true":"false");return;}
    if(v.isDouble()){const double n=static_cast<double>(v);if(std::isfinite(n))out<<std::setprecision(17)<<n;else out<<"null";return;}
    if(v.isInt()||v.isInt64()){out<<static_cast<juce::int64>(v);return;}
    if(v.isString()){out<<juce::JSON::toString(v,true).toStdString();return;}
    if(const auto* a=v.getArray())
    {out<<'[';bool first=true;for(const auto& x:*a){if(!first)out<<',';first=false;jsonValue(out,x,depth+1);}out<<']';return;}
    if(const auto* o=v.getDynamicObject())
    {
        out<<'{';const auto& props=o->getProperties();for(int i=0;i<props.size();++i)
        {if(i)out<<',';out<<'\n'<<std::string(static_cast<size_t>(depth+1)*2,' ')<<juce::JSON::toString(props.getName(i).toString(),true).toStdString()<<": ";jsonValue(out,props.getValueAt(i),depth+1);}
        out<<'\n'<<std::string(static_cast<size_t>(depth)*2,' ')<<'}';return;
    }
    throw std::runtime_error("Unsupported report value");
}
}
juce::var report(const Request& r,const Candidate* c,bool publishing)
{
    if(c&&!c->matchesAudio(r))throw std::runtime_error("Artifact report does not match captured audio/runtime");
    auto o=object();set(o,"schema","levelscope-measurement-report-v5");
    set(o,"kind",c?"verified-artifact-measurement":"measurement-only-no-exported-audio");
    set(o,"publication",publishing?"published-verified-pcm":c?"NOT EXPORTED - verified staging artifact":"NO EXPORTED AUDIO");
    set(o,"scope","Selected whole-programme loudness, true-peak and any explicit job LRA requirement only. No complete broadcaster delivery, BWF/BW64, metadata, synchronization or future encoded-media certification.");
    set(o,"document_id",r.documentId);set(o,"programme_id",r.programmeId);set(o,"document_revision",juce::String(r.documentRevision));
    set(o,"session_unsaved",r.unsaved);set(o,"session_basename",r.sessionFile.getFileName());set(o,"runtime_id",r.current.sessionId);
    set(o,"app","LevelScope Standalone");set(o,"build",r.current.buildVersion);set(o,"platform",r.current.platformAndFloatingPointMode);
    set(o,"cpu",juce::SystemStats::getCpuModel());set(o,"processing_version",r.current.processing.engineVersion);
    set(o,"processing_gain_db",r.current.processing.staticGainDb);set(o,"processing_gain_hex",exact(r.current.processing.staticGainDb));
    auto plan=object();set(plan,"engine",r.current.processing.engineVersion);set(plan,"global_gain_db",r.current.processing.staticGainDb);set(plan,"global_gain_hex",exact(r.current.processing.staticGainDb));
    set(plan,"order","original programme/composition -> section gain -> static-gain-v1 -> monitor or final export quantization");
    juce::Array<juce::var> gains,ramps;
    if(r.current.processing.sectionGains)
    {
        const auto g=r.current.processing.sectionGains->canonical();set(plan,"input_identity",g.inputIdentity);set(plan,"audio_key",g.audioKey());
        set(plan,"transition","section-gain-v1: centred 100 ms maximum, linear dB, shortened to half adjacent spans, half-open audio-frame intervals");
        for(const auto& row:g.intervals){auto v=object();set(v,"begin_frame",juce::String(row.begin));set(v,"end_frame",juce::String(row.end));set(v,"gain_db",row.db);set(v,"gain_hex",exact(row.db));gains.add(v);}
        for(size_t n=1;n<g.intervals.size();++n){auto v=object();const auto b=g.intervals[n].begin,rad=g.radius(n);set(v,"begin_frame",juce::String(b-rad));set(v,"end_frame",juce::String(b+rad));ramps.add(v);}
    }
    set(plan,"applied_intervals",gains);set(plan,"transitions",ramps);set(o,"current_processing_plan",plan);
    set(o,"programme_frames",juce::String(r.current.media.frames));set(o,"sample_rate",r.current.media.sampleRate);
    set(o,"duration_rational_seconds",juce::String(r.current.media.frames)+"/"+juce::String(r.current.media.sampleRate,0));
    set(o,"positions","Audio-source frames; true-peak positions in eighth-frames. Video offsets and measurement EOF drain do not extend programme duration.");
    auto profile=object();const auto& p=r.current.profile;
    set(profile,"id",p.id);set(profile,"definition_version",p.definitionVersion);set(profile,"digest",p.baseProfileDigest.isEmpty()?p.digest():p.baseProfileDigest);set(profile,"evaluated_requirements_digest",p.digest());
    set(profile,"edition",p.standardEdition);set(profile,"reference",p.referenceUrl);set(profile,"applicability",static_cast<int>(r.applicability));
    set(profile,"applicability_scope","Pre-produced whole-programme linear-audio QC only; excludes live, deliberate lower targets, short-form and distribution-specific requirements.");
    set(profile,"meter_tolerance_documentary_only",p.meterToleranceNote);set(o,"profile",profile);set(o,"job_requirements",r.jobRequirements.encode());set(o,"job_requirements_status",r.jobRequirements.text());
    juce::Array<juce::var> members;
    if(r.current.composition)
    {
        const auto& recipe=*r.current.composition;
        auto composition=object();set(composition,"schema",ProgrammeComposition::schemaVersion);set(composition,"engine",ProgrammeComposition::engineVersion);
        set(composition,"definition_digest",recipe.definitionIdentity());set(composition,"confirmation",recipe.confirmation);
        set(composition,"confirmation_meaning","Complementary, correctly balanced stems intended to be combined; exact sample-zero alignment confirmed.");
        set(composition,"preset",static_cast<int>(recipe.preset));juce::Array<juce::var> components;
        for(const auto* c:{&recipe.me,&recipe.dialogue})
        {
            auto component=object();set(component,"id",c->id);set(component,"purpose",c->purpose);set(component,"language",c->language);set(component,"version_label",c->versionLabel);
            set(component,"audio_identity",c->audioIdentity());set(component,"component_gain_db",0);set(component,"sample_zero_offset_frames","0");
            juce::Array<juce::var> files,roles,mapping;
            for(const auto& m:c->files())
            {auto f=object();set(f,"slot",c->id+"/"+m.id);set(f,"basename",m.file.getFileName());set(f,"sha256",m.media.sha256);files.add(f);members.add(f);}
            for(int n=0;n<c->layout.count;++n)roles.add(juce::AudioChannelSet::getAbbreviatedChannelTypeName(c->layout.roles[n]));
            for(auto n:recipe.map(c==&recipe.dialogue))mapping.add(n<0?juce::var():juce::var(n+1));
            set(component,"files",files);set(component,"ordered_input_roles",roles);set(component,"output_to_input_channels",mapping);components.add(component);
        }
        set(composition,"components",components);set(o,"composition",composition);
    }
    else if(r.current.assembly)
    {
        for(const auto& m:r.current.assembly->members)
        {auto row=object();set(row,"basename",m.file.getFileName());set(row,"sha256",m.media.sha256);set(row,"role",juce::AudioChannelSet::getAbbreviatedChannelTypeName(m.role));
         set(row,"sample_zero_offset_frames",juce::String(m.offsetFrames));members.add(row);}
        set(o,"assembly_definition",r.current.assembly->definitionIdentity());set(o,"alignment_confirmation",r.current.assembly->alignmentConfirmation);
    }
    else if(r.current.media.sha256.isNotEmpty())
    {auto row=object();set(row,"basename",r.current.sourceFile.getFileName());set(row,"sha256",r.current.media.sha256);members.add(row);}
    set(o,"source_members",members);
    juce::Array<juce::var> roles;for(int n=0;n<r.current.fileLayout.count;++n)roles.add(juce::AudioChannelSet::getAbbreviatedChannelTypeName(r.current.fileLayout.roles[n]));set(o,"ordered_input_roles",roles);
    set(o,"source",snapshot(r.sourceAnalysis.get(),r.source,MeasurementDomain::source,p,r.applicability,r.mediaIssue));
    set(o,"current",snapshot(r.currentAnalysis.get(),r.current,MeasurementDomain::currentState,p,r.applicability,r.mediaIssue));
    if(c)
    {
        auto artifact=snapshot(&c->analysis(),c->analysis().context,MeasurementDomain::renderedOutput,p,r.applicability);
        set(artifact,"sha256",c->hash());set(artifact,"basename","programme.wav");set(artifact,"container",c->spec().format.rf64?"RF64":"RIFF/WAVE");
        set(artifact,"bits",24);set(artifact,"channel_mask",juce::String::toHexString(static_cast<int>(c->spec().format.mask)));
        set(artifact,"file_bytes",juce::String(c->spec().format.fileBytes));set(artifact,"format_identity",c->spec().identity);
        set(artifact,"quantizer",pcm24::version);set(artifact,"dither_seed_hex",c->spec().seedHex);set(artifact,"dither_applications",1);
        set(artifact,"dither_note","TPDF applied once at final PCM24 quantization, outside Current State. Even zero-gain export is not a bit-identical copy.");
        set(artifact,"verification_gain_db",0);juce::Array<juce::var> map;for(auto n:c->spec().outputToInput)map.add(n+1);set(artifact,"output_to_input_channels",map);set(o,"artifact",artifact);
    }
    juce::Array<juce::var> sections;
    if(r.sections.sections.size()>10000)throw std::runtime_error("Too many report annotations");
    for(const auto& s:r.sections.sections){auto row=object();set(row,"id",juce::String(s.id));set(row,"name",s.name);set(row,"begin_frame",juce::String(s.begin));set(row,"end_frame",juce::String(s.end));sections.add(row);}
    set(o,"sections_annotations_only",sections);
    set(o,"proposal_preferences",r.proposalPolicy.encode());set(o,"proposal_preferences_digest",r.proposalPolicy.digest());
    set(o,"proposal_inspection_not_current_measurement",r.proposalInspection);
    set(o,"last_applied_proposal_policy_historical",r.lastAppliedProposalPolicy?r.lastAppliedProposalPolicy->encode():juce::var());
    set(o,"last_applied_proposal_binding_historical",r.lastAppliedProposalBinding);
    set(o,"last_applied_proposal_inspection_transient",r.lastAppliedProposalInspection);
    set(o,"proposal_provenance_scope","Preferences are preservation choices, not delivery rules. Last applied policy is historical and may differ from current processing/preferences. Full transient proposal inspection is reportable before closing; sessions retain policy/binding and the exact applied envelope, not private assessment authority.");
    set(o,"marker_import_provenance",r.markerImport?proposal::markers::encode(*r.markerImport):juce::var());
    set(o,"annotations_match_last_import",r.markerImport&&r.sections.valid()&&r.markerImport->imported.digest()==r.sections.digest());auto video=object();set(video,"basename",r.video.file.getFileName());set(video,"sha256",r.video.sha256);
    set(video,"offset_rational",juce::String(r.video.startOnAudio.text()));set(video,"is_sync_evidence",false);set(o,"reference_picture_context",video);
    return o;
}
juce::String reportJson(const juce::var& value)
{std::ostringstream out;out.imbue(std::locale::classic());jsonValue(out,value,0);out<<'\n';return out.str();}
juce::String reportHtml(const juce::var& model)
{
    const bool artifact=model.hasProperty("artifact");
    juce::String html="<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>LevelScope measurement report</title>"
        "<style>body{font:16px system-ui,sans-serif;margin:2em;max-width:80em;color:#18212b}table{border-collapse:collapse}td,th{padding:.6em;border:1px solid #bbc}pre{white-space:pre-wrap;overflow-wrap:anywhere}h1{font-size:1.5em}</style><body><h1>";
    html+=artifact?"VERIFIED ARTIFACT MEASUREMENT":"MEASUREMENT-ONLY - NO EXPORTED AUDIO";
    html+="</h1><p>"+escaped(model["publication"].toString())+"</p><p>"+escaped(model["scope"].toString())+"</p>";
    html+="<p>Programme: "+escaped(model["programme_id"].toString())+" | Processing gain: "+escaped(model["processing_gain_db"].toString())+" dB | "+(static_cast<bool>(model["session_unsaved"])?"Unsaved captured session revision":"Saved session revision")+"</p>";
    if(model.hasProperty("composition"))html+="<p><strong>Source: original confirmed M&amp;E + dialogue composition.</strong> Unity component gains; Current State applies the recorded section and global processing plan. Component purpose, language/version, mapping and individual file hashes are recorded below. This is not dialogue analysis or verification of creative balance.</p>";
    html+="<h2>Selected job requirements</h2><pre>"+escaped(model["job_requirements_status"].toString())+"</pre>";
    html+="<table><tr><th>Measured state</th><th>Status / selected rules</th><th>Integrated LUFS</th><th>True peak dBTP</th><th>LRA LU</th></tr>";
    for(const auto& key:artifact?juce::StringArray{"source","current","artifact"}:juce::StringArray{"source","current"})
    {
        const auto state=model[juce::Identifier(key)];html+="<tr><th>"+escaped(state["domain"].toString())+"</th><td>"+escaped(state["status"].toString())+" / "+escaped(state["evaluation"]["verdict"].toString())+"</td>";
        for(const auto& metricName:{"integrated","true_peak","lra"})
        {
            const auto v=state["metrics"][metricName];const auto n=v["value"];
            html+="<td>"+(n.isDouble()||n.isInt()?escaped(juce::String(static_cast<double>(n),6)):escaped(v.hasProperty("status")?v["status"].toString():"unavailable"))+"</td>";
        }
        html+="</tr>";
    }
    html+="</table><p>Display rounding does not determine a verdict. LRA is informative unless an explicit job ceiling is selected; required unstable or unavailable LRA cannot pass. Stability and individual rule reasons are recorded below.</p>";
    if(artifact)html+="<p>Actual PCM SHA-256: <code>"+escaped(model["artifact"]["sha256"].toString())+"</code></p><p>"+escaped(model["artifact"]["dither_note"].toString())+"</p>";
    html+="<h2>Exact measurements, reasons and provenance</h2><pre>"+escaped(reportJson(model))+"</pre></body></html>";return html;
}
}

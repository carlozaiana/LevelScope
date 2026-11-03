#include "AppliedProposalExplanation.h"
#include "StandaloneSections.h"
#include <juce_cryptography/juce_cryptography.h>
#include <cstring>
namespace levelscope::standalone
{
juce::String proposalPlanIdentity(const ProcessingPlan& p)
{
    uint64_t gain;std::memcpy(&gain,&p.staticGainDb,sizeof(gain));
    const auto key=juce::String(p.schemaVersion)+";"+p.engineVersion+";"+juce::String::toHexString(juce::int64(gain))+";"+p.graphJson+";"+(p.sectionGains?p.sectionGains->audioKey():juce::String());
    return juce::SHA256(key.toRawUTF8(),key.getNumBytesAsUTF8()).toHexString();
}
ProposalMeasurementRecord proposalMeasurementRecord(const AnalysisSnapshot& a,const proposal::TargetProfile& target,proposal::Applicability applicability)
{
    ProposalMeasurementRecord r;r.identity=a.identity;r.content=a.measuredContentHash;r.method=a.measurementVersion;
    r.build=a.context.buildVersion;r.platform=a.context.platformAndFloatingPointMode;r.domain=a.domain;
    r.originSession=a.originSessionId.isEmpty()?a.context.sessionId:a.originSessionId;
    r.originJob=a.originSessionId.isEmpty()?a.context.ticket.id:a.originJobId;
    r.originRevision=a.originSessionId.isEmpty()?a.context.ticket.sessionRevision:a.originSessionRevision;
    r.identicalSource=a.identicalSourceMeasurementIdentity;
    const auto& d=a.detail;r.metrics={d.integrated,d.truePeak,d.lra,d.momentaryMax,d.shortTermMax,d.samplePeak};
    r.evaluation=evaluationText(evaluateSnapshot(&a,a.context,a.domain,target,applicability));return r;
}
juce::String proposalMeasurementText(const ProposalMeasurementRecord& r,const juce::String& label)
{
    juce::String text=label+" (measured)\n";
    const char* names[]={"IL LUFS","True peak dBTP","LRA LU","Max M LUFS","Max S LUFS","Sample peak dBFS"};
    const char* states[]={"available","silence","insufficient duration","below gate"};
    for(size_t n=0;n<r.metrics.size();++n)text+=juce::String(names[n])+" = "+(r.metrics[n].value?juce::String(*r.metrics[n].value,12):juce::String(states[int(r.metrics[n].status)]))+"\n";
    return text+r.evaluation+"\nMeasurement "+r.identity+"\nMethod "+r.method+" / build "+r.build+" / "+r.platform
        +"\nOrigin session "+r.originSession+" / job "+juce::String(r.originJob)+" / revision "+juce::String(r.originRevision)
        +(r.identicalSource.isNotEmpty()?"\nInitial Current from identical Source measurement "+r.identicalSource:juce::String())+"\n";
}
AppliedProposalExplanation captureProposalExplanation(const SectionProposalRequest& r,const proposal::JobRequirements& requirements,
    const juce::String& programmeId,const proposal::TrimProposal& global,const SectionProposalResult* result)
{
    if(!r.source||!r.current)throw std::runtime_error("Measured Source and Current baseline required for an applied explanation");
    AppliedProposalExplanation h;h.programmeId=programmeId;h.inputIdentity=sectionInputIdentity(r.currentContext);
    h.binding=result?result->identity():global.identity;h.method=result?proposal::sectionTrimMethodVersion:proposal::trimMethodVersion;
    h.before=r.currentContext.processing;h.after=h.before;
    h.roles=r.currentContext.fileLayout;h.sampleRate=int(r.currentContext.media.sampleRate);h.frames=r.currentContext.media.frames;
    if(result&&result->status()==SectionProposalStatus::ready)h.after=result->plan();
    else if(global.totalGainDb)h.after.staticGainDb=*global.totalGainDb;
    h.beforeIdentity=proposalPlanIdentity(h.before);h.afterIdentity=proposalPlanIdentity(h.after);
    h.target=r.currentContext.profile;h.requirements=requirements;h.preferences=r.policy;h.applicability=r.applicability;
    h.annotations=result?result->originalAnnotations():r.annotations;h.global=global;h.selection.policy=r.policy;
    h.source=proposalMeasurementRecord(*r.source,h.target,h.applicability);h.baseline=proposalMeasurementRecord(*r.current,h.target,h.applicability);
    if(result)
    {
        h.narrative=result->narrativeText();h.selection=result->selection();h.passes=result->passes();
        if(result->assessment())h.candidate=proposalMeasurementRecord(*result->assessment(),h.target,h.applicability);
    }
    else h.narrative=global.explanation;
    return h;
}
juce::String proposalTargetText(const proposal::TargetProfile& p)
{
    const auto& l=p.limits;
    return "Target "+p.displayName+" / "+p.id+" / "+p.definitionVersion+" / "+p.digest()+"\n"
        +"IL "+(l.hasIntegratedLoudnessTarget?juce::String(l.integratedLufs,9)+" LUFS +/- "+juce::String(l.integratedToleranceLu,9)+" LU":juce::String("No Limit"))
        +" / True peak "+(l.hasTruePeakCeiling?juce::String(l.truePeakCeilingDbtp,9)+" dBTP":juce::String("No Limit"))
        +" / LRA "+(l.hasMaxLoudnessRange?juce::String(l.maxLoudnessRangeLu,9)+" LU":juce::String("No Limit"))
        +" / Max M "+(l.hasMaxMomentary?juce::String(l.maxMomentaryLufs,9)+" LUFS":juce::String("No Limit"))
        +" / Max S "+(l.hasMaxShortTerm?juce::String(l.maxShortTermLufs,9)+" LUFS":juce::String("No Limit"))+"\n"
        +p.description+"\n"+p.standardEdition+" / "+p.referenceUrl+"\n"
        +"Requirement provenance\nIL "+p.integratedProvenance+"\nTrue peak "+p.truePeakProvenance+"\nLRA "+p.lraProvenance
        +"\nMax M "+p.momentaryProvenance+"\nMax S "+p.shortTermProvenance+"\n"+p.meterToleranceNote+"\n";
}
juce::String AppliedProposalExplanation::text() const
{
    juce::MemoryOutputStream t;t << juce::String("Last successfully applied proposal — historical explanation\n")+method+"\n"+narrative+"\n"
        +"Programme "+programmeId+" / input "+inputIdentity+"\nProposal binding "+binding
        +"\nBefore plan "+beforeIdentity+"\nAfter plan "+afterIdentity+"\n"
        +proposalTargetText(target)
        +requirements.text()+"\n"+preferences.text()+"\n"
        +proposalMeasurementText(source,"Source baseline")+"\n"+proposalMeasurementText(baseline,"Current baseline before Apply")+"\n"
        +proposal::trimProposalText(global)+"\nApplied global gain "+juce::String(after.staticGainDb,12)+" dB (before "+juce::String(before.staticGainDb,12)+" dB).\n";
    t << "Programme " << sampleRate << " Hz / " << frames << " frames. Confirmed file order";
    for(int n=0;n<roles.count;++n)t << " " << n+1 << "=" << juce::AudioChannelSet::getAbbreviatedChannelTypeName(roles.roles[size_t(n)]);
    t << "\n";
    if(!selection.corrections.empty())t << proposal::sectionSeedText(selection)+"\n";
    t << "Original annotation names and intervals at application\n";
    for(const auto& a:annotations.sections)t << "Section "+juce::String(a.id)+" ["+juce::String(a.begin)+", "+juce::String(a.end)+") "+a.name+"\n";
    t << "Applied section plateaus (transitions use the saved processing engine)\n";
    if(after.sectionGains)for(const auto& i:after.sectionGains->intervals)t << "["+juce::String(i.begin)+", "+juce::String(i.end)+") local "+juce::String(i.db,12)+" dB / local + global "+juce::String(i.db+after.staticGainDb,12)+" dB\n";
    else t << "Local gain 0 dB throughout. Combined gain "+juce::String(after.staticGainDb,12)+" dB.\n";
    t << "Private passes "+juce::String(passes)+" / 3\n";
    t << (candidate?proposalMeasurementText(*candidate,"Private candidate at proposal time"):juce::String("No private candidate was measured. Global values above are predictions.\n"));
    t << "\nHistorical evidence does not verify current processing. Current verification is shown independently in the main window.\n";return t.toUTF8();
}
}

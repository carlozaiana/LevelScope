#include "StandaloneRegionProposal.h"
#include "AppliedProposalExplanation.h"
#include <stdexcept>

namespace levelscope::standalone
{
namespace
{
void require(bool b,const char* why){if(!b)throw std::invalid_argument(why);}
juce::String hash(const juce::String& s){return juce::SHA256(s.toRawUTF8(),s.getNumBytesAsUTF8()).toHexString();}
}
RegionAssessmentReceipt::RegionAssessmentReceipt(const AnalysisSnapshot& actual,const SectionProposalRequest& r)
{
    auto expected=r.currentContext;expected.processing=actual.context.processing;
    require(actual.domain==MeasurementDomain::privateCandidate&&!actual.restoredFromSession
        &&OfflineAnalysisService::hasCompleteCoverage(actual)&&measurementMatchesContext(actual,expected,MeasurementDomain::privateCandidate),"Incomplete, restored or mismatched worker assessment");
    require(actual.detail.channelTruePeaks.size()==size_t(expected.fileLayout.count)&&actual.detail.channelSamplePeaks.size()==size_t(expected.fileLayout.count)
        &&actual.detail.truePeakEighthFrames.size()==size_t(expected.fileLayout.count)&&actual.detail.samplePeakFrames.size()==size_t(expected.fileLayout.count),"Incomplete candidate channel peaks");
    measuredContext=actual.context;snapshot=actual.identity;binding=r.binding();planIdentity=proposalPlanIdentity(actual.context.processing);
    profileDigest=r.currentContext.profile.digest();applicability=r.applicability;measurementMethod=actual.measurementVersion;measuredContent=actual.measuredContentHash;
    programme=measuredProgrammeForSnapshot(&actual,expected,MeasurementDomain::privateCandidate,applicability);
    require(programme.currentIdentityMatches,"Worker metric/context identity rejected");
    const auto& d=actual.detail;detail.sampleRate=d.sampleRate;detail.frames=d.frames;
    detail.integrated=d.integrated;detail.lra=d.lra;detail.momentaryMax=d.momentaryMax;detail.shortTermMax=d.shortTermMax;detail.samplePeak=d.samplePeak;detail.truePeak=d.truePeak;
    detail.integratedRelativeGate=d.integratedRelativeGate;detail.lraRelativeGate=d.lraRelativeGate;detail.p10=d.p10;detail.p95=d.p95;
    detail.integratedIncluded=d.integratedIncluded;detail.lraIncluded=d.lraIncluded;detail.lraStable=d.lraStable;
    detail.channelSamplePeaks=d.channelSamplePeaks;detail.channelTruePeaks=d.channelTruePeaks;detail.samplePeakFrames=d.samplePeakFrames;detail.truePeakEighthFrames=d.truePeakEighthFrames;
    id=hash("region-assessment-receipt-v1;"+binding+";"+planIdentity+";"+profileDigest+";"+programme.snapshotIdentity);
}
juce::String RegionAssessmentReceipt::validationError(const SectionProposalRequest& r,const ProcessingPlan& plan)const
{
    if(!r.regionSettings||!(proposal::isEventMethod(r.regionSettings->methodVersion())||(r.regionSettings->schema==1&&r.regionSettings->methodVersion()==proposal::regionSearchVersion))||r.binding()!=binding||r.pendingManualDraft)return "Stale region assessment request or pending manual draft";
    auto expected=r.currentContext;expected.processing=plan;
    if(proposalPlanIdentity(plan)!=planIdentity||r.currentContext.profile.digest()!=profileDigest||r.applicability!=applicability)return "Assessment belongs to another exact plan/requirement definition";
    if(!r.source||!r.current||!r.sourceContext.sourceFilesMatch()||!r.currentContext.sourceFilesMatch())return "Programme members changed after assessment";
    const auto scalar=r.scalar();if(!scalar.source.currentIdentityMatches||!scalar.current.currentIdentityMatches||scalar.sourceInputIdentity!=scalar.currentInputIdentity)return "Source and Current no longer match measured inputs";
    if(measurementMethod!=measurement::version||measuredContent!=expected.media.sha256||measuredContext.sessionId!=expected.sessionId
        ||OfflineAnalysisService::validateRequest(expected,MeasurementDomain::privateCandidate).isNotEmpty()
        ||snapshot!=OfflineAnalysisService::identityFor(expected,MeasurementDomain::privateCandidate,expected.media.sha256))return "Stale source, roles, range, build or processing receipt";
    return {};
}
proposal::Evaluation RegionAssessmentReceipt::evaluation(const SectionProposalRequest& r,const ProcessingPlan& plan)const
{
    auto input=programme;input.currentIdentityMatches=validationError(r,plan).isEmpty();
    return proposal::evaluateMeasuredProgramme(r.currentContext.profile,input);
}
size_t RegionAssessmentReceipt::storageBytes()const
{
    size_t n=sizeof(*this)+(detail.channelSamplePeaks.capacity()+detail.channelTruePeaks.capacity())*sizeof(double)
        +(detail.samplePeakFrames.capacity()+detail.truePeakEighthFrames.capacity())*sizeof(int64_t);
    for(const auto* s:{&snapshot,&id,&binding,&planIdentity,&profileDigest,&measurementMethod,&measuredContent})n+=s->getNumBytesAsUTF8()*4;
    const auto stringBytes=[](const juce::String& s){return 64+4*s.getNumBytesAsUTF8();};
    const auto mediaBytes=[&](const MediaInfo& m)
    {
        size_t bytes=0;for(const auto* s:{&m.format,&m.sha256,&m.layoutProvenance,&m.measurementRestriction})bytes+=stringBytes(*s);
        for(const auto& s:m.tags.getAllKeys())bytes+=sizeof(juce::String)*2+stringBytes(s);
        for(const auto& s:m.tags.getAllValues())bytes+=sizeof(juce::String)*2+stringBytes(s);
        return bytes;
    };
    const auto assemblyBytes=[&](const ProgrammeAssembly& a)
    {
        size_t bytes=a.members.capacity()*sizeof(ProgrammeMember)+stringBytes(a.alignmentConfirmation);
        for(const auto& member:a.members)bytes+=stringBytes(member.id)+stringBytes(member.file.getFullPathName())+mediaBytes(member.media);
        return bytes;
    };
    const auto& c=measuredContext;n+=mediaBytes(c.media);
    for(const auto* s:{&c.sessionId,&c.profileVersion,&c.buildVersion,&c.platformAndFloatingPointMode,&c.processing.graphJson,&c.processing.engineVersion})n+=stringBytes(*s);
    n+=stringBytes(c.sourceFile.getFullPathName())+stringBytes(proposalTargetText(c.profile));
    if(c.processing.sectionGains)n+=c.processing.sectionGains->intervals.capacity()*sizeof(processing::GainInterval)+stringBytes(c.processing.sectionGains->inputIdentity);
    if(c.assembly)n+=assemblyBytes(*c.assembly);
    if(c.composition)
    {
        n+=stringBytes(c.composition->confirmation);
        for(const auto* part:{&c.composition->me,&c.composition->dialogue})
        {
            for(const auto* s:{&part->id,&part->purpose,&part->language,&part->versionLabel})n+=stringBytes(*s);
            n+=mediaBytes(part->media)+stringBytes(part->file.getFullPathName());if(part->assembly)n+=assemblyBytes(*part->assembly);
        }
    }
    return n+16384; // fixed metric/provenance strings and allocator overhead
}
}

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
    if(!r.regionSettings||!proposal::isEventMethod(r.regionSettings->methodVersion())||r.binding()!=binding||r.pendingManualDraft)return "Stale event assessment request or pending manual draft";
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
    // Context/programme strings and shared input recipes are conservatively charged.
    return n+16384;
}
}

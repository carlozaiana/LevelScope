#include "StandaloneMarkerImport.h"
#include "SessionDocument.h"
namespace levelscope::standalone
{
bool preflightMarkerImport(const StandaloneSessionModel& s,const proposal::markers::Preview& p,juce::String& error)
{
    try
    {
        if(s.sessionId!=p.runtime||!s.sectionDocument.map().valid()||s.sectionDocument.map().digest()!=p.previousMap
           ||s.sourceContentInvalid||p.receipt.imported.inputIdentity!=sectionInputIdentity(s.measurementContext(MeasurementDomain::source)))
        {error="Stale marker preview: programme, annotation map or runtime changed. Preview again.";return false;}
        const auto checked=proposal::markers::decode(proposal::markers::encode(p.receipt));
        auto next=s;next.sectionTrimDraft.reset();next.sectionProposal.reset();
        if(!next.sectionDocument.restore(checked.imported,s.sectionDocument.map().inputIdentity))return false;
        next.markerImport=checked;auto document=sessionfile::capture(next);return sessionfile::coreFits(document,error);
    }
    catch(const std::exception& e){error=e.what();return false;}
}
bool commitMarkerImport(StandaloneSessionModel& s,const proposal::markers::Preview& p,const juce::String& hash,juce::String& error)
{
    if(s.audioEditsBlocked()||s.job.isActive()){error="Stop playback and finish/cancel the active job before confirming markers.";return false;}
    if(hash!=p.receipt.inputHash){error="Marker bytes changed since preview; preview again.";return false;}
    if(!s.sourceDocument.matchesFileStamps()){error="Programme source changed; annotations unchanged.";return false;}
    if(!preflightMarkerImport(s,p,error))return false;
    if(s.sectionDocument.map().digest()==p.receipt.imported.digest()&&s.markerImport&&s.markerImport->identity()==p.receipt.identity())
    {error="NO CHANGE: annotations, processing, measurements and revisions preserved.";return true;}
    if(s.documentRevision>=UINT64_MAX-1){error="Document revision exhausted";return false;}
    // All validation and allocation precede the single atomic annotation edit.
    auto receipt=p.receipt;proposal::SectionDocument document;
    if(!document.restore(receipt.imported,s.sectionDocument.map().inputIdentity)){error="Invalid marker replacement";return false;}
    s.sectionDocument=std::move(document);s.markerImport=std::move(receipt);s.documentChanged();
    error="Annotations replaced. Audio, applied absolute gain intervals and measurements unchanged. Geometry-dependent drafts/proposals require matching boundaries.";
    return true;
}
}

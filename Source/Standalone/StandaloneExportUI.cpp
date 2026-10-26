#include "StandaloneMainComponent.h"

namespace levelscope::standalone
{
namespace
{
// UI eligibility is a hint over accepted immutable session snapshots. Workers
// recompute complete snapshot bindings and both measured gates before publication.
bool exportUiReady(const StandaloneSessionModel& session,const exporting::Request& request)
{
    exporting::Specification specification;
    return request.mediaIssue.isEmpty()&&request.sourceAnalysis&&request.currentAnalysis
        &&exporting::specification(request,specification).isEmpty()
        &&session.evaluation(MeasurementDomain::currentState).verdict==proposal::Verdict::pass;
}
}

void StandaloneMainComponent::configureExportControls()
{
    for (auto* b : {&renderPcmButton,&publishPcmButton,&measurementReportButton,&artifactReportButton,&discardExportButton}) addAndMakeVisible(b);
    renderPcmButton.onClick=[this]{setPage(WorkflowPage::exportResults);chooseExportDestination(false);};
    measurementReportButton.onClick=[this]{setPage(WorkflowPage::exportResults);chooseExportDestination(true);};
    artifactReportButton.onClick=[this]{setPage(WorkflowPage::exportResults);chooseExportDestination(true,true);};
    publishPcmButton.onClick=[this]{publishExport();};
    discardExportButton.onClick=[this]
    {if(!session.job.isActive()){exportCandidate.reset();exportNotice="Rendered staging artifact discarded; Source/Current unchanged.";refreshFromSession();}};
}
void StandaloneMainComponent::refreshExportControls()
{
    const bool idle=canReplaceSession()&&!audition.connected();
    const auto request=exporting::capture(session);
    if(exportSourceCopy&&!exportSourceCopy->matches(request.current.media))exportSourceCopy.reset();
    if(exportCandidate&&!session.job.isActive()&&!exportCandidate->matchesAudio(request))
    {exportCandidate.reset();exportNotice="Previous rendered artifact is no longer current; no export approval retained.";}
    const bool ready=idle&&exportUiReady(session,request);
    renderPcmButton.setEnabled(ready);
    renderPcmButton.setTooltip(ready?"Render this Current State and measure the actual finalized PCM":idle?"Matching measured Source/Current and passing selected Current rules are required":"Stop playback, close the device and finish other jobs first.");
    measurementReportButton.setEnabled(idle);
    measurementReportButton.setTooltip("Source/Current measurements only, including failures or unavailable states. No exported audio or artifact verdict.");
    publishPcmButton.setEnabled(ready&&exportCandidate&&exporting::artifactEvaluation(*exportCandidate,request).verdict==proposal::Verdict::pass);
    discardExportButton.setEnabled(idle&&exportCandidate!=nullptr);
    artifactReportButton.setEnabled(idle&&exportCandidate!=nullptr);
}
void StandaloneMainComponent::chooseExportDestination(bool reportOnly,bool artifactOnly)
{
    if(!canReplaceSession()||audition.connected())return;
    if(!reportOnly&&!exportUiReady(session,exporting::capture(session)))return;
    const auto base=session.documentFile!=juce::File()?session.documentFile.getParentDirectory():juce::File::getSpecialLocation(juce::File::userDocumentsDirectory);
    exportChooser=std::make_unique<juce::FileChooser>(reportOnly?"New measurement-report folder (no audio)":"New PCM export folder (created only after explicit Publish)",
        base.getChildFile(reportOnly?"LevelScope-measurements":"LevelScope-export"));
    juce::Component::SafePointer<StandaloneMainComponent> safe(this);
    exportChooser->launchAsync(juce::FileBrowserComponent::saveMode|juce::FileBrowserComponent::canSelectFiles,
        [safe,reportOnly,artifactOnly](const auto& picker)
        {
            if(!safe||(artifactOnly&&!safe->exportCandidate)||picker.getResult()==juce::File()||!safe->canReplaceSession()||safe->audition.connected())return;
            const auto dest=picker.getResult();
            if(dest.exists()||dest.isSymbolicLink()){safe->exportNotice="Choose a new folder name; no existing destination will be overwritten.";safe->refreshFromSession();return;}
            safe->exportRequest=exporting::capture(safe->session);
            if(!reportOnly&&!exportUiReady(safe->session,*safe->exportRequest))return;
            if(!safe->session.job.begin(reportOnly?JobStage::report:JobStage::render,safe->session.revision))return;
            safe->exportTicket=safe->session.job.ticket;safe->session.persistenceLoading=true;
            safe->importControl=std::make_shared<JobControl>();safe->publishingExport=false;safe->savingArtifactReport=artifactOnly;
            if(reportOnly)
            {
                safe->session.job.message=artifactOnly?"Writing measured artifact report - NOT EXPORTED":"Writing measurement-only report; no audio or artifact verdict";
                if(artifactOnly)safe->exportWriteResult=std::async(std::launch::async,[r=*safe->exportRequest,dest,c=safe->importControl,candidate=safe->exportCandidate]{return exporting::artifactReport(*candidate,r,dest,*c);});
                else safe->exportWriteResult=std::async(std::launch::async,[r=*safe->exportRequest,dest,c=safe->importControl]{return exporting::measurementReport(r,dest,*c);});
            }
            else
            {
                safe->exportCandidate.reset();safe->exportDestination=dest;
                if(!safe->exportSourceCopy){const auto auditionCopy=safe->audition.sharedPreparedCopy();if(auditionCopy&&auditionCopy->matches(safe->exportRequest->current.media))safe->exportSourceCopy=auditionCopy;}
                safe->session.job.message="Rendering PCM24 and verifying the finalized output; nothing published";
                safe->renderResult=std::async(std::launch::async,[r=*safe->exportRequest,parent=dest.getParentDirectory(),c=safe->importControl,copy=safe->exportSourceCopy]{return exporting::render(r,parent,*c,{},nullptr,copy);});
            }
            safe->refreshFromSession();
        });
}
void StandaloneMainComponent::publishExport()
{
    if(!canReplaceSession()||audition.connected()||!exportCandidate)return;
    const auto fresh=exporting::capture(session);
    if(!exportUiReady(session,fresh)||exporting::artifactEvaluation(*exportCandidate,fresh).verdict!=proposal::Verdict::pass)return;
    if(!session.job.begin(JobStage::render,session.revision))return;
    exportTicket=session.job.ticket;exportRequest=fresh;session.persistenceLoading=true;publishingExport=true;
    importControl=std::make_shared<JobControl>();session.job.message="Rechecking exact bytes and publishing complete audio/report folder";
    exportWriteResult=std::async(std::launch::async,[candidate=exportCandidate,fresh,dest=exportDestination,c=importControl]
        {const auto approval=exporting::approvalIdentity(*candidate,fresh);return exporting::publish(*candidate,fresh,dest,approval,*c);});refreshFromSession();
}
void StandaloneMainComponent::pollExport()
{
    if(renderResult.valid()&&renderResult.wait_for(std::chrono::seconds(0))==std::future_status::ready)
    {
        auto result=renderResult.get();session.persistenceLoading=false;
        const bool current=exportRequest&&exportRequest->current.sessionId==session.sessionId
            &&exportRequest->documentRevision==session.documentRevision;
        if(session.job.finish(exportTicket,session.revision,current?result.error:ErrorCode::staleRevision,result.message)&&result.value)
        {
            exportCandidate=*result.value;exportSourceCopy=exportCandidate->verifiedSource();
            exportNotice="NOT PUBLISHED.\nSHA-256: "+exportCandidate->hash()
                +"\nDestination: "+exportDestination.getFullPathName();
            session.job.message="Rendered artifact measured. Review actual results and explicitly Publish if eligible.";
        }
        else exportNotice="No audio published: "+session.job.message+". Current State processing and measurements retained.";
        importControl.reset();exportRequest.reset();refreshFromSession();
    }
    if(exportWriteResult.valid()&&exportWriteResult.wait_for(std::chrono::seconds(0))==std::future_status::ready)
    {
        auto result=exportWriteResult.get();session.persistenceLoading=false;
        // The exclusive rename is a completed transaction. A later Cancel must
        // not falsely describe an already committed complete folder as absent.
        if(result.value&&session.job.state==JobState::cancelling)session.job.state=JobState::running;
        const bool current=exportRequest&&exportRequest->current.sessionId==session.sessionId&&exportRequest->documentRevision==session.documentRevision;
        const bool accepted=session.job.finish(exportTicket,session.revision,current?result.error:ErrorCode::staleRevision,result.message);
        if(accepted&&result.value)
        {
            exportNotice=(publishingExport?"Published verified PCM and reports: ":savingArtifactReport?"Saved artifact report - NOT EXPORTED: ":"Saved MEASUREMENT-ONLY report - NO EXPORTED AUDIO: ")+result.value->getFullPathName();
            if(publishingExport){exportNotice+="\nPCM SHA-256: "+exportCandidate->hash();exportCandidate.reset();}
        }
        else exportNotice="Publication/report failed: "+session.job.message+". No successful new set published.";
        importControl.reset();exportRequest.reset();refreshFromSession();
    }
}
juce::String StandaloneMainComponent::exportDetails() const
{
    const auto request=exporting::capture(session);
    juce::String text=exporting::preview(request,false)+"\n\n";
    // Applicability/profile changes reuse encoded audio, but its displayed verdict
    // must always describe the current selection rather than render completion.
    const bool currentArtifact=!session.job.isActive()&&exportCandidate&&exportCandidate->matchesAudio(request);
    if(currentArtifact)text+=evaluationText(exporting::artifactEvaluation(*exportCandidate,request))+"\n\n";
    text+=exportNotice+"\n\n";
    auto measuredLine=[&](const AnalysisSnapshot* snapshot,const ServiceContext& context,MeasurementDomain domain,const char* label)
    {
        const bool bound=snapshot&&snapshot->domain==domain&&snapshot->context.sessionId==context.sessionId
            &&snapshot->identity==OfflineAnalysisService::identityFor(context,domain,context.media.sha256);
        if(!bound||request.mediaIssue.isNotEmpty())return juce::String(label)+": measurement unavailable/stale\n";
        auto value=[](const measurement::Metric& m){return m.value&&std::isfinite(*m.value)?juce::String(*m.value,6):juce::String("unavailable");};
        return juce::String(label)+": I "+value(snapshot->detail.integrated)+" LUFS | TP "+value(snapshot->detail.truePeak)+" dBTP | LRA "+value(snapshot->detail.lra)+" LU\n";
    };
    text+=measuredLine(request.sourceAnalysis.get(),request.source,MeasurementDomain::source,"Source");
    text+=measuredLine(request.currentAnalysis.get(),request.current,MeasurementDomain::currentState,"Current State");
    if(currentArtifact)
    {
        text+=measuredLine(&exportCandidate->analysis(),exportCandidate->analysis().context,MeasurementDomain::renderedOutput,"Exported artifact");
    }
    text+="Save measurement report is independent: Source/Current only, including failures/unavailable states. It cannot publish PCM.\n"
          "Output is a new local folder with programme.wav, report.json, report.html and manifest.json. No complete broadcaster/BWF/BW64 delivery claim.";
    if(audition.connected())text+="\nClose the audio device before offline export/report work.";
    return text;
}
}

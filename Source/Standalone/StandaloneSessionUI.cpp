#include "StandaloneMainComponent.h"

namespace levelscope::standalone
{
bool StandaloneMainComponent::canReplaceSession() const
{ return !audition.active()&&!session.job.isActive()&&!session.persistenceLoading&&!documentResult.valid()&&!importResult.valid()&&!analysisResult.valid()&&!proposalResult.valid(); }
void StandaloneMainComponent::configureSessionControls()
{
    for(auto* c:std::initializer_list<juce::Component*>{&newSessionButton,&openSessionButton,&saveSessionButton,&saveAsSessionButton,&relinkAudioButton,&relinkVideoButton,&cancelSessionButton,&sessionStatus})addAndMakeVisible(*c);
    sessionStatus.setName("Session document status");
    newSessionButton.onClick=[this]
    {
        if(!canReplaceSession())return;
        confirmDiscard([this]{audition.disconnect();session.newDocument();audioPanel.sourceChanged();videoPanel.restoreSessionReference();documentNotice="New session; device closed and all routes muted";refreshFromSession();});
    };
    openSessionButton.onClick=[this]
    {
        if(!canReplaceSession())return;
        confirmDiscard([this]
        {
            sessionChooser=std::make_unique<juce::FileChooser>("Open LevelScope session",session.documentFile,"*.levelscope");
            juce::Component::SafePointer<StandaloneMainComponent> safe(this);
            sessionChooser->launchAsync(juce::FileBrowserComponent::openMode|juce::FileBrowserComponent::canSelectFiles,[safe](const auto& p){if(safe&&p.getResult()!=juce::File())safe->beginOpenSession(p.getResult());});
        });
    };
    saveSessionButton.onClick=[this]{beginSaveSession(false);};saveAsSessionButton.onClick=[this]{beginSaveSession(true);};
    relinkAudioButton.onClick=[this]{chooseRelink(false);};relinkVideoButton.onClick=[this]{chooseRelink(true);};
    cancelSessionButton.onClick=[this]{if(documentControl)documentControl->cancel.store(true);};
}
void StandaloneMainComponent::refreshSessionControls()
{
    const bool busy=documentResult.valid();const bool replace=canReplaceSession();
    newSessionButton.setEnabled(replace);openSessionButton.setEnabled(replace);saveSessionButton.setEnabled(!busy&&!session.persistenceLoading&&!audition.active());saveAsSessionButton.setEnabled(!busy&&!session.persistenceLoading&&!audition.active());
    relinkAudioButton.setEnabled(replace&&session.sourceDocument.hasSource);relinkVideoButton.setEnabled(replace&&session.referenceVideo.file!=juce::File());cancelSessionButton.setEnabled(busy);
    juce::String text=(session.documentFile==juce::File()?"Untitled":session.documentFile.getFileName())+(session.documentDirty()?" — UNSAVED":" — saved/unchanged");
    if(busy&&documentControl)text+=" | "+juce::String(fileOperation==FileOperation::save?"Saving snapshot":"Opening/verifying media")+" "+juce::String(juce::roundToInt(100*documentControl->progress.load()))+"%";
    else if(documentNotice.isNotEmpty())text+=" | "+documentNotice.replace("\n","; ");
    sessionStatus.setText(text,juce::dontSendNotification);sessionStatus.setTooltip(documentNotice);
}
void StandaloneMainComponent::confirmDiscard(std::function<void()> action)
{
    if(!session.documentDirty()){action();return;}
    juce::Component::SafePointer<StandaloneMainComponent> safe(this);
    juce::AlertWindow::showAsync(juce::MessageBoxOptions().withIconType(juce::MessageBoxIconType::QuestionIcon)
        .withTitle("Unsaved session").withMessage("Save your session changes before continuing?").withButton("Save").withButton("Discard").withButton("Cancel"),
        [safe,action](int result){if(!safe)return;if(result==1)safe->beginSaveSession(false,action);else if(result==2)action();});
}
void StandaloneMainComponent::requestClose(std::function<void()> close)
{
    if(documentResult.valid()||session.persistenceLoading){documentNotice="Cancel or finish the session file operation before closing";refreshSessionControls();return;}
    audition.stop();audition.disconnect();
    confirmDiscard(close);
}
void StandaloneMainComponent::beginSaveSession(bool saveAs,std::function<void()> continuation)
{
    if(documentResult.valid()||session.persistenceLoading||audition.active()){documentNotice="Stop audition and finish the file operation before saving";return;}
    if(saveAs||session.documentFile==juce::File())
    {
        sessionChooser=std::make_unique<juce::FileChooser>("Save LevelScope session",session.documentFile==juce::File()?juce::File::getSpecialLocation(juce::File::userDocumentsDirectory).getChildFile("Untitled.levelscope"):session.documentFile,"*.levelscope");
        juce::Component::SafePointer<StandaloneMainComponent> safe(this);
        sessionChooser->launchAsync(juce::FileBrowserComponent::saveMode|juce::FileBrowserComponent::canSelectFiles|juce::FileBrowserComponent::warnAboutOverwriting,[safe,continuation](const auto& picker)
        {
            if(!safe||picker.getResult()==juce::File())return;
            // Setting the destination is not a successful save; restore it until completion.
            auto original=safe->session.documentFile;safe->session.documentFile=picker.getResult();
            safe->beginSaveSession(false,continuation);safe->session.documentFile=original;
        });return;
    }
    try { capturedSave=sessionfile::capture(session); }
    catch(const std::exception& e){documentNotice="Cannot capture session: "+juce::String(e.what());refreshSessionControls();return;}
    pendingDocumentFile=session.documentFile;afterSave=std::move(continuation);fileOperation=FileOperation::save;documentControl=std::make_shared<JobControl>();
    documentResult=std::async(std::launch::async,[snapshot=*capturedSave,file=pendingDocumentFile,c=documentControl]{return sessionfile::save(snapshot,file,*c);});refreshSessionControls();
}
void StandaloneMainComponent::beginOpenSession(const juce::File& file)
{
    if(!canReplaceSession())return;
    audition.disconnect();session.persistenceLoading=true;pendingDocumentFile=file;fileOperation=FileOperation::open;documentControl=std::make_shared<JobControl>();documentNotice.clear();
    documentResult=std::async(std::launch::async,[file,c=documentControl]
    {
        auto r=sessionfile::read(file,*c);if(!r.document)return r;return sessionfile::verify(std::move(*r.document),*c);
    });refreshFromSession();
}
void StandaloneMainComponent::chooseRelink(bool video)
{
    if(!canReplaceSession())return;
    sessionChooser=std::make_unique<juce::FileChooser>(video?"Relink identical reference picture":"Relink identical source audio",juce::File(),video?"*.mov;*.mp4;*.m4v":getSourceFileWildcard());
    juce::Component::SafePointer<StandaloneMainComponent> safe(this);
    sessionChooser->launchAsync(juce::FileBrowserComponent::openMode|juce::FileBrowserComponent::canSelectFiles,[safe,video](const auto& picker)
    {
        if(!safe||picker.getResult()==juce::File()||!safe->canReplaceSession())return;
        auto d=sessionfile::capture(safe->session);if(video)sessionfile::relinkVideo(d,picker.getResult());else sessionfile::relinkAudio(d,picker.getResult());
        safe->audition.disconnect();safe->session.persistenceLoading=true;safe->relinkingVideo=video;safe->pendingDocumentFile=safe->session.documentFile;safe->fileOperation=FileOperation::relink;safe->documentControl=std::make_shared<JobControl>();
        safe->documentResult=std::async(std::launch::async,[d=std::move(d),c=safe->documentControl]{return sessionfile::verify(d,*c);});safe->refreshFromSession();
    });
}
void StandaloneMainComponent::installDocument(sessionfile::Document d,bool relinking,std::shared_ptr<sessionfile::Prepared> prepared)
{
    audition.disconnect();juce::String message;
    if(!prepared||!sessionfile::adopt(d,*prepared,session,message)){documentNotice=message;refreshSessionControls();return;}
    session.documentFile=pendingDocumentFile;if(relinking)session.documentChanged();
    audioPanel.sourceChanged();videoPanel.restoreSessionReference();documentNotice="Opened stopped, device closed, routes muted";
    if(message.isNotEmpty())documentNotice+="\n"+message;
    refreshFromSession();
}
void StandaloneMainComponent::pollSessionOperation()
{
    if(!documentResult.valid()||documentResult.wait_for(std::chrono::seconds(0))!=std::future_status::ready)return;
    sessionfile::Outcome result;
    try { result=documentResult.get(); }
    catch(const std::exception& e){result.error=ErrorCode::invalidMetadata;result.message=e.what();}
    const bool wasCancelled=documentControl&&documentControl->isCancelled();documentControl.reset();auto operation=fileOperation;fileOperation=FileOperation::none;session.persistenceLoading=false;
    if(operation==FileOperation::save)
    {
        // A successful atomic replacement is committed even if Cancel arrived just afterwards.
        if(result.succeeded()&&capturedSave&&sessionfile::saved(session,*capturedSave,pendingDocumentFile))
        {
            documentNotice=session.documentDirty()?"Captured revision saved; newer changes remain UNSAVED":"Session saved";
            auto action=std::move(afterSave);capturedSave.reset();if(action&&!session.documentDirty())action();
        }
        else{documentNotice=result.message;afterSave={};capturedSave.reset();}
    }
    else if(!result.document||wasCancelled)documentNotice=wasCancelled?"Session opening cancelled; previous session retained":result.message;
    else
    {
        auto d=std::move(*result.document);const bool unresolved=!d.audioResolved||!d.videoResolved;
        if(operation==FileOperation::relink&&!(relinkingVideo?d.videoResolved:d.audioResolved)){documentNotice=d.notices.joinIntoString("\n");}
        else if(unresolved && operation != FileOperation::relink)
        {
            session.persistenceLoading=true;
            juce::Component::SafePointer<StandaloneMainComponent> safe(this);
            juce::AlertWindow::showAsync(juce::MessageBoxOptions().withIconType(juce::MessageBoxIconType::WarningIcon).withTitle("Session media unresolved")
                .withMessage(d.notices.joinIntoString("\n")+"\nOpen settings and sections with unresolved media?").withButton("Open unresolved").withButton("Keep previous session"),
                [safe,d,prepared=result.prepared](int choice){if(!safe)return;safe->session.persistenceLoading=false;if(choice==1&&safe->canReplaceSession())safe->installDocument(d,false,prepared);safe->refreshFromSession();});
        }
        else installDocument(std::move(d),operation==FileOperation::relink,result.prepared);
    }
    refreshFromSession();
}
}

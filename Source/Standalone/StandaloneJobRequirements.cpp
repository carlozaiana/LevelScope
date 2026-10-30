#include "StandaloneJobRequirements.h"
#include <stdexcept>

namespace levelscope::standalone
{
std::optional<proposal::JobRequirements> loadRequirementsPreset(const juce::File& file,juce::String& error)
{
    try
    {
        const auto size=file.getSize();
        if(!file.existsAsFile()||size<1||size>proposal::JobRequirements::maximumBytes)throw std::runtime_error("Preset missing, empty or exceeds 64 KiB");
        auto stream=file.createInputStream();if(!stream||stream->failedToOpen())throw std::runtime_error("Preset unreadable");
        juce::MemoryBlock data;stream->readIntoMemoryBlock(data,size);
        if(data.getSize()!=size||!stream->isExhausted())throw std::runtime_error("Preset changed or could not be read completely");
        if(std::memchr(data.getData(),0,data.getSize())||!juce::CharPointer_UTF8::isValidString(static_cast<const char*>(data.getData()),int(data.getSize())))throw std::runtime_error("Preset must contain valid UTF-8 JSON without NUL bytes");
        auto result=proposal::JobRequirements::fromJson(juce::String::fromUTF8(static_cast<const char*>(data.getData()),int(data.getSize())));
        error.clear();return result;
    }
    catch(const std::exception& e){error=e.what();return {};}
}
bool saveRequirementsPreset(const proposal::JobRequirements& requirements,const juce::File& file,juce::String& error,
                            const std::function<bool(const char*)>& fail)
{
    error=requirements.validationError();if(error.isNotEmpty())return false;
    if(file.isSymbolicLink()||!file.getParentDirectory().isDirectory()){error="Select an existing folder and a regular preset destination.";return false;}
    const auto fault=[&](const char* stage){return fail&&fail(stage);};
    juce::TemporaryFile temporary(file);
    auto stream=fault("open")?nullptr:temporary.getFile().createOutputStream();
    if(!stream||stream->failedToOpen()){error="Cannot create temporary preset; previous file preserved.";return false;}
    const auto text=juce::JSON::toString(requirements.encode(),false)+"\n";
    if(text.getNumBytesAsUTF8()>proposal::JobRequirements::maximumBytes||fault("write")||!stream->write(text.toRawUTF8(),text.getNumBytesAsUTF8()))
    {error="Preset write failed; previous file preserved.";return false;}
    stream->flush();const bool ok=stream->getStatus().wasOk();stream.reset();
    if(!ok||fault("replace")||!temporary.overwriteTargetFileWithTemporary()){error="Preset replacement failed; previous file preserved.";return false;}
    error.clear();return true;
}
JobRequirementsPanel::JobRequirementsPanel(StandaloneSessionModel& s,std::function<void()> c):session(s),changed(std::move(c))
{
    units.setText("LU",juce::dontSendNotification);nameLabel.setText("Name",juce::dontSendNotification);
    provenanceLabel.setText("User requirement reference (optional)",juce::dontSendNotification);
    ceiling.setName("Draft maximum LRA in LU (0.0 to 100.0)");ceiling.setTextToShowWhenEmpty("Unset",juce::Colours::grey);
    ceiling.setInputRestrictions(6);name.setName("Job requirements name");name.setInputRestrictions(80);
    provenance.setName("User-specified requirement provenance");provenance.setInputRestrictions(2048);
    enabled.onClick=[this]{ceiling.setEnabled(enabled.getToggleState());message("Draft only. Use EBU requirements to confirm; audio unchanged.");};
    for(auto* field:{&ceiling,&name,&provenance})field->onTextChange=[this]{message("Draft only. Use EBU requirements to confirm; active requirements unchanged.");};
    for(auto* component:std::initializer_list<juce::Component*>{&enabled,&ceiling,&units,&nameLabel,&name,&provenanceLabel,&provenance,&use,&save,&load,&status})addAndMakeVisible(component);
    use.onClick=[this]{confirm();};save.onClick=[this]{choose(true);};load.onClick=[this]{choose(false);};refresh();
}
void JobRequirementsPanel::message(const juce::String& s){status.setText(s,juce::dontSendNotification);status.setTooltip(s);}
void JobRequirementsPanel::refresh()
{
    const bool available=!session.audioEditsBlocked()&&!session.job.isActive();
    for(auto* c:std::initializer_list<juce::Component*>{&enabled,&name,&provenance,&use,&load})c->setEnabled(available);
    ceiling.setEnabled(available&&enabled.getToggleState());save.setEnabled(available&&session.jobRequirements.validationError().isEmpty());
    if(runtime==session.sessionId&&activeDigest==session.jobRequirements.digest())return;
    runtime=session.sessionId;activeDigest=session.jobRequirements.digest();const auto& r=session.jobRequirements;
    enabled.setToggleState(bool(r.maximumLraTenths),juce::dontSendNotification);
    ceiling.setText(r.maximumLraTenths?juce::String(*r.maximumLraTenths/10.,1):juce::String(),false);
    name.setText(r.name,false);provenance.setText(r.provenance,false);ceiling.setEnabled(available&&bool(r.maximumLraTenths));
    message(r.validationError().isNotEmpty()?r.validationError():r.maximumLraTenths?"Active job ceiling: "+juce::String(*r.maximumLraTenths/10.,1)+" LU. User-specified; short LRA requires stable measurement.":"No LRA acceptance limit specified. Base EBU rules unchanged.");
}
void JobRequirementsPanel::confirm()
{
    if(runtime!=session.sessionId||activeDigest!=session.jobRequirements.digest()){refresh();message("Session changed; review current requirements.");return;}
    auto r=session.jobRequirements.validationError().isEmpty()?session.jobRequirements:proposal::JobRequirements{};
    r.maximumLraTenths.reset();
    if(enabled.getToggleState())
    {r.maximumLraTenths=proposal::JobRequirements::parseCeiling(ceiling.getText());if(!r.maximumLraTenths){message("Enter an explicit 0.0 to 100.0 LU value in exact tenths. Zero is a strict ceiling; no suggested value is supplied.");return;}}
    r.name=name.getText();r.provenance=provenance.getText();
    if(r.digest()!=session.jobRequirements.digest())
    {
        if(r.id==proposal::JobRequirements{}.id){r.id=juce::Uuid().toString();r.revision=1;}
        else ++r.revision;
    }
    juce::String error;if(!session.setJobRequirements(r,error)){message(error);return;}
    session.setSelectedTargetProfileIndex(1);if(changed)changed();activeDigest.clear();refresh();
}
void JobRequirementsPanel::choose(bool saving)
{
    if(session.audioEditsBlocked()||session.job.isActive())return;
    const auto opening=session.sessionId,binding=session.jobRequirements.digest();const auto captured=session.jobRequirements;
    chooser=std::make_unique<juce::FileChooser>(saving?"Save active requirements preset (draft edits excluded)":"Load named job requirements",juce::File{},"*.levelscope-target");
    chooser->launchAsync((saving?juce::FileBrowserComponent::saveMode|juce::FileBrowserComponent::warnAboutOverwriting:juce::FileBrowserComponent::openMode)|juce::FileBrowserComponent::canSelectFiles,
        [safe=juce::Component::SafePointer<JobRequirementsPanel>(this),saving,opening,binding,captured](const juce::FileChooser& fc)
        {
            if(!safe)return;auto file=fc.getResult();if(file==juce::File()){safe->message("Preset operation cancelled; active definition unchanged.");return;}
            auto& self=*safe;juce::String error;
            if(opening!=self.session.sessionId||binding!=self.session.jobRequirements.digest()||self.session.audioEditsBlocked()||self.session.job.isActive()){self.message("Session changed; preset operation discarded.");return;}
            if(saving)
            {
                if(!file.hasFileExtension("levelscope-target"))file=file.withFileExtension("levelscope-target");
                self.message(saveRequirementsPreset(captured,file,error)?"Active requirements saved; session definition unchanged.":error);
            }
            else
            {
                const auto r=loadRequirementsPreset(file,error);
                if(!r||!self.session.setJobRequirements(*r,error)){self.message(error);return;}
                self.session.setSelectedTargetProfileIndex(1);if(self.changed)self.changed();self.activeDigest.clear();self.refresh();
            }
        });
}
void JobRequirementsPanel::resized()
{
    auto r=getLocalBounds();auto a=r.removeFromTop(30);enabled.setBounds(a.removeFromLeft(260));ceiling.setBounds(a.removeFromLeft(75).reduced(2));units.setBounds(a.removeFromLeft(28));nameLabel.setBounds(a.removeFromLeft(52));use.setBounds(a.removeFromRight(175).reduced(2));name.setBounds(a.reduced(3));
    auto b=r.removeFromTop(30);provenanceLabel.setBounds(b.removeFromLeft(260));load.setBounds(b.removeFromRight(125).reduced(2));save.setBounds(b.removeFromRight(125).reduced(2));provenance.setBounds(b.reduced(3));status.setBounds(r);
}
}

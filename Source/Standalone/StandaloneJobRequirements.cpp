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

namespace
{
juce::String targetText(std::optional<int> value)
{return value?juce::String(*value/10.,1):juce::String("No Limit");}
juce::String measuredText(const measurement::Metric& m)
{return m.status==measurement::Status::available && m.value?juce::String(*m.value,1):juce::String("--");}
}
JobRequirementsPanel::JobRequirementsPanel(StandaloneSessionModel& s,std::function<void()> c,juce::File folder)
    :session(s),changed(std::move(c)),presetDirectory(folder==juce::File{}?
      juce::File::getSpecialLocation(juce::File::userApplicationDataDirectory).getChildFile("LevelScope/Target Presets"):folder)
{
    const char* names[]={"IL | LUFS", "Max Peak | dBTP", "LRA | LU", "Max M | LUFS", "Max S | LUFS"};
    for(size_t n=0;n<5;++n)
    {
        headers[n].setText(names[n],juce::dontSendNotification);headers[n].setJustificationType(juce::Justification::centred);
        addAndMakeVisible(headers[n]);auto& f=fields[n];f.setName(names[n]);f.setFont(16.f);f.setJustification(juce::Justification::centred);
        f.setInputRestrictions(32);f.setTooltip(n==0?"IL target and tolerance, e.g. -23.0 +/- 0.2, or No Limit. Enter to apply.":"Exact tenths, or No Limit. Enter to apply.");
        f.onTextChange=[this]{if(updating)return;editing=true;profiles.setText("Unsaved preset",juce::dontSendNotification);message("Editing Target -- Enter or leave the field to apply.");};
        f.onReturnKey=[this]{commitFields();};f.onFocusLost=[this]{if(editing)commitFields();};
        f.onEscapeKey=[this]{editing=false;activeDigest.clear();refresh();};addAndMakeVisible(f);
    }
    for(auto& label:measured){label.setFont(16.f);label.setJustificationType(juce::Justification::centred);addAndMakeVisible(label);}
    for(auto& label:rowNames){label.setFont(16.f);addAndMakeVisible(label);}
    profiles.setName("Target profile");profiles.onChange=[this]{if(!updating)selectProfile();};addAndMakeVisible(profiles);
    save.onClick=[this]{savePreset();};addAndMakeVisible(save);addAndMakeVisible(status);
    rebuildProfiles();refresh();
}
void JobRequirementsPanel::message(const juce::String& text){status.setText(text,juce::dontSendNotification);status.setTooltip(text);}
void JobRequirementsPanel::rebuildProfiles()
{
    updating=true;profiles.clear(juce::dontSendNotification);profiles.addItem("No target profile",1);profiles.addItem("EBU R128",2);
    profiles.addItem("ATSC A/85 -- planned",3);profiles.setItemEnabled(3,false);
    profiles.addItem("Unsaved preset",4);profiles.setItemEnabled(4,false);profiles.addSeparator();presets.clear();
    const auto files=presetDirectory.findChildFiles(juce::File::findFiles,false,"*.levelscope-target");
    for(int n=0;n<juce::jmin(files.size(),256);++n)
    {juce::String error;auto r=loadRequirementsPreset(files[n],error);if(r){profiles.addItem(r->name,100+int(presets.size()));presets.push_back(*r);}}
    profiles.addSeparator();profiles.addItem("Load preset...",10);updating=false;
}
void JobRequirementsPanel::refreshAvailability()
{
    const bool available=!session.audioEditsBlocked()&&!session.job.isActive();
    profiles.setEnabled(available);save.setEnabled(available);
    for(auto& f:fields)f.setEnabled(available);
    const auto reason=session.audioEditsBlocked()?"Stop playback/preparation or finish loading to edit Target.":session.job.isActive()?"Finish or cancel the active job to edit Target.":"Choose a target profile; changing rules invalidates dependent drafts, not audio measurements.";
    profiles.setTooltip(reason);save.setTooltip(available?"Save all Target values, including No Limit, as a named preset.":reason);
}
void JobRequirementsPanel::refresh()
{
    refreshAvailability();
    for(int row=0;row<2;++row)
    {
        const auto& summary=row==0?session.source:session.currentState;
        const auto domain=row==0?MeasurementDomain::source:MeasurementDomain::currentState;
        const auto result=session.evaluation(domain);
        rowNames[size_t(row)].setText(juce::String(row==0?"Original":"Current")+"  |  "+(summary.hasMeasurement?"Ready":row==1&&session.currentStateNeedsRemeasurement()?"Unverified":"Pending")+(result.verdict==proposal::Verdict::pass?" / PASS":result.verdict==proposal::Verdict::fail?" / FAIL":""),juce::dontSendNotification);
        auto origin=summary.statusText;
        if(summary.analysis)
        {
            const auto& a=*summary.analysis;
            origin += a.identicalSourceMeasurementIdentity.isNotEmpty()?"\nInitialised from identical Source measurement "+a.identicalSourceMeasurementIdentity:"\nMeasured directly. "+a.identity;
            origin += "\nOrigin job "+juce::String(a.originSessionId.isEmpty()?a.context.ticket.id:a.originJobId)+"; runtime "+(a.originSessionId.isEmpty()?a.context.sessionId:a.originSessionId);
        }
        rowNames[size_t(row)].setTooltip(origin+"\n"+evaluationText(result));
        for(size_t col=0;col<5;++col)
        {
            auto& label=measured[size_t(row)*5+col];juce::String text="--";
            if(summary.analysis)
            {
                const auto& d=summary.analysis->detail;
                const measurement::Metric* metrics[]={&d.integrated,&d.truePeak,&d.lra,&d.momentaryMax,&d.shortTermMax};
                text=measuredText(*metrics[col]);if(col==2&&d.lra.value&&!d.lraStable)text+=" *";
            }
            label.setText(text,juce::dontSendNotification);label.setTooltip(origin+"\n"+evaluationText(result));
            const char* ids[]={"programme-loudness","true-peak","loudness-range","maximum-momentary","maximum-short-term"};
            auto colour=juce::Colour(0xffdce4ed);
            for(const auto& rule:result.rules)if(rule.id==ids[col])
            {if(rule.verdict==proposal::Verdict::fail)colour=juce::Colour(0xffffb27d);else if(rule.verdict==proposal::Verdict::pass)colour=juce::Colour(0xff88d8bd);}
            label.setColour(juce::Label::textColourId,colour);
        }
    }
    rowNames[2].setText("Target",juce::dontSendNotification);
    if(runtime==session.sessionId&&activeDigest==session.jobRequirements.digest()&&selectedIndex==session.selectedTargetProfileIndex)return;
    runtime=session.sessionId;activeDigest=session.jobRequirements.digest();selectedIndex=session.selectedTargetProfileIndex;
    auto r=proposal::JobRequirements::editable(session.getSelectedTargetProfile());
    updating=true;editing=false;
    fields[0].setText(r.integratedTenths?targetText(r.integratedTenths)+" +/- "+juce::String(r.integratedToleranceTenths/10.,1):"No Limit",false);
    const std::optional<int> values[]={r.truePeakTenths,r.maximumLraTenths,r.maximumMomentaryTenths,r.maximumShortTermTenths};
    for(size_t n=1;n<5;++n)fields[n].setText(targetText(values[n-1]),false);
    profiles.setSelectedId(selectedIndex==0?1:selectedIndex==1&&session.jobRequirements.schema==1&&!session.jobRequirements.maximumLraTenths?2:4,juce::dontSendNotification);
    if(profiles.getSelectedId()==4)profiles.setText(session.jobRequirements.name==proposal::JobRequirements{}.name?"Unsaved preset":session.jobRequirements.name,juce::dontSendNotification);
    updating=false;message(session.jobRequirements.validationError().isNotEmpty()?session.jobRequirements.validationError():"Edit Target values directly. No Limit disables that rule. IL includes its tolerance. * LRA under 60 s is unstable.");
}
bool JobRequirementsPanel::commitFields()
{
    if(!editing)return true;
    if(runtime!=session.sessionId||activeDigest!=session.jobRequirements.digest()||selectedIndex!=session.selectedTargetProfileIndex)
    {activeDigest.clear();refresh();message("Session changed; review the current Target.");return false;}
    auto r=session.jobRequirements.schema==2?session.jobRequirements:proposal::JobRequirements::editable(session.getSelectedTargetProfile());
    const auto read=[&](juce::String text,std::optional<int>& value,bool lra=false)
    {
        if(text.trim().equalsIgnoreCase("No Limit")){value.reset();return true;}
        value=lra?proposal::JobRequirements::parseCeiling(text):proposal::JobRequirements::parseSignedTenths(text);return bool(value);
    };
    auto il=fields[0].getText().replace(juce::String::charToString(0x00b1),"+/-");const int split=il.indexOf("+/-");
    if(!read(split<0?il:il.substring(0,split),r.integratedTenths)) {message("Enter IL in exact tenths or No Limit.");return false;}
    if(split>=0)
    {auto tolerance=proposal::JobRequirements::parseCeiling(il.substring(split+3));if(!tolerance){message("Enter an IL tolerance from 0.0 to 100.0 LU.");return false;}r.integratedToleranceTenths=*tolerance;}
    if(!read(fields[1].getText(),r.truePeakTenths)||!read(fields[2].getText(),r.maximumLraTenths,true)
        ||!read(fields[3].getText(),r.maximumMomentaryTenths)||!read(fields[4].getText(),r.maximumShortTermTenths))
    {message("Enter a value in exact tenths or No Limit. LRA range is 0.0-100.0 LU; other limits are -100.0-24.0.");return false;}
    // Leaving an unchanged field must preserve a saved profile and its identity.
    if(session.jobRequirements.schema==2&&r.digest()==session.jobRequirements.digest())
    {editing=false;activeDigest.clear();refresh();return true;}
    r.schema=2;r.name="Unsaved preset";++r.revision;
    juce::String error;if(!session.setJobRequirements(r,error)){message(error);return false;}
    session.setSelectedTargetProfileIndex(3);editing=false;activeDigest.clear();refresh();if(changed)changed();return true;
}
void JobRequirementsPanel::selectProfile()
{
    const int id=profiles.getSelectedId();if(id==10){loadPreset();return;}
    if(id==4)return;
    juce::String error;proposal::JobRequirements r;
    if(id>=100&&id<100+int(presets.size()))r=presets[size_t(id-100)];
    else if(id!=1&&id!=2)return;
    if(!session.setJobRequirements(r,error)){message(error);return;}
    session.setSelectedTargetProfileIndex(id==1?0:id==2?1:r.schema==2?3:1);
    editing=false;activeDigest.clear();refresh();if(changed)changed();
}
void JobRequirementsPanel::savePreset()
{
    if(saveDialog||!commitFields()||session.audioEditsBlocked()||session.job.isActive())return;
    auto captured=proposal::JobRequirements::editable(session.getSelectedTargetProfile());
    captured.name=session.jobRequirements.schema==2?session.jobRequirements.name:profiles.getText();
    captured.provenance=session.jobRequirements.provenance;
    const auto opening=session.sessionId,binding=session.jobRequirements.digest();const int index=session.selectedTargetProfileIndex;
    saveDialog=std::make_unique<juce::AlertWindow>("Save target preset",captured.text(),juce::MessageBoxIconType::NoIcon);
    saveDialog->addTextEditor("name",captured.name=="Unsaved preset"?juce::String():captured.name,"Preset name");
    saveDialog->addTextEditor("reference",captured.provenance,"Delivery reference (optional)");
    saveDialog->addButton("Save",1,juce::KeyPress(juce::KeyPress::returnKey));saveDialog->addButton("Cancel",0,juce::KeyPress(juce::KeyPress::escapeKey));
    saveDialog->enterModalState(true,juce::ModalCallbackFunction::create([safe=juce::Component::SafePointer<JobRequirementsPanel>(this),captured,opening,binding,index](int result) mutable
    {
        if(!safe)return;
        auto& self=*safe;
        auto dialog=std::move(self.saveDialog);
        dialog->setVisible(false);
        if(result!=1)return;
        if(opening!=self.session.sessionId||binding!=self.session.jobRequirements.digest()||index!=self.session.selectedTargetProfileIndex||self.session.audioEditsBlocked()||self.session.job.isActive())
        {self.message("Session changed; preset save cancelled.");return;}
        captured.name=dialog->getTextEditorContents("name").trim();captured.provenance=dialog->getTextEditorContents("reference");
        juce::String error=captured.validationError();if(error.isNotEmpty()){self.message(error);return;}
        if(self.presetDirectory.createDirectory().failed()){self.message("Cannot create the preset folder.");return;}
        if(!saveRequirementsPreset(captured,self.presetDirectory.getChildFile(captured.id+".levelscope-target"),error)){self.message(error);return;}
        if(!self.session.setJobRequirements(captured,error)){self.message(error);return;}
        self.session.setSelectedTargetProfileIndex(3);self.rebuildProfiles();self.activeDigest.clear();self.refresh();if(self.changed)self.changed();
    }),false);
}
void JobRequirementsPanel::loadPreset()
{
    const auto opening=session.sessionId,binding=session.jobRequirements.digest();const int index=session.selectedTargetProfileIndex;
    chooser=std::make_unique<juce::FileChooser>("Load target preset",juce::File{},"*.levelscope-target");
    chooser->launchAsync(juce::FileBrowserComponent::openMode|juce::FileBrowserComponent::canSelectFiles,
        [safe=juce::Component::SafePointer<JobRequirementsPanel>(this),opening,binding,index](const juce::FileChooser& fc)
    {
        if(!safe)return;auto& self=*safe;auto file=fc.getResult();self.activeDigest.clear();self.refresh();if(file==juce::File{})return;
        if(opening!=self.session.sessionId||binding!=self.session.jobRequirements.digest()||index!=self.session.selectedTargetProfileIndex){self.message("Session changed; preset load cancelled.");return;}
        juce::String error;const auto r=loadRequirementsPreset(file,error);
        if(!r||!self.session.setJobRequirements(*r,error)){self.message(error);return;}
        self.session.setSelectedTargetProfileIndex(r->schema==2?3:1);self.activeDigest.clear();self.refresh();if(self.changed)self.changed();
    });
}
int JobRequirementsPanel::rowHeight() const
{
    const int available=getHeight()-(status.getParentComponent()==this?22:0);
    return juce::jlimit(22,54,(available-22-(available<112?0:24))/3);
}
void JobRequirementsPanel::resized()
{
    auto r=getLocalBounds();const int labelWidth=230,colWidth=(r.getWidth()-labelWidth)/5;
    const bool compact=getHeight()-(status.getParentComponent()==this?22:0)<112;
    auto heading=r.removeFromTop(22);auto leading=heading.removeFromLeft(labelWidth);
    for(auto& h:headers)h.setBounds(heading.removeFromLeft(colWidth));
    const int height=rowHeight();
    for(size_t row=0;row<3;++row)
    {
        auto line=r.removeFromTop(height);rowNames[row].setBounds(line.removeFromLeft(labelWidth));
        for(size_t col=0;col<5;++col)
        {auto cell=line.removeFromLeft(colWidth).reduced(4,std::min(3,std::max(0,(height-22)/2)));if(row==2)fields[col].setBounds(cell);else measured[row*5+col].setBounds(cell);}
    }
    auto choices=compact?leading:r.removeFromTop(24);
    if(!compact)choices.removeFromLeft(labelWidth);
    profiles.setBounds(choices.removeFromLeft(compact?123:280).reduced(2,0));save.setBounds(choices.removeFromLeft(compact?107:160).reduced(2,0));
    if(status.getParentComponent()==this)status.setBounds(r);
}
void JobRequirementsPanel::paint(juce::Graphics& g)
{
    g.fillAll(juce::Colour(0xff202730));g.setColour(juce::Colour(0xff394451));
    for(int n=0;n<=3;++n)g.drawHorizontalLine(22+n*rowHeight(),0.f,float(getWidth()));
}
}

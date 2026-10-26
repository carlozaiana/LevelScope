#include "StandaloneCompositionImport.h"

namespace levelscope::standalone
{
namespace
{
class GroupWindow final:public juce::DocumentWindow
{
public:
    explicit GroupWindow(ProgrammeImportPanel& p):DocumentWindow("Component channel files",juce::Colour(0xff181b20),closeButton),panel(p)
    {setUsingNativeTitleBar(true);setContentNonOwned(&p,true);centreWithSize(p.getWidth(),p.getHeight());setVisible(true);}
    void closeButtonPressed() override{panel.cancel();}
private:ProgrammeImportPanel& panel;
};
}
CompositionImportPanel::CompositionImportPanel(std::optional<ProgrammeComposition> existing,Completion done):complete(std::move(done))
{
    if(existing)recipe=*existing;
    else{recipe.me.id="me";recipe.me.purpose="music-effects";recipe.dialogue.id="dialogue";recipe.dialogue.purpose="dialogue";}
    confirmation.setToggleState(existing&&recipe.validate().isEmpty(),juce::dontSendNotification);
    preset.addItem("Stereo M&E + stereo dialogue",1);preset.addItem("5.1 M&E + mono dialogue to Centre",2);preset.addItem("5.1 M&E + 5.1 dialogue",3);
    preset.setSelectedId(static_cast<int>(recipe.preset)+1,juce::dontSendNotification);
    preset.onChange=[this]{recipe.preset=static_cast<CompositionPreset>(preset.getSelectedId()-1);changed();refresh();};
    for(auto* c:std::initializer_list<juce::Component*>{&preset,&explanation,&notice,&accept,&cancelButton,&confirmation})addAndMakeVisible(*c);
    explanation.setText("Original Source = confirmed M&E + dialogue at unity gain. Current State applies the explicit global trim.\nMatching rates/lengths do not prove alignment or balance. Changing composition audio replaces programme-bound measurements and section annotations; save your previous session first.",juce::dontSendNotification);
    for(int k=0;k<2;++k)
    {
        choose[k].setButtonText(k==0?"Choose M&E file...":"Choose dialogue file...");groups[k].setButtonText("Use channel files...");canonical[k].setButtonText("Assign stated file order");
        languages[k].setText(part(k).language,false);versions[k].setText(part(k).versionLabel,false);
        languages[k].setTextToShowWhenEmpty("Language label (optional)",juce::Colours::grey);versions[k].setTextToShowWhenEmpty("Version label (optional)",juce::Colours::grey);
        languages[k].setInputRestrictions(128);versions[k].setInputRestrictions(256);
        languages[k].onTextChange=[this,k]{part(k).language=languages[k].getText();};versions[k].onTextChange=[this,k]{part(k).versionLabel=versions[k].getText();};
        for(auto* c:std::initializer_list<juce::Component*>{&choose[k],&groups[k],&canonical[k],&descriptions[k],&languages[k],&versions[k]})addAndMakeVisible(*c);
        choose[k].onClick=[this,k]{pickFile(k);};groups[k].onClick=[this,k]{pickGroup(k);};
        canonical[k].onClick=[this,k]
        {
            auto& p=part(k);const auto expected=expectedRoles(k);if(p.assembly||p.media.channels!=static_cast<int>(expected.size()))return;
            std::copy(expected.begin(),expected.end(),p.layout.roles.begin());p.layout.manuallyInterpreted=true;changed();refresh();
        };
        for(size_t channel=0;channel<roles[k].size();++channel)
        {
            addAndMakeVisible(channelNames[k][channel]);addAndMakeVisible(roles[k][channel]);roles[k][channel].setTextWhenNothingSelected("Confirm channel role");
            roles[k][channel].onChange=[this,k,channel]{part(k).layout.roles[channel]=static_cast<ChannelRole>(roles[k][channel].getSelectedId());part(k).layout.manuallyInterpreted=true;changed();refresh();};
        }
    }
    confirmation.onClick=[this]{refresh();};
    accept.onClick=[this]
    {
        if(!confirmation.getToggleState()||pending.valid()||groupPanel)return;
        auto candidate=recipe;candidate.me.layout.confirmed=candidate.dialogue.layout.confirmed=true;
        if(!candidate.confirmComplementaryBalancedAlignment()||!candidate.matchesFileStamps())
        {notice.setText("Composition incomplete or changed; inspect files, roles and confirmation again",juce::dontSendNotification);return;}
        complete(std::move(candidate));
    };
    cancelButton.onClick=[this]{cancel();};setSize(1000,700);refresh();startTimerHz(10);
}
CompositionImportPanel::~CompositionImportPanel()
{stopTimer();if(control)control->cancel.store(true);if(pending.valid())pending.wait();groupWindow.reset();groupPanel.reset();}
std::vector<ChannelRole> CompositionImportPanel::expectedRoles(int k) const
{return k==1&&recipe.preset==CompositionPreset::surroundMonoCentre?std::vector<ChannelRole>{juce::AudioChannelSet::centre}:assemblyRoles(recipe.preset==CompositionPreset::stereoStereo?2:6);}
void CompositionImportPanel::changed(){recipe.confirmation.clear();operationNotice.clear();confirmation.setToggleState(false,juce::dontSendNotification);}
void CompositionImportPanel::pickFile(int k)
{
    if(pending.valid()||groupPanel)return;
    chooser=std::make_unique<juce::FileChooser>(k==0?"Select M&E component":"Select dialogue component",juce::File(),"*.wav;*.wave;*.aif;*.aiff;*.flac");
    juce::Component::SafePointer<CompositionImportPanel> safe(this);
    chooser->launchAsync(juce::FileBrowserComponent::openMode|juce::FileBrowserComponent::canSelectFiles,[safe,k](const auto& p)
    {
        if(!safe||p.getResult()==juce::File())return;
        safe->pendingFile=p.getResult();safe->pendingComponent=k;safe->control=std::make_shared<JobControl>();
        safe->pending=std::async(std::launch::async,[file=safe->pendingFile,c=safe->control]{return validateMedia(file,*c);});safe->changed();safe->refresh();
    });
}
void CompositionImportPanel::pickGroup(int k)
{
    if(pending.valid()||groupPanel||expectedRoles(k).size()==1)return;
    juce::Component::SafePointer<CompositionImportPanel> safe(this);
    groupPanel=std::make_unique<ProgrammeImportPanel>([safe,k](std::optional<ProgrammeAssembly> group)
    {
        juce::MessageManager::callAsync([safe,k,group=std::move(group)]
        {
            if(!safe)return;safe->groupWindow.reset();safe->groupPanel.reset();
            if(safe->closing){safe->complete({});return;}
            if(group){auto& p=safe->part(k);p.file=juce::File();p.assembly=*group;p.media=group->summary();p.layout=group->layout();safe->changed();}
            safe->refresh();
        });
    });
    groupWindow=std::make_unique<GroupWindow>(*groupPanel);refresh();
}
void CompositionImportPanel::cancel()
{closing=true;if(control)control->cancel.store(true);if(groupPanel){groupPanel->cancel();return;}if(!pending.valid())complete({});else refresh();}
void CompositionImportPanel::timerCallback()
{
    if(!pending.valid())return;
    if(pending.wait_for(std::chrono::seconds(0))!=std::future_status::ready)
    {notice.setText(closing?"Cancelling; previous session retained":"Verifying component: "+juce::String(juce::roundToInt(100*control->progress.load()))+"%",juce::dontSendNotification);return;}
    MediaValidationResult result;
    try{result=pending.get();}catch(const std::exception& e){result.error=ErrorCode::decodeFailed;result.message=e.what();}
    control.reset();if(closing){complete({});return;}
    if(result.succeeded())
    {
        auto& p=part(pendingComponent);p.file=pendingFile;p.assembly.reset();p.media=result.media;p.layout=result.media.suggestedLayout;
        p.layout.confirmed=false;
    }
    operationNotice=result.succeeded()?juce::String():result.message+"; previous component retained";refresh();
}
void CompositionImportPanel::refresh()
{
    const bool busy=pending.valid()||groupPanel;const bool editable=!busy&&!closing;
    preset.setEnabled(editable);confirmation.setEnabled(editable);
    auto probe=recipe;probe.me.layout.confirmed=probe.dialogue.layout.confirmed=true;
    const auto error=probe.validate(false);accept.setEnabled(editable&&error.isEmpty()&&confirmation.getToggleState());
    if(!busy)notice.setText(operationNotice.isNotEmpty()?operationNotice:error.isEmpty()?"Review each component, role mapping and confirmation. No component normalization is applied.":error,juce::dontSendNotification);
    for(int k=0;k<2;++k)
    {
        auto& p=part(k);const auto expected=expectedRoles(k);
        choose[k].setEnabled(editable);groups[k].setEnabled(editable&&expected.size()>1);canonical[k].setEnabled(editable&&!p.assembly&&p.media.channels==static_cast<int>(expected.size()));
        juce::String order;for(auto role:expected){if(order.isNotEmpty())order+=" ";order+=juce::AudioChannelSet::getAbbreviatedChannelTypeName(role);}
        canonical[k].setButtonText("Assign file order: "+order);
        languages[k].setEnabled(editable);versions[k].setEnabled(editable);
        juce::String description=k==0?"M&E":"Dialogue";
        if(p.media.sha256.isNotEmpty())
        {description+=" | "+(p.assembly?juce::String(p.assembly->members.size())+" channel files":p.file.getFileName())+"\n"+juce::String(p.media.sampleRate,0)+" Hz | "+juce::String(p.media.frames)+" frames | unity gain";}
        else description+=" | No component selected";
        descriptions[k].setText(description,juce::dontSendNotification);
        for(size_t n=0;n<roles[k].size();++n)
        {
            auto& box=roles[k][n];box.clear(juce::dontSendNotification);for(auto role:expected)box.addItem(juce::AudioChannelSet::getChannelTypeName(role),static_cast<int>(role));
            const auto name=p.assembly&&n<p.assembly->members.size()?p.assembly->members[n].file.getFileName():"Channel "+juce::String(n+1);
            channelNames[k][n].setText(name,juce::dontSendNotification);channelNames[k][n].setTooltip(name);channelNames[k][n].setVisible(n<expected.size());
            box.setVisible(n<expected.size());box.setEnabled(editable&&!p.assembly&&n<static_cast<size_t>(p.media.channels));
            box.setSelectedId(static_cast<int>(p.layout.roles[n]),juce::dontSendNotification);
            box.setTooltip(p.assembly&&n<p.assembly->members.size()?p.assembly->members[n].file.getFileName():"File Channel "+juce::String(n+1));
        }
    }
}
void CompositionImportPanel::resized()
{
    auto b=getLocalBounds().reduced(12);preset.setBounds(b.removeFromTop(32).removeFromLeft(500));explanation.setBounds(b.removeFromTop(64));
    auto componentArea=b.removeFromTop(350);
    for(int k=0;k<2;++k)
    {
        auto column=componentArea.removeFromLeft(componentArea.getWidth()/(2-k)).reduced(5);auto row=column.removeFromTop(30);
        choose[k].setBounds(row.removeFromLeft(235));groups[k].setBounds(row);descriptions[k].setBounds(column.removeFromTop(45));canonical[k].setBounds(column.removeFromTop(30));
        for(size_t n=0;n<roles[k].size();++n){auto line=column.removeFromTop(29).reduced(1);channelNames[k][n].setBounds(line.removeFromLeft(190));roles[k][n].setBounds(line);}
        row=column.removeFromTop(30);languages[k].setBounds(row.removeFromLeft(row.getWidth()/2).reduced(1));versions[k].setBounds(row.reduced(1));
    }
    confirmation.setBounds(b.removeFromTop(48));notice.setBounds(b.removeFromTop(35));auto buttons=b.removeFromTop(32);accept.setBounds(buttons.removeFromLeft(290));cancelButton.setBounds(buttons.removeFromLeft(120));
}
}

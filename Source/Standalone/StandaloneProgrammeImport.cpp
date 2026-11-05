#include "StandaloneProgrammeImport.h"

namespace levelscope::standalone
{
ProgrammeImportPanel::ProgrammeImportPanel(Completion callback,juce::File roleDirectory,bool analyse):complete(std::move(callback)),roleMemory(roleDirectory)
{
    if(!analyse)accept.setButtonText("Confirm component roles");
    layout.addItem("Stereo - two mono files",2);layout.addItem("5.1 - six mono files",6);layout.setSelectedId(2,juce::dontSendNotification);
    layout.onChange=[this]{group.reset();alignment.setToggleState(false,juce::dontSendNotification);refresh();};
    for(auto* c:std::initializer_list<juce::Component*>{&layout,&choose,&accept,&cancelButton,&alignment,&notice,&explanation,&suggestion,&editOrder})addAndMakeVisible(*c);
    explanation.setText("Assign each file to one channel. No resampling, padding, trimming or mixing.\nMetadata/timecode offsets are not applied; matching lengths do not prove alignment.",juce::dontSendNotification);
    for(size_t i=0;i<roles.size();++i)
    {
        addAndMakeVisible(labels[i]);addAndMakeVisible(roles[i]);roles[i].onChange=[this,i]
        {
            if(group&&i<group->members.size())group->members[i].role=static_cast<ChannelRole>(roles[i].getSelectedId());
            alignment.setToggleState(false,juce::dontSendNotification);if(group)group->alignmentConfirmation.clear();refresh();
        };
    }
    editOrder.onClick=[this]{editingOrder=true;refresh();};
    choose.onClick=[this]
    {
        chooser=std::make_unique<juce::FileChooser>("Select all mono files for this programme",juce::File(),"*.wav;*.wave;*.aif;*.aiff;*.flac");
        juce::Component::SafePointer<ProgrammeImportPanel> safe(this);
        chooser->launchAsync(juce::FileBrowserComponent::openMode|juce::FileBrowserComponent::canSelectFiles|juce::FileBrowserComponent::canSelectMultipleItems,[safe](const auto& p)
        {
            if(safe)safe->inspectFiles(p.getResults());
        });
    };
    alignment.onClick=[this]{refresh();};
    accept.onClick=[this]
    {
        if(closing || !group || !alignment.getToggleState())return;
        auto ordered=*group;
        if(!ordered.orderByRoles()||!ordered.confirmAlignment()||!ordered.matchesFileStamps()){notice.setText("Group changed or roles incomplete; inspect and confirm again",juce::dontSendNotification);return;}
        closing=true;setEnabled(false);
        complete(std::move(ordered));
    };
    cancelButton.onClick=[this]{cancel();};setSize(820,560);refresh();startTimerHz(10);
}
void ProgrammeImportPanel::inspectFiles(const juce::Array<juce::File>& selection)
{
    if(selection.isEmpty() || pending.valid() || closing)return;
    if(selection.size()!=layout.getSelectedId())
    {
        notice.setText("Select exactly "+juce::String(layout.getSelectedId())+" mono files",juce::dontSendNotification);
        return;
    }
    std::vector<juce::File> files(selection.begin(),selection.end());
    group.reset();alignment.setToggleState(false,juce::dontSendNotification);
    control=std::make_shared<JobControl>();control->feedback=std::make_shared<JobFeedback>();
    pending=std::async(std::launch::async,[files,c=control]
    {
        try{return inspectProgrammeFiles(files,*c);}
        catch(const std::exception& e){return ProgrammeInspection{{},ErrorCode::decodeFailed,e.what()};}
    });
    refresh();
}
ProgrammeImportPanel::~ProgrammeImportPanel(){stopTimer();if(control)control->cancel.store(true);if(pending.valid())pending.wait();}
std::shared_ptr<JobControl> ProgrammeImportPanel::preparationControl() const { return pending.valid() ? control : nullptr; }
void ProgrammeImportPanel::cancel(){if(closing)return;closing=true;if(control)control->cancel.store(true);if(!pending.valid())complete({});else refresh();}
void ProgrammeImportPanel::timerCallback()
{
    if(!pending.valid())return;
    if(pending.wait_for(std::chrono::seconds(0))!=std::future_status::ready)
    {notice.setText(closing?"Cancelling; previous session retained":"Verifying every member: "+juce::String(juce::roundToInt(control->progress.load()*100))+"%",juce::dontSendNotification);return;}
    auto result=pending.get();control.reset();
    if(closing){complete({});return;}
    if(result.error==ErrorCode::none)
    {
        group=std::move(result.group);
        editingOrder=false;
        const auto proposal=roleMemory.suggest(*group);
        for(size_t i=0;i<group->members.size();++i)group->members[i].role=proposal.order.roles[i];
        suggestion.setText(proposal.explanation,juce::dontSendNotification);
        suggestion.setTooltip(proposal.explanation);
    }
    notice.setText(result.message,juce::dontSendNotification);refresh();
}
ChannelLayout ProgrammeImportPanel::confirmedFileOrder() const
{
    ChannelLayout order;
    if(group && alignment.getToggleState())
    {
        order.count=int(group->members.size());
        for(int i=0;i<order.count;++i)order.roles[i]=group->members[size_t(i)].role;
        order.confirmed=order.isValid();
    }
    return order;
}
void ProgrammeImportPanel::refresh()
{
    const bool busy=pending.valid();choose.setEnabled(!busy);layout.setEnabled(!busy);alignment.setEnabled(group.has_value()&&!busy);
    auto proposed=group;const bool valid=proposed&&proposed->orderByRoles()&&proposed->validate(false).isEmpty();
    accept.setEnabled(valid&&alignment.getToggleState()&&!busy);
    editOrder.setEnabled(group.has_value() && !busy);
    suggestion.setVisible(group.has_value());
    const auto available=assemblyRoles(layout.getSelectedId());
    for(size_t i=0;i<roles.size();++i)
    {
        const bool visible=i<available.size();roles[i].setVisible(visible);labels[i].setVisible(visible);roles[i].setEnabled(group.has_value()&&!busy&&editingOrder);
        const int selected=group&&i<group->members.size()?static_cast<int>(group->members[i].role):0;
        roles[i].clear(juce::dontSendNotification);roles[i].setTextWhenNothingSelected("Assign channel role");
        for(auto role:available)roles[i].addItem(juce::AudioChannelSet::getChannelTypeName(role),static_cast<int>(role));
        roles[i].setSelectedId(selected,juce::dontSendNotification);
        juce::String text="No file selected";
        if(group&&i<group->members.size()){const auto& m=group->members[i];text=m.file.getFileName()+" | "+m.media.format+" "+juce::String(m.media.bitsPerSample)+" bit | "+juce::String(m.media.sampleRate,0)+" Hz | "+juce::String(m.media.frames)+" frames";if(!editingOrder && m.media.hasFileRoleLabels)text+=" | File label";
            labels[i].setTooltip(m.file.getFullPathName());}
        labels[i].setText(text,juce::dontSendNotification);
    }
}
void ProgrammeImportPanel::resized()
{
    auto b=getLocalBounds().reduced(14);auto row=b.removeFromTop(30);layout.setBounds(row.removeFromLeft(260));choose.setBounds(row.removeFromLeft(200));
    explanation.setBounds(b.removeFromTop(48));
    auto review=b.removeFromTop(40);editOrder.setBounds(review.removeFromRight(100).reduced(2));suggestion.setBounds(review);
    for(size_t i=0;i<roles.size();++i){auto line=b.removeFromTop(43);roles[i].setBounds(line.removeFromRight(190).reduced(2));labels[i].setBounds(line);}
    alignment.setBounds(b.removeFromTop(32));notice.setBounds(b.removeFromTop(42));auto buttons=b.removeFromTop(30);accept.setBounds(buttons.removeFromLeft(250));cancelButton.setBounds(buttons.removeFromLeft(120));
}
}

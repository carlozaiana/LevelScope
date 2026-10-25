#include "StandaloneProgrammeImport.h"

namespace levelscope::standalone
{
ProgrammeImportPanel::ProgrammeImportPanel(Completion callback):complete(std::move(callback))
{
    layout.addItem("Stereo - two mono files",2);layout.addItem("5.1 - six mono files",6);layout.setSelectedId(2);
    layout.onChange=[this]{group.reset();alignment.setToggleState(false,juce::dontSendNotification);refresh();};
    for(auto* c:std::initializer_list<juce::Component*>{&layout,&choose,&accept,&cancelButton,&alignment,&notice,&explanation})addAndMakeVisible(*c);
    explanation.setText("Assign each file to one channel. No resampling, padding, trimming or mixing.\nMetadata/timecode offsets are not applied; matching lengths do not prove alignment.",juce::dontSendNotification);
    for(size_t i=0;i<roles.size();++i)
    {
        addAndMakeVisible(labels[i]);addAndMakeVisible(roles[i]);roles[i].onChange=[this,i]
        {
            if(group&&i<group->members.size())group->members[i].role=static_cast<ChannelRole>(roles[i].getSelectedId());
            alignment.setToggleState(false,juce::dontSendNotification);if(group)group->alignmentConfirmation.clear();refresh();
        };
    }
    choose.onClick=[this]
    {
        chooser=std::make_unique<juce::FileChooser>("Select all mono files for this programme",juce::File(),"*.wav;*.wave;*.aif;*.aiff;*.flac");
        juce::Component::SafePointer<ProgrammeImportPanel> safe(this);
        chooser->launchAsync(juce::FileBrowserComponent::openMode|juce::FileBrowserComponent::canSelectFiles|juce::FileBrowserComponent::canSelectMultipleItems,[safe](const auto& p)
        {
            if(!safe||p.getResults().isEmpty())return;
            if(p.getResults().size()!=safe->layout.getSelectedId()){safe->notice.setText("Select exactly "+juce::String(safe->layout.getSelectedId())+" mono files",juce::dontSendNotification);return;}
            std::vector<juce::File> files;for(const auto& f:p.getResults())files.push_back(f);
            safe->group.reset();safe->alignment.setToggleState(false,juce::dontSendNotification);safe->control=std::make_shared<JobControl>();
            safe->pending=std::async(std::launch::async,[files,c=safe->control]{try{return inspectProgrammeFiles(files,*c);}catch(const std::exception& e){return ProgrammeInspection{{},ErrorCode::decodeFailed,e.what()};}});
            safe->refresh();
        });
    };
    alignment.onClick=[this]{refresh();};
    accept.onClick=[this]
    {
        if(!group||!alignment.getToggleState())return;
        auto ordered=*group;
        if(!ordered.orderByRoles()||!ordered.confirmAlignment()||!ordered.matchesFileStamps()){notice.setText("Group changed or roles incomplete; inspect and confirm again",juce::dontSendNotification);return;}
        complete(std::move(ordered));
    };
    cancelButton.onClick=[this]{cancel();};setSize(820,520);refresh();startTimerHz(10);
}
ProgrammeImportPanel::~ProgrammeImportPanel(){stopTimer();if(control)control->cancel.store(true);if(pending.valid())pending.wait();}
void ProgrammeImportPanel::cancel(){closing=true;if(control)control->cancel.store(true);if(!pending.valid())complete({});else refresh();}
void ProgrammeImportPanel::timerCallback()
{
    if(!pending.valid())return;
    if(pending.wait_for(std::chrono::seconds(0))!=std::future_status::ready)
    {notice.setText(closing?"Cancelling; previous session retained":"Verifying every member: "+juce::String(juce::roundToInt(control->progress.load()*100))+"%",juce::dontSendNotification);return;}
    auto result=pending.get();control.reset();
    if(closing){complete({});return;}
    if(result.error==ErrorCode::none)group=std::move(result.group);
    notice.setText(result.message,juce::dontSendNotification);refresh();
}
void ProgrammeImportPanel::refresh()
{
    const bool busy=pending.valid();choose.setEnabled(!busy);layout.setEnabled(!busy);alignment.setEnabled(group.has_value()&&!busy);
    auto proposed=group;const bool valid=proposed&&proposed->orderByRoles()&&proposed->validate(false).isEmpty();
    accept.setEnabled(valid&&alignment.getToggleState()&&!busy);
    const auto available=assemblyRoles(layout.getSelectedId());
    for(size_t i=0;i<roles.size();++i)
    {
        const bool visible=i<available.size();roles[i].setVisible(visible);labels[i].setVisible(visible);roles[i].setEnabled(group.has_value()&&!busy);
        const int selected=group&&i<group->members.size()?static_cast<int>(group->members[i].role):0;
        roles[i].clear(juce::dontSendNotification);roles[i].setTextWhenNothingSelected("Assign channel role");
        for(auto role:available)roles[i].addItem(juce::AudioChannelSet::getChannelTypeName(role),static_cast<int>(role));
        roles[i].setSelectedId(selected,juce::dontSendNotification);
        juce::String text="No file selected";
        if(group&&i<group->members.size()){const auto& m=group->members[i];text=m.file.getFileName()+" | "+m.media.format+" "+juce::String(m.media.bitsPerSample)+" bit | "+juce::String(m.media.sampleRate,0)+" Hz | "+juce::String(m.media.frames)+" frames";labels[i].setTooltip(m.file.getFullPathName());}
        labels[i].setText(text,juce::dontSendNotification);
    }
}
void ProgrammeImportPanel::resized()
{
    auto b=getLocalBounds().reduced(14);auto row=b.removeFromTop(30);layout.setBounds(row.removeFromLeft(260));choose.setBounds(row.removeFromLeft(200));
    explanation.setBounds(b.removeFromTop(58));
    for(size_t i=0;i<roles.size();++i){auto line=b.removeFromTop(43);roles[i].setBounds(line.removeFromRight(190).reduced(2));labels[i].setBounds(line);}
    alignment.setBounds(b.removeFromTop(32));notice.setBounds(b.removeFromTop(42));auto buttons=b.removeFromTop(30);accept.setBounds(buttons.removeFromLeft(250));cancelButton.setBounds(buttons.removeFromLeft(120));
}
}

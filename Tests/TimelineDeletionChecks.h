void deletionChecks(const juce::File& root)
{
    auto req=fixture(root.getChildFile("delete.wav"),48000,2,240013);StandaloneSessionModel s;
    expect(s.setValidatedSource(req.context.sourceFile,req.context.media),"deletion fixture imported");s.sourceDocument.layout=req.context.fileLayout;s.initializeCurrentStateFromSource();s.synchronizeSections();
    juce::String error;auto& doc=s.sectionDocument;doc.split(doc.map().sections.front().id,96000,error);doc.split(doc.map().sections.back().id,192000,error);
    const auto first=doc.map().sections[0].id,middle=doc.map().sections[1].id,last=doc.map().sections[2].id;
    s.setCurrentGainDb(4);expect(s.setManualSectionGains({first},-4)&&s.setManualSectionGains({middle},-6)&&s.setManualSectionGains({last},-2),"different neighbouring applied gains");
    const auto baseline=s.manualState();TimelineGainLane lane(s);MeasurementTimeline ui;ui.setSize(1150,640);ui.attachGainLane(lane);ui.setProgramme(240013,48000,"delete");lane.synchronize("delete");
    ui.addToDesktop(juce::ComponentPeer::windowHasTitleBar);ui.setVisible(true);ui.focusForTransport();
    const auto select=[&]{lane.synchronize("delete");lane.selectAt(120000,lane.plotBounds().getCentreY());ui.focusForTransport();};
    const auto answer=[&](int choice){expect(lane.editPrompt().window()!=nullptr,"different neighbour gains offer explicit choices");lane.editPrompt().window()->exitModalState(choice);const auto start=juce::Time::getMillisecondCounterHiRes();while(lane.editPrompt().active()&&juce::Time::getMillisecondCounterHiRes()-start<1500)juce::MessageManager::getInstance()->runDispatchLoopUntil(10);expect(!lane.editPrompt().active(),"choice closes after asynchronous callback");};
    select();expect(ui.keyPressed(juce::KeyPress(juce::KeyPress::backspaceKey)),"Backspace routes from timeline");imageOf(*lane.editPrompt().window(),root.getChildFile("delete-neighbour-choices.png"));answer(0);expect(s.manualState()==baseline,"Cancel preserves sections and gain");
    select();ui.keyPressed(juce::KeyPress(juce::KeyPress::deleteKey));answer(1);
    expect(s.currentStateDocument.sectionGains->gainAt(120000)==0&&s.currentStateDocument.staticGainDb==4&&doc.map().sections.size()==3&&!s.sectionTrimDraft,"Set to zero applies immediately, retaining necessary boundaries and global gain");
    expect(s.undoManual()&&s.manualState()==baseline,"Undo restores reset plus original annotation IDs");
    select();lane.deleteSelected();answer(2);expect(doc.map().sections.size()==2&&doc.map().sections[0].end==192000&&s.currentStateDocument.sectionGains->gainAt(120000)==-4,"Merge left immediately extends left gain through selected passage");s.undoManual();
    select();lane.deleteSelected();answer(3);expect(doc.map().sections.size()==2&&doc.map().sections[1].begin==96000&&s.currentStateDocument.sectionGains->gainAt(120000)==-2,"Merge right preserves left neighbour and extends right gain");s.undoManual();
    expect(s.setManualSectionGains({first},-2),"equal nonzero neighbours prepared");select();lane.deleteSelected();
    expect(button(*lane.editPrompt().window(),"Merge into neighbours (-2.00 dB)")&&!button(*lane.editPrompt().window(),"Merge right (-2.00 dB)"),"equal neighbours offer one merge choice");answer(2);expect(doc.map().sections.size()==1&&s.currentStateDocument.sectionGains->gainAt(120000)==-2,"merge into equal neighbours removes both redundant boundaries");s.undoManual();
    s.setManualSectionGains({first,last},0);select();lane.deleteSelected();
    expect(!lane.editPrompt().active()&&!s.currentStateDocument.sectionGains&&doc.map().sections.size()==1&&s.currentStateDocument.staticGainDb==4,"zero neighbours reset and merge in one action without a prompt");const auto noOpCount=s.manualHistory->undoCount();select();lane.deleteSelected();expect(s.manualHistory->undoCount()==noOpCount,"deleting the remaining unity programme creates no spurious Undo entry");expect(s.undoManual(),"direct reset is undoable");
    select();laneEditor(lane)->grabKeyboardFocus();expect(!ui.keyPressed(juce::KeyPress(juce::KeyPress::backspaceKey)),"Backspace in numeric editor never deletes sections");
    select();s.setManualSectionGains({first},-4);select();lane.deleteSelected();doc.rename(doc.map().sections.front().id,"New annotation",error);const auto changed=doc.map().digest();answer(2);expect(doc.map().digest()==changed,"stale modal answer cannot overwrite a newer annotation");
    ui.removeFromDesktop();
}

#include "StandaloneAudioPanel.h"
#include "StandaloneMeasurementLayout.h"

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-AUDIO-PANEL-IMPLEMENTATION]
StandaloneAudioPanel::StandaloneAudioPanel (StandaloneSessionModel& model, StandaloneAudition& engine,
                                           std::function<void()> changed,juce::File roleDirectory,
                                           std::shared_ptr<OutputDeviceCatalog> catalog, OutputDeviceSettings settings)
    : session (model), audio (engine), onChanged (std::move (changed)),roleMemory(roleDirectory),
      monitoring(settings.file()),connection(engine,std::move(catalog),std::move(settings))
{
    session.auditionBusy = [&engine] { return engine.active() && engine.getState() != AuditionState::cancelling; };
    for (auto* component : std::initializer_list<juce::Component*> {
        &deviceBox, &automaticReconnect, &scanButton, &openButton, &closeButton, &playButton,
        &stopButton, &restartButton, &seekButton, &seekPosition, &confirmButton, &gain, &statusLabel, &heading, &viewport, &domainBox, &prepareButton, &auditionLabel, &layoutPreset,&suggestionLabel,&editOrder,&monitoringStatus,&rememberRoutingButton })
        addAndMakeVisible (*component);
    monitoringStatus.setFont(juce::Font(13.0f));monitoringStatus.setJustificationType(juce::Justification::topLeft);
    monitoringStatus.setComponentID("monitoring-memory-status");
    rememberRoutingButton.setTooltip("Remember the displayed role assignments for this Mac, interface and layout, including an all-muted configuration.");
    rememberRoutingButton.onClick=[this]{rememberRouting();refresh();};
    deviceBox.setName("Output interface");
    deviceBox.setTextWhenNothingSelected ("Choose an output interface (offline analysis is available)");
    automaticReconnect.setToggleState(connection.automaticReconnect(),juce::dontSendNotification);
    automaticReconnect.setTooltip("Off by default. Reconnect this Mac's last successful output when the programme rate is known. Playback stays stopped. Disconnect pauses automatic reconnection until LevelScope restarts.");
    automaticReconnect.onClick = [this]
    {
        if(session.persistenceLoading || session.job.isActive() || audio.active())return;
        connection.setAutomaticReconnect(automaticReconnect.getToggleState());
        connectionChanged(); refresh(); onChanged();
    };
    scanButton.onClick = [this] { scan(); };
    deviceBox.onChange = [this]
    {
        if (session.persistenceLoading || session.job.isActive() || audio.active()) return;
        const auto index=deviceBox.getSelectedId()-1;
        if(juce::isPositiveAndBelow(index,int(connection.choices().size())))
        {
            const auto& choice=connection.choices()[size_t(index)].identity;
            if(connection.selection() && choice==*connection.selection())return; // Ignore a queued notification after internal restoration.
            connection.select(choice);
        }
        actionError.clear();connectionChanged();refresh();onChanged();
    };
    openButton.onClick = [this]
    { if(session.persistenceLoading || session.job.isActive() || audio.active())return;connection.connect();actionError.clear();connectionChanged();refresh();onChanged(); };
    closeButton.onClick = [this] { connection.disconnect();actionError.clear();connectionChanged();refresh();onChanged(); };
    domainBox.addItem ("Source", 1); domainBox.addItem ("Current State", 2);
    domainBox.setSelectedId (1, juce::dontSendNotification);
    domainBox.onChange = [this]
    {
        session.selectAuditionDomain (domainBox.getSelectedId() == 2 ? MeasurementDomain::currentState : MeasurementDomain::source);
        actionError.clear(); synchronize(); onChanged();
    };
    prepareButton.onClick = [this]
    {
        if (session.persistenceLoading) return;
        synchronize(); actionError = audio.prepare (session.auditionRequest()); refresh(); onChanged();
    };
    playButton.onClick = [this]
    {
        if (session.persistenceLoading) return;
        synchronize(); actionError = audio.play (session.auditionRequest(), session.sourceDocument.routing);
        refresh(); onChanged();
    };
    stopButton.onClick = [this] { audio.stop(); actionError.clear(); refresh(); onChanged(); };
    restartButton.onClick = [this]
    {
        if (session.persistenceLoading) return;
        synchronize();
        actionError = audio.play (session.auditionRequest(), session.sourceDocument.routing, true); refresh(); onChanged();
    };
    seekPosition.setText ("0"); seekPosition.setTooltip ("Audio timeline in seconds or an exact fraction; seek rounds down to a source sample.");
    seekButton.onClick = [this]
    {
        if (session.persistenceLoading) return;
        try
        {
            const auto& s = session.sourceDocument;
            const auto frame = video::Time::parse (seekPosition.getText().trim().toStdString()).floorTicks (static_cast<int64_t> (s.media.sampleRate));
            actionError = audio.seek (frame, session.auditionRequest(), s.routing);
        }
        catch (const std::exception& e) { actionError = e.what(); }
        refresh();
    };
    seekPosition.onReturnKey = [this] { seekButton.triggerClick(); };
    editOrder.onClick=[this]
    {
        editingOrder=true;refresh();
        if(roleBoxes[0].isShowing())roleBoxes[0].grabKeyboardFocus();
    };
    confirmButton.onClick=[this]
    {
        if(audio.active() || session.job.isActive() || !session.sourceDocument.matchesFileStamps() || !session.confirmLayout())return;
        const auto saved=roleMemory.remember(session.sourceDocument.layout);
        if(saved.failed())actionError=saved.getErrorMessage();
        const bool analyse=reviewing;
        if(analyse)session.sourceDocument.importStatus=editingOrder?"Explicitly confirmed edited file order":suggestionLabel.getText();
        const auto analysisError=OfflineAnalysisService::validateRequest(session.measurementContext(MeasurementDomain::source),MeasurementDomain::source);
        if(analyse && analysisError.isNotEmpty())actionError=analysisError;
        reviewing=false;
        suggestionLabel.setText("Confirmed channel order for this programme.",juce::dontSendNotification);
        refresh();resized();onChanged();
        if(analyse && analysisError.isEmpty() && onAnalyse)onAnalyse();
    };
    layoutPreset.setTextWhenNothingSelected ("Explicit System I order preset (then confirm)");
    layoutPreset.addItem ("L R C LFE Lrs Rrs Lss Rss — WAVE / FLAC order", 1);
    layoutPreset.addItem ("L R C LFE Lss Rss Lrs Rrs — side before rear", 2);
    layoutPreset.onChange = [this]
    {
        if (session.sourceDocument.layout.count != 8 || layoutPreset.getSelectedId() == 0) return;
        const auto layout = systemILayout (layoutPreset.getSelectedId() == 2);
        for (int i = 0; i < layout.count; ++i) session.setChannelRole (i, layout.roles[i]);
        layoutPreset.setSelectedId (0, juce::dontSendNotification);
        onChanged();
    };
    gain.setRange (-60, 0, 0.5);
    gain.setSliderStyle (juce::Slider::LinearHorizontal);
    gain.setTextBoxStyle (juce::Slider::TextBoxRight, false, 80, 24);
    gain.setTextValueSuffix (" dB");
    gain.setTooltip ("Monitor gain only. Adjust smoothly during playback. Processing and measurements are unaffected.");
    gain.setValue (-30, juce::dontSendNotification);
    gain.explicitChange = [this]
    {
        audio.setMonitorGainDb(static_cast<float>(gain.getValue()));
        // Only explicit slider/text/keyboard notifications persist a user value.
        // Initialisation and every internal restoration use dontSendNotification.
        if(const auto& id=connection.selection())monitoring.rememberVolume(*id,gain.getValue());
    };
    heading.setText ("File channel / Speaker role / Device output (no downmix)", juce::dontSendNotification);
    heading.setFont (juce::Font (12.0f));
    statusLabel.setFont (juce::Font (13.0f));
    statusLabel.setJustificationType (juce::Justification::topLeft);
    viewport.setViewedComponent (&rows, false);
    for (int channel = 0; channel < maxChannels; ++channel)
    {
        rows.addAndMakeVisible (channelLabels[channel]);
        rows.addAndMakeVisible (roleBoxes[channel]);
        rows.addAndMakeVisible (outputBoxes[channel]);
        channelLabels[channel].setText ("Channel " + juce::String (channel + 1), juce::dontSendNotification);
        roleBoxes[channel].setName("Role for file channel "+juce::String(channel+1));
        roleBoxes[channel].addItem ("Unassigned", 1);
        for (size_t role = 0; role < supportedRoles.size(); ++role)
            roleBoxes[channel].addItem (juce::AudioChannelSet::getChannelTypeName (supportedRoles[role]), static_cast<int> (role) + 2);
        roleBoxes[channel].onChange = [this, channel]
        {
            const auto index = roleBoxes[channel].getSelectedId() - 2;
            session.setChannelRole (channel, juce::isPositiveAndBelow (index, static_cast<int> (supportedRoles.size()))
                ? supportedRoles[static_cast<size_t> (index)] : session.sourceDocument.layout.count == 8
                ? systemIRoleForAlias (roleBoxes[channel].getText()) : juce::AudioChannelSet::unknown);
            refresh();onChanged();
        };
        outputBoxes[channel].setName("Output for file channel "+juce::String(channel+1));
        outputBoxes[channel].onChange = [this, channel]
        {
            if(audio.active()||session.job.isActive()||session.persistenceLoading||!session.sourceDocument.layout.confirmed)return;
            const auto output=outputBoxes[channel].getSelectedId()-2;
            if(output==session.sourceDocument.routing.fileToDevice[size_t(channel)])return; // Internal/queued notification, not a new route edit.
            if(session.setOutputRoute(channel,output))rememberRouting();
            actionError.clear();refresh();
        };
    }
    rebuildOutputs();
    refreshDeviceChoices();
    sourceChanged();
}

void StandaloneAudioPanel::scan()
{
    if(session.persistenceLoading || session.job.isActive() || audio.active())return;
    connection.rescan();actionError.clear();connectionChanged();refresh();onChanged();
}

void StandaloneAudioPanel::refreshDeviceChoices()
{
    deviceBox.clear(juce::dontSendNotification);
    int chosen=0;
    const auto& wanted=connection.selection();const auto& entries=connection.choices();
    for(size_t i=0;i<entries.size();++i)
    {
        const auto& item=entries[i];
        deviceBox.addItem(item.identity.name+" ("+item.identity.type+")",int(i)+1);
        if(wanted && item.identity==*wanted)chosen=int(i)+1;
    }
    if(wanted && chosen==0)
    {
        chosen=int(entries.size())+1;
        deviceBox.addItem(wanted->name+" (selected; unavailable)",chosen);
        deviceBox.setItemEnabled(chosen,false);
    }
    deviceBox.setSelectedId(chosen,juce::dontSendNotification);
    deviceBox.setTooltip(connection.status());
}

void StandaloneAudioPanel::connectionChanged()
{
    const bool changed=connectionRevision!=connection.deviceRevision();
    connectionRevision=connection.deviceRevision();
    refreshMonitoring(true,changed);
    refreshDeviceChoices();
}

void StandaloneAudioPanel::refreshMonitoring(bool discover,bool resetVolume)
{
    const auto& selected=connection.selection();
    const auto volumeKey=selected?MonitoringSettings::identityKey(*selected):juce::String();
    if(volumeDevice!=volumeKey || resetVolume)
    {
        volumeDevice=volumeKey;
        const auto saved=selected?monitoring.volume(*selected):std::optional<double>{};
        gain.setValue(saved.value_or(-30),juce::dontSendNotification);
        audio.setMonitorGainDb(static_cast<float>(gain.getValue()));
    }
    auto outputs=monitoringOutputs;
    if(discover)
        outputs=audio.connected()?audio.outputConfiguration():selected?connection.outputConfiguration():audio.outputConfiguration();
    else if(audio.connected())outputs=audio.outputConfiguration();
    const auto identity=connection.monitoringIdentity();const auto& source=session.sourceDocument;
    juce::String key=volumeKey+"|"+(identity?"identified":"unavailable")+"|";
    key<<int(source.hasSource)<<int(source.layout.confirmed)<<int(session.sourceContentInvalid)<<int(session.persistenceLoading)<<"|";
    for(int n=0;n<source.layout.count;++n)key<<int(source.layout.roles[n])<<",";
    for(const auto& name:outputs)key<<"|"<<name.length()<<":"<<name;
    if(routingContext==key)return; // Routine refresh never rewrites the user's live adjustments.
    if(audio.active())audio.stop(); // Restore/reconfigure only while transport is stopped.
    routingContext=key;monitoringOutputs=outputs;
    session.sourceDocument.routing=OutputRouting();
    if(!source.hasSource || !source.layout.confirmed || session.sourceContentInvalid || session.persistenceLoading)
        routingNotice="Confirm this programme's speaker roles before restoring output routing.";
    else if(!identity)
        routingNotice="No uniquely identified interface available; routes muted. Select or scan an output. Offline analysis is available.";
    else
    {
        const auto restored=monitoring.restore(*identity,monitoringOutputs,source.layout);
        session.sourceDocument.routing=restored.routing;routingNotice=restored.explanation;
    }
    rebuildOutputs();
}
void StandaloneAudioPanel::rememberRouting()
{
    const auto id=connection.monitoringIdentity();const auto& source=session.sourceDocument;
    if(audio.active() || session.job.isActive() || session.persistenceLoading || !id || !source.layout.confirmed)return;
    if(MonitoringSettings::layoutKey(source.layout).isEmpty())return; // Existing other layouts stay manually usable.
    if(source.routing.validate(source.layout.count,std::min(maxChannels,monitoringOutputs.size()),true).isNotEmpty())
    {routingNotice="Resolve duplicate or unavailable outputs before this routing can be remembered.";return;}
    if(monitoring.rememberRouting(*id,monitoringOutputs,source.layout,source.routing))
        routingNotice="Remembered "+MonitoringSettings::layoutKey(source.layout)+" routing for "+id->name+", including muted roles.";
}

void StandaloneAudioPanel::rebuildOutputs()
{
    auto names = monitoringOutputs;
    while(names.size()>maxChannels)names.remove(names.size()-1);
    for (int channel=0;channel<maxChannels;++channel)
    {
        auto& box=outputBoxes[size_t(channel)];
        box.clear (juce::dontSendNotification);
        box.addItem ("Muted", 1);
        for (int i = 0; i < names.size(); ++i) box.addItem (juce::String (i + 1) + ": " + names[i], i + 2);
        box.setSelectedId(session.sourceDocument.routing.fileToDevice[size_t(channel)]+2,juce::dontSendNotification);
    }
}

void StandaloneAudioPanel::sourceChanged(bool freshImport)
{
    auto& source=session.sourceDocument;
    reviewing=freshImport && source.hasSource && !source.assembly && !source.composition
        && (source.layout.count==2 || source.layout.count==6);
    editingOrder=false;
    if(reviewing)
    {
        const auto proposal=roleMemory.suggest(source.media);
        source.layout=proposal.order;
        suggestionLabel.setText(proposal.explanation,juce::dontSendNotification);
        suggestionLabel.setTooltip(proposal.explanation);
    }
    else suggestionLabel.setText(source.layout.confirmed?"Project channel assignments":source.media.layoutProvenance,juce::dontSendNotification);
    routingContext.clear(); // A new/opened programme must map its own confirmed file order.
    audio.sourceChanged (session.sourceDocument.media);
    connection.programmeChanged(source.hasSource && !session.sourceContentInvalid && !session.persistenceLoading ? source.media.sampleRate : 0);
    connectionChanged();
    actionError.clear();
    refresh();
    resized();
}

void StandaloneAudioPanel::synchronize()
{
    refreshMonitoring();
    const bool previouslyInvalid = session.sourceContentInvalid;
    if ((session.sourceDocument.assembly || session.sourceDocument.composition) && !session.sourceContentInvalid && !session.persistenceLoading
        && !session.sourceDocument.matchesFileStamps())
        session.measurementUnsuccessful(session.auditionDomain, ErrorCode::sourceChanged);
    if (audio.failureCode() == ErrorCode::sourceChanged && ! session.sourceContentInvalid)
        session.measurementUnsuccessful (session.auditionDomain, ErrorCode::sourceChanged);
    audio.synchronize (session.auditionRequest(), session.sourceDocument.routing);
    if (!previouslyInvalid && session.sourceContentInvalid) onChanged();
}

void StandaloneAudioPanel::refresh()
{
    synchronize();
    const auto& source = session.sourceDocument;
    const bool playing = audio.active();
    const bool editable = ! playing && ! session.job.isActive() && !session.persistenceLoading;
    const auto request = session.auditionRequest();
    const auto requestError = request.validationError();
    domainBox.setEnabled (editable);
    domainBox.setSelectedId (session.auditionDomain == MeasurementDomain::source ? 1 : 2, juce::dontSendNotification);
    prepareButton.setEnabled (editable && requestError.isEmpty());
    prepareButton.setTooltip (requestError);
    const bool current = session.auditionDomain == MeasurementDomain::currentState;
    const auto& summary = current ? session.currentState : session.source;
    const auto expected = requestError.isEmpty() ? OfflineAnalysisService::identityFor (request.context, request.domain, request.context.media.sha256) : juce::String();
    const bool measured = summary.hasMeasurement && summary.analysis && ! session.sourceContentInvalid
        && summary.analysis->identity == expected && (! current || ! session.currentStateNeedsRemeasurement());
    juce::String measurementStatus = measured ? "Measurements current" : summary.hasMeasurement
        ? "UNVERIFIED: measurements stale" : current && session.currentStateNeedsRemeasurement()
        ? "UNVERIFIED: measurements unavailable; remeasure Current State" : "UNVERIFIED: measurements unavailable";
    auditionLabel.setText ((current ? "Current State" : "Source") + juce::String (session.sourceDocument.composition ? " | composed programme | Processing " : " | Processing ")
        + juce::String (request.context.processing.staticGainDb, 12) + " dB global"
        + (request.context.processing.sectionGains ? " + applied section gains" : "") + " | " + measurementStatus
        + "\n" + (requestError.isNotEmpty() ? requestError : "Playback does not measure or establish compliance."), juce::dontSendNotification);
    gain.setEnabled (true); // Listening volume remains live, even after an overload fault.
    automaticReconnect.setEnabled(editable);
    deviceBox.setEnabled (editable);
    scanButton.setEnabled (editable);
    openButton.setEnabled (editable && connection.selection().has_value());
    closeButton.setEnabled (audio.connected() || connection.selection().has_value());
    const auto routeError = source.routing.validate (source.layout.count, audio.outputNames().size());
    playButton.setEnabled (editable && audio.connected() && requestError.isEmpty() && routeError.isEmpty());
    stopButton.setEnabled (playing);
    restartButton.setEnabled (!session.persistenceLoading && !session.job.isActive() && ! audio.preparing() && audio.connected() && requestError.isEmpty() && routeError.isEmpty());
    seekButton.setEnabled (!session.persistenceLoading && !session.job.isActive() && ! audio.preparing() && requestError.isEmpty());
    seekPosition.setEnabled (seekButton.isEnabled());
    editOrder.setVisible(reviewing);
    editOrder.setEnabled(editable);
    confirmButton.setButtonText(reviewing?"Confirm roles and analyse":"Confirm roles");
    confirmButton.setEnabled (editable && source.hasSource && source.layout.isValid() && ! source.layout.confirmed);
    layoutPreset.setEnabled (editable && source.hasSource && source.layout.count == 8);
    layoutPreset.setTooltip (source.media.layoutProvenance);
    for (int channel = 0; channel < maxChannels; ++channel)
    {
        const bool visible = channel < source.layout.count;
        channelLabels[channel].setVisible (visible);
        roleBoxes[channel].setVisible (visible);
        outputBoxes[channel].setVisible (visible);
        roleBoxes[channel].setEnabled (editable && !source.assembly && !source.composition && (!reviewing || editingOrder));
        channelLabels[channel].setText(source.channelName(channel),juce::dontSendNotification);
        channelLabels[channel].setTooltip(source.channelName(channel));
        roleBoxes[channel].setEditableText (source.layout.count == 8);
        roleBoxes[channel].setTooltip (source.layout.count == 8
            ? "Select a role or enter a documented System I alias (L/R/C/LFE/Lss/Rss/Lrs/Rrs). Confirm roles afterwards."
            : "Select the original file-channel speaker role, then confirm.");
        outputBoxes[channel].setEnabled(editable && source.layout.confirmed && !monitoringOutputs.isEmpty());
        int selected = 1;
        for (size_t role = 0; role < supportedRoles.size(); ++role)
            if (source.layout.roles[channel] == supportedRoles[role]) selected = static_cast<int> (role) + 2;
        // A timer refresh must not erase an alias while the user is typing it.
        if ((reviewing && !editingOrder) || !roleBoxes[channel].hasKeyboardFocus (true))
            roleBoxes[channel].setSelectedId (selected, juce::dontSendNotification);
    }
    const auto monitorIdentity=connection.monitoringIdentity();
    rememberRoutingButton.setEnabled(editable && monitorIdentity.has_value() && source.layout.confirmed
        && MonitoringSettings::layoutKey(source.layout).isNotEmpty() && !monitoringOutputs.isEmpty()
        && source.routing.validate(source.layout.count,std::min(maxChannels,monitoringOutputs.size()),true).isEmpty());
    auto monitorText=routingNotice;
    const auto savedVolume=connection.selection()?monitoring.volume(*connection.selection()):std::optional<double>{};
    monitorText+="\nMonitor: "+juce::String(gain.getValue(),1)+" dB"+(savedVolume?" (last user setting on this Mac/interface).":" (no saved value; initial default is -30 dB).");
    if(monitoring.warning().isNotEmpty())monitorText+="\n"+monitoring.warning();
    monitoringStatus.setText(monitorText,juce::dontSendNotification);monitoringStatus.setTooltip(monitorText);
    juce::String status = actionError.isNotEmpty() ? actionError : audio.status();
    status << "\n" << connection.status();
    if(connection.settingsWarning().isNotEmpty())status << "\n" << connection.settingsWarning();
    if (audio.preparing()) status << " | " << juce::String (juce::roundToInt (100 * audio.preparationProgress())) << "%";
    if (playing || audio.getState() == AuditionState::finished)
        status << "  |  " << juce::String (audio.positionSeconds(), 1) << " s  |  underruns " << juce::String (audio.underruns());
    status << "\n" << (source.layout.confirmed ? "Roles confirmed. " : "Review the proposed file-channel order and confirm. ");
    if (routeError.isNotEmpty()) status << routeError;
    status << "\n" << (isSystemI (source.layout) ? "BS.2051 System I 7.1 | " : "") << monitoringCoverage (source.layout, source.routing);
    statusLabel.setTooltip (source.media.layoutProvenance + (source.layout.manuallyInterpreted ? "\nExplicit manual interpretation selected." : ""));
    statusLabel.setText (status, juce::dontSendNotification);
}

void StandaloneAudioPanel::resized()
{
    auto area = getLocalBounds();
    auto row = area.removeFromTop (28);
    domainBox.setBounds (row.removeFromLeft (165).reduced (2));
    prepareButton.setBounds (row.removeFromLeft (150).reduced (2));
    auditionLabel.setBounds (area.removeFromTop (45).reduced (2));
    row = area.removeFromTop (28);
    scanButton.setBounds (row.removeFromRight (110).reduced (2));
    deviceBox.setBounds (row.reduced (2));
    automaticReconnect.setBounds(area.removeFromTop(28).reduced(2));
    row = area.removeFromTop (30);
    openButton.setBounds (row.removeFromLeft (105).reduced (2));
    closeButton.setBounds (row.removeFromLeft (105).reduced (2));
    playButton.setBounds (row.removeFromLeft (105).reduced (2));
    stopButton.setBounds (row.removeFromLeft (65).reduced (2));
    gain.setBounds (row.reduced (2));
    row = area.removeFromTop (28);
    restartButton.setBounds (row.removeFromLeft (90).reduced (2));
    seekButton.setBounds (row.removeFromRight (65).reduced (2));
    seekPosition.setBounds (row.reduced (2));
    statusLabel.setBounds (area.removeFromTop (112).reduced (2));
    layoutPreset.setVisible(session.sourceDocument.layout.count==8);
    if(layoutPreset.isVisible())layoutPreset.setBounds(area.removeFromTop(28).reduced(2));
    suggestionLabel.setBounds(area.removeFromTop(40).reduced(2));
    suggestionLabel.setJustificationType(juce::Justification::centredLeft);
    monitoringStatus.setBounds(area.removeFromTop(64).reduced(2));
    row = area.removeFromTop (30);
    rememberRoutingButton.setBounds(row.removeFromRight(150).reduced(2));
    confirmButton.setBounds (row.removeFromRight (205).reduced (2));
    if(reviewing)editOrder.setBounds(row.removeFromRight(100).reduced(2));
    heading.setBounds (row.reduced (2));
    viewport.setBounds (area);
    rows.setSize (juce::jmax (0, area.getWidth() - 16), session.sourceDocument.layout.count * 28);
    for (int channel = 0; channel < maxChannels; ++channel)
    {
        auto bounds = juce::Rectangle<int> (0, channel * 28, rows.getWidth(), 28).reduced (2);
        channelLabels[channel].setBounds (bounds.removeFromLeft ((session.sourceDocument.assembly || session.sourceDocument.composition) ? 240 : 80));
        roleBoxes[channel].setBounds (bounds.removeFromLeft (210).reduced (2, 0));
        outputBoxes[channel].setBounds (bounds.reduced (2, 0));
    }
}
// [END LS-STANDALONE-AUDIO-PANEL-IMPLEMENTATION]
}

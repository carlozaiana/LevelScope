#include "StandaloneAudioPanel.h"
#include "StandaloneMeasurementLayout.h"

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-AUDIO-PANEL-IMPLEMENTATION]
StandaloneAudioPanel::StandaloneAudioPanel (StandaloneSessionModel& model, StandaloneAudition& engine,
                                           std::function<void()> changed)
    : session (model), audio (engine), onChanged (std::move (changed))
{
    session.auditionBusy = [&engine] { return engine.active() && engine.getState() != AuditionState::cancelling; };
    for (auto* component : std::initializer_list<juce::Component*> {
        &typeBox, &deviceBox, &scanButton, &openButton, &closeButton, &playButton,
        &stopButton, &restartButton, &seekButton, &seekPosition, &confirmButton, &gain, &statusLabel, &heading, &viewport, &domainBox, &prepareButton, &auditionLabel, &layoutPreset })
        addAndMakeVisible (*component);
    const auto& types = audio.deviceTypes();
    for (int i = 0; i < types.size(); ++i) typeBox.addItem (types[i]->getTypeName(), i + 1);
    if (! types.isEmpty()) typeBox.setSelectedId (1, juce::dontSendNotification);
    deviceBox.setTextWhenNothingSelected ("Scan, then select an output device");
    typeBox.onChange = [this] { scan(); };
    scanButton.onClick = [this] { scan(); };
    deviceBox.onChange = [this]
    {
        if (session.persistenceLoading) return;
        actionError = audio.selectDevice (typeBox.getSelectedId() - 1, deviceBox.getText());
        gain.setValue (-30, juce::sendNotificationSync);
        rebuildOutputs();
        refresh();
    };
    openButton.onClick = [this] { if (session.persistenceLoading) return; actionError = audio.connect (session.sourceDocument.media.sampleRate); refresh(); };
    closeButton.onClick = [this] { audio.disconnect(); actionError.clear(); refresh(); onChanged(); };
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
    confirmButton.onClick = [this] { session.confirmLayout(); onChanged(); };
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
    gain.onValueChange = [this] { audio.setMonitorGainDb (static_cast<float> (gain.getValue())); };
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
        roleBoxes[channel].addItem ("Unassigned", 1);
        for (size_t role = 0; role < supportedRoles.size(); ++role)
            roleBoxes[channel].addItem (juce::AudioChannelSet::getChannelTypeName (supportedRoles[role]), static_cast<int> (role) + 2);
        roleBoxes[channel].onChange = [this, channel]
        {
            const auto index = roleBoxes[channel].getSelectedId() - 2;
            session.setChannelRole (channel, juce::isPositiveAndBelow (index, static_cast<int> (supportedRoles.size()))
                ? supportedRoles[static_cast<size_t> (index)] : session.sourceDocument.layout.count == 8
                ? systemIRoleForAlias (roleBoxes[channel].getText()) : juce::AudioChannelSet::unknown);
            onChanged();
        };
        outputBoxes[channel].onChange = [this, channel]
        {
            session.setOutputRoute (channel, outputBoxes[channel].getSelectedId() - 2);
            actionError.clear();
            refresh();
        };
    }
    rebuildOutputs();
    sourceChanged();
}

void StandaloneAudioPanel::scan()
{
    if (session.persistenceLoading) return;
    audio.disconnect();
    deviceBox.clear (juce::dontSendNotification);
    const auto devices = audio.scanOutputs (typeBox.getSelectedId() - 1);
    for (int i = 0; i < devices.size(); ++i) deviceBox.addItem (devices[i], i + 1);
    audio.selectDevice (-1, {});
    rebuildOutputs();
    actionError = devices.isEmpty() ? "No output devices found" : "Select a device; no device is open";
    refresh();
}

void StandaloneAudioPanel::rebuildOutputs()
{
    session.sourceDocument.routing = OutputRouting();
    const auto names = audio.outputNames();
    for (auto& box : outputBoxes)
    {
        box.clear (juce::dontSendNotification);
        box.addItem ("Muted", 1);
        for (int i = 0; i < names.size(); ++i) box.addItem (juce::String (i + 1) + ": " + names[i], i + 2);
        box.setSelectedId (1, juce::dontSendNotification);
    }
}

void StandaloneAudioPanel::sourceChanged()
{
    audio.sourceChanged (session.sourceDocument.media);
    actionError.clear();
    rebuildOutputs();
    refresh();
    resized();
}

void StandaloneAudioPanel::synchronize()
{
    const bool previouslyInvalid = session.sourceContentInvalid;
    if (session.sourceDocument.assembly && !session.sourceContentInvalid && !session.persistenceLoading
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
    auditionLabel.setText ((current ? "Current State" : "Source") + juce::String (" | Processing ")
        + juce::String (request.context.processing.staticGainDb, 12) + " dB | " + measurementStatus
        + "\n" + (requestError.isNotEmpty() ? requestError : "Playback does not measure or establish compliance."), juce::dontSendNotification);
    gain.setEnabled (true); // Listening volume remains live, even after an overload fault.
    typeBox.setEnabled (editable);
    deviceBox.setEnabled (editable);
    scanButton.setEnabled (editable);
    openButton.setEnabled (editable && source.hasSource && deviceBox.getSelectedId() > 0);
    closeButton.setEnabled (audio.connected());
    const auto routeError = source.routing.validate (source.layout.count, audio.outputNames().size());
    playButton.setEnabled (editable && audio.connected() && requestError.isEmpty() && routeError.isEmpty());
    stopButton.setEnabled (playing);
    restartButton.setEnabled (!session.persistenceLoading && !session.job.isActive() && ! audio.preparing() && audio.connected() && requestError.isEmpty() && routeError.isEmpty());
    seekButton.setEnabled (!session.persistenceLoading && !session.job.isActive() && ! audio.preparing() && requestError.isEmpty());
    seekPosition.setEnabled (seekButton.isEnabled());
    confirmButton.setEnabled (editable && source.hasSource && source.layout.isValid() && ! source.layout.confirmed);
    layoutPreset.setEnabled (editable && source.hasSource && source.layout.count == 8);
    layoutPreset.setTooltip (source.media.layoutProvenance);
    for (int channel = 0; channel < maxChannels; ++channel)
    {
        const bool visible = channel < source.layout.count;
        channelLabels[channel].setVisible (visible);
        roleBoxes[channel].setVisible (visible);
        outputBoxes[channel].setVisible (visible);
        roleBoxes[channel].setEnabled (editable && !source.assembly);
        channelLabels[channel].setText(source.channelName(channel),juce::dontSendNotification);
        channelLabels[channel].setTooltip(source.channelName(channel));
        roleBoxes[channel].setEditableText (source.layout.count == 8);
        roleBoxes[channel].setTooltip (source.layout.count == 8
            ? "Select a role or enter a documented System I alias (L/R/C/LFE/Lss/Rss/Lrs/Rrs). Confirm roles afterwards."
            : "Select the original file-channel speaker role, then confirm.");
        outputBoxes[channel].setEnabled (editable);
        int selected = 1;
        for (size_t role = 0; role < supportedRoles.size(); ++role)
            if (source.layout.roles[channel] == supportedRoles[role]) selected = static_cast<int> (role) + 2;
        // A timer refresh must not erase an alias while the user is typing it.
        if (! roleBoxes[channel].hasKeyboardFocus (true))
            roleBoxes[channel].setSelectedId (selected, juce::dontSendNotification);
    }
    juce::String status = actionError.isNotEmpty() ? actionError : audio.status();
    if (audio.preparing()) status << " | " << juce::String (juce::roundToInt (100 * audio.preparationProgress())) << "%";
    if (playing || audio.getState() == AuditionState::finished)
        status << "  |  " << juce::String (audio.positionSeconds(), 1) << " s  |  underruns " << juce::String (audio.underruns());
    status << "\n" << (source.layout.confirmed ? "Roles confirmed. " : "Review decoder-suggested roles; then confirm. ");
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
    typeBox.setBounds (row.removeFromLeft (115).reduced (2));
    scanButton.setBounds (row.removeFromRight (110).reduced (2));
    deviceBox.setBounds (row.reduced (2));
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
    statusLabel.setBounds (area.removeFromTop (72).reduced (2));
    layoutPreset.setBounds (area.removeFromTop (28).reduced (2));
    row = area.removeFromTop (30);
    confirmButton.setBounds (row.removeFromRight (115).reduced (2));
    heading.setBounds (row.reduced (2));
    viewport.setBounds (area);
    rows.setSize (juce::jmax (0, area.getWidth() - 16), session.sourceDocument.layout.count * 28);
    for (int channel = 0; channel < maxChannels; ++channel)
    {
        auto bounds = juce::Rectangle<int> (0, channel * 28, rows.getWidth(), 28).reduced (2);
        channelLabels[channel].setBounds (bounds.removeFromLeft (session.sourceDocument.assembly ? 240 : 80));
        roleBoxes[channel].setBounds (bounds.removeFromLeft (210).reduced (2, 0));
        outputBoxes[channel].setBounds (bounds.reduced (2, 0));
    }
}
// [END LS-STANDALONE-AUDIO-PANEL-IMPLEMENTATION]
}

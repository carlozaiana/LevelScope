#include "StandaloneAudioPanel.h"

namespace levelscope::standalone
{
// [BEGIN LS-STANDALONE-AUDIO-PANEL-IMPLEMENTATION]
StandaloneAudioPanel::StandaloneAudioPanel (StandaloneSessionModel& model, StandaloneAudition& engine,
                                           std::function<void()> changed)
    : session (model), audio (engine), onChanged (std::move (changed))
{
    for (auto* component : std::initializer_list<juce::Component*> {
        &typeBox, &deviceBox, &scanButton, &openButton, &closeButton, &playButton,
        &stopButton, &confirmButton, &gain, &statusLabel, &heading, &viewport })
        addAndMakeVisible (*component);
    const auto& types = audio.deviceTypes();
    for (int i = 0; i < types.size(); ++i) typeBox.addItem (types[i]->getTypeName(), i + 1);
    if (! types.isEmpty()) typeBox.setSelectedId (1, juce::dontSendNotification);
    deviceBox.setTextWhenNothingSelected ("Scan, then select an output device");
    typeBox.onChange = [this] { scan(); };
    scanButton.onClick = [this] { scan(); };
    deviceBox.onChange = [this]
    {
        actionError = audio.selectDevice (typeBox.getSelectedId() - 1, deviceBox.getText());
        gain.setValue (-30, juce::sendNotificationSync);
        rebuildOutputs();
        refresh();
    };
    openButton.onClick = [this] { actionError = audio.connect (session.sourceDocument.media.sampleRate); refresh(); };
    closeButton.onClick = [this] { audio.disconnect(); actionError.clear(); refresh(); };
    playButton.onClick = [this]
    {
        const auto& source = session.sourceDocument;
        actionError = audio.play (source.file, source.media, source.layout, source.routing);
        refresh();
    };
    stopButton.onClick = [this] { audio.stop(); actionError.clear(); refresh(); };
    confirmButton.onClick = [this] { session.confirmLayout(); onChanged(); };
    gain.setRange (-60, 0, 0.5);
    gain.setSliderStyle (juce::Slider::LinearHorizontal);
    gain.setTextBoxStyle (juce::Slider::TextBoxRight, false, 80, 24);
    gain.setTextValueSuffix (" dB");
    gain.setTooltip ("Monitor gain only. Original source and future render settings are unaffected.");
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
        channelLabels[channel].setText ("File " + juce::String (channel + 1), juce::dontSendNotification);
        roleBoxes[channel].addItem ("Unassigned", 1);
        for (size_t role = 0; role < supportedRoles.size(); ++role)
            roleBoxes[channel].addItem (juce::AudioChannelSet::getChannelTypeName (supportedRoles[role]), static_cast<int> (role) + 2);
        roleBoxes[channel].onChange = [this, channel]
        {
            const auto index = roleBoxes[channel].getSelectedId() - 2;
            session.setChannelRole (channel, juce::isPositiveAndBelow (index, static_cast<int> (supportedRoles.size()))
                ? supportedRoles[static_cast<size_t> (index)] : juce::AudioChannelSet::unknown);
            onChanged();
        };
        outputBoxes[channel].onChange = [this, channel]
        {
            session.sourceDocument.routing.fileToDevice[channel] = outputBoxes[channel].getSelectedId() - 2;
            actionError.clear();
            refresh();
        };
    }
    rebuildOutputs();
    sourceChanged();
}

void StandaloneAudioPanel::scan()
{
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
    audio.stop();
    actionError.clear();
    rebuildOutputs();
    refresh();
    resized();
}

void StandaloneAudioPanel::refresh()
{
    const auto& source = session.sourceDocument;
    const bool playing = audio.active();
    const bool editable = ! playing && ! session.job.isActive();
    typeBox.setEnabled (editable);
    deviceBox.setEnabled (editable);
    scanButton.setEnabled (editable);
    openButton.setEnabled (editable && source.hasSource && deviceBox.getSelectedId() > 0);
    closeButton.setEnabled (audio.connected());
    const auto routeError = source.routing.validate (source.layout.count, audio.outputNames().size());
    playButton.setEnabled (editable && audio.connected() && source.hasSource && source.layout.confirmed && routeError.isEmpty());
    stopButton.setEnabled (playing);
    confirmButton.setEnabled (editable && source.hasSource && source.layout.isValid() && ! source.layout.confirmed);
    for (int channel = 0; channel < maxChannels; ++channel)
    {
        const bool visible = channel < source.layout.count;
        channelLabels[channel].setVisible (visible);
        roleBoxes[channel].setVisible (visible);
        outputBoxes[channel].setVisible (visible);
        roleBoxes[channel].setEnabled (editable);
        outputBoxes[channel].setEnabled (editable);
        int selected = 1;
        for (size_t role = 0; role < supportedRoles.size(); ++role)
            if (source.layout.roles[channel] == supportedRoles[role]) selected = static_cast<int> (role) + 2;
        roleBoxes[channel].setSelectedId (selected, juce::dontSendNotification);
    }
    juce::String status = actionError.isNotEmpty() ? actionError : audio.status();
    if (playing || audio.getState() == AuditionState::finished)
        status << "  |  " << juce::String (audio.positionSeconds(), 1) << " s  |  underruns " << juce::String (audio.underruns());
    status << "\n" << (source.layout.confirmed ? "Roles confirmed. " : "Review decoder-suggested roles; then confirm. ");
    if (routeError.isNotEmpty()) status << routeError;
    statusLabel.setText (status, juce::dontSendNotification);
}

void StandaloneAudioPanel::resized()
{
    auto area = getLocalBounds();
    auto row = area.removeFromTop (28);
    typeBox.setBounds (row.removeFromLeft (115).reduced (2));
    scanButton.setBounds (row.removeFromRight (110).reduced (2));
    deviceBox.setBounds (row.reduced (2));
    row = area.removeFromTop (30);
    openButton.setBounds (row.removeFromLeft (105).reduced (2));
    closeButton.setBounds (row.removeFromLeft (105).reduced (2));
    playButton.setBounds (row.removeFromLeft (105).reduced (2));
    stopButton.setBounds (row.removeFromLeft (65).reduced (2));
    gain.setBounds (row.reduced (2));
    statusLabel.setBounds (area.removeFromTop (55).reduced (2));
    row = area.removeFromTop (30);
    confirmButton.setBounds (row.removeFromRight (115).reduced (2));
    heading.setBounds (row.reduced (2));
    viewport.setBounds (area);
    rows.setSize (juce::jmax (0, area.getWidth() - 16), session.sourceDocument.layout.count * 28);
    for (int channel = 0; channel < maxChannels; ++channel)
    {
        auto bounds = juce::Rectangle<int> (0, channel * 28, rows.getWidth(), 28).reduced (2);
        channelLabels[channel].setBounds (bounds.removeFromLeft (60));
        roleBoxes[channel].setBounds (bounds.removeFromLeft (210).reduced (2, 0));
        outputBoxes[channel].setBounds (bounds.reduced (2, 0));
    }
}
// [END LS-STANDALONE-AUDIO-PANEL-IMPLEMENTATION]
}

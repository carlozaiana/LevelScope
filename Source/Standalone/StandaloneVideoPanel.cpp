#include "StandaloneVideoPanel.h"

namespace levelscope::standalone
{
StandaloneVideoPanel::StandaloneVideoPanel (StandaloneSessionModel& s, StandaloneAudition& a) : session (s), audio (a)
{
    for (auto* component : std::initializer_list<juce::Component*> { &attach, &clear, &apply, &streams, &offset, &offsetLabel, &status, &metadata })
        addAndMakeVisible (*component);
    metadata.setMultiLine (true); metadata.setReadOnly (true); metadata.setCaretVisible (false);
    offset.setText ("0");
    offset.setTooltip ("Signed seconds or exact fraction, e.g. -0.250 or 1001/30000. Positive starts picture later.");
    offsetLabel.setText ("Video starts on audio (seconds):", juce::dontSendNotification);
    streams.setTextWhenNothingSelected ("Video stream");
    attach.onClick = [this]
    {
        chooser = std::make_unique<juce::FileChooser> ("Attach a local reference picture", juce::File(), "*.mov;*.mp4;*.m4v");
        juce::Component::SafePointer<StandaloneVideoPanel> safe (this);
        chooser->launchAsync (juce::FileBrowserComponent::openMode | juce::FileBrowserComponent::canSelectFiles,
            [safe] (const juce::FileChooser& picker)
            {
                if (!safe || picker.getResult() == juce::File()) return;
                safe->session.referenceVideo.file = picker.getResult();
                safe->session.referenceVideo.sha256.clear(); safe->session.referenceVideo.streamId = 0;
                safe->actionError.clear(); safe->displayedMetadata = nullptr;
                safe->metadata.clear(); safe->streams.clear (juce::dontSendNotification);
                safe->player.attach (picker.getResult());
            });
    };
    clear.onClick = [this]
    {
        player.clear(); session.referenceVideo = {}; offset.setText ("0");
        streams.clear (juce::dontSendNotification); metadata.clear(); displayedMetadata = nullptr; actionError.clear();
    };
    apply.onClick = [this]
    {
        try
        {
            const auto time = video::Time::parse (offset.getText().trim().toStdString());
            // Validate the current mapping before accepting the edit.
            const auto rate = session.sourceDocument.media.sampleRate;
            (void) (video::Time (audio.positionFrames(), rate > 0 ? static_cast<int64_t> (rate) : 48000) - time);
            session.referenceVideo.startOnAudio = time; player.invalidate(); actionError.clear();
        }
        catch (const std::exception& error) { actionError = error.what(); }
    };
    offset.onReturnKey = [this] { apply.triggerClick(); };
    streams.onChange = [this]
    {
        const int id = streams.getSelectedId(); if (id <= 0) return;
        session.referenceVideo.streamId = id; displayedMetadata = nullptr; actionError.clear(); metadata.clear();
        player.attach (session.referenceVideo.file, id);
    };
    setSize (700, 600); startTimerHz (120);
}
StandaloneVideoPanel::~StandaloneVideoPanel() { stopTimer(); }
void StandaloneVideoPanel::timerCallback()
{
    try
    {
        const auto rate = session.sourceDocument.media.sampleRate;
        const auto time = video::Time (audio.positionFrames(), rate > 0 ? static_cast<int64_t> (rate) : 48000)
                        - session.referenceVideo.startOnAudio;
        view = player.update (time, audio.transportGeneration());
        if (view.asset && displayedMetadata != view.asset.get())
        {
            displayedMetadata = view.asset.get();
            const auto& asset = *view.asset;
            session.referenceVideo.sha256 = asset.sha256; session.referenceVideo.streamId = asset.selectedStream;
            streams.clear (juce::dontSendNotification);
            for (const auto& stream : asset.streams)
            {
                streams.addItem ("Stream " + juce::String (stream.id) + " " + stream.codec + " "
                    + juce::String (stream.width) + "x" + juce::String (stream.height), stream.id);
                streams.setItemEnabled (stream.id, stream.supported);
            }
            streams.setSelectedId (asset.selectedStream, juce::dontSendNotification);
            const auto selected = std::find_if (asset.streams.begin(), asset.streams.end(), [&asset] (const auto& item) { return item.id == asset.selectedStream; });
            juce::String text;
            text << asset.file.getFullPathName() << "\n" << asset.backend << "\nSHA-256 " << asset.sha256
                 << "\n" << (asset.variableRate ? "VFR" : "CFR") << " | " << asset.frames.size() << " frames"
                 << " | nominal/minimum-duration rate " << selected->nominalRate.text() << " fps | timebase 1/" << selected->timeScale
                 << "\nCoverage " << asset.frames.front().pts.text() << " to " << asset.frames.back().end.text() << " s"
                 << "\n" << asset.timecodeStatus << "\nEmbedded audio tracks: " << asset.embeddedAudioTracks << " (excluded from playback)";
            metadata.setText (text, false);
        }
        juce::String text = actionError.isNotEmpty() ? actionError : view.message;
        if (view.picture && view.asset)
        {
            text << " | frame " << view.picture->index << " | PTS " << view.asset->frames[view.picture->index].pts.text();
            if (view.asset->timecode)
            {
                const auto& tc = *view.asset->timecode;
                const auto number = video::Time (tc.startFrame, 1) + (time - tc.origin) * tc.rate;
                text << " | TC " << video::timecodeLabel (number.floorTicks (1), tc.rate, tc.dropFrame, tc.wraps24Hours);
            }
        }
        text << "\nLate polls " << view.latePolls << " | skipped frames " << view.droppedFrames
             << " | output latency " << audio.outputLatencyFrames() << " samples"
             << (audio.hasHardwareClock() ? " | HAL output timestamp" : " | callback-time estimate / stopped");
        status.setText (text, juce::dontSendNotification); repaint (pictureBounds);
    }
    catch (const std::exception& error)
    {
        view.picture.reset(); status.setText ("Video timing error: " + juce::String (error.what()), juce::dontSendNotification);
        repaint (pictureBounds);
    }
}
void StandaloneVideoPanel::paint (juce::Graphics& g)
{
    g.fillAll (juce::Colour (0xff181b20)); g.setColour (juce::Colours::black); g.fillRect (pictureBounds);
    // A timer snapshot can straddle a VFR boundary before the next paint.
    // Select from the bounded decoded queue against the clock at paint time.
    try
    {
        const auto rate = session.sourceDocument.media.sampleRate;
        const auto time = video::Time (audio.positionFrames(), rate > 0 ? static_cast<int64_t> (rate) : 48000)
                        - session.referenceVideo.startOnAudio;
        const auto painted = player.update (time, audio.transportGeneration());
        if (painted.picture) g.drawImageWithin (painted.picture->image, pictureBounds.getX(), pictureBounds.getY(),
            pictureBounds.getWidth(), pictureBounds.getHeight(), juce::RectanglePlacement::centred);
        else { g.setColour (juce::Colours::lightgrey); g.drawText (painted.message, pictureBounds.reduced (12), juce::Justification::centred); }
    }
    catch (const std::exception& error)
    { g.setColour (juce::Colours::lightgrey); g.drawText (error.what(), pictureBounds.reduced (12), juce::Justification::centred); }
}
void StandaloneVideoPanel::resized()
{
    auto area = getLocalBounds().reduced (12);
    auto row = area.removeFromTop (30);
    attach.setBounds (row.removeFromLeft (135).reduced (2)); clear.setBounds (row.removeFromLeft (100).reduced (2)); streams.setBounds (row.reduced (2));
    row = area.removeFromTop (32);
    offsetLabel.setBounds (row.removeFromLeft (235)); apply.setBounds (row.removeFromRight (105).reduced (2)); offset.setBounds (row.reduced (2));
    status.setBounds (area.removeFromBottom (58)); metadata.setBounds (area.removeFromBottom (150));
    pictureBounds = area.reduced (0, 8);
}
}

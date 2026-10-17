// Explicit opt-in hardware/display qualification. Never enumerates or opens a
// default device. No installed app/plugin or system default is modified.
#include "Standalone/StandaloneVideo.h"
#include "Standalone/StandaloneAudition.h"
#include <juce_gui_basics/juce_gui_basics.h>
#include <iostream>

using namespace levelscope::standalone;
namespace
{
class Run final : public juce::Component, private juce::Timer
{
public:
    Run (juce::StringArray args)
    {
        auto argument = [&] (const juce::String& name)
        {
            const auto i = args.indexOf (name);
            return i >= 0 && i + 1 < args.size() ? args[i + 1] : juce::String();
        };
        if (!args.contains ("--allow-audio")) throw std::runtime_error ("Requires explicit --allow-audio after user coordination");
        audioFile = juce::File (argument ("--audio")); videoFile = juce::File (argument ("--video"));
        reportFile = juce::File (argument ("--report")); deviceName = argument ("--device");
        duration = argument ("--seconds").getDoubleValue();
        if (deviceName.isEmpty() || duration <= 0 || duration > 7200 || reportFile.exists())
            throw std::runtime_error ("Specify --device, --seconds (1..7200) and a new --report file");
        JobControl job; auto result = validateMedia (audioFile, job);
        if (!result.succeeded()) throw std::runtime_error (result.message.toStdString());
        media = result.media;
        if (media.channels != 2 || (media.sampleRate != 48000 && media.sampleRate != 44100) || media.durationSeconds() < duration)
            throw std::runtime_error ("Qualification requires a sufficiently long 44.1 or 48 kHz stereo source");
        layout.count = 2; layout.roles[0] = juce::AudioChannelSet::left; layout.roles[1] = juce::AudioChannelSet::right; layout.confirmed = true;
        routing.fileToDevice[0] = 0; routing.fileToDevice[1] = 1;
        muted = args.contains ("--muted");
        manualStart = args.contains ("--manual-start");
        startButton.setButtonText ("Start " + juce::String (duration, 0) + "-second check");
        addChildComponent (startButton); startButton.setVisible (manualStart); startButton.setBounds (10, 10, 195, 30);
        startButton.onClick = [this] { requestedStart = true; startButton.setVisible (false); startedNs = StandaloneAudition::hostTimeNs(); };
        if (muted) routing = OutputRouting();
        audio.sourceChanged (media);
        // Scan/selection only here; open waits for successful video preroll.
        bool found = false;
        for (int type = 0; type < audio.deviceTypes().size(); ++type)
            if (audio.scanOutputs (type).contains (deviceName))
            {
                if (found) throw std::runtime_error ("Ambiguous device name");
                const auto error = audio.selectDevice (type, deviceName);
                if (error.isNotEmpty()) throw std::runtime_error (error.toStdString()); found = true;
            }
        if (!found) throw std::runtime_error ("Requested device was not found");
        if (audio.outputNames().size() != 2) throw std::runtime_error ("Use an explicitly coordinated two-output ordinary device");
        audio.setMonitorGainDb (-30);
        player.attach (videoFile); startedNs = StandaloneAudition::hostTimeNs();
        setSize (800, 520); startTimerHz (120);
    }
    ~Run() override { stopTimer(); audio.disconnect(); }
    void paint (juce::Graphics& g) override
    {
        const auto paintStart = StandaloneAudition::hostTimeNs();
        const video::Time paintTime (audio.positionFrames(), static_cast<int64_t> (media.sampleRate));
        view = player.update (paintTime, audio.transportGeneration());
        g.fillAll (juce::Colours::black);
        if (view.picture)
            g.drawImageWithin (view.picture->image, 0, 0, getWidth(), getHeight() - 65, juce::RectanglePlacement::centred);
        if (audio.getState() == AuditionState::playing && view.asset && audio.positionSeconds() >= 2)
        {
            ++paintCount;
            const auto time = paintTime;
            const auto expected = video::frameAt (view.asset->frames, time);
            if (expected && view.picture)
            {
                const auto& wanted = view.asset->frames[*expected];
                const auto difference = (view.asset->frames[view.picture->index].pts - wanted.pts).seconds();
                const auto bound = (wanted.end - wanted.pts).seconds();
                const auto paintMs = (StandaloneAudition::hostTimeNs() - paintStart) / 1e6;
                maximumPaintMs = std::max (maximumPaintMs, paintMs);
                if (paintMs > bound * 1000) ++paintDeadlineFailures;
                maximumSkew = std::max (maximumSkew, std::abs (difference));
                if (std::abs (difference) > bound + 1e-9) ++boundFailures;
                skewSum += difference; ++skewCount;
                if (audio.positionSeconds() < 62) { firstSum += difference; ++firstCount; }
                if (audio.positionSeconds() > duration - 60) { lastSum += difference; ++lastCount; }
            }
            else if (expected) ++blankPaints;
        }
        g.setColour (juce::Colours::white);
        g.drawFittedText ("Coordinated audio/display qualification | " + deviceName + " | monitor -30 dB\n"
            + juce::String (audio.positionSeconds(), 2) + " / " + juce::String (duration, 0)
            + " s | " + view.message + "\nPaint skew max " + juce::String (maximumSkew * 1000, 3)
            + " ms | bound failures " + juce::String (boundFailures) + " | blank paints " + juce::String (blankPaints)
            + (muted ? " | ALL ROUTES MUTED" : ""),
            getLocalBounds().removeFromBottom (65), juce::Justification::centred, 3);
    }
private:
    void finish (juce::String error)
    {
        stopTimer();
        const auto seconds = audio.positionSeconds(); const auto underruns = audio.underruns();
        const bool hardware = audio.hasHardwareClock(); const auto latency = audio.outputLatencyFrames();
        audio.disconnect();
        juce::DynamicObject::Ptr report = new juce::DynamicObject();
        auto set = [&] (const char* key, juce::var value) { report->setProperty (key, std::move (value)); };
        set ("scope", "Internal audio-clock and UI paint selection; no optical/acoustic measurement");
        set ("passed", error.isEmpty() && boundFailures == 0 && paintDeadlineFailures == 0 && blankPaints == 0 && underruns == 0 && seconds >= duration && skewCount > 0);
        set ("error", error); set ("os", juce::SystemStats::getOperatingSystemName());
#if defined(__arm64__)
        set ("architecture", "arm64");
#else
        set ("architecture", "x86_64");
#endif
        set ("sample_rate", media.sampleRate); set ("device", deviceName); set ("audio_sha256", media.sha256);
        set ("video_sha256", view.asset ? view.asset->sha256 : juce::String());
        set ("backend", view.asset ? view.asset->backend : juce::String());
        set ("source_seconds", seconds); set ("requested_seconds", duration); set ("latency_frames", latency);
        set ("hal_output_timestamp", hardware); set ("paints", static_cast<juce::int64> (paintCount));
        set ("blank_paints", static_cast<juce::int64> (blankPaints)); set ("bound_failures", static_cast<juce::int64> (boundFailures));
        set ("maximum_frame_pts_skew_seconds", maximumSkew);
        set ("maximum_cpu_paint_ms", maximumPaintMs);
        set ("paint_deadline_failures", static_cast<juce::int64> (paintDeadlineFailures));
        set ("mean_skew_seconds", skewCount ? skewSum / skewCount : 0);
        set ("first_minute_mean_skew_seconds", firstCount ? firstSum / firstCount : 0);
        set ("last_minute_mean_skew_seconds", lastCount ? lastSum / lastCount : 0);
        set ("underrun_blocks", static_cast<juce::int64> (underruns));
        set ("late_polls", static_cast<juce::int64> (view.latePolls)); set ("skipped_frames", static_cast<juce::int64> (view.droppedFrames));
        set ("embedded_audio_outputs", 0);
        set ("all_routes_muted", muted);
        const auto json = juce::JSON::toString (juce::var (report.get()));
        const bool saved = reportFile.replaceWithText (json + "\n"); std::cout << json << '\n';
        juce::JUCEApplication::getInstance()->setApplicationReturnValue (saved && static_cast<bool> (report->getProperty ("passed")) ? 0 : 1);
        juce::JUCEApplication::quit();
    }
    void timerCallback() override
    {
        try
        {
            view = player.update (video::Time (audio.positionFrames(), static_cast<int64_t> (media.sampleRate)), audio.transportGeneration());
            if (view.state == video::Player::State::failed) { finish (view.message); return; }
            if (!opened && view.picture && view.queued >= 3 && (!manualStart || requestedStart))
            {
                if (view.asset->frames.back().end < video::Time (static_cast<int64_t> (duration), 1))
                    throw std::runtime_error ("Video does not cover the requested run");
                auto error = audio.connect (media.sampleRate);
                if (error.isEmpty()) error = audio.play (audioFile, media, layout, routing, true, muted);
                if (error.isNotEmpty()) throw std::runtime_error (error.toStdString()); opened = true;
                view = player.update (video::Time (audio.positionFrames(), static_cast<int64_t> (media.sampleRate)), audio.transportGeneration());
            }
            audio.poll (player.prerollReady()); repaint();
            if (audio.getState() == AuditionState::failed) { finish (audio.status()); return; }
            if (opened && audio.positionSeconds() >= duration) { finish ({}); return; }
            if ((!manualStart || requestedStart) && (StandaloneAudition::hostTimeNs() - startedNs) / 1e9 > duration + 120)
                finish ("Qualification timed out");
        }
        catch (const std::exception& error) { finish (error.what()); }
    }
    StandaloneAudition audio;
    video::Player player;
    video::Player::View view;
    MediaInfo media; ChannelLayout layout; OutputRouting routing;
    juce::File audioFile, videoFile, reportFile;
    juce::String deviceName;
    juce::TextButton startButton;
    bool opened = false, muted = false, manualStart = false, requestedStart = false;
    uint64_t startedNs = 0, paintCount = 0, boundFailures = 0, blankPaints = 0, paintDeadlineFailures = 0;
    double duration = 0, maximumPaintMs = 0, maximumSkew = 0, skewSum = 0, firstSum = 0, lastSum = 0;
    int64_t skewCount = 0, firstCount = 0, lastCount = 0;
};
class Window final : public juce::DocumentWindow
{
public:
    explicit Window (Run* component) : DocumentWindow ("LevelScope synchronization qualification", juce::Colours::black, allButtons)
    { setUsingNativeTitleBar (true); setContentOwned (component, true); centreWithSize (800, 520); setVisible (true); toFront (true);
      juce::Process::makeForegroundProcess(); }
    void closeButtonPressed() override { juce::JUCEApplication::quit(); }
};
class Application final : public juce::JUCEApplication
{
public:
    const juce::String getApplicationName() override { return "LevelScope Video Qualification"; }
    const juce::String getApplicationVersion() override { return "0.1"; }
    void initialise (const juce::String& command) override
    {
        try { juce::StringArray args; args.addTokens (command, true); args.removeEmptyStrings(); args = unquote (args);
              window = std::make_unique<Window> (new Run (args)); }
        catch (const std::exception& e) { std::cerr << e.what() << '\n'; setApplicationReturnValue (1); quit(); }
    }
    void shutdown() override { window.reset(); }
private:
    static juce::StringArray unquote (juce::StringArray args) { for (auto& arg : args) arg = arg.unquoted(); return args; }
    std::unique_ptr<Window> window;
};
}
START_JUCE_APPLICATION (Application)

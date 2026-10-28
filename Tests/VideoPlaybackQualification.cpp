// Explicit opt-in hardware/display qualification. Never enumerates or opens a
// default device. No installed app/plugin or system default is modified.
#include "Standalone/StandaloneVideo.h"
#include "Standalone/StandaloneAudition.h"
#include "Standalone/StandaloneMeasurementLayout.h"
#include "Standalone/SessionDocument.h"
#include <juce_gui_basics/juce_gui_basics.h>
#include <iostream>
#if LEVELSCOPE_VIDEO_PAINT_DIAGNOSTICS
#include "VideoPaint/Trace.h"
namespace paintdiag = levelscope::paintdiagnostic;
#endif
#if JUCE_MAC
#include <mach/mach.h>
#include <sys/resource.h>
#endif

using namespace levelscope::standalone;
namespace
{
struct ProcessingPresence final : ProcessedAudioObserver
{
    // Bounded callback evidence; full sample parity is qualified by M7's injected
    // capture suite. Sample both ends of each callback on every original channel.
    std::atomic<int64_t> frames { 0 };
    std::atomic<unsigned> seen { 0 };
    std::atomic<bool> valid { true };
    int expected = 0;
    void observe (const float* const* data, int channels, int count, int64_t) noexcept override
    {
        if (channels != expected || count <= 0) { valid.store (false); return; }
        unsigned present = seen.load();
        for (int c = 0; c < channels; ++c)
            if (data[c][0] != 0 || data[c][count - 1] != 0) present |= (1u << c);
        seen.store (present); frames.fetch_add (count);
    }
};
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
#if LEVELSCOPE_VIDEO_PAINT_DIAGNOSTICS
        traceFile=juce::File(argument("--trace"));
        if (!args.contains("--muted") || argument("--trace").isEmpty() || traceFile.exists()
            || duration>100 || args.contains("--manual-start"))
            throw std::runtime_error("Diagnostic requires muted routes, <=100 seconds and a new --trace file");
#endif
        JobControl job;
        std::optional<ProcessingPlan> savedProcessing;
        const bool sessionProcessing=args.contains("--session-processing");
        if(sessionProcessing&&(argument("--programme").isEmpty()||argument("--audition-state")!="current"||argument("--processing-gain").isNotEmpty()))
            throw std::runtime_error("Saved processing requires --programme and current audition without a gain override");
        if (argument("--programme").isNotEmpty())
        {
            if (argument("--audio").isNotEmpty()) throw std::runtime_error("Choose either --audio or --programme");
            auto loaded=sessionfile::read(juce::File(argument("--programme")),job);
            if(!loaded.document)throw std::runtime_error(loaded.message.toStdString());
            auto verified=sessionfile::verify(*loaded.document,job);
            if(!verified.document||!verified.document->audioResolved||!verified.prepared)
                throw std::runtime_error("A complete confirmed programme is required");
            if(sessionProcessing)
            {
                if(!verified.prepared->model.currentStateDocument.isInitialized)throw std::runtime_error("Saved Current State required");
                savedProcessing=verified.prepared->model.currentStateDocument.processingPlan();
            }
            request=verified.prepared->model.auditionRequest();request.domain=MeasurementDomain::source;request.context.processing=ProcessingPlan{};
            media=request.context.media;layout=request.context.fileLayout;
        }
        else {
        auto result = validateMedia (audioFile, job);
        if (!result.succeeded()) throw std::runtime_error (result.message.toStdString());
        media = result.media;
        if ((media.channels != 2 && ! isSystemI (media.suggestedLayout))
            || (media.sampleRate != 48000 && media.sampleRate != 44100) || media.durationSeconds() < duration)
            throw std::runtime_error ("Qualification requires a sufficiently long 44.1 or 48 kHz stereo or explicit System I source");
        layout = media.suggestedLayout;
        if (media.channels == 2) { layout.roles[0] = juce::AudioChannelSet::left; layout.roles[1] = juce::AudioChannelSet::right; }
        layout.confirmed = true;
        request = AuditionRequest::source (audioFile, media, layout);
        }
        if ((media.sampleRate != 48000 && media.sampleRate != 44100) || media.durationSeconds() < duration)
            throw std::runtime_error("Qualification programme must cover the requested duration at 44.1 or 48 kHz");
        request.context.buildVersion = LEVELSCOPE_MEASUREMENT_BUILD;
        const auto auditionState = argument ("--audition-state");
        if (auditionState.isNotEmpty() && auditionState != "source" && auditionState != "current")
            throw std::runtime_error ("Audition state must be source or current");
        if (auditionState == "current")
        {
            request.domain = MeasurementDomain::currentState; request.currentInitialized = true;
            request.context.ticket.sessionRevision = 1;
            if(savedProcessing)request.context.processing=*savedProcessing;
            else {
            const auto gainText = argument ("--processing-gain");
            if (gainText.isEmpty()) throw std::runtime_error ("Current qualification requires explicit --processing-gain");
            size_t used = 0; request.context.processing.staticGainDb = std::stod (gainText.toStdString(), &used);
            if (used != gainText.toStdString().size()) throw std::runtime_error ("Invalid processing gain");
            }
        }
        else if (argument ("--processing-gain").isNotEmpty()) throw std::runtime_error ("Source has identity gain only");
        const auto requestError = request.validationError();
        if (requestError.isNotEmpty()) throw std::runtime_error (requestError.toStdString());
        routing.fileToDevice[0] = 0; routing.fileToDevice[1] = 1;
        muted = args.contains ("--muted");
        if ((media.channels == 8 || request.context.assembly) && ! muted) throw std::runtime_error ("Multichannel video qualification requires all monitor routes muted");
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
        audio.synchronize (request, routing);
        processingPresence.expected = media.channels;
        if (! audio.setCaptureObserver (&processingPresence)) throw std::runtime_error ("Cannot attach processing-presence observer");
        const auto preparationError = audio.prepare (request);
        if (preparationError.isNotEmpty()) throw std::runtime_error (preparationError.toStdString());
        player.attach (videoFile); startedNs = attachmentStartedNs = StandaloneAudition::hostTimeNs();
        setSize (800, 520); startTimerHz (120);
        rssSampler = std::thread ([this]
        {
            while (!samplingClosed.load())
            {
                sampleMemory();
                std::this_thread::sleep_for (std::chrono::milliseconds (50));
            }
        });
    }
    ~Run() override { stopTimer(); stopSampling(); audio.disconnect();
#if LEVELSCOPE_VIDEO_PAINT_DIAGNOSTICS
        trace.stop();
#endif
    }
    void paint (juce::Graphics& g) override
    {
#if LEVELSCOPE_VIDEO_PAINT_DIAGNOSTICS
        paintdiag::Trace::Span paintSpan(trace,paintdiag::Kind::paintEnter,paintdiag::Kind::paintExit);
#endif
        const auto paintStart = StandaloneAudition::hostTimeNs();
        const video::Time paintTime (audio.positionFrames(), static_cast<int64_t> (media.sampleRate));
        view = player.update (paintTime, audio.transportGeneration());
#if LEVELSCOPE_VIDEO_PAINT_DIAGNOSTICS
        trace.selection(player,view,paintTime);
#endif
        g.fillAll (juce::Colours::black);
        if (view.picture)
            g.drawImageWithin (view.picture->image, 0, 0, getWidth(), getHeight() - 65, juce::RectanglePlacement::centred);
        if (audio.getState() == AuditionState::playing && view.asset && audio.positionSeconds() >= 2)
        {
            ++paintCount;
            if (view.picture)
            {
                const auto index = view.picture->index;
                if (lastPainted && index > *lastPainted + 1)
                {
                    paintedSkippedFrames += index - *lastPainted - 1;
                    if (skipEvents.size() < 256)
                    {
                        juce::DynamicObject::Ptr event = new juce::DynamicObject();
                        event->setProperty ("previous_frame", static_cast<juce::int64> (*lastPainted));
                        event->setProperty ("frame", static_cast<juce::int64> (index));
                        event->setProperty ("audio_seconds", paintTime.seconds());
                        event->setProperty ("previous_paint_audio_seconds", previousPaintAudioSeconds);
                        event->setProperty ("paint_interval_ms", previousPaintNs ? (paintStart - previousPaintNs) / 1e6 : 0);
                        event->setProperty ("queued_pictures", static_cast<int> (view.queued));
                        event->setProperty ("worker_selection_skips", static_cast<juce::int64> (view.droppedFrames));
                        skipEvents.add (juce::var (event.get()));
                    }
                }
                if (previousPaintNs) maximumPaintIntervalMs = std::max (maximumPaintIntervalMs, (paintStart - previousPaintNs) / 1e6);
                previousPaintNs = paintStart;
                previousPaintAudioSeconds = paintTime.seconds();
                if (!lastPainted || index != *lastPainted) ++uniquePaintedFrames;
                lastPainted = index;
            }
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
    void stopSampling()
    {
        samplingClosed.store (true);
        if (rssSampler.joinable()) rssSampler.join();
    }
    void sampleMemory()
    {
#if JUCE_MAC
            const auto second = static_cast<int> (observedSeconds.load());
            if (second >= 2 && second != lastRssSecond)
            {
                lastRssSecond = second;
                const auto start = StandaloneAudition::hostTimeNs();
                // Read one coherent VM snapshot. Reusable pages remain in RSS but
                // macOS can reclaim them; retain raw RSS and footprint independently.
                task_vm_info_data_t info {}; mach_msg_type_number_t size = TASK_VM_INFO_REV1_COUNT;
                const auto result = task_info (mach_task_self(), TASK_VM_INFO,
                    reinterpret_cast<task_info_t> (&info), &size);
                maximumMemoryQueryMs = std::max (maximumMemoryQueryMs, (StandaloneAudition::hostTimeNs() - start) / 1e6);
                if (result != KERN_SUCCESS || size < TASK_VM_INFO_REV1_COUNT
                    || info.resident_size == 0 || info.phys_footprint == 0 || info.reusable > info.resident_size)
                {
                    ++memoryQueryFailures;
                    return;
                }
                const auto rss = static_cast<double> (info.resident_size);
                const auto reusable = static_cast<double> (info.reusable);
                const auto retained = rss - reusable;
                const auto footprint = static_cast<double> (info.phys_footprint);
                maxPlaybackRss = std::max (maxPlaybackRss, rss);
                maxRetainedRss = std::max (maxRetainedRss, retained);
                maxFootprint = std::max (maxFootprint, footprint);
                if (second < 62)
                {
                    firstRssSum += rss; firstRetainedSum += retained;
                    firstFootprintSum += footprint; ++firstRssCount;
                }
                if (second > duration - 60)
                {
                    lastRssSum += rss; lastRetainedSum += retained;
                    lastFootprintSum += footprint; ++lastRssCount;
                }
                juce::DynamicObject::Ptr sample = new juce::DynamicObject();
                sample->setProperty ("seconds", second);
                sample->setProperty ("rss_bytes", rss);
                sample->setProperty ("reusable_bytes", reusable);
                sample->setProperty ("non_reusable_resident_bytes", retained);
                sample->setProperty ("physical_footprint_bytes", footprint);
                memorySamples.add (juce::var (sample.get()));
            }
#endif
    }
    void finish (juce::String error)
    {
        stopTimer(); stopSampling();
        const auto seconds = audio.positionSeconds(); const auto underruns = audio.underruns();
        const auto firstRssMean = firstRssCount ? firstRssSum / firstRssCount : 0;
        const auto lastRssMean = lastRssCount ? lastRssSum / lastRssCount : 0;
        const auto firstRetainedMean = firstRssCount ? firstRetainedSum / firstRssCount : 0;
        const auto lastRetainedMean = lastRssCount ? lastRetainedSum / lastRssCount : 0;
        const auto firstFootprintMean = firstRssCount ? firstFootprintSum / firstRssCount : 0;
        const auto lastFootprintMean = lastRssCount ? lastFootprintSum / lastRssCount : 0;
        const bool memoryPassed = memoryQueryFailures == 0 && maxPlaybackRss > 0
            && maxPlaybackRss <= 768.0 * 1024 * 1024
            && (duration < 120 || (firstRssCount >= 55 && lastRssCount >= 55
                && lastRetainedMean - firstRetainedMean <= 32.0 * 1024 * 1024
                && lastFootprintMean - firstFootprintMean <= 32.0 * 1024 * 1024));
        const bool hardware = audio.hasHardwareClock(); const auto latency = audio.outputLatencyFrames();
        audio.disconnect();
#if LEVELSCOPE_VIDEO_PAINT_DIAGNOSTICS
        trace.stop();
        const auto traceError=trace.save(traceFile,media.sampleRate);
        if (traceError.isNotEmpty()) error += " Diagnostic: " + traceError;
#endif
        juce::DynamicObject::Ptr report = new juce::DynamicObject();
        auto set = [&] (const char* key, juce::var value) { report->setProperty (key, std::move (value)); };
        set ("scope", "Internal audio-clock and UI paint selection; no optical/acoustic measurement");
        set ("passed", error.isEmpty() && memoryPassed && attachmentPassed && boundFailures == 0 && paintedSkippedFrames == 0 && paintDeadlineFailures == 0 && blankPaints == 0 && underruns == 0 && seconds >= duration && skewCount > 0
            && processingPresence.valid.load() && processingPresence.frames.load() > 0
            && (media.channels != 8 || processingPresence.seen.load() == 255));
        set ("error", error); set ("os", juce::SystemStats::getOperatingSystemName());
#if defined(__arm64__)
        set ("architecture", "arm64");
#else
        set ("architecture", "x86_64");
#endif
        set ("sample_rate", media.sampleRate); set ("device", deviceName); set ("audio_sha256", media.sha256);
        set ("audition_state", request.domain == MeasurementDomain::source ? "Source" : "Current State");
        set ("processing_gain_db", request.context.processing.staticGainDb);
        set ("section_plan_key", request.context.processing.sectionGains?request.context.processing.sectionGains->audioKey():juce::String());
        set ("audition_identity", request.identity());
        set ("measurement_build", LEVELSCOPE_MEASUREMENT_BUILD);
        set ("measurement_status", "UNVERIFIED; timing qualification does not measure compliance");
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
        set ("skipped_frame_limit", 0);
        set ("maximum_paint_interval_ms", maximumPaintIntervalMs);
        set ("skip_events", skipEvents);
        set ("actual_paint_skipped_frames", static_cast<juce::int64> (paintedSkippedFrames));
        set ("unique_painted_frames", static_cast<juce::int64> (uniquePaintedFrames));
        set ("hash_ms", view.asset ? view.asset->hashMilliseconds : 0);
        set ("index_ms", view.asset ? view.asset->indexMilliseconds : 0);
        set ("embedded_audio_policy", "Excluded by video-only backend; independent nonzero-sink audit is a separate test");
        set ("all_routes_muted", muted);
        set ("selected_audio_channels", media.channels);
        set ("programme_assembly", request.context.composition ? "aligned-two-component-sum-v1" : request.context.assembly ? "mono-channel-assembly-v2" : "interleaved");
        if(request.context.composition)set("composition_identity",request.context.composition->definitionIdentity());
        if(request.context.assembly)set("alignment_confirmation",request.context.assembly->alignmentConfirmation);
        set ("measurement_build", LEVELSCOPE_MEASUREMENT_BUILD);
        set ("layout_capability", isSystemI (layout) ? systemICapability : request.context.assembly ? "mono-channel-assembly-v2" : "stereo");
        set ("pre_monitor_frames_observed", static_cast<juce::int64> (processingPresence.frames.load()));
        set ("pre_monitor_nonzero_channel_mask", static_cast<int> (processingPresence.seen.load()));
        set ("pre_monitor_channel_count_valid", processingPresence.valid.load());
        set ("attachment_first_picture_s", attachmentSeconds);
        set ("attachment_budget_s", attachmentBudget);
        set ("attachment_within_budget", attachmentPassed);
        set ("peak_sampled_playback_rss_bytes", maxPlaybackRss);
        set ("first_minute_mean_rss_bytes", firstRssMean);
        set ("last_minute_mean_rss_bytes", lastRssMean);
        set ("memory_growth_assessed", duration >= 120);
        set ("memory_passed", memoryPassed);
        set ("memory_sampling", "Background sampler; no Mach memory query on the paint/timer thread");
        set ("maximum_memory_query_ms", maximumMemoryQueryMs);
#if LEVELSCOPE_VIDEO_PAINT_DIAGNOSTICS
        set ("qualification_eligible", false);
        set ("diagnostic_only", true);
        set ("diagnostic_trace_valid", trace.valid() && traceError.isEmpty());
        set ("diagnostic_trace_error", trace.error());
        set ("diagnostic_trace_file", traceFile.getFullPathName());
#else
        set ("qualification_eligible", true);
#endif
        set ("memory_method", "prores-memory-v2-user-approved-2026-10-20");
        set ("allocator_relief_requested", false);
        set ("memory_query_failures", memoryQueryFailures);
        set ("memory_samples", memorySamples);
        set ("first_minute_memory_sample_count", firstRssCount);
        set ("last_minute_memory_sample_count", lastRssCount);
        set ("peak_non_reusable_resident_bytes", maxRetainedRss);
        set ("peak_physical_footprint_bytes", maxFootprint);
        set ("first_minute_mean_non_reusable_resident_bytes", firstRetainedMean);
        set ("last_minute_mean_non_reusable_resident_bytes", lastRetainedMean);
        set ("first_minute_mean_physical_footprint_bytes", firstFootprintMean);
        set ("last_minute_mean_physical_footprint_bytes", lastFootprintMean);
#if JUCE_MAC
        struct rusage usage {};
        getrusage (RUSAGE_SELF, &usage);
        set ("process_peak_rss_including_load_bytes", static_cast<juce::int64> (usage.ru_maxrss));
#endif
        const auto json = juce::JSON::toString (juce::var (report.get()));
        const bool saved = reportFile.replaceWithText (json + "\n"); std::cout << json << '\n';
        juce::JUCEApplication::getInstance()->setApplicationReturnValue (saved && static_cast<bool> (report->getProperty ("passed")) ? 0 : 1);
        juce::JUCEApplication::quit();
    }
    void timerCallback() override
    {
#if LEVELSCOPE_VIDEO_PAINT_DIAGNOSTICS
        paintdiag::Trace::Span timerSpan(trace,paintdiag::Kind::timerEnter,paintdiag::Kind::timerExit);
#endif
        try
        {
#if LEVELSCOPE_VIDEO_PAINT_DIAGNOSTICS
            if (!trace.start(*this)) { finish(trace.error()); return; }
#endif
            const video::Time timerTime(audio.positionFrames(),static_cast<int64_t>(media.sampleRate));
            view = player.update (timerTime, audio.transportGeneration());
#if LEVELSCOPE_VIDEO_PAINT_DIAGNOSTICS
            trace.selection(player,view,timerTime);
#endif
            if (view.state == video::Player::State::failed) { finish (view.message); return; }
            if (view.picture && attachmentSeconds == 0)
            {
                attachmentSeconds = (StandaloneAudition::hostTimeNs() - attachmentStartedNs) / 1e9;
                attachmentBudget = 30.0 + double (view.asset->bytes) / (64 * 1024 * 1024)
                    + 10.0 * view.asset->frames.size() / 100000.0;
                attachmentPassed = attachmentSeconds <= attachmentBudget;
            }
            if (!opened && !audio.preparing() && audio.preparedCopy() && view.picture && view.queued >= 3 && (!manualStart || requestedStart))
            {
                if (view.asset->frames.back().end < video::Time (static_cast<int64_t> (duration), 1))
                    throw std::runtime_error ("Video does not cover the requested run");
                auto error = audio.connect (media.sampleRate);
                if (error.isEmpty()) error = audio.play (request, routing, true, muted);
                if (error.isNotEmpty()) throw std::runtime_error (error.toStdString()); opened = true;
                playbackStartedNs = StandaloneAudition::hostTimeNs();
                const video::Time startTime(audio.positionFrames(),static_cast<int64_t>(media.sampleRate));
                view = player.update (startTime, audio.transportGeneration());
#if LEVELSCOPE_VIDEO_PAINT_DIAGNOSTICS
                trace.selection(player,view,startTime);
#endif
            }
            audio.poll (player.prerollReady());
#if LEVELSCOPE_VIDEO_PAINT_DIAGNOSTICS
            trace.mark(paintdiag::Kind::repaintRequest);
#endif
            repaint();
            observedSeconds.store (audio.positionSeconds());
            if (audio.getState() == AuditionState::failed) { finish (audio.status()); return; }
            if (opened && audio.positionSeconds() >= duration) { finish ({}); return; }
            if ((!manualStart || requestedStart) && (StandaloneAudition::hostTimeNs() - (opened ? playbackStartedNs : startedNs)) / 1e9
                > (opened ? duration + 120 : 230.0 + double (videoFile.getSize()) / (64 * 1024 * 1024)))
                finish ("Qualification timed out");
        }
        catch (const std::exception& error) { finish (error.what()); }
    }
#if LEVELSCOPE_VIDEO_PAINT_DIAGNOSTICS
    paintdiag::Trace trace;
    juce::File traceFile;
#endif
    std::thread rssSampler;
    std::atomic<bool> samplingClosed { false };
    std::atomic<double> observedSeconds { 0 };
    double maximumMemoryQueryMs = 0;
    int memoryQueryFailures = 0;
    double maxRetainedRss = 0, maxFootprint = 0;
    double firstRetainedSum = 0, lastRetainedSum = 0, firstFootprintSum = 0, lastFootprintSum = 0;
    juce::Array<juce::var> memorySamples;
    ProcessingPresence processingPresence; // Outlives audio and its callback.
    StandaloneAudition audio;
    AuditionRequest request;
    video::Player player;
    video::Player::View view;
    MediaInfo media; ChannelLayout layout; OutputRouting routing;
    juce::File audioFile, videoFile, reportFile;
    juce::String deviceName;
    juce::TextButton startButton;
    bool opened = false, muted = false, manualStart = false, requestedStart = false;
    uint64_t startedNs = 0, paintCount = 0, boundFailures = 0, blankPaints = 0, paintDeadlineFailures = 0;
    double attachmentSeconds = 0, attachmentBudget = 0, maxPlaybackRss = 0, firstRssSum = 0, lastRssSum = 0;
    bool attachmentPassed = false;
    uint64_t playbackStartedNs = 0, attachmentStartedNs = 0;
    int lastRssSecond = -1, firstRssCount = 0, lastRssCount = 0;
    juce::Array<juce::var> skipEvents;
    uint64_t previousPaintNs = 0;
    double maximumPaintIntervalMs = 0, previousPaintAudioSeconds = 0;
    std::optional<size_t> lastPainted;
    uint64_t paintedSkippedFrames = 0, uniquePaintedFrames = 0;
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
#if LEVELSCOPE_VIDEO_PAINT_DIAGNOSTICS
              const auto offline=args.indexOf("--diagnostic-offline-check");
              if (offline>=0)
              {
                  if (offline+1>=args.size()) throw std::runtime_error("Offline check needs a new trace path");
                  paintdiag::Trace trace;
                  trace.mark(paintdiag::Kind::timerEnter);
                  trace.mark(paintdiag::Kind::timerExit);
                  const auto failure=trace.save(juce::File(args[offline+1]),48000);
                  if (failure.isNotEmpty()) throw std::runtime_error(failure.toStdString());
                  std::cout<<"PASS diagnostic startup and trace export; no window or audio device opened\n";
                  setApplicationReturnValue(0); quit(); return;
              }
#endif
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

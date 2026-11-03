# M18 short end-to-end acceptance on Sequoia

Both jobs passed on 2 November 2026 on the development M1 Max, macOS 15.7.9. The same two 60-second, 48 kHz interleaved PCM24 files each ran three times in fresh processes and session directories. The stereo job exercised a predicted global proposal. The 5.1 job exercised a measured section proposal with a crest on LFE. The three repetitions alternated layout order. Files were on local storage with normal OS caching.

The native driver linked the existing Release standalone objects, replacing only the application launcher. It hosted the actual main component in a native JUCE window, used its 10 Hz timer, and exercised the same role review, timeline, gain lane, proposal inspector, session, transport and export handlers as the development app. The driver supplied chooser results, a private role-preferences directory and an injected audio output. It did not automate OS file-panel navigation or measure a hardware output path. No test changed the user's media, sessions or role preferences. The SHA-256 of the existing development executable matched its enabled build manifest. [Build and measured ranges](summary.json) · [Driver and reproduction](../../../Profiling/ShortAcceptance/README.md)

## Measured time

Seconds, median of three. Import starts when the selected audio filename is handed to the app. A fresh job reviews roles and configures one-to-one playback outputs. The driver confirms the profile applicability and roles without a human pause. Play is requested as soon as both measurements are ready. It triggers first-use audition preparation.

| Elapsed time | Stereo | 5.1 |
|---|---:|---:|
| Import to role review | 0.19 | 0.39 |
| Import to first painted partial waveform | 0.61 | 1.23 |
| Import to at least 16 waveform columns containing signal | 0.81 | 1.44 |
| Role confirmation to first painted partial waveform | 0.26 | 0.67 |
| Import to Original and untouched Current ready | 1.44 | 2.17 |
| Play to first non-zero routed sample | 0.30 | 0.71 |
| Import to first non-zero routed sample | 1.74 | 2.88 |
| Propose | 0.019 | 4.00 |
| Independent Current remeasurement | 1.01 | 1.50 |
| Render and independently measure finalized PCM | 1.33 | 2.64 |

The waveform clock stops after the native timeline can paint actual sample extrema. The first painted prefix can be only a few pixels wide, so a second observation requires at least 16 columns with signal. Both programmes begin at a zero crossing. After the normal transport started the device, the driver captured its first 256-frame output callback at 48 kHz. The first non-zero post-monitor sample was frame 1 on every routed channel. The reported time includes its one-frame offset. The device reports zero hardware latency. Analysis and first audition each prepare a verified source, and the UI polls jobs on its normal timer. [Stereo capture](accept-01-stereo/progressive-waveform.png) · [5.1 capture](accept-01-surround/progressive-waveform.png)

## Outcomes

| Result | Stereo | 5.1 |
|---|---:|---:|
| Source IL / true peak | −24.43 LUFS / −10.46 dBTP | −24.72 LUFS / −1.94 dBTP |
| Global gain | +1.431 dB | +1.781 dB |
| Local gain | None | −1.287 dB in the 30–33 s section |
| Private measurement passes during Propose | 0 | 3 |
| Independently verified Current IL / true peak | −23.00 LUFS / −9.03 dBTP | −23.00 LUFS / −1.44 dBTP |
| Measured exported PCM IL / true peak | −23.00 LUFS / −9.03 dBTP | −23.00 LUFS / −1.44 dBTP |
| Current and exported audio, selected EBU R128 rules | PASS | PASS |

Each job confirmed these steps in order:

1. Fresh import showed proposed channel roles. Confirm roles and analyse started the Source pass. The timeline painted incomplete samples before either completed measurement was available. Untouched Current then referenced the exact accepted Source measurement.
2. Play started Original on the explicitly routed output and Space stopped playback. Switching the inspection domain left the audition domain unchanged. Dragging in the main waveform created a section, and dragging the same boundary in Sections adjusted it. Stopped ruler seeking worked in Sections. Zooming and selecting a gain left the playhead in place. Entering and discarding a numerical gain draft preserved both measurements.
3. Propose and Review left Current intact. Apply installed exactly the reviewed processing and disabled verified export until independent remeasurement. Current passed. A stopped main-ruler seek near the section boundary followed by Current playback produced 256 samples per channel identical to the expected applied processing, including the 5.1 gain transition.
4. Saving and reopening restored the same section map, roles, processing plan, complete last-applied explanation and independent Current measurement identity. Playback reopened stopped with the output disconnected and routes muted. The historical explanation opened from the main window.
5. Rendering measured the finalized PCM24 artifact before Publish became available. The explicit Publish action wrote the WAV, manifest and JSON/HTML reports. A separate Python check read the WAV header and checked six or two ordered channels, full frame count, original file hash, bundle hashes, independent Current and artifact verdicts. It regenerated selected PCM24 values from the source, gains, transition positions and recorded dither seed. All 126 stereo and 486 surround sample comparisons matched exactly on every run. [Stereo audit](accept-01-stereo/export-audit.json) · [5.1 audit](accept-01-surround/export-audit.json)

Each stereo run passed 61 driver checks. Each 5.1 run passed 77. Evidence includes all six native and independent audit results, plus one set of artifact reports and captures per layout. Local build outputs retain the playable WAVs and saved sessions in the corresponding acceptance directories.

Two stereo setup pilots stopped early. The first asserted Play button availability before the app's next normal UI refresh. The second confirmed applicability before replacing Source, which intentionally resets applicability for the new programme. A subsequent exploratory pair completed both workflows. Its diagnostic read the baseline global field when recording the stereo prediction. The final runs use the same preview plan as the UI and check exact equality with the plan installed by Apply. The product required no change. [Pilot results](pilots)

This passes the short development-Mac workflow. Actual High Sierra, cross-Mac exchange, selected hardware/video and three-hour/full-duration resource qualification remain pending.

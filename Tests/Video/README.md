# Reference-picture qualification

The backend is a video-track-only AVAssetReader adapter. Tests open no embedded
audio decoder or audio output. Hardware playback is a separate, explicit opt-in
application and is never registered with CTest.

Use `/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation` as the source
root. Configure a fresh `build/video-dev` inside that worktree, following the pinned
external JUCE prerequisite and named-target commands in
[the development workflow](../../Docs/ProjectIndex/19_DevelopmentWorkflow.md#fresh-builds-after-relocation).
The VST3 target installs automatically, so build only the specified standalone/test
targets. The original tested binaries, evidence and old cache are preserved outside
this source root; the old cache is not reusable after relocation. Historical result
records keep their original paths. The commands below are future fixture-generation
examples to run from the new source root, never over the preserved artifacts.

`LevelScopeVideoTests` without arguments checks rational timing, six frame rates
over simulated hour-long timelines, VFR gaps, drop-frame labels, overflow, the
bounded audio presentation clock, transport transitions with an injected device,
and measurement/routing independence. No real audio device is opened.

Optional positional arguments are generated 120-frame movies. The native checks
decode each frame's binary burn-in and compare it to its indexed PTS, seek in both
directions, switch streams, cancel/reject stale generations, clear and close during
work, and check queue bounds, coverage and missing-media errors. `timecode.mov`
additionally requires the generated 01:00:00;00 origin and 30000/1001 rate.
`half-hour.mov` instead checks all 54,000 exact PTS/coverage pairs and random
decoded burn-ins across the long fixture. Compressed indexing uses Core Media's
**output** timing, including container edit-list offsets; raw media timestamps
are not interchangeable with the decoded presentation timeline.

## Generated fixtures

Build the native fixture tool on macOS (universal when transferring to High Sierra):

```sh
xcrun clang++ -std=c++17 -fobjc-arc -arch arm64 -arch x86_64 -mmacosx-version-min=10.13 -Wno-deprecated-declarations -framework AVFoundation -framework AppKit -framework QuartzCore -framework CoreMedia -framework CoreVideo Tests/VideoBackendProbe.mm -o build/LevelScopeVideoFixtureBuilder
build/LevelScopeVideoFixtureBuilder generate build/cfr-exact.mov
build/LevelScopeVideoFixtureBuilder generate-vfr build/vfr-exact.mov
build/LevelScopeVideoFixtureBuilder generate-tc build/timecode.mov
build/LevelScopeVideoFixtureBuilder generate-tc build/half-hour.mov 54000
build/LevelScopeVideoFixtureBuilder remux-mp4 build/cfr-exact.mov build/cfr-exact.mp4
build/LevelScopeVideoFixtureBuilder remux-mp4 build/cfr-exact.mov build/cfr-exact.m4v
python3 Tests/generate_video_sync_audio.py build/half-hour-clicks.wav
python3 Tests/generate_video_sync_audio.py build/half-hour-clicks-44100.wav --sample-rate 44100
```

The 320x180 H.264 movie uses a 30000-unit timescale, 1001 units per CFR frame,
54,000 frames (1801.8 seconds), a readable frame/PTS burn-in and binary frame bits.
Every thirtieth frame has a bright red/orange background flash. Paired stereo WAVs
have a 5 ms, 1 kHz pulse at each flash, with peak amplitude 0.1 before monitor gain.
Clicks are placed from exact rational event times, rounded to the nearest source
sample; rounding does not accumulate. The 48 kHz interval is exactly 48,048 samples.
The VFR fixture doubles every third frame interval. The writer assigns its last
sample one base-frame duration, leaving a final gap before the session end; the
adapter preserves this gap instead of extending the last picture. Timecode fixture generation
uses 29.97 drop-frame only. Generated media and binaries stay outside Git.

`mux input.mov input.wav output.mov` creates a fixture with two identical video
tracks and one audio track; the production decoder only reads the selected video
track. Identical streams test selection and indexing, not perceptual discrimination.

## Coordinated playback

After coordinating the exact ordinary device, rate, volume and whether the routes
are muted with the user, run the qualification app's executable with absolute paths:

```sh
"LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification" --allow-audio --muted --device "Built-in Output" --audio /absolute/half-hour-clicks-44100.wav --video /absolute/half-hour.mov --seconds 1800 --report /absolute/new-report.json
```

Omit `--muted` only for explicitly coordinated audible playback. The harness uses
no inputs, opens the named two-output device at the fixture's rate, sets its own
monitor gain to -30 dB and disconnects at the end. It does not choose a default
device. Match the fixture to the device's current rate to avoid changing it.
The normal application's all-muted playback protection remains enabled; only the
harness opts into silent transport with all routes at -1. A headphone check must
be observed on the remote Mac's physical display, not its Screen Sharing image.
Add `--manual-start` for a listening check: the picture is shown while audio stays
closed until the Start button is clicked. Bring the window to the intended monitor
and coordinate headphone readiness before clicking. The audio generator optionally
accepts `--peak-amplitude 0.4`, 12.04 dB above its default 0.1; coordinate this change
before playback as well. The test app does not alter system volume.

The report compares the frame PTS selected for each UI paint with the frame expected
from the predicted audio presentation clock, after two seconds of preroll. Its
one-frame bound uses the expected frame's actual duration, including VFR. It records
maximum skew, early/late means, blank paints, skipped frames and audio underruns.
Selection reads the clock at paint time, avoiding a stale timer snapshot across
VFR boundaries. CPU drawing duration is measured against the expected frame's
duration, with deadline misses counted separately. The expected-frame calculation
shares the tested rational mapping; native burn-in tests independently verify that
decoded image content agrees with that index. Zero selection skew is not evidence
of zero physical display delay. Late polls include preroll and generation changes
and count polls without the desired decoded picture; they are not dropped paints.
These are internal scheduler measurements. They do not measure scanout, display
processing, DAC/acoustic delay, or physical lip sync. Closing the window early
ends the run without a pass report. A 30-minute run must complete in real time.

Keep failed and superseded attempts identifiable. Bind successful reports to
binary, source/build and fixture hashes and actual OS/architecture. Preserve the
M1 evidence and compare fresh official ITU/EBU results against its pinned reports.

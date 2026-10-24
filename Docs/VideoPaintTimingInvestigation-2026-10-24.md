# Video paint timing investigation — 24 October 2026

This is a read-only investigation and a proposed next diagnostic, not a fix or a
new qualification result. M7 candidate-02 remains frozen, unqualified and
unpublished. No playback, build, installation, source integration or audio-device
opening was performed. The prior qualification records remain historical records
of their exact runs; this analysis does not replace them.

The evidence supports an intermittent gap between successful picture selection
and UI painting on Sequoia. It does not yet identify which scheduling stage or
session condition caused that gap. M6 stereo reproduced the failure, so M7's
eight-channel processing is not necessary to trigger it. Neither rebooting nor
the later passing comparisons establishes a permanent fix.

## Repository and identities

- Video worktree: `/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation`,
  branch `codex/video-implementation`, unchanged at
  `a7fec0d6e7e8ef5da7f56819722d73d96dd10f78`.
- Published main inspected after fetching:
  `b57a1ed23cbe7d74afe98ed16ec8fd615db9d171`, qualified M6. Main is two commits
  beyond the video worktree. Inspection used Git objects and preserved package
  sources; no branch was advanced. The accepted workflow-file rename is preserved.
- Original checkout: `/Users/fuzz/Developer/LevelScope`, branch
  `codex/logic-module`, HEAD `b57a1ed23cbe7d74afe98ed16ec8fd615db9d171`, with
  uncommitted M7 work and evidence. Nothing there was edited.
- Pinned dependency: JUCE 7.0.12, commit
  `4f43011b96eb0636104cb3e433894cda98243626`, at
  `/Users/fuzz/Developer/LevelScope-dependencies/JUCE-7.0.12`; unchanged.

| Identity | Preserved M6 candidate-01 | Frozen M7 candidate-02 |
| --- | --- | --- |
| Package | `/Users/fuzz/Developer/LevelScope/build/m6-qualification-2026-10-23/candidate-01` | `/Users/fuzz/Developer/LevelScope/build/m7-qualification-2026-10-24/candidate-02` |
| Source digest | `0702768309e7328549fb86b677ec295fdc197d0382c45f0487b9060053b9a286` | `c3df888f567a69adb83cbb6c7d2c7a1dcdf1da9498d43d688ffe890786c0e43e` |
| identity.json SHA-256 | `ab7440f4d06275ea55055f59eec98c4f5872a775342e5b31777cefdf1dfdb4dd` | `ca3eda67d2ab4505237075a0f9ff41968ff88f0d5e615652891f886c7a94d59b` |
| Video qualification executable SHA-256 | `c8efc1a6d0b760206c8f99105e553542a4d32fdb9722dc3e19b2172f60bb36e9` | `7b0113714b90fa582dd6ae725652630e7bd4f601fecc608fcb4dd15acd0d512e` |
| Manifest artifacts checked | 308 | 339 |

The video executable is `LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification`
inside each package. Golem's existing M7 package at
`/Users/fuzz/LevelScope-M7-qualification-2026-10-24/candidate-02` was checked over
SSH without running it: its identity and harness hashes match those above. Golem
reported Intel x86_64, High Sierra 10.13.6 build 17G66. This remote check covers
those two files only, not a fresh full-package or playback qualification.

The actual reference movie remains
`/Users/fuzz/Developer/LevelScope/Tests/51_SyncTest.mov`, SHA-256
`0870d5ab93d1492f23524ed1fff9319af1a6a9cf71882543f2190929afeb8e24`:
1920×1080 ProRes 422, progressive, 25 fps, 2772 frames, 110.88 seconds,
375714116 bytes. The stereo selected-audio fixture is
`/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/prores-20261020/fixtures/selected-clock-48000.wav`.
The M7 eight-channel comparison uses its package's
`m7-fixtures/clock-71-48000.flac`. Exact fixture hashes for every run are retained
in the raw reports and the companion analysis JSON. The movie is not added to Git.

## Reconciled evidence

Evidence root:
`/Users/fuzz/Developer/LevelScope/Tests/Reference/results/m7-qualification-2026-10-24/candidate-02`.
The following are nine playback groups, not an inventory of all M7 tests.
Only full raw reports are counted; `*-outcome.json`, `runs.json` and summaries
are not additional runs. There are 73 raw reports, 68 passing and five failing,
with seven missed painted frames in seven events.

| Group relative to evidence root | Runs | Passed | Missed painted frames | Largest paint-entry gap (ms) |
| --- | ---: | ---: | ---: | ---: |
| `sequoia-muted-video-first-attempt` | 16 | 14 | 3 | 51.340 |
| `paint-investigation-01` | 4 | 4 | 0 | 27.840 |
| `sequoia-muted-video-controlled` | 16 | 15 | 1 | 52.120 |
| `paint-reference-investigation-01` | 2 | 0 | 3 | 51.125 |
| `high-sierra-muted-video` | 16 | 16 | 0 | 21.494 |
| `post-reboot-local-six-01` | 6 | 6 | 0 | 38.317 |
| `screen-sharing-local-six-01` | 6 | 6 | 0 | 36.520 |
| `screen-sharing-idle-ssh-local-six-01` | 6 | 6 | 0 | 38.280 |
| `settings-change-m7-single-01` | 1 | 1 | 0 | 25.227 |

The first failures were Proxy 29.97 (two misses) and Proxy 30 (one). The four
short M6/M7 stereo/eight-channel comparisons then passed. In the controlled
matrix, all fifteen generated files passed but the actual reference movie missed
one frame. The two 100-second stereo reference runs failed with both preserved M6
(two misses) and M7 (one). High Sierra's sixteen-file matrix passed. Its device
ran at 44.1 kHz while the Sequoia failures used Mac Studio Speakers at 48 kHz;
the hosts therefore are not a comparison that changes only the OS.

After reboot, each six-run set alternated three M6 and three M7 stereo runs of
100 seconds. All eighteen passed: disconnected; Screen Sharing connected to idle
Golem; and Screen Sharing plus idle SSH. These are six runs per group, not twelve.
Idle connections alone are not supported as the cause under these conditions.
The passing groups are not replacements for the earlier failures.

One subsequent M7 run passed after the user changed settings. The later,
separately preserved clarification says Firefox Local Network permission was
granted and browser tabs were opened. It does not establish that Firefox was
absent before, what network/browser/GPU workload ran during any test, or why the
permission reportedly changed after reboot. No causal browser or permission
diagnosis follows from this single pass.

Post-reboot snapshots agree on the same boot session, PB248 main display at
1920×1200/60 Hz, no mirroring, and recorded power settings. They contain selected
connection/process snapshots, not continuous WindowServer, browser, GPU, UI
occlusion, display-link or event-loop histories. The earlier environment record
identifies OS/hardware/binaries but lacks an equivalent complete session timeline.
Present-day inspection cannot reconstruct those missing observations.

## What the code and counters establish

1. **Successful selection can occur without a paint.**
   `Source/Standalone/StandaloneVideo.cpp:31–59` implements `Player::update`.
   Both the harness timer and paint callback call it. The function updates
   `lastDisplayed` when a desired decoded picture is found, and increments
   `view.droppedFrames` if that successful selection jumps over an index.
   The worker does not update `lastDisplayed`. The report field named
   `worker_selection_skips` is therefore a shared selection counter, not a
   decoder-worker skip counter.

   Every missed-paint event has that selection counter at zero. In these
   continuous, single-start runs, this strongly implies that a timer update
   selected the intermediate frame even though no qualified paint recorded it.
   It narrows the failure to the path after an intervening successful selection
   and before the next paint. It is stronger evidence than the queue depth of
   seven or eight at the resumed paint alone. Full generation/view-epoch and
   per-update traces are still needed to exclude unrecorded reset histories
   directly and prove each handoff.

2. **The skipped frames coincide with approximately 49–51 ms paint-entry gaps.**
   Across the seven events, selected-audio advancement differs from the measured
   host-time gap by at most 4.542 microseconds, less than one 48 kHz sample
   (20.833 microseconds). This supports a real paint-callback gap rather than a
   sudden audio-position jump. It is an internal timing comparison, not a
   physical output-latency or optical/acoustic synchronization measurement.

   Three late actual-movie events cluster at source times 56.327, 56.842 and
   57.726 seconds, but at different frame indices. There are also early failures
   around 2.51, 3.57, 3.60 and 5.55 seconds. Neither one bad source frame nor a
   periodic 60-second cause is established.

3. **The CPU-paint metric is narrower than the complete drawing path.**
   In preserved M6 `Source/Tests/VideoPlaybackQualification.cpp:88–149`,
   `maximum_cpu_paint_ms` is sampled before the status text is drawn and before
   control returns through JUCE/AppKit/CoreAnimation. The M7 equivalent has the
   same measurement boundary. Zero CPU-paint deadline failures therefore does
   not exclude slow status drawing, native drawing/submission or delayed native
   callbacks. Memory sampling is already off the UI thread; it should not be
   blamed on synchronous UI memory queries. Allocation of a skip-event object
   happens after the triggering interval and cannot explain the first gap.

4. **JUCE already uses display-link-driven repaint dispatch.**
   The pinned `modules/juce_gui_basics/native/juce_NSViewComponentPeer_mac.mm`
   queues dirty rectangles in `repaint` (around line 1068), receives CVDisplayLink
   notifications through an AsyncUpdater on the main thread (around line 1755),
   and calls `setNeedsDisplayInRect` from `onVBlank` (around line 1083). Async
   notifications can be coalesced. Live resize can throttle requests, but the
   evidence does not record live resize and does not prove it occurred.

   The application's 120 Hz JUCE timer is not a presentation guarantee.
   `VBlankAttachment` would receive notifications after the same main-thread
   handoff. Replacing the timer with that listener alone would not bypass all
   existing scheduling or native draw delays.

5. **JUCE's optional Metal-layer renderer is absent from these harnesses.**
   The M6/M7 build flags do not enable
   `JUCE_COREGRAPHICS_RENDER_WITH_MULTIPLE_PAINT_CALLS`. Both preserved arm64
   binaries have no `CoreGraphicsMetalLayerRenderer` symbols and their
   `NSViewComponentPeer::displayLayer` implementations are a single `ret`
   instruction (M6 address `0x1000092c8`, M7 `0x1000096c0`). The compile-time
   guard and binary inspection exclude that optional renderer's command-buffer
   deferral path as the cause here. This does not exclude operating-system
   compositor/GPU scheduling below the asynchronous CoreGraphics layer.

The common video implementation is byte-identical in the M6 package, M7 package
and inspected main: `StandaloneVideo.cpp`, `StandaloneVideoPanel.cpp`,
`StandaloneVideo_mac.mm`, and `AudioPresentationClock.h`. The M7 audio and harness
have additional multichannel/processing observations, but M6's reproduction
shows those additions are not required. The product panel shares timer/paint
selection logic while doing different UI work; harness callback observations do
not prove the product's physical presentation cadence.

## Next discriminating diagnostic

Prepare a fresh diagnostic candidate from qualified main/M6
`b57a1ed23cbe7d74afe98ed16ec8fd615db9d171` in the video worktree, after comparing
its implementation inputs to preserved M6. Preserve the original executable as
the control. Give the diagnostic source, binary and trace schema new identities.
Do not import M7 or overwrite any existing candidate/build. Build only the named
`LevelScopeVideoQualification` target in a fresh build directory.

Start with a preallocated, bounded event trace recording:

- timer entry/exit, successful selection, repaint request, paint entry and full
  callback completion, all using the same monotonic clock;
- audio position read on the UI thread, selected/painted frame, generation and
  view epoch, desired frame readiness and queue indices;
- independently timestamped display-link notifications and the corresponding
  JUCE main-thread vblank delivery, with counts to expose coalescing;
- display identity/rate, window foreground, visibility, occlusion and live-resize
  transitions, plus full callback duration including status text.

Keep trace recording allocation-free and bounded, outside the audio callback;
serialize after playback closes. Report overflow explicitly so an incomplete
trace cannot masquerade as an explanation. Retain existing counters unchanged
and add accurately named counters rather than rewriting old reports. If these
observations leave the delay between vblank handling and paint unresolved, add a
second narrowly scoped native draw/invalidation hook in the diagnostic build.
Any dependency patch must be reviewable and isolated; never edit shared JUCE or
the preserved packages. Native draw completion still must not be described as
actual display scanout.

The proposed first playback comparison is a fixed six-run Sequoia sequence,
alternating preserved M6 stereo and the fresh diagnostic candidate, three runs
each, 100 seconds per run, with the same reference movie, selected-audio fixture,
Mac Studio Speakers at the recorded 48 kHz, muted routes, display geometry and
hands-off interval. Coordinate device use and that interval before starting.
Record session/workload conditions instead of assuming they match. Do not reboot
or change settings during the comparison. A separately coordinated Golem
control/diagnostic comparison should retain its existing device configuration and
record its actual rate; do not change the rate merely to match the Studio.

Interpret the trace before choosing a correction:

| Trace observation | Next bounded action |
| --- | --- |
| Display ticks continue, main-thread delivery stalls | Identify the long callback/event-loop work; bound or move that work off the UI thread. |
| Main-thread vblank and repaint requests continue, native paint does not | Trace invalidation/coalescing, visibility and native draw entry; evaluate a window-specific scheduling change only if supported. |
| A paint starts promptly but full callback/native return is long | Separate picture draw, status draw and native submission costs; correct the measured slow stage. |
| Desired frame was unavailable or generations changed unexpectedly | Revisit decode readiness/lifecycle with the recorded history; do not infer readiness from the resumed queue alone. |
| No failure reproduces | Retain the pass as a bounded diagnostic sample. Do not claim a fix; predefine the next single workload/session change. |

A later source-time versus elapsed-time comparison can start near the reference
movie's 50-second point to investigate the 56–58-second cluster, if it remains
relevant. Normal interaction stress should be separately labelled and preserve
audio continuity and recovery. Neither is a substitute for the unchanged
controlled zero-skipped-frame gate. Do not force late frames to be displayed by
delaying audio or accumulating picture lag.

## Correction acceptance and disposition

No production correction is selected yet. Evidence is insufficient to justify a
timer-frequency increase, a blanket display-link replacement, a JUCE upgrade,
or disabling a renderer that these binaries do not use.

After a justified correction, qualify a new candidate on both actual Macs:
all accepted 1080p ProRes profiles/rates, the actual movie, H.264 CFR/VFR,
Source/Current, lifecycle/seek/offset/frame identity, measurement and selected-audio
isolation, product UI, and three 30-minute ProRes runs per Mac (one per profile).
Retain failed attempts and exact source/binary/fixture identities. Require zero
actual painted-frame skips after the existing two-second preroll, zero blanks,
frame-bound violations, CPU-paint deadline failures and audio underruns. Preserve
the existing attachment/cancellation/queue gates and approved memory gate:
768 MiB total playback RSS ceiling, with non-reusable RSS and physical-footprint
last-minus-first minute mean growth each at most 32 MiB. Add more complete timing
observations without silently replacing or loosening existing limits.

Preserve exact rational timestamps, signed offsets, audio-authoritative
selection, transport generations, source identity, measurement behaviour and
monitor separation. Embedded movie audio remains undecoded and unrouted.
Physical optical/acoustic synchronization remains unmeasured. No VST3, installed
application, Pro40, driver, service or system-audio changes belong to this work.

M7's prior passing device-closed/native/reference evidence remains preserved;
its outstanding physical video groups are still incomplete. No M7 checkpoint,
integration or publication follows from this investigation. If a correction is
later accepted, integrate it deliberately into a newly identified M7 candidate
and rerun affected checks on both platforms. Physical 7.1 listening remains a
separate pending observation.

The companion file
`Tests/Video/results/paint-investigation-2026-10-24/read-only-analysis.json`
records all 73 raw-report paths/hashes, individual skip arithmetic, manifests,
binary-renderer observations, package verification, and preservation checks.
It is derived evidence only. The original historical paths and records have not
been rewritten. No message was sent to another chat.

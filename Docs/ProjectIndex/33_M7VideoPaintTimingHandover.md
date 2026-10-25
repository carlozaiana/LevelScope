# Video handover — reproduced UI paint skips; M7 remains frozen

24 October 2026. The user explicitly chose **prepare a Video Implementation
handover and keep M7 frozen**. No display-scheduling change has been made here.
Do not publish, merge, rebuild over, or relabel M7 candidate-02 as qualified.

## Starting state and ownership

- Published main and logic-branch HEAD:
  `b57a1ed23cbe7d74afe98ed16ec8fd615db9d171`, qualified M6.
- M7 implementation is uncommitted in the established logic working copy
  `/Users/fuzz/Developer/LevelScope`, branch `codex/logic-module`. Preserve it.
- Registered video working copy:
  `/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation`, branch
  `codex/video-implementation`, observed HEAD
  `a7fec0d6e7e8ef5da7f56819722d73d96dd10f78`. Inspect its own state and newer remote
  main before choosing a development base. Do not reset it or copy M7's unrelated
  changes into it. This chat has not modified that working copy.
- All ten existing checkpoint tag objects and target commits were verified
  unchanged locally and remotely. No M7 checkpoint exists.

Read [32 — M7 qualification](32_M7QualificationAndHandover.md),
[30 — M6 handover](30_M6QualificationAndHandover.md), and the accepted video/ProRes
contracts in documents 18, 25 and 27. M7's approved PCM scope remains in [31](31_M7ChannelBasedPcm71Proposal.md).

Frozen M7 candidate-02 source digest:
`c3df888f567a69adb83cbb6c7d2c7a1dcdf1da9498d43d688ffe890786c0e43e`.
Its exact 155 input hashes and 339 artifact hashes are in
[`identity.json`](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/identity.json).
The source-manifest ancestry commit is not the uncommitted candidate's source identity.

Preserved packages:

- Sequoia: `build/m7-qualification-2026-10-24/candidate-02` under the logic root.
- High Sierra: `/Users/fuzz/LevelScope-M7-qualification-2026-10-24/candidate-02`.
- Baseline M6 on Sequoia: `build/m6-qualification-2026-10-23/candidate-01`.

## Reproduced problem

The strict gate is **zero skipped UI-painted frames**, independently of audio
continuity and frame-selection error. The selected-audio clock remained
authoritative in every failed case: zero audio underruns, blank paints, CPU-paint
deadline failures and frame-selection bound failures. These are internal clock/UI
measurements, not optical/acoustic synchronization measurements.

On arm64 Sequoia 15.7.9 (24G830), Apple M1 Max / Mac13,1, Mac Studio Speakers at
48 kHz, all output routes muted:

| Run | Result |
| --- | --- |
| First M7 eight-channel Source ProRes matrix | 14/16 passed; Proxy 30000/1001 skipped 2 frames, Proxy 30 skipped 1 |
| One controlled repeat, with explicit hands-off notice | All 15 generated cases passed; `51_SyncTest.mov` skipped 1 frame during its 100-second run |
| Bounded Proxy comparison | M6 stereo, M7 stereo, M7 eight-channel 30000/1001 and M7 eight-channel 30 each passed once |
| Same reference movie, preserved M6 binary with stereo audio | Failed: 2 skipped UI frames |
| Same reference movie, M7 binary with stereo audio | Failed: 1 skipped UI frame |

Thus eight-channel processing is not required to reproduce the issue. The cause
is **not established**. The user was unsure about activity during the first
matrix; do not attribute it to user interference. The controlled repeat and both
reference comparisons retain the same gates and unchanged binaries. Original
failed reports have not been overwritten or replaced by passing comparisons.
Earlier M6 qualification records are also unchanged; these are new observations.

Actual Intel High Sierra 10.13.6 (17G66), using Built-in Output at 44.1 kHz with
all routes muted, completed all sixteen Source matrix cases successfully,
including the reference movie, with zero skipped frames. Current-state matrix,
H.264 physical checks and three 30-minute endurance runs per Mac remain deferred.
All running qualification processes have exited and their devices are closed.

Useful trace points from the reference-movie runs:

- M7 eight-channel: frame 1419 → 1421 at audio 56.842292 s; 50.023709 ms between
  paints; seven decoded pictures queued; no worker-selection skips.
- M6 stereo: frame 87 → 89 at 3.566063 s, and 1441 → 1443 at 57.726083 s;
  51.125417 / 48.753792 ms gaps; seven pictures queued.
- M7 stereo: frame 1406 → 1408 at 56.327146 s; 50.289375 ms gap;
  eight pictures queued.

The late events cluster around programme time 56–58 seconds across these three
runs. Investigate this correlation; it does not establish a content-dependent
cause. Queue readiness, measured CPU paint duration, OS/UI dispatch and actual
presentation are different measurements.

Both Macs' existing display-sleep and system-sleep settings were zero; they were
read, not modified. The qualification window calls JUCE's foreground activation.
Neither observation proves that every dispatch/paint was timely.

Evidence is preserved under
[`candidate-02`](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/):
`sequoia-muted-video-first-attempt`, `sequoia-muted-video-controlled`,
`paint-investigation-01`, `paint-reference-investigation-01`, and
[`video-paint-investigation.md`](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/video-paint-investigation.md).
Each run records exact commands, media identities, binary hashes and events.
The final [qualification summary](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/qualification-summary.json)
links the independent passing evidence, failed video audit and final integrity
records. It explicitly records M7 as unqualified and frozen.

## Minimal reproduction

### Later isolation result: six passes after reboot

At the user's request, six fixed local stereo reference tests ran on 24 October
2026, approximately 17:06–17:17 Europe/Rome, after the Mac Studio was rebooted and
Golem disconnected. Screen Sharing was closed; local preflight/final checks found
no Screen Sharing or SSH client process. No remote connection was made.
Preserved M6 and frozen M7 alternated for three runs each, 100 seconds per run.
**All six passed with zero skipped UI frames, underruns, blank paints, clock-bound
failures or CPU-paint deadline failures.** The total elapsed test time was
10 minutes 38 seconds. Source/binary/package integrity and local tags remained
unchanged, and all test processes exited with devices closed.

The [post-reboot addendum](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/post-reboot-local-six-01/README.md)
contains exact commands, environment, hashes and every result. This is evidence
that the failure did not reproduce in that fixed sample. Reboot and removal of
remote-session workload changed multiple conditions together: the cause is not
isolated, and earlier failures remain valid observations. Keep M7 frozen; its
deferred physical matrix/endurance checks and overall qualification are unchanged.

### Connected Screen Sharing comparison: another six passes

The user then connected Screen Sharing to Golem's idle desktop, kept SSH
disconnected and left both Macs untouched. A second fixed six-run set completed
approximately 17:27–17:38 Europe/Rome in the same Mac Studio boot session. The
same display/power configuration, binaries, source/audio/video identities,
alternating run order, 100-second duration, muted 48 kHz output and limits were
retained. **All three M6 and three M7 runs passed**, again with zero skipped UI
frames, underruns, blank paints, clock-bound failures or CPU-paint deadline
failures. Summed executable elapsed time was 10 minutes 38 seconds.

The [connected comparison addendum](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/screen-sharing-local-six-01/README.md)
records the six results and matched comparison with the disconnected baseline.
Fourteen local boundary observations found an established Screen Sharing TCP
connection and no SSH client; no remote command or continuous network capture was
used. Final integrity passed and all test apps exited with their device closed.
Screen Sharing was left connected for the user.

Both fixed samples passed 6/6. This does not support blaming an idle Screen
Sharing connection alone. Active remote-screen content, SSH and concurrent
remote workloads were not exercised in this comparison. The cause remains
unresolved; preserve these passes alongside the original failures and keep M7
frozen. No additional qualification group was started by that comparison.

### Screen Sharing plus idle SSH: another six passes

At the user's request, a third fixed six-run set completed approximately
17:49–18:00 Europe/Rome in the same boot session, with Screen Sharing and an idle
SSH shell connected to Golem. **All three M6 and three M7 runs passed**, with zero
skipped frames, underruns, blank paints, clock-bound failures or CPU-paint
deadline failures. The same recorded configuration, files, run order and gates
were retained. Summed executable elapsed time was 10 minutes 38 seconds.

The [idle-SSH comparison addendum](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/screen-sharing-idle-ssh-local-six-01/README.md)
contains raw results and a matched three-condition comparison. Fourteen local
boundary observations confirmed established Screen Sharing and SSH connections
to matching destination hosts. Idle remote activity is user-confirmed; no remote
commands were issued and no continuous network/activity trace was taken. Test
apps exited with their device closed, both connections were left open, and final
source/package/tag integrity passed.

All eighteen post-reboot 100-second stereo runs passed: six disconnected, six
with Screen Sharing only and six with Screen Sharing plus idle SSH. These samples
do not support idle connections alone as the explanation under the tested
conditions. Active remote work, changing remote-screen content and other
session/load conditions remain untested. The original cause is unresolved;
earlier failures and M7's frozen, overall unqualified status remain unchanged.

### One further M7 run after user-reported settings changes

The user then requested one test after changing Mac Studio settings.
Frozen M7 completed another 100-second muted stereo reference run, approximately
18:23:52–18:25:39 Europe/Rome, in 106.928 seconds. It passed with zero skipped
frames, underruns, blank paints, clock-bound failures or CPU-paint deadline
failures. Source/package/tag integrity passed and the test device closed.

The [single-run addendum](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/settings-change-m7-single-01/README.md)
retains current observations separately. The recorded display/power configuration
and boot session matched the preceding condition; those snapshots do not cover
all settings. In a later [clarification](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/settings-change-m7-single-01/user-settings-clarification-2026-10-24.md),
the user identified granting Firefox Local Network access and opening browser
tabs before that run. They report that this access becomes blocked again after
each reboot and that they need it for local router/receiver signal checks. That
permission-reset behavior is user-reported, not independently diagnosed here.
Screen Sharing and SSH remained connected. Exact browser/network workload was
not measured; do not assume Firefox was absent from earlier runs or infer a
cause/fix from this single pass. It is not another completed six-run condition.
M7 remains frozen and unqualified overall. Original records are preserved with
the clarification added separately.

### Command used in the reference comparison

The reference movie SHA-256 is
`0870d5ab93d1492f23524ed1fff9319af1a6a9cf71882543f2190929afeb8e24`.
It remains at `/Users/fuzz/Developer/LevelScope/Tests/51_SyncTest.mov`; exclude it
from Git. The stereo clock comes from the preserved ProRes fixture folder.

Preserved M6 video-qualification executable SHA-256:
`c8efc1a6d0b760206c8f99105e553542a4d32fdb9722dc3e19b2172f60bb36e9`.
M7 candidate-02 executable SHA-256:
`7b0113714b90fa582dd6ae725652630e7bd4f601fecc608fcb4dd15acd0d512e`.

After coordinating a new run and giving a clear hands-off notice, use a fresh
report path. This example opens only the named output and keeps every route muted:

```sh
"/Users/fuzz/Developer/LevelScope/build/m6-qualification-2026-10-23/candidate-01/LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification" \
  --allow-audio --muted --device "Mac Studio Speakers" --seconds 100 \
  --audio "/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/prores-20261020/fixtures/selected-clock-48000.wav" \
  --video "/Users/fuzz/Developer/LevelScope/Tests/51_SyncTest.mov" \
  --report "/absolute/new-test-folder/m6-reference.json" --audition-state source
```

Use the M7 executable above for the same stereo comparison. To reproduce its
eight-channel run, use candidate-02's `m7-fixtures/clock-71-48000.flac` instead.
Do not modify preserved packages or run repeated trials until one happens to pass.
New instrumentation belongs in a fresh build/package with its own identity.

## Recommended bounded investigation

Inspect `StandaloneVideoPanel.cpp`, `StandaloneVideo.cpp`,
`StandaloneVideo_mac.mm`, and `Tests/VideoPlaybackQualification.cpp` in the
appropriate working copy. Product and harness currently poll through a 120 Hz
JUCE message-thread Timer and request repaint; selection is refreshed at paint.
Pinned JUCE 7.0.12's `juce_Timer.h` documents approximate timing and message-queue
delays. This makes timer/paint dispatch a hypothesis to instrument, not a proven
root cause or a reason to weaken the gate.

1. Capture bounded timer, repaint-request, paint-start/end, foreground/occlusion,
   display-refresh, queue and host-time observations. Separate decoder readiness,
   UI scheduling, draw CPU time and compositor presentation. Keep instrumentation
   off the audio callback and avoid perturbing the measurements with synchronous
   logging or memory queries on the UI thread.
2. Account for all three six-run post-reboot samples: disconnected, Screen Sharing
   only, and Screen Sharing plus idle SSH, all passed. Avoid attributing the
   failure to idle connections alone. Consider a
   separately coordinated fixed comparison that changes one additional
   session/workload condition at a time, with the same binaries and gates.
   If instrumentation is needed, compare preserved M6 stereo with a fresh
   instrumented candidate. Investigate the 56–58 s interval and control display
   conditions explicitly. Add a separate normal-interaction stress case: it must
   preserve audio continuity and clock selection/recovery. Record visual misses
   honestly; do not merge its interpretation with the unchanged controlled gate.
3. Propose a bounded correction after identifying evidence for the cause. A
   display-refresh-driven path is an option to evaluate, not an adopted design.
   Keep actual High Sierra support and the pinned dependency in view.
4. After a fix, qualify a new video candidate on both actual Macs, retaining
   source hashes, binary hashes, failed attempts and passing results. Include
   H.264 CFR/VFR, all accepted ProRes rates/profiles, Source/Current, the reference
   movie, and the three 30-minute ProRes runs per Mac with existing paint/clock,
   underrun, attachment, queue and memory limits. Do not silently change limits.

Preserve audio presentation time, rational seeks/offsets, transport generations,
stop/resume/restart/EOF, programme duration and source identity. Embedded video
audio stays undecoded and unrouted. Preserve measurement identities, monitor
separation and full-precision processing. No installed app/plugin, Pro40, driver,
service or system-audio setting changes. Build only explicitly named standalone
or test targets; no VST3/default/all build or installation.

## Return path to M7

M7's device-closed offline, injected eight-channel routing/sample parity, actual
three-hour runs, official references, M1–M6 regressions, native video and silent
UI evidence have passed. Full physical 7.1 listening remains separately pending
and does not cause this blocker. Required physical video groups are incomplete.

Keep candidate-02 and its reports immutable. Return the diagnosis, exact tested
source/build identities, evidence and any justified video correction to the logic
work. If code changes, integrate deliberately into a new M7 candidate, compare
its implementation inputs with candidate-02, and rerun affected qualification on
both actual platforms before considering an M7 checkpoint. If code is unchanged,
the remaining required checks still need completion and the intermittent failure
needs an explicit disposition. No merge/publication is authorized by this handover.

## Starter message for the Video Implementation chat

> Please investigate the reproduced UI paint skips documented in ProjectIndex
> 33, keeping M7 candidate-02 and the logic working copy frozen. The same reference
> movie previously failed the strict paint-skip gate with preserved M6 stereo and
> M7 stereo/eight-channel playback, while audio continuity and frame-selection
> bounds passed. After rebooting the Studio and leaving Golem, Screen Sharing and
> SSH disconnected, a fixed six-run comparison passed all three M6 and three M7
> stereo runs. A second fixed six-run set in the same boot session also passed
> with Screen Sharing connected to idle Golem and SSH disconnected. A third fixed
> six-run set also passed with Screen Sharing and idle SSH connected. All eighteen
> post-reboot runs passed; idle connections alone are not supported as the cause
> under these conditions, and no permanent fix is established.
> One further M7 run also passed after the user granted Firefox Local Network
> access and opened browser tabs; details and observation limits are in document 33.
> Investigate the remaining session/workload differences and preserve earlier
> failures alongside these passes.
> Inspect the existing video working copy and latest main, preserve all pending
> work/tags/packages, and identify the cause before proposing a bounded correction.
> Preserve the audio-authoritative timing contract, actual High Sierra/Sequoia
> support and unchanged acceptance limits. Coordinate physical runs and announce
> hands-off intervals; keep outputs muted unless separately authorized. Do not
> change M7 source, installations, plugins, Pro40 or system audio configuration.

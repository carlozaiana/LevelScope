# Background playback preparation — UI handover, 23 November 2026

This UI development continues from published main
`18adeb520006b7fe7e0a22300cc26abd28eec51c`, retaining shared-source service commits
`93df14e1466e96ff073c5fd1440cc84291b3119e` and
`c9b624e688b982526e60d661b8917b312871e83a` and the final export lease guard
`e54ee83be3d0fb9d697e0b4a40e9a4a1d15d6a91` as ancestors. Logic owns the combined
universal freeze and final actual-platform qualification. These development
receipts are not a substitute for that frozen-package run. Main, tags, frozen
packages, manuals, original programmes and earlier evidence remain unchanged.

## Product behavior

After initial Source analysis and its progressive overview complete, the main
window schedules device-free playback preparation. It prepares one independent
reader, the exact selected Source/Current processor and at most 4096 processed
frames (98,304 bytes at the supported 6-channel layout; the engine can also use
12 channels, 196,608 bytes). Scratch decoding uses 2048 frames; FLAC seek pre-roll
is bounded to 65536 discarded frames and is not retained as PCM. No full decoded
programme is cached by this UI change. Shared verified bytes remain governed by
the service store. Gain transitions continue through the same processor after
the prefix; monitor gain, smoothing and output routing remain downstream.

Preparation opens no output and cannot authorize playback. The main feedback
shows its stage, progress, cancellation and separate Playback ready state.
Explicit Connect is allowed while background work remains stopped. Explicit
Play joins the preparation and becomes cancellable delayed transport intent.
Stop, Cancel, Disconnect, stale processing/source/roles/routing, New and Open
revoke that intent. Picture preroll still gates audio output. Stopped seeking
cancels the old-position prefix and schedules the new one. Inspection zoom does
not seek or play. Analysis, proposal and export verification remain independent. The shell no longer
retains a separate export source lease: rendering obtains its reader through the
active session context, while each staged artifact retains its own required
source and final-artifact checks.

Reader trust binds the request to both session runtime and weak store owner,
plus original and sealed-copy identity/stamps. An identical-byte relink can
replace the owner within the same runtime; both automatic scheduling and the
engine account for that. Cancellation/failure is attempted only once for the
same request/position/owner; explicit Play retries. New inputs or a changed
position/processing/owner may schedule new preparation. Cancelled or expired
workers cannot restore a former position or start sound.

## Build and run

Configure the existing JUCE 7.0.12 build with deployment target 10.13 and:

```text
LEVELSCOPE_BUILD_REGION_UI_TESTS=ON
LEVELSCOPE_BUILD_AUDITION_TESTS=ON
LEVELSCOPE_BUILD_DEVICE_SETTINGS_TESTS=ON
LEVELSCOPE_BUILD_TIMELINE_PROGRESS_TESTS=ON
LEVELSCOPE_BUILD_PLAYBACK_TIMING_TESTS=ON
```

Required targets and device-free commands, each using a fresh absolute output:

```sh
LevelScopeRegionUITests ABS_NEW_PLAYBACK_DIRECTORY --playback-ui
LevelScopeAuditionTests --report ABS_NEW_AUDITION_JSON
LevelScopeDeviceSettingsTests ABS_NEW_DEVICE_DIRECTORY
LevelScopeTimelineProgressTests --root ABS_NEW_TIMELINE_DIRECTORY --output ABS_NEW_TIMELINE_JSON
LevelScopeRegionUITests ABS_NEW_UI_DIRECTORY
LevelScopePlaybackTiming ABS_WAV ABS_NEW_TIMING_JSON
```

The immutable-candidate runner `Tests/Qualification/M20/run_ui.py` also accepts
`--case playback-ui`, `--case devices`, `--case timeline` and `--case ui`.
Audition must use its ordinary no-long-source path for this affected regression.
No supplied programme is required for these native/regression groups. They
synthesize owned fixtures and use injected outputs. Logic separately includes
no-argument VideoTests for short clock/preroll continuity, not movie qualification.

## Timing scope

`PlaybackPreparationTiming.cpp` compiles against both the qualified frozen
baseline and integrated source. Baseline alone installs a counted diagnostic
writer. Integrated uses default storage and records store counts after accepted
Source analysis, after warmup and after first/subsequent Play. It reports copy
and byte deltas as well as totals. A custom storage writer would bypass the
normal session cache and is deliberately absent from the integrated path.

These are device-free engine timings with approximately 1 ms polling, measured
to the first nonzero sample delivered to the injected app output. They are not
physical-output measurements or full native UI click latency (the normal main
poller runs at 10 Hz). Source validation/analysis and automatic preparation are
reported separately. Do not add those phases to first Play or hide them.

Unchanged inputs for the coordinated stereo and 5.1 comparison:

- `/Users/fuzz/N15/Mixdown/MusicInBackAndForegroundSomeDialog.wav`
  (`7ba6404041c94c909df290a6829a0fe92ce594329a1d575490c2bb4a5b6cc8f0`)
- `/Users/fuzz/N15/Mixdown/5.1_Busy.wav`
  (`4fa30bfd9a8fbb64e68fc00eaf11f8adf81ead8502d05ce639b22cc07645d7b6`)

Baseline receipts and original-link-input identity are preserved under
`build/background-playback-preparation-2026-11-23/baseline-02` in the UI worktree.
The baseline is implementation `79b55eb177a01e56fdbe86f843e03140df356bcd`,
measurement build `9b379f0479f69a6e2ded23dc778ad5c50c5e92a20e4076a5936e06b5a82698fd`.
Stereo first/subsequent Play: 0.912446833 / 0.002609916 s, one extra 79,375,076-byte
copy. 5.1: 1.844278667 / 0.002630458 s, one extra 160,617,696-byte copy.
Final integrated comparison belongs to the combined package/timing slot.

## Development evidence and limits

Receipts in `results/development-sequoia-2026-11-23` identify the exact development
measurement build and affected results. Earlier native development passed 629
checks; the intermediate expanded run passed 657. Final-source development passed
**7,934 assertions** on Sequoia arm64: 668 playback/native, 4,069 audition, 450
device settings, 2,580 timeline, and 167 existing native proposal/UI. Implementation
`b618ce9dfa77a89b56ba50fd209e70e6eca0305b`, Source tree
`6a099bb9325326cf4b6136f3854e54dbc258ce77`, measurement build
`67009a8ea02c5cb442e2fd9a209de9c868fc932c1e63488dbb17fcf1ba46b5b2`.
See `final-summary.json` for executable hashes and exact receipts. The original
stopped/repaint assertion passes with the Clear fix. Final native screenshots
were visually checked for both readiness and preparation. An early portion of
the existing UI run may overlap Logic compilation; all unchanged timing gates
passed. The timeline run started after Logic released the Mac. The first ordinary audition rerun
stopped at the old assertion forbidding even explicit Connect during background
preparation. That expectation was updated to the approved contract: no implicit
open; explicit Connect stays stopped; pending explicit Play still blocks device
replacement. The failed receipt is retained. No production fix was made to
satisfy that obsolete assertion. The ordinary native UI receipt also retained its
first failure: its test helper classified disconnected cancellation of a
background worker as active transport. Review then found a real Clear Source
conflict: that same cancelling state could
make the model reject the requested clear. Transport activity now ends when
Stop quiesces callbacks and revokes Play; preparation lifetime remains a separate
guard. The original stopped helper and immediate-repaint assertions are retained,
and the new native test explicitly clears a programme during background work.
A source-identity test label now correctly
names the inode/ctime protection rather than claiming hashing was necessary.

Final universal package and actual High Sierra execution remain the combined
owner's work. No physical device was opened and no new listening acceptance is
claimed. Broad video, full-film and three-hour endurance remain deferred, with
no automatic 90-minute substitute. The existing completion pause is separate.
Limited acceptance of occasional dropped video frames does not waive audio
interruptions, accumulating drift or sustained freezes.

## Short manual check after the combined build is available

1. Import stereo or 5.1, confirm roles and let initial analysis finish. Observe
   playback preparation and Playback ready without opening Audio setup.
2. With your explicitly chosen output/routes and safe monitor volume, press
   Play. If preparation is still underway, Cancel playback or Stop must prevent
   later sound. Press Play again to retry deliberately.
3. Stop, change position, switch Original/Current, or edit a gain/transition.
   Playback uses that exact selection; Current still needs remeasurement after
   a real processing edit. Inspection zoom must not move the transport.
4. Open another session or choose New during background preparation. The old
   programme must not begin playing or restore its former position. Remembered
   connection preference and explicit Disconnect behavior remain in force.

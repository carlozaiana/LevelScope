# Video implementation handover

Prepared 17 October 2026. This document is a handover specification; no video
implementation is included in the measurement milestone.

## Integration continuation (17 October 2026)

Measurement has been integrated into main. The final-candidate 30-minute gate now
passes on both actual systems without source changes or rebuilding. See
[the current acceptance and main handover](20_MainIntegrationAcceptance.md).
Earlier checkpoint-only instructions below describe their historical stage; the
user has now authorized staged main integration. The fixed video tag stays intact.

## Current video checkpoint and source location (17 October 2026)

Continue video work at
`/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation`, on
`codex/video-implementation`. Git moved this registered worktree out of the build
area; `Source/`, `Tests/`, `Docs/` and all untracked files moved intact. The original
`/Users/fuzz/Developer/LevelScope` remains on `codex/logic-module`.

The published [video branch](https://github.com/carlozaiana/LevelScope/tree/codex/video-implementation)
has fixed annotated tag
[`video-checkpoint-2026-10-17`](https://github.com/carlozaiana/LevelScope/releases/tag/video-checkpoint-2026-10-17)
at `700b9e77e439a8ca9bb7583adf826fac01dbe924`. Its implementation is unchanged from
`54e2aa9d2773229dd666bb36df410ff932e9f9bf`. Location/workflow documentation follows
in a separate commit without moving the tag.

Tested products, generated fixtures, raw evidence and old CMake caches remain under
`/Users/fuzz/Developer/LevelScope/build/video-implementation/{build,probe,evidence}`.
A verified binary/evidence snapshot is in that root's
`preserved/video-checkpoint-2026-10-17/`. The old cache is historical and cannot be
reused after relocation; future builds use a fresh directory inside the new source
root. Read [the development workflow](19_DevelopmentWorkflow.md) for exact artifact
locations, hashes, verification, dependency setup and future build commands.

Qualification remains bounded: final candidate 06 has the numeric/native tests and
short real-device runs; preceding harness 05 has the 30-minute and listening runs.
Physical lip-sync latency is not measured. Embedded soundtrack playback remains
unimplemented. Preserve the historical test records and their recorded paths.

## Video development continuation (17 October 2026)

Implementation and qualification evidence are recorded in
[18_VideoQualificationResults.md](18_VideoQualificationResults.md). The specification
below remains the original milestone contract; the continuation describes the
implemented scope and its measurement limits.

The measurement tag was verified locally and with `git ls-remote origin`: annotated
tag object `76affb667b6dab066cb36b2a2c741525b9b6eeac`, peeled commit
`b29bc7c183743f35c6271d2d7ec8fb5ba1464591`. Development branch
`codex/video-implementation` started at that commit in the former isolated checkout
`build/video-implementation/checkout`. That was the source path during qualification;
it has since moved to the location above. Historical build, generated media and raw
evidence remain in `build/video-implementation/{build,probe,evidence}` under the
original project. The original checkout remains on `codex/logic-module`; M1 evidence
is preserved.

Backend spike: JUCE is still 7.0.12 at
`4f43011b96eb0636104cb3e433894cda98243626`. Its `VideoComponent` wraps `AVPlayer`,
exposes seconds as doubles, and does not expose the required indexed PTS, explicit
video-track choice or generation-tagged seek completion. The chosen macOS adapter
uses `AVAssetReaderTrackOutput` for exactly one video track. It creates no player,
audio reader/output, audio engine or audio device. Decoded BGRA pixels are copied
on the video worker into immutable software images; the message thread displays them.
No third-party codec or runtime installation is required.

The universal `Tests/VideoBackendProbe.mm` spike ran on actual arm64 Sequoia 15.7.9
and Intel High Sierra 10.13.6. Both decoded all 120 CFR and VFR H.264 frames, passed
random-start decoding, and submitted 120 frames to a native window layer. Probe
SHA-256: `e3d1ab8b5e5e7035f7e4692ba0814ee918e4186df872a7085344dbd11e1a456b`.
CFR fixture: `ca0f64031788361e61c5c80eea7226aa87861fef2b52a18de4f0ae6f243cfad8`;
VFR: `c8e48b93056d2848b8cc5d36b851ba4724c821ad895a8fa1359d697b57e16712`.
These are API/decode/presentation-submission observations, **not** playback-sync
qualification or optical/acoustic measurements. Later playback and regression
results are recorded separately in document 18.

The spike found that default AVAssetWriter timescales quantize 30000/1001 material;
the generator therefore pins its track/movie timescale to 30000. AVAssetReader also
clips the first returned PTS when a time-range seek falls inside a frame. The adapter
indexes compressed **output** PTS/durations (sorted from decode order), then seeks to
an exact indexed PTS. It never infers VFR positions from nominal rate or frame count.
Output timing is essential: a passthrough-remuxed MP4 shifted raw compressed PTS
by two frames while its decoded presentation started at zero. Empty marker buffers
are skipped and batched samples are expanded using their output timing.
See Apple's [track-output settings](https://developer.apple.com/documentation/avfoundation/avassetreadertrackoutput/outputsettings)
and [timecode reader/writer note](https://developer.apple.com/library/archive/technotes/tn2310/_index.html).

Initial format scope: local MOV/MP4/M4V, SDR H.264, square pixels/full aperture,
maximum 3840x2160 decoded raster; at most two million indexed frames and eight
queued pictures. Other codecs, HDR and unsupported aperture formats fail explicitly.
VFR uses half-open PTS coverage with no frame before/after coverage or in an indexed
gap. Single-origin 32/64-bit timecode tracks carry their exact rate, start and drop
flag; unavailable or multiple-origin metadata is labelled, never invented.

The audio clock uses callback-submitted source spans, not decoder position. The
pinned JUCE CoreAudio callback passes HAL's output `mHostTime` (nanoseconds). That
timestamp already includes buffering/scheduling; only read-only queried device and
stream hardware latency is added. When that mapping is unavailable, callback entry
time plus JUCE's total reported output latency is used and labelled an estimate.
The latter is not added to the HAL timestamp. The bounded span history preserves
source position during inserted underrun silence and drains presentation at EOF.
Stop retains the predicted presented source frame; resume starts there. Seek/restart
quiesce the callback, flush queues, change generation and preroll before resuming.
Audio and device latency values are internal predictions, not physical measurements.

Reference attachment, stream choice and exact signed offset are independent session
state. Positive offset means `videoTime = audioTime - offset`; time quantization
floors towards negative infinity. Rational operations reject overflow. An optional
embedded-soundtrack source switch was discussed as a follow-up; this milestone keeps
all embedded audio excluded. For physical testing, pair burned-in flash markers with
selected-audio clicks and measure screen/speaker output independently.

For subsequent builds, use the fresh-directory commands in
[the development workflow](19_DevelopmentWorkflow.md#fresh-builds-after-relocation)
from `/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation`. Build only the
named standalone and test targets, including video tests and the separate playback
harness. Historical build commands and cache paths are retained with the evidence;
do not execute them against the preserved products.

The separate Golem SSH handle did not authenticate noninteractively. The verified
existing `/tmp/pro40-high-sierra-ssh` handle was used without modifying or closing it.
Remote probes live under `~/LevelScope-video-2026-10-17/probe`; installed applications
and the M1 qualification bundle are untouched. The user authorized silent GUI checks
on both machines. Subsequent explicitly coordinated tests used Mac Studio Speakers
at 48 kHz and Golem's Built-in Output at its user-selected 44.1 kHz headphone rate.
Golem's endurance run has all routes muted. A 20-second headphone check was heard
and seen by the user, who reported that clicks and flashes seemed aligned.

## Starting checkpoint and qualification gate

Start from annotated tag **`measurement-milestone-2026-10-17`** on
`codex/logic-module`. Resolve its exact starting commit with
`git rev-parse measurement-milestone-2026-10-17^{commit}`. Its tested implementation
is `b3fe9dece6b9d4cd2d8f6a3aebca69b8d0534e9c`; subsequent checkpoint changes contain
qualification evidence and documentation only. Actual Intel High Sierra 10.13.6
(17G66) and arm64 Sequoia 15.7.9 (24G830) passed all M1 numeric and silent UI checks.
Read [17_M1QualificationResults.md](17_M1QualificationResults.md) for the complete
scope, artifact hashes, comparisons and observed checklists.
See [the measurement contract](14_OfflineMeasurementContract.md),
[the milestone handover](13_StandaloneScaffoldHandover.md), and
`Tests/Reference/README.md` for evidence and remaining product limitations.

Preserve `main` at `556d387ad772dc5e604ead95f1525f60e66fd5f9` and the annotated tags
`high-sierra-pre-scaffold-2026-10-15` and `standalone-scaffold-2026-10-15`.
Do not change the Pro40 repository, installed drivers/services, system audio
configuration or the installed app/plugin. Use a separate checkout/build/test folder.

## First video milestone: synchronized reference picture

Attach a local reference-video file, select its video stream, show its metadata and
picture, and play it in synchronization with the selected audio. Provide signed
offset adjustment with the convention below. Selected audio remains authoritative;
embedded video audio remains muted throughout loading, seeking and playback.
The first milestone includes reliable stop/restart/seek behavior and explicit missing,
unsupported or failed-video states. It excludes editing, video rendering/export,
embedded-audio mixing, delivery profiles and polished UI.

## Timing and transport contract

Carry forward section “Video and transport contract” in `14_OfflineMeasurementContract.md`:

- Keep session UUID and separate audio/video content identities, stream ID and
  transport generation. Attaching video or editing its offset never changes audio
  measurement identity, Current State measurement or file-to-speaker/device routing.
- Audio timeline positions use signed 64-bit source frames. Video PTS/timebase and
  nominal frame rate use reduced integer rationals. Preserve VFR timestamps; never
  derive VFR presentation time from frame count.
- `videoTime = audioTimelineTime - videoStartOnAudio`. Positive offset starts video
  later on the audio timeline. Use checked rational conversion and specify rounding.
- Timecode carries its start origin and drop-frame flag. Drop-frame changes labels,
  not playback speed. Support signed offsets without overflow or truncation drift.
- Drive picture from the audio presentation clock, accounting for output latency.
  Audit whether a counter represents submitted or presented frames before applying
  latency; do not compensate twice. Worker decode position is not a playback clock.
- Stop retains position; restart seeks to zero. Seek stops/flushes old queued media,
  increments generation, prerolls and resumes only if transport was already playing.
  Audio EOF stops transport; display no frame outside the video's coverage.
- Video decode and UI display stay outside the real-time audio callback. Reject
  stale completions after seek, replacement, clear or close. Buffering/error policy
  must preserve the audio clock and expose dropped/late video frames.

## Existing interfaces and integration boundaries

- `Source/Standalone/StandaloneSessionModel.h`: source/current documents, session UUID,
  revisions, processing gain, measurement identity and invalidation. Add reference
  video as independent session state. Session persistence has not been implemented.
- `StandaloneAudition.h/.cpp`: current output-only audio transport, decoder worker,
  bounded SPSC buffer, sample counters, lifecycle and device ownership. It currently
  supports play-from-start/stop/EOF, without pause/seek/resampling or presented-clock
  telemetry. Add a transport abstraction and latency contract before synchronization.
- `StandaloneAudioPanel.h/.cpp`: explicit device opening, role confirmation, muted
  routes and monitor gain. Preserve explicit selection and hardware independence.
- `StandaloneServices.h`, `StandaloneJob.h`, `StandaloneMeasurementService.h/.cpp`:
  worker-only value requests, cancellation, ticket/revision checks and measured truth.
  Video jobs need equivalent cancellation/generation discipline without changing M1.
- `StandaloneMainComponent.h/.cpp`: temporary workflow UI and worker polling. Keep
  the reference-picture view and transport controls separate from measurement state.
- `Source/Core/Measurement/OfflineMeasurement.h/.cpp`: qualified offline measurement;
  no host/UI/video dependencies. Video must not alter its numerical input or lifecycle.

## Compatibility and backend decision

Support actual Intel macOS High Sierra 10.13 and arm64 macOS Sequoia 15. Preserve the
universal build and deployment minimums. Inspect the pinned JUCE 7.0.12 facilities and
Apple APIs available on both systems before selecting JUCE video or a platform adapter.
Prove a small decode/presentation spike on both operating systems. Record supported
containers/codecs explicitly; fail unsupported media clearly. No minimum-OS increase,
codec/runtime installation or system audio change is implicit in this milestone.

## Synchronization acceptance tests

Use deterministic frame-number/burn-in reference movies paired with known audio-frame
events. Record backend, codec, rate, timebase, hardware/display, audio latency and hashes.

1. Rational-clock tests: 24, 25, 30, 24000/1001, 30000/1001 and 60000/1001 rates;
   drop-frame minute/ten-minute boundaries; VFR; signed offsets; large timestamps;
   rounding and overflow; at least one-hour simulated timelines without accumulating drift.
2. Offline frame selection: expected PTS/frame at exact and adjacent boundaries,
   offset changes and coverage limits. Decode/preroll never publishes an old generation.
3. Real playback on both OS versions: after preroll, picture selection stays within
   one source-frame duration of the predicted audio presentation time; record measured
   skew and drift over a 30-minute run, including fractional-rate material. For VFR,
   use the duration of the expected frame. State separately what display scheduling
   was observed and what acoustic/optical output was measured; do not infer physical
   lip-sync accuracy from internal timestamps alone.
4. Seek while stopped/playing, repeated seeks, stop/restart, positive/negative offsets,
   source/video replacement, mismatched durations, EOF, decode failure and close during
   load/decode. No stale picture, stuck worker, unbounded buffering or audio discontinuity
   caused by video work. Quantify dropped/late frames under load.
5. Embedded audio is inaudible/absent from the audio graph; attaching video leaves
   selected-audio samples, routing, gain and M1 measurements unchanged. Re-run M1 suites.
6. Coordinate any audible/hardware synchronization checks with the user first. Use
   ordinary output devices; Pro40 is not a dependency. Silent tests keep device closed,
   routes muted and transport stopped.

The one-frame bound above is a proposed product acceptance requirement, not a claim
of an EBU/ITU video-sync standard. Agree and document any tighter delivery requirement
before treating it as a compliance criterion.

## Build and tests

Use the reproducible build/test commands in `13_StandaloneScaffoldHandover.md` and
`Tests/Reference/README.md`. Build only `LevelScopeStandalone`,
`LevelScopeMeasurementTests` and `LevelScopeScaffoldTests`; the VST3 target has automatic
installation enabled. Downloaded reference audio and build products stay outside Git.
Preserve the pinned manifests. Run actual High Sierra tests, not just the Intel slice
under Sequoia. VST3 measurement parity is a separate milestone in
[16_Vst3MeasurementParityMilestone.md](16_Vst3MeasurementParityMilestone.md).

## Starter message for the dedicated chat

> Continue LevelScope from `measurement-milestone-2026-10-17` on `codex/logic-module`.
> Read ProjectIndex documents 14, 15 and 17 before changing code. M1 measurement is
> qualified on actual Intel High Sierra 10.13.6 and arm64 Sequoia 15.7.9. Implement the
> first video milestone: synchronized reference-picture playback and signed offset
> adjustment, with selected audio authoritative and embedded video audio muted.
> Carry forward the rational timing, audio presentation-clock and transport-generation
> contract. First inspect backend compatibility on both OS versions and define the
> measurable synchronization tests in document 15. Preserve M1 results, the original
> main branch and baseline tags, installed apps/plugins, and all Pro40/system audio
> configuration. Build only named standalone/test targets. Keep VST3 parity separate
> (document 16), and coordinate any audible tests before opening an audio device.

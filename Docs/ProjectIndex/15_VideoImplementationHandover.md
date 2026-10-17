# Video implementation handover

Prepared 17 October 2026. This document is a handover specification; no video
implementation is included in the measurement milestone.

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

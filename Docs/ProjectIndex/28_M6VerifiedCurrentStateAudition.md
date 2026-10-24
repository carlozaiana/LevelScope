# M6 proposal — verified Current State audition

Prepared 22 October 2026. **Approved 23 October; candidate 01 implemented and qualified
on both actual platforms, published and integrated into main on 24 October 2026.**
The user authorized implementation, qualification and subsequent checkpoint publication.
Physical-device opening and audible testing require advance coordination.
This document also records the requested future channel-label
clarification and the user's separate informal playback observation.

Candidate 01 implementation, exact source/binary identity and current evidence are
recorded in [30_M6QualificationAndHandover.md](30_M6QualificationAndHandover.md).
Both actual-platform automated suites, interactive checks and coordinated physical
runs passed. The linked handover records complete evidence and the annotated checkpoint.

## Baseline and dependencies

Reviewed published main `6d1c22a23244b4e421dd2773862451ac8a863992`, containing qualified
M5 `ae88e187d288c53dbca1287e8341b74b9b5b2619`, the fixed ProRes checkpoint
`dcc34aabba6b5b3ca4805e2965d21ffebf2e4d0a`, and accepted ProRes integration
`af7ed751ab869f7829c65e7736939e45425c0d54`. The established logic branch was
fast-forwarded to this main without conflicts; the video worktree stays untouched.
The latest main commit only renames `.github/workflows/build.yml` to `buildNOTyml`.
That file no longer has an Actions workflow extension; this plan preserves that
published change and does not assume automatic CI qualification.

Required contracts: [14](14_OfflineMeasurementContract.md),
[24](24_M5ProgrammeSections.md), [27](27_ProResMainIntegration.md), the
[ProRes acceptance criteria](25_ProResAcceptance.md), and
[ProposalEngineSpec_v2.md](../ProposalEngineSpec_v2.md). M1 owns measurement truth,
M2 owns measured-rule evaluation, M3 owns explicit global-trim application, and
M4/M5 own read-only diagnostics and in-memory section annotations.

At the reviewed baseline, `StandaloneAudition` provided bounded streaming, explicit
device opening, exact-rate playback, monitor gain/routing, transport generations
and the audio presentation clock. It played Source only, reopening the original
file and checking size/timestamp. `StandaloneMeasurementService` instead verified
a private source snapshot by SHA-256 and applies the typed `static-gain-v1` plan.
M6 closed that content-identity gap and added the shared gain operation.

## Bounded scope and user-visible meaning

Add an explicit Source / Current State audition selector. Display the selected
state, processing trim (Source always 0 dB; Current uses the full stored double),
separate monitor gain, preparation/playback status, and measurement freshness:
current, unavailable/not measured, or stale/unverified after a processing edit.
Missing individual metrics keep their existing reasons, including silence and
short-programme limitations; they do not become zero. M2's actual verdict stays
separate and keeps its narrow wording. Playback never implies a compliance pass.

Approved policy: a valid, initialized Current State may be auditioned before
remeasurement, visibly **UNVERIFIED**. The request must still bind the current
source and processing plan. A stale playback request is always rejected; an old
passing measurement or an unapplied M3 prediction cannot authorize or describe the
new audio. Audition does not measure, generate/apply proposals, normalize, retry a
failed verification or roll back gain. Explicit existing measurement remains the
way to establish measured Current State and M2 results.

The qualified Current State matrix is M1's mono C/L/R, stereo L/R, conventional
5.0/5.1 with accepted surround aliases, at 44100/48000/96000/192000 Hz; preserve
ordered file channels, including LFE. Reject unsupported Current processing plans,
layouts and rates explicitly. Preserve broader existing Source audition support,
with measurement-unavailable wording outside M1's matrix; do not widen measurement
claims or silently substitute Source for an unavailable Current State.

## Approved architecture

Implementation locations: `Source/Standalone/StandaloneProcessing.h` owns the
sample operation; `VerifiedSource.h/.cpp` owns cancellable copy/hash/payload
validation; `StandaloneAuditionRequest.h` validates and binds domain/source/roles/
plan/revision; `StandaloneAudition` prepares and streams; `StandaloneSessionModel`
enforces stopped-only edits; `StandaloneAudioPanel` exposes selection, preparation,
trim and measurement freshness. The label displays 12 decimal places while the
stored double and sample operation retain all precision; labels never feed the plan.
`Tests/AuditionTests.cpp` exercises
the actual consumed stream, injected device lifecycle, storage failures and product
control callbacks. The separate video harness records audition state and gain.

Source hashing uses bounded 64 KiB input buffers. The algorithm, full-content hash,
payload checks and cancellation remain unchanged; buffering avoids one filesystem
read per hash block on long programmes. M1 reference suites are rerun because
measurement uses this shared utility. The source-sized disk copy is separate from
the duration-independent memory budget. The current device/plan settings and valid
session are retained on preparation failure, with playback stopped.

1. A small standalone processing primitive validates the existing typed plan and
   computes `g = pow(10.0, gainDb / 20.0)` in double precision once per request.
   For every channel/sample it produces `float(double(x) * g)`, exactly matching
   the qualified measurement operation. Source uses identity processing. Retain
   the existing -60..+24 dB range, full gain precision, no-fast-math/FP-contract-off
   build policy and non-finite rejection. Do not route this through plugin DSP,
   loudness channel weights, dynamics masks or monitor controls. LFE gets the same
   trim even though its loudness energy weight is zero.
2. Extract the current bounded-copy/hash/payload-validation logic into an owned,
   cancellable source-snapshot utility, preserving measurement behavior. Playback
   decodes only the verified private copy. Retain one source copy for the active
   imported content so stop/resume, seek and gain/domain changes do not recopy a
   three-hour file. Rebind its consumer identity for ordered-role/range changes.
   The copy is ephemeral, read-only to consumers and released on source teardown;
   no persistent cache or processed-file render is introduced. Failure/cancel
   cleans owned partial files; abrupt process termination retains the existing
   temporary-file-cleanup limitation. Record preparation progress and disk cost.
3. Create an immutable audition request binding session UUID, content hash,
   metadata, ordered confirmed roles, full programme range, Source/Current domain,
   exact typed plan and relevant processing revision. Bind each prepared/queued
   result to its request and transport generation. Profile, applicability, monitor
   gain, routing, video offset and section annotation revision stay outside the
   audio processing identity. Routing has a separate stopped-only configuration.
4. Keep the existing decoder worker, bounded FIFO and presentation clock. The
   worker decodes a fixed-size block, applies the shared operation, checks finite
   samples and feeds the FIFO. A test-only preallocated capture observes frames
   actually consumed immediately before monitor smoothing and physical routing.
   Production callback remains free of allocation, file I/O, hashing, locks,
   decoding, measurement and notifications. There is no duration-sized PCM cache.

Signal order is: verified source bytes -> float32 decode -> identity/global trim
-> bounded FIFO -> test capture point -> monitor ramp/gain -> one-to-one routing
-> explicit device output. The callback's consumed source spans remain the clock.
Monitoring can change what reaches speakers without changing the processing plan,
offline samples, measurements or evaluations.

An owned cached copy represents the imported content, not a claim that a mutable
external path still contains those bytes. Same-size/timestamp mutation before
preparation must fail the content hash. Later external writes cannot enter the
private copy; a detected source-stamp change stops audition and requires reimport.

## Lifecycle and timing

- Preparation is possible with the device closed. Opening a device and starting
  playback remain explicit actions. No device opens on selection, import, prepare,
  measurement, application launch or recovery. Preparation/cancellation have their
  own generation and cannot publish after source replacement or shutdown.
- Require Stop before changing domain, roles, processing gain, applying an M3
  proposal, reinitializing/clearing Current State or editing routing. Disable UI
  actions and enforce the same rule in the action/controller boundary; also reject
  stale engine requests. An unexpected programmatic invalidation stops output and
  flushes the old generation before new samples can be accepted. Re-selecting the
  same state or setting an identical gain is a no-op.
- Switching Source/Current while stopped retains the held source-frame position,
  flushes/rebinds playback and requires explicit Play. It does not change gain,
  processing revision, valid measurements, M2 verdicts or M3 drafts. Gain edits
  retain position and invalidate affected Current measurements as today. No change
  is made to the selected Source samples.
- Import/replace/clear first stops and cancels audition. Successful source
  replacement resets coordinates and routes as today; role changes require
  confirmation and invalidate dependent Current State. Failed import retains the
  previous valid source under existing transactional behavior, but stays stopped.
  Device error/disconnect stops, flushes and does not automatically reopen/resume.
- Seek retains exact source-frame rounding, flushes audio/video generations,
  prerolls and resumes only if previously playing. Stop retains presented position;
  Play resumes; Restart begins at zero. Natural EOF drains presented audio and
  stops; seek-to-EOF and subsequent Play retain the existing explicit behavior.
  Underflow silence does not advance the source clock beyond consumed samples.
- Source and Current have identical frame counts, source sample rate and duration.
  Static gain adds no latency, lookahead or tail. Measurement FIR/LRA padding is
  never played. Preserve rational PTS/timecode/VFR and
  `video time = audio time - videoStartOnAudio`, including negative offsets,
  outside-coverage blanking and audio-authoritative EOF. Embedded video audio
  remains undecoded and unrouted, including in ProRes MOV files.

## Approved over-full-scale policy

Finite processed float values above absolute 1 are valid pre-monitor samples and
must remain unchanged for capture/measurement parity. Do not clamp, limit, dither,
normalize, reduce the plan's gain or add a compliance allowance. Such a Current
State may truthfully fail its measured true-peak rule.

Use a separate monitor-output fault: if the actual routed, monitor-scaled
callback would contain a non-finite value or a sample strictly outside [-1,1],
output silence for that callback, latch a clear monitor-overload error and stop.
Do not present the suppressed span as successfully played or automatically retry.
Report the offending channel/source frame and leave processing/monitor settings
unchanged. The user explicitly reduces monitor gain or edits the plan before
restarting. Sufficient manual monitor attenuation permits audition of finite
over-full-scale processing while preserving its pre-monitor samples.

This is an output fault policy, not processing correction or true-peak protection.
A sample-safe output does not establish analog/intersample-peak safety, complete
delivery compliance or calibrated listening level. Tests must exercise the actual
monitor ramp, including transitions from a previously louder setting.

## Acceptance tests

| Area | Required acceptance |
| --- | --- |
| Sample parity | Capture actual FIFO-consumed samples before monitor controls. Zero gain is bit-identical to decoded Source, including initial/final samples. Nonzero audition samples are bit-identical to the offline processing stream on the same platform; compare independently generated scalar/high-precision reference cases as well. No display rounding enters either path. |
| Gain and channels | Test 0, positive, negative, -60/+24 dB and non-rounded M3-derived gain, distinct channel signals, LFE-only, role permutations, silence, short files, impulses and near/over-full-scale float input. Every channel gets exactly one gain multiplication; no mix/downmix or LFE exclusion. Invalid/non-finite/out-of-range plans fail without mutation. |
| Formats/rates | WAV/AIFF/FLAC at all four qualified rates, supported PCM/float depths, mono/stereo/5.0/5.1. Device-rate mismatch fails with no resampling/fallback. Broader Source layouts keep existing behavior and unavailable measurement status; unsupported Current combinations fail explicitly. |
| Streaming/timing | Vary worker/callback block sizes, partial final block, FIFO wrap, underrun, rapid seek generations, source frame zero/final frame, Stop/Play/Restart, seek while stopped/playing and seek-to-EOF. Captured intervals match offline slices with no lost/duplicated samples; both domains submit identical frame/time spans for identical transport schedules. |
| Identity and UI | Wrong source hash, role order, range, domain, plan or revision cannot play. Cancelled/superseded preparation cannot publish. Gain changes invalidate Current measurements/evaluations; Source survives. Monitor/routes/profile/video/sections do not invalidate audio measurements. Test current/unavailable/stale UI states, measured failures, unapplied predictions, and M3 already-passing no-op behavior. |
| Source integrity/lifecycle | Same-size/timestamp mutation, mutation during snapshot, changed payload/metadata, missing/truncated/invalid files, cancellation during copy/hash/validation/buffering, repeated prepare/close/reopen, replacement/clear, closing during work, device open failure/loss/rate change and absent device. No unexpected open/resume, stale callback or leaked owned temporary file. |
| Output faults | Capture finite values above 1 unchanged before monitoring. With sufficient monitor attenuation, output follows the existing ramp/routing. At an actual post-monitor excess, the entire offending callback is silent and transport stops with an error; no clamp, altered plan, hidden attenuation or automatic retry. Unrouted channels stay silent. |
| Measurement bridge | Feed full captured streams to the qualified kernel and compare actual I/M/S/LRA/SP/TP and statuses with ordinary Source/Current offline analysis. Test gate-crossing and short/silent cases; do not infer loudness solely by adding trim. M2 receives only real, current offline measurements. |

Same-platform processing parity is exact, not a loose loudness-only comparison.
Across platforms, zero gain stays bit-exact; propose at most one float32 ULP for
finite nonzero processed samples, with exact frame/channel/status identity.
An ULP is one representable float step. Report maximum errors rather than rounding
them away. Retain M1 cross-platform metric budgets (0.01 LU/dB, 0.1 LU LRA) and all
official per-case limits; none is added to M2's profile acceptance allowances.

Test actual three-hour stereo 48 kHz input through preparation, random seeks and
final-second inspection on both Macs; use streamed multichannel/high-rate fixtures
for the remaining matrix. Proposed M6-owned additional buffer/workspace budget is
16 MiB, independent of programme duration (excluding existing video/decoder memory
and the source-sized temporary disk copy). UI action/refresh p95 should remain
<=100 ms. Cooperative preparation cancellation should complete <=2 s on responsive
local storage; no hard OS-read timeout is claimed. Record preparation time and
temporary disk usage separately from steady playback and video attachment time.

## Qualification sequence

1. After approval, implement/stabilize in the established logic working copy using
   a fresh `build/logic-m6` directory, external pinned JUCE 7.0.12 and minimum macOS
   10.13. Add a named audition test target; explicitly build that target plus
   Standalone, Measurement, Scaffold, Compliance, Proposal, Diagnostics, Sections,
   VideoTests and VideoQualification. Never default/all, VST3 or install targets.
2. Run M6 capture/identity/lifecycle tests, M1–M5/scaffold/video tests and runner
   checks on native arm64 Sequoia. Since shared measurement snapshot/gain code is
   extracted, rerun both pinned official suites and compare old/new measurement
   behavior; do not rely only on unchanged measurement-kernel source.
3. Run H.264 and native ProRes regressions, including movie attachment/seek/stream
   changes/clear against nonzero selected audio. The captured selected audio must
   remain unchanged; observe actual video-reader outputs and zero embedded-audio
   decoder outputs. Preserve the existing codec/metadata rejection matrix.
4. Freeze source manifest, build configuration, executable hashes, fixture hashes
   and evidence. Run that same universal candidate natively on actual Intel High
   Sierra in a new isolated folder after inspecting its existing runtime. An Intel
   slice on Sequoia is not a substitute. Any fix creates a new candidate and reruns
   affected checks on both platforms; historical failures remain recorded.
5. Complete silent product UI checks on both Macs with devices closed, routes
   muted and playback stopped, using injected devices only in automated tests.
   Coordinate physical-device runs separately before opening any device, including
   nominally silent timing runs. Agree device, routing, rate, monitor level and
   duration before short audible Source/Current comparisons; no auto-discovery
   action changes system defaults or the Pro40 configuration.
6. For changed audition/clock integration, rerun the existing short ProRes
   profile/rate matrix in Source and nonzero-gain Current modes, plus H.264 CFR/VFR
   cases, signed offsets, seeks, stop/resume/restart and EOF. Rerun the six existing
   30-minute profile/platform endurance cases with Current nonzero gain. Source
   controls and exact captured timing-span comparisons supplement these; old
   binaries' endurance evidence is not relabelled as a new candidate's run.
   Preserve document 25's zero interior actual-paint skips, zero blanks/bound
   failures/deadline failures/underruns, queue <=8, peak RSS <=768 MiB and both
   non-reusable-resident/physical-footprint growth <=32 MiB criteria. No gate is
   relaxed because processing is inaudible or synchronization appears maintained.
7. Qualify M6 only after both platforms and all required UI/device runs pass.
   If access, device approval or a required test is unavailable, record that exact
   gap; keep M6 unqualified. Internal clock/paint evidence and listening observations
   remain distinct from physical optical/acoustic synchronization measurement.

## User observation and future UI clarification

On 22 October the user reported successful manual tests on both systems with one
or two observed dropped video frames while synchronization appeared maintained.
This is informal user feedback with no new instrumented timing result. Investigation
belongs to the Video Implementation chat. Preserve formal video regression coverage
and existing acceptance gates during M6; do not reclassify the historical evidence
or treat this observation as proof of synchronization accuracy.

For interleaved audio, replace **File 1, File 2, ...** with **Channel 1, Channel 2,
...** in a future UI clarification. These rows are channels within one file.
Future separate-file import should show filenames alongside channel roles.
Separate-file import itself remains deferred. This planning task records the wording
request only; it does not change the UI or expand M6 into import work.

## Approval clarifications and exclusions

Approved clarifications (23 October): monitor listening volume remains adjustable
and smoothly ramped during playback. Only processing gain is stopped-only; monitor
changes never invalidate processing identity, measurements or evaluations.
After monitor overload, reset FIFO and volume-ramp state for explicit Play (resume
held position) or Restart (zero). Silence the whole offending callback and never
advance presentation for the suppressed span. Lowering monitor volume alone must
not resume playback; stale buffers/ramp must not cause a repeat fault.

Verified-source preparation must check available temporary space, report progress,
handle write/flush failures, clean owned partial copies and leave the previous valid
session stopped. Use deterministic injected storage failures in tests; do not fill
a real volume. Retain the verified copy across seeks and stopped domain/gain changes.

Keep processing, FIFO, identity and routing generic through the existing 12-channel
capacity and ordered height roles. Only the current measurement-support gate may
restrict Current State to M1 layouts. See
[the staged multichannel follow-up](29_MultichannelPcmMilestone.md).

Verified audition means a verified
source/plan/sample path, not an automatic measurement or compliance guarantee.

Section trims, dynamics, automatic correction, seamless A/B switching, persistence,
public export, separate-file imports, extra video features, embedded soundtrack and
VST3 parity are excluded. Sections remain in-memory annotations. Preserve all nine
existing checkpoint tags, both worktree locations, earlier tested packages and
installations. The separately authorized M6 publication and integration are recorded
in the qualification handover.

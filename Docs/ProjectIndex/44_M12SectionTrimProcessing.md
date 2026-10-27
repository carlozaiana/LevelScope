# M12 proposal — explicit section-trim processing and measured verification

Status: **approved for implementation**, 27 October 2026.

M12 approved Sequoia phase completed on candidate-03, 27 October 2026; see
[46](46_M12QualificationAndHandover.md). Full qualification remains deferred.
Required development-checkpoint label:
**“Sequoia tested; High Sierra and video qualification pending.”**

Base: published M11 `cfb1033aa7fb65a2655db276bee696161196bdeb`, annotated tag
`m11-aligned-stem-composition-2026-10-26`. Direct GitHub inspection confirms that
remote main identified this commit when proposed. Keep main at exactly this M11
checkpoint throughout M12/M13. Continue in the established
`/Users/fuzz/Developer/LevelScope` working copy on `codex/logic-module`. Preserve unrelated work and every earlier
checkpoint/package; do not integrate M12/M13 into main before combined qualification.

Read with [the proposal-engine specification](../ProposalEngineSpec_v2.md),
[M5 sections](24_M5ProgrammeSections.md), [M11 handover](43_M11QualificationAndHandover.md)
and [the measurement contract](14_OfflineMeasurementContract.md).

## Recommendation and dependency order

Make the next intervention after global trim available through the complete
measured workflow: **user-selected section gains, explicit Apply, whole-programme
remeasurement, verified audition and verified PCM export**.

The proposal specification places section/scene trims after global trim and before
dynamic processing. M5 supplies exact programme sections and descriptive analysis;
M3 supplies a measured global-trim proposal; M6/M8/M10 supply audition, sessions and
export. At the M11 baseline, code executes only `static-gain-v1`, and M5 summaries are
not independent section loudness/peak measurements or a causal edit simulator.
The existing plugin DSP does not establish qualified standalone processing parity.

M12 should therefore establish the shared section-processing and verification path
before asking an automatic method to choose section gains. It is useful for an
operator correcting a deliberately identified passage and supplies a concrete plan
representation for later proposals. It does **not** implement automatic section
levelling, promise a solution to every loudness/peak conflict, or normalize every
section to the programme target. A scene remains a user annotation, not a detected
semantic category.

Following M12, propose bounded automatic section-trim selection with explicit
alteration budgets and measured candidate assessment. Subsequent whole-programme
dynamics and residual peak control need separate proposals. Dialogue-relative DSP
follows that broader workflow, using explicitly defined analysis/acceptance methods.

## Approved scope and defaults

- New section processing supports whole-programme **stereo and conventional 5.1**,
  at the existing 44.1/48/96/192 kHz rates, up to three hours. Inputs may be qualified
  interleaved files, M9 mono groups or any of the three M11 compositions.
- Create an editable gain draft from the current M5 section boundaries. The user
  chooses each gain; default 0 dB, range **-12 to +12 dB**, stored at full
  double precision. This range is a product boundary, not an EBU requirement.
- One section coefficient acts on **every programme channel, including LFE**.
  For compositions it follows the confirmed M&E/dialogue sum. It changes the whole
  mix, not the dialogue-to-M&E balance. No component gain, mute, solo or DSP.
- Use a versioned, nominal **100 ms total transition**, centred at each internal
  change of gain. Short adjacent spans shorten the transition deterministically;
  impossible nonzero transitions are rejected rather than becoming hard steps.
- Preserve the existing global trim as a separate final processing stage with its
  existing range and numerical operation. Monitor volume/routing remain separate.
- Show draft and applied curves distinctly. **Apply section trims** is explicit
  and stopped-only. **Re-measure Current State** remains a separate explicit action.
  A structurally valid applied plan may be auditioned as **UNVERIFIED**, as in M6.
- Source stays original. Nothing here changes sample rate, duration, frame zero,
  channel order, alignment, programme coordinates or video offset semantics.
- Existing mono/5.0/7.1 static-gain workflows remain covered by regressions. The new
  section processor rejects unsupported layouts; it never falls back, drops
  channels or downmixes. Shared interfaces retain twelve-channel capacity.

## Keep annotations and audio processing distinct

M5 remains an annotation document. Merely splitting, moving, merging, renaming or
selecting or importing an annotation must not change audio or invalidate unchanged measurements.

A draft records the exact input identity, Current revision and annotation geometry
from which it was made. Apply installs an immutable processing envelope with its
own frame intervals and gain values. Later annotation edits do not move an applied
gain transition. The UI shows when annotations and applied processing boundaries
differ; updating processing requires an explicit new draft and Apply.

Geometry changes make a dependent draft stale. Refreshing the draft requires review;
do not silently transfer gains to newly split, merged or moved intervals. Cosmetic
names stay outside audible identity. Applied interval/gain data and its original
annotation provenance remain inspectable independently of the current annotations.

Drafts are not saved as applied settings. Save/New/Open/Close with unapplied edits
must offer an explicit resolution; never silently apply them or claim they were
saved. A bounded first UI can require Apply or Discard before Save, and Discard or
Cancel before leaving the draft. A history browser, arbitrary curve drawing and
general undo/redo are outside this milestone.

## Follow-up import and mono milestones

Use the existing manual section editor in M12. Follow with marker import, then
bounded EDL import: explicit programme timecode origin/frame-rate mapping, exact
audio-frame boundaries, preview and confirmation. Imported edits are proposed
annotations only, never automatic processing changes. Neither importer is M12 work.

Preserve existing mono measurement/global trim. Keep channel-driven processing
capable of mono, but reject mono section processing in this release. Complete mono
section-processing and verified PCM export support is a near-term separately
qualified extension, not part of M12 stereo/5.1 claims.

## Shared processing and numerical contract

Introduce a typed, versioned whole-programme gain-plan model/evaluator in Core,
with the standalone worker adapter. Do not create an arbitrary module graph or
reuse mutable plugin parameter state. Keep the existing static-gain implementation
and plugin processing behavior unchanged.

Processing order is:

`original verified programme/composition → section gain → qualified global trim`

The section evaluator takes an absolute audio-source frame index. Its result must
be independent of block size, callback boundaries, seek history or preceding reads.
It adds no lookahead, latency, filter state, resampling or full-length render cache.
Measurement, audition worker and export use this exact same processing entry point.
Audition's FLAC decoding preroll must not shift the processing frame coordinate.

Approved `section-gain-v1` convention:

1. Validate finite gains, supported input, contiguous nonempty half-open intervals
   covering `[0, frames)`, bounded count and exact integer frame positions.
2. Canonicalize adjacent intervals with identical gains before constructing ramps.
   IDs, labels and redundant equal-gain boundaries do not alter the audio digest.
3. At a changing boundary `b`, set radius `r` to the minimum of rounded 50 ms at the
   source rate and half of each adjacent canonical span, rounded down. The ramp
   covers `[b-r,b+r)`. This prevents adjacent ramps from overlapping.
4. If `r == 0`, reject the differing-gain transition. Tiny annotations with equal
   gains remain valid. Do not silently introduce a discontinuity.
5. Across the `2r` samples, interpolate gain linearly in dB with
   `t = (frame-(b-r))/(2r-1)`. Endpoints equal the adjoining plateau gains. Outside
   ramps use the relevant plateau. There is no implicit programme-start/end fade.
6. Convert dB to a double coefficient once per frame and multiply each float32
   channel sample, rounding to float32. Then execute the unchanged global-gain
   operation with its existing separate rounding boundary. Do not fuse stages.

Zero section gain bypasses its multiplication. An entirely zero envelope uses the
legacy static-gain plan/path, preserving its samples and identity. A draft whose
canonical plan equals the currently applied plan is an audio no-op, including
equal-gain repartitioning; no new revision or remeasurement just for labels.
Returning from a different active plan to a historical plan is an audio change
and must follow the normal revision and invalidation rules.

Finite above-full-scale output remains measurable without attenuation or clipping.
Non-finite inputs/results fail explicitly. Retain monitor-overload stop/restart
behavior and M10's full-scale/export gates. Dither occurs once, at final PCM24
quantization only; section processing and Current measurement are not dithered.

This convention is an engineering choice to qualify, not a standard-prescribed
fade or a guarantee that every transition is perceptually inaudible.

## Measurement, proposals and identities

Source/Current must still describe the same content, ordered roles, composition
and programme range. Current adds the entire canonical processing plan: operation
versions/order, exact gain bits, interval frames and ramp policy. Bind worker jobs,
snapshots, diagnostics, prepared audition and export approvals to that identity.

| Action | Required effect |
| --- | --- |
| Draft edit, selection or annotation-only edit | No audio change; retain valid measurements/evaluation |
| Apply different section processing or global gain | Preserve Source; invalidate Current measurements/evaluation, proposals, prepared Current playback and export approval |
| Cancel/fail remeasurement after Apply | Retain the applied plan; Current stays clearly unverified; no rollback, retry or hidden correction |
| Complete failing measurement | Show the actual measured failure; never expose an older pass as current |
| Compatible profile/applicability change | Reevaluate valid measurements; no audio remeasurement |
| Monitor routing/volume or video attachment/offset | Preserve audio measurement/processing identity; maintain existing report-context rules |
| Member bytes, roles, mapping or alignment change | Invalidate affected audio; do not silently rebind the applied envelope to different input |
| Identical verified relink | Preserve audible plan identity, subject to existing cache compatibility |

Global-trim proposal support is included for a **freshly measured section-processed
Current State**. Keep the section plan fixed and propose only a new final global
gain, with the existing deterministic decision policy and narrow explanation.
Version the method/adapter and bind the exact fixed upstream plan into the draft;
do not merely bypass M3's current `static-gain-v1` check. Applying that proposal
must preserve every section setting. An already passing Current makes no change.
If no global-trim solution exists, explain that limit; M12 does not search for new
section gains or add a limiter. Prediction never satisfies a measured rule.

M4/M5 diagnostics are rebuilt from actual Current measurements. Keep programme
gates, short-programme warnings, missing values and descriptive LRA semantics.
Do not label a section median or contained-window share as section integrated
loudness, independent peak safety or the predicted effect of a trim.

The versioned EBU profile/digest remains unchanged. Whole-programme evaluation
retains the -23 LUFS target, explicit +/-0.2 LU QC allowance and strict -1 dBTP
ceiling. The standard's separate true-peak measurement-tolerance provision does
not enlarge that ceiling. No universal LRA ceiling or per-section normalization
rule is added. [EBU R128-2023](https://tech.ebu.ch/docs/r/r128.pdf) calls for whole-
programme measurement; [BS.1770-5](https://www.itu.int/rec/R-REC-BS.1770-5-202311-I/en)
remains the existing measurement-method reference. Processing defaults are not
new delivery requirements or intelligibility claims.

## Sessions, reports and future extension points

Recommend session schema 4 with a typed Current plan and readers for schemas 1–3.
Old documents stay untouched until explicit Save; warn that the new format needs
an M12 reader and offer Save As. Preserve exact frame/gain/ramp values, independent
annotation data, composition confirmation, profile/applicability and rational video
alignment. Reopen stopped, device-closed and muted. Reject unknown processing
versions rather than silently playing unity. Corrupt/incompatible caches remain
discardable without losing a valid core document; reconstruct evaluations from
validated measurements and preserve original provenance.

Save a coherent revision; later edits remain dirty. Runtime job identities must
continue to reject late completions after New/Open, even for the same document ID.
Reports should identify the exact applied plan and transition intervals, using a
new report schema version. Retain separate Source/Current/artifact domains and
measurement-only reports. Public audio export still requires passing measured
Current and passing finalized, reopened artifact at zero extra gain. No automatic
attenuation/retry when an encoded artifact fails. No new export format or layout.

Preserve separately accessible original M&E/dialogue components. The new envelope
is whole-programme processing, not a permanent substitute for a component stage.
Future component DSP must precede composition and enter Current identity. Original
Source remains the confirmed composition. Analysis-only references remain
non-audible with independently invalidated dependent analysis. Dialogue loudness,
speech-gated mix measurements, dialogue/background balance and intelligibility
remain separate future methods requiring explicit targets and verification of both
the intended improvement and the resulting complete programme.

## Acceptance tests

1. **Independent sample oracle:** stereo/5.1, all four rates, interleaved/grouped/
   composed inputs; zero/positive/negative global and section gains; shuffled
   confirmed roles; LFE and explicitly muted monitoring channels. Compare every
   sample in short offline, pre-monitor and pre-quantization captures. Same-build
   consumers must agree exactly. Cross-platform processing uses the existing
   at-most-one-float32-ULP bound; never add this allowance to compliance thresholds.
2. **Envelope boundaries:** first/last sample, exact half-open endpoints, adjacent
   and short spans, shortened ramps, minimum representable ramp, unequal/equal
   neighbours, 1,000/10,000 spans, full-precision gains and invalid/non-finite data.
   Independent oracle verifies the formula, clipping absence and constant/zero
   behavior. A no-op keeps the audio revision and valid results unchanged.
3. **Random access and lifecycle:** compare full sequential processing with seeks
   before/inside/after every kind of transition, different chunks, restart/resume,
   EOF and FLAC preroll. Preserve exact duration/positions. Test cancellation,
   close/source replacement, missing members, stale jobs, changed roles, disk-full
   preparation and write failures. Preserve prior valid sessions/owned-copy cleanup.
4. **Measured workflow:** apply a locally selected correction, observe unverified
   Current, remeasure and evaluate actual output; draft/apply M3 global trim while
   retaining the envelope; verify already-pass no-op, no solution and prediction/
   measurement disagreement. Cancelled or failing verification never revives an
   old pass. Monitor-muted LFE still determines peak constraints.
5. **Annotation separation/UI:** navigation never plays or edits audio. Pending
   drafts and applied curves are distinct; annotation changes do not alter the
   active envelope; stale drafts cannot apply. Test all stopped-only processing
   changes, live monitor volume, overload recovery, explicit discard and save/open
   with unapplied edits. Exercise the silent checklist on Sequoia now; actual High Sierra is deferred until after M13.
6. **Session/report/export:** exact plan round-trip (cross-Mac exchange deferred until after M13); older
   schema migration and old-reader rejection; damaged cache recovery; original
   confirmation retention; stale profile/report/export rejection. Independent
   decoder/quantizer oracle verifies PCM24, one final TPDF stage, unchanged rate/
   frames/mapping, strict Current/artifact gates and measurement-only reports.
7. **Compatibility and timing:** focused M1–M11 cases for legacy static gain,
   interleaved mono/stereo/5.1/7.1, mono groups and all M11 presets. New envelopes
   explicitly reject unqualified layouts. Pinned EBU/ITU checks run because the
   shared measurement dispatch changes. Audio-only injected seeking, frame-position, duration, EOF and continuity
   tests remain required. All tests that load, decode or play video—including
   injected H.264/ProRes tests—are deferred until after M13. Preserve video code
   and rational alignment data; saved alignment is tested without decoding media.
8. **One consolidated Sequoia long workflow:** three-hour 48 kHz 5.1 composition,
   1,000 gain intervals, Source/Current measurements, M3 final trim and remeasurement,
   first/middle/final/boundary captures, diagnostics, save/reopen and verified
   PCM24 RF64 export. Reuse verified immutable long inputs where available. No
   separate repetitions of equivalent older three-hour workflows.

## Impact-selected qualification and resource plan

Before implementation, record an explicit test selection covering direct edits and
callers/build dependencies. Stabilize on Sequoia, freeze a new universal candidate
and complete the required Sequoia audio-only phase. Build x86_64 locally using the
established 10.13 deployment target/dependencies; compilation is not actual High
Sierra execution. Do not access Golem or run any video-loading/decoding/playback
test during routine M12 qualification. After M13, qualify the combined candidate
on actual High Sierra and run selected deferred video checks before requesting
main integration. Earlier M11 results never qualify the new M12/M13 features.

Expected affected areas are the typed Current plan, shared processing dispatch,
measurement identity/worker, audition worker frame positions, M3 adapter, sessions,
reports/export and section-edit UI. Reuse unchanged measurement-kernel, decoder,
real-time callback and broad video evidence only with exact original identities
and documented dependency equivalence. Do not silently invoke recursive M1–M11
runners. The M12 runner must enumerate required gates and fail incomplete manifests,
missing required Sequoia tests or UI results. Keep fresh, justified reused and
explicitly deferred evidence separate; deferral is not a pass. Rerun affected checks after a
fix; expand coverage only for observed failures or uncertain dependencies.

Approved resource budgets for the Sequoia phase:

- New envelope evaluator/plan storage: <=4 MiB at 10,000 intervals, bounded block
  work, no full-programme float render/cache. Full workflow additional working RSS
  <=256 MiB above existing measured-result storage.
- UI p95 <=100 ms; cancellation acknowledgement <=1 s; save and core/cache parsing
  <=5 s; session <=128 MiB; JSON+HTML <=16 MiB.
- Consolidated long workflow <=180 minutes on Sequoia; planning estimate 60–90 min
  using M11's 56.54-minute observation as a baseline plus an additional trim cycle.
  Short/reference regressions approximately 20–40 min; silent UI approximately
  20–40 min. Estimates are not evidence. Replan combined High Sierra work after M13.
- Reused long inputs do not need another permanent copy. Allow approximately
  65 GB free working headroom for verified private copies, output, tests and reserve;
  measure actual peak usage and clean only owned temporaries. Follow the compact
  M11 packaging approach, not repeated copies of all historical M9/M10 media.

The prior `build/logic-m*` build trees have been cleaned; their configurations are
archived under `build/sable-cleanup-2026-10-26`. Use a fresh M12 build directory and
only named standalone/test targets, never the default all/VST3/install target.
Frozen M9/M10 packages now resolve through their original symlinks to
`/Volumes/Tucan/LevelScope-qualification`; require the archive volume when those
assets are selected. Do not overwrite packages, regenerate missing assets in place,
or treat unavailable archived evidence as a pass. Do not access Golem in this phase;
its preserved ProRes path and all tested packages remain unchanged.

Required M12 playback checks use injected audio-only devices and pre-monitor
captures; silent UI checks leave physical devices closed. No physical device is
automatically opened. Any optional audible check requires separate coordination.
All actual High Sierra and all video-loading/decoding/playback checks remain
deferred until after M13. If changes to timing, buffering, presentation coordinates
or video handling become necessary, or a concrete failure warrants earlier checks,
explain that reason to the user before revising this deferral. Passing an existing
frame index to the new stateless gain evaluator does not itself change presentation
coordinates; regression-test that audio-only contract explicitly.

Retain occasional picture skips as an accepted development limitation for the
later video phase; drift, sustained freezes, audio interruptions and identity/
lifecycle/memory failures remain blockers. Full physical multichannel listening
remains separate.

## Exclusions and approval boundary

No automatic section-gain selection or optimizer, section-by-section normalization,
LRA target/steering, semantic scene/speech detection, compressor, limiter, per-stem
processing, freezes, general automation editor or plan-history browser. No dialogue
analysis/intelligibility score, 7.1 section-processing expansion, 7.1.4/Atmos,
resampling, new containers/formats, additional video features or VST3 parity.

The user approved the manual-first boundary, gain range, transition convention and
annotation/audio separation, with the revised Sequoia-only phase above. Stop once
that phase is complete for readiness review, before checkpoint publication or M13.
Report exact source/binary identities, fresh/reused evidence, required checks and
deferred High Sierra/video checks. No main integration is authorized.

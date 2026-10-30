# M17 proposal — explicit job LRA requirements and reusable presets

30 October 2026. **Approved for implementation and impact-selected Sequoia qualification, 30 October 2026.** Publication and main integration require later review. Candidate-04 completes the approved Sequoia phase; [62 — qualification/handover](62_M17QualificationAndHandover.md) records exact identity, results and pending combined qualification.

Continue from published M16 candidate-03, commit
`73a7b25e71f87299fa2d93c1b1106c118c55beec`, annotated tag
[`m16-proposal-preferences-sequoia-2026-10-30`](https://github.com/carlozaiana/LevelScope/tree/m16-proposal-preferences-sequoia-2026-10-30).
The documentation-only publication follow-up does not change that implementation.
Use the established logic working copy and `codex/logic-module`. Main stays at
qualified M13 `a9be070ab9eb5cb4f591ff9d991d8999c3a016d0`.

Retain **Sequoia tested; High Sierra, video and three-hour checks pending** for
development checkpoints. Combined qualification follows M18, no later than M19,
before separately approved main integration. Preserve every checkpoint, package,
original media file, installation and the video worktree.

## Outcome and dependencies

Let a user select the existing EBU programme profile, explicitly add a job-specific
maximum LRA, and save/load a named reusable requirements preset. Evaluate Source,
Current and actual exported audio against the selected rules. Make an unmet or
unavailable additional requirement visible throughout the workflow, without
inventing an automatic LRA correction strategy.

Dependencies already exist: M2 measured rules and profile provenance; M4 final
LRA and interpretability; M8–M11 coherent sessions and export approvals; M12 manual
processing and independent remeasurement; M13 proposal identity; M16 preferences,
complete inspection and stale-result handling. Read ProjectIndex 14, 21, 39, 47,
56, 58 and 59 and ProposalEngineSpec_v2 before implementation.

At the published M16 baseline, optional LRA evaluation existed in `TargetProfile`/`VerificationEngine`,
including unavailable and unstable handling. The authoritative EBU catalogue entry
has no LRA ceiling. The old `RuleBasedMethod` rejected profiles containing an LRA
rule before its no-change branch; `SectionTrimMethod` depended on that decision.
Session schema 6 persisted a catalogue profile reference. These shared
callers require explicit integration; merely exposing an existing numeric field is
insufficient. No new audio analysis or arithmetic is required by this proposal.

## Defined rules and interpretation

Keep the authoritative `ebu-r128-2023-programme-linear-qc` definition/version and
digest unchanged:
`ef6f5cc071a8c9d5a347da0dd61fcc46e9f24da0aed787bc9f3ea48cb3db7e66`.
The intended −23 LUFS target, ±0.2 LU QC allowance, strict −1 dBTP ceiling and
whole-programme linear-audio applicability remain unchanged. EBU R128 does not
supply a universal maximum LRA; its short-programme caution must remain visible.
[EBU R128 v5 (2023)](https://tech.ebu.ch/docs/r/r128.pdf)

An added rule is **measured whole-programme LRA ≤ the explicitly selected job
ceiling**, in LU. Use the existing qualified EBU Tech 3342 method and method/build
compatibility checks. LRA is a gated statistical loudness-range measure, distinct
from peak-to-loudness range and dialogue-gated metrics; no rolling or partial value
can satisfy the final rule.
[EBU Tech 3342 v4 (2023)](https://tech.ebu.ch/files/live/sites/tech/files/shared/tech/tech3342.pdf)

Preserve the current stability indicator: a result from a programme under 60 seconds
cannot satisfy a required LRA rule. The number may remain available for descriptive
inspection with its limitation. A longer programme still needs a valid, available,
finite metric. Silence, excluded windows and unavailable results never become zero.
No selected LRA rule means unavailable LRA cannot block the existing programme rules.

Evaluate full-precision measured values. Use an inclusive ceiling without epsilon,
rounding or measurement-engine test tolerance added to the requirement. Display
rounding cannot change its verdict. Result wording identifies Source, Current or
measured export and states only that the selected programme/job rules pass. It does
not certify delivery metadata, a broadcaster contract or an unmeasured future file.

## User-facing scope and defaults

Use an inline **Job requires a maximum LRA** control beside the target profile;
default unchecked, with **No LRA acceptance limit specified** displayed. Selecting
it requires an explicit value; do not prefill a supposedly universal 7, 10 or 18 LU
limit or silently infer one from genre, the cinema example or proposal preferences.

Approved input range: **0.0–100.0 LU in exact 0.1 LU steps**, stored
as integer tenths. This is an input/engineering bound, not a broadcast recommendation.
Zero is a real strict ceiling, distinct from unset. Reject malformed, negative,
non-finite, unsupported-precision or out-of-range values without clamping. The
measurement being compared retains its full precision.

Show the base profile and added rule separately. The added rule is labelled
**user-specified job requirement**, with a name and optional requirement-document
reference/edition. User-entered provenance is not independently verified and cannot
be labelled an official EBU rule. No arbitrary rule editor or modification of the
base IL/TP allowances is included.

Provide **Save requirements preset…** and **Load requirements preset…**, using one
small versioned JSON file per named preset (`.levelscope-target`). Include stable
preset ID/revision, definition digest, base profile reference, method/applicability,
the optional ceiling and provenance. An edit creates a new revision; loading copies
an immutable definition into the session. Saving a preset must not silently change
the active definition. No installation-wide preset database or hidden settings are
needed in this milestone. A user can organize those files in an ordinary folder.

The preset contains requirements, not gain settings, media paths, measured results,
verdicts or M16 processing preferences. A session stores the complete selected
definition, so renaming/deleting/modifying an external preset cannot change it.
Unknown IDs, unsupported methods, digest mismatches and malformed presets fail
transactionally without changing the active session. No preset file can request
network access, execute code or acquire authority by setting an `authoritative` flag.

## Identity, session, evaluation and export integration

Introduce a small versioned job-requirements definition around the immutable base
profile. Keep base-profile identity, added-rule identity and the combined requirements
digest explicit. A definition includes its method, applicability and provenance;
a new definition invalidates approvals, even when its numeric thresholds are equal.
Identical confirmation is a genuine no-op. Keep M16's policy digest separate.

Edit/load requirements while stopped with no conflicting active job. Use the existing
coherent session revision rules. A change marks the session unsaved, preserves
unchanged Source/Current audio identity and valid measurements, and immediately
reevaluates them. Incompatible measurements remain unavailable until remeasurement.
Previously generated proposals, private assessments and export-publication approvals
bound to older requirements become stale. Late completions cannot update a newly
opened session. Monitor settings and video alignment do not alter audio identity.

The session schema should advance from 6 to 7. Schemas 1–6 open with their exact
media, applied processing, annotations and M16 preferences, with no added LRA rule.
Do not overwrite the original on open; explain the M17-reader requirement on explicit
Save and offer Save As. Keep caches replaceable. An independently identifiable
invalid requirements record may recover otherwise valid settings/media with
requirements unavailable and all dependent approvals disabled until explicit
replacement; do not silently fall back to plain EBU. Invalid core structure still
fails safely. Reconstruct evaluations from validated measurements, never saved
verdicts, and preserve original measurement provenance.

For export, both Current and finalized/reopened artifact must pass **all** selected
rules. Source need not pass. Keep exact artifact/audio processing identity separate
from approval/report metadata identity: a requirements edit can reuse a valid,
unchanged rendered artifact and its compatible measurement, but must reevaluate it
and refresh its approval/report. It cannot publish under an old approval. Preserve
separate Render and verify / Publish export, final-only dither and exact artifact
hashes. A measurement-only report remains available for failed/unavailable results
without claiming an artifact measurement, hash or verdict.

JSON/HTML reports show base IL/TP and additional LRA results and provenance separately,
the combined requirements identity, full-precision values, stability/unavailability
reasons, and the evaluated domain. Version the report extension. Escape all names
and provenance; a preset is data, including in HTML.

## Approved automatic-proposal boundary

**Approved bounded M17 behavior:** first validate matching measured inputs and
the applied plan, then return a genuine no-change result when **all selected rules**
already pass. Preserve the exact gain/envelope, measurement/evaluation objects,
audio revision and job count. Do not move a satisfactory result towards a preferred
number. An IL/TP pass cannot hide failed, missing or unstable required LRA.

For this milestone, keep changing automatic proposals available for the existing
programme-only profile. When an additional LRA rule is selected and the complete
Current result does not pass, explain **Automatic correction for this selected set
of requirements is not supported yet**. This applies even when LRA alone passes
and IL/TP needs correction: the existing scalar prediction cannot promise preserved
gating/LRA. Show the actual failing/indeterminate rules, not a claim that no solution
exists. Manual M12 processing, explicit Apply and independent remeasurement remain
available; completed results are evaluated against the complete selected rule set.

This deliberately keeps M17 to requirements/presets and trustworthy decisions. It
does not strip LRA internally to enable an apparently compliant proposal. Supporting
changing IL/TP proposals under an additional LRA constraint would need exact private
candidate assessment on the global-only path as well as the section path; that is
a separately reviewable extension. LRA-directed region selection/correction remains
a later bounded milestone. No new strategy, pass budget, retry, attenuation or DSP
is introduced here. Selecting a less restrictive requirement is an explicit user
decision recorded with its identity, never an automatic workaround.

## Acceptance tests

1. **Rules and numerical boundaries:** unchanged base digest/IL/TP decisions; no LRA
   rule versus explicit zero; exact ceiling and `nextafter` either side; misleading
   display rounding; silence, missing/non-finite/negative metrics and every unavailable
   status. Stable/unstable results, under/exactly/over 60-second programmes, and
   informative-but-unacceptable short LRA. Required unavailable LRA never passes;
   optional unavailable LRA never introduces a failure.
2. **Presets and validation:** named round trips, stable revisions/digests, duplicate
   names distinguished by ID, unsupported method/edition/base digest, tampering,
   malformed/oversized input, numeric bounds, Unicode and HTML escaping. Cancelled
   dialogs or write failures preserve the active definition and prior file. No
   silent authority elevation, media access or external-file dependency after load.
3. **Identity and reuse:** changing only the ceiling reevaluates cached measurements
   without an audio worker; processing/audio/mapping changes retain their normal
   invalidation. Stale proposal Apply/export approval/late completion is rejected.
   Identical selection is a no-op. M16 preference changes remain independent of
   job-rule evaluation. Test video-offset metadata identity without loading video.
4. **Proposal behavior:** all selected rules pass → zero audio/revision/job changes;
   IL/TP pass plus LRA fail/indeterminate → never already compliant. Changing
   proposals with selected LRA are clearly unavailable, with manual workflow intact.
   Plain-EBU global-first/section, policy limits/full inspection, Apply and measured
   failure/cancellation behavior remain unchanged. Predictions never satisfy LRA.
5. **Sessions and reports:** schema 1–6 migration defaults unset, exact processing/
   markers retained; schema 7 self-contained round trip; old M16 reader rejects it
   safely; incompatible/damaged cache recovery; invalid requirements recovery without
   fabricated pass; edits during save remain unsaved. All report domains, unknown
   metrics, external provenance and content/hash checks are exercised.
6. **Export/manual verification:** real short mono/stereo/5.1 fixtures exercise
   manual processing and independent Current/artifact measurement. A Current failure
   blocks rendering/publication as specified; an artifact-only LRA failure blocks
   publication despite a passing Current. Use injected finalized-result boundary
   cases alongside real decoder/worker tests to cover tiny dither-induced differences
   deterministically. No artifact hash/verdict is invented in measurement-only mode.
   Reevaluate an unchanged artifact after requirement changes without rerendering.
7. **Silent UI and resources:** inline unset/explicit requirement, provenance, preset
   save/load/error/cancel, all-pass no-op, unsupported automatic correction, manual
   verification, session recovery and reports; remain stopped/device-closed/muted.
   Enforce the budgets below and missing-required-test detection in the selected runner.

## Qualification plan and resources

Before implementation, record an explicit impact manifest. Fresh Sequoia selection:
new requirements/preset tests; measured evaluation; global/section proposal and M16
preferences tests; short session migration/recovery; short export/report and manual
processing regressions; representative interleaved, assembled and composed programmes;
silent UI and the existing universal/Intel 10.13 cross-build. Build only named
standalone/test targets. Do not invoke recursive historical runners.

Use short synthetic audio, normally 65–120 seconds for stable-LRA integration, plus
short boundary fixtures. Reuse documented M13–M16 lower-level arithmetic/decoder/
audition evidence only if unchanged, with original identities and explicit rationale.
No automatic full cinema or sustained workflow rerun. If shared measurement arithmetic,
buffering, clocks or dependencies change, reassess impact and explain the concrete
need before revising the agreed deferrals.

Proposed resource gates: preset ≤64 KiB, name ≤80 Unicode characters, provenance
≤2 KiB UTF-8; no unbounded directory scan. Cached three-rule reevaluation p95 ≤10 ms;
control feedback p95 ≤100 ms; added retained requirements/preset state ≤1 MiB;
explicit file operations off the audio callback. Preserve existing session/report
size and resource caps; typical short-session Save/reopen ≤5 seconds. Transactional
failure and cancellation preserve the previous valid document. No new audio copies
or private measurement passes solely for a requirements change.

Planning estimate, not measured qualification: 5–15 minutes for the selected frozen
automated run, 15–30 minutes for silent UI, ≤2 GiB disposable fixtures/results,
and roughly 1–3 GiB for a new frozen universal app/test package. Incremental build
time is additional and host-dependent; check space and reuse the established build
tree without deleting prior candidates. Keep the no-three-hour phase boundary.

Actual High Sierra/cross-Mac exchange, video load/decode/play and three-hour workflows
remain pending until the combined post-M18 qualification, no later than M19.
Local Intel compilation is not runtime evidence. Coordinate physical-device use
separately. Accepted occasional picture skips do not excuse audio discontinuity,
accumulating drift, sustained freezes or identity/lifecycle failures when video
testing resumes. Full physical 7.1 listening remains independently pending.

## Exclusions and completion boundary

No LRA correction method, coordinated neighbour adjustment, automatic segmentation,
selective raising, limiter/compressor, new dialogue metric, general target editor,
new codec/layout, VST3 parity or physical multichannel claim. Automatic analysis,
progressive displays and timeline editing remain the following separately bounded
proposal; do not make them prerequisites for this requirements layer.

The user approved the unset default and 0.0–100.0 LU input range, portable named
preset files/session snapshots, and the conservative automatic-proposal boundary.
M17's Sequoia phase stops for readiness review before publication or merging.

The upcoming timeline must eventually restore applied-proposal explanations along
with exact processing: reasons, relevant baseline and requirement/policy provenance
are historical information. Do not restore private candidate measurements or old
verdicts as current verification. Mark explanations whose processing was subsequently
edited. This is architectural preparation, not complete explanation persistence in
M17. Keep the timeline/progressive-analysis slice bounded so it does not unnecessarily
delay region suggestions and selective LRA correction. See the impact plan in 61.

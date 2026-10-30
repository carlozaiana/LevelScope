# M16 proposal — bounded proposal preferences and configurable limits

30 October 2026. **Approved revisions implemented in candidate-03. Sequoia tested; High Sierra, video and three-hour checks pending.** See [59](59_M16QualificationAndHandover.md) for exact identity/evidence; published as [`m16-proposal-preferences-sequoia-2026-10-30`](https://github.com/carlozaiana/LevelScope/tree/m16-proposal-preferences-sequoia-2026-10-30), commit `73a7b25e71f87299fa2d93c1b1106c118c55beec`; main integration remains deferred.

Start from published M15 candidate-05, commit
`1b56804ecfb469aef3276b2d4ac96039ce3e4cc1`, annotated tag
[`m15-marker-import-sequoia-2026-10-29`](https://github.com/carlozaiana/LevelScope/tree/m15-marker-import-sequoia-2026-10-29).
Its source manifest is
`38676a45d47ddba6eb54f7ab299f52f622f6be26554256781b8960749c609564`.
Use the established `/Users/fuzz/Developer/LevelScope` working copy on
`codex/logic-module`. Main remains qualified M13
`a9be070ab9eb5cb4f591ff9d991d8999c3a016d0` until combined qualification and subsequent
integration approval. Preserve all checkpoint objects/targets, original media,
frozen packages, qualification receipts, installations and the video working copy.

## Outcome and dependency boundary

Let the user explicitly set a small number of proposal limits, see exactly which
constraints prevented an attempt, and reproduce the same decision from the same
measured baseline and policy. These preferences bound an existing section-led
method; they are not delivery rules or evidence of perceptual quality.

Dependencies already present: M2 versioned profiles/evaluation; M3 measured
Source/Current identity and no-op/global-first decisions; M8–M11 programme/session/
export identity; M12 canonical section gains and shared offline/audition/export
processing; M13 attribution and one-envelope candidate assessment; M14 mono and
RF64 input support; M15 marker mapping/provenance. Read ProjectIndex 14, 47, 49,
52, 53 and 55 and ProposalEngineSpec_v2 section 17 before implementation.

Retain the M13 decision sequence: valid measured matching Source/Current, genuine
already-pass no-op, global-only correction first, then the existing peak-attributed
section attenuation method if needed. One proposed envelope, at most three private
measurement/attribution passes, explicit stopped-only Apply, independent Current
remeasurement and actual outcome display. A failed/cancelled verification retains
the applied settings as unverified; no rollback, retry or second strategy. Additional
objectives need their own separately validated methods/acceptance criteria; none
are inferred or added by M16.

## Proposed controls and defaults

Use a session-scoped **Proposal limits** panel with Edit, Use these limits and Cancel.
There is no general preference editor or global installation preference. Defaults
are conservative continuations of existing behavior, not values chosen for the
cinema file. All supported controls are hard proposal constraints.

| Control | Proposed default | Supported range / representation | Reason |
| --- | --- | --- | --- |
| Maximum changed processing regions | No region-count restriction | Optional explicit integer 1–10,000, within structural processing capacity | Count is a preservation preference, not an engineering budget. No ceiling of 64 and no truncation to fit a limit. |
| Maximum newly affected local-envelope duration | 20% | 0–100%, 0.1 percentage-point steps; exact integer thousandths of programme length | Keep the existing default. Wider values explicitly allow broader local processing, with preview of seconds/percentage. Zero disables changing the local envelope. |
| Maximum additional local attenuation per proposal | 3 dB | 0–6 dB, 0.1 dB steps | Preserve the current default; allow a bounded reviewable extension, still limited by the cumulative and processing bounds. Zero disables local attenuation. |
| Maximum resulting cumulative local attenuation | 12 dB | 0–12 dB, 0.1 dB steps | The existing M12 lower plateau bound is −12 dB. Making that ceiling explicit preserves current reach; users may choose a stricter value, e.g. 3 dB. It is not a recommended amount to apply. |
| Requested planning reserve | 0.5 dB | 0–1 dB, 0.1 dB steps, not above the additional cap | Preserve the existing reserve inside the correction cap. Zero requests no reserve; partial achievement is disclosed. |

Store gain preferences as integer tenths of a dB and duration as integer thousandths
of programme length. These numeric preferences have explicit stored units; do not store locale-formatted
strings or round evaluated audio to the controls' step size. Existing gains remain
full precision. Reject non-finite/out-of-range values and incompatible settings;
do not silently clamp a saved policy. A zero additional cap requires reserve zero
(or the user explicitly chooses both); show the conflict without modifying settings.
An explicit optional count limit is distinct from the always-enforced 10,000-interval structural capacity. Cumulative gain remains bounded. Duration 100% is an
explicit whole-programme local allowance, not bypass of any other gate.

Keep non-configurable: target profile and its exact rules; local upward proposal
gain zero; M12 manual plateau range ±12 dB; global range −60 to +24 dB; nominal
100 ms transition/shortening/rounding convention; candidate/pass/resource limits;
all-channel processing including LFE; explicit application and verification gates.
Changing the selected profile remains the existing separate workflow. Do not add
measurement test tolerances to programme allowances or relax true-peak ceilings.

### Additional, cumulative and whole-chain gain

Freeze the measured Current plan as baseline for each request. Additional reduction
is the positive difference between its old and proposed local plateau gains.
Cumulative local attenuation is the resulting local envelope relative to unity
section processing of the original confirmed Source, including existing manual
settings. It is not a sum of historical automatic requests, and not global trim.
No attempt is made to infer unavailable manual-versus-automatic edit history.

Example: an existing −2.5 dB plateau followed by another −2.5 dB correction is −5 dB
cumulative attenuation. It satisfies a 3 dB additional cap but fails a user-selected
3 dB cumulative cap. Repeated requests and save/reopen cannot reset that comparison.
A positive existing manual plateau can be attenuated without first being negative;
report the old/new values and delta explicitly rather than calling all deltas
attenuation relative to Source. Report local, global and total static-plus-local
resulting gains separately; a local cap is never a whole-chain attenuation promise.

If an existing manual envelope already exceeds a newly selected cumulative limit,
leave it unchanged. Already-pass/no-op and global-only paths remain available under
the measured programme rules. A new local proposal is unavailable until the user
explicitly chooses compatible limits or manually changes/reverifies that envelope.
Do not erase manual processing, quietly grandfather an over-limit local candidate,
or turn a preference violation into a measurement/compliance failure.

For a positive peak deficit D, additional cap A, cumulative ceiling C and old local
gain g, available reduction headroom is `min(A, g + C)` (also subject to valid M12
bounds). After rejecting an already incompatible baseline, choose the representable
reduction up to `min(D + reserve, headroom)` using the existing inward rounding rule.
Report requested/actual reduction, achieved reserve `min(reserve, max(0, actual−D))`,
and uncovered deficit `max(0, D−actual)`. Reserve consumes both caps; it never shifts
the profile ceiling. Candidate measurement, not this planning calculation, decides
whether the complete proposed programme passes. Partial reserve or an uncovered
planning deficit must remain visible even if the exact result differs.

## Marker density and changed audio

Frozen M13 code counts `corrections.size()`: changed annotation rows. It also
unions changed interval/transition support for the 20% check. M16 should version a
new accounting policy, **not silently reinterpret the old eight-section rule**.

Display both **changed annotation rows** (traceability) and **changed processing
regions** (the configurable count gate). A region is a connected interval of newly
changed local-envelope support relative to the frozen Current envelope, after
canonicalizing redundant equal-gain boundaries under `section-gain-v1`. Join touching
support intervals; do not bridge an unchanged audio-frame gap. Include complete
transition support from before and after, conservatively counting boundary endpoints
as required by the versioned sample-coordinate convention. Clip to [0, programme
frames); never extend programme duration for filter padding. The duration gate is
the union of that support, not a sum of row lengths or a percentage of marker count.

The implementation partitions both canonical envelopes into plateaus/ramps, merges their event coordinates, and counts spans with different plateau values or ramp definitions, including complete ramp endpoints. Touching spans join. This conservative versioned support may include an equal ramp endpoint; it is not a sample-by-sample count of unequal floating-point values. Independently sampled tests must show that every unequal processing sample is covered. Identical
canonical envelopes must have zero change. An arbitrary rename, ID change or redundant
equal-gain split must not increase measured support or region count when the actual
before/after processing is identical. Two separated modified spans remain two regions;
adjacent changes can be one region even with different plateau levels. Also report
canonical plateau/transition count, maximum delta, resulting cumulative attenuation,
and largest neighbour difference/slew so one connected region cannot conceal a
complex envelope. The existing 10,000-interval hard resource bound remains in force.

This is processing-support accounting, not a listening-quality or intelligibility
score. Silence is not substituted as a zero loudness measurement; a gain-envelope
change over silence still belongs to its processing support. Global trim affects
the whole programme and remains separately disclosed. Show cumulative non-unity
local support separately from newly changed support; M16 does not add a configurable
cumulative-duration gate or pretend the per-request 20% limit is one.

More markers can legitimately change section peak attribution, selected gains and
transition shortening. M16 cannot promise identical proposals for arbitrary splits.
It promises no count penalty from annotation density **when processing is identical**,
and a clear explanation when processing actually differs. Existing annotation/applied-
envelope geometry safeguards remain: no automatic remapping of applied intervals or
bypass of stale geometry. Optional marker import never applies gains.

This departure from raw section counting is limited to accounting/acceptance. It
does not merge editorial annotations, choose a larger correction region, propagate
gain to neighbouring sections or search for a smoother envelope. Independent
section selection and existing transitions stay unchanged. Coordinated neighbours
and alternative segmentation remain separately reviewed work.

### Cinema example

The recorded test selected four changed sections covering approximately **89.529560%**
of the programme, including transitions. Its **20% duration gate**, not its eight-
changed-section gate, blocked that attempt. Keep this historical result unchanged.
The new default 20% still rejects the same proposed support. Increasing only the
count allowance cannot resolve that duration failure. A deliberately larger duration
allowance may allow candidate assessment, but never guarantees measured compliance
or acceptable creative impact. Do not create a cinema-specific preset or tune gains,
profile tolerances or defaults to make that file pass. Use the existing evidence and
a short synthetic equivalent for M16 regression; no repeated full cinema workflow
is proposed merely to demonstrate a less restrictive preference.

## Versioning, persistence and stale results

Introduce an immutable versioned policy record separate from `TargetProfile` and
`ProcessingPlan`. Include policy/method/accounting versions, explicit values/units,
capability flags, baseline definitions and a canonical digest. Bind generation and
Apply to that digest plus measured Source/Current identities, ordered roles/range,
current plan, annotation geometry, profile definition/applicability and runtime job
identity. The result retains the policy snapshot that actually generated it.

- Policy edits alone change document metadata and invalidate pending proposals and
  their Apply eligibility. They do not change audio revision, valid measurements,
  compliance evaluations, audition identity, or rendered audio. Reuse compatible
  baseline measurements/attribution; rerun exact candidate assessment for a newly
  generated changed proposal. Late completion from an old policy cannot replace a
  new result. Identical confirmed policy values are a genuine metadata no-op.
- Editing remains explicit and stopped/no active proposal or file transaction for
  the first UI. Show draft versus active policy; Cancel preserves the old policy.
  Nothing regenerates, applies, measures or starts playback automatically.
- New sessions and pre-M16 sessions (schemas 1–5) use configurable policy v2 defaults automatically. Preserve media, annotations and exact applied processing; discard transient old proposals. No legacy executable method or adoption prompt. Historical packages/evidence retain old decisions. Opening never rewrites the original file.
- Propose session schema 6 for the policy record. Old formats open without rewriting;
  Save As preserves an M15-readable original. Saving v2 requires an M16 reader.
  Unknown/invalid mandatory policy fields block proposal generation, without trusting
  a saved private assessment. If the otherwise valid core can be safely recovered,
  preserve media, exact applied processing, sections and usable measurements while
  marking the unsupported policy unavailable; require explicit replacement. Structural
  corruption still fails safely. Never silently substitute default limits.
- Save a coherent document/policy revision. Policy edits during saving remain unsaved;
  job/runtime identity remains separate. Keep the existing conservative cache policy.
  A new build may require remeasurement without losing settings or policy provenance.
- Reports identify active limits and, where applicable, the distinct policy associated
  with an applied proposal. Updating current policy does not relabel an earlier
  application. Refresh metadata/report approval without rerendering unchanged audio;
  actual Current and exported-artifact measured rules remain the strict M10 gates.
  Do not persist private candidate measurements/verdicts as Current authority.

No rules are invented for LRA, dialogue or metadata. Narrow selected-programme-rule
wording stays unchanged. Prediction/planning, private measured assessment and independent
Current verification remain visibly different stages.

## Acceptance tests and impact-selected qualification

Before implementation, write an explicit manifest of affected tests, shared callers,
reuse rationale and budgets. Stabilize on Sequoia, freeze one candidate and require
all listed evidence. Do not recursively launch historical milestone runners.

1. **Policy validation and identity:** every default, min/max and immediately outside
   boundary; finite values, exact units, disabled-local zero cases, inconsistent
   reserve/caps, unknown versions/fields, canonical digest, identical-edit no-op,
   deterministic request/result text and reproduction after session round trip.
2. **Gain accounting:** zero/new/manual negative/positive baselines, additional versus
   cumulative caps, repeated requests and reopen, pre-existing over-limit settings,
   full/partial/no reserve, uncovered deficit, global/local/whole-chain disclosure.
   No undocumented epsilon; display rounding never changes a gate or verdict.
3. **Extent and density:** independently computed sample-support oracle for short
   envelopes at 44.1/48/96/192 kHz; transitions/short sections/EOF; union with no double
   counting; exact duration boundary; separate versus touching regions; counts above 64, disabled/restricted count and
   structural min/max. Equal canonical processing with renamed/reordered metadata
   or redundant splits has equal support/count; truly different gains/transitions may
   differ and are explained. Inspect 10,000-row inputs without per-sample three-hour scans.
4. **Method behavior:** pre-M16 migration adopts new defaults; v2 default and
   explicit alternate policies use the same global-first/one-envelope sequence;
   already-pass preserves gain/revisions/measurements and requests no new pass.
   No truncated subset chosen merely to fit count/time budgets. Cinema-equivalent
   four-section/89.53% support still fails 20%; count-only relaxation cannot pass it.
   Passing an extent gate never fabricates a measured pass.
5. **Apply/verification and isolation:** stale policy/source/roles/range/processing/
   geometry/profile/runtime rejection; valid measurements reused after preference
   edits; no policy-induced audio identity changes. Exact candidate and independent
   Current verification; cancellations/failures keep applied gain with honest status;
   no automatic second correction. Include mono/stereo/5.1, interleaved/separate-file/
   composition, LFE and existing 7.1 capability rejection boundaries.
6. **Shared audio/session/report callers:** short shared processing/offline versus
   pre-monitor capture and export parity, seeks across ramps, EOF/duration/continuity;
   schema 1–5 migration and schema-6 round trip, legacy-reader rejection, corrupt/
   unknown policy recovery, coherent saves/late jobs; report provenance and strict
   artifact gates, single dither stage and unchanged rendered-audio reuse. Preserve
   M15 exact marker mapping and annotation-only invalidation. No playback device.
7. **Silent UI:** draft/edit/confirm/cancel/reset, automatic pre-M16 migration, disabled Apply
   on stale proposal, requested/achieved limits and failure explanations, saved policy
   and reopening stopped/device-closed/routes-muted. Proposal-limit edits do not show
   a new compliance verdict or automatically regenerate audio.

Primary affected code is SectionTrimMethod, request/result binding, session/UI and
report metadata. Measurement arithmetic, decoder/input and processing/transport
algorithms should remain byte-identical. Reuse pinned EBU/ITU, mono/RF64 and earlier
long evidence only after a dependency/hash audit; rerun affected short checks if
shared callers/build changes invalidate that rationale. New policy features are not
qualified on High Sierra by old M13 results or a local Intel-slice run.

### Resource estimates and acceptance budgets

Planning estimate, excluding implementation/fix cycles: universal build **5–15 min**,
selected automated Sequoia checks **10–25 min**, silent UI **15–25 min**. Use short
fixtures with **<=2 GiB disposable media/evidence**, plus approximately **4–8 GiB**
free build/frozen-package headroom; inspect actual free space before creating data.
Do not move or delete older packages to meet these estimates without authorization.

Keep inherited limits: one envelope/three private passes; proposal/attribution/plan
state <=16 MiB at 10,000 sections, extra working RSS <=256 MiB beyond measured-result
storage, heartbeat p95 <=100 ms, cancellation acknowledgement <=1 s, session <=128
MiB, JSON+HTML <=16 MiB, save/parse <=5 s. For new policy/support accounting alone,
require <=0.5 s at 10,000 intervals and <=32 MiB added transient memory on Sequoia;
use canonical interval/event operations, not programme-length sample enumeration.
Record observed timings/RSS and cleanup on error, not estimates as passes. Larger
count/duration settings must not increase candidate/pass limits or retain a full render.

No new 90-minute, three-hour or full cinema run is automatic. Actual High Sierra,
cross-Mac exchange, and every test loading/decoding/playing video remain deferred.
Continue explicitly named standalone/test universal builds with the established Intel
10.13 target/dependencies; never default-all, install or VST3. Combined qualification
is planned **after M18 and no later than M19**, on one frozen combined candidate and
an impact-selected manifest. This cadence is not approval of intervening milestones.

If a concrete change to macOS APIs/dependencies, measurement arithmetic, buffering,
clocks/latency or a failure warrants earlier deferred checks, explain the specific
reason before revising the plan. Coordinate any physical audio use separately.
Keep occasional picture skips accepted when deferred tests run, recording counts;
audio interruptions, accumulating drift, sustained freezes and identity/lifecycle
failures remain blockers. Full physical 7.1 listening remains separately pending.

## Complete inspection and exclusions

Every correction remains available in the bounded scrollable proposal view and JSON/HTML reporting, including accepted proposals with more than 32 or 64 rows. A shortened summary is never labelled a rejected proposal. The 10,000-interval structural bound, 16 MiB retained-state budget and report budgets remain independent safeguards. No silent subset selection.

No coordinated neighbouring adjustment, new transition rules, automatic segmentation, new strategies, limiter/compressor integration, new marker formats/EDL, component DSP or intelligibility analysis inside M16. See the separately staged UI/product priorities in ProjectIndex 58. High Sierra/Golem, video and three-hour tests remain deferred after M18, no later than M19. Main stays at M13. Stop for review before publication or integration.

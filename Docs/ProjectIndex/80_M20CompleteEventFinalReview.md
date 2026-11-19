# M20 complete-event preference — final Logic review

19 November 2026. **Design review only. No implementation or audio experiments
performed. Wait for Strategy's explicit implementation instruction.**

## Outcome and authority

Recommend proceeding to bounded implementation and focused non-audible checks after
Strategy accepts this review. No unresolved design blocker remains under the
contracts below. Detector effectiveness, memory high-water use and listening quality
remain implementation/qualification gates, not results of this review. No new DSP,
expanded search budget or relaxed delivery rule is proposed.

Reviewed the accepted amendment in [70](70_M20RegionFormationAndPredictiveGainPlanning.md),
SHA-256 `ec57b412e6eed225f9775929d5b7e1f6f19c7cf15c3d80f3566307eebe7d5364`.
Preserve [79](79_M20BoundaryAndDurationLogicReview.md), SHA-256
`b5eefd9f2bc350a121a816788ea2d1f34ead9e8721b898286c42218688311079`, as the initial
review. This final review supersedes its neutral decay preference and its
duration-first comparison between different intervention signatures. Neither
document is edited by this review.

Direct human authorization was verified in Strategy chat
`01a1ac62-28bd-7910-a131-59c3b1cea6d9`, user message
`01a1f539-8a6e-7482-866a-217548d216fa`: final review, report back to Strategy, then
Strategy may explicitly authorize implementation if it agrees. Reporting back is
authorized; this review does not authorize itself to start implementation.

## Exact starting point and preservation

Read-only local ancestry inspection establishes:

| Item | Revision / status |
| --- | --- |
| Logic checkout / branch | `/Users/fuzz/Developer/LevelScope`, `codex/m20-region-planning` |
| Current Logic HEAD | `294f9b52abc1ae21754d5cf42c879b9de674244e` |
| Accepted published UI continuation | `607a3f54a24474ebcb6bf7cb5daa59a4e64e4bf7`, `codex/timeline-selection-boundaries-2026-11-11` |
| Merge base | Exactly Logic HEAD; there are no Logic commits absent from that UI continuation |
| UI tested implementation | `c4b4620de88fb74f44016832fde8d9cba61160b6`; UI records 2,462 focused checks on each actual Mac |
| Main / locally recorded origin/main | `22961f03de7919e1de631f8b13eb4b852108ac0b`; unchanged |

Use **607a3f54… as the implementation source base**, including all M20 Logic and
accepted UI/manual-edit/session history. No competing checkout is needed. After
authorization, preserve the pending design files in a deliberate documentation
checkpoint, integrate that UI continuation into the existing Logic branch, resolve
documentation/index differences, and record the resulting combined commit before
production edits. Its source/build/test files must initially equal the UI base;
documentation changes alone must not be described as newly qualified software.
The combined commit does not exist yet. Refresh ancestry before integration if
another accepted continuation appears; never discard it or rewrite history.

Read UI documents 73–78 in the registered UI working copy and its
[11 November handover](/Users/fuzz/.codex/worktrees/a411/LevelScope/Docs/TimelineSelectionBoundaries-2026-11-11.md).
Preserve the fixed remeasurement-preview lifetime, keyboard/Undo focus, separate
selection, shared pending boundary batches, exact transition shading, narrative
read bound, and device preferences. The UI Source tree is
`e3cd3ab769eb0e49821e016f3f28e40a3a3ae202`; current Logic's older Source tree is
`ab0819224ef725f3ebc482c7cc5356213b678feb`. Do not start from the older tree alone.

This review performs no fetch, branch integration, source edit, build, audio decode,
device opening, Golem access, experiment, UI dispatch, commit or publication. All
25 local checkpoint refs are preserved; the sorted `for-each-ref` snapshot of
`refname objectname peeled-objectname` has SHA-256
`4f7c31708f633fc13f09a30b4f331cffc782973559e6fa0bd2b972a76a97c6d9`.
Remote freshness was not re-established by this design-only review.

## 1. Evidence and coverage contracts

Keep three different inputs: frozen valid **Current measurements** establish need;
verified **Source PCM after assembly/composition, before whole-programme processing**
supplies event evidence; the unchanged Current prediction sidecar estimates complete
candidate outcomes. Keep a matching Current energy envelope to expose old automation.
Never infer Source by dividing processed samples by gain.

Freeze event groups and their need associations before scoring any plan. Each group
has a stable ID, Source onset bracket, last-attack/decay-end bracket, channel witnesses,
confidence, contributing measured-window IDs and any missing/capped context. An
intervention signature identifies the need-bearing groups actually changed by the
candidate's **additional local envelope**. Global gain alone does not create a local
signature. Groups touched incidentally by an entry, return or internal change are
also inspected and disclosed; they cannot disappear merely because their loudness
window centre lies outside the changed support. Unexamined incidental evidence
prevents a supported-complete claim for that plan.

For frame `f`, use `D(f) = beforeLocalDb(f) - candidateLocalDb(f)` from the exact
canonical processors. The coverage test concerns additional correction, not the
candidate's absolute local gain or global normalization. A pre-existing fade cannot
satisfy a newly proposed correction's coverage obligation.

The **intended correction** is the selected positive gain alternative for this
candidate/part, after gain caps and deterministic representable-plan construction.
It is not the uncapped statistical request or unachieved planning reserve. A capped
3 dB treatment can cover an event completely while leaving estimated loudness deficit;
the complete programme measurement decides whether that plan satisfies requirements.
Do not choose an epsilon floor after construction merely to call a late ramp complete.
Persist the selected intended correction schedule before evaluating its coverage.

A supported complete group requires all of the following:

- Corroborated, complete event evidence with a bounded onset and supported decay end.
- Its intended nonzero correction is reached no later than the earliest plausible
  onset frame, and is maintained through the upper end of the supported decay bracket.
  A ramp entirely after the first hit fails this check, even without crossing an attack.
- Return toward the preceding gain starts only after that bracket. An in-decay return
  is a known partial treatment when supported evidence establishes that the decay
  continues; an uncertain decay is not automatically known truncation.
- Sustained measured need can justify internal changes between positive gains.
  Their fixed ramps must follow the recorded intended schedule and cannot hide a
  zero-gain hole inside one continuing event. A justified spared middle requires
  separately frozen event/need structure, not post-hoc redefinition by the candidate.
- Exact additional-envelope support respects permissions, pins, protected intervals,
  gain/duration/separation limits and fixed transition behavior. Validation uses all
  relevant old/new ramp breakpoints and frame endpoints, not plot pixels or labels.

Do not add a floating-point allowance to these geometry checks. Derive representable
reference gains with the same canonical construction and compare exact stored
values/processor behavior. Recheck shorter-radius and old-ramp interactions after
canonicalization. No gain/coverage test relaxes measured loudness or true-peak limits.

### Per-group and whole-plan states

Store independent `knownPartial` and `evidenceIncomplete` flags plus reasons; do not
discard one when both occur. A group can have a proven missed entrance and an unknown
ending. Complete means neither flag and all affirmative checks above succeeded.

| Groups touched by a local candidate | Whole-plan coverage / selection treatment |
| --- | --- |
| Nonempty need-bearing signature; every implicated and incidental group supported complete | **Supported complete**; eligible for preferred recommendation if measured rules pass |
| At least one known partial, no incomplete evidence | **Known partial**; explicit-review alternative only |
| No known partial, some incomplete/uncertain evidence | **Uncertain**; explicit-review alternative only |
| At least one known partial and some incomplete evidence, in the same or different groups | **Known partial + uncertain**; preserve both; explicit-review alternative only |
| Empty signature, zero correction, or no supported association with correction need | **Not an event-treatment candidate**; cannot earn completeness or enter the event preference ranking |

Untouched need groups remain visible in the overall report but are not added to a
plan's signature. There is no reward for the number of groups corrected. Extra local
processing without a supported need association is not a way to improve the score.
Missing feature capacity is uncertain evidence, never a favourable zero-risk value.
Already-passing no-op and measured sufficient global-only paths keep separate labels
and authority; they need no fictitious event classification.

## 2. Deterministic whole-set selection

First validate legality and separate actual measured passes, actual failures and
unavailable/incomplete assessments. Forecasts cannot enter the passing set.

Across **all measured-passing local plans**, collect the supported-complete subset.
If nonempty, select only from it, regardless of intervention signature. Thus a known
partial 3 s treatment cannot defeat a supported complete 6 s treatment through the
old duration-first comparison. For that subset use the lexicographic order:

1. Number of newly introduced audible internal transitions unsupported by a frozen,
   sustained correction-need change. Do not count annotation dividers or untouched
   existing transitions. Unknown justification remains explicit, not zero unsupported
   transitions with complete evidence.
2. Actual affected-frame union, including old/new transition differences.
3. Integral of absolute additional local dB over frames.
4. Maximum additional local attenuation.
5. Number of new audible internal transitions, including justified ones.
6. Absolute global-gain delta.
7. Lexicographic source-frame signature and canonical plan ID as final stable ties.

Use complete-set stages and exact comparisons, not pairwise epsilon comparators.
Signature is explanatory and a final tie-break; it cannot partition away the global
complete-before-partial rule. No event-count or artistic-quality score is introduced.
Do not continue optimizing LRA downward once every selected rule passes.

Keep passing partial, uncertain and mixed alternatives in separate labelled lists.
Within a list use documented preservation ordering for stable presentation, showing
unknown transition evidence; list order is not an artistic recommendation. Do not
prefer uncertainty over known partial treatment or vice versa through a hidden ordinal.
If no supported-complete plan passes, **recommendedAttempt is empty**. A diagnostic
displayed attempt may exist, but Apply is unavailable until explicit selection below.

The existing anchored 0.01 LU failed-result shortlist and remaining-rule ranking
remain unchanged. Coverage never promotes a measured failure, widens that band or
changes failure ranking. Report the best failure separately from passing alternatives.

## 3. Explicit fallback selection, identity and Apply

Current code only permits the retained attempt to Apply; the review-panel attempt
menu is inspection-only. A new arbitrary array index must not bypass that safeguard.
The following small interface extension is required before exposing fallback:

- Separate immutable `recommendedAttempt`, optional `selectedAttempt`, and
  `inspectedAttempt`. Keep the complete immutable attempt list and its result ID.
- Add outcome **Review required — passing alternatives, complete treatment not
  established**. It is not `ready` or a measured failure. Existing `noChange`,
  global-only and legacy method behavior stay separate.
- Inspection and timeline navigation never select Apply authority. Add an explicit
  **Use this measured alternative** action for an actually passing attempt. Show
  the coverage category, missing/partial reasons, exact gains, full measured result
  and displaced alternatives at that action. A passing complete alternative can
  also be explicitly chosen; no failed or unavailable attempt can be selected.
- Selection performs no processing, playback or measurement. It asks the service
  to create a new immutable selected result bound to the current request and the
  exact measured plan. Only then can the ordinary stopped-only Apply be enabled.
  Selection is the explicit acknowledgement; a second generic modal is unnecessary.
- A fallback remains labelled **Selected passing alternative — known partial /
  uncertain coverage** in main and timeline. Never rename it a successful complete-
  event recommendation. Retain the reason in applied history.

### Preserve measured authority without retaining three window populations

Use a service-owned, immutable compact `RegionAssessmentReceipt` for v2 attempts.
Create it only after a complete qualified-worker result has passed context/content/
domain/metric validation. It records the exact plan ID, input/range/ordered roles,
request binding, worker snapshot/method/build/platform/origin identities, all relevant
metric statuses/values, and selected profile/requirements/applicability digests.
Include the compact per-channel sample/true-peak values and exact locations, plus
gate/percentile/stability scalar provenance needed by existing comparison reports.
Keep the current all-channel, full-precision independent-remeasurement assertions;
adapting a v2 consumer must not simply remove them because window arrays are released.
Reevaluate those measured scalar metrics under the **same exact bound requirements**
at selection and Apply; reject any binding change rather than silently refreshing a
stale proposal. Actual full analysis is how the receipt is created, not how it is
reconstructed from saved text or UI-supplied values.

The receipt has no public constructor/decoder/fabrication path. A selected-result
factory accepts an attempt ID belonging to that service result, not caller-supplied
metrics or a mutable `RegionAttempt` row. Historical scalar records remain documentary
and cannot be converted into this runtime authority. This requires a version-specific
M20 v2 validation path; it must not weaken legacy `evaluateSnapshot` checks.

Release private full window populations after producing bounded attempt evidence;
do not keep another population for each selectable alternative. Adapt v2 consumers
of `SectionProposalResult::assessment()` (Apply and history capture) explicitly;
legacy methods retain their existing snapshot path. The v2 receipt cannot be installed
as Current, diagnostics, a measurement cache, or exported-artifact verification.

The selection token binds result/request identity, attempt/plan/receipt identity,
coverage/evidence digest, profile/requirements, selection origin and runtime session/
job generation. Proposal settings, source/member identity, roles/range, processing,
requirements, pending manual edits or New/Open invalidate it. No late job may restore
it. Display inspection alone does not invalidate it, but Apply continues to guard
against displaying a different diagnostic plan as if selected. Device routing,
monitor gain and video offset stay outside programme measurement identity.

Apply rechecks stale state, receipt authority, legality and exact chosen plan; builds
the complete bounded history before mutation; then applies once and marks Current
UNVERIFIED. Independent Current measurement and M10 finalized-artifact checks remain
mandatory. No automatic rollback, retry, attenuation or publication is added.

This is a required implementation dependency, not a claim that today's UI can already
select fallback. Until service and UI guards pass focused tests, fallback remains
inspection-only and cannot Apply.

## 4. Adaptive three-assessment schedule

The new explicit preference uses comparative planning; it must not inherit the v1
`compare=false` early stop on the first partial pass. Require effective comparison
enabled for schema-2 complete-event mode and show it; do not change old settings.
Keep original-boundary, supported-complete and justified tighter/spared-middle families
in the existing finite predictive pool, with existing round-robin and deduplication.

Deterministic scheduling after predictions, always at most three actual attempts:

1. If a legal supported-complete candidate forecasts IL/LRA and any selected M/S
   within their bounds, with those estimates finite and available, assess the best
   such candidate first. This is screening eligibility, not a predicted verdict.
   Unknown true peak remains unknown;
   it is not a forecast failure. Otherwise assess the best original-boundary candidate
   if available, then the best legal remaining candidate if that family is unavailable.
2. After a failed complete-family measurement, an eligible exact global-feedback
   correction under the existing measured IL/TP/M/S intersection gets **at most one
   priority retry**. Preserve its local plan; recompute coverage and binding anyway.
   If no feasible scalar correction exists, do not manufacture one. Other local/gain
   alternatives remain ordinary finite-pool candidates, not unlimited feedback passes.
3. Without that priority retry: if no supported-complete pass exists, next assess the
   best distinct unmeasured supported-complete candidate, if any. Do not discard all
   such candidates solely because an approximation predicts failure; use the existing
   anchored failed forecast ranking to choose among them when necessary.
4. If a complete pass exists, reserve the next available slot for an unmeasured original
   challenger, then a distinct justified tighter/spared-middle/other candidate. If
   no complete candidate exists or remains, use the unmeasured original then the best
   remaining finite-pool candidate. Stop on exhausted candidates, cancellation or budget.

The second and third slots apply the same rules. Never stop at a partial pass while
a legal unmeasured complete candidate and a slot remain. Never promise all three
families plus a repair. For example, complete-fail / corrected-complete-pass / original
uses the whole budget and displaces the tighter comparison. Record that displacement.
If all three measured plans fail, none can Apply. If only partial/uncertain plans
pass, enter explicit review rather than silently recommending them. Cancellation,
incomplete reads or watchdog termination provide no new Apply authority.

Within a scheduler-selected family, order screening-eligible forecasts by the
coverage/preservation tuple above. Otherwise retain the existing anchored forecast
failure ordering and canonical plan-ID ties. Missing required estimates are not
screening passes: keep a deterministic legal challenger when all are unavailable,
label that choice, and measure it if its family receives a slot. A loose true-peak
upper bound above the ceiling remains inconclusive, not a measured failure.

At most 8 roots, 12 nonmandatory location nominations per root, 256 constructions
per root, 16 retained geometries per root, beam width 4 and 64 complete proxy plans
remain the caps. All complete, original, gain/global and repair variants share them.
No fresh budget per signature, child, family or user-inspected attempt.

## 5. Feature definition and bounded extraction

Accept document 79's explicit detector as the **initial evaluation hypothesis**:
continuous auxiliary first-order 20 Hz HP, 160 Hz LP low-band branch; 20 ms broadband
and 80 ms low-band energy windows; 10 ms hop with controlled 5/20 ms comparisons;
-80 dBFS RMS activity floor; 250 ms preceding median; 6 dB broadband rise/20 ms
persistence and 3 dB low-band rise/40 ms persistence; two kinds of corroboration
within an initial 80 ms bracket; 1 s grouping/lower-need persistence; 100 ms last-
attack reference; 6 dB fall/200 ms in-decay nomination and within 3 dB of background/
300 ms after-decay nomination. Exact formulas, silence cases and channel rules are
in 79 and are part of the serialized feature-definition digest, not file-specific knobs.

Define bracket coordinates conservatively: onset includes the beginning of the
supporting fast-energy window through the confirming cell; the supported after-decay
end includes the full confirmation interval. Use actual source sample frames and
partial-cell lengths, with hop edges `floor(k*Fs*hopMs/1000)`. The 80 ms agreement
test is not onset accuracy. Return-after-decay uses the upper bracket end; in-decay
is a separate candidate. Clipped context, continued activity at EOF, weak/conflicting
channels and missing detail stay uncertain. No padded measurement tail creates an event.

Retain per-channel energy before aggregation, non-LFE role-labelled mean and contrast
witnesses, independent LFE context, and no waveform summation. Finite over-full-scale
Source samples remain analyzable. LFE evidence cannot create correction need, but
linked processing/peak measurement continue across all channels. Keep capacity 12;
qualify this method only for supported stereo/conventional 5.1. No new loudness weight.

Use the existing auxiliary verified-Source read: collect event evidence before
`ProgrammeGainProcessor`, then Current context/prediction after it. Stream all roots
and context with continuous filter state; no extra PCM file or decoder pass. Freeze
coarse need nominations before that pass. Use bounded rolling histories for new
feature nominations; do not promise a post-hoc high-resolution neighborhood already
discarded. If a new nomination lacks the required history/detail, mark it uncertain
and report the omission. No hidden second pass or cadence reduction.

Memory feasibility is plausible, **not yet measured**. Maximum detail is 14 six-second
neighborhoods per root for eight roots (672 s total). At 5 ms and 64-byte compact cells,
this is about 8.21 MiB. A **10 MiB sub-budget inside 64 MiB total** must also include
ring/history capacity, headers, per-channel witnesses and up to 64 group summaries
per root. Whole-programme high-rate feature arrays are prohibited. Compact/stream
deterministically; incomplete evidence must never earn complete status.

Existing three-hour-sized analysis plus sidecar evidence is 35,846,608 bytes, before
all prediction scratch/search/metadata. It does not prove combined fit. Preflight and
track actual reserved-capacity/high-water costs for features, diagnostics, bins,
transformed windows/prefixes, finite gain alternatives, beam/pool, IDs/strings, plans,
receipts and explanations, including simultaneous lifetimes and reallocations. Count
allocations once when shared; bound them before allocation, not after copying. Shared
pre-existing authoritative Source/Current caches and one in-flight worker result are
reported separately from additional planning memory; do not duplicate them to evade
the cap. Retained candidate windows are not free baseline data.

Keep cached planning targets <=10/20 s on Sequoia/High Sierra, 30/60 s watchdog,
<=1 s cancellation and native p95 <=100 ms responsiveness. Measure extraction time
separately: the new per-sample filters add real work even without another decode.
If the accepted budget cannot fit, fail clearly and return for review rather than
raise the limit or hide missing evidence. No performance improvement is claimed yet.

## 6. Exact version and session contract

The accepted UI base still uses session schema 9, programme schema 9, region settings
schema 1, history attachments 1/2 and the v1 search. The earlier allocations remain
available; allocate explicitly as follows:

| Definition | New explicit version / preserved behavior |
| --- | --- |
| Region settings | schema 2; explicit `prefer-complete-event` preference, feature definition/digest, automatic-margin permission and comparative scheduling |
| Constraints / formation / search | `region-constraints-v2` / `region-formation-v2` / `region-gain-search-v2` |
| Event features / coverage | `event-boundary-features-v1` / `event-treatment-coverage-v1` |
| Passing comparison / scheduling | `candidate-review-rank-v3` / `event-assessment-schedule-v1` |
| Runtime receipt / selection binding | `region-assessment-receipt-v1` / `region-selection-v1`; not session caches |
| Session document | schema 10 for the new writer; programme representation stays schema 9 |
| Applied explanation | attachment 3 for v2; read legacy attachments 1/2 unchanged |
| Unchanged | `window-energy-proxy-v1`, qualified meter/profile digests, `section-gain-v1`, proposal-policy-v2 accounting and existing audio/cache formats |

Use an explicit version registry/dispatch, not a global replacement of the current
`regionSearchVersion` constant. Schema-1 encode/digest/decode and v1 generation must
retain exact legacy method strings, inward-only automatic permissions and ranking.
History validation currently compares against that single constant; update it to
validate each supported history/settings combination without rewriting old digests.
The UI's large-narrative read-bound fix must survive attachment-3 support.

Schema-2 definitions embed the full feature constants and accepted enums; reject
unknown/mismatched mandatory fields. All parameters, including fixed filter/grouping
hypotheses and scheduler, participate in the request/evidence identity. No invisible
per-file tuning. Altering a hypothesis requires a new definition digest and fresh
affected evidence, even during development.

Read session schemas 1–9 through their existing migrations; add the explicit legal
pair `(document=10, programme=9)`. Do not equate these counters or widen all old pairs.
Opening an old session does not opt it into new settings or rewrite its file. Saving
with the new writer requires a new reader; keep explicit Save/Save As messaging and
preserve the original until Save. Unknown mandatory core/settings fail safely;
independently damaged/incompatible measurement caches remain replaceable. Historical
explanation damage may discard history with a clear notice, never the applied audio.

Persist exact applied settings and compact historical evidence: event brackets,
feature/requirement/policy versions, coverage category, chosen alternative, explicit
fallback reason and baseline/result provenance. Preserve the total 8 MiB history,
4 MiB core and 128 MiB session bounds. Never persist runtime selection authorization,
private candidate window arrays or verdicts as Current verification. Reopening keeps
the existing stopped/device-closed/muted contract and conservative cache compatibility.

The 60% experimental configuration uses the existing policy schema, with its full
definition embedded. It is neither an automatic mode default nor a settings migration.
Keep 20% conservative; additional 3 dB, cumulative 12 dB and reserve 0.5 dB unchanged.
Source/Current audio measurements can be reused after planning-preference changes;
proposals and their selection tokens become stale. Export gates are never relaxed.

## 7. Focused implementation and qualification sequence after dispatch

Logic owns feature collection, formation/coverage/ranking, measured receipts and
selection validation, versioned settings/history and focused tests. UI owns new
controls, review selection and presentation after a focused-tested interface
checkpoint. Do not concurrently edit the UI worktree. Until integrated UI fallback
tests pass, no product control may activate fallback selection. A later UI handover
must include the result/selection state table above, not just new enum names.

| Fresh tests required | Decisive checks |
| --- | --- |
| Coverage | Early/full onset, ramp entirely after attack, return before last attack, in-decay versus after-decay, shortened ramps, old automation, internal nonzero schedule, zero holes, incidental events; exact half-open frame boundaries and no epsilon pass |
| Whole-set ranking | Complete wins over shorter partial across different signatures; complete+unknown, complete+partial and partial+unknown mixtures; no empty/zero/unrelated reward; permutation/tie determinism; no uncertainty ordinal; legacy failed ranking unchanged |
| Selection/Apply authority | Inspection cannot Apply; explicit passing fallback can; forged/copied documentary receipts and wrong plan/index rejected; stale source/roles/requirements/settings/manual batch/New rejected; no saved receipt authority; applied state remains unverified until independent measurement |
| Scheduler | Complete fail/repair/original, original/complete/other, no complete candidate, only uncertain passes, exhausted budget, cancelled third pass; all variants count; no first-partial-pass early success; displaced comparison explained |
| Detector | Synthetic known entrances/tails, silence/DC/noise, quiet/gradual/sustained/dialogue, repeated hits, antiphase stereo, reordered 5.1/LFE, finite over-full-scale composition; selected rate/cadence cells and chunk invariance |
| Geometry/limits | Frozen midpoint ownership even after pruning, no borrowing/stacking, explicit margins versus context, pins/protection, per-original-root split cap, exact old/new duration union, reserve/caps and separated global gain |
| Existing workflows | Focused legacy M19/M20, manual editing/Undo/boundary batch, remeasurement-preview lifetime, session/history read/write/migration, measurement sample parity and applicable export/selection invalidation |
| Resource | Worst three-hour-sized model at 5 ms, eight broad roots, dense nominations/summary overflow, all reserved capacities, receipt lifetimes, cancellation/watchdog/UI heartbeat; no three-hour media decode by default |
| Real examples | Current listening cases are development material; evaluate first-hit/decay and coherent later passage/spared middle without hard-coded coordinates; existing simple/heterogeneous/5.1/gating and bounded cinema failure; fresh holdouts only after freeze |

Run targeted Sequoia checks first and the established universal/Intel deployment
build. Freeze the combined Logic/UI candidate for selected actual High Sierra and
cross-platform checks; do not equate cross-build success with runtime qualification.
Use explicit test groups; never recursively invoke all historical milestone runners.
Reuse unchanged official-reference/processor/video/long-media evidence only after
source/dependency inspection. Repeat affected references if authoritative arithmetic
changes; selected video or sustained PCM runs require a concrete timing/buffering/
dependency change or observed problem. Physical device opening/listening is separately
coordinated. Occasional picture skips retain their accepted conditions; audio loss,
drift, sustained freezes and lifecycle/identity faults remain blockers.

Separate evaluation variables: short 5/10/20 ms detector comparison; old/new geometry
with identical zero margins and 20%; new-mode zero/1 s permission comparison; then
20/40/60% from the same baseline/permissions. Predeclare a finite run list instead of
all combinations. No requirement that the cinema excerpt reach 15 LU. Freeze methods
before fresh holdouts/listening, and retain all historical results unchanged.

Provisional resources remain those in 79: 20–45 min named build with existing caches,
10–25 min targeted development checks, roughly 30–60 min frozen offline checks per
Mac excluding native/manual/listening time; 5–10 GiB new build/package and <=2 GiB short
fixtures. These are unmeasured planning estimates, not promised qualification times.

## Final conditions and report-back

No material departure from the accepted amendment is requested. This review settles
the previous open selection/scheduling questions; the compact runtime measured receipt
is the required bridge between selectable alternatives and the existing memory/Apply
contract. Its authority isolation and history compatibility must pass before fallback
is exposed. The detector, one-pass bounded retention and resource fit remain empirical
gates; if they fail, report the specific issue rather than quietly weakening coverage.

Strategy should inspect this exact review and, if it agrees, explicitly dispatch the
bounded implementation. Until then, leave production files and branches unchanged.
Review/documentation verification is not qualification; no new tests, audio runs or
listening judgments are claimed here. Publication/main integration remain outside
the implementation-start authority.

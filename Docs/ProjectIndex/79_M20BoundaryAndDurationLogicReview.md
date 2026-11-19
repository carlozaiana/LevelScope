# M20 boundary and duration amendment — Logic review

19 November 2026. **Design review only. No implementation, changed defaults, new
audio experiments, publication or merge is authorized.** Proposed constants below
are evaluation hypotheses, not validated detector performance or product defaults.

Reviewed Strategy's [document 70 amendment](70_M20RegionFormationAndPredictiveGainPlanning.md#boundary-and-duration-amendment-for-logic-review),
SHA-256 `bd5537044d56f9eeaf6e861a0d2a57660370b247c63dd441a70b806ba22f64dd`.
That document is preserved unchanged by this review. The Logic checkout remains at
`294f9b52abc1ae21754d5cf42c879b9de674244e`; main/origin/main remain
`22961f03de7919e1de631f8b13eb4b852108ac0b`. Strategy's three pending documentation
edits are preserved. This review adds this document and an index pointer only.

## Recommendation and necessary qualifications

Proceed to a separately approved, bounded boundary-method experiment, with these
design corrections. Separating correction need, inferred event extent and actual
gain geometry is the right direction. A 60% experimental duration ceiling is
compatible with the existing accounting; it does not solve boundary placement.

1. **An event detector cannot identify the editorially correct phrase.** Energy
   rises, attacks and decays can nominate alternatives. They cannot establish that
   successive hits belong to one musical idea, identify the last *relevant* hit
   unambiguously, or certify an unobjectionable edit. Use evidence brackets and
   explicit uncertainty, plus listening and user constraints.
2. **The automatic permission change must be explicit and versioned.** Current
   `RegionPlanningSettings::validationError()` rejects automatic margins. Do not
   reinterpret schema-1 settings, or equate reading context with permission to
   process it. Outward permission remains off unless explicitly confirmed.
3. **Ranking must not reward processing more events.** Compare boundary treatment
   for the same implicated event groups before comparing preservation costs. Do
   not build a global count of events covered, or let low confidence appear as
   zero boundary risk. A deterministic proposal for this is below.
4. **Three exact assessments cannot test every combination.** Boundary, gain,
   duration and global-trim variants compete for the same budget. Reserve diverse
   challengers, disclose omissions and accept an inconclusive comparison. No
   hidden fourth measurement or assertion of optimality.
5. **The listening examples are now development material.** The prior independent
   review remains historically valid, but these files cannot independently validate
   a method tuned to their entrances/exits. Freeze the new method before another
   holdout review. Do not make the cinema excerpt reaching 15 LU an acceptance gate.
6. **Integrate against the accepted UI continuation before future coding.** The
   UI checkout is at `607a3f54a24474ebcb6bf7cb5daa59a4e64e4bf7`, on
   `codex/timeline-selection-boundaries-2026-11-11`, with later selection/manual-edit
   and session fixes absent from this Logic branch. Read UI documents 73–78 and
   preserve that history; agree the combined starting commit after implementation
   approval. This review does not merge or switch either checkout.

## Existing evidence that constrains the design

The exact `listening-gate-04` sessions use the same 275.6 s stereo programme,
48 kHz, with content SHA-256
`7ba6404041c94c909df290a6829a0fe92ce594329a1d575490c2bb4a5b6cc8f0`.
Their caches and historical narratives were inspected read-only in this Logic chat.
The baseline is IL -20.523327351947 LUFS, TP -0.035707558335 dBTP,
LRA 19.042928894183 LU, P10 -35.830809388641 LUFS. The selected 18 LU
ceiling gives the original planning threshold -17.830809388641 LUFS.

- The 12.7–15.7 s window is -17.894106343 LUFS; the 12.8–15.8 s window
  is -17.709045292 LUFS. That crossing creates the first automatic root at
  12.8 s. M20's actual entry ramp is 12.8–12.9 s, nominal cut 12.85 s;
  M19's is 12.75–12.85 s. The user's approximately 12.283 s first hit lies
  outside both. A 1 s left margin can contain a 100 ms ramp ending before that
  landmark, but the exact attack bracket remains to be measured and reviewed.
- The first root's 30.1 s endpoint is a loudness-window endpoint. There is no
  newly established exact last-hit/decay endpoint. Do not invent one from the
  user's observation that the correction returns too soon.
- The initial 49.9–60.2 s root is a union of high windows. M20 retained
  50.0–53.1 s support, including ramps, through a bounded numerical search.
  The user's 52.615–57.762 s landmarks are an evaluation hypothesis. They must
  not become hard-coded detector answers, gain boundaries or acceptance tolerances.
- At 20%, the duration budget is 55.12 s. M19's 51.5 s plus the omitted
  175.1–180.1 s correction and ramps would be about 56.6 s. M20's retained
  43.9 s leaves enough capacity for another five seconds; its omission cannot
  be explained solely by the final duration budget. The search stopped at a pass.
- M20's first complete candidate failed IL by about 0.000143 LU, despite a
  passing forecast. The second passed. Retain this strict full-precision example
  when finer boundaries increase prediction sensitivity.

The [M19 record](/Users/fuzz/.codex/worktrees/a411/LevelScope/build/ui-inspection-2026-11-07/listening-recovery-02/MusicInBackAndForegroundSomeDialog/m19/proposal.txt)
and [M20 record](/Users/fuzz/.codex/worktrees/a411/LevelScope/build/ui-inspection-2026-11-07/listening-recovery-02/MusicInBackAndForegroundSomeDialog/m20/proposal.txt)
match the saved historical narratives. This review performs no new decode,
measurement, rendering, playback or listening experiment.

## Proposed detector domains and inputs

Maintain three independent authorities:

| Purpose | Input and treatment |
| --- | --- |
| Correction need and final acceptance | Existing validated frozen Current M/S/LRA populations, gate membership, profile/job requirements and authoritative measurement worker. No new loudness metric or gate. |
| Event evidence | Verified original programme PCM, **after channel assembly/composition but before whole-programme processing**. Preserve channel order; include a matching inexpensive frozen-Current energy envelope for comparison. An existing processing ramp alone cannot supply source-event evidence. |
| Candidate prediction | Existing 100 ms Current energy sidecar and full-window overlap model. Finer event features do not silently change `window-energy-proxy-v1`, measurement arithmetic or its accuracy claims. |

Use the existing single auxiliary streaming read. Feed the dry programme features
before `ProgrammeGainProcessor`, then the unchanged prediction collector after it.
Original here means the confirmed composition when stems are present, not a different
unmixed file. Do not reconstruct Source by dividing processed samples by gain: that
would lose original float-stage details and is unnecessary with the verified copy.

Compute per-channel energies before aggregation, never a waveform sum. Preserve
antiphase stereo and reordered 5.1 interpretation. For loudness-related context keep
the qualified weights and LFE exclusion. For non-loudness attack evidence use
role-labelled per-channel activity and the unweighted mean of non-LFE channel
energies; retain the maximum eligible per-channel contrast as a separate witness.
Do not include below-floor channels in the contrast witnesses or change the mean's
channel denominator from cell to cell. An isolated channel rise may be visible even
when the other channels mask it. LFE is separate evidence:
it may flag a transient that a whole-mix ramp would cross, but cannot create an LRA
root or justify a new correction by itself. Applied gains and all-channel peaks
continue to include LFE. Arrays remain count/role driven, with capacity 12; this
method's qualification stays stereo/conventional 5.1.

### Initial feature definition for evaluation

These are proposed, explicit heuristic settings. They are not standards, proven
onset accuracy, or new user-adjustable DSP parameters.

- Accumulate in double precision with finite/overflow checks. Use a continuously
  running, auxiliary first-order 20 Hz high-pass per source channel to suppress DC;
  use its output for broadband energy. Add a first-order 160 Hz low-pass on that
  output for low-band energy. Specify coefficients by
  `a = exp(-2*pi*fc/sampleRate)`, HP `y[n]=a*(y[n-1]+x[n]-x[n-1])`, LP
  `y[n]=(1-a)*x[n]+a*y[n-1]`. This is deliberately a modest heuristic band, not a
  calibrated bass-drum detector. Record the coefficients and filter version.
- Process from programme sample zero without resets at roots or chunk boundaries.
  Zero initial state is part of feature identity. No padded EOF samples count as
  events. The filters can have tails; do not subtract a supposed constant delay
  and label the result an exact onset.
- Proposed fast broadband energy window: 20 ms; low-band energy window: 80 ms.
  Keep existing measured 400 ms momentary windows as sustained context with their
  real supports: Source for corroborating Source events, frozen Current for linking
  correction need. Do not reinterpret their 100 ms sampling as onset precision.
  Retain a raw-sample maximum/channel label to corroborate activity, not to place
  every gain edge at a peak or zero crossing.
- Primary hop is 10 ms. Compare 5 and 20 ms with identical window/filter/threshold
  definitions. Boundaries are `floor(k * sampleRate * hopMs / 1000)` from source
  zero, not repeated rounded increments: 5 ms at 44.1 kHz is 220.5 samples.
  Carry actual first/end frame and sample count for every partial cell/window.
- Planning values are RMS/energy ratios, not LUFS. Exact silence is unavailable
  for log ratios; never replace it with a measured zero dB value. A proposed
  numerical activity floor of -80 dBFS RMS, combined with local contrast, prevents
  near-zero ratios from looking significant. Its influence on very quiet passages
  must be tested; it does not change measurement availability or acceptance.
  A genuinely silent pre-window uses a separately labelled onset-from-silence test
  (activity above the floor with the same persistence), not an infinite contrast.
  A silent background's return test is sustained energy below the floor; its ratio
  remains unavailable. Keep these heuristic classifications distinct from measured
  silence/below-gate statuses.

### Entry, exit, grouping and confidence

Propose this reproducible state machine for review:

1. Freeze coarse roots and need anchors from Current before feature collection.
   Scan their source intervals plus at most 3 s context per edge. Also preserve
   source evidence around existing internal nominations. Continuous filtering
   outside retained contexts is history preparation, not expanded processing.
2. Nominate a broadband rise when 20 ms energy rises by at least 6 dB relative
   to the median of the preceding 250 ms feature cells and remains elevated for
   at least 20 ms. Independently nominate a low-band rise of at least 3 dB with
   40 ms persistence. The active signal must clear the activity floor. Freeze
   both domain/role baselines before each test; no per-file threshold tuning.
3. A **corroborated** entrance needs source evidence, association with an existing
   need anchor, and agreement of two kinds of evidence: broadband/low-band rise,
   or one rise plus sustained Source momentary/fast-energy elevation. A raw peak alone
   is insufficient. Agreement within an initial 80 ms bracket is a test hypothesis,
   not a claimed timing error. LFE-only evidence remains contextual. A rise present
   only in Current and coincident with an old processing ramp is flagged as
   processing-associated, not promoted to corroborated Source evidence.
4. Record the first qualifying rise cell and its preceding window support as an
   onset bracket, rather than one supposedly exact attack sample. Retain earliest
   activity and raw peak positions within that bracket. Values close to detector
   thresholds, clipped context, conflicting channels or weak corroboration produce
   **tentative** evidence. Missing/capped features produce **unknown**, not tentative
   success. Confidence is categorical evidence provenance, not a probability.
5. Group successive attacks when a gap is shorter than 1 s, or the sustained
   context has not returned toward its preceding background. Do not let every
   inter-hit dip reset a group. A quieter interior becomes a split candidate only
   after both the existing correction-need rule and at least 1 s of sustained
   lower-energy context support it. Retain both complete-group and spared-middle
   alternatives; do not equate grouping with a musical phrase.
6. For the last corroborated attack, freeze the median energy over the following
   100 ms as its local reference (mark truncated/missing context). Nominate an
   **in-decay** return when energy has fallen at least 6 dB below that reference
   for 200 ms; nominate an **after-decay/background** return after energy remains
   within 3 dB of the pre-group background for 300 ms. Apply checks to broadband
   and an active low-band channel; an ongoing tail in either is disclosed. A new
   grouped attack resets the last-attack reference. These alternatives require
   listening; neither has automatic artistic priority. Never extend beyond U to
   satisfy them. If no return is found, report an open/uncertain ending.
7. Sustained material can have no distinct attack or decay. Use sustained context
   brackets and the original geometry challenger with tentative evidence. Do not
   force a peak-aligned edit or claim a complete event was detected. EOF clips an
   open group and is labelled; it is not an inferred decay to silence.

Association with correction need is required throughout. An event group is relevant
when it overlaps measured above-threshold window centres in H, or is the immediately
preceding/following activity connected to that group within the bounded context.
No chains of events outside H extend the search indefinitely. Whether the first/last
hit should receive the same gain as the entire passage remains a listening judgment.
The requested alternative can express that intention; the detector does not prove it.

## Processing geometry and permission ownership

For an alternative intended to fully attenuate the first attack, the **additional
correction** must reach its requested value no later than the onset bracket's
earliest plausible frame. For a normal 100 ms ramp this places its start roughly
100 ms earlier and its nominal cut 50 ms earlier. Use the actual processor radius,
sample endpoints and before/after difference, not this arithmetic approximation,
to validate the final canonical envelope. Near programme zero or existing short
intervals, an appropriate ramp may be impossible; disclose that instead of padding
audio, moving a pin or changing the processor.

Inspect returns and internal changes by the same rule. A return before the last
corroborated attack has completed is a boundary conflict; returning during a decay
versus afterward is an explicit comparison, not automatically a conflict. Existing
ramps can change radius after canonicalization: re-evaluate the complete old/new
envelope and event treatment after construction. Do not judge geometry solely by
the proposed nominal interval labels.

**Recommend fixed ownership of overlapping automatic margins:**

1. Sort frozen, disjoint H roots by source begin/end, with stable root ID last.
2. For adjacent H intervals, form a separator at
   `floor((left.H.end + right.H.begin)/2)` using overflow-safe arithmetic.
   Half-open ownership assigns the separator sample to the right side.
3. Intersect each requested outward U with programme bounds, explicit strict
   permissions and its neighbouring separators. H is never shrunk or enlarged
   by this ownership rule. The requested/effective margins and clipping reason
   are shown independently. A short or zero gap may provide no useful margin.
4. Do not reassign ownership after skipping a root or creating descendants. A
   correction may not borrow an unused neighbour's area, reset split allowances,
   stack with another root, or create a joined-root gain assignment.
5. Then apply exact pins, local protections and full changed-support/separation
   validation. Adjacent permissions do not guarantee a valid pair of changes;
   the unchanged separation after ramps still applies. Reject a conflicting
   combination or consider another already-nominated legal candidate, without
   silently sliding an accepted edge or flattening an existing gain transition.

The selected-section/selection mode still has one explicit H/U; its established
strict-area and margin semantics remain. Refinement-off keeps original edge pins.
Protection with holes remains a prohibition checked on actual support, not permission
to bridge the hole. If one inferred event spans two roots, report the limitation:
coordinated cross-root processing is outside this amendment.

## Enumeration and ranking proposal

Keep the existing bounded totals. Before construction reserve distinct nominations
for the earliest supported entrance, complete ending and in-decay ending; use at
most four of the twelve nonmandatory location slots for the outer alternatives.
Reserve two locations for an eligible spared-middle pair. Fill unused slots with
the existing need-change nominations, ordered by evidence tier, strength, source
frame and stable root ID. Mandatory exact pins retain their existing treatment;
they do not become detector suggestions. Dedupe exact locations and canonical audio.

Within 16 geometries per root, reserve the legal original, complete-event, tighter
and spared-middle alternatives before generic combinations. Failed legality does
not spend a fictitious measured pass; report why an alternative could not be formed.
Use round-robin root/geometry visits so a long rhythmic root does not consume the
whole search. All gain choices, complete root combinations and global variants
still share 64 proxy evaluations and beam width four. Full-event coverage does not
authorize extra attenuation, a new root, or a fourth full measurement.

For an explicitly comparative request, reserve an original-geometry complete-plan
challenger and an event-based complete-plan challenger for actual measurement,
subject to legality/availability. The third slot is the best distinct tighter or
spared-middle challenger, or a feedback-corrected complete plan needed after an
actual failed rule. Record the scheduling choice and omitted comparison. Comparing
three boundary families does not also guarantee three normalisation retries. A
previously saved passing candidate is not a fresh measurement of changed audio.

**Recommended passing-plan comparison, for Strategy review:**

- Only actual all-selected-rule passes enter final passing selection. Unavailable
  evidence or event confidence cannot convert a measured failure to a pass.
- Freeze a set of relevant event groups before evaluating any plan. Assign each
  plan an intervention signature of groups whose need-bearing material it changes.
  A plan cannot redefine the evidence population to hide its own clipped entrance.
- Compare event treatment **within the same intervention signature**. If the shared
  evidence is complete and corroborated, prefer: (a) no added transition cutting a
  corroborated attack; (b) no return before the last corroborated attack completes;
  (c) fewer new internal transitions unsupported by sustained correction-need
  change; then the existing preservation order. Use Boolean conflict flags for
  (a)/(b), not the number of beats or events covered. Count only new audible
  transitions, not annotation dividers or pre-existing untouched transitions.
- If evidence is tentative, contradictory or incomplete, do not award a zero-risk
  score. Use preservation ordering within that signature, retain the alternate for
  review, and label the event comparison inconclusive. During-decay versus
  after-decay returns are also incomparable on artistic grounds without listening.
- Select one representative from each signature, then use existing preservation
  order across signatures. Expose distinct representatives as trade-offs; do not
  claim that covering more events is better. This constrained use of event
  preference avoids favouring unrelated processing merely to improve coverage.
- Evaluate the complete set using deterministic stages and source-frame/plan-ID
  tie-breaks, never a pairwise comparator with non-transitive tolerance bands.
  Preserve **exactly** the existing anchored 0.01 LU failed-candidate shortlist
  and its remaining-rule ordering. The event preference applies only to passing
  candidates and heuristic exploration, not to widening that failed shortlist.

This is a deliberate qualification of the amendment's broad phrase "event coverage
before duration." A single total artistic-quality score would be unjustified.
The proposed partitioned comparison can improve two treatments of the same entrance
while being honest about different sets of corrected material. Strategy should
approve this rule or specify another explicit comparison before implementation.

## Data, memory and one-pass feasibility

Proposed immutable records, with bounded capacities:

- `EventFeatureDefinition`: filter/window/hop/threshold/grouping definitions,
  role policy, numeric mode and version; no new measurement-method label.
- `EventEvidenceBinding`: verified original programme content/assembly/composition,
  ordered roles/range, frozen Current snapshot/plan, feature definition and context
  definition. Correction associations additionally bind requirements and H/U.
- `EventGroupEvidence`: root/group IDs, onset/last-attack/decay brackets, feature
  values and channel witnesses, supporting measured window IDs, confidence/reasons,
  clipped/missing-context flags and nominated alternatives.
- `BoundaryTreatment`: requested intent, exact entry/full-gain/return/end frames,
  before/after ramp differences, permission conflicts, considered/omitted alternatives
  and evidence version. Keep private prediction/measurement authority separate.

Do not allocate 5/10 ms features for an entire three-hour programme. At 10 ms,
1,080,000 records of even 64 bytes alone would exceed 64 MiB before existing analysis.
Instead compute continuously and retain bounded summaries plus detailed feature
neighbourhoods. Proposed maximum detail: two original edges plus twelve nomination
neighbourhoods per root, each up to six seconds, for eight roots. Their union is at
most 672 s; at 5 ms and a 64-byte compact aggregate record this is approximately
8.21 MiB. Record coordinates implicitly in immutable segment headers; do not keep
per-candidate copies. Deduplicated overlapping contexts only reduce that bound.

Use a **10 MiB additional event-feature sub-budget within the existing total 64 MiB**,
including detail, rings, capacity, summaries and witnesses. Retain at most 64 group
summaries per root (512 total); coalesce successive grouped attacks into first/last
and sustained-activity summaries instead of retaining every hit. Cap detail and
summary selection deterministically by association with frozen need anchors, then
evidence tier and source frame. Any omitted required group evidence marks that
root's event comparison incomplete; it cannot earn corroborated coverage. Do not
silently lower cadence, omit a critical first/last event, or expand the memory limit.

Keep filters/running windows per channel and reusable small PCM scratch. Full ordered
role witnesses can be retained for nominated events; full per-channel high-rate
arrays for the whole programme are unnecessary. Capture the designated original
edge/coarse-nomination neighbourhoods during the existing pass; newly discovered
events use streaming rings and bounded summaries. If an event needs unavailable
retrospective detail, mark it tentative rather than launching another PCM pass.

The existing three-hour-sized Logic receipt reports 35,846,608 bytes for analysis
plus prediction sidecar, **before all model/search scratch**. The service separately
accounts for transformed windows, prefixes, ranked plans and text. Preflight the
combined allocation, including the new sub-budget, and fail clearly if it cannot
fit. The earlier whole-process footprint/RSS is not proof of additional-storage
compliance. Verify actual reserved capacity/high-water use during future tests.

Expected added work is O(samples × channels) in the existing auxiliary pass, with
small filters and streaming sums; it adds no decoder pass or retained PCM file.
No speed multiplier is claimed before measurement. Preserve cached planning targets
<=10/20 s (Sequoia/High Sierra), 30/60 s watchdog, <=1 s cancellation and native
p95 <=100 ms responsiveness. Report auxiliary extraction time separately from
cached search and complete candidate measurement.

## Identity, sessions and UI explanations

Propose settings schema 2 and explicit new IDs for permission, event formation,
search and passing ranking: `region-constraints-v2`, `event-boundary-features-v1`,
`region-formation-v2`, `region-gain-search-v2`, `candidate-review-rank-v3`.
Retain `window-energy-proxy-v1` only while its arithmetic is unchanged. Retain the
measurement method, profile definition and `section-gain-v1` unchanged. Fixed
threshold changes participate in the feature definition/digest, not an invisible
build-only heuristic.

Read old region settings as their original v1 method. Reproduce their identities
and inward-only semantics; do not auto-enable event analysis, margins or 60%.
The duration value 60% already fits the existing versioned proposal preferences;
an experimental named configuration must embed its full definition. Keep the
conservative default at 20% until a later explicit product decision.

Recommend session document schema 10 for sessions using the new settings/history,
while retaining the existing programme representation version (the media/assembly
model has not changed). Separate those version counters rather than bumping both
because they currently happen to be 9. Historical explanation attachment v3 adds
feature/evidence/ranking provenance and compact boundary reasons; keep v1/v2 reading.
Validate mandatory new settings, and retain replaceable measurement caches. Saving
the new format requires a new reader; explicit Save/Save As preserves the old original.
Old readers must reject unknown mandatory permissions safely, not ignore them.

Do not persist fine feature arrays, private candidate measurements or verdicts as
Current verification. Persist the exact selected processing plan plus compact
historical evidence and reasons under the existing 8 MiB explanation limit. Source
or processing edits invalidate features/proposals appropriately; changing requirements,
duration or margins may reuse valid measurements but stales affected plans and export
approval. Independent runtime job identity rejects late completions after New/Open.
Video offset, routing and monitor volume remain outside programme audio identity.

UI fields should answer, for every applied/proposed region:

1. Which measured windows established need, their gate status and planning threshold.
2. Which event features nominated the entrance/exit, which channels supported them,
   their actual window supports, confidence and unavailable/omitted evidence.
3. H versus requested/effective U, ownership clipping, exact pins/protection, and
   why a wanted boundary was rejected or never examined.
4. Ramp start, full-gain point, return start/end and internal transitions; show
   event brackets and nominal cuts separately. Do not display 10 ms cadence as
   10 ms proven onset accuracy or label events as beats/dialogue automatically.
5. Requested/achieved local attenuation, reserve/uncovered estimate and separate
   global gain; 20/40/60% selected ceiling, actual newly changed union, cumulative
   support, split/count/capacity and search budgets independently.
6. Actual versus predicted outcomes, event-comparison conflicts/inconclusiveness,
   retained alternatives and reason for stopping. Historical reasons remain
   historical after subsequent editing; current verification remains independent.

## Impact-selected qualification proposal after approval

No tests in this table were run for this review. First implement/stabilize on
Sequoia, then freeze one combined candidate for affected actual High Sierra and
cross-platform checks. Build only named standalone/test targets. Preserve the
existing M20/UI candidates, session pairs and evidence unchanged.

| Area | Fresh checks and acceptance |
| --- | --- |
| Detector and channel/domain rules | Known synthetic impulses/bursts and gradual changes; DC, silence, quiet content, noise, low-frequency attacks, antiphase stereo, reordered 5.1, isolated LFE, finite over-full-scale sums. Compare 5/10/20 ms at 44.1/48/96/192 kHz; verify rational cadence, chunk invariance, explicit uncertainty and no drift. No measured loudness changes from observation. |
| Existing processing | Source versus frozen Current; gain changes without source events, old ramps, unity, compositions and separate-file callers. Do not generate event evidence solely from the prior gain automation. Exact independent offline/audition sample parity for unchanged plans. |
| First/last music entrance | Include the user's ~12.283 s landmark and establish the actual last-hit bracket through new evidence/listening. With confirmed permission, enumerate a ramp completing before the corroborated first attack and alternatives covering the last. With permission off, explain the conflict without widening U. No hard-coded times or automatic claim of musical correctness. |
| 52.615–57.762 s investigation | Compare original, coherent complete passage and tighter/spared-middle plans; inspect the ~55 s dip. Preserve the possibility that the landmarks need adjustment. No individual-peak chopping as the only available challenger. Measure whole plans; numerical pass alone does not pass listening. |
| Sustained and dialogue/ambience | Gradual rise/decay, background noise, pauses, no clear onset, conflicting channels and open endings. Uncertain detector results preserve a legal original challenger and do not become artistic-quality claims. |
| Simple success / heterogeneous / 5.1 | Retain prior useful correction and a spared-middle challenger. Exercise margins, exact pins, old/new ramp accounting, adjacent roots, protection holes, root split limits and unavailable confidence. Include the existing 5.1 Busy material with honest TP failures if they remain. |
| Gate crossing / cinema | Retain full gate recomputation, anchored failed ranking and measured-versus-predicted failures. Bounded existing eight-minute cinema excerpt only; 15 LU success is not required. No repeat of the full-film workflow. |
| Search/ranking | Permutation-independent whole-set comparison, equal-cost tie-breaks, signatures with different corrected groups, tentative/corroborated/capped evidence, neutral during/after-decay comparison, exact budget exhaustion and no fourth measurement. No pass from confidence or prediction. |
| Duration control | From identical frozen baseline and permissions, compare 20/40/60% with the same gains/search budget. Include overlap union, old ramps, metadata-only splits, global gain independence and a case where a larger cap does not improve the result. |
| Session/lifecycle/UI | Old sessions preserve v1 semantics; explicit new selection; settings/history round trip, older-reader rejection, recoverable caches, coherent saves, stale jobs/cancel/source replacement; silence/device isolation and visible uncertainty. Preserve recent manual timeline editing and connection preference behavior. |
| Resource/platform | Worst-case three-hour-sized feature/population model without three-hour PCM decode; all eight broad roots, rapid many-group input, maximum pins/detail at 5 ms. Enforce 64 MiB total model/search budget and feature sub-budget, cancellation, cached targets and native responsiveness. Run selected complete numerical cases on actual High Sierra; cross-build alone is insufficient. |

Use a staged comparison, not a Cartesian sweep: (A) detector cadence on short
fixtures with all other constants fixed; (B) v1 versus new boundaries at 20% with
the same zero outward permission, followed by a separate new-method zero-versus-1 s
margin contrast (v1 cannot accept outward margins); (C) 20/40/60% with the chosen
boundary method, fixed permissions and same baseline;
(D) freeze definition/settings and evaluate new independent material/listening.
Each comparative request retains its own three-assessment ceiling. Predeclare the
finite run list, and preserve policy-limited failures instead of retrying until pass.

Reuse the unchanged qualified measurement arithmetic/official-reference, processor,
video and long-media evidence only with source/dependency equivalence. Fresh integration
tests must show that the modified auxiliary caller has no side effects on these paths.
Repeat official references if authoritative meter code/arithmetic changes. Repeat
selected video/long execution only for concrete changed timing/buffering/dependencies
or an observed failure. No automatic codec matrix, endurance or full cinema rerun.
Coordinate any physical-device opening/listening separately. Occasional picture skips
remain accepted under the existing continuity, drift and lifecycle conditions.

Provisional resource envelope, to refine after approved implementation: 20–45 minutes
for a named universal build when caches are available; 10–25 minutes short development
checks; approximately 30–60 minutes selected frozen offline checks per Mac excluding
native/manual/listening time. These are scheduling estimates, not measured timings.
Reserve 5–10 GiB for the new build/package and up to 2 GiB temporary short fixtures;
reuse existing immutable media and avoid new three-hour PCM. Do not discard any frozen
candidate to meet these estimates.

## Decisions to resolve before implementation

Strategy/user review should settle these substantive choices:

- Accept the energy-based event heuristic and its stated limitations; the numeric
  detector hypotheses need controlled evaluation rather than presumed approval.
- Accept fixed midpoint ownership for automatic outward margins, with conflicts
  disclosed instead of cross-root coordination; confirm 1 s per edge as an
  explicit experimental setting, not an automatic default.
- Accept comparison within matching intervention signatures and neutral handling
  of uncertain/decay evidence, or specify a different deterministic passing rule.
- Confirm the 60% named experimental configuration while 20% remains unchanged,
  and reserve fresh holdout material after using the current comparisons to develop
  the method. No implementation should be tuned only to make this programme pass.
- Accept the proposed feature sub-budget, evidence-completeness fallback and new
  method/settings/session/history identities. No budget or supported layout expansion.

**Review outcome:** conceptually suitable with the explicit qualifications above;
not yet an implementation authorization or a validated detector design. No production
files/defaults/sessions/media/packages were changed, and no new audio experiments,
builds, tests, hardware actions, messages to UI, commits or publication were performed.

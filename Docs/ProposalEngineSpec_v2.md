# LevelScope — Proposal Engine Specification v2

Status: Design clarified; bounded M13 candidate-02 implemented and qualified, 28 October 2026.
[ProjectIndex 49](ProjectIndex/49_M13QualificationAndHandover.md) records actual-platform
evidence and limitations. The user subsequently authorized checkpoint publication and main integration.

Scope: Standalone-first; shared Core contracts; plugin parity remains separate.
Depends on: `Docs/ArchitectureSnapshot.md`, `Docs/ComplianceMatrix.md`, `Docs/ProjectIndex/*`

This revision reconciles the earlier broad design with the implemented M1–M12
contracts and the agreed three-strategy direction. The user subsequently approved
bounded M13 implementation and qualification, then authorized publication and main
integration after readiness review. The combined [M13 checkpoint](https://github.com/carlozaiana/LevelScope/tree/m13-bounded-section-proposals-2026-10-28)
fast-forwards main from preserved M11 `cfb1033aa7fb65a2655db276bee696161196bdeb`. Published M12 candidate-03 remains
**“Sequoia tested; High Sierra and video qualification pending.”**
[ProjectIndex 47](ProjectIndex/47_M13BoundedSectionTrimProposals.md) defines the first
bounded section-led method; future capabilities below are not M13 scope.

---

## 1. Purpose

Generate explainable processing plans for a selected, authoritative profile while
preserving the programme's artistic intent. An honest failure within a bounded
policy is a valid method outcome. It does not prove that another strategy cannot
succeed, and must not be disguised as compliance or solved by relaxing a profile.

For the broader architecture, no change is appropriate when **all applicable,
explicitly selected objectives are satisfied**. Passing basic programme rules must
not suppress an explicitly requested validated dialogue or other objective. Every
additional objective needs an identified method, applicability and acceptance
criteria; preferences never invent objectives or new compliance rules. The objective
set and its versions belong in proposal identity and comparison-baseline bindings.
M13 supports only the selected programme rules, so its existing already-pass no-op
remains unchanged. Consider sufficient global-only correction before other processing. Beyond that common prefix, intervention preferences
are strategy-specific. Every added processor needs a reason and qualified behavior;
being available never means it must be active. Apply and verification stay explicit.

---

## 2. Product philosophy and terminology

Separate four independent concepts:

| Concept | Meaning |
| --- | --- |
| Delivery profile | Versioned measurement method, applicability and rules. Defines what a measured pass means. |
| Processing strategy | Permitted interventions: dynamics-led, section-led or mixed. |
| Selection algorithm / method | How settings are chosen within a strategy; multiple algorithms may later serve each strategy. |
| Processing preferences / policy | Versioned preservation priorities, activity budgets and search/resource limits, independent of delivery rules. |

The older `HybridMethod` combined rule-based diagnosis, distribution analysis,
proxy simulation and bounded refinement. That is **combined-analysis selection**,
not the mixed processing strategy. Use `CombinedAnalysisMethod` in future design
names; do not rename existing source files or migrate code as part of this revision.
It is neither a mandated production default nor a prerequisite for M13.

Reports identify strategy, algorithm and policy separately. They explain the
constraint, evidence, exact change, expected and measured effects, processing
activity, limitations and next available action. Artistic preservation and any
future dialogue objective require listening/method-specific validation beyond
scalar programme loudness and peak checks.

---

## 3. Authoritative compliance model

Use the same versioned profile definition and evaluator for every strategy. The
implemented EBU R128-2023 whole-programme linear-audio QC profile remains unchanged:
intended programme target **-23 LUFS**, explicit **±0.2 LU** allowance, strict
**-1 dBTP** ceiling, and informative LRA with **no universal LRA ceiling**. This is
the existing M2 contract, not a new interpretation of the standard; see
[21](ProjectIndex/21_M2MeasuredCompliance.md) and
[14](ProjectIndex/14_OfflineMeasurementContract.md).

The separate documented meter-tolerance provision never becomes an evaluation
allowance. Compare full precision; display rounding, engine test tolerances,
planning reserves and user preferences cannot widen accepted limits. Remove the
older illustrative ±0.5 LU default and “near-hard” loudness interpretation.
Every required profile rule must pass; no preservation mode turns a failure into
an acceptable delivery result. Invalid profiles, unsupported applicability, stale
snapshots, missing required metrics or predictions cannot yield an overall pass.
Unavailable LRA does not block this profile when it has no LRA rule.

Use narrow domain-specific wording: “Passes selected programme loudness and
true-peak rules.” Do not imply complete broadcaster delivery, metadata correctness
or future encoded-export verification. Live exceptions, short-form supplements,
deliberate lower targets and distribution-specific requirements remain outside the
current profile. ATSC/dialogue-gated methods and intelligibility targets require
separately defined methods, applicability and validation; they are not interchangeable
with level-gated whole-mix loudness.

Preferences may govern permitted alteration, costs and future optional objectives.
They cannot edit a delivery profile or silently create an LRA/intelligibility rule.

---

## 4. Truth domains, application and export

| Domain | Authority and required display |
| --- | --- |
| Source | Actual measurement of the original confirmed programme, including any original composition. |
| Target | Selected profile rules and applicability; not a measurement. |
| Prediction | Estimate for a proposed plan, labelled predicted; never a measured verdict. |
| Private measured candidate | Real worker measurement bound to that exact candidate plan, independently labelled; not installed as Current or an export approval. |
| Current State | Actual measurement of the explicitly applied plan. Audio edits immediately make it UNVERIFIED until independently remeasured. |
| Verified exported audio | Measurement of the finalized, reopened artifact, with its own content hash and format identity, no extra processing gain. |

Generate/review/cancel never changes Current. Explicit stopped-only Apply installs
one immutable plan atomically, invalidates affected Current truth and export
approvals, and retains Source. Explicit independent Current remeasurement remains
required even after a privately measured candidate passes. A prediction, private
candidate or previous passing Current result cannot fill that gap.

Verification disagreement displays actual values and verdict and stops. Cancelled
or unsuccessful measurement leaves applied settings unverified. Completed failure
remains a measured failure. Retain settings; no implicit rollback, retry, correction
or automatic strategy fallback. Already-pass and identical-plan outcomes are true
no-ops: no revision increment, invalidation or worker request.

Render and verify / Publish export remain separate M10 actions. Both fresh Current
and actual artifact must pass the selected rules before local audio publication;
Source need not pass. Apply deterministic dither once at final PCM24 quantization,
not in Current processing. Preserve separate measurement-only reports, unavailable
statuses and missing artifact identity; do not invent artifact hashes or verdicts.
The earlier generic `ExportReady` from Current alone is insufficient.

---

## 5. Gate Visibility

The app should expose gate values numerically and visually.

Recommended visible gates:

- Integrated Loudness absolute gate
- Integrated Loudness relative gate
- LRA effective lower inclusion threshold / gate

## 5.1 UI representation
On the loudness timeline:
- show gates as horizontal lines

On histogram/distribution views:
- show gates as vertical markers on the loudness axis

Purpose:
- explain what content is counted
- explain why some sections contribute to IL and LRA and others do not
- increase trust and reduce trial-and-error

---

## 6. Boundaries, time scales and optional markers

| Entity | Meaning and provenance |
| --- | --- |
| Editorial scene/edit | A user/editor-supplied narrative or cut annotation. An audio algorithm cannot certify it as a cinematic scene. |
| Processing section | Exact source-frame interval used by an explicitly applied gain plan. It can differ from editorial annotations. |
| Loudness region | Method-derived interval of similar measured behavior, with method/version/confidence. It is descriptive, not a verified scene or an applied change. |
| Peak event | Future grouped peak excursions with an explicitly defined threshold, duration, merge-gap, channel and interpolation policy. A maximum is not an event count. |

Marker import is optional for the eventual application. Dynamics-led operation
needs no section map. Section-led methods may use manual sections, imported marker
annotations or future audio-derived regions. M13 uses the existing manual editor;
no marker/audio-region importer is required. It may explain that a whole-programme
section offers insufficient localization instead of inventing cinematic scenes.

Keep source coordinates, half-open intervals, origin, author/importer or analysis
version and user confirmation visible. Marker/EDL timecode requires explicit
programme origin and rational frame-rate mapping, exact audio-frame conversion,
preview and confirmation. An imported boundary never automatically moves an
applied transition. User protection instructions take precedence over automatic
suggestions; unsupported mandatory protection must reject a changing proposal
rather than be silently ignored. A genuine no-change result needs no intervention
and remains a no-op.

### Measurement gates, loudness regions and DSP thresholds

Measurement gates select populations for authoritative metrics and must follow
the measurement method. They are not edit boundaries or compressor thresholds.
Loudness-region boundaries come from a separate versioned segmentation algorithm
with detector/window, hysteresis, minimum duration and merge rules. Selecting such
regions for processing requires an explicit mapping/confirmation policy.

A DSP threshold must identify detector, domain/units, reference point, window,
channel linking/weighting, gain position and mapping method. A LUFS gate, a region's
median or a histogram percentile cannot be copied into a dBFS threshold as though
they meant the same thing. A future loudness-controlled compressor may consume a
compatible LUFS detector only through a defined, tested transfer rule. These
mappings require qualification; M13 adds none of them.

---

## 7. Three processing strategies and intervention preferences

Common prefix: validate measured baseline/profile/explicit objective set/protection;
preserve Current when every applicable selected objective is satisfied; consider sufficient global-only correction with upstream processing
fixed. A global-only prediction still follows explicit Apply and measured verification.

| Strategy | Permitted future processing | Preferences after the common prefix |
| --- | --- | --- |
| Dynamics-led | Global gain, downward and upward compression, final true-peak limiting | Diagnose sustained downward/upward needs and residual peak behavior. Prefer the sufficient qualified dynamics stage(s) within activity budgets; use limiting only when justified. No sections or markers required. |
| Section-led | Global gain, section gains and transitions, final true-peak limiting | Prefer reviewable local corrections within section/transition/protection budgets. Consider the cost to satisfactory neighbours. Future limiting may handle residual events when justified; no compression in this strategy. |
| Mixed | A justified combination of section adjustments, dynamics, global gain and final limiting | Demonstrate why the combination is warranted against the same baseline and permitted simpler plans. Account for interacting stages and cumulative alteration; do not enable both merely because they are available. |

Intervention choice order is distinct from sample-processing stage order. Every
plan declares exact active stage order and versions. M12's qualified order remains
composition → section gain → global trim; M13 reuses it unchanged. Future qualified
plans put final peak limiting after all gain-changing programme stages, before
export quantization and outside monitor controls. New latency/stateful stages need
explicit timing, seek/preroll, audition/offline/export parity qualification.

M13 is the first **bounded section-led algorithm**, with global/section gains only.
Its declared capability subset has no qualified limiter or dynamics. It can report
no solution; it cannot invoke unqualified DSP or silently change strategy. Full
section-led limiting and the other strategies arrive through separate milestones.

---

## 8. One shared foundation

### 8.1 Shared services, separate selection algorithms

All strategies use the same authoritative profile catalog/evaluator, verified
programme inputs, measurement kernel, descriptive analysis, processor implementation,
identity and invalidation rules, job lifecycle, candidate assessment, Current
verification and artifact export gates. Alternative algorithms consume these
contracts; they do not fork metering, DSP or compliance logic.

Responsibilities remain separable:

- **Proposal orchestration:** selected strategy, algorithm/capability and frozen
  preference snapshot; bounded jobs and results; no automatic Apply.
- **Measurement/analysis/attribution:** immutable measured snapshots; profile-independent
  descriptive data; profile/policy-dependent annotations bound separately. Peaks and
  energy shares do not by themselves establish the result of an edit.
- **Policy:** permitted interventions, hard user limits/protection, soft preservation
  costs and compute budgets. It never changes delivery thresholds.
- **Prediction and bounded refinement:** optional algorithm-owned private planning,
  explicitly versioned and resource-limited. M13 has one envelope and no search loop.
- **Processing:** shared exact immutable plan executed by offline measurement,
  audition and export. No plugin-state shortcuts or hidden limiting/attenuation.
- **Verification/export:** private candidate assessment, independent post-Apply Current
  measurement, then finalized-artifact measurement as distinct authority levels.
- **Explanation/activity:** reasons, exact settings, constraints, measured activity
  and limitations. Future protection is a shared constraint, not a special algorithm.

Existing DSP source in Core/the plugin is implementation material, not standalone
qualification. No compliance claim follows until its shared standalone path,
state/reset/seek behavior, channel policy, resource bounds and actual-platform
sample parity are qualified. Do not install or modify the plugin for this work.

### 8.2 Immutable contracts and identity

A request binds source and Current snapshot identities, programme/member content,
ordered roles, alignment/range, current plan and revision, method/build/platform,
section geometry/provenance if relevant, selected profile digest/applicability,
strategy ID/version, algorithm ID/version, policy ID/version and preference digest.
Runtime session/job identity is separate from persistent document identity.

A reviewable plan/report contains exact ordered active stages and full-precision
settings, unchanged stages, baseline/additional/cumulative corrections, boundary
provenance, reasons, prediction limitations, private candidate measurements with
domain labels, processing activity, budgets consumed, failed constraints and next
action. Inactive stages stay inactive. Unknown/unsupported versions fail closed.

Processing changes invalidate affected Current measurements, diagnostics, proposals
and export approvals. Programme/member/role/alignment changes additionally invalidate
Source as applicable. Compatible profile changes reuse audio measurements but rebuild
evaluations and dependent proposals. Preference/algorithm changes invalidate proposals,
not unchanged audio measurements. Monitor routing/volume and video attachment/offset
stay outside audio identity. Annotation edits affect dependent analysis/proposals,
never silently modify applied audio. Late jobs cannot enter another runtime session.

Original components remain separately accessible; Source is the original confirmed
composition. Future component DSP participates in Current identities before composition;
whole-programme processing follows the sum. Analysis-only references are independently
identified and non-audible: changing one invalidates its dependent dialogue analysis,
not ordinary measurements of unchanged programme audio.

### 8.3 Strategy comparisons and resources

Generate only the selected strategy by default. Reuse compatible verified input,
measurements and analysis, checking identity, geometry, method and provenance. A
private candidate is not interchangeable with baseline analysis.

Comparisons are explicit and begin from the **same frozen measured baseline**, not
from each other's applied output. Record the same profile and common preferences;
disclose strategy-specific differences. Before a comparison, bound selected methods,
candidate count, measurement passes, time, memory and disk; cancellation removes
owned temporaries and leaves Current intact. Do not fan out every strategy or hide
unbounded retries. M13 implements no multi-strategy comparison or selection UI.

---

## 9. Workflow states and invalid transitions

Use explicit distinctions rather than a single generic “Verified” state:

`MeasuredBaseline → ValidatedRequest → Draft/Prediction → PrivateCandidateMeasured
→ UserReview → ExplicitApply → CurrentUnverified → ExplicitCurrentRemeasurement
→ CurrentMeasuredPass/Fail → RenderAndVerify → ArtifactMeasuredPass/Fail
→ ExplicitPublish`

Global-only M3 retains its existing prediction/Apply/remeasure path; the new M13
section path requires private measured candidate success before offering Apply.
No-change exits before new jobs. Unavailable, cancelled, stale and no-solution are
terminal outcomes of that request, not instructions to retry. Failed post-Apply
verification stops with the applied plan retained. Publication requires both fresh
Current and the exact finalized artifact to pass; each failed gate stays visible.

Future protection edits alter proposal constraints, not audio by themselves.
Changing the audible plan is always explicit. Sessions retain applied settings and
validated replaceable caches; saved predictions/verdict strings are not authority.
Reopening remains stopped, device-closed and routes-muted.

---

## 10. Strategies, algorithms and preference presets

Strategy IDs describe allowed interventions (`dynamics-led`, `section-led`, `mixed`).
Algorithm IDs identify selection behavior within them; an algorithm declares supported
strategies and processor capabilities. Version them separately from processing-engine
and profile versions. M13's proposed algorithm is `section-trim-proposal-v1` within
a bounded `section-led` capability subset.

Future presets such as preservation-focused or balanced may set processing budgets
and costs. They never make a required rule soft. The older `StrictCompliance` name
must not suggest other modes may publish noncompliant audio. User-defined preference
controls/preset editing are not M13 implementation scope.

Rule-based, distribution-aware, combined-analysis and later optimizer methods are
alternative **selection algorithms**, not the three processing strategies. No
algorithm, including combined-analysis, is automatically the production default.

---

## 11. Measurement, analysis and activity outputs

Reuse the authoritative integrated, momentary, short-term, LRA, sample-peak and
true-peak results with availability and LRA-stability status intact. Do not convert
silence, below-gate, insufficient-duration or missing values to measured zero.

Descriptive analysis includes full-resolution windows, gates, distributions,
rankings, exact audio-source positions and optional section summaries. Display
reduction preserves extrema/gate membership but never drives planning calculations.
Profile annotations bind separately. Region classification and peak-event collection
are future capabilities, not requirements for ordinary measurement or dynamics-led
operation. LRA distributions remain separate from the programme timeline.

For M13, activity means changed intervals/transition support, requested and applied
gains, reserve use/capping and affected duration. Section maxima identify affected
sections and locations, **not the number or density of distinct peak events**.
Future dynamics activity must distinguish measured gain reduction/boost, affected
channel/time fraction, sustained duration, distribution/percentiles and prediction.
Specify the activity detector threshold and collection method; absent activity is
unavailable, not a fabricated zero or proof of inaudibility.

---

## 12. Future dynamics threshold design (T0–T3)

T0–T3 may describe future support/preserve/control regions for a particular
qualified detector and transfer curve. They are processor settings, not measurement
gates, editorial cuts or universal loudness boundaries. They cannot identify normal
dialogue or intelligibility from level alone.

Any mapping from descriptive loudness statistics to these settings must declare
units, detector/window, weighting, channel policy, reference stage, hysteresis,
attack/release, ratio/knee, range and calibration. Section 6 defines the separation.
Upward processing must address silence/noise/ambience explicitly; downward processing
must address sustained and transient activity. Neither is authorized by M13.

A future mapping needs synthetic and real-programme tests, actual-platform
measurement/audition/export consistency and listening review. Availability of existing
plugin curves/detectors does not establish this qualification.

---

## 13. Deterministic first-pass planning

For identical measured inputs and versioned preferences, use a fixed, explainable
reasoning chain with stable tie handling:

1. Validate the common baseline, identities, profile applicability and constraints.
2. Return no-change if all applicable explicitly selected objectives are satisfied
   (programme rules only in M13); otherwise consider sufficient global-only correction.
3. Select only interventions allowed by the chosen strategy and qualified capability.
4. Build an immutable intent/plan with reasons and explicit preservation costs.
5. Optionally predict/refine within that algorithm's declared private budget.
6. Measure the private candidate through the shared worker where required, keeping
   it separate from Current. Report success, unavailable or the limiting constraint.
7. Offer review and explicit Apply; independently remeasure Current afterward.

No universal section-before-dynamics hierarchy applies after step 2. No LRA
correction is implied without an explicitly selected and validated objective.
[47](ProjectIndex/47_M13BoundedSectionTrimProposals.md) owns the exact M13 decision
sequence and three-pass bound. This broad design does not expand it.

---

## 14. Selection algorithms

A deterministic rule-based algorithm is the first bounded section-led method.
Subsequent strategies may each offer alternative algorithms with a shared request
and result contract. They must retain the same truth, identity and export gates.

Distribution-aware selection may later use descriptive statistics to guide a
specified objective; it cannot turn LRA into an unrequested compliance rule.
`CombinedAnalysisMethod` (formerly the design's `HybridMethod`) combines analytical
techniques within an explicitly selected processing strategy. It is not synonymous
with `mixed`, and neither method is implemented by M13.

Any future local optimizer or learned initializer has explicit search/resource
budgets, deterministic/reproducible provenance and independent measured verification.
Algorithm innovation cannot redefine profiles or install predictions into Current.

---

## 15. Attribution Maps

Reuse descriptive measured attribution where available, keeping profile-dependent
annotations separate. These maps do not prove the causal result of removing or
changing a passage; that requires candidate measurement.

### `ILContributionMap`
Shows included measured-block energy shares, with overlap and gate context intact.

### `LRALowerTailContributionMap`
Shows which windows drive the lower tail.

### `LRAUpperTailContributionMap`
Shows which windows drive the upper tail.

### `PeakContributionMap`
Shows measured peak locations/section maxima and separately bound ceiling annotations.
It is not a distinct-event counter.

### Later: `IntelligibilityContributionMap`
A future method-defined diagnostic; no validated intelligibility claim exists yet.
Keep dialogue-reference and final-mix evidence distinct.

This is the core of explainability and efficient correction.

---

## 16. Failure explanations and future limiter decisions

“No solution within this policy” must identify the strategy/method, baseline and
profile, limiting constraint, units and configured bound, observed/measured values,
required versus allowed change and affected intervals/channels. Identify whether a
candidate was measured, which rule failed, what remains unknown and whether a
reserve was capped. Distinguish unavailable data, unsupported processing, budget
exhaustion, protected-passage conflict and actual measured failure. Do not imply
that no other strategy could succeed; do not promise that an untested one will.

For M13, report section counts, durations and maxima only. Distinct peak-event
counts require a separately versioned collector/grouping definition: threshold,
minimum duration, gap merging, cross-channel grouping and interpolation/EOF handling.
One channel maximum can hide many events; six channel maxima can belong to one event.

Future limiter decisions consider **excess magnitude, duration and distribution**,
including sustained activity, temporal/channel concentration and the cost of
attenuating entire passages or satisfactory neighbours. Count alone cannot choose
between limiting and broad level changes. Predict and then measure limiter activity
(maximum/percentile gain reduction, affected duration and consecutive activity),
verify true peak afterward and review audibility. Limiting is optional, not a
silent safety net, and must be qualified in the standalone path before any claim.

Future dialogue-stem loudness, speech-gated final-mix measurement, dialogue/background
diagnostics and intelligibility scores remain distinct. Each requires its own method,
validity conditions and acceptance criteria; a scalar programme pass proves none.

---

## 17. Configurable processing-preference design

Define a versioned, immutable `ProcessingPreferences` design record, independent of
`TargetProfile`. Include schema/policy version, units, explicit values versus disabled
or unsupported fields, constraint severity, scope/reference baseline and a digest.
Selection algorithms consume the same structure within supported capabilities.
Unknown mandatory limits must fail closed. Controls, persistence migration and
coordinated neighbour adjustment are future implementation work, not M13 scope.

| Preference family | Definition | Initial M13 effective default |
| --- | --- | --- |
| Additional section correction | Signed delta from frozen measured Current; separate upward/downward magnitudes | Upward 0 dB; downward at most 3 dB per section, reserve included |
| Cumulative section correction | Resulting section envelope relative to original confirmed Source, distinct from the new delta and global gain | Existing M12 plateau range -12 to +12 dB; retain prior manual settings, show baseline/new/resulting and whole-chain effects |
| Neighbour differences | Absolute difference between adjacent resulting plateau gains; optionally also a separate limit on adjacent new deltas | Report the difference; no new independent hard neighbour cap or coordinated adjustment in M13 |
| Transitions / slew | Versioned duration/shape, permitted shortening and optionally maximum dB/s or jump; constraints must be mutually consistent | Existing `section-gain-v1` nominal 100 ms, exact shortening/rounding; report derived slew; no new rate cap or transition control |
| Affected extent | Count and union of actual changed intervals plus transitions, with per-request and cumulative extents distinguished | At most eight changed sections and 20% programme duration per request; show cumulative envelope separately |
| Planning reserve | Additional requested local attenuation, never a compliance tolerance | 0.5 dB included inside the 3 dB cap; disclose requested/achieved reserve and uncovered deficit |
| Compression activity | Per-stage max boost/reduction, percentile/mean activity, affected time/channel fraction and maximum sustained duration, with defined detector threshold | Upward/downward processors disabled; permitted added activity zero until separately qualified and approved |
| Limiter activity | Maximum/percentile reduction, affected duration/fraction, consecutive activity and event distribution with defined measurement | Disabled; permitted added activity zero in M13 |
| Protected passages | Explicit intervals, stage/scope and strength, transition/filter support, author and provenance | No new protection editor in M13; never override a supplied mandatory instruction or pretend an unsupported protection is enforced |
| Private work | Candidate count, full measurement passes, wall-time, RSS, disk and cancellation bounds | One envelope, three measurement/attribution passes; ProjectIndex 47 budgets |

“Additional” is not reset into “cumulative” by repeatedly generating proposals. Each
explicit new request may use a newly measured baseline but must still show and
validate the resulting cumulative settings. Global gain and future dynamic gain
activity are reported separately and in total; a section-only bound must not be
advertised as a whole-chain attenuation guarantee. For time-varying future DSP,
activity is measured, not inferred from static plateau sums. M13's 20% limit is on
the changed **local section envelope**, including transitions; it excludes uniform
global trim. Report global trim's whole-programme extent separately so the local
budget is never presented as a 20% limit on all altered samples.

Example: applied corrections **0, -2.5, -5, -2.5, 0 dB** have 5 dB maximum cumulative
section attenuation and 2.5 dB maximum neighbouring difference. They satisfy a
2.5 dB neighbour limit but fail a 3 dB cumulative attenuation limit. A new -2.5 dB
change on an existing -2.5 dB section is a -5 dB result; M13's 3 dB *additional* cap
must not be reported as a 3 dB cumulative cap. These are gain-envelope constraints,
not instructions to equalize neighbouring scenes' measured loudness.

For an M13 peak deficit D > 0 at the intended global correction, requested reduction
is D + 0.5 dB; actual reduction A is capped at 3 dB (and must meet the existing
plan validity checks). Achieved reserve is min(0.5, max(0, A - D)); uncovered deficit
is max(0, D - A). For D=2.8, A=3.0: only 0.2 dB reserve remains. For D=3.2, A=3.0:
reserve is zero and deficit remains 0.2 dB. Capping is visible, never extra hidden
attenuation; neither formula guarantees a pass after gates/transitions/global gain.
Exact candidate assessment still decides.

Future coordinated neighbouring adjustments must price the collateral change to
otherwise satisfactory material: count newly affected duration, gain magnitude,
slew and protected content. Smoother gain steps are not free if they attenuate
neighbours that did not need correction. Do not normalize their loudness or move
boundaries silently. No coordinated adjustment is implemented in M13; retain its
fixed ramps and disclose any transition-region effect on neighbours.

Hard safety/user limits cannot be traded against cost. Soft preservation weights
rank feasible alternatives without relaxing a profile. Future preference editing
must distinguish them and disclose unsupported combinations. For example, an
immutable boundary, a gain-step limit and a slew/duration bound may conflict: report
that conflict instead of extending a transition through protected material.

---

## 18. Future protection and region provenance

Manual sections and optional marker annotations are usable without automatic scene
recognition. Audio-derived loudness regions remain labelled as such, with method
and confidence; they are not verified cinematic scenes. Editorial, processing,
loudness-region and peak-event boundaries retain separate identities and provenance.

User protection instructions take precedence over automatic suggestions. A future
protection model must distinguish no *additional local processing*, fixed processing
settings and exact audible-output preservation. The latter also constrains global
gain, linked dynamics and transition/filter/lookahead support. Defaults cannot
silently permit a global change through a fully protected passage. Conflicting
profile/protection constraints produce an explanation, not an override.

Suggested future intents include preserving dramatic impact, ambience, dialogue
perspective and music phrasing; these labels are not measurable compliance rules.
Unknown or unavailable protection semantics must be disclosed and generation
blocked where mandatory. Protection controls and coordinated smoothing remain
separately reviewable work, not M13 implementation.

---

## 19. Explicit final correction and limiting

M12/M13 retain section gain followed by the separate global trim. M13 may derive
that global setting within the private candidate calculation, then offer the entire
plan for explicit Apply. It does not create a hidden post-verification correction.

Future final limiting, when qualified and justified, follows every programme gain
stage and precedes final export quantization. Any later gain change creates a new
plan and invalidates dependent truth/approvals. A measured miss stops and is shown;
there is no automatic micro-correction, retry, rollback or profile relaxation.
No second leveller is mandated. Existing monitor-overload protection remains
separate and must not silently change programme samples/settings.

---

## 20. Standalone-First Strategy

This workflow should be implemented first in the standalone application.

Reasons:
- whole-program analysis
- bounded private candidate assessment
- current-state re-measurement
- later protection-aware planning under separately qualified semantics
- export readiness workflow

The plugin may later reuse:
- measurement pieces
- proposal preview logic
- explanation snippets
- selected module planning components

But the full proposal/export workflow belongs first in the standalone product.

---

## 21. Edge cases and local evaluation

Acceptance must cover already-pass/no-change, global-only success, localized and
sustained peak constraints, LFE/muted-channel constraints, silence and short inputs,
large cumulative versus small additional changes, boundary/transition peaks,
collateral neighbour alteration, unsupported protection, exhausted budgets, stale
identities, private-candidate versus actual verification disagreement and cancelled
jobs. Synthetic fixtures must demonstrate real success as well as honest bounded
failure. Future dynamics tests additionally cover noise/ambience lift, linked-channel
behavior, sustained limiting, state/reset/seek and artistic listening trade-offs.

Planned local evaluation: the supplied cinema 5.1 mix and user-reported 28 October 2026 Nuendo statistics
listed in [47](ProjectIndex/47_M13BoundedSectionTrimProposals.md#planned-local-cinema-evaluation).
Preserve the original, keep raw media/screenshots out of GitHub, and run no evaluation
before its planned bounded Sequoia evaluation. The historical screenshot is Cubase 10
supplementary evidence, not the intended Nuendo comparison. External statistics are contextual observations,
not LevelScope measurement authority or proof of exact current content/settings.
The film need not pass M13: an evidenced failure inside the bounded method is valid.
Retain synthetic success cases and eventual coordinated listening review.

---

## 22. Staged development roadmap

Keep future stages small and separately reviewable; assign no detailed milestone
numbers beyond the existing M13 proposal:

1. **M13: first bounded section-led algorithm.** Existing manual sections, fixed
   conservative preferences, no limiter/dynamics; private candidate measurement,
   explicit Apply and independent Current verification.
2. **Mono section processing and verified export.** Extend/qualify the shared channel
   path and format gates without bundling new automatic methods.
3. **Basic optional marker import.** Explicit time origin/rate and exact frame mapping,
   preview/confirmation and provenance; annotations only. Bounded EDL import follows
   as a separate proposal and is not a dependency for useful proposal generation.
4. **Standalone DSP qualification, in separate small steps:** first final true-peak
   limiter, then downward processing, then upward processing, each with exact shared
   offline/audition/export behavior, identity/activity, state/seek/timing/resource
   evidence on both Macs and listening review. Existing plugin code alone is insufficient.
5. **Dynamics-led proposals.** Qualified processors and bounded versioned activity
   preferences; useful without markers/sections. Explicit method selection/review.
6. **Enhanced section-led proposals.** Qualified optional final limiter, configurable
   section/transition/neighbour preferences, provenance/protection and later audio-
   derived regions as individually approved increments; no forced scene equalization.
7. **Mixed strategy.** Justified interaction of qualified section and dynamics stages,
   comparable against the same baseline with explicit bounded comparison costs.
8. **Later component-relative dialogue/M&E processing and validated objectives.**
   Keep dialogue-stem loudness, speech-gated mix measurement, dialogue/background
   diagnostics and intelligibility distinct. Verify the intended improvement and
   resulting complete programme; analysis-only references remain non-audible.

Generate the selected strategy only; alternatives/algorithm comparisons require an
explicit bounded request. Automatic region generation, broad optimization, general
freeze UI and ML remain separate future proposals. Further immersive expansion and
VST3 parity stay deferred.

The agreed M13 sequence is unchanged: Sequoia audio development first, freeze the
combined M12/M13 candidate, then actual Intel High Sierra and selected deferred video
checks. Reuse evidence only with justified compatibility, never as new-feature runtime
qualification. Accepted occasional picture skips retain recorded counts/timing;
audio interruptions, accumulating drift, sustained freezes and identity/lifecycle
failures remain blockers. Main stays at M11 until combined qualification and
subsequent explicit integration approval. The approved impact plan selects the work.

---

## 23. Recommended File Placement

Recommended documentation placement:

- main spec file: `Docs/ProposalEngineSpec_v2.md`
- add an index entry in `Docs/ProjectIndex/README.md`
- add a summary line in `Docs/ArchitectureSnapshot.md` under planned proposal/policy hooks
- later add code-level design notes under `Docs/ProjectIndex/` if the implementation splits into multiple tracked milestones

Reason:
- this keeps one canonical proposal-engine spec at the top `Docs` level
- keeps `ProjectIndex` as the navigation/index layer
- keeps `ArchitectureSnapshot.md` as the contract summary, not the full design document

---

## 24. Recommended Code Namespace / Folder Direction

Suggested future code grouping in Core:

- `Source/Core/Proposal/`
  - `ProposalEngineFacade.h/.cpp`
  - `MeasurementEngineBridge.h/.cpp`
  - `AnalysisEngine.h/.cpp`
  - `AttributionEngine.h/.cpp`
  - `PolicyEngine.h/.cpp`
  - `ProxySimulationEngine.h/.cpp`
  - `RefinementEngine.h/.cpp`
  - `VerificationEngine.h/.cpp`
  - `FreezeEngine.h/.cpp`
  - `ExplanationEngine.h/.cpp`

Suggested separate methods:

- `Source/Core/Proposal/Methods/RuleBasedMethod.h/.cpp`
- `Source/Core/Proposal/Methods/DistributionReshapingMethod.h/.cpp`
- `Source/Core/Proposal/Methods/CombinedAnalysisMethod.h/.cpp` (future design name; formerly `HybridMethod`, unrelated to the mixed processing strategy)

Suggested shared model types:

- `Source/Core/Proposal/Model/ProgramAnalysis.h`
- `Source/Core/Proposal/Model/SceneAnalysis.h`
- `Source/Core/Proposal/Model/ContributionMap.h`
- `Source/Core/Proposal/Model/ProposalIntent.h`
- `Source/Core/Proposal/Model/ModulePlan.h`
- `Source/Core/Proposal/Model/ProposalReport.h`
- future strategy/capability and versioned processing-preference records (design only)

These are design locations, not instructions to create or rename files now. Keep
strategy capabilities, selection algorithms and processing preferences separate
while preserving one shared truth model.

---

## 25. Design status

Three processing strategies share one authoritative foundation and may each host
multiple selection algorithms. M13 remains a bounded section-led subset with the
existing conservative defaults. Configurable preferences, neighbouring adjustment,
qualified standalone dynamics/limiting and comparisons are staged future work.
Bounded M13 implementation is now approved; broader strategies/controls remain
future work. Readiness review is complete and checkpoint publication/integration
was subsequently authorized; no broader strategy implementation is included.

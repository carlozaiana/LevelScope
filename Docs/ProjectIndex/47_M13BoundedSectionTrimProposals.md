# M13 proposal — first bounded section-led method and measured verification

Publication update, 28 October 2026: combined candidate-02 completed actual-platform
qualification and the user subsequently authorized commits, an annotated tag and
main integration. The [M13 checkpoint](https://github.com/carlozaiana/LevelScope/tree/m13-bounded-section-proposals-2026-10-28) records the accepted source.
The original development-phase restrictions below are retained as historical
plan/authorization boundaries; the later approval permits this publication.

Status: **approved for implementation and qualification**, 28 October 2026.
The two clarifications below preserve M13's bounded scope and conservative defaults.
Stop for review before checkpoint publication or integration into main.

Start from M12 candidate-03, checkpoint
[`m12-section-trim-sequoia-2026-10-27`](https://github.com/carlozaiana/LevelScope/tree/m12-section-trim-sequoia-2026-10-27),
commit `2fb0e3cfee22a527bacaa09c4ddace23129442fb`, in the established
`codex/logic-module` working copy. Main remains qualified M11
`cfb1033aa7fb65a2655db276bee696161196bdeb`. Preserve all checkpoint tags and packages.
The M12 checkpoint is **“Sequoia tested; High Sierra and video qualification pending.”**
M13 approval authorizes development and the agreed qualification, not a main merge
or checkpoint publication. The pre-implementation change-impact plan is [48](48_M13ImpactQualificationPlan.md).

Implementation/qualification update, 28 October: candidate-02 has completed the
approved combined actual-platform checks. [49](49_M13QualificationAndHandover.md)
records its exact identity, fresh/reused evidence and limitations. Publication and
main integration were subsequently approved after readiness review.

Read with [ProposalEngineSpec_v2](../ProposalEngineSpec_v2.md),
[M12 contract](44_M12SectionTrimProcessing.md), [M12 handover](46_M12QualificationAndHandover.md)
and [measurement contract](14_OfflineMeasurementContract.md).

## Recommended outcome

Offer one explained, deterministic section-gain proposal when a measured programme
cannot satisfy the selected loudness/true-peak rules with global trim alone.
Privately measure the proposed complete programme before offering Apply. Applying
it remains an explicit stopped-only operation; Current State then requires its own
remeasurement. Generation, preview and private assessment never change audible audio.

This first automatic method addresses **localized peak constraints**. It reduces
the level of explicitly annotated passages and compensates with final global trim
where feasible. It does not promise a solution, equalize scenes, infer dialogue,
steer LRA, compress samples or substitute a limiter. The operator reviews the
artistic balance change. M12 manual trims remain available when this bounded method
cannot help. The specification's broader optimizer and dialogue targets remain future work.

### Position within the three-strategy design

The shared [specification](../ProposalEngineSpec_v2.md#7-three-processing-strategies-and-intervention-preferences)
now distinguishes **dynamics-led** (global gain, upward/downward compression, final
true-peak limiting), **section-led** (global gain, section gains/transitions, final
true-peak limiting), and **mixed** (a justified combination). Available processors
need not be active. The broader no-change rule means **all applicable explicitly
selected objectives are satisfied**. A basic programme-rule pass must not suppress
a future explicitly requested validated dialogue/other objective. Such objectives
need defined methods, applicability and acceptance criteria, with identity/version
bindings; preferences cannot infer them or add compliance rules. M13's only objectives
are the currently supported selected programme rules: its already-pass no-op and
global-only-first behavior remain unchanged.

M13 is the first selection algorithm within the section-led strategy, with a
strictly smaller capability set: **global and section gain only; no limiter**.
Standalone limiter/dynamics qualification precedes their use in later proposals;
existing plugin DSP is not evidence for standalone compliance. A bounded-method
failure does not silently invoke another processor or strategy.

Strategy describes permitted interventions; algorithm describes how settings are
selected. `section-trim-proposal-v1` identifies this bounded algorithm/policy version,
separately from the strategy, profile and processing-engine versions. The old
specification's “hybrid” meant combined analysis techniques; use “combined-analysis
method” for that idea, reserving **mixed** for section-plus-dynamics processing.
Each strategy may later have multiple algorithms on the same foundation.

Generate only the selected method. Reuse compatible analysis by exact identity.
Future explicit strategy comparisons begin from the same frozen measured baseline,
not sequentially applied results, and declare their candidate/pass/time/memory/disk
budgets before work. M13 adds no strategy selector or multi-strategy comparison UI.

## Existing implementation and dependencies

| Foundation | Reuse / required addition |
| --- | --- |
| M2 profile and measured evaluator | Keep the authoritative EBU definition/digest and applicability. Obtain all rule limits from it; no new LRA rule or compliance allowance. |
| M3 `draftGlobalTrim` / `global-trim-m12.1` | Keep the global-first/no-change decision and fixed upstream section identity. A global-only solution remains the existing M3 workflow. |
| M4/M5 analysis and manual sections | Reuse full-resolution window data, source-frame intervals and annotation editor. `SectionAnalysis` is descriptive: its energy shares are not causal predictions, and its peaks contain only programme maxima. |
| M12 `SectionGainPlan` / shared `ProgrammeGainProcessor` | Reuse the exact 100 ms convention, float32 stage order, all-channel processing, audition/export path, sessions and measured verification. No new DSP stage or transport behavior. |
| M11 verified programme inputs | Reuse original component identities and deterministic composition. Section changes follow the sum; they never change dialogue relative to M&E. |
| New bounded proposal service | Add proposal-specific peak attribution, immutable candidate plans, private whole-programme assessment, explicit atomic Apply, cancellation and explanation UI. |

Every strategy must share profiles, measurement, analysis, exact processing,
identity/invalidation, independent Current verification and export gates. Each
reviewable plan must identify strategy/algorithm/policy versions, baseline, exact
settings, reasons, candidate measurements, processing activity and limitations.
M13 activity reports changed section/transition extent and baseline/additional/
cumulative gains; it must not report invented compressor/limiter activity.

No existing section statistic establishes the maximum peak of every section.
M13 needs an optional, bounded per-section peak collector using the **same continuous
true-peak interpolation and sample stream** as the qualified measurement engine.
Do not restart a meter at each section or substitute sampled plot points. Keep
this attribution separate from ordinary measured-result authority and do not
alter numerical accumulation, filter history, gates or EOF drain behavior.
If introducing the collector changes shared measurement code, the pinned EBU/ITU
checks and numerical equivalence tests are required on both actual platforms.
Section maxima identify affected sections and peak locations, not the exact number
of distinct peak events. Event collection/grouping would need its own threshold,
duration, gap, channel and interpolation/EOF rules. It is not part of M13.

## Exact scope and proposed policy defaults

- Whole-programme stereo and conventional 5.1, including LFE, at 44.1/48/96/192 kHz,
  up to three hours. Support qualified interleaved files, mono-file assemblies and
  all three M11 composition mappings. Preserve legacy mono/5.0/7.1 global-only behavior.
- Use existing manually confirmed annotations; no imported marker track is required.
  Do not create, move, merge or infer section boundaries automatically. Editorial
  scenes/edits, processing sections, loudness regions and peak events are distinct.
  Keep boundary provenance visible; an annotation or audio-derived region is not
  a verified cinematic scene. Future dynamics-led proposals need no sections;
  future section-led methods may use manual, imported or audio-derived boundaries.
- Require valid, current measured Source and Current State describing the exact same
  programme content, ordered roles, range and alignment, with Current matching the
  applied plan. Source need not pass. Confirm profile applicability.
- Require unapplied manual drafts to be resolved explicitly first. If an applied
  non-unity section plan exists, its annotation geometry must still match the
  proposed geometry. Otherwise stop with an explanation; do not remap, reset or
  overwrite manual gains. Cosmetic names stay outside audio identity.
- A versioned fixed policy, initially `section-trim-proposal-v1`: additional local
  attenuation only, **at most 3 dB per section**, **at most eight changed sections**,
  and **at most 20% of programme duration affected**, counting the union of changed
  sections and transition regions. Final plateau gains must also remain within
  M12's -12 to +12 dB range. Existing global limits remain unchanged.
- Include a visible **0.5 dB planning reserve** in requested local attenuation to
  allow for the loudness lost when reducing a passage. This is a product heuristic,
  not a standard tolerance, guaranteed headroom, or relaxed evaluation ceiling.
  The reserve is **inside the 3 dB total additional-reduction cap**, never 3 dB
  plus another 0.5 dB. Show requested reduction, applied reduction, achieved reserve
  and any uncovered deficit when the cap prevents the full reserve.
- One proposed section envelope per Generate action; no combinatorial search,
  automatic retries or background pursuit of a passing result. At most **three
  additional full-programme measurement/attribution passes**: attribution if needed, local-plan measurement,
  then measurement of the complete plan with its proposed final global trim.
  Reuse compatible attribution and skip an identical final pass. Explicit later
  Current remeasurement and export verification are separate user actions.
  Verified-copy preparation and content hashing are additional I/O, included in
  disk/runtime budgets rather than hidden within the three measurement passes.

These alteration budgets are deliberately conservative defaults for review, not
EBU requirements. They limit the change **relative to measured Current State**;
show existing, additional and resulting gains, and the cumulative applied envelope.
They do not prove minimum perceptual alteration. No automatic upward section gain,
per-section normalization to -23 LUFS, or new LRA target is introduced.

### Preference structure and neighbouring material

The [shared preference design](../ProposalEngineSpec_v2.md#17-configurable-processing-preference-design)
separates additional upward/downward correction, cumulative correction, neighbouring
plateau differences, transition duration/rate, affected extent, future compression/
limiter activity and protected-passage constraints from profile requirements.
Define that versioned structure now; do not implement configurable controls or
coordinated neighbouring-section adjustment in M13.

M13 retains zero additional upward gain, 3 dB maximum additional downward gain,
M12's cumulative plateau range -12 to +12 dB, eight changed sections, 20% affected
programme duration, fixed versioned 100 ms transitions and three private passes.
Neighbour differences and derived slew are descriptive; no new independent cap is
introduced in this revision. Compression and limiting remain disabled. Show both
per-request and cumulative settings/extent; repeated requests cannot conceal the
total applied correction. The 20% duration limit covers changed local section gain
and transitions, excluding uniform global trim; show global trim's whole-programme
extent separately. Do not claim that only 20% of samples change when global gain
also changes.

For example, **0, -2.5, -5, -2.5, 0 dB** has 5 dB total section attenuation but only
2.5 dB differences between neighbours. These are separate preference limits on
applied gain, not goals for equal measured scene loudness. Future smoothing must
account for changing otherwise satisfactory neighbouring material, not only smaller
steps. M13 keeps its existing transition rule and counts affected transition support;
it does not propagate gain adjustments into neighbours merely to smooth the curve.
Direct peak attribution to both sides of a boundary remains distinct from such
coordinated adjustment.

Future user protection instructions take priority over automatic suggestions.
Protection must specify whether it permits global gain and how transition/filter
support is handled. M13 adds no protection editor; if a supplied mandatory constraint
cannot be represented/enforced, explain that generation is unavailable rather than
ignore it when proposing changes. A true no-change outcome remains a no-op.
Measurement gates, future loudness-region boundaries and DSP thresholds
remain separate concepts with explicit future mapping rules, not interchangeable
numbers or automatic section cuts.

## Deterministic decision sequence

1. Validate snapshots, profile, applicability, required finite metrics, programme
   identity and job/session freshness. Missing metrics, stale measurements,
   unsupported inputs and predictions cannot enter the method as measured facts.
   Unavailable LRA is not a blocker when the profile has no LRA rule.
2. If Current already passes, return **NO CHANGE**: preserve all gains, revisions,
   valid measurements and evaluations; start no new measurement/attribution job.
3. Try M3 global trim with the existing section envelope fixed. If it offers a
   global-only draft, return that outcome. Do not apply section changes merely to
   move an acceptable result closer to exactly -23 LUFS. Enter the new path only
   for an identified loudness/peak-headroom conflict, not arbitrary invalid input.
4. Obtain a complete per-section/all-channel peak map for measured Current using
   a private verified-source pass if none matches its exact identity and geometry.
   Include LFE and muted channels. Record sample and true peaks, exact source
   frame/eighth-frame coordinates, and explicit EOF filter-padding attribution.
   Padding does not lengthen the displayed programme.
5. At the profile's intended global correction, identify every section whose
   mapped peak would exceed the strict profile ceiling. Attribute peaks within
   transition neighborhoods and the interpolation support conservatively to both
   affected sides. At an exact boundary use a documented deterministic rule;
   leading/trailing interpolation padding belongs to the first/last interval.
   Never hide a crossing peak by clipping its coordinate to a display axis.
6. For each affected section request the peak deficit plus the explicit planning
   reserve, capped at 3 dB; leave other section plateaus unchanged. Build the
   complete M12 envelope, validate resulting transitions and all budgets. An
   over-budget changed set is rejected, not silently truncated to eight sections
   or remapped. Impossible very short transitions are rejected. A capped candidate
   can still be assessed, but carries no predicted-pass claim. For deficit D > 0,
   request D + 0.5 dB and apply A = min(3, D + 0.5) subject to M12 validity.
   Achieved reserve = min(0.5, max(0, A - D)); uncovered deficit = max(0, D - A).
   Thus D=2.8 uses A=3.0 with 0.2 dB reserve; D=3.2 uses A=3.0 with no reserve
   and 0.2 dB still uncovered. Disclose each, without changing the strict ceiling.
7. Measure that local envelope at the existing final global gain through the
   qualified worker. From this **private measured candidate**, derive a final
   global correction using the existing M3 numerical rules. If no feasible global
   correction exists, finish with **NO SOLUTION FOUND WITHIN THIS POLICY**. This is
   not proof that no other plan could work.
8. Measure the final combined section/global candidate, including all transition
   samples, gating and every channel. Strict M2 evaluation determines whether it
   can be offered for Apply. A failed or unavailable assessment is shown as such;
   never promote an estimate to a pass. Do not try a second section envelope.
9. Present changed intervals, baseline/resulting gains, affected duration, method
   and profile identities, why global trim was insufficient, candidate measurements
   and limitations. Rank/explanation uses full-resolution data; plot reduction
   cannot affect the decision. Report section count/extent and maxima, not a count
   of distinct peak events. Equal inputs and policy give equal plans; stable
   frame/role ordering resolves ties. Do not claim a general optimum.

The peak map and reserve are planning aids. Exact whole-programme assessment is
necessary because filtering, transitions, gate membership and the final global
gain can change the result. Do not add measurement-engine test tolerances to rules.

### Failure explanation contract

Every no-solution/unavailable result identifies the strategy and method, baseline
and profile identities, the limiting rule or budget, units and configured bound,
observed supporting measurements, requested versus allowed correction, affected
sections/channels and whether a private candidate was actually measured. State
when the reserve was capped, attribution is unavailable, geometry/protection is
unsupported or the work budget was exhausted. An unmeasured estimate is never
reported as a measured failure/pass. Show actual post-candidate/post-Apply failures.

Examples are insufficient peak headroom for the allowable global gain, more than
eight sections, affected extent over 20%, a cumulative gain or transition conflict,
and a final measured loudness/true-peak miss. “No solution found within this policy”
does not mean no dynamics-led, enhanced section-led or mixed strategy could succeed;
none is automatically attempted or promised to succeed.

Future limiter choices will require excess magnitude, duration and distribution,
plus measured activity and listening evidence. Section maxima alone cannot decide
whether limiting many brief peaks is preferable to reducing a long passage; even
an event count would be insufficient on its own.

## Separate candidate, applied and measured state

Use an immutable proposal request containing source/current snapshot identities,
programme/member identities and alignment, ordered roles, range, base processing
revision/digest, annotation geometry, profile digest/applicability, policy/method
version, selected strategy/capability and algorithm identity, preference digest,
build/platform provenance and runtime session/job token.

Private measurements live in a distinct **candidate-assessment domain**. They may
say “Candidate passes selected programme loudness and true-peak rules”; they are
neither Current State measurements nor verified exported audio. They cannot satisfy
Current freshness or M10 publication gates and are not installed into session caches.
Keep predictions explicitly labelled as predictions. Candidate cancellation, read
failure or exhausted budgets leave the previous valid session and audio untouched.

**Apply proposal** installs the proposed section envelope and final global gain
atomically, once, while stopped and with all identities rechecked. It invalidates
affected Current measurements/evaluations, other proposals and export approvals;
Source stays unchanged. Current is immediately **UNVERIFIED**. Require a separate
explicit **Re-measure Current State** through the qualified worker and M2 evaluation,
even if the identical private candidate was measured successfully.

If actual verification disagrees, show the actual values and verdict and stop.
Cancellation/failure leaves the applied plan unverified; completed failure remains
a measured failure. Retain applied settings; no rollback, retry or extra correction.
An identical plan is a genuine no-op and preserves revisions/results.

Source/member/role/processing changes reject dependent drafts and late jobs. An
annotation geometry edit invalidates its proposal/attribution but not unchanged
applied audio or measurements. Names may update display labels only. A compatible
profile change reuses valid audio measurements but requires fresh policy evaluation
and proposal binding. Changed processing preferences or algorithm selection
invalidate dependent proposals, not unchanged audio measurements. Monitor
volume/routing and video attachment/offset preserve
audio identity; their established playback restrictions still apply.

Save/reopen persists the applied M12 plan and exact settings using schema 4; it
does not make a transient proposal or private assessment authoritative. Reopened
sessions stay stopped, device-closed and routes-muted. Follow M12's explicit
Apply/Discard/Cancel handling for outstanding drafts. No proposal-history browser
or general undo system is required. Existing reports describe actual Source,
Current and exported audio; strict Current-plus-artifact publication gates and
one final TPDF quantization remain unchanged.

## Minimal UI

Add Generate section proposal / Cancel, a concise interval-and-gain diff, separately
labelled candidate measurements, and explicit Apply. Reuse the existing timeline
for interval highlighting and draft/applied curves. Navigation never starts audio.
Show why generation is unavailable or no bounded solution was found. A candidate
that has not passed its private assessment cannot be applied via this automatic
workflow; the independent manual M12 editor remains available.

Do not redesign the application, add candidate audition, make device opening
automatic, or implicitly generate a proposal when another job completes.

## Acceptance tests

| Area | Required evidence |
| --- | --- |
| Decision hierarchy | Already-pass preserves settings/revision/results with zero worker calls; feasible global-only case keeps every section unchanged; localized-peak fixture yields an explained successful section proposal; distributed peaks/budget exhaustion produce honest no-solution without mutation. |
| Numeric policy | Full-precision profile comparisons and gains; exact/adjacent threshold values; reserve included inside 3 dB, full/partial/missing reserve and uncovered deficit disclosed; additional versus cumulative and neighbouring-difference examples; eight-section, 20% and M12 range boundaries; no epsilon/rounded-display verdict changes; deterministic ties and repeated identical requests. |
| Attribution | All-channel per-section sample/true peaks against an independent short-fixture oracle; multiple peaks, near-boundary peaks, shortened ramps, first/last-frame impulses, EOF drain, silence and non-finite rejection. Collector on/off must preserve every ordinary measurement result; section maxima must not be labelled event counts. |
| Input coverage | Four rates × stereo/5.1 × interleaved/grouped/composed representations, all three composition mappings; LFE or a muted channel alone limits the plan. No channel removal or component rebalance. Unsupported layouts fail clearly. |
| Identity | Mismatched Source/Current content, roles, alignment, range or plan rejected; stale map/profile/request/late completion rejected; rename and video/monitor-only edits preserve audio; annotation movement changes no applied samples. |
| Manual baseline | Existing nonzero gains survive; geometry mismatch never resets them; no pending draft overwritten; budgets count complete intervals and ramp unions; applying section plus global is one transaction. |
| Truth domains | Prediction/private-candidate pass cannot become Current/artifact pass; explicit Apply invalidates old Current truth; independent remeasurement succeeds or displays actual failure; cancel/fail retains applied settings unverified with no automatic action. |
| Audio and lifecycle | Independent sample oracle, offline/pre-monitor/export parity including all channels and transitions; random seeks/block sizes, duration/frame position/EOF, overload recovery, interrupted jobs and owned-copy/write/disk failures. Device remains closed for injected tests. |
| Sessions and export | Exact applied-plan round trip, schema 1–4 compatibility and cache recovery; coherent saving/late jobs; narrow report domains/hashes; final PCM24 single dither and actual-artifact gates, including a deliberate failing Current/artifact. |
| UI | Generate/Cancel, explanation and candidate labels, draft/Apply/UNVERIFIED/actual result, no-change, no-solution, stale actions, reasoned limiting-constraint failure, boundary provenance, operation without imported markers, save/discard, small-window reachability; no accidental playback or device opening. |

Use known synthetic programmes with successful, infeasible and misleading proxy
cases; include a deliberately injected disagreement between private assessment and
post-Apply verification to prove it stops rather than silently correcting. Frozen
fixtures must exercise the actual success path, not only graceful rejection.
Cross-platform tests use the established measurement tolerances, separate from
strict rule thresholds. Determinism means identical measured inputs/policy produce
identical plans; near-threshold cross-platform verdict differences must be reported,
never hidden by increasing the profile allowance.

## Planned local cinema evaluation

Local original, to preserve unchanged and keep out of GitHub:
`/Volumes/Tucan/GuerraDiCesare/LGdC_Masters/LGdC__51_Master_V2.wav`.
The single bounded evaluation has now completed; [49](49_M13QualificationAndHandover.md)
records the PCM-identical container rewrap, actual measurements and honest
whole-programme duration-budget failure.
No external observation below is attributed to LevelScope. The user confirmed on
28 October that this test should retain **one whole-programme section**; no automatic
scene boundaries are to be created. A duration-budget failure is an acceptable outcome.

The intended external comparison is the statistics supplied by the user on
**28 October 2026**, reported as measured with the latest Nuendo version available
to them. Do not invent an exact version. The user identifies the file above as the
source of these **external reference observations**:

| Field | User-reported Nuendo observation, 28 October 2026 |
| --- | --- |
| Channels | L, R, C, LFE, Ls, Rs |
| Sample rate / bit depth | 48 kHz / 24 bit |
| Integrated loudness | -26.0 LUFS |
| Loudness range | 19.9 LU |
| Maximum true peak | -0.85 dBTP |
| Maximum momentary / short-term | -5.4 / -9.8 LUFS |
| Dialogue-gated loudness / range | -28.3 LUFS / 17.4 LU |
| Detected dialogue | 28% |
| Per-channel true peaks, stated order | -0.85, -1.43, -1.00, -1.98, -1.96, -1.99 dBTP |

The historical screenshot
`/Volumes/Tucan/GuerraDiCesare/LGdC_Masters/Screen Shot 2024-04-03 at 21.24.08.png`
was produced in **Cubase 10**, not Nuendo. Retain it separately as supplementary
historical evidence. It displays -25.98 LUFS, 19.94 LU LRA and -0.85 dBTP for
`LGdC_51_Master_V2`; do not substitute those rounded readings for the supplied
Nuendo export. The other historical stereo/Insight screenshots remain distinct.

Confirm exact Nuendo version, relevant settings (including dialogue gating/detection),
measurement range, content hash and ordered channel interpretation when needed for
the later comparison, without blocking implementation. These are external observations,
not LevelScope results or exact reference vectors. Do not treat differences as engine
errors until provenance/compatibility is established. Dialogue-gated values and the
28% dialogue figure remain contextual future-work references, not M13 objectives.

For illustration only, the supplied whole-programme values suggest a scalar +3.0 dB
move toward -23 LUFS and +2.15 dBTP predicted peak, 3.15 dB over the strict ceiling.
The requested local deficit plus reserve is 3.65 dB before the 3 dB cap: the initial
estimate leaves no reserve and 0.15 dB uncovered. This is not a final section-plan
prediction: gates, local extent, transitions and final global correction require
actual measurements. It supplies no event count, clipping claim, LRA failure or
proof that every strategy is infeasible.

For the approved single bounded Sequoia evaluation, verify actual content, roles/range and the external
statistics provenance, measure Source/Current with the candidate, retain one bounded
proposal attempt and its exact success/failure explanation, and review the artistic
trade-off through a separately coordinated listening check. Do not alter the original,
relax the profile, broaden the policy or repeatedly retry to force this film to pass.
An honestly explained bounded-method failure is a valid evaluation outcome; synthetic
fixtures must still demonstrate successful proposal/application/verification.
Only text evidence, measured identities and explanations belong in version control;
raw media and the local screenshots remain outside GitHub.

## Qualification sequence and resource limits

Before implementation, record an explicit change-impact table covering shared
callers, selected executable modes, reused identities and estimates. Short targeted
development tests first; no recursive historical qualification runners. Build only
named standalone/test targets with the established universal/10.13 dependencies.

1. Stabilize M13 on Sequoia with audio-only automated and silent UI checks. Preserve
   the M12 candidate and keep main at M11. No Golem or video loading/decoding/playback
   during this development stage, per the existing deferral. Explain any concrete
   timing/buffering/video change before requesting an earlier change to that plan.
2. Freeze a new combined M12/M13 candidate, recording source and executable hashes.
   Then complete actual Intel High Sierra runtime and cross-Mac session checks on
   that same candidate, including all new manual and automatic section behavior.
3. Run pinned EBU/ITU and affected M1–M12 checks on both Macs, with collector on/off
   equivalence and focused legacy mono/5.0/7.1, group/composition and schema coverage.
   Reuse historical evidence only with a documented dependency/behavior equivalence;
   M11 results never qualify new M12/M13 behavior on High Sierra.
4. One consolidated three-hour 48 kHz 5.1 composition workflow per Mac, with 1,000
   annotations and localized peak constraints: attribution, bounded proposal,
   explicit Apply/verification, seek/capture, diagnostics, sessions and real verified
   RF64 export. Do not separately rerun equivalent older long workflows. Exercise
   10,000-interval storage/boundary budgets in a short synthetic test.
5. After the freeze, perform the deferred combined video checks on both actual Macs:
   H.264 CFR and VFR, and ProRes exact- and fractional-rate fixtures, each in Source
   and Current (**eight short cases per Mac**). Include a nonzero section envelope
   plus global trim, positive/negative rational offsets, seeks across transitions,
   pause/resume/restart and EOF. Use injected audio first; coordinate a representative
   muted physical-device subset before opening devices. Preserve the audio clock
   and embedded-video-audio exclusion. No automatic full codec/colour matrix or
   30-minute endurance suite unless changed dependencies or a concrete failure
   justify it. Record selected coverage and the original evidence reused.

Planning estimate: short/reference checks 30–60 minutes per Mac, silent UI 20–40
minutes, and short video checks 15–30 minutes, excluding builds/fixes/operator delays.
M12's Sequoia consolidated workflow took 70 minutes; added exact candidate scans
make the proposed long-workflow estimate **90–150 minutes on Sequoia, 120–210 minutes
on High Sierra**, to be measured, not assumed. Proposed hard limits are **240 minutes
on Sequoia and 300 minutes on High Sierra**, including the combined long workflow;
missed limits remain qualification failures, not silently widened budgets.

Keep the generation pass cap at three. Reuse one verified source preparation
through the private passes; stream samples without retaining full-length renders.
Require at least **65 GiB free working headroom per long run**, checking actual
input/copy/export sizes before starting and avoiding duplicate original packages.
Proposal/attribution/plan data <=16 MiB at 10,000 sections; additional working RSS
<=256 MiB beyond existing measured-result storage; UI heartbeat p95 <=100 ms;
cancellation acknowledgement <=1 second; save/parse <=5 seconds; session <=128 MiB;
JSON+HTML report <=16 MiB. Track peak disk/RSS and owned cleanup on success/failure.

Occasional picture skips remain accepted development limitations, with counts and
timing retained. Audio interruptions, accumulating drift, sustained freezes,
identity errors and memory/lifecycle failures remain blockers. Do not repeat a run
merely to seek zero skips. Full physical multichannel listening remains separately
pending; offline/automated qualification does not require eight outputs. Coordinate
audible tests separately; preserve installations, plugins, Pro40 and system audio.

The completion validator must explicitly require each selected test/receipt and
candidate identity; missing required evidence cannot pass. Until actual High Sierra
and selected video checks are complete, retain **“Sequoia tested; High Sierra and
video qualification pending.”** Report combined readiness for review; do not merge
main or publish a new checkpoint without subsequent authorization.

## Exclusions and following order

No dynamics/limiting, lower-tail or LRA steering, per-component gains/DSP, dialogue
analysis, semantic segmentation, optimization search, section freezes, general
policy editor, marker/EDL import, export-format expansion, 7.1.4/Atmos or VST3 parity.
Mono section processing and mono verified export remain a near-term separate
extension; their existing measurement/global-trim behavior must keep working.

Following the combined milestone, use separately reviewed stages, without assigning
future milestone numbers prematurely:

1. Mono section processing and verified export.
2. Basic optional marker import; a later bounded EDL proposal adds explicit timecode/
   frame mapping and preview. Importers create annotations, never automatic gains.
3. Standalone DSP qualification in small steps: final true-peak limiter, downward
   processing, then upward processing, including activity, identity, state/seek and
   shared measurement/audition/export parity. Existing plugin DSP is insufficient.
4. Dynamics-led proposals, with useful marker-free operation and bounded activity.
5. Enhanced section-led proposals: qualified optional final limiting, configurable
   preferences and later protection/coordinated-neighbour/audio-region work as
   separately approved increments.
6. Mixed strategy, with justified combinations and explicit resource-bounded
   comparisons from the same baseline.
7. Later component-relative dialogue/M&E processing and validated dialogue/intelligibility
   objectives. Verify intended improvement and the resulting complete programme;
   analysis-only references remain non-audible and independently invalidated.

The shared specification defines this roadmap. Preserve original accessible
components, original Source composition and extensible versioned Current plans.
Useful proposals never depend on having imported editorial markers. No future
strategy, controls, neighbour optimizer, DSP or importer is implemented by this
revision or added to M13's bounded scope.

## Approval boundary and future choices

M13 implementation and resource-aware qualification are complete. The later
28 October authorization permits publication and integration; no next milestone
is started by this checkpoint.
Preserve the published M12 package and all earlier evidence. No configurable controls,
coordinated neighbour adjustment, automatic segmentation, marker import, limiting,
compression or multi-strategy comparison is included.

Future neighbour/slew preferences, dynamics activity budgets, protection semantics,
event grouping and detector mappings remain separate proposals. Exact external
Nuendo settings/version will be confirmed for comparison without blocking development.
No dialogue objective is inferred from the supplied contextual statistics.

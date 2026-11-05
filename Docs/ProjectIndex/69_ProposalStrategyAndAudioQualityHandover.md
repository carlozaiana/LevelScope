# Proposal Strategy and Audio Quality — post-M19 design handover

5 November 2026. **Design and analysis only.** The user approved the development
approach below and this documentation handover. M19 publication and main integration
were subsequently authorized and are being prepared as `m19-lra-passages-2026-11-05`.
That authorization does not extend to new DSP, strategy implementation or test work.
The next design chat's first deliverable remains a bounded proposal for review.

## Starting point and ownership

Use the existing logic checkout `/Users/fuzz/Developer/LevelScope`, branch
`codex/m19-lra-passages`, for reading these documents and implementation. Do not
start from the older main and mistake missing M19 code for unfinished implementation.
No new checkout is necessary for this design assignment. Do not edit the separate
UI or video worktrees. Coordinate shared specification edits rather than concurrently
modifying implementation owned by another chat.

| Identity | Exact revision/status |
| --- | --- |
| Qualified published main | `51e47f300fcc38a386a5a92a7277029175cfab0c` — M18, unchanged |
| M18 checkpoint | `m18-interactive-timeline-2026-11-03`, target `a9550d9d43ed70fde4c73f90e44081a706d30a67` |
| Published Intel optimisation | `intel-multicore-true-peak-2026-11-05`, target `c4eea6f1c4569ec4f31725694ddf34db8340f15e`; receipt `09c8b8f6345c7c889f21f5399a38284ce1d353b2` |
| M19 candidate-01 product | `2da3ac4ccb03a75976403c03bbc3461d00639752` — original numerical qualification |
| M19 candidate-02 tested product | `cf275afb133118dd2a06e7bca59f82f98898b3cd` — corrected preview legend, same numerical implementation |
| M19 qualification/review documentation | `33e2248cb119290aaa7f16e03bb31faffbbbd0f0` — final evidence, supplemental helpers, completed listening review; no later product change |
| Original design handover | `b475134e9b7e8e2afe0f41ea3edca7c27a196aa6` — documentation-only continuation, not a new tested binary |
| Publication boundary | User authorized M19 publication and main integration after readiness review; preparation in progress. Preserve all 23 pre-existing checkpoint tag objects and targets. |

Frozen candidate-02 source-manifest SHA-256:
`0cc78f242abc6b1b799c04c2154fa25fb58ca57be295af2e047866c7d0a431e7`.
Measurement build: `993f26a9295201a87fedea54a60bc8ca57877e93bb8bebc11729ab5bf5d2af1a`.
Universal executable SHA-256: `4bab04604bce5e7137dfb103ef880dd127a6cb402844281473307775d1c42536`.

- Sequoia app: `/Users/fuzz/Developer/LevelScope/build/m19-qualification-2026-11-05/candidate-02/LevelScope Standalone.app`.
- High Sierra app: `/Users/fuzz/LevelScope-M19-qualification-2026-11-05/candidate-02/LevelScope Standalone.app`.
- Both `candidate-01` siblings remain preserved, with their own identities and receipts.
  Do not rebuild over, modify or relabel either package.

The advisory chat independently reran only the existing read-only completion audit
to `/tmp/levelscope-m19-advisory-review.json` and reported a pass. That is readiness
review, not fresh execution qualification or publication authorization. Original
evidence remains under `Tests/Reference/results/m19-qualification-2026-11-05`.

## Reading order and current versus historical status

1. [Shared proposal architecture](../ProposalEngineSpec_v2.md), especially truth domains,
   strategies, activity/protection, current section 8.4 and current roadmap in section 22.
2. [58 — selective-LRA/timeline priorities](58_LraWorkflowAndInteractiveTimelineDesign.md)
   and [64 — implemented M18 UI/timeline handover](64_M18UITimelineHandover.md).
3. [56 — current preference design](56_M16ProposalPreferencesAndLimits.md),
   [60 — job requirements](60_M17JobLraRequirementsAndPresets.md), and
   [14 — measurement contract](14_OfflineMeasurementContract.md).
4. [66 — exact M19 method and approved deterministic contract](66_M19LraPassagesAndMeasuredCandidates.md).
5. [67 — impact-selected test plan](67_M19ImpactQualificationPlan.md),
   [68 — qualification/results/limitations](68_M19QualificationAndHandover.md), and
   [evidence index](../../Tests/Reference/results/m19-qualification-2026-11-05/README.md).

The old roadmap statements “main at M13”, “M18 not implemented”, the M17 blanket
LRA-proposal restriction, fixed eight-section defaults, and the pending combined
High Sierra/video/three-hour qualification are historical. M18 completed that combined
qualification; M19 adds the bounded method below and has actual-platform qualification.
This does not mean every earlier future feature is implemented. General region editing,
freehand automation, rolling LRA/cumulative IL estimators, freezes, neighbour strategies,
standalone dynamics/limiting, multi-strategy comparison and dialogue analysis remain proposals.

## Approved development approach

Develop a processing capability with its first useful selection method through small
successive slices. Do not wait until every DSP module exists to design rules. Keep one
shared design separating exact processors, selection methods, hard delivery requirements,
user preservation preferences/resource budgets, and measured/perceptual assessment.

For each capability: identify a concrete mixing problem and success/failure examples;
qualify its processor with sufficient manual/test control; expose full settings and
actual activity; implement one bounded proposal method; measure and listen; revise the
rules from evidence with versioning and affected requalification. No general processing
graph refactor or complete dynamics suite is a prerequisite.

| Owner | Responsibility |
| --- | --- |
| Logic Module chat | Processor/algorithm implementation, shared identity/session contracts and qualification |
| UI/Timeline chat | Reviewable controls, timeline/activity displays, domain/state distinctions and explanations |
| Proposal Strategy and Audio Quality chat | Evidence-led rule design, preservation trade-offs, examples, decision records and bounded implementation handovers |

Keep section-led, dynamics-led and mixed choices distinct. Version algorithms separately
from processors, requirements and preferences. A more capable algorithm should reuse
qualified measurement, execution, Apply/history, UI and export rather than replace them.

Numerical compliance and preservation of mix quality are separate assessments. Record
the user's mixing judgments as examples; convert them into measurable constraints only
with an explicit definition and evidence. A scalar alteration score is neither artistic
transparency nor intelligibility. Neither can be inferred from a narrowed histogram.

**Analytical and predictive design is included in “evidence-led”.** Develop calculations
for a proposed region change together with surrounding material, overlapping windows,
transitions, existing processing and programme-wide gates, forecasting IL/LRA/true peak
and other selected metrics. Model-based selection, bounded simulation/search and
coordinated-neighbour methods are open design directions. This clarifies design scope;
it does not expand M19's implemented attenuation-only method during publication.

For every forecast, identify required data and identity, assumptions, units, stage and
resource cost, and distinguish exact conditional results from approximations/bounds.
Fixed-gain scaling can be exact in an ideal linear, fixed-population model; real gate
changes and float-stage rounding still need assessment. Regional window shifts require
overlap and boundary/filter context; whole-programme LRA needs the full population and
recomputed gates/percentiles. Arbitrary DSP/time-varying true-peak forecasts need sample/
interpolation/state context or a justified conservative bound. Calibrate predictions
against measured examples and show errors/failure cases. See the explicit model/data
table in [specification section 8.4](../ProposalEngineSpec_v2.md#analytical-and-predictive-design-is-evidence-led-design).
Actual complete candidate assessment, independent Current verification and finalized-
artifact verification remain acceptance authority. No model restores predictions as
Current, relaxes profile rules or eliminates those verification steps.

## What M19 actually does

Changing LRA proposals support confirmed stereo/conventional 5.1 at 44.1/48/96/192 kHz,
including existing interleaved, separate-file and confirmed M&E/dialogue compositions.
Every programme channel, including LFE and monitor-muted channels, participates in
processing and relevant peak constraints. Existing mono/7.1 workflows remain; this new
automatic LRA method is not qualified for them. No new DSP is active.
The authoritative EBU base remains −23 LUFS with its explicit ±0.2 LU QC allowance
and strict −1 dBTP ceiling. Job LRA is optional/unset by default, separately identified
and never an implicit EBU rule. There is no three-strategy selection UI yet: selecting
an LRA requirement dispatches the supported M19 method; the existing peak-led workflow
remains when LRA is unset.

`lra-passages-v1` consumes the complete validated LRA-window population, retaining
availability, stability, inclusion and EOF-padding information. Low passages contain
included windows at/below P10. High passages contain included windows strictly above
P95 when descriptive, or above `P10 + selected maximum LRA` when a requirement exists.
LRA has no universal EBU ceiling. These are passages to investigate, not additive LRA
contributions, defective scenes or proof of what an edit will accomplish.

Windows have actual three-second source-time support. Merge overlapping/touching
supports; keep genuine gaps separate; clip to programme EOF while identifying padding.
Rank high regions by sum of positive loudness excess above the search level, in original
window order (units LU × windows), then decreasing window count, then earlier begin/end.
Stable percentile ordering uses original window index for ties. Marker names/boundaries
are context; adding redundant annotations cannot change the default audible envelope.

`lra-selective-trim-v1` requires valid matching measured Source/Current content, ordered
roles/range, current applied plan, confirmed applicability and usable selected IL/TP/LRA
rules. Required metrics must be available and LRA stable; additional selected M/S ceilings
also apply. Pending manual drafts must be resolved. Its sequence is:

1. If all selected requirements pass, preserve settings, revision and measurements:
   genuine no-op. Otherwise, when LRA already meets and a scalar gain intersects the
   IL/TP/M/S bounds, actually assess that global-only plan first.
2. Otherwise request each high region's selected-window P95 minus the search level,
   plus the planning reserve. Cap attenuation by both additional headroom from measured
   Current and least cumulative headroom across intersecting applied plateaus. Disclose
   capped reserve and uncovered overshoot. Full precision applies, including representable
   inward adjustments at limits; no acceptance allowance is added.
3. Visit ranked regions greedily. Keep a complete region only if the canonical change
   support, including transitions, fits count/duration/capacity limits. Do not shorten it
   to fit or coordinate its neighbours. Partition prior plateaus, then validate exact
   linear-ramp endpoint differences: no additional local boost or excess attenuation.
   Incompatible existing transition geometry is a reported method limitation.
4. Use at most **three complete private measurements and one verified source copy**.
   Measure the first envelope at baseline global gain. If needed, one refinement deepens
   only those selected regions by the remaining measured LRA deficit, within original
   caps. Reserve an assessment for the retained local plan plus a changed final global
   gain. No unmeasured final plan can become Apply authority.
5. Rank measured attempts lexicographically: all-selected-rule pass; smaller LRA excess;
   fewer other failed rules; smaller sum of other positive bound violations; smaller
   affected local duration; smaller maximum additional attenuation; earlier attempt.
   Missing required metrics cannot win. Deeper attenuation is not assumed better.
6. The final global estimate is nearest the IL target inside IL/TP/M/S intersections
   and total −60…+24 dB. Measure the complete combined plan. Positive global gain raises
   the whole mix, including quiet material: it is not selective quiet-passage boost.
   Local duration/gain budgets must not be advertised as whole-chain alteration limits.

Current defaults are `section-preferences-v2` / `connected-envelope-support-v2`:

| Preference | Default | Supported control range / distinction |
| --- | ---: | --- |
| Changed-region count | No count restriction | Optional 1–10,000; underlying processing capacity is still 10,000 intervals |
| Affected local-envelope duration, including transitions | 20% | 0–100%, 0.1 percentage-point steps; global trim's whole-programme extent is separate |
| Additional attenuation from frozen Current | 3 dB | 0–6 dB, 0.1 dB steps |
| Resulting cumulative local attenuation | 12 dB | 0–12 dB, 0.1 dB steps; not reset by repeated proposals |
| Requested planning reserve | 0.5 dB | 0–1 dB, 0.1 dB steps, no greater than additional cap; inside both caps |
| Transition | `section-gain-v1`, nominal 100 ms | Fixed linear-in-dB convention, shortened symmetrically by adjacent interval lengths; not a new UI control |

The existing processor permits −12…+12 dB local plateaus; M19 adds attenuation only
relative to the frozen applied envelope. Canonical equal-gain boundaries disappear.
Transition radius is the minimum of rounded 50 ms and half each adjoining interval;
endpoint/sample rounding and float32 stage order are versioned, not aesthetic defaults
the planner may silently replace.

## Truth, identity and application safeguards

Generation/inspection/cancellation leaves Current unchanged. A changing result binds
both measured snapshots, exact programme/member identities, ordered roles/range,
current processing/revision, requirement/profile/applicability, method and policy.
Source replacement or stale jobs discard authority; reopened runtime identity rejects
late completions. Compatible requirement changes reuse measurements but reevaluate and
invalidate proposals/export approvals. Preferences invalidate proposals, not unchanged
measurements. Video offsets and monitoring stay outside programme audio identity.

Explicit stopped-only Apply is available only for the exact complete measured passing
candidate. It makes Current unverified; independently remeasure Current. A cancelled or
failed verification retains applied settings, with no rollback/retry. Source remains
unchanged. Historical explanations retain method, baseline, requirements/policy, exact
before/after plans and reasons; documentary private metrics never become reopened Current
verification. Subsequent edits mark the history's association accurately.

Export still renders exact Current, dithers once at final PCM24 quantization, reopens
the finalized artifact and measures it without additional gain. Both Current and artifact
must pass all selected requirements before explicit local publication. Preview, private
candidate, monitor output or historical verdicts cannot satisfy either gate.

## Compact measured evidence to discuss

Numerical receipts below are original candidate-01 evidence, reused for candidate-02's
legend-only change with documented equivalence. Fresh candidate-02 native UI and extra
5.1 rates passed on actual Sequoia/Intel High Sierra. Values here are rounded for reading;
selection/verdicts use full precision. Requirement ceilings are test choices, not EBU defaults.

| Example and evidence | Actual outcome | Design implication |
| --- | --- | --- |
| [65-second stereo success](../../Tests/Reference/results/m19-qualification-2026-11-05/sequoia-01/stereo/proposal.txt), 48 kHz, LRA ceiling 4.9 LU | LRA **5.910634 → 4.290034 LU**; one local region, **1.620600 dB** attenuation, **12.153846%** affected; final global **+6.531787 dB**; two private passes, all selected rules pass, then independent Current/artifact checks pass | Useful bounded correction exists. The global lift materially changes listening level, so local attenuation alone does not describe the audible change. |
| [Cinema 27:00–35:00 excerpt](../../Tests/Reference/results/m19-qualification-2026-11-05/sequoia-01/real/proposal.txt), 5.1/48 kHz, LRA ceiling 15 LU | **23.435333 → 22.385113 LU**, still fails; **19.0625%** affected, **3 dB** local cap, other complete regions exceed remaining duration allowance; two private passes, final global **+0.143883 dB**; no Apply, Current unchanged | Gain/duration restrictions and this method's search limit are visible. Do not enlarge defaults solely to make this film pass or assume limiting will resolve its LRA. |
| [Ten-minute gate-crossing fixture](../../Tests/Reference/results/m19-qualification-2026-11-05/sequoia-01/gating/proposal.txt), LRA ceiling 6.7 LU | Baseline **10.000000 LU**; first attempt **7.000000 LU**; deeper refinement **25.989999 LU** after gate/population change; final global-adjusted attempt about **7.000000 LU** also fails. Exact rank retains attempt 1; three passes, no Apply | Histogram intuition and deeper attenuation are insufficient. Preserve the better measured attempt; discuss priorities among failed attempts without weakening full-precision pass rules. |

Cinema provenance: original
`/Volumes/Tucan/GuerraDiCesare/LGdC_Masters/LGdC__51_Master_V2.wav`, SHA-256
`1fb7f93a33d600bf0a62d63e642e5703f03aec8e6920f8a6e510cf65f61fdeaf`.
Excerpt original frames `[77760000,100800000)`, byte-identical PCM24, ordered
L/R/C/LFE/Ls/Rs; excerpt SHA-256
`cfd57bea9a108db850b4e371de0f299f3396fbc1730e58e357380042a3064dd5`.
Selected region is excerpt `[8784000,13171200)` (183.0–274.4 s); transitions extend
the affected union. Initial requested reduction 9.300182 dB includes reserve, but only
3 dB is permitted; achieved reserve is zero and initial uncovered overshoot 5.800182 dB.
These are excerpt findings, not full-film conclusions or comparison to Nuendo's external
observations. “No acceptable candidate in this bounded attempt” is not mathematical
impossibility. The earlier whole-film/marker examples retain their distinct identities.

Listening is also evidence, separately from those metrics. Coordinated Saffire playback
used outputs 1/2, first −5 dB then a requested 0 dB monitor repeat. Synthetic sine
transitions/endings had subtle audible effects the user judged normal. Cinema had no
reported clicks, but quiet Left/Right content made the uncued sequence ambiguous.
The subsequent independent supplied-file inspection concluded:
**“I checked all 4 stereo and 5.1 files and they are flawless.”** This completes the
bounded M19 boundary review. Inspection method/hardware were not specified; full physical
multichannel routing/listening qualification remains separate. Do not relabel partial
Left/Right playback as 5.1 listening, or generalize the clip judgment to all programmes.

## Implementation interfaces and reusable material

Paths below are relative to the established checkout; they are reading entry points,
not instructions to activate new processing.

| Entry point | Current contract / reuse boundary |
| --- | --- |
| `Source/Standalone/StandaloneServices.h` | `ServiceContext`, half-open source `FrameRange`, immutable `ProcessingPlan`. `graphJson` is reserved; the qualified gain path rejects nonempty graphs. New DSP requires a versioned supported plan and migrations, not arbitrary plugin-state injection. |
| `Source/Standalone/StandaloneProcessing.h`; `Source/Core/Processing/SectionGainPlan.h` | `ProgrammeGainProcessor::process(..., firstFrame)` shares exact composition → local envelope → global trim sample operations across callers; all channels/LFE, double coefficients and float32 stage rounding, no clipping/monitor gain. Current envelope is stateless and seek/chunk invariant. |
| `Source/Standalone/VerifiedSource.*`, `ProgrammeSource.*`, `ProgrammeComposition.*` | Verified programme/member content and confirmed assembly/composition; `prepareVerifiedSource`, owned-copy lifecycle, failure/cancellation and readers. One-copy reuse within M19 assessment is not session-wide sharing. |
| `Source/Core/Proposal/LraPassages.*`, `Source/Standalone/StandaloneLraProposal.*` | Full-resolution grouping/selection, `validateLraEnvelope`, `lraFinalGlobalGain`, `lraAttemptRank`, actual candidate worker. Keep method-version and measured-domain boundaries. |
| `Source/Standalone/StandaloneSectionProposal.*`, `StandaloneSessionModel.h` | `SectionProposalRequest::binding`, immutable result, `generateLraProposal`, `sectionProposalRequest`/acceptance/staleness and `applySectionProposal`. Source/Current mismatch, pending drafts, runtime changes and cancellation reject changing authority. |
| `Source/Core/Proposal/VerificationEngine.*`, `Model/ProposalPolicy.h`, `Model/JobRequirements.*` | Authoritative full-precision evaluations and separately versioned rules/preferences. New methods must reuse these gates. |
| `Source/Standalone/AppliedProposalExplanation.*`, `AppliedProposalSerialization.cpp`, `SessionDocument.*` | Last applied historical explanation, exact before/after plan identities, baseline and method/requirement/policy provenance; bounded separate attachment (8 MiB), replaceable measurement caches, no private-candidate authority restoration. |
| `Source/Standalone/LraPassagesPanel.*`, `TimelineGainLane.*`, `MeasurementTimeline.*`, `ProposalWorkflowPanel.*` | Temporary estimated overlay, actual assessed plan/reasons, draft/applied/verified distinctions. `onAssess` is explicit; selection highlights/inspects without seek/play. Ruler seek remains explicit and stopped-only. Inspection domain and audition domain differ. |
| `Source/Standalone/StandaloneAudition.*`, `StandaloneExport.*` | `ProcessedAudioObserver` supports bounded pre-monitor capture; monitor/routing and overload protection are separate. Final artifact measurement and strict publication gates remain independent. |

Existing plugin/Core reuse candidates:

- `Core/Processing/DSP/LookaheadLimiter.*`: prepare/reset/process, lookahead, attack/release,
  drive, detector-only 2×/4× FIR oversampling, effective latency query and block gain meters.
  The header retains a sample-peak-v0 description and implementation uses an oversampled
  detector. Neither the name nor detector proves a qualified true-peak output ceiling.
  Review the actual algorithm, overshoot, delay/tail, channel linking and seek/state first.
- `BroadbandDownwardCompressor.*`, `BroadbandUpwardCompressor.*`,
  `SpectralUpwardCompressor.*`, `UpwardGainLaw.h`, `BS1770MomentaryLufsDetector.h`:
  transfer/detector/smoothing and activity starting points, not qualified standalone
  methods. Audit detector units/weights and LFE policies; legacy defaults must not silently
  replace programme channel policy. Spectral processing introduces latency/state issues.
- `Core/Processing/Modules/MultiThresholdDynamicsModule.*` and `LevelerModule.*`:
  plugin parameter binding, per-stage atomic current/peak/hold meters, latency and stage
  wiring are useful examples. APVTS/DAW control, live-history assumptions and plugin
  parameter/state compatibility stay separate. Some old skeleton comments are stale;
  inspect actual implementations rather than infer capability from headers.
- `Core/Processing/IAudioModule.h`, `ProcessContext.h`, `ProcessorCore.h`:
  lifecycle/RT module-host material already exists. Its mutable module graph is not the
  qualified immutable standalone processing plan; adopting/refactoring a general graph
  is not a prerequisite or authorization for the next slice.
- `VolumeHistoryComponent` family and `DynamicsCurveComponent.*`: reusable activity/curve
  presentation ideas under UI ownership. Block peak/hold summaries cannot supply exact
  event counts, affected duration or a time/channel trace by themselves.

A future planner needs a processor capability/version contract and full settings:
input/output domains and units; channel detector/link/application masks including LFE;
stage order and bypass; gain reduction/boost over source time and channels; activity
threshold/resolution, extrema/percentiles, affected union, longest sustained activity
and event/repetition definition. Label estimated versus measured activity. Collect it
boundedly without affecting DSP; missing activity stays unavailable.

It also needs latency at each rate/setting, detector delay, tail/flush and exact output
length, state reset, required preroll, deterministic seek/block behaviour, finite/overload
policy, cancellation/resource bounds and parity across offline/audition/export. Distinguish
processor latency from authoritative audio-source coordinates and monitor/device delay.
Qualify compensating trim/preroll/EOF behaviour before preserving the existing duration
and video contract can be claimed. Keep final dither outside the shared float processing.

## Design questions for the new chat

1. Among candidates passing all selected objectives, what should be preserved: dramatic
   contrast, dialogue perspective, ambience/noise floor, music phrasing, transients or
   spatial balance? Which judgments are explicit protections, soft preferences or merely
   annotations? Which examples demonstrate the difference?
2. When is broad scene/region gain less intrusive than selective compression or limiting?
   Compare equal delivered requirements from the same frozen baseline. Report local and
   global activity separately; distinguish delivery-level audition from any explicitly
   level-matched comparison so a louder candidate is not automatically preferred.
3. How should amount, total duration, longest run, repetition/density, source location,
   transitions and effects on satisfactory neighbours influence choice? Define measurable
   activity before assigning costs. Do not conflate rank score with audible damage.
4. Which M19 limitation should be addressed first: whole-region greedy selection, policy
   budgets, fixed transition geometry, neighbour coordination, or residual peaks requiring
   a new qualified processor? A limiter can help peak constraints but cannot be assumed
   to solve the cinema excerpt's LRA deficit.
5. When is one explanation enough, and when should the user explicitly request alternatives?
   Bound candidate count/pass/runtime/disk costs. Compare from the same measured Current,
   not successive applied outputs. Revisit ranking among failed attempts separately from
   exact all-rule pass gates; never hide failed measurements behind a preferred score.
6. Define future freeze/protection precisely: fixed settings, no additional local change,
   or preserved audible output? Does it constrain global gain, cumulative prior changes,
   linked channels, detector/lookahead/transition support or neighbouring regions? What
   tolerance and priority applies when protections conflict with delivery rules? Unsupported
   mandatory protection must not be advertised as enforced.
7. Which user control maps to which measured quantity, units, scope, baseline and hard/soft
   rule? How are unset, zero, unavailable and unsupported distinguished? How do method/
   policy/protection changes invalidate proposals while compatible audio measurements remain?

## Evaluation collection and first assignment

Start with the preserved success, gate-crossing and cinema cases above, then propose
a small varied collection before any new rendering/tests. Include stereo/5.1 quiet
ambience/dialogue, sustained loud music/effects and brief peaks; keep independent examples
for checking generalisation. For each: content/member/role/range identity, baseline,
requirements/policy, concrete mixing concern, acceptable and unacceptable outcomes,
activity to inspect, measurement gates, listening context and user judgment. A dialogue
example is not authorization for speech detection, stem processing or intelligibility scoring.

**First assignment:** review these documents, M19 evidence and actual code; deliver a
short decision matrix comparing the next useful slices, then recommend one bounded
design/implementation sequence. Include the processor and first proposal method together,
manual/test controls, settings/activity needed, UI ownership, failure examples, acceptance
tests, dependencies, resource costs and explicit exclusions. Explain why the chosen slice
is preferable now and what would change that decision. Leave the result at proposal stage.

A qualified true-peak limiter and one section-gain-plus-limiter method are promising,
but not designated the next milestone. Improved region/gain/transition/neighbour work
must be considered against the actual M19 cases. Do not build all DSP first, silently
expand M19, edit production code, run new tests/DSP experiments, access physical audio
or start qualification during this design assignment.

Future implementation needs an impact-selected plan before editing: targeted processor
and shared-caller tests, same frozen candidate on relevant actual platforms, exact
offline/audition/export parity, state/seek/EOF and latency/timing checks, bounded resources
and coordinated listening. Reuse unchanged evidence with exact identity/rationale;
no recursive historical suite or automatic repeated long cinema workflow. Latency,
buffering/clock/macOS dependency changes warrant focused timing/video tests rather than
assuming M19's video reuse still applies. Physical opening/audible use needs coordination.
Occasional picture skips remain accepted with counters/timing retained; audio interruption,
accumulating drift, sustained freezes and identity/lifecycle failures remain blockers.

## Independent performance/UI work

These are separate work items, not strategy prerequisites or new approval:

- Remembered device selection, opening on explicit selection and optional reconnection
  are separate planned device/UI work. Define the user's connection intent and lifecycle
  before implementation; none of these behaviours changes in M19 publication.
- Clearly show initial import/verified-source preparation and progressive analysis states,
  progress/cancellation and pending versus authoritative results.
- Design session-wide reuse of a verified source across analysis, proposals, audition
  and export where identities permit. M19 already shares one copy within its private
  passes; it does not supply session-wide ownership, sharing or persistent reuse.
  Define exact member/assembly/composition identity, reference lifetimes, concurrent
  jobs, replacement/cleanup and disk/write failures without weakening verification.
- Consider optional project-owned audio with user-visible ownership/relink/cleanup and
  disk estimate. PCM24/48 kHz is about 0.97 GiB/hour stereo or 2.90 GiB/hour 5.1;
  three-hour 5.1 about 8.69 GiB. Float32 5.1 is about 3.86 GiB/hour. Add actual
  members/components, headers, staging/candidate/export working space; avoid duplicate
  copies. Reopening validates content/assembly and cache compatibility, never filenames
  or timestamps alone. This is not approved implementation of “copy audio into project”.
- Profile the inherited Intel completion pause independently; the Intel optimisation
  accepted its approximately 0.8-second legacy pause, not a fix. Preserve its historical
  strict-threshold results. M19's bounded inspection responsiveness figures are not proof
  the inherited pause disappeared.
- The user deferred the separate first-playback feedback task while shared-source reuse
  is planned. Do not revive it incidentally in this strategy chat.

Preserve installations, Pro40 code/drivers/services/system audio, original media, all
checkpoints, both M19 packages and unrelated UI/video work. Full physical multichannel
routing, immersive expansion, VST3 parity and new dialogue objectives stay separate.

## Pasteable starter message

> Please begin the dedicated **Proposal Strategy and Audio Quality** design chat.
> Read `/Users/fuzz/Developer/LevelScope/Docs/ProjectIndex/69_ProposalStrategyAndAudioQualityHandover.md`
> and its reading order. M19 candidate-02 is qualified within its approved scope:
> tested product `cf275afb133118dd2a06e7bca59f82f98898b3cd`, qualification/review commit
> `33e2248cb119290aaa7f16e03bb31faffbbbd0f0`, branch `codex/m19-lra-passages`.
> Use the subsequent documentation-only handover revision on that branch; do not
> treat it as a new tested binary. Main remains M18 `51e47f300fcc38a386a5a92a7277029175cfab0c`;
> publication/main integration of M19 remains a separate decision.
>
> Review the actual M19 success, policy-limited cinema and gate-crossing results.
> Develop each processing capability alongside its first useful proposal method in
> bounded successive slices. Keep section-led, dynamics-led and mixed strategies,
> exact processors, delivery rules, preservation preferences and measured/perceptual
> assessment distinct. Logic owns implementation/qualification; UI owns controls and
> activity/timeline displays; this chat owns evidence-led strategy rules/trade-offs
> recorded in the shared proposal specification.
> Evidence-led includes analytical/predictive models, surrounding/overlapping-window
> and transition/gating calculations, bounded simulation/search and neighbour-aware
> design. State assumptions, exact versus approximate/bounded predictions and required
> data; calibrate against measured examples. Complete candidate, independent Current
> and finalized-artifact measurements remain acceptance authority.
>
> Propose a bounded next design/implementation sequence with alternatives, evaluation
> cases, activity information, success/failure criteria, costs and ownership. Compare
> qualified limiting plus a first section-gain-and-limiter method against better
> region/gain/transition/neighbour rules; no option is already approved. Preserve exact
> identity, explicit Apply, independent Current and finalized-artifact verification.
> Keep preparation/source-reuse/completion-pause work separate. Begin with design and
> analysis only: no production edits, new tests/DSP work, devices, publication or merge.
> Preserve media/packages/tags/installations and stop with a proposal for review.

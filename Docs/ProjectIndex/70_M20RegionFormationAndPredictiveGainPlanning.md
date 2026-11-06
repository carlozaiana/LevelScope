# M20 region formation and predictive gain planning

6 November 2026. **Approved for Logic implementation and focused non-audible checks
following advisory review.** UI implementation follows a committed, focused-tested
Logic interface/core checkpoint; full combined qualification follows UI integration.
M20 publication and merge into main are not authorized. Include bounded splitting
and boundary refinement in the first useful
M20 slice, together with predictive gain planning. Reuse the qualified
`section-gain-v1` processor, its fixed transitions and the shared measurement,
Apply/history and export contracts. Better selection of unchanged M19 regions alone
would retain the geometry limitation this milestone is intended to address.

Implementation checkout: `/Users/fuzz/Developer/LevelScope`,
`codex/m20-region-planning`, based on published main/origin/main
`22961f03de7919e1de631f8b13eb4b852108ac0b`. This includes all three qualified UI
convenience slices, checkpoint `ui-convenience-2026-11-06` at
`b7350128d7b65507a5db434e38ecd73b8f2c4cfd`; read the
[published UI handover](../UIConveniencePublication-2026-11-06.md) and preserve its
device, preparation and monitoring behavior. The original design was prepared on
`codex/m19-lra-passages` at `5fc8c6ee95e705077cf996cecbc6543dbc20d257`;
that branch and all tested packages/checkpoints remain preserved.
Qualified M19 checkpoint: `a14c88e8a94412eb2e2536860aa7cbd79de0dd39`; tested
candidate-02: `cf275afb133118dd2a06e7bca59f82f98898b3cd`. Read
[69](69_ProposalStrategyAndAudioQualityHandover.md) for preserved identities and
ownership, [66](66_M19LraPassagesAndMeasuredCandidates.md) for the unchanged M19
method, and [the shared specification](../ProposalEngineSpec_v2.md) for authority.
This proposal supersedes the earlier chat recommendation to defer all splitting.

## Approval and implementation handover

The user authorized review and implementation dispatch in advisory chat
`01a13c8f-90a9-7150-813f-cb057626ef7e`: “The Strategy chat refined its proposal.
Please check it and if it is ok, tell the chat to send the necessary messages to
the logic module chat to start implementation.” The advisory review accepted this
design on 6 November and directed Strategy to commit the handover on updated main
before dispatching Logic. This supersedes the earlier design-only execution boundary
for M20's approved scope; it does not authorize new DSP or milestone publication.

Logic Module chat `01a13fdd-ecc3-7090-a004-bc921e58124f` on host `local` owns
constraints/accounting, formation/prediction, candidate measurement, selection/ranking,
session persistence and explanation contracts. Strategy owns rule-design review and
evaluation trade-offs. UI owns advanced controls, processing-region/transition
overlays and readable explanations, after the stable Logic handover. Strategy finishes
repository mutations before dispatch; Logic then owns this shared checkout.

Logic's first deliverable is a focused impact/test plan and concrete immutable
interface/identity contract, recorded before production edits. Proceed through the
owned core/service/session implementation and focused checks to a committed, tested,
integrated core checkpoint. Then prepare a concrete UI handover with exact commit,
API/data examples, callback/invalidation rules, explanation/activity fields,
session changes, test commands/results and remaining combined qualification.
Stop there when UI integration is the next dependency; final milestone qualification
is not a prerequisite to handing UI a stable interface. Do not start competing UI
edits or dispatch UI implementation now.

Initial controls, budgets and the 0.01 LU diagnostic shortlist remain evaluation
hypotheses. Freeze evaluated ordering/settings with provenance; retain the
heterogeneous-root usefulness gate and independent real-material review. Existing
media may be read without modification and non-audible offline implementation
checks may proceed. Audible/device use requires prior coordination with the user.

## Decision and scope

| Option | Useful outcome | Limitation and cost | M20 decision |
| --- | --- | --- | --- |
| Select different whole M19 regions | Improves allocation among existing groups | Cannot spare a satisfactory interior or give different parts different gains; low geometry cost | Retain as a comparison seed, insufficient on its own |
| Automatic region formation without markers | Subdivide measured high-population groupings and refine processing edges using correction-need changes | Can mistake a musical or dramatic change for an intervention opportunity; medium analysis/search cost | Include, with frozen analysis seeds and explicit search truncation |
| Refine inside one selected section | Engineer supplies context; find smaller useful intervals without editing the section | Strict containment may prevent a transition or exclude useful context from processing; small incremental scope/UI cost | Include, strict local area by default |
| Refine with an explicitly permitted surrounding margin | Move a correction's onset/return into a suitable nearby location | Changes neighbouring material and consumes duration; needs clear authorization geometry | Include through the same engine and advanced numeric controls; zero margin by default |
| Adjustable ramp shape or coordinated attenuation of additional neighbours | Could improve transitions beyond what moving existing fixed ramps achieves | New processing/version/qualification and more collateral change | Defer; M20 coordinates boundary choices and existing ramps only |
| Qualified limiter plus first section-gain-and-limiter method | Can address brief peaks obstructing otherwise acceptable delivery | New state, latency, output-ceiling and activity qualification; not the demonstrated cinema LRA remedy | Separate later capability/method decision |

The first implementation supports one selected section or temporary selection at
a time, or automatic programme search. It does not combine overlapping selected
sections. Selected-section mode reuses the same formation/planning engine with
different constraints; the margin is an explicit expansion of its permitted area,
not a separate DSP capability. Multiple selected sections and a general region
editor can follow without delaying this useful slice.

Eligibility remains M19's stable measured LRA, selected usable IL/TP/LRA rules,
all additional selected M/S rules, matching measured Source/Current and confirmed
stereo/conventional 5.1 at 44.1/48/96/192 kHz. Processing includes every programme
channel, including LFE. Resolve pending manual drafts first. All selected rules
already passing remains a genuine no-op. M19 and the LRA-unset workflow retain
their existing versions; M20 is an explicitly selected experimental method.

## Terminology and independent objects

| Term | Meaning |
| --- | --- |
| Section | User-defined or imported organisational span. Its annotation, name, bounds and provenance are preserved. |
| Analysis region | Measured grouping with a method/version and source-time support. It is neither a semantic scene nor proof that correction is needed. |
| Processing region | Interval proposed or used for a correction. One section or analysis region may contain several processing regions. |
| Planning root | Internal provenance object tying all descendants to the original selected section/selection or frozen automatic analysis seed. It is not a new user section. |
| Permitted local area U | Union of half-open source-frame intervals in which additional local-envelope change is allowed. It includes any explicitly granted margin. |
| Changed support D | Canonical support of the difference between frozen Current's local envelope and the complete proposed local envelope, including transitions under existing v2 accounting. |

Use these terms in new design text and new explanations. Existing API names such
as `SectionGainPlan`, strategy ID `section-led`, saved schemas and historical
reports remain intact. In particular, the old specification's “processing section”
corresponds to a processing region; no broad rename or migration follows from the
terminology correction.

## Selection intent and boundary permissions

Store the selection hint H separately from U. Reading or predicting surrounding
audio never authorizes processing it. Every candidate must satisfy `D subset U`
and must have no overlap between D and protected local material P. Validate the
complete canonical before/after envelopes, not just the proposed plateau rows.

| Mode | H and U | Boundary behavior |
| --- | --- | --- |
| Automatic | No user hint; U is the programme minus explicit local protections | Frozen M19 analysis regions are initial roots, not compulsory correction spans. First-version search refines inward within each seed; it does not extend seeds or join separate roots. |
| Inside selected section | H is the original section; U is exactly its bounds | Find shorter intervals and splits inside U. All transition support must fit; nothing outside is silently included. |
| Selection with permitted margin | H retains the original section; U is `[H.begin-leftMargin, H.end+rightMargin)`, intersected with programme bounds | Outer edges can move inside the explicitly displayed U. A changing candidate must intersect H and address a need identified in H; the margin is not a search for unrelated nearby corrections. |

A rough selection means its endpoints are hints, not mandatory correction edges.
It still has a finite, explicit U. Never interpret “hint” as permission to process
the whole programme. “Strict area” limits where processing may occur; it does not
require all audio within the area to change. Exact-boundary pins are separate.
Changing the requested margin updates a draft U; assessment uses only confirmed
values. At EOF/start, show requested and effective clipped bounds before assessment.

The UI's **processing start/end** are the outer edges of D, including ramps, in
source frames. Store a closed permitted range for each edge, intersected with U
and ordering/minimum-duration constraints. The advanced controls initially offer
the full legal inward range; numeric edits can narrow it. A pin sets min=max.
Do not label a plateau boundary as the start of all audible processing. Display
the actual plateau cuts/ramp support separately in the inspector.

If “Allow boundary refinement” is off, pin the outer changed-support edges to the
original selected section endpoints (or automatic seed endpoints). A granted
margin alone does not move these pins. Enabling refinement permits shortening or
movement only within the two stated ranges. Pinning the outer edges does not prohibit an allowed
internal split. Internal proposed cuts each retain their own permitted range;
the first UI needs a numeric range/pin in a row inspector, not draggable automation.
Reassessment starts from frozen Current, never from an unaccepted preview.

Pins constrain a changing root; they do not require a correction or defeat the
all-rule no-op. An internal pinned divider may remain in explanatory metadata
when equal corrections canonicalize to identical audio; report that it has no
audible boundary. A mandatory unchanged gap is represented as protected material,
not as a marker or an equal-gain divider.

An exact support edge may require moving a plateau cut inward by the fixed ramp
radius. Construct a legal `section-gain-v1` plan and validate its resulting support;
do not assume the radius is always 50 ms or round the user's pin to the search grid.
If exact edges cannot coexist with existing geometry, reject that candidate with
a method limitation. No ramp truncation, new curve or hidden pin movement is allowed.

## Proposed advanced controls

These initial values are hypotheses for evaluation, not EBU requirements or claims
of transparency. Store seconds entered in integer milliseconds, convert outward
permissions inward to source frames (lower bounds ceil, upper bounds floor), and
retain explicit frame inputs exactly. Reject empty/contradictory ranges. Apply
gain and requirement comparisons at their existing full precision.

| Control | Initial value and range | Meaning and reason |
| --- | --- | --- |
| Allow splitting | On in the explicitly chosen experimental method | Enables different correction decisions within a broad root. Off requires no new interior correction cuts. |
| Maximum splits per original selected section | 2; integer 0–4 | Counts new internal correction boundaries, not resulting regions. Two permits three parts, including a spared middle. Automatic mode labels the same limit “per analysis seed”. |
| Allow boundary refinement | On | Needed to evaluate approximate selections; pins/ranges remain authoritative. |
| Left/right surrounding margin | 0 s; independently 0–3 s in 0.1 s steps | Outward movement is opt-in. The initial maximum is one short-term window, a bounded context allowance rather than a preferred transition length. |
| Minimum processing interval | 3 s; 0.5–10 s in 0.1 s steps | Minimum length of a newly assigned nonzero-correction interval, between its nominal cuts, before fixed ramps. Starts at one LRA window; shorter choices are experimental, not statistically independent LRA measurements. |
| Minimum unchanged separation | 1 s; 0–3 s in 0.1 s steps | Minimum gap between distinct components of D after ramps. Avoids rapid repeated intervention; does not require a gap between adjacent different-gain intervals. Zero still cannot erase protected material. |
| Consider filling a short gap | Off | If enabled, a merged correction is a separately predicted candidate. Never bridge protection, a strict area gap, another root or an exact separator. |
| Need-change threshold / persistence | 1 LU / 1 s; 0.5–3 LU / 0.3–3 s | Candidate-boundary discovery controls. Suppress tiny fluctuating changes without claiming scene detection. Both use full-resolution measured windows. |
| Allow global trim of whole programme | On in automatic mode; off initially in either selected-section mode | An explicit independent permission. With it off, global gain equals frozen Current. With it on, show that audio outside U can change globally. |

Keep M16's gain/reserve/duration/count controls in their existing policy panel:
default additional attenuation 3 dB, cumulative 12 dB, affected duration 20%, no
user count restriction, reserve 0.5 dB inside caps. No new local upward correction.
The split limit, changed-region count, 10,000-interval engine capacity, affected
duration and search budget are independent controls/limits.

Evaluate min-duration values 1/3/6 s, split counts 0/1/2/4, margins 0/0.5/1.5/3 s,
and threshold/persistence variants on selected representative cases, not their
Cartesian product. Freeze choices before two independent real holdouts. Change
defaults only for a documented improvement in usefulness and preservation, not to
make the cinema example pass.

## Split accounting and geometry enforcement

Count splits in the final correction partition of each original planning root,
after merging truly redundant correction assignments. Outer start/end refinement
does not count as an internal split. A new distinction in additional correction
across an existing applied boundary does count; it is not free because a sample
coordinate already existed. Existing baseline boundaries alone do not count.

Examples relative to frozen Current:

- `-2 / -1 / -2 dB` inside the refined outer edges uses two splits. It may be one
  connected changed region for M16 accounting.
- `-2 / unchanged / -2 dB` uses two splits and normally two changed regions.
  The unchanged gap is retained in the partition; it cannot be omitted to undercount.
- `-2 / -2 dB` with an annotation-only subdivision uses no split when canonical
  processing is identical. Keep the annotation/provenance; do not create new DSP.
- Shortening both outer edges of one correction uses no internal split but must
  obey boundary permissions, minimum interval length and complete support accounting.

Thus K splits permits at most K+1 decision parts within the refined outer span,
including unchanged internal gaps. The minimum processing-interval duration applies
to nonzero assignments; unchanged interior gaps obey the separation rule instead.
All descendants, including margin portions, retain the same original root ID and
share its cap. Repeated split operations cannot reset it. Exploratory alternatives
consume search work, not the final split allowance. Repeated explicit applied
proposals still obey cumulative gain limits; the split cap is not lifetime history.

Do not merge just because two labels or desired deltas are equal if removing the
cut changes an existing ramp. Compare canonical resulting envelopes; baseline cuts
needed to reproduce Current remain. Exact linear-segment endpoint checks establish
additional gain bounds and zero change outside U/P. `changedEnvelopeSupport` defines
the conservative charged union, including its endpoint convention. Charge all of
D once across the whole programme, with the original programme length as denominator;
never use the selected section or child duration as denominator. Count connected
components separately from split count and resulting plan interval count.

Touching identical corrections canonicalize into one. Adjacent unequal corrections
may remain with the existing ramp and consume a split. Separate corrections whose
post-ramp unchanged gap is too short are rejected, or compete against a permitted
merged candidate. A merge pays for the entire filled gap and is rechecked/predicted;
it is not a free smoothing operation. Never merge across protected material.

Imported markers are organisational context, not implicit pins or split candidates.
Redundant marker changes cannot alter automatic output. Selecting an imported
section makes its original bounds/provenance the root; making a marker exact is a
separate explicit pin. Pins/section bounds never rewrite annotations.

Existing gain transitions belong to frozen Current. Candidate geometry must account
for both old and new ramps; a newly inserted cut can shorten a neighbouring ramp.
Reject any resulting extra boost, cap violation, unauthorized support or incompatible
pin rather than moving/flattening old processing. Show such rejections as method
limitations, not proof that every possible envelope is impossible.

## Protection supported by M20

Support mandatory **no additional local-envelope change** on explicit source-frame
intervals, for all programme channels. Global-trim permission remains separate.
With frozen global gain and the current stateless gain chain, identical local
processing outside U must preserve pre-monitor samples there; qualification checks
that parity. This narrow guarantee is not a generic whole-chain freeze facility.

Locking named settings and preserving rendered audio under arbitrary global,
neighbouring or stateful processing remain different semantics. A mandatory
unsupported protection blocks changing generation; do not substitute local
protection. Nonzero global trim cannot be offered as preserving a protected audible
result. If delivery rules and supported protections conflict, show the failed rules
and constraint conflict, leaving Apply unavailable. Protection edits themselves
do not alter audio or invalidate otherwise compatible baseline measurements.

## Required analysis and prediction

Existing `ProgramAnalysis` supplies full-resolution window bounds, energies,
inclusion/padding/status, gates, percentiles, ordered peaks and identities. Add an
immutable request-owned prediction sidecar: weighted K-filtered energy sums and
sample counts in 100 ms source-aligned bins, per-channel absolute sample maxima,
and exact baseline/envelope/role/method provenance. It is neither a display overview
nor a new authoritative measurement. Collect it in at most one auxiliary streaming
pass through the same verified source and exact frozen Current processing.
Measurement arithmetic/order must remain unchanged; use a non-mutating collector
or separately qualified auxiliary analysis. Do not claim the existing window-only
observer already supplies these bin energies or sample-resolution maxima.

Formation uses short-term correction-need differences, momentary context, gate
proximity, existing transitions and protected neighbourhoods. Define provisional
need at each available LRA window as `q=max(0, L-P10-selectedLraCeiling)`; use all
available windows for context, retaining included/excluded status. The original
P10 is a seed reference, not a fixed target during prediction. A window's centre
can nominate a boundary hypothesis; its actual three-second support still drives
energy overlap. Do not equate that centre with an audible event onset.

For a complete candidate, approximate transformed window energy by summing baseline
bin energies weighted by the candidate-to-baseline squared gain ratio, integrated
over each bin/window intersection. Use the actual piecewise dB ramps and global
gain, not one gain for every window touching a region. Uniform energy within a bin,
gain/filter noncommutation, boundary memory and float-stage rounding are model
limitations. Constant-gain scaling is exact only in its declared ideal linear model.
Baseline bins do not reconstruct arbitrary new filtered samples exactly.

Anchor each forecast to the measured baseline window energy and add the modelled
energy delta; the identical-plan proxy must reproduce the measured baseline rather
than accumulate a fresh bin-summation error. Retain separate EOF-drain bin/state
provenance needed for LRA; it never extends processing or duration accounting.
Changes near EOF make the transformed filter tail uncertain and must be flagged,
not modelled as a fabricated measured silence or a precisely scaled tail.

Recompute absolute/relative gates and percentiles from the full predicted programme,
including formerly excluded windows and EOF padding. Recompute predicted IL, LRA,
P10/P95, membership and M/S curves jointly. No “local LRA contribution” is assigned
to a region. Gate-sensitive predictions retain that warning even when proxy LRA
improves. Reading context extends beyond U, including window support and filter
history; it never widens D. A finite three-second halo is not an exact IIR-state bound.

TP estimates require all channels and interpolation support. A conservative bound
may use the qualified interpolation coefficients and upper sample-magnitude/gain
bounds through `abs(sum(h*x)) <= sum(abs(h)*abs(x))`, with rounding handled. Such a
bound can be loose: an upper bound above the ceiling is inconclusive, not evidence
that the actual candidate fails. M/S curves sampled at 100 ms cannot certify
sample-resolution maxima. Mark unavailable/uncertain predictions explicitly and
let complete measurement decide. No limiter prediction is included in M20.

Calibrate prediction errors against actual candidates by metric, gate-membership
changes, plateau interiors and boundaries. Retain maximum and signed errors and
material ranking reversals. Empirical uncertainty is not a proven bound. A model
miss may change which pre-budgeted candidate is measured next, but cannot relax a
requirement or silently authorize a fourth assessment.

## Deterministic bounded formation and planning

Proposed independent IDs: `region-formation-v1`, `region-gain-search-v1`,
`region-constraints-v1`, `window-energy-proxy-v1`, `candidate-review-rank-v2`.
They do not rename `lra-passages-v1`, `lra-selective-trim-v1`, the processor, delivery
requirements or M16 policy. Persist exact effective settings, not just these names.

1. Validate identities, capabilities, scope, pins/protections and controls. Return
   genuine all-rule no-op first. Consider sufficient global-only processing only
   if explicitly permitted; it still needs complete candidate measurement.
2. Freeze roots once. Selected mode has one original root, including its granted
   margin. Automatic mode takes at most eight M19 high analysis seeds: alternate
   the existing excess-score order and excess-score-per-support-second order,
   skipping duplicates, with window count then earliest frame as tie-breaks.
   Report how many roots were not explored. This cap is a search budget, not a
   changed-region limit. Empty automatic seeds yield a specific method limitation.
3. Nominate internal cuts on the existing 100 ms source grid. At each location,
   compare median q on the preceding/following three seconds, without zero-filling
   unavailable data. Require each side's persistence interval to remain on its
   corresponding side of the two medians' midpoint and their difference to meet
   the configured need threshold. Rank by absolute median difference, then earlier
   frame. For unchanged-gap candidates also nominate the ends of runs with
   `q <= needThreshold/2` for at least the persistence duration. Rank these by the
   neighbouring median contrast, then earlier frame. Expand each nomination in
   order to its centre, minus 1.5 s, plus 1.5 s; deduplicate and discard out-of-range
   coordinates, retaining up to 12 nonmandatory locations per root. This explores
   window smearing without claiming event detection. Preserve
   exact off-grid pins separately; pins must fit the split cap and legal geometry.
4. Form at most 16 geometry alternatives per root: legal unsplit original first;
   unsplit outer refinement second, using the highest combined onset/return
   nomination strength. Then alternate the allowed split counts 1 through K,
   considering their combinations by decreasing summed cut strength and
   lexicographic source-frame order. Reserve a legal two-cut unchanged-middle
   alternative when available; one-cut variations must not consume all slots.
   Inspect at most 256 topology/edge combinations per root, round robin across
   split counts, retaining up to 16 legal distinct geometries and reporting
   truncation. Outer edges are drawn from root limits, edge-range endpoints,
   pins and nominated cut locations. Rank equal-cut alternatives by total edge movement then earlier
   edges. For each proposed support edge, map to a legal fixed-ramp plateau cut and
   verify the resulting exact support; skip unrepresentable hypotheses with reasons.
   Keep minimum length, separation, caps, pins and provenance throughout. No silent
   extension to fit a minimum duration. Never infer that omitted geometry is infeasible.
5. For each decision part, let h be the least permitted additional/cumulative
   headroom across its baseline. Use included high windows whose centres fall in
   the part only to seed gains; full-support overlap still determines the proxy.
   Let r be their provisional P95 overshoot plus reserve, capped to h. Without
   positive need the default choice is zero only; a separately permitted merged-gap
   candidate derives a common gain from its combined span. Candidate reductions
   are the unique feasible values `{0, r/2, r, min(h,r+0.5 dB), h}`. All numbers stay
   full precision, with inward representable adjustment at caps. Zero keeps an
   interior part unchanged and pays the corresponding split/separation accounting.
   Reserve stays inside caps and requested/achieved amounts remain visible.
6. Explore geometry/gain combinations with a deterministic beam of at most four
   complete plans. Start with the legal M19 recipe and legal unsplit roots; nominate
   one-root add/drop/replace, one-cut add/remove, boundary move and shallower/deeper
   changes. Visit roots in frozen order and alternatives in the order above; round
   robin root expansions so a long first root cannot consume every evaluation.
   Cache by exact canonical audio plan identity. Validate all hard geometry/gain
   constraints before proxy evaluation; identical audio has one evaluation.
7. Include global trim in proxy plans if allowed. Candidate global choices include
   frozen gain, predicted IL-target correction and predicted feasible interval
   endpoints; clamp only to declared parameter limits, not to a relaxed delivery
   threshold. Recompute populations for each distinct gain; do not assume IL shifts
   linearly through a gate crossing. Uncertain TP/M/S bounds cannot certify feasibility.
   Every distinct complete plan/global choice consumes the proxy budget.
8. Spend at most 32 proxy evaluations before complete assessment 1, then up to 16
   after its feedback and 16 after assessment 2. Rank proxies with predicted rule
   deficits, the failed-review rule below and preservation activity; uncertainty
   accompanies the rank and never becomes measured authority. Keep a distinct
   shallower or different-geometry challenger when one exists. Measure the most
   promising complete plans, at most three; each already includes its chosen global
   gain. A final global adjustment consumes an assessment, never an unmeasured patch.
   Stop at a measured passing recommendation by default; explicit comparison may
   use the remaining budget to assess a materially different alternative.

Unknown TP or sample-resolution M/S predictions do not disqualify a plan from
complete measurement. Proxy ranking uses available predicted deficits and displays
unknowns; the missing-required-metric rejection below applies to actual baseline/
candidate measurements. In particular, do not turn a loose TP upper bound into
a false measured failure or discard every region because such a bound is above
the ceiling. If refined measured candidates fail, retain a legal unsplit/M19-recipe
challenger for the remaining assessment instead of automatically deepening again.

The finite proposal ordering is part of the method version. If implementation
changes discovery, tie-breaking or enumeration to improve coverage, record the
change and affected checks. All unexamined roots/geometries, unavailable predictions
and exhausted budgets are disclosed. None establishes mathematical infeasibility.

Every proposed processing region explains its root, baseline/window need, why its
edges/cuts were nominated, the effect predicted for neighbouring windows, predicted
gate changes, intended gain, caps/protection and actual complete-plan result when
measured. Say “this combined plan changed LRA by ...”; do not claim that its regions'
effects add independently without evidence.

## Comparing measured results

Apply eligibility is binary and exact: complete available metrics, every selected
rule passing at full precision, supported protections and valid immutable binding.
All-pass beats any failure. Among passing candidates prefer a non-dominated choice
over affected duration, attenuation magnitude/integral, longest changed run, split
and transition complexity, neighbour activity and global change. If neither
dominates, show the trade-off; use a documented default lexicographic order of
affected frames, integrated absolute local dB-change over frames, maximum added
attenuation, split count, absolute global delta and stable candidate ID for the
initial recommendation. This engineering order is not an artistic-quality score.
Do not reward lower LRA once its requirement passes.

For a failed-result representative, first exclude missing-required-metric attempts
from numerical ranking. Let e_min be the smallest actual nonnegative LRA excess
among complete measured failures. Form an anchored shortlist `e <= e_min + 0.01 LU`.
Within it prefer fewer failed selected rules; then smaller maximum and summed other
positive violations normalized by explicitly recorded 1 LU (IL/M/S) and 1 dB (TP)
reference scales; then preservation order; then exact LRA excess and stable ID.
This avoids a non-transitive pairwise “approximately equal” comparator. The 0.01 LU
band is a proposed failure-ranking resolution, never an acceptance tolerance or
permission to round measurements. Evaluate 0/0.001/0.01/0.05 LU on failure examples
before fixing v2. Retain every measured failure for inspection, including trade-offs
outside this headline shortlist. A user may prioritize another failure diagnostically;
no failed choice can Apply.

The M19 gate fixture makes the reason concrete: LRA 6.999999882961 with failed IL
currently outranks 6.999999887254 with passing IL. Both fail the 6.7 LU ceiling.
The proposed ranking presents the latter as the better failed compromise while
retaining the former and the 25.989999164220 LU gate-crossing failure. All-rule
thresholds, numerical reports and historical M19 decisions remain unchanged.

## UI and persistence responsibilities

UI/Timeline owns an advanced Region planning panel attached to the existing proposal
workflow: mode and selected section, strict area and margins, edge ranges/pins,
split/minimum/separation controls, local protections, global-trim permission and
method/search settings. Logic owns validation and returns effective geometry,
reasons and rejected alternatives. A temporary selected time span need not create
an annotation. Provide numeric read/edit controls and an assessment action; defer
general drawing tools, freehand automation and semantic scene labels.

Show sections, analysis regions, proposed processing regions and complete transition
support as distinct layers. Display all counts independently: new splits per root,
decision parts, connected changed regions, plan interval capacity, affected seconds/
percentage, cumulative local support, whole-programme global delta and work consumed.
Show strict areas/protections and rejected out-of-area support. A metadata-only
division must not appear as an audible change or trigger a new assessment.
Inspection does not seek/play or switch audition domain. Private previews do not
become Current. Any future private-candidate audition needs its own explicit contract;
M20 does not silently add it. Listening qualification can use identified diagnostic
captures; product Apply and independent Current verification remain explicit.

Bind region requests/results to measured source/current identities, roles/range,
baseline plan/revision, selected original section identity/bounds, H/U/P, pins,
controls, global permission, strategy/method/model/rank versions, requirements,
policy and runtime job. Name-only marker edits refresh display; changed bounds of
a referenced section invalidate its request without moving applied audio. Unselected
marker edits do not alter automatic geometry. Metadata changes retain compatible
audio measurements but invalidate dependent proposals and approval metadata.

Persist settings/provenance and the last applied explanation through a versioned
extension coordinated with Logic's session owner (schema 9 from current schema 8
if still free at implementation). Old sessions retain their exact audio and M19
method until M20 is explicitly selected. Unknown mandatory constraints fail closed;
do not fall back to wider permissions. Candidate predictions/measurements never
restore Current verification. The sidecar is initially request-owned, not a new
session-wide source/cache ownership system.

## Resource budgets proposed for qualification

| Resource | Bound or acceptance target |
| --- | --- |
| User scope | One selected section/temporary selection, or automatic programme mode |
| Automatic root shortlist | At most 8; frozen before subdivision; disclose omitted roots |
| Internal discovery | 12 nonmandatory cut locations, at most 256 topology/edge constructions and 16 retained geometries per root; exact pins retained within split/structural limits |
| Search | Beam width 4; 64 distinct complete proxy plans total, including global variants, partitioned 32/16/16; no hidden Cartesian search |
| Private audio work | At most 1 auxiliary analysis pass and 3 complete candidate measurements, using 1 verified source copy; no persistent candidate PCM |
| Memory | At most 64 MiB additional region/model/search storage, separate from existing measurement/waveform storage; stream bins and reuse working arrays |
| Cached planning performance | Target <=10 s Sequoia / <=20 s actual High Sierra for three-hour-sized synthetic populations; measure before claiming support |
| Watchdog | Cached planning stops at 30/60 s respectively as resource-exhausted, with no partial Apply authority; host timing must not silently select different audio |
| Responsiveness | Cancellation <=1 s; UI heartbeat p95 <=100 ms; no planning/audio I/O on UI or audio callback |
| Retained explanation | Existing 8 MiB attachment and session/report bounds; no silent truncation presented as complete |

The sidecar pass is additional to M19 and must be timed separately. Existing
Source/Current acquisition, post-Apply independent Current measurement, export render
and finalized-artifact measurement are not hidden inside the three private passes.
Report their costs separately. Do not require session-wide verified-source reuse.
Do not store one full window population per proxy; retain summaries and at most the
bounded active working populations. On storage/capacity exhaustion report unavailable.
When automatic root/cut/geometry shortlists omit possibilities, a bounded result
is still possible but must state the explored subset; watchdog/cancellation ends
the request without Apply. No unbounded retries or automatic method fallback.

## Acceptance collection and impact selection

Record immutable content/member/role/range identities, requirements, frozen Current,
policy/permissions, exact plans, predictions/errors, actual measurements/activity
and separate listening context/judgments for each example. Use original M19 receipts
as evidence of the baseline, not as fresh M20 execution.

| Case | Required distinction and outcome |
| --- | --- |
| Existing stereo success | Preserve a simple unsplit candidate in the search; splitting must justify a preservation benefit. M19 achieved 5.910633862481 to 4.290033958906 LU with independent Current/artifact passes. No requirement that M20 choose identical gains, but it must not lose all passing choices on this fixture. |
| Existing gate-crossing fixture | Predict/reveal population entry/exit, compare shallower/deeper gains, retain the better measured attempts and exercise the new failed-result headline. Never accept 6.999999887254 against 6.7 LU. |
| Cinema 27:00–35:00 excerpt | Examine heterogeneous correction needs within the 183.0–274.4 s seed and alternatives to consuming 19.0625% at 3 dB. Evaluate less collateral change and useful prediction even if 15 LU remains unmet. Baseline 23.435332665255, M19 result 22.385113193953 LU; no full-film conclusion. |
| Synthetic broad root with two sustained high parts and a satisfactory middle | Two cuts can spare the middle; one cut cannot produce the same three-part decision. Demonstrate both adjacent unequal gains and two disjoint changed supports. Compare equal delivery requirements. |
| Independent stereo dialogue/ambience and sustained music/effects | Evaluate perspective, ambience continuity, phrasing and dramatic contrast. Keep at least one real example of each as holdout; no dialogue or intelligibility score. |
| Independent 5.1 scene with surrounds and LFE peaks | Enforce all-channel processing/TP, monitor independence and local protections; assess spatial and transient preservation with coordinated listening. No inference from partial L/R playback. |
| Boundary/permission fixtures | Exact off-grid pins, 0/positive margins, near EOF/start, shorter-than-minimum roots, gap shorter than separation, protected gaps, prior ramps, all split limits and one-frame-outside support. Reject impossible combinations without expansion. |
| Already compliant and unchanged-envelope cases | No-op preserves plan/revision/measurement; equal-gain splits and redundant markers neither change samples nor spend candidate assessments. |

The fixed-population ideal model bounds an attenuation-only window shift of at most
3 LU to at most 3 LU LRA reduction. For this cinema excerpt that leaves at least
20.435332665255 LU before considering the duration restriction. This conditional
bound explains why 15 LU is not an M20 success criterion. Gating/filtering can
invalidate the assumptions; it is not proof that the actual programme is impossible.

The usefulness gate is a passing heterogeneous-root example where subdivision
either finds a legal complete solution missed by the whole-region method or spares
more satisfactory material at the same selected requirements and policy. Retain
the simple stereo passing workflow. Independent real examples must receive separate
engineering/listening review of the proposed benefit and any collateral change;
do not declare preservation from the synthetic result alone. If these benefits do
not materialize, revisit formation/search rules before calling the slice complete.

Before production edits, Logic records an impact manifest selecting:

1. Pure permission/partition/accounting checks, independent sampled-envelope oracles
   and metamorphic marker/identical-audio tests; exact limits and existing-ramp cases.
2. Model checks against analytical constant-gain/overlap cases and actual complete
   measurements; gate-sensitive cases, prediction failures and deterministic budget
   termination, including no causal attribution from a selected window alone.
3. Short integrated stereo/5.1 requests at supported rates, existing assembly/
   composition callers, stale/cancel/reopen handling, strict all-rule Apply gates,
   unchanged Source, independent Current and finalized PCM24 verification.
4. All-channel offline/audition/export parity at new geometry/seek/EOF points; sample
   preservation outside strict local areas with global fixed; final-only dither.
5. Native advanced controls/overlays/reasons/history and bounded population/resource
   checks on the same frozen candidate on actual Sequoia and Intel High Sierra.
6. Coordinated listening of representative changed intervals, including surrounding
   context and proposed splits, with delivery-level and any explicitly level-matched
   comparisons labelled separately. Numerical success cannot replace this judgment.

No recursive historical suite, automatic full film/three-hour decode or qualification
by local Intel compilation. Reuse unchanged evidence with exact identity/rationale.
Any change reaching measurement arithmetic, buffering, latency, clock or dependencies
requires explicit impact reassessment and focused timing/video checks. Non-audible
offline implementation checks and bounded evaluation renders are now authorized
within the focused impact plan. Coordinate audible/device tests with the user
before opening hardware; preserve original media and frozen packages.

## Implementation sequence and ownership

1. **Contract and impact plan:** Logic records affected production/shared callers,
   immutable region constraints, support-edge meaning, identity and explanation
   fields, tests/reuse and resources before production edits. Strategy supports
   rule review. Freeze finite ordering and the initial evaluation collection with
   provenance. No broad UI/code rename.
2. **Logic core checkpoint:** Implement constraints, bounded formation including
   splits and inward/refined/margin boundaries, sidecar/model, gain search, measured
   selection/ranking, session persistence and explanation contracts using existing
   DSP. Integrate the core services, run focused checks and commit a stable interface.
   Prepare the concrete UI handover described above and stop at that dependency.
3. **UI integration:** UI uses that exact Logic checkpoint to add advanced numeric
   controls, separate proposed geometry/transition overlays and readable reasons.
   Stage shared-file work explicitly; do not edit those UI files concurrently.
   Together these stages form the first useful M20 slice; a whole-region-only
   planner is not completed M20.
4. **Combined qualification and preservation review:** After UI integration, qualify
   the same frozen combined candidate on the specified actual platforms, exercise
   targeted shared callers and the varied cases, and coordinate listening. Record
   prediction errors and user judgments before rule revisions. Earlier focused
   Logic checks are development evidence, not full milestone qualification.
5. **Following scope decisions:** Consider multiple selected sections, broader
   automatic boundary search, adjustable ramps/coordinated neighbour gains or a
   qualified limiter with its first method only when examples justify them.

M20 adds neither a limiter/compressor nor selective boosting, component/dialogue
processing, new delivery rules, immersive layouts, VST3 parity or general output
freeze. Preserve the published device, preparation and monitoring convenience slices;
source reuse and project-owned media remain separate work. No M20 publication,
checkpoint tag or merge into main is authorized by this implementation approval.

## Evidence and primary measurement references

- [M19 stereo proposal](../../Tests/Reference/results/m19-qualification-2026-11-05/sequoia-01/stereo/proposal.txt),
  [cinema proposal](../../Tests/Reference/results/m19-qualification-2026-11-05/sequoia-01/real/proposal.txt),
  [gate-crossing proposal](../../Tests/Reference/results/m19-qualification-2026-11-05/sequoia-01/gating/proposal.txt)
  and [qualification scope](68_M19QualificationAndHandover.md).
- [Measurement contract](14_OfflineMeasurementContract.md),
  [current preference accounting](56_M16ProposalPreferencesAndLimits.md),
  `Source/Core/Proposal/LraPassages.*`, `Source/Standalone/StandaloneLraProposal.cpp`,
  `Source/Core/Processing/SectionGainPlan.h`, `Source/Standalone/StandaloneProcessing.h`
  and `Source/Core/Proposal/Model/ProgramAnalysis.h` are the implementation anchors.
- [EBU Tech 3342 v4, November 2023](https://tech.ebu.ch/docs/tech/tech3342.pdf)
  defines the LRA window population, gates and percentiles;
  [ITU-R BS.1770-5](https://www.itu.int/rec/R-REC-BS.1770-5-202311-I/en)
  defines the loudness/true-peak measurement basis. Neither specifies these region
  boundaries, split limits, search budgets or preservation defaults; those are
  proposed LevelScope rules requiring the evaluation above.

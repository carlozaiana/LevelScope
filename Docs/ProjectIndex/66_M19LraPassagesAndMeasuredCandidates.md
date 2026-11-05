# M19 proposal — LRA passages and measured selective-gain candidates

5 November 2026. **Approved for implementation and focused qualification.**
The user approved the scope and the review clarifications via the standalone
coordination chat. Work starts on `codex/m19-lra-passages` from the published Intel
checkpoint plus publication receipt `09c8b8f6345c7c889f21f5399a38284ce1d353b2`.
M19 publication and main integration remain for later review.
Start after the qualified M18 implementation and review of the separate Intel
true-peak continuation. The latter improves measurement throughput and introduces
no LRA correction method. This document supersedes the old roadmap's pending-M18
status, not its historical test receipts.

## Useful first outcome

Explain which passages occupy the programme's low and high loudness populations,
with or without imported markers. Offer a bounded, explicitly requested trial
of attenuation of selected high passages. Show what its complete measured result
actually changes before offering Apply. An unsuccessful bounded attempt means
that this method/budget did not find an acceptable candidate, not that no solution
exists.

Plain EBU R128 still has no universal maximum LRA. “Excessive” means above an
explicitly selected job LRA ceiling with its provenance. Without that requirement,
the same passage view is descriptive and makes no correction recommendation.
Keep existing measured Source/Current status, LRA stability and all selected
IL/TP/M/S rules visible. No programme is described as passing because a preview
histogram looks narrower.

## Statistical basis and limits of interpretation

Reuse the qualified complete LRA-window population, shared gate predicates and
percentile implementation. LRA comes from the gated distribution of three-second
loudness windows, including its 10th and 95th percentiles; it is not an additive
quantity that can be allocated uniquely to individual scenes. A short peak alone
need not drive it. This interpretation follows [EBU Tech 3342, November 2023,
sections 3–3.1](https://tech.ebu.ch/docs/tech/tech3342.pdf).

Show gate membership, window counts, population share, low/high percentile
neighbourhoods and their exact source-time support. Label rows “passages to
investigate”, not defective scenes or proven causes. A high programme LRA can
reflect sustained quiet material as well as loud passages. This first method
does not raise quiet material or infer speech, dialogue importance or artistic
intent. Explain when those boundaries prevent a useful result.

Do not rank only the loudest five percent: changing isolated extremes can leave
the 95th percentile unchanged. The suggested high population must include the
neighbourhood determining that percentile and the proposed reduction, with a
deterministic full-resolution ordering. Any inexpensive modified-distribution
estimate is a search aid only. Recomputed gates, overlapping windows, transitions
and programme normalization require measurement of the actual candidate PCM.

## Bounded scope

- Stereo and conventional 5.1 programme audio, including existing interleaved,
  separate-file and M&E/dialogue composition sources. Processing acts on the
  complete mix and every channel, including LFE. Preserve mono/7.1 measurement
  and all existing workflows; their new automatic LRA method remains unqualified
  until separately covered.
- Inspect Source or Current. A changing candidate must use valid, matching
  Source/Current inputs and the current applied processing plan. Source remains
  immutable. Use Current's complete measured LRA population as the baseline.
- Generate analysis regions independently of annotation density. Group adjacent
  supported windows on the existing 100 ms analysis grid; use the union of their
  actual three-second support, bounded to programme frames. Overlap merges;
  a genuine gap remains separate. Do not label these regions editorial scenes.
  Preserve EOF-padding flags and never extend programme duration.
- Imported/manual markers provide names and navigation context. Their absence
  must not reduce the programme to one unusably large processing region. Merely
  adding redundant markers must not change ranking or a default candidate's
  audible envelope. An explicit “use these annotation boundaries” alternative is
  deferred, to avoid two strategies in the first slice.
- Suggested regions are temporary analysis objects. Selecting or viewing them
  does not alter annotations, seek, play, processing or device state. Existing
  ruler seeking and audition/inspection-domain separation remain authoritative.
- One attenuation-only proposal method, with the existing versioned 100 ms
  section-gain transition convention. Existing applied gains are the baseline;
  retain the approved additional/cumulative gain caps and global-gain bounds.
  No limiter, compressor, automatic makeup gain outside the explicit candidate,
  component processing or quiet-passage boost.
- Use the active versioned M16 preference values for changed-region count,
  affected-duration union including transitions, gain limits and planning
  reserve. Do not reset user limits or silently loosen defaults to achieve LRA.
  Count actual nonzero processing regions, not marker rows. Peak planning reserve
  stays inside the configured gain budget; it is not an LRA acceptance tolerance.
- First changing-proposal eligibility requires an available, stable programme
  LRA and selected IL, TP and LRA requirements. Additional selected M/S ceilings
  are also checked for every candidate; this method does not optimise them
  separately. Missing requirements/metrics or unsupported configurations give a
  specific unavailable explanation. All selected requirements already passing
  returns the existing genuine no-op, preserving gain, revisions and results.

This explicitly proposes one departure from the earlier “automatic segmentation
deferred” boundary: deterministic **analysis-region grouping**, needed for the
marker-free workflow. It neither changes editorial annotations automatically nor
adds semantic scene detection. Coordinated neighbouring-section strategies,
multiple competing algorithms and a general region editor remain outside M19.

## Measured candidate sequence

1. Build suggestions from validated baseline data. Explain the selected support,
   provisional gains, affected-duration union and each limiting policy. Offer a
   preview overlay distinguishable from annotations and applied processing.
2. On explicit **Assess candidate**, freeze source/member identities, ordered
   roles, range, current plan, requirements, method/version and policy. Construct
   one deterministic envelope. Retain a hard budget of **three full-programme
   private assessments**, using one verified source copy.
   Recommended initial search rule: use measured P10 plus the selected maximum
   LRA as a provisional upper-population level. Group included windows exceeding
   that level, rank their union regions by full-resolution excess/population
   support with an earliest-frame tie-break, and select within the current
   duration/count limits. For each selected region, derive attenuation from its
   high-percentile overshoot, capped by both additional and cumulative policy
   limits. Expose this provisional rule and resulting full-precision gains in
   the preview; it is not a new measurement or acceptance formula.
3. Measure a first local envelope through the ordinary processing/measurement
   worker. If appropriate, allow one deterministic refinement within the same
   selected high regions. Reserve one assessment for the final local-plus-global
   plan when its gain differs. Never measure partial regions and infer a complete
   programme verdict. Cancel, source/plan/requirements changes or exhausted
   resources stop the attempt without Apply.
   The optional refinement may only deepen attenuation of those same regions
   using the remaining measured LRA deficit, within unchanged caps. It cannot
   expand into neighbouring regions or begin a second strategy automatically.
4. Recompute complete IL, LRA, TP, M/S maxima, gate membership and percentiles
   from each actual candidate. A local LRA improvement alone is insufficient.
   Apply availability requires the exact final plan's measurement to pass **all**
   selected requirements. No rounding allowance or measurement-engine test
   tolerance enters acceptance. If the pass budget cannot verify the final plan,
   it remains unavailable for Apply.
5. Explain measured baseline-to-candidate changes separately from estimates.
   Do not attribute the combined delta causally to each region without separate
   measurements. Persist the eventual applied explanation with its baseline,
   requirement/policy/method and plan provenance, using M18's historical model.
6. **Apply** remains explicit and stopped-only. It changes Current and invalidates
   its previous measurements. Independently remeasure Current. Cancellation or
   failure leaves the applied plan unverified; keep its gains with no rollback,
   retry or hidden second correction. Show any actual failure/disagreement and
   stop. Verified export still independently measures the finalized rendered file
   and requires both Current and artifact to pass.

The suggested three-assessment cap limits search quality intentionally. If no
candidate passes, show the best actually measured attempt, unmet rules, policy
constraints and the scope of the search. Never promise arbitrary LRA attainment.

## Interfaces and dependencies

- `Core/Measurement/OfflineMeasurement.h` and `GatePopulation.h`: immutable full
  windows, gates/percentiles, metrics/status and exact sample/eighth-frame peaks.
- `StandaloneDiagnostics.*`, `MeasurementTimeline.*`, `TimelineGainLane.*`: share
  source-frame coordinates and full-resolution calculations. Display reduction
  never decides region scores or acceptance.
- `Core/Proposal/Methods/SectionTrimMethod.*`, `SectionPeakMap.*`,
  `StandaloneSectionProposal.*`: preserve current peak-led method; add a separate
  named/versioned LRA method, not a silent reinterpretation of M13 results.
- `StandaloneProcessing.*`, `Core/Processing/SectionGainPlan.*`: exact existing
  whole-programme envelope and global trim; no new DSP.
- `StandaloneSessionModel.h`, `StandaloneJobRequirements.*`,
  `AppliedProposalExplanation.*`, `SessionDocument.*`: immutable identities,
  stale/cancel handling, historical reasons and replaceable caches.
- UI/Timeline owner reviews the small suggestion-list/overlay interaction and
  shared interfaces. A comprehensive UI refactor is not a dependency.

Bind descriptive suggestions to measured snapshot plus region-method version.
Bind candidate authority additionally to processing, selected requirements and
policy. Marker-label-only changes can refresh presentation without remeasurement;
an explicitly changed processing geometry must invalidate the corresponding
candidate. Video offset/routing/monitor volume never change audio analysis identity.

## Acceptance and qualification proposal

Use deterministic two-/multi-level programmes, smooth ramps, brief high transients,
quiet excluded regions, silence, boundary/EOF transitions, ties, short/unstable
programmes and LFE-only high peaks. Test both a feasible LRA reduction and an honest
unachievable-within-budget outcome. Include a case where a visually plausible
histogram shift fails when gates and the complete candidate are remeasured.

Require exact source-frame intervals and padding labels; identical suggestions
and processing under no markers versus redundant/dense markers; unchanged
annotations until an explicit separate action; and accurate overlap/transition
union accounting. Test limits at, immediately inside and outside their bounds,
including cumulative existing corrections and zero/unset LRA differences.

Require no-op preservation when all selected rules pass, independent requirements
reevaluation, stable deterministic proposals, cancellation/source replacement,
policy/requirements/plan staleness, historical explanation association after
reopen/edit, measured private rejection, explicit Apply and independent Current
verification. Capture processed samples for offline/audition parity and retain
strict export gates with one final dither stage. No candidate cache may become
Current or exported-artifact authority.

Stabilize with short Sequoia tests, freeze one candidate, then run the affected
checks on actual High Sierra. M18 has completed the former combined qualification
deferral; do not carry that obsolete pending status forward. Reuse unchanged
M18/Intel evidence with exact identity/rationale. No recursive M1–M18 runners or
automatic repeat of the cinema/full three-hour workflow. Select video/physical
checks only for affected timing/buffering/UI contracts or a concrete problem;
coordinate physical audio use separately and retain the accepted picture-skip
limitation.

Proposed budgets: suggestions from a three-hour-sized full-resolution population
in <=2 s Sequoia / <=5 s High Sierra, <=64 MiB incremental analysis storage,
page at most 200 visible suggestions, and <=100 ms UI heartbeat p95 with an
appropriate 20 ms test timer. Cancel <=1 s; at most three full private passes,
one verified copy, no persistent candidate audio. Use synthetic population stress
plus one bounded 10-minute real stereo/5.1 workflow per Mac, rather than another
three-hour decode. Estimate 1–2 hours of selected qualification per Mac and
5 GiB scratch excluding installed tools and preserved fixtures. Measure actual
times/memory and record any unsupported product-duration claim honestly.

## Approved deterministic method contract

Version `lra-passages-v1` describes populations; `lra-selective-trim-v1` selects
and assesses processing. Use validated full-resolution M1 LRA windows, never
display bins. Stable sorting breaks equal loudness values by original window
index. Window inclusion, P10/P95 and rounding remain the qualified definitions.

The descriptive lower population includes gated windows at or below P10. The
upper population uses strictly greater than P95 without a job LRA ceiling, or
strictly greater than `P10 + selected LRA ceiling` with one. This latter level
is a provisional search threshold, not a compliance formula. Group overlapping
or touching three-second source-frame supports; clip displayed support at EOF
and retain padding counts. For each region, sum its selected windows' positive
LUFS excess above the threshold, in original time order. This score has units
LU × windows, not an allocated amount of programme LRA or unique audio energy.
Order by decreasing score, then decreasing window count, then increasing begin
frame, then increasing end frame. Markers only supply optional names/context.

For a high region, initial requested attenuation is its selected-window P95
minus the search threshold, plus the active planning reserve. The reserve is
planning-only headroom inside both gain caps; it is not an LRA allowance or a
guarantee about the resulting TP. Clip reduction to the additional cap and the
least cumulative headroom across the existing applied plateaus intersecting the
region. Report requested/achieved reserve and uncovered overshoot separately.
Use full precision and move a representable result inward when subtraction
would exceed the exact additional cap.

Visit ranked regions greedily. Tentatively add each complete region; retain it
only if the actual canonical before/after envelope support meets count/duration
preferences and the existing processing capacity. Do not shorten regions to fit
a budget or change adjacent unselected regions. Partition the existing plateaus
at selected region boundaries. Validate actual sample-position gain differences
at all linear-ramp segment endpoints: no additional local boost and no excess
additional/cumulative attenuation are allowed. If interacting existing ramps
make that impossible, explain the method limitation and skip the region.
This guard preserves the existing processing convention without claiming that
an arbitrary new partition preserves every prior transition automatically.

Assess the first selected envelope with the baseline global gain. If its measured
LRA still exceeds the ceiling, one optional refinement may deepen only the same
selected regions by that measured deficit, constrained by the original caps and
support budget. A refinement is never assumed to improve LRA. Compare actually
measured attempts lexicographically: an all-selected-rules pass first, then lower
measured LRA excess, fewer other failing rules, smaller sum of their positive
bound violations, smaller changed-envelope duration, smaller maximum additional
attenuation, and finally the earlier assessment index. Missing required metrics
cannot win or be applied. Retain the better measured attempt if refinement is
worse; show both actual outcomes and explain the choice.

Reserve one of the three full-programme assessments for the retained local plan
plus explicit final global trim. Calculate a scalar trim nearest the IL target
within IL/TP/M/S limits and -60 to +24 dB total global gain. Its estimates confer
no authority. Global gain may be positive after local attenuation; the complete
programme, including quiet passages, is then raised by that global amount.
Selective quiet-passage boosts remain deferred. If the retained exact plan
already passes every rule, preserve it; otherwise measure the exact combined
plan before Apply can become available. No more than three assessments, no
automatic expansion of regions, rollback, second strategy or relaxed preference.

An unsuccessful outcome distinguishes **measured rule failure**, **policy
restriction** (such as duration/count/gain caps), **method limitation** (such as
attenuation-only geometry), and **search budget**. It means this attempt did not
find an acceptable candidate under these conditions, not that no solution exists.
Only the exact measured passing complete plan can be applied. Historical reasons
and baseline/requirements/policy/method survive via M18 history; candidate snapshots
never become reopened Current verification.

## Coordinated follow-ups, outside M19

Track completion-pause profiling and scaling separately from this algorithm.
A later UI task should prominently expose initial import/source preparation and
analysis progress. A separate first-playback progress change is deferred while
session-wide reuse of one verified source copy is designed as coordinated
performance work. M19 retains existing preparation, cleanup and integrity rules.

Discuss an optional persistent “Copy audio into project” design before implementing
it. Show an explicit disk estimate and free-space check: copying encoded PCM24 at
48 kHz costs approximately 0.97 GiB/hour stereo or 2.90 GiB/hour 5.1, plus headers;
three-hour 5.1 is about 8.69 GiB. A float32 decoded copy would instead require
about 3.86 GiB/hour 5.1. Multiple members/components add their actual sizes.
Reopening must validate content/assembly identity and cache compatibility; a
project copy cannot be trusted merely because its filename or timestamp matches.
Relinking, partial-copy cleanup, ownership and user-visible deletion need design.
These are planning notes, not approval for bundling or changed preparation.

Raising quiet passages, freeze/preservation semantics, coordinated neighbour
strategies, dialogue analysis, new DSP, immersive expansion and VST3 parity remain
separately reviewed work.

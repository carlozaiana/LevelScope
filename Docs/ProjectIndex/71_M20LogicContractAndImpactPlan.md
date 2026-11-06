# M20 Logic contract and impact-selected development plan

6 November 2026. Recorded **before production edits** on
`codex/m20-region-planning`, handover `b4ccb66b70ea4dcad227d106925fc04db493881e`,
published main `22961f03de7919e1de631f8b13eb4b852108ac0b`.
[70](70_M20RegionFormationAndPredictiveGainPlanning.md) is the approved rule contract.
This stage ends with a committed, focused-tested core/interface for UI integration;
it is not completed M20 qualification or authorization to publish/merge.

## Immutable interface and identity

Add `Core/Proposal/RegionPlanning` value contracts, an auxiliary prediction collector
and model, and `StandaloneRegionProposal` as a worker service. Keep the existing
`SectionProposalRequest/Result` dispatch/Apply boundary. Optional, explicitly selected
`RegionPlanningSettings` activate M20; absent settings retain M19 and the LRA-unset
workflow with their original identities. Session schema 9 stores an independent
settings definition. No UI controls or automatic opt-in are added in this stage.

The UI edits a private draft and commits it through a validating model setter. Each
request copies settings, immutable measured Source/Current snapshots, processing,
requirements, preferences and selected section binding. Workers return shared const
results. Completion and Apply compare a freshly reconstructed binding; existing runtime
session/job validity independently rejects late completions after New/Open.

| Value contract | Contents and authority |
| --- | --- |
| `RegionPlanningSettings` | Automatic/selected-section/temporary-selection mode; section ID or frame span; margins; support-edge ranges/pins; internal nominal-cut ranges/pins; mandatory local-envelope protections; global permission; controls and comparison intent. Integer milliseconds and exact explicit frame coordinates. Automatic defaults allow global; selected defaults do not. |
| `RegionConstraints` | Resolved input/rate/frame binding, H/U/P, frozen original roots and effective permissions. Unknown mandatory protection, invalid ranges and unresolved selections fail closed. Validation does not mutate audio or annotations. |
| `RegionGeometry` / assignments | Original root, nominal cuts/parts, support edges, discovery provenance, reduction and reserve. Descendants share their root split cap. Zero interior parts count. Identical canonical audio removes redundant divisions. |
| `RegionActivity` | Full changed support including old/new ramps, programme-frame denominator, connected components, splits per root, interval capacity, added/cumulative attenuation, absolute dB-frame integral, longest run, margin/transition activity and global delta. Descriptive costs are not quality/compliance scores. |
| `PredictionSidecar` | Immutable request-owned 100 ms source-aligned weighted K-filtered energy sums/counts, per-channel sample maxima and separate EOF-drain bins; baseline snapshot/plan/ordered-role/method identity. One streaming auxiliary pass over the verified source with exact frozen processing. |
| `RegionPrediction` | Approximate IL/LRA/P10/P95, sampled M/S, population entry/exit, per-metric availability and uncertainty. Anchor model deltas to baseline windows, integrate actual squared gain ramps over bin intersections, and recompute complete gates. TP may be unknown/one-sided; loose bounds and sampled M/S cannot certify compliance or veto all candidates. |
| `RegionProposalReview` | Effective settings/binding, frozen root order, bounded rejection/omission reasons, proposed geometry/activity, all measured attempt summaries, exact retained plan, forecast errors and gate changes, resource/pass counts and ranking shortlist. Only complete measured all-rule passes enable Apply. |

An isolated auxiliary filter uses the qualified engine's coefficient expressions and
order, independently checked against its window energies. Do not alter authoritative
measurement arithmetic or its observer layout merely to obtain bins. Model limitations
include uniform within-bin energy, filter/gain noncommutation, boundary memory and
float-stage rounding. Near-EOF changes flag transformed-tail uncertainty; padded
support never extends programme duration or permissions. Identity proxies reproduce
the measured baseline. No candidate PCM or population per every proxy is retained.

M20 binding includes exact measured Source/Current, ordered roles/range, baseline plan,
requirements/policy, selected root ID/bounds, H/U/P, pins, controls/global permission,
method/model/rank IDs and runtime session revision. Exclude unselected annotations
and names; selected bound changes stale the request. Keep legacy binding bytes when
M20 is absent. Settings/requirements/protection changes reuse compatible measurements
and invalidate dependent proposals/approvals. Monitor/video metadata stays independent.

Apply rechecks actual envelope permissions and exact measured acceptance before using
existing stopped-only mutation and independent Current remeasurement. Cancellation or
failure retains applied gain and leaves Current unverified. Version the historical
explanation attachment to retain M20 settings/provenance/scalar results; never restore
private candidate windows, predictions or verdicts as current verification. Schema 8
sessions migrate in memory with M20 absent and identical audio; explicit saving writes
schema 9 requiring an M20 reader. Preserve Save As/original-file safeguards.

## Frozen evaluation hypotheses

IDs: `region-constraints-v1`, `region-formation-v1`, `region-gain-search-v1`,
`window-energy-proxy-v1`, `candidate-review-rank-v2`. Initial defaults, full-precision
comparisons, deterministic finite enumeration and tie-breaks follow document 70.
Reserve two-cut unchanged-middle alternatives, round-robin roots, legal unsplit/M19
challengers, shallower/deeper gains and exact canonical plan deduplication. Bound
roots/cuts/constructions/geometries/beam/proxies to 8/12/256/16/4/64, with proxy stages
32/16/16 and at most three complete candidate measurements. The 0.01 LU anchored
failed-result shortlist never enters acceptance thresholds. Record justified
enumeration/default changes and affected checks instead of silently substituting
whole-region-only planning.

## Explicit focused impact/test selection

| Area and shared callers | Fresh Logic-stage checks |
| --- | --- |
| Constraints, formation, search | Every control/range; splits versus components/capacity; exact off-grid pins; margins/clipping; protections and one-frame-outside U; independent sampled old/new-envelope oracle; redundant markers/equal gains; deterministic bounded enumeration, cancellation and omitted-alternative explanations |
| Auxiliary collector/model | Analytic constant gain/ramp overlaps; identical-plan identity; ordered all-channel weights and LFE maxima; bin energy versus qualified windows; gate entry/exit and EOF tail; unknown predictions; measured forecast errors and ranking reversals |
| Worker/binding/Apply | Matching Source/Current identities; no-op; strict all-selected-rule gates; stale/cancel/new-session completion; one verified copy; no fourth pass; retained measured failures and exact anchored shortlist |
| Audio/audition/export | Short stereo/5.1 at 44.1/48/96/192 kHz; representative reordered/assembly/composition callers; all-channel pre-monitor capture including LFE; seek/chunk/EOF/outside-U sample parity with fixed global; independent Current and finalized PCM24/final-only dither |
| Session/history | Schema 8 migration; schema 9 exact round trip; unknown constraints fail closed; recoverable caches; stopped/device-closed/muted reopening; historical explanations not authority; coherent saves and runtime identity |
| Existing methods | Explicit M19 core/stereo/gating/failures and representative surround/input groups; focused section-proposal, requirements, session/history and injected audition checks selected by touched callers. No recursive historical runners. |
| Build | Named standalone/new/focused test targets only; established JUCE and universal arm64/x86_64 build with Intel 10.13 minimum. Cross-build is not High Sierra execution. |

Use a new owned `build/m20-logic-2026-11-06` directory for commands, outputs, timing
and candidate records. Run short development groups, fix failures, then record the
stable core checkpoint with relevant textual evidence. Do not rebuild over packages.

Estimated resources: 5–10 GiB new universal build tree using existing dependency
sources; under 2 GiB short generated audio/session/export fixtures; roughly 15–30
minutes first build and 10–25 minutes selected tests, subject to measured host load.
Sidecar/model/search additional storage <=64 MiB, no full population per proxy.
Three-hour-sized synthetic populations test cached planning, not three-hour decoding:
<=10/20 s Sequoia/High Sierra targets, 30/60 s watchdog, <=1 s cancellation and native
p95 <=100 ms after UI integration. Watchdog/cancel provides no partial Apply authority.
Report source preparation, auxiliary pass, assessments, independent Current and export
separately. Clean up only task-owned partial data on storage failures.

## Reused evidence and combined qualification still required

M19 candidate-02 `cf275afb133118dd2a06e7bca59f82f98898b3cd` supplies baseline evidence,
including stereo success, gate-crossing ranking counterexample and honest cinema
failure. Use the existing eight-minute excerpt for bounded new offline evaluation;
do not rerun the full film. Reuse M18/Intel arithmetic/processor/video evidence only
where unchanged. UI convenience product `96159ebd79f762bea147b41ad696fd80f449dfc9`
supplies unchanged device/preparation/monitoring evidence; repeat affected shared
callers if touched. Earlier results are never fresh M20 execution evidence.

After Logic hands over the stable interface, UI owns advanced controls/overlays.
The frozen combined candidate needs actual Sequoia/High Sierra numerical/native
checks, targeted cross-platform/resource tests, the heterogeneous-root usefulness
gate and independent real-material engineering/listening review. Evaluate selected
control variants rather than a Cartesian product; freeze defaults before two real
holdouts. Coordinate audible/device tests separately. No automatic video/full-film/
three-hour workflow; reassess for concrete timing, latency, buffering, authoritative
arithmetic or dependency changes. Main remains published UI convenience/M19.
No M20 publication, checkpoint tag or merge is authorized.


## Development observations (6 November; not combined qualification)

The first integrated stereo workflow passed a complete candidate, explicit Apply,
independent Current measurement and independently measured PCM24 artifact. The
initial broad-root experiment did not find a passing candidate: nearly identical
boundary strengths crowded the twelve cut slots around one of two high passages,
and a low *average* need could hide loud material inside a proposed unchanged gap.
Historical failed development receipts remain under the owned build directory.

The finite pre-publication ordering in document 70 now reserves one paired interior
low-need-run hypothesis before the cut cap, uses maximum available need for the
spared-gap reservation, retains an inward window-support alternative, and explicitly
seeds its middle decision at zero. This changes discovery coverage, not permissions,
processing, defaults, measurement arithmetic or strict rule acceptance. Required
reruns: formation/core, broad-root split-versus-unsplit comparison, simple stereo/5.1,
full-resolution population proxy, cancellation/identity, and affected legacy callers.
A positive usefulness finding and focused-tested interface are required before UI
handover; these observations are not a claim of M20 qualification.

Completed Logic-stage receipts and UI contracts are in [72](72_M20LogicHandover.md).
The finite discovery correction passed the broad-root split-versus-unsplit usefulness
case; no production measurement arithmetic, gain processing or transport was changed.

# Logic implementation continuation and M20 boundary investigation

20 November 2026. The user requested a new chat named **Logic_Implementation_2**,
using **Astra 6.2 / Extra High**, to take over ongoing Logic responsibility.
Its immediate assignment is to finish the M20 boundary investigation and propose
general improvements for review. Production rule changes are not yet authorized
for this new slice. The previous M20 implementation and UI integration are complete
within their recorded scopes; their remaining musical-quality problem is the task.

Read the short [current task record](../LogicImplementationCurrentTask.md) first.
After an interruption or context compaction, reconcile that record with the latest
human instruction before resuming. Historical instructions to begin implementation
or send the UI handover describe earlier, completed assignments.

## Ownership and working state

Use the existing local checkout `/Users/fuzz/Developer/LevelScope`, branch
`codex/m20-region-planning`. No new worktree is needed for this investigation.
Logic_Implementation_2 becomes the active Logic owner; the previous Logic chat is
an idle historical reference. Strategy owns rule review and preservation trade-offs;
UI owns interaction/timeline work. Preserve the separate UI and video checkouts.

| Item | Verified state at handover |
| --- | --- |
| Logic HEAD | `8e50899d5487cef0466a8bd7ba30559047f972be` |
| Logic implementation | `4cbcca7f6b63b792b9a8ccbe73dbd5702ee311c5` |
| UI checkout | `/Users/fuzz/.codex/worktrees/a411/LevelScope` |
| UI branch | `codex/m20-complete-event-ui-2026-11-20` |
| UI HEAD and handover | `832952db5cbe89090869107dcfb9846bfe7340df` |
| Tested UI implementation | `8b28348dee77179b395253728d71818404bd5022` |
| Published main | `22961f03de7919e1de631f8b13eb4b852108ac0b` |

The UI continuation contains the Logic checkpoint. The Logic checkout does not yet
contain the new UI continuation. Read UI evidence at its absolute location; do not
mistake its absence in this checkout for unfinished UI implementation. Later
implementation must deliberately preserve/include that UI continuation in its base.
No checkout integration or source change is needed to complete the present report.

Tracked files were clean before Strategy added this handover, the current-task
record and its index entry. Those documentation changes are intentionally local.
Existing untracked `.DS_Store` files and `Tests/BusyPass2.levelscope` are unrelated;
preserve them. Preserve main, all existing tags, original media/sessions, frozen
packages, installed apps/plugins and system audio/Pro40 configuration.

Chat references, for bounded history lookup only:

- Strategy: `01a1ac62-28bd-7910-a131-59c3b1cea6d9`.
- Previous Logic: `01a13fdd-ecc3-7090-a004-bc921e58124f`.
- UI: `01a18e68-1087-7283-b605-7183e089a038`.

## Read the current contracts rather than all historical milestones

1. [70 M20 design](70_M20RegionFormationAndPredictiveGainPlanning.md), especially
   the boundary/duration amendment and accepted refinements.
2. [80 accepted final review](80_M20CompleteEventFinalReview.md), the authoritative
   implemented event-coverage, ranking, permission and budget contract.
3. [81 Logic implementation and evidence](81_M20EventImplementationAndImpact.md)
   and [82 interface handover](82_M20CompleteEventLogicHandover.md).
4. [83 completed UI integration](/Users/fuzz/.codex/worktrees/a411/LevelScope/Docs/ProjectIndex/83_M20CompleteEventUIIntegration.md).
   Its newer completion status supersedes statements in 81/82 that UI is pending.
5. [69 shared development approach](69_ProposalStrategyAndAudioQualityHandover.md)
   and relevant sections of [ProposalEngineSpec_v2](../ProposalEngineSpec_v2.md)
   for architecture, Source/Current/artifact truth, predictive planning and future DSP.
   Their old milestone status and design-only dispatches are historical.

Use [79 initial review](79_M20BoundaryAndDurationLogicReview.md) where it explains
a decision; 70/80 contain the accepted refinements. Do not reread the entire chat
or all historical test runners to establish the current assignment.

## The actual problem and the user's judgement

Primary development file:
`/Users/fuzz/N15/Mixdown/MusicInBackAndForegroundSomeDialog.wav`.

- SHA-256: `7ba6404041c94c909df290a6829a0fe92ce594329a1d575490c2bb4a5b6cc8f0`.
- Recorded source: stereo, 48 kHz, 13,228,800 frames, 275.6 seconds.
- Original IL approximately -20.523327 LUFS, TP -0.035708 dBTP, LRA 19.042929 LU.

The user listened to the preserved M19/M20 projects at
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/ui-inspection-2026-11-07/listening-gate-04/`,
under `m19/inspection.levelscope` and `m20/inspection.levelscope`. Both pass their
selected requirements, but the user sees no useful improvement in the boundaries.

First passage: the user sees a region beginning around 12.827 s, but hears the
musical entrance/first bass-drum hit shortly after approximately 12.283 s. Correction
should enter before that first relevant hit; the old ending also misses a relevant
final loud hit. Distinguish approximate UI readings, exact stored frames, detector
decision time and proposed transition geometry.

Later passage: the M19 region around 49.906–60.196 s, attenuated by 3 dB, sounds too
long at both ends. Investigate the louder activity around 52.615–57.762 s and its
surrounding material. These times are listening landmarks, not constants or exact
endpoints to force. Relevant decay and transition placement still require evidence.

Develop and demonstrate a general rule on this file first. Additional material is
a later generalisation check; it must not postpone solving this development case.
Numerical compliance alone does not establish useful event formation or mix quality.

## Salvage and verify the investigation already performed

Existing diagnostics are in
`/Users/fuzz/Developer/LevelScope/build/m20-boundary-investigation-2026-11-20/`:

- `identity-before.json`, `diagnostic.json`, `compile-commands.json`, build/run logs.
- `trace.csv`: detector decisions, contrasts, backgrounds and active-group state.
- `groups.csv`, `windows.csv`, `geometry.csv`: grouping, measured windows and nominations.
- `saved-plan-summary.json`: exact local intervals/global gains from the older sessions.
- `prepare_trace.py`, `trace.cpp`, `DiagnosticEventFeatures.*`: diagnostic-only
  instrumented collector and harness, with the production collector used for comparison.

`diagnostic.json` reports matching original/instrumented event identities, seven
groups, zero supported groups, no new correction assessments and no device opening.
It records schema-2 automatic planning, explicit 1 s margins and the feature settings.
Treat the equality assertion as a useful check, not proof that every diagnostic
interpretation is correct. Preserve these files; review/reuse them before rerunning.
Record useful evidence and provenance in durable repository documentation/results
when producing the report, rather than leaving the only explanation in this chat.

Established code mechanisms to trace:

- `Source/Core/Proposal/EventFeatures.cpp` creates a group only when no group is
  active. Later attacks update its last attack and decay reference without creating
  a competing entrance. Original-onset corroboration has an 80 ms upgrade window.
- Closure depends on return to the background captured at group start and the
  persistence/last-attack conditions. The background-return test permits values
  below that background; it is not an absolute plus/minus 3 dB band.
- `Source/Core/Proposal/EventCoverage.cpp::formEventRegions` nominates preferred
  event geometry only from need-associated, corroborated groups with supported
  endings and no incomplete flag. Existing loudness-based alternatives remain.
- The saved music run has a tentative first relevant onset around 10.05–10.16 s;
  all seven groups are tentative/incomplete. This explains why event reporting
  could change without introducing the desired event-based geometry.

The previous Logic chat's preliminary trace interpretation says the first strong
hit is corroborated by about 12.37 s but absorbed into the 10.05 s group. It also
reports a strong final hit around 30.26 s after the old correction returns to zero.
Verify those interpretations and source-frame mapping against the saved trace;
12.37 s is a detection/corroboration time, not a prescribed boundary. The later
passage still needs its own complete explanation. No final improvement proposal
was delivered: after compaction the old chat repeated the completed UI handover.

For comparison, original reports are in
`Tests/Qualification/M20Events/results/logic-01/reports/`. The 60% margin run's
selected numerical pass affects about 17.74% including transitions, so its duration
cap is not binding. The approximately 12.283 s entrance with a 100 ms lead-in fits
the first root's recorded 1 s permission. Simply enlarging allowances is not a fix.

## Required report and proposal

Complete the investigation, then deliver one concrete recommendation for review:

1. Explain each passage through detection, corroboration, grouping, association
   with overlapping measurement windows, candidate formation, predictive ranking,
   exact full-plan assessment and selection. Identify where useful alternatives
   disappear; distinguish observed facts, inference and missing evidence.
2. Compare a distinct entrance within continuing activity; a small bounded set of
   competing start/end/group hypotheses; and tentative boundary nomination separate
   from the stronger evidence needed to claim supported-complete treatment. Analyse
   endings/releases as carefully as starts. Include momentary/short-term loudness,
   energy and peaks where useful; a peak alone is not proof of correction need.
3. Recommend general rules with required data, explicit parameter meanings,
   defaults with reasoning, failure cases, resource/search costs, version/session
   effects, UI responsibilities and a small implementation sequence. Identify any
   change to the accepted contract rather than silently weakening its guarantees.
4. Specify acceptance on this file: first relevant hit covered, justified treatment
   through the relevant ending, tighter later geometry where supported, appropriate
   transitions and an intelligible explanation of each proposed processing region.
   Compare full plans under unchanged requirements and plan coordinated listening.
5. Preserve shallower/deeper gain and region-combination comparisons. After useful
   results here, validate the frozen rules on additional material. Include success,
   gate-crossing and cinema checks without requiring the cinema excerpt to reach
   15 LU; any remaining failure must remain honestly reported.

## Contract essentials to preserve or explicitly propose changing

- A section is an organisational annotation, an analysis region is measured, and
  a processing region is proposed/applied correction. Preserve annotations and show
  processing geometry separately; no broad code/UI rename.
- Separate correction need from entry/exit placement. Account for surrounding
  material, overlapping windows, whole-programme gates, existing transitions and
  protected material. Reading context grants no processing permission.
- Strict allowed areas must not silently expand. Approximate hints, explicitly
  permitted margins, pinned/user/imported boundaries and supported protection
  semantics remain distinct; preserve frozen ownership between roots.
- Two splits means two new internal boundaries per original root, not two resulting
  regions. Keep this separate from changed-root limits, exact affected-duration
  union including transitions, gain caps and search limits. Dividing an identical
  envelope is not a new audible operation; splitting cannot evade restrictions.
- Current defaults include 20% affected duration, 3 dB additional/12 dB cumulative
  attenuation and reserve inside the cap. 60% and outward margins are explicit
  experiments, not silently changed defaults. Read exact applicability in 70/80.
- Preserve exact compliance thresholds. The failed-candidate LRA near-tie shortlist
  is a ranking device, not relaxed compliance. Retain the distinction between
  predictions and actual complete-plan measurements, independent Current and
  finalized-artifact verification.
- Supported-complete preference, partial/uncertain labels and explicit measured
  alternative selection are separate. Inspection cannot grant Apply authority;
  Apply is explicit and stopped-only, and leaves Current unverified. History cannot
  restore live candidate authority. Preserve identity, cancellation and stale-work guards.
- Existing bounds include 64 proxy plans, three full assessments, the 10 MiB feature
  reservation inside 64 MiB planning, and the formation caps in 80. Quantify any
  proposal's cost before increasing a budget; do not silently widen search.

## Qualification and authorization

Logic evidence records 26 groups/3,919 assertions on native Sequoia and an Intel
10.13 universal cross-build. UI evidence records 2,767 focused assertions on
Sequoia, including the real selection/Apply workflow. This does not establish
musical quality or actual High Sierra execution of the new combined method.
Broader detector controls, independent material, actual High Sierra/new-method
checks and coordinated listening remain pending. The gate and cinema examples
remain measured failures; music has measured passes with partial/uncertain coverage.

The first assignment authorizes code/evidence reading, bounded unchanged-rule
offline diagnostics if needed, and documentation/report writing. Proceed without
another permission request for those actions. No production detector/UI/DSP edits,
altered-rule correction experiments, publication, tags, main merge or physical audio
devices are part of this assignment. Present the improvement proposal before those
implementation steps. Do not reopen the completed UI dispatch or send unsolicited
instructions to another chat. Preserve known video limitations and audio/drift guards.

## Longer term Logic responsibility

LevelScope's proposal engine extends beyond region gain. The approved approach is
to develop each processing capability with a first useful selection method, keeping
processors, selection algorithms, requirements and preservation preferences distinct.
Future limiter/compressor work needs explicit processing-order, state/latency,
identity, persistence, resource, prediction and measured-verification contracts.
Existing plugin modules are not proof of qualified standalone integration.

Consult the DSP inventory and parameter/measurement contracts when that milestone
begins. Do not make a general processing-graph rewrite or the complete dynamics
suite a prerequisite for fixing M20. Maintain durable decisions and evidence so
successive focused chats can continue the Logic responsibility reliably.

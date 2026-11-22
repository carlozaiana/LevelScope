# M20 combined focused qualification — 22 November 2026

Read with [88](88_M20CompetingBoundaryUIIntegration.md), [87](87_M20InnerBoundaryLogicHandover.md)
and the unchanged definition/limits in 85–86. The user authorized remaining silent
qualification, coordination with Logic_Implementation_2 and a separately identified
replacement only if a concrete defect required it. No publication or main integration.
Both manuals, installations, original media/sessions and listening materials are protected.

The finite [plan and amendments](../../Tests/Qualification/M20Events/results/combined-2026-11-22/plan.md)
were recorded before the corresponding outcomes. The two independent inputs and
criteria were frozen by Logic before either measurement. Golem was reached through
the user's existing SSH master and verified as actual High Sierra 10.13.6 / 17G66,
x86_64. Sequoia is 15.7.9 arm64. No physical output was opened by these checks.

## Candidate identities and concrete UI correction

Original candidate 01 remains byte intact: implementation
`ce31c4679782c69d02c29ddc675fee086b89395c`, handover
`2ab8b6616c582e3544a1203ad11d2615c08337cf`, package
`8f1497d12d4d5a1c3a0e8cd4aa3e864cc70a540023a399041d5dbd808ac93e79`,
measurement build `aa2ef42883b323fbf6a2dc84dcc40d6fafdcc1ec5db8cc1ac1e74fdead997531`.

Its 281 competing-boundary native checks passed on Golem. However, the existing
event-review heartbeat exceeded the unchanged 100 ms p95 limit twice: 108.190 and
107.571 ms. This was not waived. A locator-copy formatting error before that run
was separately corrected in test tooling; the product correctly rejected the
noncanonical JSON, and no original session was changed.

Separate instrumentation localized recurring long gaps to switching, with proposal
heartbeat p95 29.180 ms and switch-loop p95 109.726 ms. The approximately 2 s maximum
spanned scripted synchronous guard manipulations; PNG captures cost 136/148 ms.
Neither explains the recurring p95 miss. Detailed callbacks and a native CPU sample
are archived with the failed receipts.

Candidate 02, commit `85ec92b224c8c6010e5dd60f9529efa6d7e110ee`, eliminated duplicate
hidden-report polling. It passed 593 affected Sequoia assertions, but Golem still
measured 107.866 ms. It is preserved as an unsuccessful performance fix, not a
qualified replacement. Native sampling then found repeated coverage-summary
TextEditor/CoreText layout: 398 of 875 samples within the switch callback. The long
report was already cached; its summary was recreated on each switch.

Candidate 03 adds a retained scrollable summary for each bounded report slot. It
updates text when status/details change, preserves all explanation access and keeps
the existing report-page state. Embedded reports have one parent refresh owner;
hidden inspection windows do not poll. Explicit refresh and all Select/Apply-time
checks remain intact. Only four UI implementation/header files differ from candidate
01. No Core, dispatcher, model authority, measurement, processing, session codec,
export, proposal method/default, gain strategy or boundary mathematics changed.

Candidate 03 implementation: `2f90dab6eaa25f9d2ce9f033341c948b2994d9da`.
Package identity: `f44c1fea69a41416c3c1dd0284a3a5b9cd52757bc484f9623cbb7fb13aa48dd3`.
Measurement build: `5cf1b904db533a8ee8379325833a6e0bf8efd4f0e8f88956a654ce9fa4b026e7`.
Source manifest: `01f8d5fbd7f4392ea83d640efab1f6c93da4ced0566c4648dfc82b8d66a5ad18`.

Exact preserved application:
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m20-competing-boundary-ui-2026-11-21/candidate-03/LevelScope Standalone.app`.
Universal arm64/x86_64; the Intel slice targets 10.13. Actual execution is recorded
separately below; a deployment target is not runtime evidence.

## Final replacement results and decision

All three affected groups pass from frozen candidate 03 on both actual platforms:

| Native group | Assertions per platform | Sequoia seconds | High Sierra seconds |
| --- | ---: | ---: | ---: |
| Event review and alternative authority | 146 | 22.128 | 44.588 |
| Competing-boundary integration and polling lifecycle | 286 | 46.674 | 98.493 |
| Generic review, manual UI and dense history | 163 | 18.325 | 42.229 |
| **Total per platform** | **595** | **87.127** | **185.310** |

The original <=100 ms p95 gate is unchanged. Legacy review switching p95 is
7.881 / 17.527 ms (Sequoia / High Sierra), heartbeat 33.797 / 78.693 ms.
Competing review switching is 2.812 / 6.053 ms, heartbeat 22.998 / 33.397 ms.
High Sierra legacy cancellation is 83.557 ms for Propose and 324.792 ms for
Current remeasurement. New polling checks cover hidden/visible standalone reports,
parent-owned embedded reports and explicit refresh; summary reuse checks retain
the exact complete/selected text. All stale, pending-edit, New/Open, selection,
matching-plan Apply and separate verification guards pass. Settings schemas 1/2/3,
independent margins, prospective Save As/history-only schema 3, supported standard
no-LRA/global-only/mono paths and dense 10,000-interval / 2,280-page access pass.

The final development session's Source/Current numeric payloads and exact processing
match both the original candidate and the other actual platform, bit for bit. The
build identity is intentionally different. No reference values or tolerances were
updated to pass. The four-file Source difference map justifies inherited engine,
session exchange, service/export, resource and independent-programme evidence.

The heartbeat is a percentile check spanning the scripted workflow, including
synchronous guard manipulations, initial report construction and snapshots. Its
full-run maxima remain recorded: High Sierra legacy 3.647 s and competing 2.976 s;
Sequoia 0.684 s and 0.592 s. These receipts are not a maximum-interaction-latency
guarantee. The inherited completion pause and sustained-media qualification remain
separate; no sustained freeze or audio interruption is waived.

**The agreed finite development qualification supports a later main merge**, with
the unsuccessful independent searches and remaining limitations below explicitly
retained. It does not claim general correction success, independent holdout listening
acceptance or video/endurance qualification. Main integration still needs separate
user authorization.

[Acceptance and preservation receipt](../../Tests/Qualification/M20Events/results/combined-2026-11-22/acceptance.json)
records the final identity, exact platform groups, Source reuse map and integrity
checks. Nine preserved packages (the seven originals plus candidates 02/03), all
26 tag objects/targets, 72 manual files, 15 listening files and other worktrees
remain intact. Main and tracking main remain
`22961f03de7919e1de631f8b13eb4b852108ac0b`. Unrelated `.DS_Store` files were preserved.
Candidate manifests remain immutable build-time records; this separate acceptance
receipt records the later qualification outcome.

Short manual review: open a copy of the supplied session in candidate 03 and allow
its build-specific preparation to finish; switch between Complete explanation,
actual attempts and shared/proxy evidence; verify status/coverage remain current
after selecting an eligible measured pass. Apply separately while stopped, then
Re-measure Current, Save As and reopen. Closing/reopening Review should retain a
responsive, current status. No physical playback is needed for this checklist.

## Original-candidate evidence retained for unchanged paths

The integrated Sequoia six-group / 3,004-assertion run and Logic 28-group /
3,556-assertion evidence remain inherited by their recorded source hashes. The
replacement's affected native UI runs supersede those UI results; unaffected
timeline/engine/control evidence is reused rather than recursively rerun.

Actual Golem additions from candidate 01:

| Group | Assertions | Seconds | Result |
| --- | ---: | ---: | --- |
| Competing native UI | 281 | 95.296 | Pass |
| Stereo service / Current / finalized export | 60 | 14.935 | Pass |
| 5.1 service / Current / finalized export | 57 | 21.360 | Pass |
| Fast synthetic three-hour-size population | 17 | 9.633 | Pass |

The population test processes synthetic arrays, not hours of audio. Its cached
64-plan calculation took 2.838 s against the 20 s Intel target; actual scoped heap
was 54,279,744 bytes, below 64 MiB. Maximum feature reservation was 9,946,592 bytes,
below 10 MiB. This is incremental observed allocation, not a whole-process RSS bound.
Original native review switching p95 was 45.698 ms; competing heartbeat p95 35.237 ms.
Cancellation and no-late-authority assertions passed without physical outputs.

Independent native development reassessments on both architectures matched in
exact processing, settings and full measurement payloads, including every measured
window, channel peak magnitude and sample position. The read-only checksum parser
records full numeric-payload SHA-256, not just rounded display values.

Schema-11/programme-9/history-4 exchange passed 35 receiving checks in each direction,
plus 34 same-platform positive-control checks. Each receiver discarded foreign
verification, independently measured Source and applied Current, and saved/reopened
its valid local caches. No live proposal authority was revived. Confirmed roles,
nominal interval boundaries/gains/global gain, settings, requirements, annotations,
policy and complete applied-history bytes matched. All recorded scalar, gate,
population, percentile and per-channel peak/position differences were zero.

The supplemental executable is separately identified by
`1161e56c07cf36c598b2e7b00d8c8797fb7b66da51d99694674a2903bcdbc67b`.
It used frozen Source and compiler/link inputs. Re-linking the original test object
first reproduced the frozen executable byte for byte. Its evidence is inherited for
unchanged engine/codec paths, not described as a run of the replacement app binary.
The Source-wide measurement-build identity changes for a UI-only edit, so old-build
caches still require independent measurement under the existing policy.

## Independent programmes and listening

Protocol SHA-256: `bb723fa6cd4eee849aee9673d655b034c4a08d7b256a6cbb2c171837c8dab853`.
Full user-supplied DialogAmbiance stereo and 5.1_Busy conventional 5.1 were selected
before outcomes, with confirmed ordered roles matching earlier records. Both are
PCM24 / 48 kHz. These are holdouts for this competing-boundary change, with prior
testing history; the music excerpt is the development comparison.

Fixed criteria: −23 LUFS ±0.2, −1 dBTP, no M/S ceilings; LRA ceiling in tenths is
floor(10 × (unrounded baseline LRA − 1)). Schema 3 automatic, explicit 1 s per side,
20% local-duration allowance, unchanged gain/reserve caps; no outcome-driven tuning.

| Programme | Baseline LRA / ceiling | Outcome | Assessment / scoped heap |
| --- | --- | --- | --- |
| DialogAmbiance | 16.660717 / 15.6 LU | Three measured failures; TP +2.627 to +2.999 dBTP | 5.778 s / 2,675,136 B |
| 5.1_Busy | 15.548069 / 14.5 LU | Three measured failures; TP +3.360 to +3.704 dBTP; two also just below −23.2 LUFS | 9.398 s / 3,487,280 B |

Both bounded searches correctly denied selection/Apply: 27 stereo and 24 surround
software assertions passed. No corrected Current or verified delivery export was
produced for these holdouts. This is honest unsuccessful correction, not evidence
that no solution exists. Stereo cancellation returned in 0.714 ms without authority.
An initial sandbox attempt could not create the ordinary verified-source cache;
the identical executable/input/settings passed with normal cache access, before any
proposal was produced. The infrastructure receipt remains preserved.

The user separately reproduced **Current – Passed** on the development programme,
then completed the requested uninterrupted listen and passage comparisons. Their
exact listening report: “I did the listen and passage comparisons and the results
are acceptable. The boundaries could and will be improved, but they are ok in this
stage of development. Overall sound is very good as expected and without any
distortion or clicks.” No independent-holdout listening or timecoded alternative
ranking is inferred. Numeric PASS and boundary/sound acceptance remain distinct.

## Scope, limitations and later direction

The inherited completion pause remains separate. This finite qualification is not
a maximum-latency, full-film, physical-output, video or endurance claim. Broad video,
full-film and three-hour endurance remain deferred; no substitute 90-minute run.
Occasional video frame drops retain their existing limited acceptance; audio
interruption, accumulating drift and sustained freezes are not waived. Development
packages remain unsigned and are not installed or notarized here.

After a separately authorized main integration: remembered Propose settings; one
clearly described action to apply an eligible measured passing plan and independently
remeasure Current; visible separation of numeric PASS from listening review; useful
passing alternatives with safe replacement/recovery; direct dragging of applied
region boundaries invalidating Current; optional advanced inspection; simpler
extension limits and later gain-strategy evaluations. None is implemented here.

Main, tracking main and every checkpoint remain unchanged. No push, new tag or main
merge is part of this assignment. Both manuals are untouched.

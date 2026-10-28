# M13 implementation impact and qualification plan

Publication update, 28 October 2026: combined candidate-02 completed actual-platform
qualification and the user subsequently authorized commits, an annotated tag and
main integration. The [M13 checkpoint](https://github.com/carlozaiana/LevelScope/tree/m13-bounded-section-proposals-2026-10-28) records the accepted source.
The original development-phase restrictions below are retained as historical
plan/authorization boundaries; the later approval permits this publication.

Recorded before source edits, 28 October 2026. Approved scope: [47](47_M13BoundedSectionTrimProposals.md).
Work stays in the existing `codex/logic-module` checkout at documentation HEAD
`0fc4c5a7f2364f57760228f5c5ecf2f6e0714197`, above frozen M12 candidate-03
`2fb0e3cfee22a527bacaa09c4ddace23129442fb`. Preserve its package and all tags.
Main remains qualified M11 `cfb1033aa7fb65a2655db276bee696161196bdeb`.
No checkpoint publication or merge is authorized by implementation approval.

## Explicit affected behavior and selected checks

| Change / shared dependencies | Required fresh evidence on the frozen combined candidate |
| --- | --- |
| Optional continuous true-peak observation in the measurement kernel | Collector on/off numeric equality; independent localized-peak/boundary/padding oracle; all channels including LFE, all four rates; M1 measurement tests and pinned EBU v05/ITU BS.2217 suites on both actual Macs |
| Pure bounded section planner and private worker | No-op/global-first, exact budgets/reserve/capping, gain and geometry validation, all 3 private-pass branches, stable ordering, complete candidate assessment, profile/application separation, explicit limiting-constraint errors; M2/M3 focused regressions |
| Shared verified-source preparation and processing callers | One verified input reused across private passes; injected read/copy/disk/cancel failures; no partial audio; unchanged processing samples and Source; interleaved/group/composition matrix; M6/M9/M11 safety modes and M12 processing regression |
| Candidate identity/model/UI | Private candidate cannot become Current/artifact truth; stale/forged/late result rejection; one atomic stopped-only Apply; independent Current remeasurement; cancel/fail preserves applied UNVERIFIED state; annotation/profile/routing/video identity rules; no-change zero worker calls; silent UI on both Macs |
| Persistence/report/export callers | Schema 1–4 compatibility, coherent save and cache recovery, no private result authority on reopen, exact applied plan; PCM24 single-dither oracle and strict Current/artifact gates; M8 short/M10 focused regressions |
| Shared UI/analysis | Existing M4/M5 inspection, rank and geometry regressions; private job progress/cancellation and small-window reachability; full-resolution planning, no event-count claims |
| Long workload | One consolidated 3-hour 48 kHz 5.1 composition per Mac, 1,000 annotations, bounded proposal/Apply/verify, captures/seek/EOF, resource/session and real RF64 export; independent bounded sample oracle and complete artifact hash. No duplicate historical long runs |
| Post-freeze transport/video integration | H.264 CFR/VFR and ProRes exact/fractional rates × Source/Current: eight short cases per Mac, section/global plan, signed rational offsets, seek across transitions, resume/restart/EOF; cross-Mac session exchange; selected coordinated muted physical subset |
| Real cinema input, Sequoia only | One bounded measured evaluation of the supplied 5.1 file with immutable original/content identity. User-reported 28 October Nuendo values are external context; Cubase 10 screenshot is separate. No forced pass; dialogue values do not become M13 objectives |

Use explicitly named standalone/test targets and selected executable modes; inspect
their call graph before running. Never invoke recursive historical milestone runners.
The completion validator enumerates required cases/receipts and rejects missing,
failed or mismatched evidence. Added tests must exercise successful automatic local
correction and honest failure, not only rejection. Short development loops precede
frozen qualification; rerun affected checks after any implementation fix.

## Reuse boundaries

M12 candidate-03 source digest
`d6501c0d234ccf9c02b90527fc4af6c613dd4d538a71cc5520fb6b7b21085c1f`
and standalone executable
`ab319348a5ca586fd9048581165f7fdefb0cb1bd7221d827c78fd3740ec2aec4`
are the preserved Sequoia-only baseline. Reuse original long media fixtures and
unmodified codec/colour evidence with recorded source/dependency equivalence.
Reuse a result only when the tested executable/path, relevant implementation,
environment and assertion remain equivalent; preserve original identities.
Shared kernel observation, model, job and build changes require affected regressions,
not a blanket claim that M12/M11 evidence qualifies M13. M12 has no actual High Sierra
section-processing evidence: run that feature's required checks in the combined phase.
Full old codec matrices and 30-minute endurance runs are conditional on changed
timing/buffering/presentation behavior or concrete failures, not automatic requirements.

## Phases, budgets and coordination

1. Sequoia audio development and silent UI; Intel slice built locally with the
   established JUCE/deployment configuration (10.13). No Golem/video execution yet.
2. Freeze a new universal combined candidate and source/artifact manifest.
3. Actual High Sierra runtime, cross-platform session/results comparison and selected
   deferred video qualification of that same candidate. Coordinate access and device
   use; an Intel slice on Sequoia never substitutes for High Sierra execution.

At preflight Sable has approximately 203 GiB free and Tucan 256 GiB. Use a new
`build/logic-m13` build and `build/m13-qualification-2026-10-28` evidence area, preserving
all prior packages. Reuse immutable verified long inputs; require 65 GiB free per
long run with exact input/copy/export size preflight. No media enters GitHub.
One proposal envelope, at most three full attribution/measurement passes; verified
copy/hash I/O is separately recorded. No unbounded search or automatic retry.

Keep approved limits: <=16 MiB retained proposal/attribution/plan data at 10,000
sections; <=256 MiB additional working RSS above existing measurement storage;
UI heartbeat p95 <=100 ms; cancel acknowledgement <=1 s; save/parse <=5 s;
session <=128 MiB; JSON+HTML <=16 MiB. Long workflow <=240 min Sequoia / <=300 min
High Sierra; estimates 90–150 / 120–210 min. Short/reference work estimate 30–60 min,
silent UI 20–40 min and short video 15–30 min per Mac, excluding builds/fixes and
operator waits. Record actuals and cleanup; do not silently widen missed budgets.

The selected physical subset is four eight-second cases per Mac: H.264 VFR and
ProRes 422, each with Source and the saved, applied Current section/global plan.
Use the approved built-in output only after coordinating this M13 block, with every
route muted. The eight-case injected video set separately includes H.264 CFR and
fractional-rate ProRes. No 30-minute endurance run is automatic. Silent UI evidence
also covers combined M12 manual drafting/Apply, full-precision entry, annotation
independence, stale drafts and cache recovery on actual High Sierra; earlier M12
Sequoia observations are not substituted for that runtime evidence.

Occasional picture skips are accepted with counts/timing retained. Audio interruption,
accumulating drift, sustained freeze, identity and memory/lifecycle failures still
block completion. Do not rerun only to obtain zero skips. Physical multichannel
listening is separate. Devices stay closed for audio-only/injected checks; any
physical or audible run is separately coordinated. Preserve installations, plugins,
Pro40, system audio and the separate video worktree.

Until the combined actual-platform/video phase passes, report exactly what is
complete and what remains pending. Stop for readiness review before publication
or main integration.

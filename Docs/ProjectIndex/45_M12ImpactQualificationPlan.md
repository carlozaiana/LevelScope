# M12 change-impact plan — Sequoia phase

Recorded 27 October 2026 before implementation. Approved contract: [44](44_M12SectionTrimProcessing.md).
Base/main remains `cfb1033aa7fb65a2655db276bee696161196bdeb`; work stays on
`codex/logic-module`. No checkpoint publication or main integration in this phase.

## Required fresh checks

| Change and shared callers | Required Sequoia evidence |
| --- | --- |
| Typed section plan and shared sample processor | Independent gain/ramp oracle, full precision, zero/no-op, invalid plans, frame/block/seek invariance, all rates and stereo/5.1 including LFE |
| Measurement/audition/export dispatch | Source remains identity; offline/pre-monitor/pre-quantization parity; audio-only seek/EOF/continuity, overload/recovery/cancel/owned-copy failure handling; pinned EBU and ITU vectors |
| Current identity, M3 and model | Draft/apply separation, exact input/plan binding, annotation independence, stale rejection, failed verification, global trim preserving sections, already-pass no-op |
| Schema 4, caches, reports and export | Exact round trips, schemas 1–3 migration, recoverable caches, coherent saves/late jobs, report plan identity and strict finalized-artifact gates, deterministic single dither |
| UI and build dependencies | Silent Sequoia manual-section workflow, draft/applied curves, stopped-only edits, unverified/actual results, save/discard; named standalone/test targets, universal architecture and 10.13 deployment inspection |
| Sustained work | One consolidated three-hour 48 kHz 5.1 composition with 1,000 intervals, measurements/trim/verification, captures, sections/diagnostics, session and real verified RF64 export; resource receipts |

Use short targeted tests while developing. The frozen runner must enumerate required
cases and fail missing Sequoia evidence; do not recursively invoke historical
qualification runners. Inspect individual tests before running them: tests that
load/decode/play video are excluded even when their target also contains audio tests.
Reuse helper code only after checking its call graph and mode selection.

## Reuse and deferral

M11 candidate-02 source manifest
`669b350652dbc3fb09d0ff724758208fec21bba40e0f6bc0188b89e6404c2c61`
and app executable
`deb189e943529dd3e9e9d63a93e1d1cd2003f36752f7ab43f1e099d6c204eae5`
remain original historical evidence for unchanged composition, measurement kernels,
container/quantizer and video components. Record final dependency/diff equivalence
before reuse; rerun affected callers. This evidence does not qualify new features.

All actual Golem/High Sierra work, cross-Mac exchange and video-loading/decoding/
playback tests are explicitly deferred until after M13. Do not access Golem. Build
the Intel slice locally with established dependencies and the 10.13 target; never
call that High Sierra runtime qualification. Preserve serialized rational video
alignment without loading video. Audio-only position/seek/duration/EOF tests remain
required. Explain any concrete timing/buffering/presentation/video change or failure
before proposing earlier deferred tests. Missing required Sequoia checks block this
phase; deferred checks are not passes.

## Budgets and execution

At preflight Sable has 213 GiB available; Tucan is mounted with 256 GiB available.
Use a new build directory and named standalone/test targets only. Reuse validated
long inputs, no equivalent historical long reruns or duplicate frozen packages.
Reserve ~65 GB working headroom. Envelope storage <=4 MiB for 10,000 intervals;
additional working RSS <=256 MiB above existing measured-result storage; UI p95
<=100 ms; cancel <=1 s; save/parse <=5 s; session <=128 MiB; reports <=16 MiB.
Long workflow <=180 min (estimate 60–90 min), short/reference 20–40 min, silent UI
20–40 min, excluding development/builds. Record actuals and owned temporary cleanup.

Physical devices remain closed for silent UI; audio capture uses injected devices.
No audible test without coordination. Keep installations, tags, video worktree and
all prior packages unchanged. After completion report the exact candidate identity
with **“Sequoia tested; High Sierra and video qualification pending.”** Stop for
review before publication or M13. Combined post-M13 qualification precedes any
request to integrate into main.

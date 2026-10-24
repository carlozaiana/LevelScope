# Fixed video paint comparison — 24 October 2026

All twelve agreed runs passed: three preserved M6 controls and three instrumented
runs on each Mac, each lasting 100 seconds, in the fixed alternating order.
Sequoia completed first, then actual Intel High Sierra. No extra playback runs or
retries were added. Both sequences are finished and their audio devices closed.

This comparison did not reproduce the earlier intermittent paint skips. It does
not establish a fix, erase earlier failures, or qualify M7. M7 remains frozen,
unqualified and unpublished. No playback scheduling correction was made.

## Diagnostic publication scope

The user authorized committing and pushing this diagnostic code, evidence and
handover to `codex/video-implementation`. This publication is limited to the
instrumented test target, trace and runner tests, evidence and documentation. It
does not merge into main, change checkpoint tags, resume M7, rebuild tested
packages or update the VST3 plugin. The production player has only a read-only
friend declaration for the probe; its playback and embedded-audio behaviour are
unchanged.

A future recurrence can provide these observations only when a diagnostic run
is explicitly coordinated and its trace is saved successfully. The normal app
and preserved control do not collect this trace automatically. A complete trace
may locate the delay between selection, repaint and paint, but does not by
itself measure compositor presentation or physical synchronization. This bounded
diagnostic allows at most 100 seconds with routes muted; trace overflow or export
failure invalidates the diagnostic evidence.

Earlier preparation and comparison records retain their original pre-publication
status. `pre-publication-handover.md` in the comparison evidence preserves the
exact document identified by `repository-save.json`; `publication-scope.json`
binds that snapshot to this updated handover. No earlier manifest is rewritten.

## Exact identities and setup

Base: qualified M6 `b57a1ed23cbe7d74afe98ed16ec8fd615db9d171`.
The video worktree remains `/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation`.
Candidate-03 source digest:
`0760e080d2e5daa8c89612029236308e134237a142c879003b53c7da1d3727cf`.
Candidate manifest SHA-256:
`35081f6e31cd02f202162525699ecd87912f4cabbc7177a2919f4c53497cad1c`.
Diagnostic executable SHA-256:
`3fd0fe58bbba684c06ac55acaf1618c318812396e3a121434e54a7dcfde1f93b`.
Frozen M6 executable SHA-256:
`c8efc1a6d0b760206c8f99105e553542a4d32fdb9722dc3e19b2172f60bb36e9`.

The actual reference movie was the unchanged 1080p/25 ProRes 422
`51_SyncTest.mov`, SHA-256
`0870d5ab93d1492f23524ed1fff9319af1a6a9cf71882543f2190929afeb8e24`.
Each host used the same movie and selected-audio fixture for its controls and
diagnostic runs. Sequoia used Mac Studio Speakers and the 48 kHz selected-clock
WAV; Golem used Built-in Output/headphones and the 44.1 kHz WAV. The original media
files were not modified. Package manifests contain exact fixture identities.

All twelve reports confirm muted routes, Source mode, the selected device/rate,
real HAL output timestamps, zero audio underruns, zero blank paints, zero
frame-bound violations and zero existing CPU-paint deadline failures. All
attachment and short-run memory gates passed. The 100-second samples do not
assess the >=120-second memory-growth gate or replace endurance qualification.
Embedded movie audio remained excluded by the unchanged video-only backend.
This was not a new independent embedded-tone isolation test.

The user explicitly approved “Start Sequoia, then Golem” after the proposed
outputs, rates, mute state and hands-off intervals were stated. Small progress
records were read during the sequences; no test-window interaction, setting
changes, builds or package changes were made during them. Sequoia traces were
analyzed locally while Golem ran. Connection/workload snapshots are retained;
this is not a controlled experiment attributing behaviour to Firefox, Screen
Sharing, SSH or a particular OS component.

## Results

Paint gaps below are the original harness's maximum qualified paint-entry
interval after its two-second preroll, not a new acceptance threshold.

| Host | Run | Role | Existing short-run gates | Missed painted frames | Largest paint gap (ms) |
| --- | ---: | --- | --- | ---: | ---: |
| sequoia | 1 | control | Pass | 0 | 39.851 |
| sequoia | 2 | diagnostic | Pass | 0 | 35.713 |
| sequoia | 3 | control | Pass | 0 | 35.321 |
| sequoia | 4 | diagnostic | Pass | 0 | 29.801 |
| sequoia | 5 | control | Pass | 0 | 30.799 |
| sequoia | 6 | diagnostic | Pass | 0 | 26.781 |
| high-sierra | 1 | control | Pass | 0 | 19.970 |
| high-sierra | 2 | diagnostic | Pass | 0 | 19.992 |
| high-sierra | 3 | control | Pass | 0 | 19.410 |
| high-sierra | 4 | diagnostic | Pass | 0 | 19.662 |
| high-sierra | 5 | control | Pass | 0 | 19.411 |
| high-sierra | 6 | diagnostic | Pass | 0 | 19.380 |

## What the traces add

All six instrumented traces are complete with zero buffer loss. Their derived
qualified paint counts and missed-frame counts exactly match the corresponding
reports. During qualified playback they record one transport generation (6) and
one view epoch (4), no missing-picture selections, no desired/selected-frame
mismatches and a queue no larger than eight pictures. Window flags remain 15:
visible, unoccluded, active and key, without minimized or live-resize flags.
The observed display remains constant within each run.

Sequoia's independent display observer had maximum qualified intervals of
17.699–18.750 ms, while main-thread vblank delivery reached 30.016–37.128 ms.
Its largest instrumented paint interval was 35.713 ms and its largest full
harness paint callback was 3.836 ms. The controls' largest interval was 39.851 ms.
Golem's largest instrumented paint interval was 19.992 ms, its largest full
paint callback was 2.089 ms, and its controls' largest interval was 19.970 ms.
Both observed displays reported nominal 60 Hz.

These observations distinguish display-refresh callbacks from UI dispatch and
paint duration. They show scheduling variation on Sequoia in this passing sample,
not the cause of the earlier 49–51 ms missed-paint events. The independent
CVDisplayLink is not a hook into JUCE's own callback; no one-to-one dispatch
latency or compositor-presentation measurement follows from it. Larger maxima
in the unfiltered trace summaries include startup/preroll and must not be
substituted for qualified-playback results. Full paint timing includes status
text but not the complete native drawing/compositor stack.

The three pairs per host, fixed order and added instrumentation do not establish
statistical equivalence or an improvement caused by tracing. Cross-host rates
and hardware differ, so this is not an experiment changing only the OS. Physical
optical/acoustic synchronization remains unmeasured; no listening judgement was
requested for these muted runs.

## Disposition and next boundary

Keep all earlier failed attempts alongside these passes. No production scheduling
change is justified by this comparison alone. M7's remaining physical video
qualification and disposition of the intermittent failure are still outstanding.
If further investigation is authorized, predeclare one additional session or
workload condition and use the same controls and diagnostic observations, or
capture an actual recurrence with a separately identified diagnostic. Do not
repeat until a preferred outcome appears. Native invalidation/draw tracing is a
possible later step if a failing trace leaves that boundary unresolved; do not
silently replace the timer, loosen limits or resume M7 publication.

## Preservation and evidence locations

Both runners verified all candidate artifacts at the end, and each recorded
unchanged before/after audio configuration. A separate final local audit verified
2235 original-checkout files, all 647 artifacts in the preserved M6/M7 packages,
all checkpoint tag objects/targets and the diagnostic source inputs unchanged.
No installed application, VST3, Pro40, driver or system setting was edited.

Raw persistent local results:
- `/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/paint-diagnostic-20261024-01/sequoia-comparison-01`
- `/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/paint-diagnostic-20261024-01/high-sierra-comparison-01`

Original Golem results:
`/Users/fuzz/LevelScope-VideoPaintDiagnostic-2026-10-24/high-sierra-comparison-01`.

The complete repository evidence is saved at
`/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/Tests/Video/results/paint-diagnostic-comparison-2026-10-24`.
It includes both sets of raw reports, losslessly compressed raw traces, exact
hashes, the analysis script, derived statistics and the final preservation audit.
All 70 original bundle files were verified against their manifest. Existing
identical evidence was retained; missing files were added without overwriting
any existing records. This handover is saved in the video worktree's Docs folder.

Two earlier automatic approval reviews timed out before the save started. After
the user explicitly requested the retry, the repository save completed. This was
an approval-service issue, not a test failure. The temporary source bundle at
`/private/tmp/levelscope-paint-comparison-20261024-results` and all raw persistent
results remain intact. Historical evidence that describes the then-pending save
is preserved verbatim; `repository-save.json` records its completion and binds
this updated handover to the original evidence manifest.


Metadata caveat: the frozen runner leaves its preflight `mode` text at its
preview label even for an approved run. Those original records are preserved
unchanged. The coordination note, actual per-run invocations, 100-second hardware
reports, logs and outcomes identify the executed runs. This label does not alter
execution or the acceptance checks.

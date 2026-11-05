# Intel multicore true-peak continuation

5 November 2026. **Published after review of the focused two-platform qualification,
with the inherited native completion pause accepted as a development limitation.**
Main integration remains for separate approval. Branch `codex/intel-multicore-true-peak` in
`/Users/fuzz/Developer/LevelScope`, based on published main
`51e47f300fcc38a386a5a92a7277029175cfab0c` and the M18 checkpoint
`m18-interactive-timeline-2026-11-03`.

## Change and candidate

The ordinary macOS standalone measurement service now requests the existing
bounded per-channel true-peak implementation on Intel as well as Apple silicon.
No interpolation coefficients, numerical operation order, loudness arithmetic,
programme input/processing, transport or video code changed. M18's earlier
corrected Intel diagnostic demonstrated equivalence but did not enable the
ordinary Intel app; this continuation addresses that remaining difference.

The existing core selects at most six lanes, limited by channel count and logical
CPU count minus one. It uses the synchronous system dispatch pool and borrows one
4096-frame input block. Each channel retains its own ordered interpolation
history and peak maxima; loudness reduction and drain remain serial. It restores
each worker's floating-point environment. Mono, 5.0, 7.1, insufficient/unknown
CPU capacity, non-concurrent observers, small blocks and failed environment
capture retain serial execution. The CMake option can disable the optimisation.

| Frozen identity | Value |
| --- | --- |
| Implementation / original test-driver commit | `f6fe99716618e94a3a3f5f9c066b5961a06a6075` |
| Candidate | `candidate-01` |
| Source-manifest SHA-256 | `3733d4be90f08cb82291872057ededc4f3c794534492be05cdd85bc4a7f728fb` |
| Measurement build | `77b7a0fa28402b43580096fa7e828af8bd08e681ddca29d4d36fe1e0183b51f6` |
| Universal app SHA-256 | `7833642804d2b4199a31411e91a02c7067f05fb6f42168072e574d19c67f59e0` |
| Local app | `build/intel-peaks-2026-11-05/candidate-01/LevelScope Standalone.app` |
| Golem app | `/Users/fuzz/LevelScope-IntelPeaks-qualification-2026-11-05/candidate-01/LevelScope Standalone.app` |
| Architectures / Intel minimum | arm64 + x86_64 / macOS 10.13 |

The immutable candidate records freeze-time `qualification_complete: false`.
Later receipts establish the outcome without rewriting it. The separate native
UI diagnostic records its own source/binary hashes; changes to that test launcher
do not replace the candidate app or its product sources. The same-build serial
oracle changes only the service compile-time execution switch, and is never
presented as the ordinary application.

## Explicit qualification selection

[Impact plan and budgets](../../Tests/Qualification/IntelPeaks/README.md).
Use the direct core parity oracle and short measurement/section-proposal suites;
then compare matched 60-second 48 kHz PCM24 stereo/5.1 workflows. Source analysis
includes waveform/M/S observation. Section Propose includes a three-pass fresh
candidate and a two-pass reuse of exact peak attribution. Explicit Apply is
followed by independent Current measurement. Seven decisions cover global-only,
section-ready, preference-blocked, measured rejection, unsupported requirements,
reuse and already-meets behavior, with stale-proposal rejection.

All measurement statuses, windows/gates, per-channel peaks and positions enter
the exact snapshot digest. Decisions, narratives, section-peak maps and candidate
plans must match within each host without numerical tolerance. Supported rates,
shuffled roles, boundary samples, ties, silence, padding, observer failures,
non-finite rejection and all four rounding modes receive fresh direct-core checks.

Native checks use the actual main component and ordinary service objects with
driver-supplied file selection; they are not claims of manual OS-picker operation.
They exercise cancel/retry for Source, section Propose and independent Current,
including preservation of applied-but-unverified gains. The final ordinary-app
stack receipt on Golem separately proves the deployed app's actual parallel path.

Reuse M18's pinned EBU/ITU, processing, source-format, export, video, cross-Mac
sessions and three-hour evidence only for unchanged behavior. The source decoder,
measurement arithmetic and sample processing are unchanged and fresh exact serial
oracles cover the scheduling change. The mono three-hour fixture cannot exercise
this optimisation and is not rerun or claimed as new multichannel performance
evidence. No physical audio device is opened in this qualification.

## Driver failures and supplemental evidence

The first Sequoia native pilot passed functional checks but its timer included
picker/setup/screenshot work during concurrent compilation (p95 170.87 ms,
maximum 433.83 ms). A second pilot scoped the timer to active jobs and passed.
The frozen native driver then exceeded its newly selected 200 ms maximum gap
on High Sierra stereo and Sequoia 5.1. These raw failed results are retained.

A separately identified diagnostic records every gap by job stage and compares
the identical native workload with the serial service. It separates actual
periodic timer delivery from the final synchronous completion tail, rather than
changing or deleting the original measurements. The final reviewed decision
includes these results and ordinary-app execution; the initial failed runner
summary alone is not a passing qualification record.

The matched comparison established an inherited completion pause, not a new
worker scheduling regression. High Sierra's largest native gap was 822/840 ms
for serial/parallel stereo and 827/803 ms for serial/parallel 5.1. Native p95
was 106.4/108.0 ms and 103.2/106.6 ms respectively. Sequoia's full native p95 was
103.5/109.6 ms stereo and 102.5/125.4 ms 5.1; maxima were 237/225 and 165/229 ms.
The final native comparison review allows one 100 ms timer period of maximum-gap
variation and 30 ms p95 variation against the same serial workload, with both
maximum and cancellation below one second. This is a documented review of the
original absolute-budget failure, not a rewrite of that test as passing.

The isolated service/timeline heartbeat gates remain unchanged and passed:
parallel p95 104.8/105.0 ms on Sequoia and 109.3/109.4 ms on Golem, maximum below
111 ms across stereo/5.1. Actual cancellation in the matched native parallel runs
was at most 397 ms on Golem; applied Current remained unverified until an
independent successful retry. No rollback or automatic Apply occurred.

## Completed results and practical timings

Both hosts passed fresh core bit-exact parity, 185 measurement checks and 366
section-proposal smoke checks. Complete measured snapshot/peak-position digests,
section maps, plans, explanations and all seven decisions match the serial oracle
exactly for both layouts. Native Source, section Propose and Current cancellation,
retry and verification checks passed. The revised UI review is described above.

These are medians of three alternating, warm-cache runs of the **same 60-second
48 kHz PCM24 input**, with fresh private verified copies and timeline/M/S
observation. They measure service jobs rather than file-picker or click latency.
The first separate seven-decision comparison is excluded from the medians.

| Job | Golem stereo serial → parallel | Golem 5.1 serial → parallel | Sequoia stereo serial → parallel | Sequoia 5.1 serial → parallel |
| --- | ---: | ---: | ---: | ---: |
| Source analysis | 2.609 → 1.555 s | 7.756 → 2.802 s | 1.393 → 0.855 s | 4.179 → 1.336 s |
| Section Propose, fresh attribution / three passes | 8.112 → 4.695 s | 23.741 → 7.819 s | 4.078 → 2.402 s | 12.337 → 4.395 s |
| Section Propose, reused attribution / two passes | 4.948 → 2.841 s | 14.435 → 4.574 s | 2.615 → 1.532 s | 7.794 → 2.085 s |
| Independent Current remeasurement | 2.662 → 1.605 s | 7.829 → 2.874 s | 1.411 → 0.871 s | 4.222 → 1.353 s |

Golem is an Intel i7-3770K, 8 logical CPUs, 32 GiB, macOS 10.13.6 build 17G66.
Sequoia is Apple M1 Max, arm64, macOS 15.7.9 build 24G830. The universal binary
has an x86_64 minimum of 10.13. Actual Golem runtime, not Rosetta, supplies the
Intel evidence. The stereo benefit is about 1.7×; 5.1 Source/Current about 2.7×,
fresh section Propose 3.0×. Other media, rates, storage and machines can differ.
Median process CPU-time/wall-time during Source was 1.75 / 3.38 cores for
stereo/5.1 on Golem and 1.71 / 3.22 on Sequoia; these whole-job averages include
the serial portions and are not a count of simultaneously active peak lanes.

Computer-driven native import in the **ordinary frozen app on Golem** confirmed
L/R and L/R/C/LFE/Ls/Rs, then completed Source and initial Current. Both retained
0 dB processing, 60-second duration, displayed channel peaks and device-closed,
muted routes. The app's sampled stacks show `StandaloneMainComponent::beginMeasurement`
through `OfflineAnalysisService::analyze` into `pushChannelPeaks`, plus that peak
callback on separate system dispatch workers for each layout. These samples are
execution proof only, not timing trials. The disposable app session was closed.

All 361 frozen artifact hashes were rechecked on each Mac. Earlier packages,
all checkpoint tag objects/targets, main, installations and video worktrees are
unchanged. No VST3 build/install, system-audio change or physical playback occurred.
No new long-duration multichannel claim is made. The pre-existing brief completion
pause and M18 cosmetic explanation-text issue remain outside this focused fix.

[Reviewed evidence and historical failures](../../Tests/Reference/results/intel-peaks-2026-11-05/README.md)
include exact commands, raw timing repetitions, both actual OS/CPU identities,
ordinary-app traces, frozen/source hashes and the separate diagnostic provenance.

## Reproduction

Configure a fresh Release tree with universal architectures, deployment target
10.13, JUCE 7.0.12 from the established dependency checkout, and only
`LEVELSCOPE_BUILD_MEASUREMENT_TESTS=ON` plus
`LEVELSCOPE_BUILD_SECTION_PROPOSAL_TESTS=ON`. Build these named targets:

```sh
cmake --build <build>/compile --target LevelScopeStandalone \
  LevelScopeMeasurementTests LevelScopeSectionProposalTests \
  LevelScopePeakChannelParity -j 8
python3 Tests/Profiling/OfflineStages/build_probe.py <build>/compile \
  <build>/drivers/peaks --compare-peaks --universal
python3 Tests/Qualification/M18/build_driver.py <build>/compile \
  Tests/Qualification/IntelPeaks/Native.cpp <build>/drivers/native --private
python3 Tests/Profiling/OfflineStages/fixtures.py <build>/fixtures
python3 Tests/Qualification/IntelPeaks/package.py <build> <build>/candidate-01
python3 Tests/Qualification/IntelPeaks/run.py <build>/candidate-01 <fresh-results>
```

The native supplemental pair is built with
`Tests/Qualification/IntelPeaks/build_native.py <build> <fresh-diagnostic-dir>`.
Its ordinary variant links the app's existing objects; its serial variant adds
the separately compiled service oracle. Use the frozen package's original
scripts for historical reproduction and the recorded supplemental commands for
the corrected diagnostic. Never build an unnamed default/all target or install
the VST3 plugin for these checks.

Compatible older sessions retain settings and references; the conservative
measurement-build binding requires remeasurement when its identity differs.
This optimisation changes throughput, not profiles, acceptance allowances or
the independence of Current/artifact verification.

[Next proposed logic milestone — LRA passages and measured candidates](66_M19LraPassagesAndMeasuredCandidates.md)
remains a proposal only.

## Publication and follow-up

Published checkpoint [intel-multicore-true-peak-2026-11-05](https://github.com/carlozaiana/LevelScope/tree/intel-multicore-true-peak-2026-11-05)
targets `c4eea6f1c4569ec4f31725694ddf34db8340f15e`. Its tested product remains
`f6fe99716618e94a3a3f5f9c066b5961a06a6075`; the later commit adds only qualification
support, evidence and documentation. The immutable candidate is unchanged.
[Publication verification](../../Tests/Reference/results/intel-peaks-publication-2026-11-05/README.md)
records the remote tag object/target, unchanged main and every earlier tag.

The user accepted the inherited approximately 0.8-second native completion pause
for this optimisation. Numerical and scheduling checks passed; the original
absolute 200 ms native gate did not. Preserve that distinction. Track a separate
performance follow-up to profile completion/acceptance/UI work, measure how it
scales with programme length and full-resolution data, and bound message-thread
work. Do not treat this checkpoint as a responsiveness fix or silently rewrite
the original failed receipts. No further long suite was run for publication.

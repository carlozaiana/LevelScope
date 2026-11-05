# Intel multicore true-peak continuation — reviewed qualification

5 November 2026. **Focused qualification complete on actual High Sierra and
arm64 Sequoia. Local continuation only; publication and main integration await
review.** The native completion-pause limitation below remains visible.

[Full report and timings](../../../../Docs/ProjectIndex/65_IntelMulticoreTruePeakQualification.md)
and [explicit impact plan](../../../Qualification/IntelPeaks/README.md).

## Frozen identity and custody

- Branch: `codex/intel-multicore-true-peak`, based on main
  `51e47f300fcc38a386a5a92a7277029175cfab0c`.
- Tested implementation / original driver commit:
  `f6fe99716618e94a3a3f5f9c066b5961a06a6075`.
- Candidate-01 source manifest:
  `3733d4be90f08cb82291872057ededc4f3c794534492be05cdd85bc4a7f728fb`.
- Measurement build:
  `77b7a0fa28402b43580096fa7e828af8bd08e681ddca29d4d36fe1e0183b51f6`.
- Universal app executable:
  `7833642804d2b4199a31411e91a02c7067f05fb6f42168072e574d19c67f59e0`.
- Local preserved app:
  `/Users/fuzz/Developer/LevelScope/build/intel-peaks-2026-11-05/candidate-01/LevelScope Standalone.app`.
- Golem preserved app:
  `/Users/fuzz/LevelScope-IntelPeaks-qualification-2026-11-05/candidate-01/LevelScope Standalone.app`.

`candidate-01-identity.json` is an unchanged copy of the original package identity,
including every source and artifact hash. Its freeze-time
`qualification_complete: false` is historical, not edited after tests.
All 361 artifacts were checked on both Macs. `golem-integrity.json` records the
remote check; `reviewed-summary.json` requires both local and remote integrity.

The later commit containing this evidence changes only qualification drivers and
documentation; it does not change product code or substitute a rebuilt binary.
`supplemental-driver-provenance.json` binds the diagnostic driver sources and its
two binaries separately. `evidence-sha256.json` hashes this archived evidence.
No raw audio, app binaries, private connection material or generated builds are
committed here. Native screenshots are test evidence.

## Fresh required evidence

`reviewed-summary.json` is the final reviewed result. It requires:

1. `sequoia-01` and `high-sierra-01`: actual arm64/x86_64 host identity, fresh core
   bit-exact parity, 185 measurement checks and 366 section-proposal smoke checks.
2. Complete serial/app comparisons for stereo and 5.1: all result bits, windows,
   gates, peaks and positions; peak attribution, proposal decisions, explanations
   and plans; independent Current after explicit Apply. Golem's missing 5.1
   comparison group is supplied by `high-sierra-surround-02`, explicitly selected
   after the original native stop. It does not replace passed earlier work.
3. Three alternating timed repetitions per path/layout on identical 60-second,
   48 kHz PCM24 fixtures, excluding the separate seven-decision comparison.
   Source, fresh/reused section Propose and Current timings are in the summary;
   all repetitions, commands, executable and report hashes remain in receipts.
4. `*-native-diagnostic-02`: same native component with app-mode and private
   serial-service controls, both layouts. Source/Propose/Current cancel and retry,
   applied-but-unverified state, independent verification and closed device.
5. `ordinary-app/*-sample.txt`: sampled native stacks from the unchanged ordinary
   candidate app on Golem, showing `beginMeasurement`, `pushChannelPeaks` and
   separate dispatch worker execution for stereo and 5.1. Thread IDs establish
   different workers over the sample, not a simultaneous lane-count claim.

Golem: Intel i7-3770K, 8 logical CPUs, 32 GiB, macOS 10.13.6 (17G66).
Sequoia: Apple M1 Max, arm64, macOS 15.7.9 (24G830). Compiler, JUCE revision,
universal architectures and Intel 10.13 minimum are in the frozen identity and
build records. Timings are warm-cache short-workflow results, not three-hour
throughput claims.

## Preserved failures and responsiveness review

The first pilot included setup/picker/screenshot work in its heartbeat; pilot-02
scoped sampling to active jobs. Both receipts are preserved.

The original `high-sierra-01/summary.json` and `sequoia-01/summary.json` remain
**failed** because their frozen native driver exceeded an absolute 200 ms gap.
They must not be cited as wholly passing suites. Supplemental matched controls
show the same completion pause in the serial path: High Sierra serial/app maxima
822/840 ms stereo and 827/803 ms 5.1. The reviewed native criterion allows at most
30 ms p95 growth and one 100 ms timer period of maximum-gap growth, with maximum
gap and cancellation below one second. All matched cases meet this review.
The original absolute budget is still not met; the pause is not fixed here.

Isolated service/timeline heartbeat gates passed unchanged (parallel maximum
under 111 ms). Native parallel cancellation on Golem was at most 397 ms.
`review.py` performs the explicit combined selection; missing required groups,
an altered package, failed functional checks or missing ordinary-app proof cannot
produce a passing reviewed result.

## Reused evidence and limits

Reuse [qualified M18](../m18-qualification-2026-11-03/README.md), product
`3665815d73fbd93c0574cdfee8d7efdd805ab7b6`, for unchanged interpolation arithmetic,
pinned EBU/ITU references, input/processing, sessions, export, video and long-work
contracts. M18's corrected Intel diagnostic was supporting evidence, not proof
that its ordinary Intel app used workers. This continuation supplies that proof
and fresh exact within-host serial comparisons. No new cross-host bit-equality
claim is made.

No video loaded/decoded/played and no physical device opened. The old three-hour
fixture is mono and does not exercise this optimisation. Mono, 5.0, 7.1 and other
unsuitable cases retain the serial fallback. Main, all pre-existing tag objects
and targets (`tags-before.txt` / `tags-after.txt`), other worktrees, installations
and frozen packages remain unchanged. M18's cosmetic explanation-text issue and
the inherited completion pause remain limitations. No new LRA method is
implemented; [M19](../../../../Docs/ProjectIndex/66_M19LraPassagesAndMeasuredCandidates.md)
is a proposal only.

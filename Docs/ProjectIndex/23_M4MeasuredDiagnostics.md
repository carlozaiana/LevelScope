# M4 — measured programme diagnostics and read-only timeline

Approved 18 October 2026; the user's subsequent extension makes **three hours** the
representative long-programme qualification case. Development uses the existing
`/Users/fuzz/Developer/LevelScope` working copy on `codex/logic-module`, based on
qualified M3/main `62ab420cb15b3ce5e1c821c963a749a9a56ae4a7`. Main was fetched before
implementation and contained no additional changes. All six earlier checkpoint
tags and the video worktree are preserved.

## Status

**Qualified for the bounded M4 scope on actual Intel High Sierra and arm64 Sequoia,
including three-hour programmes.** Automated, regression, official-reference,
performance and silent UI qualification is complete. The
[qualification report](../../Tests/Reference/results/m4-qualification-2026-10-18/README.md)
records exact coverage and evidence. Manual rebuilds completed before Cancel could
be clicked; cancellation is instead verified by deterministic tests of the actual
native panel actions and separate active-worker tests on both platforms.

Final candidate-03 source digest:
`59bc7a5f17cbcce63c86254e72dbba21ae280d475a19b801bc45c9089d20d9cc`.
Standalone executable SHA-256:
`aeb88e145945a305fec9be3be962ea888d127baae63f37c0cb8b2e7e3efb96d3`.
Candidate-02 interactive evidence applies to the byte-identical candidate-03 app;
candidate-03 adds test coverage only. The source manifest identifies the tested
working tree; the M3 base commit alone does not identify the M4 implementation.

The published annotated checkpoint is
[`m4-measured-diagnostics-2026-10-18`](https://github.com/carlozaiana/LevelScope/tree/m4-measured-diagnostics-2026-10-18).
Its peeled commit is the exact accepted M4 revision and the main revision at
publication; the annotation records the full SHA and tested source/binary hashes.
Main incorporates this checkpoint by fast-forward from qualified M3, without
integration changes to the tested implementation. The established logic working
copy remains on `codex/logic-module`. All six earlier checkpoint tags and the
registered video working copy are preserved.

The qualification manifests describe the pre-commit candidate. Their base commit
and working-tree descriptions remain historical evidence; publication adds only
documentation to those tested source files. Resolve the exact published revision
with `git rev-parse 'm4-measured-diagnostics-2026-10-18^{commit}'`.

## Scope and interpretation

M4 derives a versioned, profile-independent `ProgramAnalysis` from an immutable
measured Source or Current State snapshot. It provides M/S timeline samples,
exact gate inclusion, integrated block-energy shares/ranking, the LRA population,
histogram/CDF, percentile-tail membership and per-channel maximum peak locations.
No PCM is decoded for diagnostics. M2 supplies separately bound profile annotations.
M3 predictions are never accepted as measurements, and diagnostics cannot modify
measurements, evaluations, proposals, gain, session jobs or transport.

Integrated shares are each included overlapping block's energy divided by the sum
of included block energies. They do not represent independent passage energy,
additive passage durations or the causal effect of removing/editing a passage.
The lower/upper LRA lists contain all included windows at/below P10 and at/above
P95. Ties are retained, so lists need not contain exactly 10%/5% of the population.
The current EBU profile has no LRA ceiling. Descriptions do not infer dialogue,
scene type, excessive dynamics or a processing recommendation.

## Time and availability contract

- All locations use signed audio-source sample frames. Ranked windows show the
  exact half-open `[start,end)` interval in frames and six-decimal source seconds.
  M spans 400 ms; S/LRA spans 3 s. Timeline points sit at the **window end**.
- The common programme axis is `[0, source frame count]`; video offset never enters
  diagnostic coordinates. Zoom, scrolling, selecting and inspection do not seek
  or play. Numeric inspection finds the nearest sampled window end, with an
  earlier-window tie break; it displays full-resolution data.
- Ordinary M/S curves are sampled at 100 ms. M1's sample-resolution M/S maxima are
  retained and shown separately, since maxima may fall between plotted samples.
- The separate LRA distribution uses M1's LRA windows, including its EOF drain.
  Rows extending beyond EOF retain their true interval and an explicit padding
  label. Highlighting clips to source duration; padding never extends the timeline.
- Sample-peak positions remain integer frames; true-peak locations retain the
  kernel's eighth-frame coordinate. Boundary interpolation outside source coverage
  is labelled; the full position remains inspectable without extending duration.
  These are one maximum per channel, not a complete event/density map.
- Metric statuses remain available/silence/insufficient-duration/below-gate.
  Zero-energy windows have no finite LUFS plot value. Excluded finite windows keep
  their measured values and exclusion status. Unavailable never means zero.
- LRA's existing under-60-second instability indication is retained separately
  from availability. It introduces no compliance rule or new delivery claim.

## Architecture and identity

`Core/Measurement/GatePopulation.h` extracts M1's existing energy-domain gate
expressions, comparison operators, sum order and percentile-index convention.
Both M1 and diagnostics use these helpers; no displayed gate is converted back
into an energy threshold. M1's algorithm/version remains unchanged. Because the
shared code was touched, both official suites were rerun and passed on both actual OSes.

`Core/Proposal/AnalysisEngine` validates full-resolution populations against the
measured scalar values/counts/percentiles. `AttributionEngine` builds deterministic
rankings and descriptions. `Model/ProgramAnalysis` binds the canonical measured
snapshot digest plus `measured-diagnostics-m4.1`. Equal energy ranks break ties by
source-window end. Descriptive histogram bins are 1 LU; percentile calculations
never use bins. CDF display values at bin edges count the exact eligible population.

`StandaloneDiagnostics` validates the exact snapshot/context through the existing
M2 bridge, then builds off the UI thread. Cancellation/generation checks reject
stale completion. A changed Current processing plan invalidates Current analysis;
Source remains valid. Profile/applicability changes replace annotations only;
video and monitor/routing changes preserve audio diagnostic identities. The
worker's lifetime holds only immutable snapshot data, never UI/session references.

The minimal standalone diagnostics window contains a reusable `MeasurementTimeline`,
state selector, zoom/scroll/numeric inspection, distinct LRA distribution and a
virtualized full-resolution ranked list. The panel receives a **const** session;
the timeline receives measured analysis only, with no transport/device reference.
Closing it cancels diagnostic work. Build diagnostics explicitly retries cancellation.
No audio device is opened by this feature.

## Long-programme performance contract

Qualification uses a generated **three-hour, six-channel, 192 kHz measurement-sized
snapshot** with 100 ms window cadence, loud spikes, excluded quiet windows and EOF
padding. This is a downstream diagnostics stress fixture. A separate generated three-hour,
48 kHz stereo PCM file was also imported and measured through the actual standalone
worker on both platforms, then inspected to its final second. Neither fixture is
normative reference audio; the six-channel/192 kHz case benchmarks downstream
diagnostics, while the stereo PCM case exercises the complete import/measurement/UI path.

On actual High Sierra and arm64 Sequoia:

| Operation / resource | Acceptance limit |
| --- | --- |
| Core analysis and complete worker including snapshot identity | Each <=2,000 ms |
| Viewport preparation, 120 varied pans/zooms | 95th percentile <=50 ms |
| Offscreen timeline paint, 40 views | 95th percentile <=100 ms |
| Additional analysis storage per snapshot | <=64 MiB |
| Additional process peak RSS during diagnostic benchmark | <=128 MiB |
| Active diagnostic cancellation | Settled within 1,000 ms |

Plot storage is limited to 4,096 columns per curve. Each column keeps separate
included/excluded/ungated minima and maxima with counts, preserving narrow extrema
and mixed inclusion. Missing samples add no point; no averages or interpolated
connections imply missing data. Detailed values, ranks and statistics always use
full resolution. List rows are virtualized. Worker input rejects malformed data
and programmes beyond a defensive 24-hour ceiling; M4 performance qualification
covers through three hours, without restricting the pre-existing measurement API.

## Acceptance and reproduction

Build only `LevelScopeStandalone`, `LevelScopeDiagnosticsTests`,
`LevelScopeProposalTests`, `LevelScopeComplianceTests`, `LevelScopeMeasurementTests`,
`LevelScopeScaffoldTests` and `LevelScopeVideoTests`. Use a fresh `build/logic-m4`,
Release, universal x86_64/arm64, minimum deployment 10.13, and pinned JUCE 7.0.12
(`4f43011b96eb0636104cb3e433894cda98243626`). Enable the corresponding
`LEVELSCOPE_BUILD_*_TESTS` CMake options including `DIAGNOSTICS`.

From the repository root (use a fresh build-directory name for a new candidate):

```sh
cmake -S . -B build/logic-m4 -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_OSX_ARCHITECTURES=x86_64;arm64' -DCMAKE_OSX_DEPLOYMENT_TARGET=10.13 \
  -DFETCHCONTENT_SOURCE_DIR_JUCE=/Users/fuzz/Developer/LevelScope/build/juce-7.0.12 \
  -DLEVELSCOPE_BUILD_DIAGNOSTICS_TESTS=ON -DLEVELSCOPE_BUILD_PROPOSAL_TESTS=ON \
  -DLEVELSCOPE_BUILD_COMPLIANCE_TESTS=ON -DLEVELSCOPE_BUILD_MEASUREMENT_TESTS=ON \
  -DLEVELSCOPE_BUILD_SCAFFOLD_TESTS=ON -DLEVELSCOPE_BUILD_VIDEO_TESTS=ON
cmake --build build/logic-m4 --target LevelScopeStandalone LevelScopeDiagnosticsTests \
  LevelScopeProposalTests LevelScopeComplianceTests LevelScopeMeasurementTests \
  LevelScopeScaffoldTests LevelScopeVideoTests -j 6
ctest --test-dir build/logic-m4 --output-on-failure
```

The local JUCE path must contain the pinned revision above. Do not use a default
all-target build or an install target. These commands build standalone/test targets
only; they do not build or install VST3.

Generate the new 120-second, short and silent UI inputs with
`python3 Tests/generate_m4_ui_fixtures.py --output <new-folder> --prior-fixtures <M3-fixture-folder>`.
Add `--three-hour` to generate the streaming-written stereo PCM UI fixture
(2,073,600,044 bytes; excluded from Git). The six-channel/192 kHz measurement-sized
stress fixture is generated directly in `DiagnosticsTests.cpp`; it also checks
end-of-programme true-peak positions beyond 32-bit eighth-frame coordinates.

`package_m4_qualification.py` pins exact source, binaries, runners, reference and UI
fixtures. `run_m4_qualification.py` guards actual OS/architecture and artifact hashes,
runs M4/M3/M2/measurement/scaffold/video/runner tests and both official suites.
`compare_m4_reports.py` retains M3/M2 comparisons and requires identical check
outcomes plus all three-hour performance limits; display rounding never widens a
profile rule. New output directories are required; preserve all previous candidates.

Run the packaged suite **from its candidate directory**, whose `identity.json`
pins the executables and inputs, using an absolute fresh output path:

```sh
cd <candidate-directory>
python3 -B Tests/run_m4_qualification.py --platform sequoia --output <new-absolute-output-directory>
```

On Golem, use `/usr/local/bin/python3` and `--platform high-sierra`. Running the
Intel slice on Sequoia does not qualify High Sierra. Compare the two resulting
report directories with `Tests/compare_m4_reports.py` (see `--help` for arguments).

Silent UI checklist on both actual OSes, device closed/routes Muted/stopped:

1. Import/confirm roles and measure Source/Current. Open diagnostics and verify the
   selector, measured statuses, gates and independent LRA distribution.
2. Zoom, scroll, inspect source seconds, select integrated/tail/excluded windows;
   verify exact intervals/highlights, stable duration and padding labels.
3. Inspect per-channel SP/TP maximum positions and profile-dependent headroom.
4. Change profile/applicability; retain descriptive data, update only annotations.
5. Change Current gain: immediate stale-view removal; remeasure and rebuild. M3
   draft/no-op and evaluations remain unchanged by diagnostic navigation.
6. Change reference-video offset/routing/monitor gain while silent: preserve
   positions, selections, measurements and diagnostic identity.
7. Silence/short programme: missing values and LRA stability remain explicit.
8. Replace/clear source and close/reopen the window; no stale publication, automatic
   playback, gain change or device opening. The actual Cancel and close action
   bindings are also exercised deterministically by `LevelScopeDiagnosticsTests`
   while completion is pending. Three-hour rebuilds finish before sequential
   manual/remote clicks can reach Cancel; a completion-before-cancel attempt is
   recorded as such, never as an observed interactive cancellation. Active-worker
   cancellation is tested separately against the one-second deadline.

## Deferred work

Scene segmentation/editing/trims, gain automation, dynamics, T0–T3 placement,
automatic correction, freezes, processed audition, persistence, public export,
embedded soundtrack and VST3 measurement parity remain outside M4. Existing
reference-picture synchronization and explicit selected-audio authority are
preserved. No new audible, device or physical synchronization claim is made.

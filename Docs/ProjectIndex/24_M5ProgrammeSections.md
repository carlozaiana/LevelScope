# M5 — explicit programme sections and measured section diagnostics

Approved 19 October 2026. Develop and stabilize on arm64 Sequoia first, then
qualify the **same candidate** on actual Intel High Sierra. M5 is not qualified
until both actual-platform automated and silent UI checks pass.

## Status and baseline

**Implemented and qualified on both actual platforms, 20 October 2026.** Sequoia
stabilization preceded High Sierra execution of the identical candidate-01 bundle.
All automated, official-reference, cross-platform and silent UI checks pass,
including full three-hour PCM measurement and final-second inspection on both Macs.
The published annotated checkpoint is
[`m5-programme-sections-2026-10-20`](https://github.com/carlozaiana/LevelScope/tree/m5-programme-sections-2026-10-20),
saved on `codex/logic-module` and fast-forwarded into main from qualified M4.
Its peeled commit is the exact M5 checkpoint and main revision at publication;
the annotation records the full commit SHA and tested source/binary hashes.
Publication and integration introduce no changes to the tested implementation.

See the [qualification report](../../Tests/Reference/results/m5-qualification-2026-10-19/README.md)
and [machine-readable sign-off](../../Tests/Reference/results/m5-qualification-2026-10-19/qualification-summary.json).
Candidate source digest:
`4973465a63bb59ec6c07667f96c6ad579856d89e57c213b5037067be8076d0f2`.
Universal standalone executable SHA-256:
`c5689afab976fc5d04ab4b47208c022b74b9882ec152a3f345441182b3fc02e9`.
The report links exact per-file source/test/binary hashes, compiler/build identity,
OS/CPU/executing architecture, raw results and UI evidence. The base commit below
is ancestry; it alone does not identify the tested M5 implementation. The original
candidate manifest describes the working source before its checkpoint commit.
Its historical base SHA and working-candidate wording are deliberately preserved.
Every one of its 130 source files matches the checkpoint's corresponding file;
publication adds documentation and evidence outside that source digest.
The [publication preflight](../../Tests/Reference/results/m5-qualification-2026-10-19/publication-preflight.json)
records source/artifact integrity and a fresh comparison of the retained platform
reports. Resolve the exact accepted revision with
`git rev-parse 'm5-programme-sections-2026-10-20^{commit}'`.

Each platform passed 2,457 M5 checks, all M1–M4/scaffold/video regressions, four
runner suites, 39 ITU cases and 70 EBU checks. The 24 measured section cases agree
exactly across architectures. Three-hour/1,000-section aggregation took 3.146 ms on
Sequoia and 6.303 ms on High Sierra; interaction p95 was 7.262/20.302 ms, active
cancellation 0.829/1.758 ms, and additional section storage 1,747,195 bytes each.
All M4/M5 budgets pass. No unresolved acceptance issues remain.

Established working copy: `/Users/fuzz/Developer/LevelScope`, `codex/logic-module`.
Fetched main/base: qualified M4 `9a22749901d9f6d4b9c194c9c22c0a231f2494f6`, tagged
`m4-measured-diagnostics-2026-10-18`. Preserve all seven existing checkpoint tags,
the registered video working copy, unrelated Finder metadata and
`Tests/51_SyncTest.mov`. Build/test bundles stay under ignored `build/`.

## Scope and edit semantics

The versioned `programme-sections-m5.1` document begins with one whole-programme
section. Explicit actions split a selected section, move its end/shared boundary,
merge with the next section or rename it. Boundaries are exact signed 64-bit
audio-source sample frames. Rows show `[start,end)` in frames and seconds; the last
end is fixed to source duration. There are no gaps, overlaps or empty sections.
Names are user annotations, not inferred semantic scenes.

Split retains the left ID/name and allocates a new right ID. Move and rename retain
IDs; merge retains the left ID/name and retires the right ID. Retired IDs are never
reused within the document. Rejected edits are transactional; no-op moves/names do
not increment revision. Numeric editing accepts integer frames, never rounded
display seconds. Names must be non-empty and at most 128 characters. A defensive
10,000-section bound limits malformed input; performance qualification uses
1,000 sections across three hours.

Open **Measured diagnostics → Programme sections…**. The window has a Source/Current
selector, virtualized section list, numerical edit controls, read-only M/S timeline
with violet boundary markers, measured summaries and a virtualized window inspector.
The LRA distribution remains in the M4 window. Selection/navigation never seeks,
plays or changes audio. Close/reopen retains the in-memory section document;
persistence is deferred.

## Interpretation

Statistics consume full-resolution M4 measured windows; they do not decode audio
or remeasure sections. Primary statistics use only windows wholly contained in a
section. Windows crossing a boundary are listed separately. Windows ending beyond
EOF form a separate padding category, excluded from primary summaries. A crossing
or padded window can overlap multiple sections; counts are not additive durations.

Each M/S/LRA population shows contained/crossing/padded counts, finite/missing counts
and sampled min/median/max when available. Median is the middle finite value or
midpoint of the two middle values. M/LRA retain original programme gate membership;
S is ungated. No local gates are calculated. Tail counts retain the programme's
inclusive P10/P95 membership and ties.

Integrated contribution adds the existing programme-denominated shares of contained
included M windows. It describes overlapping-block energy, not independent section
energy or the causal effect of editing it. Without an eligible integrated programme
population the share is **unavailable**. Zero share is valid for a section with no
included blocks when a programme denominator exists; it is not zero loudness.

Silence and sections without complete windows remain unavailable. Excluded finite
windows keep their values and status. Programme I/LRA and the existing LRA stability
indicator stay labelled as programme results. No section I, section LRA, section
compliance verdict, LRA ceiling or processing recommendation is introduced.

Per-channel maxima are labelled **programme maximum located here**. SP retains
integer frames; TP retains eighth-frame coordinates. A maximum at an internal
boundary belongs to the following section. Out-of-source interpolated TP is not
assigned to a section. Local maxima are not measured; absence of a programme
maximum cannot establish peak safety. Window intervals and padding remain
inspectable without extending programme duration.

## Architecture and identity

- `Core/Proposal/Model/ProgrammeSections.h`: annotation document, versioned map,
  stable IDs and immutable summary model.
- `Core/Proposal/SectionAnalysis`: pure aggregation using indexed full-resolution
  window ranges. It retains spans rather than copies of windows, and has no
  profile, transport or processing dependency.
- `StandaloneSections`: content/role/range binding and independent cancellable
  worker, consuming M4's authoritative measured-snapshot bridge.
- `StandaloneSectionsPanel`: owns the in-memory document and observes a **const**
  session. The reusable timeline only displays boundary annotations.

The map binds source hash, sample rate/frame count, ordered confirmed roles and
full programme range. Source/Current reuse it only for matching inputs. Invalid,
unconfirmed, replaced or cleared input clears the bound document; coordinates
never silently transfer to different media. Summary identity binds the exact map,
including revision/IDs/names, M4 analysis identity and M5 method version. Stale
worker generations cannot publish after edits or selection changes.

Section edits invalidate only section summaries. Audio revisions, M1 measurements,
M2 evaluations, M3 drafts and gain remain intact. Processing changes invalidate
affected Current diagnostics while retaining source-bound sections and Source
truth. Profile/applicability changes update separate annotations only. Video offset,
monitor gain and routing never enter audio coordinates or identities. Predictions
cannot populate measured sections.

## Reproduction and acceptance

Pinned JUCE 7.0.12: `4f43011b96eb0636104cb3e433894cda98243626`.
Build only explicitly named standalone/test targets; never default/all or install.

```sh
cmake -S . -B build/logic-m5 -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_OSX_ARCHITECTURES=x86_64;arm64' -DCMAKE_OSX_DEPLOYMENT_TARGET=10.13 \
  -DFETCHCONTENT_SOURCE_DIR_JUCE=/Users/fuzz/Developer/LevelScope/build/juce-7.0.12 \
  -DLEVELSCOPE_BUILD_SECTIONS_TESTS=ON -DLEVELSCOPE_BUILD_DIAGNOSTICS_TESTS=ON \
  -DLEVELSCOPE_BUILD_PROPOSAL_TESTS=ON -DLEVELSCOPE_BUILD_COMPLIANCE_TESTS=ON \
  -DLEVELSCOPE_BUILD_MEASUREMENT_TESTS=ON -DLEVELSCOPE_BUILD_SCAFFOLD_TESTS=ON \
  -DLEVELSCOPE_BUILD_VIDEO_TESTS=ON
cmake --build build/logic-m5 --target LevelScopeSectionsTests LevelScopeStandalone \
  LevelScopeDiagnosticsTests LevelScopeProposalTests LevelScopeComplianceTests \
  LevelScopeMeasurementTests LevelScopeScaffoldTests LevelScopeVideoTests -j 6
ctest --test-dir build/logic-m5 --output-on-failure
```

Package a fresh candidate with `Tests/package_m5_qualification.py`. From its root,
run `python3 -B Tests/run_m5_qualification.py --platform sequoia --output <fresh-absolute-directory>`.
After stabilization, run the same bundle on Golem with its inspected existing
Python runtime and `--platform high-sierra`. `compare_m5_reports.py` compares both
actual platforms and enforces M4/M5 budgets. Exact source/binary hashes identify
the tested implementation; its base commit alone does not identify working changes.

The retained local fixtures are in
`build/m4-qualification-2026-10-18/candidate-03`. With that directory available,
the isolated bundle can be reproduced with a new output directory:

```sh
python3 -B Tests/package_m5_qualification.py \
  --build build/logic-m5 --output build/m5-reproduction/candidate-01 \
  --m1-bundle build/m4-qualification-2026-10-18/candidate-03 \
  --video-fixtures build/m4-qualification-2026-10-18/candidate-03/video-fixtures \
  --ui-fixtures build/m4-qualification-2026-10-18/candidate-03/ui-fixtures
```

The packager verifies the official vectors against their repository manifests.
Downloaded reference audio and generated fixtures stay outside Git. A fresh build
may have different binary hashes; qualify its own manifest rather than relabelling
it as the accepted binary.

M5 budgets for three-hour, six-channel, 192 kHz measurement-sized input with
1,000 sections: aggregation <=2,000 ms; boundary-edit/selection/pan/paint p95 <=100 ms;
active cancellation <=1,000 ms; additional section output/map/median workspace
<=16 MiB. Retain M4 budgets. This does not claim new PCM throughput or physical
synchronization performance.

Tests cover independent interval oracles, exact endpoints and one-frame sections,
stable IDs/no-ops, transactional errors, ties, finite/missing/gated data,
short/silent real audio at every qualified rate, EOF padding, maximum ownership,
identity reuse/invalidation, stale/cancelled work, native controls and long-programme
performance. Run existing M1–M4/scaffold/video regressions. Measurement code is
unchanged; any later shared-engine change requires affected official suites to
rerun. The packaged runner can rerun both pinned official suites as well.

Cross-platform section comparisons require exact interval, membership and
availability agreement for 24 real measured section cases (four rates, three
signal types, two sections). Sampled min/median/max use a numerical comparison
budget of 1e-9 LU; programme-denominated shares use 1e-12. These compare computation
results, not profile acceptance rules. Both sides also pass independent interval
and sorted-value oracles. The native panel checks all controls fit its minimum
1000 × 820 content area; summaries and window lists scroll independently.

Silent UI checklist on both actual OSes, Device closed/routes Muted/stopped:

1. Import/confirm/measure Source and Current; verify initial whole programme,
   shared boundaries, measured labels and exact audio-source time.
2. Split, move by one frame, rename, merge; inspect highlights and reject invalid
   or empty boundaries. Preserve gain, measurements, verdicts and drafts.
3. Inspect contained/crossing/padded windows, maximum wording, silence, short
   sections and unavailable populations; retain programme LRA interpretation.
4. Change Current gain: clear stale truth, remeasure and rebuild. Preserve Source
   and section coordinates.
5. Change applicability/profile and video offset/monitor settings while silent;
   retain descriptive results and positions.
6. Cancel/retry, close/reopen, replace/clear. Distinguish observed manual behaviour
   from deterministic native-action tests if rebuilding finishes before Cancel.
7. Inspect a three-hour source through its final second, including boundary
   navigation, detailed values and unchanged programme duration.

## Deferred work

Section trims, automatic segmentation/classification, protection enforcement,
freezes, dynamics, T0–T3, proxy simulation/refinement, processed audition,
persistence, export, embedded soundtrack and VST3 parity remain outside M5.
Selected audio remains authoritative; video audio stays muted. Preserve installed
apps/plugins, Golem's prior copies, Pro40 and system audio configuration.
Large PCM imports and measurement snapshot verification can take minutes on High
Sierra. The qualified worker retains its full-content verification; M5's response
budgets concern diagnostics after measurement, not import or measurement throughput.

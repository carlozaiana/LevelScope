# M15 qualification and handover

Status, 29 October 2026: **Sequoia tested; High Sierra, video and three-hour checks pending.**
Candidate-05 is published on `codex/logic-module`. The exact
checkpoint and remote verification are recorded in the [publication receipt](../../Tests/Reference/results/m15-publication-2026-10-29/README.md); main stays at M13.
Nuendo 15.0.30’s actual six-point CSV and 25 fps are confirmed; a second actual export
confirms its two-second origin. The user also confirmed Pro Tools 10.3.10 at 24 fps,
project/audio/trailer start at zero, and the first two visible marker positions. Both
actual export setups are verified within the bounded reader scope.
Published checkpoint: [`m15-marker-import-sequoia-2026-10-29`](https://github.com/carlozaiana/LevelScope/tree/m15-marker-import-sequoia-2026-10-29),
commit `1b56804ecfb469aef3276b2d4ac96039ce3e4cc1`. Main remains qualified M13
`a9be070ab9eb5cb4f591ff9d991d8999c3a016d0`.

## Implementation

[53](53_M15MarkerImportAndSourceTimeMapping.md) and
[54](54_M15ImpactQualificationPlan.md) define the approved boundary. The pure
`Core/Proposal/MarkerImport` module reads Nuendo CSV, the dedicated Pro Tools
Session Info marker table, or pasted positions. It converts explicit units/rate/
origin with checked rational arithmetic and floor-to-sample quantization; it never
loads audio/video. Actual versioned Nuendo exports and a minimal verbatim Pro Tools
10.3.10 marker table supplement the synthetic cases. Audio-frame selection uses
TIME REFERENCE / Samples only with a matching sample-rate header; other units
use LOCATION, with a disclosure of timecode subframe precision loss.
[Export instructions and dialects](../../Tests/Markers/README.md) identify their
bounded support, supplied versions and independently confirmed programme mappings.

The section editor opens an asynchronous, cancellable import preview. Its virtual
list exposes original rows, exact mapped positions/quantization, the full proposed
map, additions/removals/renames and unused-column/duplicate disclosures. Explicit
confirmation atomically replaces annotations after checking runtime, programme,
annotation map and current marker bytes. Active playback or file/measurement/
proposal/export jobs block confirmation. No implicit stop, seek, remeasurement,
proposal generation, processing Apply or device action occurs.

Ordered roles/content/range bind the annotation map. Annotation edits preserve
Source/Current audio identity, applied section intervals and full-precision gain,
compatible measurements, selected-rule evaluations and rendered audio. Geometry-
dependent drafts/proposals retain their explanation but cannot apply when stale.
Reports use current annotations/provenance and renewed report approval, while an
unchanged rendered artifact is reusable. A replaced programme cannot inherit the
previous programme’s annotation receipt in a new report.

Session schema 5 stores the latest receipt (input basename/hash, original selected
rows, mapping options/version, exact samples, imported map/IDs and disclosures).
Manual annotation changes preserve it as historical provenance and disclose map
divergence. Schemas 1–4 remain readable; opening never rewrites the original. Save
As preserves the older-readable copy. Schema 5 requires an M15 reader. Measurement
cache format/algorithm are unchanged; the conservative build/platform compatibility
policy remains. New source hashes change the build identity, so older-build caches
are discardable, with settings/media/sections retained and remeasurement required.

## Current candidate-05 evidence

The supplied Pro Tools 10.3.10 file required its actual `MARKERS LISTING` heading
and exact sample-reference selection. Only `MarkerImport.cpp` and the marker-panel
field-policy label changed in production code from candidate-04. All earlier
packages and receipts remain unchanged. No measurement, input decoder, processing,
transport, timing, video, session or report implementation changed in this follow-up.

- Source manifest: `38676a45d47ddba6eb54f7ab299f52f622f6be26554256781b8960749c609564`.
- Measurement/build identity: `8a7a2866d48f7ba6d4db1c7f532eaffbfce03fa78fbca02211d3f6b9fb79ef47`.
- Universal executable SHA-256: `2fc20d4cb9f1db21bc875ba2df9b9ebbeb048976093771affca13c7e25e7dd3c`.
- Preserved app: `build/m15-qualification-2026-10-29/candidate-05/LevelScope Standalone.app`.
- Base commit remains `be4bf23c8bb9f661ae1ede30f3b6e79a067aea37`; it is not the
  frozen M15 implementation. Checkpoint `1b56804ecfb469aef3276b2d4ac96039ce3e4cc1` now records
  all 262 manifest paths exactly; no rebuild or qualification rerun was needed.

All eight explicitly selected short commands ran fresh and passed in **94.28 s**
on actual arm64 Sequoia 15.7.9 (24G830), Apple M1 Max. Marker tests passed **535
checks in ten required groups**, including both actual Nuendo exports, the minimal
verbatim Pro Tools fixture, and the full original 672,421-byte Pro Tools export.
The new actual-Pro-Tools group is mandatory in runner completeness validation.
Missing, duplicate or failed required evidence cannot produce a short-suite pass.
All seven app/test executables contain arm64/x86_64; Intel minimum remains 10.13.
This does not establish actual High Sierra runtime success.

All 45 marker sample references are preserved exactly, including 587,888 and
1,610,400 for the first two points. At explicit 24 fps/origin zero the LOCATION
labels map to 586,000 and 1,610,000 instead. The UI discloses the subframe precision
loss and unused fields; sample selection rejects mismatched or missing rate, mixed
reference units, malformed counts and ambiguous tables. No origin is inferred from
“Start Audio.” Full-file and minimal-fixture selected content agree exactly.

The actual-file UI check used a disposable 900-second silent mono coordinate
fixture, not the user’s programme. It verified both explicit field selections,
preview, disclosures, cancellation, annotation-only replacement, Save As and safe
reopening. The saved core retains all 45 original rows, full-export hash, chosen
units/origin and 46 resulting sections. Reopening remains stopped/device-closed/
muted, with measurements unavailable and Source UNVERIFIED. No audio measurement,
playback, video or physical-device operation was performed in this UI check.
The unchanged earlier Nuendo/paste/stale-draft/report/cache UI cases retain their
original candidate identities and the explicit reuse rationale in `ui-candidate-05.json`.

Resource results: 10,000-row parse/map **0.0445 s**; approximately 8 MiB parse
**0.5645 s**, **51.1 MiB** incremental peak RSS; UI heartbeat p95 **6.38 ms**,
maximum **133.52 ms**. These meet the specified p95/parse/memory budgets; they
are not long-duration audio evidence. The extra silent UI media uses 86.4 MB.

Fresh records: `sequoia-candidate-05/summary.json`, `candidate-05-cross-build.json`,
`candidate-05-final-integrity.json`, `pro-tools-interoperability.json` and
`ui-candidate-05.json`, with text-only copies in the repository evidence folder.
The original user export remains unchanged and outside Git; only its timing
headers and complete marker table are retained as a regression fixture.

On 29 October the user independently confirmed the Pro Tools project is **24 fps**,
project and actual audio start at **00:00:00:00**, and a trailer begins at that
origin. The first two DAW-visible markers are **00:00:12:05 / 00:00:33:13**, matching
the export. “Start Audio” is an internal marker, not programme zero. This confirms
the exact zero-origin mapping already tested on candidate-05; no source change,
rebuild or repeat test is necessary. The reference-sample fields retain subframe
precision that is absent from the displayed LOCATION labels. No general Pro Tools
version/export-dialect or universal Nuendo-offset claim is made.

The approved **Sequoia phase is complete**. The user subsequently authorized and
published this **Sequoia-tested development checkpoint**, with actual High Sierra,
cross-Mac exchange, all video execution and three-hour qualification explicitly
pending. Main integration remains outside that authorization.
Main remains fully qualified M13. The frozen identity’s pending-at-freeze status is
preserved; current acceptance is in `status.json` and this handover, not a rewritten
package. All prior checkpoints and tested packages remain unchanged.

## Historical evidence, candidates 01–04

Historical candidate-04's required short automated selection passed on actual arm64
Sequoia 15.7.9 (24G830), Apple M1 Max. This is **not yet complete M15 Sequoia
qualification**: at that point the required actual Pro Tools example and comparison
were pending. Candidate-04 includes the corrected paste-mode label. Its package and
receipts remain unchanged as candidate-05 adds the actual Pro Tools adapter coverage.

- Candidate-04 source manifest:
  `0cfd128c38e78b74d890eccb6bc78ad52ae3527fda970ec06d9eab4d8e7f7462`.
- Shared measurement/build identity:
  `830e191e8e8431393702e6f4b7e7df5a1b7c89677c72288cfcff82672a2b3608`.
- Universal standalone executable SHA-256:
  `37e83225a3555fc9d484e1bf0cf2f7e62e24d57a6c04d19a2d1f1a32a9fb0fad`.
- Preserved package: `build/m15-qualification-2026-10-29/candidate-04/`.
  Candidates 01, 02 and 03 remain unchanged. The base Git commit is the M14 publication
  receipt, not the uncommitted M15 implementation; the source manifest identifies it.

The explicit selection contains M15 marker, short session, section, M12 processing,
M13 proposal and short export tests, frozen M14 reader rejection, and runner
completeness tests. Candidate-01 ran all eight commands in approximately 97 seconds.
Candidate-02 reran marker/legacy-rejection/runner checks in approximately 9.4 seconds,
reusing the five regression receipts only after verifying identical production
source, build identity and each executable. Those receipts retain candidate-01's
source digest `f140731dde637c51d92426d0c353f73bc5a45c4619b870da3dbf063b4c902278`
and original paths/results. No recursive milestone runner is invoked. The frozen
M14 reader rejects schema 5 with its existing generic `Session integer outside
supported range` message; M15 explicitly warns that schema 5 requires an M15 reader.
All seven frozen app/test executables contain arm64/x86_64 slices, with Intel
minimum macOS 10.13. This is cross-build evidence only.

Candidate-02 (manifest
`8d0f5c7a50d43317fff19874cf9895a363e96fd25dd518d2ff419986903b9acf`)
first received a separate test-only two-second-origin supplement (370 checks),
preserved under `nuendo-origin-check/`. Candidate-03 then froze the named
programme-zero export and complete actual-file regression: **381 checks in nine
named groups passed**. Its fresh commands took approximately 9.3 seconds; the
same five hash-identical regression executables reuse candidate-01 receipts.
No production source, CMake or app binary changed between candidates 01/02/03.
Candidate-03's source and executable hashes are not substituted into earlier records.

After the display-only fix, candidate-04 ran **all eight commands fresh**, with no
regression-receipt substitution, in 96.23 seconds. All passed, including 381 marker
checks in nine named groups. Seven explicitly built app/test targets contain both
architectures and retain Intel minimum 10.13. Measurement/input/processing source
remains unchanged; the conservative source/build identity nevertheless changes.

Candidate-04 measured 10,000-row parse/map at 0.056 seconds, an approximately
8 MiB input at 0.586 seconds with 53.8 MiB incremental peak RSS, and UI heartbeat
p95 6.3 ms (maximum 133.0 ms). These meet the explicit p95, parse and memory budgets;
they are short import/resource checks, not three-hour audio execution.

Native UI on the byte-identical candidate-02 app confirmed safe reopening, pasted
25 fps preview/detail, cancellation, out-of-duration rejection, explicit annotation
replacement, stale draft/no Apply, save blocking until draft discard, Save As and
reopening with the original applied gain/current measurements. Independent saved
session comparison changed only document revision, sections and marker provenance;
measurement-cache payload bytes are identical. The app stayed stopped, device-closed
and muted. The user selected NDF units; keyboard popup navigation then completed
format/rate/map/diff selection. The follow-up UI covered CSV timecode/name import,
full replacement map/diff, dedicated Pro Tools synthetic-table handling and rejection
of CSV as Pro Tools text, fractional-rate quantization and cancellation, followed by
Save As with byte-identical measurement caches. These records retain their original
candidate identities: `ui-candidate-02-partial.json` and
`ui-candidate-03-file-reader-supplement.json`. Synthetic Pro Tools checks do not
establish actual DAW interoperability.

The UI revealed that switching to pasted positions left the previously selected
filename visible. The provenance was correct, but the display was misleading.
Candidate-04's one-line visibility fix hides that label in paste mode; no parser or
audio-processing behavior changed. Fresh UI inspection verified the fix, preview,
explicit confirmation, unsaved-close cancellation, Save As and reopening. The
saved core again changes only revision/sections/provenance, with cache bytes
unchanged and `Pasted positions` correctly recorded as input. Current measurements
and full gain precision remain current; device stays closed and all routes muted.
Unchanged earlier UI cases retain their original candidate receipts with the narrow
reuse rationale in `ui-candidate-04.json` and `ui-label-fix-impact.json`. The corrected
executable was not substituted into the preserved candidate-01/02/03 packages.

All local receipts are under `build/m15-qualification-2026-10-29/`, including
`sequoia-candidate-01/summary.json`, `sequoia-candidate-02/summary.json`,
`sequoia-candidate-03/summary.json`, `candidate-03-cross-build.json`,
`candidate-03-equivalence.json`, `nuendo-origin-check/marker-tests.json`,
`nuendo-interoperability.json` and the UI/supplement records above.
Candidate-04's fresh records add `sequoia-candidate-04/summary.json`,
`candidate-04-cross-build.json`, `candidate-04-final-integrity.json` and
`ui-candidate-04.json`. Bounded text-only evidence is retained in
[the M15 evidence folder](../../Tests/Reference/results/m15-qualification-2026-10-29/),
with a SHA-256 evidence manifest; media, binaries and session files remain local.
The runner rejects missing/duplicated/failed checks and mismatched build identities.
[Current status](../../Tests/Reference/results/m15-qualification-2026-10-29/status.json)
records the completed Sequoia phase after the actual Pro Tools programme mapping
was independently confirmed; deferred qualification remains pending. That original
readiness receipt is deliberately unchanged; the separate publication receipt records
the subsequent checkpoint. No main merge has been performed. All 18 prior tag objects/targets, the frozen M14 package and
M15 packages, original Nuendo files and video branch remain unchanged.

Reused measurement/reference, input/RF64, processing and timing evidence originates
from published M14 candidate-02, source manifest
`9dd22a598281e83431039f210202a8a24fe49517294cc6bb8d83032c4975e013`.
The final candidate records exact unchanged source hashes and reuse rationale.
This is historical evidence for unchanged behavior, not actual High Sierra/video
or three-hour qualification of M15. Fresh tests cover changed session/report/model
callers and marker import; no old result is relabelled as a fresh pass.

## Reproducing the bounded checks

Use the established working copy and pinned local JUCE 7.0.12 dependency. The
Release build uses `arm64;x86_64`, Intel deployment target 10.13 and the enabled
marker, section, session, section-processing, section-proposal and export test
options. Only these explicit targets are built; never run the default all/install
or VST3 target:

```sh
cmake --build build/logic-m15 --target LevelScopeStandalone LevelScopeMarkerTests LevelScopeSectionsTests LevelScopeSessionTests LevelScopeSectionProcessingTests LevelScopeSectionProposalTests LevelScopeExportTests -j 6
```

Run a preserved candidate's own explicit runner into a **new** output directory:

```sh
python3 -B build/m15-qualification-2026-10-29/candidate-05/Tests/run_m15_qualification.py --candidate build/m15-qualification-2026-10-29/candidate-05 --output build/m15-qualification-2026-10-29/recheck-new --pro-tools-original /Users/fuzz/N15/TUTTI_UCCIDONO_V2.txt
```

The runner verifies package hashes before execution and runs only the eight named
short commands. It performs no video, physical-device or network/Golem operations.
JUCE needs permission to create its owned temporary source snapshots. Optional
`--reuse-candidate` / `--reuse-output` requires identical production source/build
and individual regression-executable hashes; it is not used across the candidate-04
or candidate-05 source changes. The optional `--pro-tools-original` checks the full
user file locally; the mandatory minimal actual table is always tested. Interactive UI evidence and actual DAW interoperability are separate
required records; an automated-runner pass never sets complete qualification true.

## Remaining coverage and product limitations

- Nuendo 15.0.30 CSV: received, SHA-256
  `24941511d781d103f801a348b022001f3eda2b33887d694b37a87393f1b835b7`.
  The user confirmed 25 fps, default UTF-8/comma export, project start two seconds
  before programme zero, and two ruler positions exactly two seconds before CSV
  labels. An explicit `00:00:02:00` export-coordinate origin reproduces both;
  a second actual export with a named programme-zero marker confirms that mapping
  and preserves all earlier rows. The Time Offset control value is not separately
  reported, but the anchor establishes the actual export-coordinate origin. Do not
  infer a universal correction or describe raw zero-origin tests as programme alignment.
  The original six-point example has empty names; the second includes a named anchor.
  Pro Tools 10.3.10 actual export and user-confirmed zero origin/24 fps/first two
  visible positions are verified for this setup. Its 45 points retain both
  observed LOCATION labels and exact Samples reference behavior. The manual review and additional-format scope are recorded in [53](53_M15MarkerImportAndSourceTimeMapping.md).
  Synthetic fixtures alone do not qualify an actual DAW version.
- Actual High Sierra, cross-Mac exchange, all video execution and three-hour
  qualification: deferred to the combined window after M18, no later than M19.
  Universal/Intel build success does not establish High Sierra runtime success.
- No new 90-minute workflow or physical-device test is required for M15.
- Point markers only; no native projects, complete track reconstruction, cycle/
  range semantics, drop-frame, inferred 24-hour wrap, VFR frame lookup, marker/EDL
  automatic gain changes, or automatic segmentation.
- Existing M13 proposal policy unchanged. Cinema test: four changed sections,
  89.53% affected envelope duration; the 20% duration gate blocked it, not eight
  changed sections. Preference review is next as a separately approved proposal,
  including marker-density sensitivity and versioned budgets. Standalone limiter
  qualification and basic reviewed measured-audio region suggestions follow, each
  separately proposed/approved. Dynamics-led use must not require markers.
- Preserve accepted occasional picture skips, while audio discontinuity,
  accumulating drift, sustained freezes and identity/lifecycle failures remain
  blockers when deferred tests run. Physical full-7.1 listening remains separate.

M15 checkpoint publication is complete. Stop for review of the separate M16 proposal;
no M16 implementation or main integration is authorized.
Preserve all tags/packages/installations, unrelated work and the video working copy.

## Publication and next review

The checkpoint status remains **“Sequoia tested; High Sierra, video and three-hour
checks pending.”** The original candidate-05 app, package identity and qualification
receipts are unchanged. The new commit ties the manifest to Git; the tested binary
still contains its original measurement/build identity rather than a substituted
post-publication build. All 18 prior tag objects/targets and qualified M13 main were
verified directly on GitHub after push.

[M16 — bounded proposal preferences and configurable limits](56_M16ProposalPreferencesAndLimits.md)
is a proposal only. It covers configurable magnitude/count/duration/reserve limits,
canonical processing-region accounting and reproducible policy identity. No new
processing strategy, neighbour optimizer, segmentation or dynamics is implemented.

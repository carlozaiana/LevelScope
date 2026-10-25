# M8 — reproducible sessions and safe reopening

Implemented on `codex/logic-module`, based on published M7 main
`cfbf3f981dece326d2363b911e5914ea61f1796c`. **Candidate-01 qualified on actual arm64
Sequoia and Intel High Sierra on 25 October 2026; published and integrated into main.**
The annotated checkpoint is
[`m8-reproducible-sessions-2026-10-25`](https://github.com/carlozaiana/LevelScope/tree/m8-reproducible-sessions-2026-10-25).
Its peeled commit is the accepted M8 checkpoint and main revision at publication.
The tag annotation records the full commit SHA, source manifest, binary hashes
and qualification evidence identity. The [publication record](../../Tests/Reference/results/m8-publication-2026-10-25/README.md)
documents that mapping and preservation checks. The approved contract and subsequent
roadmap are in [35](35_PostM7BroadcastWorkflowProposal.md).

## Implemented boundary

Explicit New / Open / Save / Save As, unsaved-change handling and per-file relink
operate on a local `.levelscope` document. Sections now belong to the session
model rather than a diagnostic window. Saved annotations remain independent of
processing, measured audio and compliance. Source and Current gain, role order,
confirmation/provenance, target definition reference, applicability and signed
rational reference-picture offset round-trip without display rounding.

Reopening never starts transport or opens a device. Device selection and physical
routes are not restored; all routes begin muted. The audition domain may be
restored, but an unmeasured/incompatible Current State remains UNVERIFIED. Monitor
gain is not processing and is not a persistent measurement setting. Saved playback
position, transient diagnostic views, proposals and predicted values are omitted.
Reference video remains subordinate to the authoritative audio programme and its
embedded audio is never decoded or routed.

`Source/Standalone/SessionDocument.{h,cpp}` owns serialization, bounded validation,
per-file verification, failure handling and coherent snapshots. The UI controller
in `StandaloneSessionUI.cpp` owns asynchronous file operations and explicit
user choices. It stages a complete replacement before committing it. Invalid core
documents, failed writes and cancelled opens preserve the previous session in a
stopped state. A final media-stamp check rejects changes after verification.

## Document and identity contract

Container version 1 starts with `LSSESS01`, an unsigned little-endian 64-bit core
length and 64 ASCII characters of SHA-256, followed by canonical UTF-8 JSON. The
core is capped at 4 MiB, nesting at 24, and the whole file at 128 MiB. Decimal
strings preserve 64-bit revisions, frames and section IDs; sixteen hexadecimal
digits preserve finite binary64 processing gain exactly. The rational time parser
retains the exact signed video alignment. Unknown core/assembly/processing versions
fail safely. The format is an application document, not an interchangeable JSON
editing interface or a certification signature.

A programme has a persistent ID independent of its file table. Assembly schema 1
records stable file IDs, content identities, locators, ordered file-channel-to-role
assignments, exact alignment and range. M8 implements only one interleaved file,
identity channel ordering and the existing whole-programme range. Multiple files,
mapping transformations, stem sums and analysis-only components are not accepted.
The future assembly identity must be versioned; it cannot inherit a single-file
measurement merely because a filename or channel count matches.

The persistent document UUID/revision is separate from the runtime opening UUID
and job ticket. Every successful open gets a new runtime identity. Previous jobs
cannot publish into it; New/Open wait for owned jobs to finish/cancel. Saving
captures an immutable revision. Subsequent edits or accepted measurement completions
remain UNSAVED when that older revision finishes saving. A late save completion
cannot mark a different opening clean.

Atomic save uses an exclusive, owned temporary in the destination directory,
checked chunked writes, flush and close, then atomic rename. Cancellation/failure
removes only its own partial. The old saved file remains intact until replacement.
Source audio and reference picture cannot be overwritten through paths, symlinks
or hard-link aliases. Exclusive-create collisions preserve the unowned file;
cleanup is enabled only after this operation successfully creates its temporary. No autosave, backup retention or power-loss recovery journal
is claimed.

## Replaceable measurement caches

Source and Current caches have independent framing, domain, size and checksum.
Missing/truncated/corrupt/incompatible caches are discarded with an explanation
when the core is intact. Settings, sections and media references remain recoverable.
Original measurement session/job/revision, method, source build and platform are
preserved; opening does not relabel old measurements as newly performed work.

Full-resolution metrics, availability/status, LRA stability, windows and per-channel
peaks/locations are serialized. Silence and missing values never become numeric
zero. No saved compliance verdict, prediction, diagnostic summary or proposal is
authoritative. Exact content, roles, range, processing, measurement method, build
and platform must match before reuse. Each required file is verified independently.
Unknown target definitions remain identified but unavailable for passing evaluation.
Unusable caches are omitted on the next save; they are not an archival history.
Known profile changes reevaluate compatible measurements without remeasuring audio.

Validation, cache binding and evaluation preparation run off the UI thread. The
model caches only newly derived evaluations of immutable accepted snapshots, keyed
by the current context, applicability and profile digest. A different snapshot or
processing/profile binding forces reevaluation. Prepared adoption repeats the
media-stamp and lifecycle guards, installs once, and is bounded on the UI thread.

Cross-platform/build sessions recover settings and references; conservative cache
policy requires remeasurement. Relinking an identical file changes its locator,
not its content identity. A same-name changed file cannot inherit measurements.
Missing audio leaves an explicitly unresolved programme; missing video leaves
audio authoritative and does not invalidate unchanged audio measurements.

## Qualification results and evidence

Candidate-01 was frozen on 25 October 2026. Its manifest identifies the complete
working source; the M7 base commit is ancestry, not an M8 tested-source claim.

| Identity | Value |
| --- | --- |
| Source manifest SHA-256 | `1c591122c03f40cdc7035d20837e4e39f9cf0a1c62de640149c0b6427d225229` |
| Measurement build | `77cf1c7a5052442d4c938ed82ff1158b2bf6557d28bf99880666ec417f4c7c05` |
| Package identity SHA-256 | `dd625b0e2241afde0e6cbd5bbe36f733d8e23bdd9afc3bba92cd62fc556e8b9f` |
| Universal standalone executable SHA-256 | `146b9a79e12387a5f1c368e1543f50e581bcd6ceea63cde2111ca4e4d9161740` |

Local package: `build/m8-qualification-2026-10-25/candidate-01`.
High Sierra ran the identical package in its own M8 test directory.
Neither package is an installation. Final verification checked all 379 packaged
artifacts on both Macs and all 174 source/test manifest entries in the working tree.
There were no post-freeze implementation changes or integration fixes.

The [qualification evidence index](../../Tests/Reference/results/m8-qualification-2026-10-25/README.md)
links the original reports, UI records, exchange results and comparisons. The
combined `qualification-summary.json` establishes completion; the raw automated
runner summaries deliberately retain `qualification_complete: false` because those
runners alone cannot establish interactive UI or bidirectional exchange.

| Required gate | Sequoia | High Sierra |
| --- | --- | --- |
| Actual platform | 15.7.9 / 24G830; Apple M1 Max; arm64 | 10.13.6 / 17G66; Intel i7-3770K; x86_64 |
| M8 session core, recovery, failure and resource checks | 2,073 assertions passed | 2,073 assertions passed |
| Restored audition and reference picture | 2,431 assertions; 48 captures; 5 video cases passed | 2,431 assertions; 48 captures; 5 video cases passed |
| Foreign-session import and verified relink | High Sierra → Sequoia passed | Sequoia → High Sierra passed |
| Silent interactive session lifecycle | Passed; `sequoia-ui.json` | Passed; `high-sierra-ui.json` |
| M1–M7 and runner regressions | Passed | Passed |
| Official reference vectors | ITU 39/39; EBU 70/70 | ITU 39/39; EBU 70/70 |

All 48 M8 restored capture cases were bit-identical across the two platforms
(maximum float32 ULP difference 0); their metrics, availability and LRA stability
also met the existing comparison rules. Fresh M6 and M7 captured streams were
likewise bit-identical. The composed M1–M5 comparison passed. No physical audio
devices were opened for M8; playback/routing tests used injected devices.

The UI records cover Save As, New/Open, coherent settings and measured-cache reuse,
unverified Current State, section edits, unsaved-change choices, wrong-content
relink rejection, exact-content relink, invalid-core rejection, independent damaged
cache recovery, active long-media cancellation and safe closure. Sequoia also
exercised missing-video recovery and restored rational video seeking interactively;
both platforms ran the five native restored-reference cases. Expected loudness FAIL
verdicts in the test fixtures are successful checks of measured failure reporting.

### Resource results

These are single recorded qualification runs, not throughput guarantees. The
three-hour persistence fixture contains full-resolution Source/Current caches and
1,000 sections; the real-media test separately verifies the complete three-hour,
eight-channel 48 kHz FLAC (518,400,000 frames).

| Measurement | Budget | Sequoia | High Sierra |
| --- | --- | --- | --- |
| Save | ≤5 s | 65.4 ms | 147.8 ms |
| Core/cache parse | ≤5 s | 325.7 ms | 740.7 ms |
| Additional persistence RSS | ≤256 MiB | 67.6 MiB | 44.6 MiB |
| UI callback interval p95 | ≤100 ms | 6.31 ms | 6.05 ms |
| Persistence cancellation | ≤1 s | 0.40 ms | 0.19 ms |
| UI prepared-state adoption | ≤100 ms | 0.21 ms | 0.25 ms |
| Document size | ≤128 MiB | 6.97 MiB | 6.97 MiB |
| Real-media reopen verification | Separately timed | 85.4 s | 113.3 s |
| Active real-media cancellation | ≤1 s | 0.65 ms | 0.27 ms |

Full-media verification remains asynchronous, with progress and cancellation. It
is not included in the five-second document parsing budget. Worker preparation
took 378 ms / 770 ms respectively; neither cache deserialization nor evaluation
preparation blocks the UI thread.

Freeze one universal Release candidate (JUCE 7.0.12, macOS 10.13 minimum, existing
floating-point policy). Stabilize on arm64 Sequoia, then run that package on actual
Intel High Sierra. Build only explicitly named standalone/test targets. No VST3 or
installed application changes are included.

New tests cover exact stereo/5.1 and supported System I 7.1 round trips, settings
and cache compatibility, missing/relinked/changed media, malformed/oversized core,
cache damage recovery, atomic I/O failures, cancellation, concurrent save edits
and late results. A synthetic full-resolution three-hour Source/Current record
with 1,000 sections measures persistence costs; a real three-hour eight-channel
FLAC separately checks verification, progress and active cancellation. The former
is not represented as a new three-hour audio measurement.

Forty-eight reopened Current State capture cases cover 44.1/48/96/192 kHz,
stereo/5.1/7.1 and zero, negative, positive and +24 dB processing gain. Every channel
is captured before monitor controls, with injected devices only, and compared
bit-for-bit with the qualified offline path on each platform. Cross-Mac capture
limits remain bit-exact at zero gain and at most one float32 ULP otherwise; metric
budgets remain 0.01 LU/dB and 0.1 LU for LRA. These are implementation-test budgets,
not changes to profile acceptance rules.

Stopped reopened H.264 CFR/VFR and ProRes Proxy/LT/422 cases retain exact signed
offsets, seek/frame selection and unchanged audio identities; native reader audit
requires zero embedded-audio outputs. M1–M7 regression/official suites and optional
diagnostic-isolation checks ran freshly. Bidirectional cross-Mac reopen/relink and
silent interactive UI checks passed as separate gates.

The M7 occasional-picture-skip development acceptance is retained. No playback
scheduling fix is claimed. Sustained freezes, drift, audio interruptions, wrong
identity, memory and lifecycle failures remain unacceptable. Full physical 7.1
listening remains separately pending and is not needed for session persistence.

The numerical engine, processing/audition worker, verified-source implementation,
native video engine and EBU profile definition are byte-identical to M7's accepted
candidate-03. The source-equivalence review explicitly scopes reuse of its actual
three-hour measurement and physical timing/endurance evidence. M8 reruns the
official references, native M1–M7 suites and affected session/UI lifecycle paths;
it adds fresh full-media verification and restored-audition/video tests. Historical
M7 timing records retain their original candidate identity and picture-skip counts;
they are not represented as fresh M8 physical-device runs.

### Reproduction

Use a fresh build directory and the pinned JUCE 7.0.12 checkout. Configure Release,
`CMAKE_OSX_ARCHITECTURES=arm64;x86_64`, `CMAKE_OSX_DEPLOYMENT_TARGET=10.13`, and
enable the `LEVELSCOPE_BUILD_*_TESTS` options for MEASUREMENT, SCAFFOLD, COMPLIANCE,
VIDEO, PROPOSAL, DIAGNOSTICS, SECTIONS, AUDITION, MULTICHANNEL and SESSION, plus
`LEVELSCOPE_BUILD_VIDEO_PAINT_DIAGNOSTIC=ON`. Preserve the existing floating-point
compiler policy. Build named targets only:

```sh
cmake --build build/logic-m8 --config Release --target \
  LevelScopeStandalone LevelScopeMeasurementTests LevelScopeScaffoldTests \
  LevelScopeComplianceTests LevelScopeVideoTests LevelScopeProposalTests \
  LevelScopeDiagnosticsTests LevelScopeSectionsTests LevelScopeAuditionTests \
  LevelScopeMultichannelTests LevelScopeSessionTests LevelScopeSessionPlaybackTests \
  LevelScopeVideoQualification LevelScopeVideoPaintDiagnostic LevelScopeVideoPaintTraceTests
```

`Tests/package_m8_qualification.py --help` describes the pinned reference, UI and
video fixture inputs. It requires a new output directory and produces the source
and artifact manifest. Run inside that preserved package, with a new report folder:

```sh
python3 -B Tests/run_m8_qualification.py --platform sequoia --output /absolute/new/sequoia
# On actual Golem, use its existing /usr/local/bin/python3 (3.8.10):
/usr/local/bin/python3 -B Tests/run_m8_qualification.py --platform high-sierra --output /absolute/new/high-sierra
```

Copy each platform's `exchange-out` folder to a new folder on the other Mac, then
run its native `bin/LevelScopeSessionTests --exchange-in /absolute/incoming-folder
--output /absolute/new/exchange-report.json`. Foreign caches must be rejected while
the exact settings/annotations and content-verified relink are retained. Compare
the two complete folders with `Tests/compare_m8_reports.py`; retain raw `.f32`
captures locally for that comparison, not in Git. The M1–M5 comparison composes
through `compare_m5_reports.py`; M6 capture and M7 capture/native-video comparisons
remain separate. Interactive UI records are an additional requirement, not implied
by any automated summary.

### Checkpoint publication

All approved M8 gates are complete for candidate-01. Publication verified all 174
source inputs against the frozen manifest, all 379 package artifacts and all 116
indexed qualification evidence files. Latest remote main remained the M7 base;
the logic checkpoint was integrated by fast-forward without implementation changes.
The original universal package and binaries were retained without rebuilding or
relabeling their original measurement-build identity. The M7 ancestry commit is
not the M8 tested-source revision. Resolve the exact saved checkpoint/main identity
with `git rev-parse 'm8-reproducible-sessions-2026-10-25^{commit}'`; its full SHA is
also recorded in the immutable annotated tag.

Qualification reports retain their original pre-publication fields. The later
publication record and annotation establish the saved state. Documentation and
publication evidence are outside the compiled source manifest. If a later change
alters that implementation, qualify the affected result before integration.

The publication audit found all eleven pre-existing tag objects and target
commits unchanged, the registered video worktree unchanged, and all seven unrelated Finder metadata files unchanged.
Earlier packages and installations were not modified. No VST3 target was built
or installed, and no driver, Pro40 or system-audio configuration was changed.

## Remaining product boundaries and next work

Sessions are opened through the app’s explicit Open control; OS document
association/double-click launching is not included. Saving starts while stopped.

No separate-file import, stems, analysis references, speech detection, dialogue
metrics, section processing, dynamics, public export/report, autosave or project
media copying is implemented in M8. 7.1.4/Atmos and VST3 parity remain separate.
Caches are replaceable performance data; sessions are not archives of their media.
An incompatible/new build may require remeasurement even if its numerical engine
is unchanged. The initial schema has no forward-core-version migration guarantee;
future readers must explicitly migrate or reject it safely.

Next: a bounded two-mono stereo / six-mono 5.1 assembly proposal, with explicit
roles, matching sample rate/frame count, confirmed alignment, filenames alongside
roles, per-file verification, grouped persistence/relinking and an incomplete
group when any file is missing. No silent resampling, padding or channel reassignment.
Then verified PCM export/artifact-bound reporting. Stem composition and dialogue-
assisted analysis follow as separate bounded proposals. [35](35_PostM7BroadcastWorkflowProposal.md)
defines their distinct audible and analysis-only identity/invalidation boundaries.

# M11 implementation, qualification and handover

Status on 26 October 2026: **candidate-02 passed the approved M11 qualification
scope and is published/integrated into main**. All required offline, injected-device,
actual-platform long-workflow, exchange and silent UI gates passed. No required
M11 check remains outstanding. Physical full-multichannel listening is separately
pending and is not claimed.

The annotated checkpoint is
[`m11-aligned-stem-composition-2026-10-26`](https://github.com/carlozaiana/LevelScope/tree/m11-aligned-stem-composition-2026-10-26).
Its peeled commit is the accepted M11 checkpoint and main revision at publication.
The annotation records the full commit SHA and original source/build/binary and
evidence hashes. The [publication record](../../Tests/Reference/results/m11-publication-2026-10-26/README.md)
records the mapping and preservation checks. Latest main remained the M10 base;
integration was a fast-forward, with no change to the tested implementation.

The authoritative [acceptance receipt](../../Tests/Reference/results/m11-qualification-2026-10-26/acceptance-candidate-02.json)
and [evidence index](../../Tests/Reference/results/m11-qualification-2026-10-26/README.md)
identify the exact tested candidate and results. The user approved publication
and integration after readiness review in the coordinating chat on 26 October.
Original qualification receipts, including their pre-publication status fields,
remain byte-for-byte intact. This handover, publication record and annotated tag
describe the subsequent saved state. No rebuild or new physical-device test was
needed solely for publication, and no next milestone is started here.

Base: published M10 `fcd4c15fe78d6352f5e96faabf3eb52f942235e1`, branch
`codex/logic-module`, existing `/Users/fuzz/Developer/LevelScope` working copy.
The approved scope is [41](41_M11AlignedStemComposition.md); the plan recorded
before implementation is [42](42_M11ImpactQualificationPlan.md).

## Implemented architecture

`ProgrammeComposition` owns two separately accessible component definitions, M&E
and dialogue. Each component retains its original interleaved file or complete M9
mono-file assembly, member hashes, ordered roles, sample-zero alignment, purpose
and optional language/version labels. Exactly three named presets are accepted:
stereo+stereo, 5.1+mono-to-Centre, and 5.1+5.1. Output order is canonical
L/R[/C/LFE/Ls/Rs]; original interleaved component order remains explicit in its
role mapping. Labels and verified locator changes do not change audible identity.

`VerifiedSource` recursively verifies and owns every required member copy. The
composition reader retains two bounded readers, not a rendered programme-sized
intermediate. `ProgrammeReader::compose` uses fixed-order double addition rounded
once to float32; mono dialogue leaves non-Centre M&E samples unchanged. All
consumers then use the existing `static-gain-v1` operation. Measurement, audition
and export share this path. The real-time callback and video decoder/scheduler are
unchanged. Monitor controls remain separate from processing identity.

Finite over-full-scale sums remain measurable. Non-finite operands or overflow
fail explicitly. Monitor overload still silences the entire offending callback,
retains presented position and requires explicit restart. Export retains both
measured Current and actual reopened-artifact gates, with dither applied once at
final PCM24 quantization. No component measurement or prediction supplies a
programme verdict. The EBU profile definition/digest is unchanged.

The import/edit window is transactional and stopped-only. It requires explicit
confirmation that the exact mapped inputs are complementary, correctly balanced
stems intended to be combined at sample zero. A new composition initializes
Current only on request. Editing its audio retains an existing global gain but
invalidates affected results and clears programme-bound section annotations.
Label-only edits retain measurements. Failure/cancellation retains the previous
valid session, stopped. Paths and filenames remain accessible in Source details
and the component editor; no source media are rewritten.

## Session and report compatibility

Session schema 3 adds the versioned composition recipe. Readers retain schema 1
and 2 support; saving upgrades only on an explicit Save/Save As. The UI explains
that M11 is required and offers Save As to preserve the original older-readable
document. Independent caches remain replaceable; damaged/incompatible caches
cannot prevent recovery of a valid core. Reopening verifies all component members,
reconstructs evaluations from compatible measurements, and remains stopped,
device-closed and routes-muted. Missing members never become silent channels.

Reports use `levelscope-measurement-report-v2`, identify the original composed
Source and measured Current, and record component file hashes, mapping, unity
gains, summing version and confirmation. Programme definition digests are labelled
separately from encoded file hashes. Actual exported audio has its own measured
domain and zero verification gain. Measurement-only reports invent no artifact.

## Future processing extension points

The architectural sequence is original components → optional versioned Current
component processing → versioned composition → versioned whole-programme
processing. M11 executes identity component processing and the existing explicit
global trim only. Schema 3 reserves component-plan and whole-programme-graph slots
and rejects unsupported nonempty plans. Introducing DSP requires a new processing
version and supported plan representation; existing Source remains the original
confirmed composition. Never reinterpret old plans under new processing semantics.

Future Current component settings must bind their component content/roles/range,
processor version and full-precision settings, then feed Current programme
identity. Changing them invalidates affected measurements, diagnostics, proposals,
prepared audition and export approvals. A separate non-audible analysis-reference
identity must invalidate only dependent analysis when programme audio is unchanged.
M11 rejects unknown purposes instead of making them audible.

Dialogue levelling/DSP relative to M&E remains a future proposal after the broader
whole-programme proposal/apply/verify workflow. Future targets must name the
analysis method and acceptance criteria. Explicit application must be followed by
measurement of the intended dialogue improvement and of the complete resulting
programme. Dialogue-stem loudness, speech-gated mix measurements, dialogue/background
diagnostics and an intelligibility score are distinct methods, not interchangeable
claims. No component DSP, speech detection or intelligibility analysis ships here.

## Qualification inventory and completed evidence

`package_m11_qualification.py` packages explicitly named standalone/test targets
and selected pinned assets. It never calls older packagers. The M11 runner has
literal required inventories for short, consolidated-long and cross-Mac exchange
phases; missing/duplicate/failed runs, wrong platform, changed binaries, missing
rate/preset cases or incomplete UI receipts prevent a complete qualification.
The validator itself has failure-injection tests.

Development evidence is under `build/m11-development-2026-10-26/`. It is not a
substitute for frozen-candidate qualification. It includes independent composition
sample comparisons, numerical/identity/failure checks, measured PCM export and
native H.264/ProRes selection with injected devices. Historical failed development
runs remain preserved, including corrected test assumptions about the monitor
startup ramp and JUCE's unsupported AIFF PCM32 writer (an independent AIFF fixture
writer now tests the supported reader).

Candidate-01 was frozen with source manifest
`4c13d106b45dfbc115a1954ff156184fadb5abc8d40b334f36390aef69cd6176`.
Its native UI exposed an incorrectly encoded separator in the composed audition
label. It is retained unchanged and is not the accepted checkpoint candidate.
The replacement uses an ASCII separator; its test receipts also enumerate every
native video selection and all six cross-platform session cases explicitly, so
missing or duplicated cases cannot satisfy qualification. No long workflow was
run on candidate-01.

The frozen identity and completed evidence are recorded below. The complete
aggregate gate passed using both actual-platform records and all required phases.

Qualified frozen candidate: **candidate-02**, published at the checkpoint identified above.

- Source manifest: `669b350652dbc3fb09d0ff724758208fec21bba40e0f6bc0188b89e6404c2c61`.
- Measurement build: `70dc29e6fabde64b0be2a0f12790528de7e088956d5817f82b6888611c723d98`.
- Universal app executable SHA-256: `deb189e943529dd3e9e9d63a93e1d1cd2003f36752f7ab43f1e099d6c204eae5`.
- Local package: `build/m11-qualification-2026-10-26/candidate-02/`.
- Golem package: `/Users/fuzz/LevelScope-M11-qualification-2026-10-26/candidate-02/`.
- The manifest covers 222 implementation, build and test-source files. Final
  qualification receipts and documentation are recorded afterward and are outside
  that manifest. The checkpoint retains those exact source hashes alongside the
  completed evidence; its Git hash differs from the recorded M10 build ancestry.
  No rebuild may replace the preserved tested universal app under this identity.
- Both actual-platform short inventories passed: 20 explicitly selected runs each,
  9.76 minutes on arm64 Sequoia 15.7.9 (24G830), Apple M1 Max, and 18.95 minutes on
  x86_64 High Sierra 10.13.6 (17G66), Intel Core i7-3770K.
- Each platform passed 24 composition cases (three presets, four rates, interleaved
  and grouped representations), 72 pre-monitor captures and 24 independently
  decoded PCM24 exports. Across Macs all 72 captures and all 24 finalized short artifacts are
  byte-identical, and the compared metric differences are zero. These comparison tolerances never enlarge
  profile acceptance allowances.
- Six cross-platform session cases passed in each direction: same audio/settings/
  confirmation survives; foreign-platform caches retain provenance but require
  local remeasurement. Eleven observed silent UI cases passed on each Mac.
- Sequoia's consolidated three-hour workflow passed in 56.54 minutes, with 1,000
  sections, 518,400,000 frames and a 9,331,200,104-byte verified RF64 export.
  Additional RSS peaked at 127.5 MiB; UI event-loop p95 was 12.57 ms; cancellation
  22.38 ms; save 91.21 ms; core/cache parse 450.99 ms; session 10,444,186 bytes;
  JSON+HTML 385,600 bytes. Observed filesystem growth peaked at 47,377,219,584
  bytes (includes concurrent small qualification/UI activity, not isolated I/O).
  Independent libsndfile 1.2.2 decoded the real RF64 and matched 98,304 PCM codes
  against the algebraic composition plus single-final-dither oracle at selected
  spans; the entire artifact hash and report/HTML/manifest hashes agree.
- High Sierra's same consolidated workflow passed in 105.64 minutes. Additional
  RSS peaked at 93.12 MiB; UI event-loop p95 was 12.25 ms; cancellation 80.08 ms;
  save 233.91 ms; core/cache parse 1,037.24 ms; session 10,443,758 bytes;
  JSON+HTML 385,672 bytes. Observed filesystem growth peaked at 46,721,753,088 bytes.
  Its independent decoder matched the same 98,304 oracle PCM codes. The final
  9,331,200,104-byte RF64 has the same SHA-256 on both Macs:
  `cdc392649f4e7e7b0cd68e99998eea71ae8d402abefcbc1ea311a0a205b4d812`.
- The complete aggregate comparison passed: zero differing float32 ULPs in all
  72 captures and zero differences in compared short and long Source/Current/
  artifact metrics. Both three-hour runs stayed inside every approved resource
  budget. Each retained long-media evidence set occupies approximately 28.0 GB;
  temporary verified copies were released by their owned lifecycle.
- Representative video integration passed six Source/Current H.264 CFR/VFR and
  ProRes cases per Mac, with 60 native frame selections each. Signed rational
  offsets, seeks/restart/EOF, uninterrupted injected audio positions, identity and
  undecoded embedded audio were checked. Recorded player skip counts were zero in
  these bounded selections; this is not a new physical-paint/endurance claim.

The selection includes fresh pinned EBU (70 files: 66 minimum-requirement and
four supplementary alignment files) and ITU (39 vector cases) checks. Individual
vector hashes are pinned and verified; no archive-verification claim is made.
The selection covers changed shared input/session/export callers; focused legacy
interleaved and separate-file compatibility; M11
four-rate/three-preset/grouped/interleaved captures; all-channel peak constraints;
and representative H.264 CFR/VFR plus ProRes Source/Current checks. Earlier real
long-media workflows are not recursively launched. The session regression uses
its short mode because real persistence/resource coverage belongs to the single
consolidated M11 long workflow.

Reused evidence remains identified as M10 candidate-03, source digest
`4ff89967e4927600d0b6defafabb28c372f4068870a0f1a9adb1ee84fe1781d8`, app SHA-256
`92c013999bdecaf904af0888ce68dd59bcd1b3004934b9ddb36ce1a1bf511e15`.
Its original full codec/colour and physical endurance records remain historical
evidence for unchanged video/callback behavior, not new M11 observations.
The final dependency audit found 55 relevant source files byte-identical to M10,
including the measurement/profile kernels, real-time callback, video timing/decoder
and PCM24 quantizer. Fresh changed-caller tests passed. The exact boundaries are in
[m10-evidence-equivalence.json](../../Tests/Reference/results/m11-qualification-2026-10-26/m10-evidence-equivalence.json).
No timing/buffering/decoder change or integration problem triggered an additional
physical block; no physical device was opened for M11. Occasional picture skips remain an accepted
development limitation; prior strict failures are not rewritten. Audio interruption,
accumulating drift, sustained freezes, identity and lifecycle/resource failures
remain blockers. Physical full-multichannel listening remains separately pending.

## Reproduction

The frozen package is the tested artifact; a later rebuild is a new candidate.
It contains `identity.json`, the exact source manifest, 13 named test executables,
the standalone app, pinned reference vectors, representative video fixtures and
the selected Python runners. The build uses JUCE 7.0.12 at
`4f43011b96eb0636104cb3e433894cda98243626`, Release, `arm64;x86_64`, deployment target
10.13. The exact compiler and build configuration are in its identity file.

For development, configure a separate build directory with
`FETCHCONTENT_SOURCE_DIR_JUCE` pointing to that checkout. Enable the explicit
`LEVELSCOPE_BUILD_{COMPOSITION,MEASUREMENT,SCAFFOLD,COMPLIANCE,PROPOSAL,DIAGNOSTICS,SECTIONS,AUDITION,MULTICHANNEL,SESSION,PROGRAMME,EXPORT,VIDEO}_TESTS`
options. This notation lists individual CMake options; it is not a literal CMake
argument. Build only `LevelScopeStandalone` and the 13 `LevelScope*Tests` targets
listed by `Tests/package_m11_qualification.py`. Do not build the default all target
or VST3, and do not install anything.

From the frozen package directory, run its `Tests/run_m11_qualification.py` with
`--platform sequoia` or `--platform high-sierra`, `--phase short`, `--output` a new
absolute evidence directory and `--decoder-library` the inspected native
libsndfile 1.2.2 library. Repeat with `--phase long` only once per candidate/platform
unless a relevant failure or implementation change requires it. The runner checks
the actual OS/CPU architecture and frozen hashes before execution. Python 3.9.6 was
used on Sequoia and 3.8.10 on High Sierra; no new runtime was installed.

The exchange phase additionally takes `--exchange` pointing to the other Mac's
six generated composition case folders. It rejects foreign-platform cache reuse,
verifies members and measures locally. `Tests/compare_m11_reports.py` takes the two
platform result roots, `--output` and `--complete`; the latter also requires both
long/exchange phases and all 11 explicit UI receipts. Raw individual phase summaries
intentionally never claim complete qualification. Sample captures remain in the
original evidence folders; the repository stores their hashes and text receipts.

## Limitations and publication boundary

No arbitrary stem count, component gains/mute/solo/DSP, nested composition,
automatic alignment, resampling, padding, truncation, pan law, stem analysis,
dialogue intelligibility, 7.1 composition/export, 7.1.4, Atmos, section processing,
or VST3 parity. Imported M&E/dialogue content and balance are user assertions,
not automatically certified. Only selected whole-programme loudness/true-peak rules
are evaluated. Existing M10 container/export/report limits continue to apply.

Implementation and the approved qualification are complete; the user subsequently
authorized publication and main integration, completed as recorded above.
Further implementation requires a new approved milestone. Preserve tested packages, all existing tags,
installations, the video working copy and Golem's relocated ProRes evidence.

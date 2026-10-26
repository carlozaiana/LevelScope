# M10 - verified PCM export and reporting: qualification and handover

Status: **candidate-03 qualified on both actual platforms and published/integrated into main, 26 October 2026.**

The annotated checkpoint is
[`m10-verified-pcm-export-2026-10-26`](https://github.com/carlozaiana/LevelScope/tree/m10-verified-pcm-export-2026-10-26).
Its peeled commit is the accepted M10 checkpoint and main revision at publication.
The annotation records the full commit SHA, original source/binary hashes and
qualification evidence identity. The
[publication record](../../Tests/Reference/results/m10-publication-2026-10-26/README.md)
documents that mapping and preservation checks.

M9 ancestry: `75843dcd017ec7990be2200454a10fc38beb44c4`, on the established
`codex/logic-module` working copy. The M9 package, all checkpoint tags and the
separate video worktree are preserved. The user authorized publication and merge
on 26 October 2026. Every required M10 gate below is complete. The aggregate
[acceptance receipt](../../Tests/Reference/results/m10-qualification-2026-10-25/acceptance-candidate-03.json)
binds the underlying reports and records `qualified: true`.

## Implemented boundary

The standalone exporter supports complete stereo and conventional 5.1 programmes,
including confirmed M9 mono groups. The shared verified programme reader and exact
all-channel static-gain processor feed a separate final PCM24 quantizer. Monitor
controls, reference-video audio and section annotations never enter the render.

The output uses the original sample rate and frame count. Canonical output order
is L/R or L/R/C/LFE/Ls/Rs; an explicit permutation preserves the confirmed input
interpretation. Ordinary extensible RIFF/WAVE is used when its complete chunk size
fits; otherwise the narrow PCM24 RF64 subset uses an exact `ds64` length/frame
count. Masks are 0x3 and 0x3f. The latter is the established conventional 5.1
back/surround interpretation, not System I 7.1 rear weighting.

`pcm24-tpdf-splitmix64-counter-v1` freezes the final quantizer. Each output channel
has a counter-addressed stream; two independent uniform 53-bit fractions are
subtracted to produce TPDF in (-1,1) PCM24 LSB. The exact implementation constants
and counter arithmetic are in `Pcm24Wave.cpp`; the independent Python oracle is
`Tests/verify_m10_pcm.py`. Ties round to even. Reject non-finite values, input
outside [-1,1), and rounded codes outside the PCM24 range. Never clip, attenuate,
retry or silently correct. The seed binds programme/processing/format identity,
not monitor settings, pathnames or profile choice. Dither happens **once**, only
at this final boundary. Zero processing gain does not mean a bit-identical copy.

Rendering holds an immutable session request and blocks conflicting edits while
running. The finalized closed output is reopened and fully measured in a distinct
artifact domain with zero additional processing gain. The original Source and
Current State measurements remain unchanged. The EBU profile definition/digest and
strict full-precision allowances remain unchanged; evaluation semantics are
versioned to `measured-programme-rules-m10.1` to admit a separately bound artifact.

**Render and verify** creates an owned staging artifact. **Publish export** is a
separate explicit action requiring both measured Current and actual artifact to
pass the selected rules. Source need not pass. Publication rechecks original
members, exact artifact bytes, coherent report approval and all report/manifest
hashes. Exclusive same-filesystem directory rename exposes a complete new bundle;
existing destinations are never overwritten. This is a local namespace-atomic
transaction, not a new power-loss recovery guarantee. Private incomplete staging
folders are never auto-published. Abrupt-process tests distinguish pre-rename
staging from a fully committed bundle and preserve unowned files. Cancellation and
ordinary failures clean owned partials. Abrupt termination can leave a hidden
`.levelscope-export-*` directory/private source copy for explicit inspection and
cleanup; M10 does not automatically scavenge or resume crash remnants.

The export UI uses the session's immutable evaluation cache and a worker-validated
artifact summary; it does not serialize the complete three-hour timeline at every
refresh. Heavy full-snapshot identity validation remains in the worker. Render and
publication jobs own their requests and block conflicting edits. Publication's
exclusive rename is its commit point; a Cancel arriving afterward cannot relabel
an already complete bundle as unpublished.

Source/Current measurement-only reports remain available when rendering is blocked,
including failing or unavailable Current measurements. They contain no artifact
fields or invented export hash/verdict. An actual measured staging artifact can
also have its own report saved without audio publication, including a failed
verification, explicitly labelled NOT EXPORTED. Reports are derived JSON/static
HTML; saved report verdicts are never imported as measurement authority.

Implementation entry points for subsequent work:

- `StandaloneExport.h/.cpp`: `capture`, `specification`, `readiness`, `render`,
  `artifactEvaluation`, `approvalIdentity` and `publish`. `Request` captures one
  session revision; `Candidate` owns the staged artifact and verified source copy.
  Worker services take explicit cancellation and return a bound result. The caller
  locks conflicting edits and rejects completions from an obsolete runtime/job.
- `Pcm24Wave.h/.cpp`: checked RIFF/RF64 sizing, explicit PCM24 writer finalization,
  bounded parser, versioned counter dither and quantization. Changing these
  semantics requires renewed sample/container qualification and a new identity.
- `StandaloneReport.cpp`: one derived model feeds `reportJson` and `reportHtml`.
  `measurementReport` and `artifactReport` remain distinct from audio publication.
- `StandaloneExportUI.cpp`: format/mapping/size preview, asynchronous progress and
  cancellation, live profile/applicability evaluation and explicit publication.
  Monitor controls remain outside the processing and approval identities.
- `StandaloneMeasurementService` and `VerificationEngine`: explicit rendered-output
  measurement domain and artifact hash binding; Source/Current rules remain intact.

The M10 service version is `verified-pcm-export-m10.1`; session schema remains 2.
No plugin processing, parameters or saved-state migration is part of these changes.

## Preserved development history

All development logs are preserved under
`build/m10-qualification-2026-10-25/`. They include early failing iterations and
must not be substituted for the frozen candidate reports below.

- Independent PCM24 oracle: initial development matrix passed 196,800 exact codes
  across 48 layout/rate/gain combinations, including deterministic chunking and
  the zero-gain sample-change distinction.
- Initial UI stereo render/reopen/explicit publication worked while device-closed,
  stopped and routes-muted. Independent libsndfile 1.2.2 decoded all 768,000 sample
  codes exactly; JSON/HTML/manifest agreed. This is development evidence only.
- Expanded development matrix: 4,688 assertions, 16 grouped/interleaved
  rate/layout artifact cases and ten native/injected video cases passed. All 16
  bundles also passed independent libsndfile decoding and report/hash comparison.
  A subsequent development run passed 4,761 assertions, including abrupt-process
  termination, shared-volume disk budgeting, verified-copy reuse, actual artifact
  failure and profile/report invalidation. These are development runs; only the
  frozen-candidate record below establishes qualification status.
- A short-programme negative test initially expected FAIL below the absolute gate;
  the correct result is INDETERMINATE. The test now uses an available out-of-target
  loudness value. No profile allowance or engine status was changed for this test.
- The initial real three-hour 5.1 RF64 run passed: 9,331,200,104 bytes,
  55.8 minutes including Source/Current prerequisites, 80.0 MiB extra RSS,
  15.1 ms heartbeat p95, 11.6 ms cancellation and 381,425 report bytes. It is
  development evidence, not the final candidate's two-platform qualification.
- Candidate-01 was frozen and preserved with source digest
  `37d67470921412a0318e76fa38e12be9e5b941cf9791d8931a8247be01302718`.
  Profiling the initial long run identified unbuffered 64-byte hash reads and
  repeated filesystem size queries. Candidate-02 buffers these reads and
  validates the final length; the SHA algorithm and PCM samples are unchanged.
  Candidate-01 is not labelled qualified. Its package and reports are retained.
- Candidate-02's Sequoia core suite passed, including M1–M9 regressions, 16 actual
  exports and independent decoding. Its interactive applicability check found a
  stale duplicate PASS in the render-completion notice: Publish correctly became
  unavailable, but the notice retained the prior evaluation. This is a **failed
  UI check**, not a qualified candidate. Candidate-03 removes the cached verdict
  and derives the displayed artifact evaluation from the current profile and
  applicability. Candidate-02's package and long-run evidence remain intact.
- The original timing fixture was rejected correctly: its normalized Current
  true peak was +0.59 dBTP. A separately generated, explicitly aligned continuous
  stereo fixture replaces it for timing preparation only. Its generator, exact
  PCM periods and file hashes are recorded under `physical-clock-fixtures`.
  Export limits were not relaxed. This change is auxiliary test media, not product
  processing or automatic correction.
- Both candidate-02 core runs passed on the actual Macs. Their sixteen exported
  artifacts have identical SHA-256s; fixed-seed PCM24 codes and pre-quantization
  captures match exactly. Maximum measured difference is
  `7.105427357601002e-15` dB. The original comparison harness erroneously compared
  entire per-channel metric objects bit-for-bit, including double/hex values.
  Candidate-03's comparator uses the existing 0.01 dB metric tolerance while
  retaining exact roles, availability, units and identical-audio peak positions.
  This fixes the harness, not the profile's strict acceptance rules. Regression
  tests reject changed positions, changed roles, out-of-budget/non-finite metrics.
- The independent decoder's sampled seek plan now bounds its pseudo-random seek
  for medium-length files; the former fixed index exceeded a 30-minute fixture.
  Short-file all-sample checks and the real three-hour first/random/middle/final
  sample plan are unchanged. This is a qualification-tool correction, with an
  explicit bounds regression test; it does not alter the product reader.
- The frozen-candidate gates below supersede the preliminary development checks:
  full M10 matrix, actual three-hour RF64 on both Macs, independent decoding,
  source/binary identities, M1–M9 regressions, comparisons/exchange, silent UI and
  the coordinated video block. Earlier development successes do not replace them.
- Sequoia development found and fixed uninitialized section-map approval identity,
  staged-report tampering checks and explanatory-text encoding. Harness setup
  errors (all-muted permission and staging-parent selection) are retained in logs.

The approved physical block is Mac Studio Speakers at 48 kHz and Golem Built-in
Output at 44.1 kHz, all routes muted: four short Source/Current H.264 VFR/ProRes 422
cases and one 30-minute Current/ProRes case per Mac. Announce before opening devices.
Stereo physical checks do not establish full multichannel listening qualification.

## Qualified candidate and results

Current frozen candidate: **candidate-03**, qualified on actual arm64 Sequoia and
Intel High Sierra. No implementation changes followed its freeze.

- Source manifest SHA-256: `4ff89967e4927600d0b6defafabb28c372f4068870a0f1a9adb1ee84fe1781d8`.
- Measurement/build identity: `f36567f5f38a47daef49bc9543252b42b3403e3bdcc2fe986b755bee1576cce9`.
- Package `identity.json` SHA-256: `0b0cd42a9f4ee917d6d62f25a903edd7ea354a7d1313d6dfa64341396453cbdd`.
- Universal app executable SHA-256: `92c013999bdecaf904af0888ce68dd59bcd1b3004934b9ddb36ce1a1bf511e15`.
- Universal export-test executable SHA-256: `496a8952947068dbc5a9cad5c9cab7a44daa412fd610ddf557011737b4f6eaab`.

Actual qualification platforms:

- Sequoia 15.7.9 (24G830), Apple M1 Max, native arm64.
- High Sierra 10.13.6 (17G66), Intel i7-3770K at 3.50 GHz, native x86_64.

Both used the same frozen universal package and pinned JUCE 7.0.12 dependency.
No Intel-under-Sequoia run substitutes for the High Sierra evidence.

The candidate-03 core suites passed on both actual platforms: each completed
4,985 M10 assertions, sixteen grouped/interleaved exports, the independent PCM24
oracle/decoder checks, ten native/injected video cases, and the relevant M1–M9
regressions with pinned EBU/ITU vectors. The sixteen exported artifact hashes match
between Macs; pre-quantization captures and fixed-seed PCM24 codes match exactly.
Maximum measured export difference is `7.105427357601002e-15` dB, within the
documented comparison budgets without widening profile acceptance.

Both silent UI records are complete. The actual UI binaries match the frozen
app hash, including the corrected live applicability verdict, monitor separation,
explicit publication, stale/failing Current reporting, cancellation/close guard
and existing-destination rejection. These checks stayed stopped, device-closed and
routes-muted. Newly saved test sessions preserve their settings and results.
The user confirmed that the candidate-02 HTML reports were readable and
correctly distinguished measurement-only reports from verified exported audio;
the report implementation is byte-for-byte unchanged in candidate-03. Candidate-03
report/HTML/manifest checks ran on both Macs. Bidirectional exchange is complete:
each recipient freshly imported/measured sixteen opposite-Mac bundles and checked
them with its independent libsndfile 1.2.2 decoder. Saved report verdicts were not
adopted as session measurements.

Both native candidate-03 real three-hour RF64 runs passed. Each produced
9,331,200,104 bytes: 518,400,000 frames, six channels, 48 kHz. Six all-channel
first/middle/final Source/Current injected captures matched exactly on each Mac.
Source, Current and artifact metrics, per-channel peak values and locations agree
exactly between platforms. Both audio artifacts have SHA-256
`b1d75b637bd8f6893109c61674ed8513a12721978eb850f1d1c5c855f3e5c2ff`.

| Three-hour test measurement | arm64 Sequoia | Intel High Sierra | Approved budget |
| --- | ---: | ---: | ---: |
| Source/Current preparation, render/verify, publish and extra RF64 audition | 44.96 min | 83.53 min | 180 min |
| Additional export RSS above existing result storage | 54.56 MiB | 1.88 MiB | 256 MiB |
| UI heartbeat p95 | 15.05 ms | 15.03 ms | 100 ms |
| UI inspection p95 | 0.384 ms | 0.882 ms | 100 ms |
| Cancellation acknowledgement | 10.06 ms | 30.35 ms | 1,000 ms |
| Combined JSON/HTML with 1,000 annotations | 381,425 B | 381,497 B | 16 MiB |

Both independent decoders verified 98,304 sampled PCM codes against the known
source/quantizer. Stage timings, exact bytes and cleanup latency are retained
in each raw `long-export.json`; the RSS row is growth, not total process memory.

Both physical blocks passed: four short Source/Current H.264 VFR/ProRes 422 cases
plus one 30-minute Current State ProRes run per Mac. All ten cases completed once;
none was repeated to obtain a zero-skip result. There were zero picture skips,
audio underruns, bounds/deadline failures or accumulating internal-clock drift.
The endurance maximum paint intervals were 33.664 ms on Sequoia and 20.382 ms on
High Sierra; both measured memory-growth checks passed. Audio device configuration
matched before and after every block. Every route remained muted. The raw reports
retain timing, memory samples and strict verdicts, and the separate
`physical-video-summary-candidate-03.json` rechecks the accepted development policy.
These are internal audio-clock/UI and device-lifecycle checks, not optical
synchronization or full physical multichannel listening qualification.

Raw sub-suite summaries deliberately retain `qualification_complete: false`:
their runners do not claim to cover the separate UI, cross-platform and physical
gates. The final M10 aggregate receipt combines those independent records;
do not rewrite the raw runner outputs to imply they independently qualify M10.

Golem's interrupted candidate-02 long run has no final result or publication
receipt. Its last log is 25 October 23:01:34 CET; the host booted again on 26 October
10:18:07 CET. It is recorded as incomplete, not a pass. New candidate-03 long runs
on both Macs use fresh result folders and retain all older packages and evidence.
Golem's candidate-03 binaries/source are independent copies; unchanged pinned media
directories refer to candidate-02 media without modifying those files. Every
manifest content hash is still checked. This avoids duplicating large reference
media while preserving space for RF64 output and injected audition copies.

Preserved candidate-02 identities (superseded, **not qualified**):

- Source manifest SHA-256: `e0eb0497e2c5c8b1e629dc5747e260cf21c8ef77ef2b0e39094294820fc6f996`.
- Measurement/build identity: `3e24a19f92e3b54738be169bc222c0a433414cca577687c42a648453e305d6ce`.
- Package `identity.json` SHA-256: `225fa11cd0274d1114c9f66dcc62e8bb8e189db77639f931f8b7559722b3d038`.
- Universal app executable SHA-256: `ab9e4c4dedcf3c2c92992483216fb0a4af5206a04cca0bd2efd99bd3b63282e9`.
- Universal export-test executable SHA-256: `0b1f0ceeb94d3a6b8eec42a1e50180ba741518c78822f667d8c2d046dd8b2ce1`.

The final per-platform reports match candidate-03 identities. Candidate-01/02
and earlier development runs remain historical evidence. Transition manifests list
every changed source and explain any component-equivalent evidence reuse without
relabeling the binary that actually ran.
The accepted checkpoint is saved on `codex/logic-module` and integrated into main
by fast-forward from M9. No integration fix or tested implementation change was
required; publication reused the unchanged candidate's complete evidence.

The complete text evidence is indexed in
[`Tests/Reference/results/m10-qualification-2026-10-25`](../../Tests/Reference/results/m10-qualification-2026-10-25/README.md).
Its manifest hashes the collected reports; original packages, binaries, media and
sample captures remain in their qualification folders. Both package integrity
checks cover all 470 recorded files. The current implementation still matches all
208 frozen source inputs. All thirteen existing tag objects/targets, registered
worktrees and unrelated files match the preservation record. Existing installations,
plugins, the separate video working copy and system audio configuration are unchanged.

Preserved tested applications:

- Sequoia: `/Users/fuzz/Developer/LevelScope/build/m10-qualification-2026-10-25/candidate-03/LevelScope Standalone.app`.
- High Sierra: `/Users/fuzz/LevelScope-M10-qualification-2026-10-25/candidate-03/LevelScope Standalone.app`.

The package's base commit remains M9 ancestry, not a commit containing M10. The
annotated M10 checkpoint establishes the relationship to this source manifest and
these original tested binaries; documentation and evidence were completed after
the freeze. There are no outstanding checks within the approved M10 scope.
Original qualification records retain their pre-publication status fields; the
publication record and tag describe the subsequent saved state. Do not substitute
a later rebuild for the preserved tested package or label a changed integration
candidate qualified without affected checks. Further implementation requires a new
approved milestone.

Each platform requires the bounded offline/export matrix, all-channel capture
parity, independent decoder and report checks, real three-hour 5.1 RF64/resource
budgets, and final silent UI. Shared measurement and reader regressions use the
pinned EBU/ITU vectors. Native/injected video checks cover reimported exported PCM,
Source/Current, H.264 CFR/VFR and ProRes Proxy/LT/422. The physical block remains a
separate gate. Earlier evidence may be reused only with explicit component-level
equivalence and original identities; no blanket prior-milestone pass is inferred.

Preserve the accepted development limitation: occasional picture skips are recorded
with timing, but are not automatic blockers if audio is continuous and video follows
the authoritative audio clock without accumulating drift. Sustained freezes,
audio interruption, identity/lifecycle errors and memory failures remain blockers.
Historical strict zero-skip failures remain historical failures.

## Reproduction and review boundaries

Use the pinned JUCE 7.0.12 dependency, Release, C++17, `arm64;x86_64`, macOS
deployment target 10.13, and the recorded strict floating-point settings. Configure
`LEVELSCOPE_BUILD_EXPORT_TESTS=ON` alongside the existing measurement, scaffold,
compliance, proposal, diagnostics, sections, audition, multichannel, session,
programme and video test options. Build only `LevelScopeStandalone`,
`LevelScopeExportTests` and the explicitly named regression/test targets recorded
in the candidate identity. Do not use the default all-target build or install a plugin.

`Tests/package_m10_qualification.py` freezes the universal app/tests, full source
manifest and pinned M9/M7/M1 assets into a **new** candidate directory. The base Git
commit is ancestry only: the source digest identifies the working implementation.
Documentation can be completed later without replacing the frozen tested binaries.

From the frozen package directory, use its `Tests/run_m10_qualification.py` with
`--platform sequoia` or `--platform high-sierra`, a fresh `--output`, `--phase core`
or `--phase long`, and the inspected `--decoder-library` path. The actual OS and
executable architecture guards must pass. The independent decoder is a qualification
tool only; it is not bundled into or required by the product.

`Tests/compare_m10_reports.py` compares core reports and fixed-seed sample captures;
`--long` compares the actual long artifact measurements. `LevelScopeExportTests
--exchange <received-bundle-folder>` reimports and remeasures all sixteen opposite-Mac
bundles instead of trusting their saved verdicts. Independent decoder checks also
validate each exchanged report/HTML/manifest and exact artifact hash. Reimported
RF64 receives bounded first/middle/final Source/Current audition captures using an
injected device, including all six channels and verified-copy reuse.

Only the separately approved `run_m10_muted_video.py --allow-audio` block may open
the specified physical device. Other qualification tools use injected devices or
remain offline. Keep all physical routes muted. A completed automated summary does
not replace the separate interactive UI, physical video and cross-platform gates.

## Remaining product boundaries

No complete broadcaster/BWF/BW64/metadata certification. No resampling, alternate
bit depths, mono/5.0/7.1 export qualification, lossy/video export, implicit downmix,
stems/dialogue processing, dynamics, automatic correction, section processing,
batch/network publication or VST3 parity. Session schema remains M9 schema 2;
active exports and approvals are not persisted. M9 sessions retain settings,
member identities, confirmed alignment, annotations and rational video offsets.
Their measurement caches follow the existing conservative build/platform policy:
opening them in this changed M10 build requires fresh measurements before export;
saved verdicts cannot satisfy either publication gate. Physical full-7.1 listening remains
separately pending. Further immersive expansion remains deferred.

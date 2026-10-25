# M7 qualification and handover — System I channel-based PCM 7.1

Current status, 25 October 2026: **candidate-03 qualified for the approved offline,
automated audition/routing and UI/video scope on both actual Macs, published and
integrated into main**. The annotated checkpoint is
[`m7-pcm-7.1-2026-10-25`](https://github.com/carlozaiana/LevelScope/tree/m7-pcm-7.1-2026-10-25).
Its peeled commit is the accepted checkpoint and main revision at publication;
the annotation records the full commit SHA and tested source/binary identities.
Full physical eight-channel listening remains separately pending.
The user lifted the
freeze and accepted occasional picture skips under the explicit bounded decision
in [34](34_M7CompletionAcceptance.md). Main's optional diagnostic tooling is now
incorporated; it is not a playback fix. Earlier packages and results below retain
their original identities and acceptance decisions.

## Combined candidate-03 — current completion record

The logic branch incorporates published main
`5e07bbd29930db34236d6dfcede7fea14983f90a`, preserving the unmerged M7 work.
The optional diagnostic target and its tests are integrated. Tracing remains
absent from the normal app and qualification executable, verified on both native
architectures. This integration does not claim a playback fix.

Frozen package `build/m7-completion-2026-10-25/candidate-03` on Sequoia and
`/Users/fuzz/LevelScope-M7-completion-2026-10-25/candidate-03` on Golem:

| Identity | Value |
| --- | --- |
| Source manifest | `7abc7587514b7f0f750720ac7853abeaf55b21c90697b25de7740cec2df034bf` |
| Measurement build | `facfe596c7126188b24627dfeab06dec907332b96442634c9ad8ee29675ea069` |
| Package identity SHA-256 | `b2cf915f2f6b644993782363b0a0c243c49b43198100ae39ff6be870eff35f9e` |
| Standalone executable SHA-256 | `14ba35d2f6ddc0bded680df267e81996f768c1537be1d6d360da87376b6889db` |
| Normal video harness SHA-256 | `dbeb25405e6d2ee7fc8f54ad425675830656a29e8103529651df92432fa8f89d` |

The frozen manifest describes the working tree before publication. Every one of
its 165 source/test/manifest inputs matches the checkpoint's corresponding file;
its original ancestry and package identities remain unchanged. Publication adds
documentation and text evidence outside the compiled source manifest, without
rebuilding or relabeling the tested binaries. Resolve the exact accepted commit
with `git rev-parse 'm7-pcm-7.1-2026-10-25^{commit}'`.
The universal Release package uses pinned JUCE
7.0.12, the existing non-fast-math FP policy and macOS 10.13 deployment support.
Only explicitly named standalone/test targets were built; no VST3 build/install.

Fresh native execution passed on arm64 Sequoia 15.7.9 (24G830), Apple M1 Max,
and Intel High Sierra 10.13.6 (17G66), i7-3770K. Both ran M1–M7 regression suites,
pinned ITU (39) and EBU (70) vectors, diagnostic isolation/runner guards and native
H.264 checks. All 218 M7 and 124 M6 captured streams were bit-identical between
platforms. M1–M5 metric comparisons and diagnostic/section resource budgets passed.
Silent integration UI checks passed on both: explicit roles, measured Source and
Current, unverified gain changes, diagnostic channel/time labels, held reference
seek/offset, stopped preparation, source clearing and normal close.

[Completion evidence](../../Tests/Reference/results/m7-completion-2026-10-25/README.md)
contains exact reports and explicit reuse reviews. Production measurement,
processing, audition, profile, diagnostics and section inputs match candidate-02;
the sole production delta is a read-only diagnostic friend declaration. Its real
three-hour measurements/long seeks and the first Source ProRes matrix per Mac are
reused with their original provenance. Earlier strict picture-skip failures remain
unchanged and are separately assessed under document 34's development decision.

The complete fixed physical programme passed on each Mac: sixteen Current ProRes,
four H.264 Source/Current and three 30-minute Current ProRes cases. All 46 fresh
cases were accepted, with zero recorded picture skips, underruns or timing-bound
failures, passing memory/lifecycle checks and unchanged audio configuration.
The largest observed paint interval was 46.262 ms; first/last mean frame-selection
skew differences were zero. Peak sampled playback RSS was 681.66 MiB on Sequoia
and 151.27 MiB on High Sierra, within the existing limits. These are internal
clock/frame-selection observations, not optical/acoustic synchronization proof.

The separately reused first Source matrices contain 32 cases, including two
original Sequoia strict failures with three skipped pictures in total. They are
accepted under the new development decision without changing their original
strict verdicts. This does not identify or fix the earlier paint-skip cause.

The [final summary](../../Tests/Reference/results/m7-completion-2026-10-25/qualification-summary.json)
and [independent physical audit](../../Tests/Reference/results/m7-completion-2026-10-25/fresh-physical-audit.json)
close all required gates. Both Macs' preserved M6 and M7 candidate-01/02/03
packages passed their identity/artifact audits. All 1,184 prior evidence files,
the video working copy, unrelated Finder metadata and all ten local/remote tag
objects and targets remain unchanged. Test processes have exited and their devices
are closed. Full physical eight-channel spatial listening remains **pending** and
does not block this approved offline/automated scope.

After M7, prioritize the practical stereo/5.1 workflow. The recommended next bounded
milestone is [reproducible sessions and safe reopening](35_PostM7BroadcastWorkflowProposal.md),
followed by verified PCM export/artifact QC and useful reports. It is a proposal,
not additional M7 implementation. Further immersive expansion remains deferred.

## Checkpoint publication — 25 October 2026

The user authorized the commit, dated annotated tag, merge and pushes. Freshly
fetched main remained `5e07bbd29930db34236d6dfcede7fea14983f90a`, already incorporated
by candidate-03. Main therefore advances by fast-forward to the accepted M7
checkpoint. There are no integration fixes or tested implementation changes.

The [publication preflight](../../Tests/Reference/results/m7-publication-2026-10-25/publication-preflight.json)
records fresh verification of all 165 source inputs, all 365 package artifacts,
the 317-file completion evidence manifest, 510 historical M7 evidence files and
all 46 fresh physical reports under the approved acceptance policy. Both actual-
platform qualification records remain complete. This read-only review reuses the
valid execution evidence; no device was reopened, binary rebuilt or test repeated
to obtain zero picture skips.

The checkpoint annotation binds its exact commit to the frozen source manifest,
binary hashes and evidence. All ten earlier tag objects and target commits are
preserved. The established logic working copy and separate video working copy
remain in place. Downloaded/generated media, raw captures, build products,
connection details, Finder metadata and unrelated files are excluded.

The qualification summaries retain their original pre-publication state, including
`published: false`; they are historical evidence, not rewritten publication receipts.
The publication record and checkpoint above establish the subsequent saved state.
Any future implementation or integration change requires affected qualification;
this checkpoint does not prequalify a changed merge result.

## Historical status through 24 October 2026

The following records describe candidate-01/02 under their original strict
zero-skip gate. They are historical results, not the current completion decision.
Offline and injected-device qualification passed on both actual Macs. They do
not substitute for the remaining required physical video checks.

Candidate-01 failed its real three-hour final-sample audition seek on both Macs.
The short automated suites and completed measurements remain recorded evidence,
but do not establish overall qualification. Replacement candidate-02 is frozen
with the bounded FLAC seek correction described below. Its complete automated
runs, including real three-hour measurements and long seeks, passed on both
actual Macs. Neither candidate is an accepted M7 checkpoint.

Candidate-02's first Sequoia muted Source ProRes matrix failed two of sixteen
cases because of intermittent skipped UI paints, while audio continuity and
frame-selection bounds passed. A bounded M6/M7 comparison passed, but did not
establish a cause or fix. A separately recorded controlled repeat passed all
fifteen generated cases but failed the reference movie with one skipped UI frame.
No implementation change or relaxed gate separates these attempts. Further
Current/H.264/endurance device groups are deferred. The already-started High
Sierra Source matrix completed normally: all sixteen cases passed. Original
failures remain preserved.
See the
[paint investigation](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/video-paint-investigation.md).

The same reference-movie failure was subsequently reproduced with the preserved
M6 binary and stereo audio (two skipped UI frames), and M7 stereo (one). This
does not require M7 eight-channel processing, but the underlying cause remains
unknown. The user chose to keep M7 frozen and hand the issue to the Video
Implementation chat. [Document 33](33_M7VideoPaintTimingHandover.md) contains exact
reproduction commands, evidence and the return path to a new M7 candidate.

A subsequent user-approved local isolation test after reboot passed all six
100-second stereo runs (three preserved M6, three frozen M7), with Golem
disconnected and no Screen Sharing or SSH client running. The
[post-reboot addendum](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/post-reboot-local-six-01/README.md)
records this additional evidence. It does not isolate the cause, replace earlier
failures or complete the remaining physical groups. M7 remains frozen and is
not ready for publication.

A second fixed six-run comparison in the same boot session also passed all
three M6 and three M7 stereo runs with Screen Sharing connected to idle Golem and
no SSH. The [connected comparison addendum](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/screen-sharing-local-six-01/README.md)
verifies matched conditions, connection observations and frozen integrity. The
two passing samples do not identify the original cause or complete the deferred
M7 qualification; all earlier evidence remains intact.

A third fixed six-run set with Screen Sharing and user-confirmed idle SSH also
passed all three M6 and three M7 runs, with matching recorded configuration and
unchanged source/packages. The [idle-SSH comparison addendum](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/screen-sharing-idle-ssh-local-six-01/README.md)
compares all eighteen passing post-reboot stereo runs. This does not identify the
original cause or resolve the deferred qualification; M7 remains frozen.

One additional frozen M7 stereo test after user-reported settings changes also
passed. The [single-run addendum](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/settings-change-m7-single-01/README.md)
records its separate scope. The later [user clarification](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/settings-change-m7-single-01/user-settings-clarification-2026-10-24.md)
identifies Firefox Local Network permission being granted and browser tabs
opened before that run. The user's report that the permission resets after reboot
was not independently diagnosed. This pass does not complete deferred
qualification or establish a cause/fix.

## Scope and preserved boundaries

The approved scope is [31](31_M7ChannelBasedPcm71Proposal.md), including the two
additional regression cases. The original M7 `origin/main` baseline was
`b57a1ed23cbe7d74afe98ed16ec8fd615db9d171`, the qualified M6 checkpoint. Work remains
in the established `codex/logic-module` working copy. Existing tags, video worktree,
installed apps/plugins, Pro40 and preserved qualification packages are unchanged.

The eight ordered roles are L/R/C/LFE, distinct M±090 sides and M±135 rears.
Channel count alone never establishes this layout. Native FLAC's normative order
and WAVE mask `0x63f` are L R C LFE Lrs Rrs Lss Rss. Unlabelled eight-channel WAV
and AIFF require explicit interpretation. Presets and bounded typed aliases do
not bypass Confirm roles. Unsupported metadata is identified; inconsistent and
malformed channel masks fail import. The original metadata and explicit manual
interpretation are retained in the in-memory session. No PCM reordering occurs.

BS.1770-5 Annex 3 coefficients are 1.00 for front and System I rear roles, 1.41
for sides, zero for LFE loudness. Every channel, including LFE and muted monitor
channels, enters sample/true peak and the unchanged static-gain operation.
Legacy WAV/FLAC 5.1 back/surround labels retain their qualified 1.41 interpretation.

Capability: `pcm-7.1-bs2051-system-i-v1`; ordered role mapping:
`ordered-pcm-roles-m7.1`; measurement method: `bs1770-5-ebu2023-m7.1`.
Measurement identities include the method/mapping/capability and ordered roles.
The EBU R128-2023 programme profile definition/version/digest is unchanged;
compatible full-programme measurements and explicit applicability are required.
No delivery-layout acceptance, metadata verification or future encoded export is
implied. Existing narrow M2 wording remains authoritative.

## Historical evidence levels — candidate-02

| Evidence | Required for approved M7 scope | Current status |
| --- | --- | --- |
| Offline import/measurement/evaluation/proposals/processing | Yes, both actual Macs, devices closed | **Passed** on candidate-02, including all official references, both three-hour runs and the cross-platform comparison. |
| Automated audition, pre-monitor sample parity and explicit routing | Yes, injected devices, both actual Macs | **Passed**. All 218 M7 and 124 M6 captures are bit-exact across platforms; all eight long-seek captures passed on each Mac. |
| Silent UI | Yes, both actual Macs | **Passed**, with documented reuse of unchanged candidate-01 UI code and candidate-02 repeats of affected preparation/seek/lifecycle checks. |
| Native and physical H.264/ProRes behavior | Yes, both actual Macs | Native suites passed; High Sierra's physical Source matrix passed 16/16. **Required physical timing remains failed/incomplete** because the controlled Sequoia repeat skipped one UI frame and later groups are deferred. |
| Physical eight-channel routing and spatial listening | Separate | **Pending; no suitable full 7.1 setup qualified** |

Two-output monitoring explicitly selects channels and mutes the others. The UI
labels this partial monitoring with no downmix. All eight channels remain in the
processing stream. Stereo or selected-channel playback is not physical 7.1
listening qualification. All-muted timing runs on the approved built-in outputs
supply device/clock/video evidence only.

## Candidate-02 and reproduction

Candidate-02 is frozen for replacement qualification, **not yet accepted**. Its
source-manifest digest is
`c3df888f567a69adb83cbb6c7d2c7a1dcdf1da9498d43d688ffe890786c0e43e`;
measurement build identity is
`4c12ce29addf7149641df1af27f136b933cf575c69757e91f61ee39c9a7e144e`.
Package `identity.json` SHA-256 is
`ca3eda67d2ab4505237075a0f9ff41968ff88f0d5e615652891f886c7a94d59b`;
the universal Standalone executable SHA-256 is
`c1686278deb2ea36b926ccb5cb9df10b169bec903594e38680fc5cbfa4817248`.
Its isolated package has the same parent location as candidate-01 below, with
the final folder named `candidate-02` on each Mac.

The [source comparison](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/source-differences.json)
records exactly four changed inputs: `StandaloneAudition.cpp`, the M7 native test,
M7 report comparator and comparator unit tests. UI, measurement, profile,
diagnostic and section implementation inputs are identical. Candidate-01's
unchanged UI observations are retained as partial reusable evidence, while
candidate-02 repeats affected preparation, seek, identity and lifecycle checks.
All automated suites, official vectors and native video suites were rerun for
candidate-02 and passed. The required physical timing qualification is incomplete
after the reproduced Sequoia failure.

### Preserved candidate-01 identity

Frozen package: **candidate-01**, created 2026-10-24. It failed the long audition
seek check and remains preserved; it is not an accepted milestone.

- Base commit: `b57a1ed23cbe7d74afe98ed16ec8fd615db9d171`.
- Source-manifest digest: `dd9f3392c3744c246951e33ca426545eb9d1aa8b97cafdf1ec0a2d21e73cad45`.
- Measurement build identity: `ab83aeed5e0845c17a63203cd1bb5581ace3aa582b6eaf77a79eb4d28f3f99f7`.
- Package `identity.json` SHA-256: `61aceec3c430e5d57abe06b0c7d9b39f0a9e6514efb4c43af8ac58c75cfa76fe`.
- Universal Standalone executable SHA-256: `7c42a0643566cd8aa35c23882e7f3cdcee39abab2e2ba36af53235c0b583cdc2`.
- Universal M7 test executable SHA-256: `ce3fb627d64a61d0914f267d1888d5faf4997e3b2158a239e5a9a2b4423eea78`.

The exact package is retained locally at
`build/m7-qualification-2026-10-24/candidate-01` and copied into Golem's new
`/Users/fuzz/LevelScope-M7-qualification-2026-10-24/candidate-01` folder.
Actual execution is arm64 on Sequoia 15.7.9 (24G830), Apple M1 Max / Mac13,1,
and x86_64 on High Sierra 10.13.6 (17G66), Intel Core i7-3770K at 3.50 GHz.
Both platform guards passed and both artifact-integrity checks reported no
mismatches. Complete artifact hashes are retained in each candidate identity;
candidate-02's final integrity and shutdown records are linked below.
A source manifest binds all implementation and runner inputs, including test-only
Objective-C reader auditing. Its ancestry commit is not a substitute for the
working candidate's exact source digest. Downloaded/generated media and binaries
stay in ignored isolated qualification folders.

Build only the explicitly named standalone/test targets in a new build directory:
`LevelScopeStandalone`, `LevelScopeMeasurementTests`, `LevelScopeScaffoldTests`,
`LevelScopeComplianceTests`, `LevelScopeProposalTests`, `LevelScopeDiagnosticsTests`,
`LevelScopeSectionsTests`, `LevelScopeAuditionTests`, `LevelScopeMultichannelTests`,
`LevelScopeVideoTests`, `LevelScopeVideoQualification`. Never use the default/all
or VST3 target. Pinned JUCE is 7.0.12; minimum deployment 10.13; universal x86_64 and
arm64; the established non-fast-math/FP-contract settings remain in force.

New tools are `Tests/package_m7_qualification.py`, `Tests/run_m7_qualification.py`,
`Tests/compare_m7_reports.py` and `Tests/run_m7_muted_video.py`. The physical runner
requires explicit coordinated `--allow-audio`, named device and muted outputs.
Invocation records are retained in each run's `runs.json`; the frozen package
contains the exact scripts used. The evidence audit under the results directory
is separate post-run verification, not an input to the frozen executable.

The following candidate-02 reproduction notes retain their historical scope.
For combined candidate-03, use the updated explicit options and thirteen-target
recipe in the [completion evidence README](../../Tests/Reference/results/m7-completion-2026-10-25/README.md).

For a fresh rebuild from the recorded candidate-02 source, use the pinned dependency
and name every target explicitly (the example creates a separate reproduction build):

```sh
cmake -S . -B build/logic-m7-repro \
  -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_OSX_ARCHITECTURES=x86_64;arm64' \
  -DCMAKE_OSX_DEPLOYMENT_TARGET=10.13 \
  -DFETCHCONTENT_SOURCE_DIR_JUCE=/Users/fuzz/Developer/LevelScope-dependencies/JUCE-7.0.12
cmake --build build/logic-m7-repro --config Release --parallel 4 --target \
  LevelScopeStandalone LevelScopeMeasurementTests LevelScopeScaffoldTests \
  LevelScopeComplianceTests LevelScopeProposalTests LevelScopeDiagnosticsTests \
  LevelScopeSectionsTests LevelScopeAuditionTests LevelScopeMultichannelTests \
  LevelScopeVideoTests LevelScopeVideoQualification
```

The frozen build used Apple clang 17.0.0 (`clang-1700.0.13.5`). A rebuild must
record its own artifact hashes; identical source does not by itself establish
binary identity or execution evidence. A future checkpoint commit should contain
all 155 source-manifest inputs byte-for-byte and add the final documentation and
evidence. The manifest's `base_commit` remains the M6 ancestry commit: do not
rewrite the preserved candidate metadata to pretend it was built from a later
checkpoint commit.

From the frozen package root, run `python3 -B Tests/run_m7_qualification.py
--platform sequoia --output /absolute/fresh/output` on Sequoia. On actual High
Sierra use its inspected `/usr/local/bin/python3` and `--platform high-sierra`.
The runner rejects the wrong OS/architecture or changed package artifacts before
testing. It runs all nine native suites, the three-hour audition regression,
both official reference suites, runner tests and the real three-hour 7.1 PCM
measurement. Keep sample captures in the ignored run folder for subsequent
`Tests/compare_m7_reports.py` comparison; the published report records their hashes.

The generated fixture manifest is
[`m7-generated-fixtures-manifest.json`](../../Tests/Reference/m7-generated-fixtures-manifest.json).
Its real three-hour source contains 518,400,000 eight-channel frames at 48 kHz,
stored as 629,476,570 bytes of native FLAC. Preparation and measurement each own
a verified encoded copy; the combined long test retains two such copies at its
peak (1,258,953,140 bytes, excluding the original fixture). PCM buffers are bounded
and do not expand to programme length. Existing free-space checks and injected
write/create/flush-failure regressions remain applicable.

## Implementation structure

`StandaloneMeasurementLayout.h` owns the bounded aliases, named file-order
presets, exact role-set capability and monitoring-coverage wording. Import keeps
original decoder metadata and records format interpretation separately from
explicit manual interpretation. The native FLAC override parser reads bounded
metadata with cancellation checks; it does not allocate arbitrary comment blocks.

`OfflineAnalysisService` validates the supported layout before measurement or
Current audition and binds the mapping/capability version into each snapshot.
The kernel adds explicit rear roles without changing its filters, gating,
true-peak reconstruction or the shared `StaticGainProcessor`. M2 and M3 consume
the expanded compatible measurement capability; their target profile and policy
remain unchanged. Diagnostics accept validated peak vectors within the existing
twelve-channel capacity. No new six-channel assumption is introduced into
processing, audition or routing.

The main window reserves a scrollable channel viewport, and interleaved imports
now use “Channel N” labels in routing and peak views. Role presets still require
confirmation. Typed aliases are limited to the documented System I vocabulary;
ambiguous generic 7.1 `Ls/Rs` entries are not silently treated as sides or rears.

### Candidate-01 failure and follow-up

Both three-hour measurements completed on each Mac, but final-sample Current audition
remained in Buffering beyond the ten-second harness deadline. An isolated reader
probe reproduced a 33.26-second direct seek to frame 518,399,999 in the pinned
JUCE 7 FLAC decoder, while the other probed seeks took under 1 ms. This was not a
measurement mismatch. Candidate-01's failure is retained and is not relabelled.

The follow-up auditions a FLAC seek only after decoding and discarding at most
65,536 preceding frames in existing 2,048-frame buffers, with cancellation checks
between reads. This exceeds the maximum FLAC block size documented in
[RFC 9639 §4.1](https://www.rfc-editor.org/rfc/rfc9639.html#section-4.1).
Discarded samples never enter processing, capture, output or the presentation
clock. The verified copy is reused, the requested first sample is unchanged, and
other source formats retain their existing path. The isolated reader probe reduced
the three tested near-EOF seeks to 5–7 ms. The focused product regression on
Sequoia passed all 51 assertions: eight Source/Current captures were exact,
prepared seek buffering took 2.6–7.7 ms, cancellation passed and the verified
copy was retained. The focused candidate-02 run on actual High Sierra also passed
51 assertions, with prepared-seek buffering of 2.45–8.89 ms. Candidate-02's complete
Sequoia long run passed 2,411 assertions: Source measurement took 937.90 s, Current
931.46 s, and all eight long-seek captures were exact with 3.07–7.58 ms buffering.
The full High Sierra run also passed 2,411 assertions: Source measurement took
1,690.26 s and Current 1,681.45 s; its eight captures were exact with 2.72–9.75 ms
prepared-seek buffering. Source and Current integrated loudness, true peak and
LRA matched Sequoia exactly. Additional owned audition workspace was 2.69 MiB on
Sequoia and 2.04 MiB on High Sierra, within the 16 MiB limit.
The expanded long regression checks all eight captures at four positions in both
Source and Current, a two-second prepared-seek limit, and stopped cancellation.

## Historical qualification checklist — candidate-02

Completed machine-readable evidence is collected under
[`m7-qualification-2026-10-24`](../../Tests/Reference/results/m7-qualification-2026-10-24/).
This directory records completed partial qualification and the outstanding video
blocker; it is not an accepted checkpoint. The full cross-platform comparison
records zero float32 ULP difference in all
218 M7 capture cases and zero difference in their reported metrics. The separate
M6 comparison records zero float32 ULP difference in its 124 regression captures
and a maximum reported metric difference of 1.78e-15 dB. M2–M5 comparisons passed
their established tolerances. The three-hour M7 metrics matched exactly.
Each final M7 native suite passed 45,377 assertions, including 40,320 ordered role
permutations and the capture-integrity assertions enabled by the qualification runner.

For candidate-02's synthetic three-hour/192 kHz, 1,000-section benchmark, Sequoia /
High Sierra diagnostic worker times were 145.38 / 320.64 ms, paint p95 5.85 /
26.73 ms, and extra peak RSS 18.95 / 47.80 MiB. Diagnostic storage was 20.06 MiB and section
workspace 1.67 MiB on each platform. These are downstream analysis benchmarks;
the completed actual three-hour PCM runs have separate reports.

- [x] Complete mapping/alias/metadata and 40,320 role-permutation checks.
- [x] Real WAV/FLAC legacy 5.1 interpretation and 1.41 numerical regression.
- [x] Rear/LFE unique highest-peak mute regression through M2/M3 and capture.
- [x] Full format/rate/gain matrix and all-eight-channel same-platform sample parity.
- [x] Offline complete workflow proves zero device opens.
- [x] M2/M3 freshness, no-change, strict boundaries and failed-verification rules.
- [x] M4/M5 eight-channel peaks, source coordinates and performance limits.
- [x] Actual three-hour 48 kHz eight-channel native FLAC on both Macs.
- [x] Synthetic three-hour/192 kHz, 1,000-section responsiveness/memory benchmark.
- [x] All M1–M6 regressions and both pinned official EBU/ITU suites.
- [x] Native H.264/ProRes source isolation, lifecycle and exact timing checks.
- [ ] Coordinated muted physical short matrix and three endurance runs per Mac: **blocked by a reproduced Sequoia Source paint-skip failure; later groups deferred**.
- [x] Silent interactive UI on Sequoia and actual High Sierra (candidate-02 affected checks plus documented reuse of unchanged candidate-01 UI code).
- [x] Full cross-platform M2–M7 comparison, including long M7 measurements.
- [x] Final candidate-integrity audit: all 155 source inputs and all 339 artifacts match; both actual-platform packages checked.
- [x] Existing ten tag objects/targets and unrelated work preserved; local and remote refs checked without mutation.

The final [qualification summary](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/qualification-summary.json)
keeps the overall result false and links the completed evidence. The final
[physical-video audit](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/muted-video-audit-final.json)
passes 31 of 32 completed reports across the controlled Sequoia and High Sierra
Source matrices, but correctly fails for the skipped frame and missing deferred
groups. First-attempt and baseline-comparison failures remain separate evidence.

[Local integrity/preservation](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/preservation-final-candidate02.json),
[High Sierra integrity/shutdown](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/high-sierra-integrity-final.json),
[Sequoia shutdown](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/sequoia-shutdown-final.json)
and [remote ref preservation](../../Tests/Reference/results/m7-qualification-2026-10-24/candidate-02/remote-preservation-final.json)
all passed. All qualification processes have exited; completed physical harnesses
closed their devices before writing reports. No further device groups are
scheduled. No M7 commit, tag, push or merge was made. The source, packages and
failed evidence remain frozen for the video handover.

## Limitations and next boundaries

The implementation keeps twelve-channel storage, height roles and count-driven
processing/routing interfaces. Only System I 7.1 is added to measured Current
support. 7.1.4 mapping/weights/profile applicability and its qualification remain
separate. There is no Atmos object/rendering support, implicit downmix, bass
management or LFE +10 dB listening gain.

RF64/BW64 and other large-file container extensions are not qualified by M7.
The real three-hour fixture uses native FLAC within its existing frame-count
limit; it is not a new large-WAV claim. Separate-file import, section processing,
dynamics, persistence, public export and VST3 parity remain deferred. Sections
remain read-only processing-wise and in-memory. Monitor volume is separate and
live; processing/domain/roles/routes remain stopped-only. Playback of unmeasured
Current State stays UNVERIFIED, and does not establish compliance.

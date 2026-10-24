# M6 qualification and handover — candidate 01

Updated 24 October 2026. **Candidate 01 is qualified on actual arm64 Sequoia and
Intel High Sierra, published and integrated into main.** All required automated,
interactive and coordinated physical checks passed. The annotated checkpoint is
[`m6-verified-current-state-audition-2026-10-24`](https://github.com/carlozaiana/LevelScope/tree/m6-verified-current-state-audition-2026-10-24).
Its peeled commit is the accepted M6 checkpoint and main revision at publication;
the annotation records the full commit SHA and frozen source/binary hashes.

## Candidate identity

Established working copy: `/Users/fuzz/Developer/LevelScope`, branch
`codex/logic-module`. Base/published main at freeze:
`6d1c22a23244b4e421dd2773862451ac8a863992`. This identifies ancestry; the M6
implementation is identified by the checkpoint above and the frozen manifest below.
The candidate manifest retains its original pre-commit working-candidate identity.

The immutable package is `build/m6-qualification-2026-10-23/candidate-01`, copied
to a new `~/LevelScope-M6-qualification-2026-10-23/candidate-01` on Golem.

| Identity | SHA-256 |
| --- | --- |
| Exact 143-file source/test/manifest definition | `0702768309e7328549fb86b677ec295fdc197d0382c45f0487b9060053b9a286` |
| Package identity.json | `ab7440f4d06275ea55055f59eec98c4f5872a775342e5b31777cefdf1dfdb4dd` |
| Compiled measurement build fingerprint | `afe72a024daea1a94ff298567bb185d2c054dcdaea61a38bcdfc36b9a1df5e06` |
| Standalone executable | `caf2b25bc13cfeee16cb59f17684227cbf1cb4822f84126bc18b78cba76024b6` |
| Audition test executable | `6e3490a15231beea07a5d7bbe9677c70837e82c76210aadd09ce041524714f2d` |
| Video qualification executable | `c8efc1a6d0b760206c8f99105e553542a4d32fdb9722dc3e19b2172f60bb36e9` |

The package identity contains all executable, fixture, pinned vector and source
hashes. Documentation updates after freeze do not change the tested implementation.
Every one of the 143 source hashes matches the checkpoint's corresponding file.
The base commit is not substituted for the source manifest, and earlier development
binaries are not relabelled as this candidate. Resolve the exact accepted revision
with `git rev-parse 'm6-verified-current-state-audition-2026-10-24^{commit}'`.

Release universal arm64/x86_64 build, minimum macOS 10.13, pinned JUCE 7.0.12
`4f43011b96eb0636104cb3e433894cda98243626`, FP contraction disabled and no fast-math.
Only the explicitly named standalone and test targets were built. No plugin or
installation targets were built or run.

## Implemented behavior

See [approved scope and acceptance criteria](28_M6VerifiedCurrentStateAudition.md).
Source uses identity processing. Current State uses the exact stored double trim,
through the same `StaticGainProcessor` used by offline measurement. Every ordered
channel, including LFE, receives one multiplication. The display's twelve decimal
places do not feed the processing plan.

An explicit device-closed preparation action copies, hashes and fully validates
the imported content before audition. The private copy is reused across seeks and
stopped domain/gain changes. Copy/hash/payload validation reports progress and
supports cancellation. Disk-space, create, write and flush failures clean owned
partials and preserve the previous valid stopped session. Same-size/timestamp
mutation fails the content hash. Playback uses the verified copy, not mutable
external file bytes.

The audition request binds session, source hash/metadata, full-programme range,
ordered confirmed roles, domain, exact processing plan and Current revision.
Domain/processing/roles/routing/M3 Apply changes are stopped-only in both UI and
controller. Monitor volume remains live and smoothly ramped, with no processing or
measurement invalidation. Profile, video offset and section annotations remain
outside audio processing identity.

Valid unmeasured Current State is allowed and labelled **UNVERIFIED**. The UI
distinguishes measured freshness from the separate M2 verdict; an actual measured
failure remains a failure even though its samples can be auditioned. Preparation,
playback and predictions do not establish a measured pass.

Finite processed samples above full scale are retained. The separate monitor
output policy prechecks the entire routed callback; an excess outside [-1,1] or a
non-finite output silences that whole callback and stops without consuming its
source span. Explicit restart resets FIFO and monitor ramp, preserving the gain
plan and held position. Lowering volume alone does not resume. This is not an
intersample/analog true-peak protection claim.

The decoder/FIFO/capture/routing interfaces retain twelve-channel capacity and
ordered roles. Current State explicitly rejects layouts/rates outside M1; broader
Source audition is retained. Future eight-channel 7.1 and twelve-channel 7.1.4
qualification is specified separately in [29](29_MultichannelPcmMilestone.md).

## Automated evidence completed

Both runs passed platform and package-integrity guards and executed native slices:

| Platform | Runtime | M6 checks | Official references and regressions |
| --- | --- | --- | --- |
| Sequoia | macOS 15.7.9, build 24G830; Apple M1 Max (Mac13,1), arm64 | 4,736 assertions and separate three-hour suite passed | M1–M5, scaffold, video, all runner suites, pinned ITU BS.2217 and EBU V05 passed |
| Golem | macOS 10.13.6, build 17G66; Intel Core i7-3770K 3.50 GHz, x86_64 | 4,736 assertions and separate three-hour suite passed | Same suites passed on actual High Sierra |

The [evidence directory](../../Tests/Reference/results/m6-qualification-2026-10-23)
contains platform reports, raw-run logs, the candidate identity, UI captures and
cross-platform comparison reports. Large raw sample captures and reference media
remain only in the preserved local/remote test folders.

Environment reports record CPU, hardware model, OS, executable architecture,
dependencies, compiler/build inputs, source identity and binary hashes. The High
Sierra archive unpacker warned about unsupported newer macOS extended metadata;
all 308 regular package artifacts passed the exact content-integrity guard before
tests ran. No old Golem package or installation was overwritten.

The 124 captured audition cases were **bit-identical across platforms**, including
nonzero trim. On each platform audition and offline sample streams match exactly;
independent long-double reference checks pass their one-float32-ULP bound.
Cross-platform metric/status/stability checks pass: maximum true-peak difference
`1.7763568394002505e-15 dB`; other reported metric differences are zero. M2–M5
cross-platform comparators also pass with their existing budgets. No test tolerance
is added to the EBU profile's acceptance allowance.

Coverage includes the four M1 rates and mono/stereo/5.0/5.1, WAV/AIFF/FLAC,
zero/positive/negative/extreme/non-rounded gain, LFE, silence, gate crossing, short
impulses, finite over-full-scale samples and non-finite rejection. Tests capture
actual consumed samples before monitoring, remeasure them, and exercise seek/FIFO
generations, identity mismatches, cancellation, source replacement, device faults,
whole-callback overload and explicit restart recovery. Injected disk failures do
not fill a real disk. Injected devices are distinct from physical-device evidence.

### Three-hour input

Stereo 48 kHz, 518,400,000 frames, 2,073,600,044 file bytes; SHA-256
`cf3d1fa528ce1798225b48496a91ced8f1696cbe17dbeb14f548d5356c86076a`.
Preparation, source-frame seeks including the final second, copy reuse and cleanup
were exercised on both Macs.

| Measurement | Sequoia | High Sierra | Required budget |
| --- | --- | --- | --- |
| Full preparation | 23.836 s | 44.046 s | Recorded, no real-time claim |
| Copy/hash/payload cancellation | 3.69 / 34.17 / 24.73 ms | 42.93 / 46.60 / 48.32 ms | Each <=2 s on responsive storage |
| UI refresh p95 | 0.151 ms | 0.199 ms | <=100 ms |
| Additional peak resident memory, including FIFO | 4,063,352 B | 4,018,296 B | <=16 MiB |

The source-sized temporary disk copy is separate from the duration-independent
memory budget. OS-stalled I/O has no hard timeout guarantee.

## Interactive and physical checks

Sequoia frozen-product UI observations, all device-closed and routes muted:

- Imported stereo Source and confirmed ordered roles; prepared verified Source.
- Measured Source; initialized and selected unmeasured Current, correctly labelled
  UNVERIFIED; prepared Current without opening a device.
- Measured Current at 0 dB; edited gain to -3.305074118944 dB, discarding old Current
  measurements while preserving Source; remeasured and displayed actual M2 failure.
- Changed monitor volume and switched audition domain without losing measurements.
- Sought to 3/2 s while stopped; attached H.264 and applied -1001/30000 s offset.
  Held picture changed from frame 44 to 45, with audio measurements unchanged.
- Drafted M3 trim without applying it; explicit stopped Apply changed gain and
  displayed UNVERIFIED; real remeasurement then passed the selected rules.
- Replaced Source with the three-hour fixture; old Current cleared. Preparation
  showed progress, cancelled to a clear stopped status, and closed cleanly during
  a subsequent preparation (process exit verified).

High Sierra frozen-product UI observations, device-closed and routes muted:
imported/confirmed stereo, prepared Source, measured Source, selected and prepared
UNVERIFIED Current using the existing copy, and measured zero-gain Current. Editing
gain to -18 dB discarded Current measurements while retaining Source. Remeasurement
showed an actual loudness failure (I -40.993358 LUFS), separately from the audition
label's current-measurement status. Monitor-volume changes preserved that result.
An unapplied M3 draft remained prediction-only; explicit Apply displayed
UNVERIFIED, and real remeasurement then passed (I -23.000000 LUFS,
TP -23.005446 dBTP). Selecting Source retained the applied Current gain and result.
On 24 October the user completed Source and zero-gain Current measurements, H.264
attachment, stopped seek to 3/2 s and offset -1001/30000 s. Inspection confirmed
held frame 45, both measurements still current, device closed and both routes
muted. Clear Source removed both Source and Current. The three-hour fixture was
then imported and its roles confirmed. Preparation displayed progress; Stop
cancelled it with an explicit stopped status. Starting preparation again and
closing the window exited the candidate cleanly (process absence verified over
SSH without forced termination). Transactional source replacement and retention
on failed preparation also pass the automated product/controller tests.

The automated suite additionally exercises actual product-panel callbacks for
fresh/unavailable/stale results, live monitor control, storage-failure retention
and source-change invalidation. Screen Sharing's blue control arrow restored mouse
control; forwarded clicks/text remain intermittent, so manual assistance is
recorded separately from automated evidence.

Both platforms' native ProRes qualification passed **53 cases each**, including supported
profiles/rates, policy rejections, VFR, selected nonzero audio isolation in both
domains, actual-file isolation and large-file loading/cancellation. Both platforms'
muted named-device ProRes matrices passed all 16 Source and all 16 nonzero Current
cases per Mac. Both H.264 CFR/VFR Source/Current matrices also passed all four cases
per Mac. These runs had zero actual-paint skips, blanks, bound/deadline failures and
underruns. All three Sequoia 30-minute Current ProRes endurance cases passed,
including the memory limits; its complete physical runner exited successfully.
All three High Sierra endurance cases also passed and its runner exited with
code 0. Every endurance case completed at least 1,800 seconds of Current playback
at -3.305074118944 dB, with zero actual-paint skips, blanks, bound/deadline failures
and underruns. Preparation and loading are additional elapsed time.

| Platform / ProRes profile | Peak playback RSS (MiB) | Non-reusable resident growth (MiB) | Physical-footprint growth (MiB) |
| --- | ---: | ---: | ---: |
| Sequoia / Proxy | 670.42 | -5.15 | +3.59 |
| Sequoia / LT | 634.14 | -6.29 | +1.15 |
| Sequoia / 422 | 681.14 | -2.23 | +0.43 |
| High Sierra / Proxy | 136.32 | +1.57 | +1.64 |
| High Sierra / LT | 143.63 | +0.40 | +0.47 |
| High Sierra / 422 | 148.94 | +0.95 | +1.09 |

All meet document 25's unchanged 768 MiB peak RSS and 32 MiB growth limits, using
the accepted `prores-memory-v2-user-approved-2026-10-20` method. Growth compares
the first and last minute means; a negative value is a decrease. Detailed reports
retain full precision. The existing native queue-bound regressions also passed.

The [physical evidence summary](../../Tests/Reference/results/m6-qualification-2026-10-23/physical-evidence-summary.json)
checks all 78 physical cases (64 ProRes matrix, eight H.264 and six endurance)
against candidate identity, domain, exact gain, device/rate, duration and acceptance
fields. It reports complete/pass with no missing reports or discrepancies.
No physical optical/acoustic synchronization measurement is claimed.

The user approved Mac Studio Speakers at 48 kHz and Golem Built-in Output at
44.1 kHz with all routes muted, including the six 30-minute ProRes runs. Separate
brief listening checks are also approved; the user confirmed readiness and changed
the requested listening monitor level to -15 dB. Eight-second fixtures bound each
play by EOF. The user confirmed Current was clearly quieter on both Macs, with no
clicks/dropouts, and remained labelled UNVERIFIED. This is user-reported listening
evidence for the requested Source/Current comparison, separate from exact sample
parity and instrumented synchronization evidence. Sequoia's final observed UI
showed Current -6.020599913280 dB, monitor -15 dB, eight-second EOF and zero
underruns. Golem's later altered levels were an additional deliberate contrast test,
as clarified by the user; they are not relabelled as the requested settings.
Golem used the initial eight-second sparse-pulse fixture. A clearer continuous-tone
fixture was prepared, but its proposed repeat was interrupted and is not claimed.
Timing-harness runs retain their fixed -30 dB monitor level and muted routes.
Device opening does not authorize changes to
system defaults, Pro40 or installed applications/plugins.

### Pause and resume record — 23 October 2026

The user requested a stop to attend to other work. The paused Sequoia orchestration
driver and its temporary sleep-prevention process were terminated; no physical
qualification runner was active on Golem. Tests did not resume automatically.
The completed reports and frozen candidate remain intact. Sequoia's test output
device was explicitly closed. Golem was observed stopped at eight-second EOF;
Screen Sharing then became unavailable during the close-device action, so remote
device closure was not confirmed at that pause. The stopped/device-closed state
was verified on the authorized resumption before further UI work.

At the pause, High Sierra UI lifecycle/video-identity checks, both H.264 matrices,
Golem's ProRes matrices and all six endurance cases remained outstanding. Those
checks were completed on 24 October. Sequoia's two earlier ProRes matrices were
reused for the unchanged candidate. The original orchestration script starts with
those matrices and uses exclusive output creation; the separate remaining-cases
helper wrote new outputs and preserved earlier evidence. No commit, tag, push or
merge was performed.

On 24 October the user resumed qualification. All 143 source files and all 308
package artifacts matched on Sequoia; all 308 preserved artifacts matched on
Golem. All nine prior tag objects/targets and both worktree registrations remained
unchanged. No qualification runner was left active from the previous day.
Sequoia's outstanding muted tests were started with the separate
`run_remaining_muted_video.py` helper and new outputs; completed matrices were not
rerun or overwritten. The user connected Saffire Pro 40 to Golem for a separate
listening repeat and confirmed outputs 1/2 and current 48 kHz. The existing 48 kHz
eight-second fixture was used at -15 dB monitor level after coordination. This
was separate from the built-in-output timing matrix; Pro40 drivers, services and
system audio configuration were not changed.

The 24 October Pro 40 supplement passed by user report: Source was louder than
Current at -6.02059991328 dB, both plays were free of clicks/dropouts, and Current
remained UNVERIFIED. The assistant inspected the final Current label
(-6.020599913280 dB), -15 dB monitor level, closed device and muted routes. The
existing 48 kHz eight-second stereo fixture was used with outputs 1/2. This remains
listening evidence, separate from instrumented parity and synchronization tests;
no agent action changed Pro40 controls, drivers, services or system settings.

## Qualification decision and publication handover

**M6 meets the approved acceptance criteria. No required qualification check is
outstanding.** The [overall status](../../Tests/Reference/results/m6-qualification-2026-10-23/status.json)
joins the two native automated runs, cross-platform comparisons, silent UI records,
coordinated listening observations and complete physical evidence. Individual
automated-run summaries retain `qualification_complete: false`: those historical
subset reports correctly say that automated suites alone do not qualify M6.

The [final preservation audit](../../Tests/Reference/results/m6-qualification-2026-10-23/final-preservation-audit.json)
passed at qualification completion, before publication. All 143 working source files match the frozen manifest; all 308 package
artifacts match on each Mac, including Golem's 143 packaged source files. Candidate
processes have exited on both systems and their test devices have been released.
At that audit all nine pre-existing annotated tag objects and target commits, the
logic HEAD and the registered worktrees were unchanged. The video working copy, installations,
plugins, Pro40 repository/drivers/services and system defaults were not modified.

## Checkpoint publication — 24 October 2026

The user separately authorized the checkpoint, annotated tag, merge and pushes.
Fetched main remained `6d1c22a23244b4e421dd2773862451ac8a863992`, the exact
qualified base. Publication adds the intended M6 source, tests, documentation and
curated evidence on `codex/logic-module`; main incorporates it by fast-forward.
There are no integration fixes or changes to the tested implementation. Existing
two-platform qualification is reused; no device is reopened or package rebuilt.

The [publication preflight](../../Tests/Reference/results/m6-qualification-2026-10-23/publication-preflight.json)
records a fresh 143-source/308-package/394-original-evidence integrity check, a
fresh comparison of all 124 retained sample captures, a fresh audit of all 78
physical reports, and the original local/remote tag objects and target commits.
The checkpoint's annotation records its exact commit and identities. All nine
earlier tags retain their original objects and commits. The video branch and
working copy remain unchanged. Downloaded/generated media, raw captures, build
products, connection details, Finder metadata and unrelated work are excluded.

The frozen identity and historical qualification audits are retained verbatim.
Publication records do not rewrite their pre-commit ancestry or claim a new binary.
Any later implementation change requires affected qualification on both Macs;
this checkpoint does not prequalify a future merge.

## Reproduction

Configure `build/logic-m6` for Release, architectures `x86_64;arm64`, deployment
target 10.13, and `FETCHCONTENT_SOURCE_DIR_JUCE` pointing at the pinned checkout.
Enable the measurement/scaffold/compliance/proposal/diagnostics/sections/audition/
video test options. Build only these targets explicitly:

```text
LevelScopeStandalone LevelScopeMeasurementTests LevelScopeScaffoldTests
LevelScopeComplianceTests LevelScopeProposalTests LevelScopeDiagnosticsTests
LevelScopeSectionsTests LevelScopeAuditionTests LevelScopeVideoTests
LevelScopeVideoQualification
```

From an intact frozen bundle, use Python 3.8+:

```sh
python3 -B Tests/run_m6_qualification.py --platform sequoia --output /new/results
# On actual High Sierra, use --platform high-sierra and its native Python runtime.
```

Every output must be a new directory. `compare_m6_reports.py LEFT RIGHT --output
REPORT` verifies raw captures and platform guards; the existing M2–M5 comparators
check their respective results. Raw float captures, fixtures and build products
are local test artifacts and must not be committed. Physical video runners require
the explicitly coordinated named device and visible test window; silent automated
suite success alone never marks M6 qualified.

## Limitations carried forward

Static global trim only; stopped domain/processing/roles/routing changes; exact
device rate, no resampler or downmix. No seamless A/B, section processing, dynamics,
automatic correction, persistence, export, embedded soundtrack or VST3 parity.
Sections remain in-memory annotations. Temporary copies are ephemeral; abrupt
termination can leave owned temp files under the pre-existing cleanup limitation.
Future interleaved labels should read “Channel 1, Channel 2…”; future separate-file
imports should show filenames beside roles. Neither import expansion nor that
wording change is included here. Informal user observations of dropped video frames
remain separate from instrumented evidence and the Video Implementation chat.

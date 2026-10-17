# LevelScope standalone scaffold handover

Prepared 14 October 2026 for a fresh development chat in `/Users/fuzz/Developer/LevelScope`.

Updated 17 October 2026: M1 qualification **passes** on actual Intel High Sierra
10.13.6 (17G66) and arm64 Sequoia 15.7.9 (24G830), including both official suites and
both silent interactive UI checklists. The checkpoint is
`measurement-milestone-2026-10-17` on `codex/logic-module`. See
[17_M1QualificationResults.md](17_M1QualificationResults.md) for final evidence and
[15_VideoImplementationHandover.md](15_VideoImplementationHandover.md) for the next chat.
Historical scaffold observations below remain scoped to their original dates.

## Current integration handover (17 October 2026)

The qualified measurement checkpoint is now on main; video has completed its final
30-minute gate on both systems. Read [the main integration handover](20_MainIntegrationAcceptance.md)
for current acceptance and the subsequent logic-module starting point. Keep VST3
measurement parity separate. The older plans and limitations below are historical.

## Logic-module continuation (15 October 2026)

The first measurement implementation is now on `codex/logic-module`, created from
the scaffold tag. See [14_OfflineMeasurementContract.md](14_OfflineMeasurementContract.md)
for accepted measurement boundaries, channel/LFE policy, identity, tolerances and
video/transport design. The sections below describe the preserved scaffold checkpoint;
this section records subsequent work, not a change to that checkpoint.

Implemented in the continuation:

- Separate worker-only `LevelScopeMeasurementCore` with published 48 kHz K coefficients,
  qualified-rate coefficient design, exact integrated block gating, independent LRA
  population/gate, M/S curves and sample-resolution maxima, all-channel sample peaks,
  and 8x/64-tap-per-phase Blackman-sinc true-peak estimation. Existing plugin analysis,
  parameter IDs, APVTS and processing behavior remain unchanged.
- Mono, stereo and conventional 5.0/5.1 measurement at 44.1/48/96/192 kHz. Unqualified
  roles/rates are refused for measurement while the broader import contract remains.
  Confirmed file order remains authoritative; no Core sorted-layout adaptation is used.
- Source and Current State measurement buttons, cancellable jobs, per-metric missing
  values/reasons, numeric gates, short-program LRA instability indication, and Current
  State static gain (-60 to +24 dB) applied to all file channels including LFE. This
  gain does not change original-source audition or monitor gain.
- A private source snapshot, SHA-256 matched to the imported bytes and revalidated
  before decoding. Same-size/same-time edits fail the hash check. Temporary snapshots
  are cleaned up on ordinary success/error/cancel. Abrupt termination may leave files.
- Session UUID and content/role/range/plan/engine/build/platform measurement identity.
  Current edits invalidate measured values; stale worker results cannot publish.
  An export-approval identity-binding primitive rejects changed revisions, plans,
  output hashes or profile editions/limits. This is not a compliance verdict or a
  public exporter; those services remain unavailable.
- An internal PCM24 render/reopen/measure harness verifies actual encoded samples,
  rejects clipping and uses only owned temporary paths. Public output path/alias
  validation is still required when public rendering is implemented.

Evidence so far: 184 new hardware-independent checks plus the original 205 scaffold
checks pass on arm64 Sequoia. All 39 downloaded ITU BS.2217 vectors pass their published
file-based integrated-loudness tolerance of ±0.1 LU; maximum observed error is about
0.04845 LU. The new tests also independently synthesize multirate filter/level signals,
Tech 3342 cases 1–4 and Tech 3341 true-peak cases 15–23, plus content-probed AIFF/FLAC
snapshot decoding. Synthesized tests are labelled separately from the original
downloaded files.

**Official EBU v05 qualification (15 October 2026):** the user supplied the archive
under `Tests/Ebu_Loudness_Testset_V05` after the automated download failed. All 70 WAVs
were measured on macOS 15.7.9 arm64. **66/66 minimum-requirements files pass** (Tech
3341 cases 1–23 and Tech 3342 cases 1–6; 72 metric rules, 1,052 checked readings).
The remaining **4/4 supplementary alignment/listening-level files pass** project
sanity checks at ±0.1 LU; that tolerance is not an additional official EBU criterion.
Eight tests of the new runner's failure detection also pass. No measurement-kernel
correction was necessary; changes for this qualification are in the test harness,
manifest, runner and evidence. The kernel remains `bs1770-5-ebu2023-m1.1`.

Maximum absolute errors across normative checks: I 0.04645 LU, M 0.04648 LU,
S 0.04646 LU, LRA 0.02339 LU, TP 0.11699 dB (an under-read). Programme LRA results
are 4.97662 LU and 14.98460 LU against 5/15 LU ±1 LU. Cases 9/12 check constant
curves; cases 11/14 check all 20 successive maxima without resetting the original
engine. The latter are offline kernel checks, not live-meter UI/scheduling qualification.

The archive CRC passes and every extracted file matches its archive bytes. Archive,
README and file SHA-256 values, metadata, roles and acceptance rules are pinned in
`Tests/Reference/ebu-v05-manifest.json`; the saved report is
`Tests/Reference/results/ebu-v05-sequoia-arm64-2026-10-15.json`. The runner records
binary/build/engine/manifest/runner identity and actual execution architecture. All
official EBU files are 48 kHz. Missing/modified input or missing readings fail the run.

The universal app and named tests build with the original deployment minimums. The
17 October qualification closes the earlier interactive-UI and actual High Sierra
evidence gaps: all named regressions, official suites, comparisons and UI checklists
passed on both systems. No physical device was opened. The earlier scaffold smoke
report is separate from this new measurement qualification.

Build and run:

```sh
cmake -S . -B build/scaffold -DLEVELSCOPE_BUILD_SCAFFOLD_TESTS=ON -DLEVELSCOPE_BUILD_MEASUREMENT_TESTS=ON -DFETCHCONTENT_SOURCE_DIR_JUCE=/Users/fuzz/Developer/LevelScope/build/juce-7.0.12
cmake --build build/scaffold --target LevelScopeStandalone LevelScopeMeasurementTests LevelScopeScaffoldTests -j 6
ctest --test-dir build/scaffold --output-on-failure
python3 Tests/run_reference_measurements.py --binary build/scaffold/LevelScopeMeasurementTests_artefacts/Release/LevelScopeMeasurementTests --vectors build/measurement-reference/itu --report build/measurement-reference/sequoia-arm64.json
python3 Tests/run_ebu_measurements.py --binary build/scaffold/LevelScopeMeasurementTests_artefacts/Release/LevelScopeMeasurementTests --vectors Tests/Ebu_Loudness_Testset_V05/ebu-loudness-test-setv05 --report build/measurement-reference/ebu-sequoia-arm64.json
python3 -B Tests/EbuReferenceRunnerTests.py
```

The tests use macOS's temporary directory; a restricted shell may require approval
for that filesystem access. No audio device is needed. Download provenance, archive
and per-file hashes and expected results are in `Tests/Reference/itu-bs2217-manifest.json`.
External ITU audio stays in ignored build storage; the supplied EBU directory is also
ignored. Neither is shipped. The reference runners
fail missing/hash-mismatched fixtures; `--compare` checks another platform's report
against the proposed 0.01 LU/dB and 0.1 LU LRA budgets. Reports include binary hashes.

Logs: `build/measurement-build.log`, `build/scaffold/Testing/Temporary/LastTest.log`,
`build/measurement-reference/itu-run.log`, `build/measurement-reference/ebu-run.log`.
M1's defined qualification now passes on both required platforms; final logs and
identities are in `Tests/Reference/results/m1-qualification-2026-10-17`.
No delivery-compliance claim is made.
Next: synchronized reference picture before polished UI.
ATSC A/85:2026-07 long-form dialogue measurement needs its own later work item.

The final qualification harness, portable packaging and generated UI fixtures are
documented in `Tests/Reference/README.md`. Future work is split into
[15_VideoImplementationHandover.md](15_VideoImplementationHandover.md) and
[16_Vst3MeasurementParityMilestone.md](16_Vst3MeasurementParityMilestone.md).
Video implementation has not started; the plugin remains unchanged/uninstalled.

## Objective and scope

Finish the existing standalone scaffold for a commercial, hardware-independent loudness-conforming application for film/TV post. The agreed sequence is **standalone scaffold → first measurement/logic milestone → synchronized reference video → full UI**, with video session/timing requirements defined during logic design. The application must work with ordinary audio devices; the Pro 40 is one test device, not a dependency.

## Starting point (14 October, before this continuation)

- Branch `main`, commit `556d387ad772dc5e604ead95f1525f60e66fd5f9`. Before this handover, only three untracked `.DS_Store` files were present; preserve unrelated work.
- CMake already defines `LevelScopeStandalone`, alongside the VST3 plugin and `LevelScopeCore`. It uses JUCE 7.0.12 and C++17; macOS defaults request Intel/Apple Silicon and deployment target 10.13. Compatibility still needs build/runtime verification.
- `Source/Standalone/` contains the app, workflow shell, session model, target-profile catalogue and readiness checks. Source selection stores file metadata and checks extensions; it does not prove the audio can be decoded.
- Source/Target/Current State navigation and Current State initialization exist. Measurement, proposal, audition/render and export are not wired. Target families are placeholders with explicitly non-authoritative numeric limits.
- The older project overview/build notes understate the existing scaffold. Inspect current source before deciding what to add. No LevelScope build or runtime test was performed during this handover.

## Read first

1. `Docs/ProjectIndex/00_ProjectOverview.md` and `11_StandaloneRequirements.md`.
2. `Docs/ArchitectureSnapshot.md` and `Docs/ProjectIndex/09_SignalFlowContract.md`.
3. `Docs/ProposalEngineSpec_v2.md` — proposed design; use it to align interfaces, not as evidence of implemented behavior.
4. `Docs/ComplianceMatrix.md`, then relevant parameter/module inventories as needed.
5. `CMakeLists.txt` and `Source/Standalone/`.

## Scaffold completion (15 October)

The existing shell was extended, not recreated. Completion criteria and implemented outcomes:

1. **Validated media:** cancellable background import probes contents as WAV/AIFF/FLAC, decodes every declared frame, rejects invalid rate/empty audio, non-finite samples and more than 12 channels, and records format, sample rate, bit depth, duration, frame count, tags and SHA-256. PCM short reads are rejected instead of accepting JUCE's zero padding. FLAC requires native STREAMINFO with a known frame count and PCM MD5; its decoded content is checked. Missing checksums, unknown lengths and lengths beyond JUCE 7's 32-bit FLAC reader limit are refused. This is audio validation, not comprehensive container/metadata conformance certification.
2. **Channel contract:** file-indexed roles remain separate from physical device indices. Decoder layout is a suggestion requiring explicit confirmation; unknown/duplicate/ambisonic roles are refused. No automatic rear/side swap or count-based routing. Routes start muted; each file channel maps to one of the first 12 device outputs or mute. Duplicate outputs and implicit summing/downmix are refused. Role changes clear measurements and Current State. Monitor-routing changes do not alter source/processing identity.
3. **Device/source audition:** explicit output selection/open/close, source Play from beginning, Stop, EOF, device/read failure handling, position and underrun count. Output-only; no input capture or system-default selection. Source/device rates must match; no resampler. Monitor gain defaults to −30 dB and is independent of render settings. A worker pre-fills a bounded SPSC buffer; the application audio callback uses only prebuilt routing, sample copies, gain smoothing and lock-free atomics. Device callback quiescence precedes teardown or configuration changes.
4. **Lifecycle:** one active import, progress, cooperative cancel, terminal success/failure/cancel states and monotonic session/job tickets. Successful replacement clears dependent state; failed/cancelled replacement retains the previous source. Clear invalidates a pending result. Close cancels and joins workers. Current State is a neutral, separate, unmeasured document, gated on confirmed roles.
5. **Future service boundaries:** `StandaloneServices.h` defines worker-only analysis/proposal/render/report APIs and value snapshots: source hash, file roles, half-open source-frame ranges, session ticket, versioned processing graph, engine/profile/build/platform identities, cancellation/progress and typed results/errors. Measured values remain optional and separate from predictions. Report verdict defaults to not evaluated. Export gates additionally require exact verification and approval for the current revision. No service implementations are registered.

Plugin DSP, parameter IDs, APVTS and state chunks were not edited. No VST3 target was built or installed. The app remains hardware-independent.

## Build and validation

Build **only named standalone/test targets**; the plugin still has automatic installation enabled. Fresh build:

```sh
cmake -S . -B build/scaffold -DLEVELSCOPE_BUILD_SCAFFOLD_TESTS=ON
cmake --build build/scaffold --target LevelScopeStandalone LevelScopeScaffoldTests -j 6
ctest --test-dir build/scaffold --output-on-failure
```

The validated build used the same pinned JUCE 7.0.12 checkout via `-DFETCHCONTENT_SOURCE_DIR_JUCE=/Users/fuzz/Developer/LevelScope/build/juce-7.0.12` (tag commit `4f43011b96eb0636104cb3e433894cda98243626`). Fresh fetches are now shallow. `.gitignore` excludes `/build/`; the original misnamed `gitignore` and unrelated `.DS_Store` files were left intact.

- **Baseline:** original commit `556d387` was extracted with `git archive HEAD` into `build/baseline-source`, configured in `build/baseline` using the pinned checkout, then built with `cmake --build build/baseline --target LevelScopeStandalone -j 6`. Passed after correcting the JUCE helper's deployment setting; log: `build/baseline-build.log`.
- **Toolchain fix:** JUCE 7 bootstraps `juceaide` with a nested CMake process. Its host-default 15.x minimum made `CGWindowListCreateImage` unavailable. Main CMake now passes the configured minimum through `MACOSX_DEPLOYMENT_TARGET` during helper configuration. For the already-created baseline helper cache, the equivalent explicit configuration was `cmake -S build/juce-7.0.12 -B build/baseline/_deps/juce-build/tools -DJUCE_BUILD_HELPER_TOOLS=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_OSX_DEPLOYMENT_TARGET=10.13`. Use a fresh build directory if an old helper cache still holds 15.x. JUCE itself was not patched or upgraded.
- **Updated app/tests:** universal Release build and CTest passed. Build log: `build/scaffold-build.log`; test log: `build/scaffold/Testing/Temporary/LastTest.log`. Remaining compiler/configure warnings are in JUCE (unused macOS variable and CMake CMP0175 compatibility).
- **Fixtures:** mono/stereo/6/12-channel PCM, AIFF/FLAC, rates 44.1/48 kHz, hash/source preservation, malformed and truncated media, empty/missing/13-channel sources, NaN samples and in-flight cancellation. Channel-role confirmation, reversed file order, route bounds/duplicates, ring-buffer wrap and concurrent producer/consumer, silence on underrun, stale job rejection, state invalidation, profile gates and frame-range validation are checked.
- **Fake device:** explicit 12-output open with zero inputs, rate rejection, prefill, playback/EOF, restart/stop during buffering and device loss. This does not qualify a physical device or OS audio backend.
- **Silent GUI smoke test:** launched on Sequoia; stereo fixture import, role confirmation, Current State initialization, replacement with a 12-channel fixture, cleared Current State, muted routing and access to file channels 10–12 verified. No audio device was opened and no playback was started.
- **User-reported High Sierra smoke test (15 October):** the user performed a short check and reported that the app works as it does on Sequoia. Record this as successful basic operation on High Sierra; detailed test scenarios, exact OS point release and audio-device coverage were not supplied. It is not full compatibility or listening qualification.

Artifact: `build/scaffold/LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app`.

## Platform distinction

- Earlier LevelScope development took place on **Mac OS X High Sierra**.
- Current development host is the user's **Mac Studio M1 on macOS Sequoia**. This run reports **macOS 15.7.9 (24G830), arm64, Xcode 16.4 / AppleClang 17, macOS SDK 15.5**.
- The configured minimum deployment target remains **10.13**, with `x86_64;arm64` builds. Mach-O inspection reports **10.13 for x86_64 and 11.0 for arm64** (the Apple Silicon platform floor). These load commands are not runtime compatibility evidence.
- The agent exercised the **arm64 Sequoia runtime**. The user subsequently reported a successful short **High Sierra / Intel runtime** check of this scaffold. Comprehensive High Sierra compatibility, Rosetta and Windows qualification remain outstanding.
- The separate experimental **Pro 40 driver was developed and tested on Sequoia**. The older High Sierra Mac served as a **vendor-driver reference system**, not the current driver development host.

## Repository checkpoints

Keep one shared codebase for High Sierra and Sequoia; do not copy changes into a separate implementation.

- Development branch: `codex/standalone-scaffold` in `carlozaiana/LevelScope`.
- Original High Sierra development baseline: annotated tag `high-sierra-pre-scaffold-2026-10-15` at `556d387ad772dc5e604ead95f1525f60e66fd5f9`. This preserves the pre-scaffold repository state; it is distinct from the new user-tested scaffold.
- Completed scaffold checkpoint: annotated tag `standalone-scaffold-2026-10-15`, including this handover and the user's High Sierra smoke-test report. Logic-module development must preserve this tag.
- `main` remains at its original baseline; publishing the development branch does not merge the scaffold into it. Build products and unrelated `.DS_Store` files are excluded from the scaffold commit.

## Video plan

Prepare the video contract during the next logic-module design; implement basic synchronized reference-picture playback after the first working measurement/logic milestone and before the polished UI.

- Attach an optional reference video to the audio session. Define its identity/path, duration, rational frame rate/timebase, starting timecode and audio/video offset. Handle frame/timecode conversion explicitly, including fractional rates and drop-frame display where applicable.
- Keep the authoritative audio source explicit. Initially, LevelScope's selected audio remains the processing/audition source and the video's embedded audio is muted; attaching a picture must not change measurement or routing silently.
- Add shared transport semantics for play, stop, seek, restart and end-of-media. Define audio-clock authority, video synchronization and playback-latency compensation; the current start-from-zero/stop audition is insufficient for synchronized picture.
- Keep video decoding and display outside the real-time audio callback and outside the measurement/proposal algorithms. Retain High Sierra/Sequoia compatibility when selecting APIs and codecs; evaluate JUCE's video component or a platform adapter with real sync tests before committing to a backend.
- First delivery: synchronized reference playback and offset adjustment with the existing audio workflow. Audio-track extraction from video, replacement-audio video export, broad codec support and polished video controls are later, separate features.
- Validate sync on both platforms using known picture/audio sync points, including seeks, offsets and long playback. Preserve the original media; do not claim frame accuracy until measured.

## Next logic-module chat

Read `StandaloneServices.h`, `StandaloneSessionModel.h`, `StandaloneJob.h`, this handover and the proposal spec. Implement exact measurement/analysis first, including shared full-program source/current/render verification semantics, explicit per-role detector/LFE policy, and reference vectors with agreed tolerances. Then connect deterministic proposal methods, authoritative versioned profiles and reproducible processing plans. Never promote predicted metrics to measured truth or use placeholder profile values for compliance.

The current interfaces are a scaffold contract, not a frozen on-disk session/report format. Define/validate the graph JSON schema, freeze model and report schema before implementing persistence or batch execution. Reorder file channels explicitly when adapting to Core's sorted `AudioChannelSet`; do not lose file-order semantics. Re-check source content identity before analysis/render and require same-revision verified output before export approval. Reject output paths/aliases that could overwrite source media.

Remaining limits: in-memory sessions only; no analysis/proposal/DSP-processed audition/render/export/report execution; no compliance claim, resampling, seek/pause, downmix, MP3/M4A, ambisonic/object workflow, video playback or polished UI. Audition detects source changes through size/mtime and decoder checks, not a continuous cryptographic identity lock; same-size/same-time external edits are not guaranteed detectable. Decoder cancellation is cooperative and shutdown may wait for the current filesystem read. The agent performed no physical device test or audible listening validation; the user's High Sierra smoke test does not establish comprehensive device or real-time qualification. Coordinate the listening setup with the user before any later audible test; software-only scaffold completion does not waive that requirement.

## Studio test context

The installed experimental Pro 40 device exposes **8 inputs / 12 outputs at 48 kHz** through CoreAudio. Nuendo lists all 12 outputs; only physical outputs 1–6 have listening validation: L, R, C, LFE, LS, RS. Outputs 7–12 are stream channels, not a claim of 12 analogue connectors. No three-hour qualification is claimed.

Keep speaker roles, file-channel order and hardware output numbers separate. Earlier QuickTime tests exposed a rear-versus-side WAV-label interoperability issue; explicitly tagged playback worked. Do not globally swap channels to accommodate one application's convention.

Leave the installed Pro 40 driver, login warm-up, keep-active service and recovery helpers untouched. Obsolete test packages were archived separately without restarting the working driver or resetting the background-item database.

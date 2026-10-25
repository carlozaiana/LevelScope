# M7 combined candidate-03 — completion evidence

Status: **approved M7 scope qualified on both actual Macs; ready for checkpoint
publication and integration of the unchanged candidate**. See
[qualification-summary.json](qualification-summary.json). Full physical 7.1
listening remains separately pending. This directory is separate from the immutable
24 October candidate-01/02 and paint-investigation records. No M7 publication has
been performed.

The user lifted the freeze and accepted occasional picture skips under
[the revised decision](revised-acceptance.md). Original strict harness results
and exit codes are retained; the separate acceptance evaluator cannot turn an
unexplained timing, audio, identity, memory or lifecycle failure into a pass.

## Candidate identity

- Integrated main: `5e07bbd29930db34236d6dfcede7fea14983f90a` plus the preserved,
  uncommitted M7 changes on `codex/logic-module`.
- Source-manifest digest: `7abc7587514b7f0f750720ac7853abeaf55b21c90697b25de7740cec2df034bf`.
- Measurement build: `facfe596c7126188b24627dfeab06dec907332b96442634c9ad8ee29675ea069`.
- Package identity SHA-256: `b2cf915f2f6b644993782363b0a0c243c49b43198100ae39ff6be870eff35f9e`.
- Standalone executable SHA-256: `14ba35d2f6ddc0bded680df267e81996f768c1537be1d6d360da87376b6889db`.
- Ordinary video qualification executable SHA-256: `dbeb25405e6d2ee7fc8f54ad425675830656a29e8103529651df92432fa8f89d`.

This is a tested working-tree identity, not a checkpoint commit. A future commit
must preserve these exact implementation/build inputs and record the mapping.
Documentation and result additions do not change the compiled source manifest.
The optional diagnostic executable is an investigation aid, not a playback fix.
No diagnostic tracing is enabled in the normal application or timing harness.

## Evidence organization

- `build-record.json`, `candidate-identity.json`: exact build and package identities.
- `sequoia/`, `high-sierra/`: fresh M1–M7 native/runner/reference suites, diagnostic
  isolation and native H.264 cases. Both automated summaries passed; their explicit
  `qualification_complete: false` means the suite alone is not the full milestone.
- `m1-m5-cross-platform.json`: documented metric tolerances and three-hour
  diagnostic/section resource budgets; passed.
- `m6-cross-platform.json`, `m7-cross-platform.json`: all 124 and 218 captured
  streams respectively were bit-identical across actual arm64 and Intel execution.
- `sequoia-ui-smoke.json`, `high-sierra-ui-smoke.json`, `ui-reuse-review.json`:
  fresh integration checks plus explicitly reused unchanged UI observations.
- `integration-equivalence-review.json`: exact production input and FP-flag review.
- `long-programme-reuse.json`: prior actual three-hour eight-channel measurement
  and seek evidence retains its original binary identity.
- `source-matrix-reuse.json`: first sixteen Source ProRes cases per Mac, including
  Sequoia's original strict failures. Those failures remain in the raw records;
  the new acceptance decision is separately recorded.
- `sequoia-unchanged-behavior.json`: old/new short M7 measurements and captures agree.
- `sequoia-prores-native/`, `high-sierra-prores-native/`: all 53 native cases per
  Mac passed, including the large-file loading checks.
- `sequoia-video-short/`, `high-sierra-video-short/`: all 20 cases per Mac accepted.
- `sequoia-video-endurance/`, `high-sierra-video-endurance/`: all three 1,800-second
  Current ProRes cases per Mac accepted. No retries were performed.
- `fresh-physical-audit.json`: all 46 fresh physical reports independently accepted;
  zero skips, underruns or timing-bound failures, memory checks passed, audio
  configuration unchanged. The reused Source matrices retain their three historical
  skips and two original strict failures under separate revised acceptance.
- `local-preservation-audit.json`, `high-sierra-preservation-audit.json`,
  `remote-preservation-check.json`: frozen packages, historical records, working-copy
  boundaries and all ten pre-existing tag objects/targets preserved.

Only text evidence is intended for Git. Binaries, reference/generated media,
raw PCM captures, temporary verified copies and connection details remain outside
this evidence directory. Their identities are recorded, and frozen packages remain
preserved in the existing local and Golem test folders.

## Evidence limits

Offline/automated support does not require eight physical outputs. Full physical
7.1 routing/listening remains pending and is not implied by muted built-in-output
tests. Video reports establish internal audio-clock/frame-selection behavior,
not optical or acoustic synchronization. Occasional picture skips remain an
accepted development limitation; sustained stalls, drift and audio interruptions
do not. No plugin, installed application, system audio setting or Pro40 change
is part of this work.

## Reproduction of the combined build

From the repository root, use a new build directory and the pinned JUCE checkout.
Do not reuse or overwrite a preserved qualification package. This records the
additional diagnostic targets and all test options explicitly:

```sh
cmake -S . -B build/logic-m7-combined-repro \
  -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_OSX_ARCHITECTURES=x86_64;arm64' \
  -DCMAKE_OSX_DEPLOYMENT_TARGET=10.13 \
  -DFETCHCONTENT_SOURCE_DIR_JUCE=/Users/fuzz/Developer/LevelScope-dependencies/JUCE-7.0.12 \
  -DLEVELSCOPE_BUILD_MEASUREMENT_TESTS=ON \
  -DLEVELSCOPE_BUILD_SCAFFOLD_TESTS=ON \
  -DLEVELSCOPE_BUILD_COMPLIANCE_TESTS=ON \
  -DLEVELSCOPE_BUILD_PROPOSAL_TESTS=ON \
  -DLEVELSCOPE_BUILD_DIAGNOSTICS_TESTS=ON \
  -DLEVELSCOPE_BUILD_SECTIONS_TESTS=ON \
  -DLEVELSCOPE_BUILD_AUDITION_TESTS=ON \
  -DLEVELSCOPE_BUILD_MULTICHANNEL_TESTS=ON \
  -DLEVELSCOPE_BUILD_VIDEO_TESTS=ON \
  -DLEVELSCOPE_BUILD_VIDEO_PAINT_DIAGNOSTIC=ON
cmake --build build/logic-m7-combined-repro --config Release --parallel 4 --target \
  LevelScopeStandalone LevelScopeMeasurementTests LevelScopeScaffoldTests \
  LevelScopeComplianceTests LevelScopeProposalTests LevelScopeDiagnosticsTests \
  LevelScopeSectionsTests LevelScopeAuditionTests LevelScopeMultichannelTests \
  LevelScopeVideoTests LevelScopeVideoQualification \
  LevelScopeVideoPaintDiagnostic LevelScopeVideoPaintTraceTests
```

`Tests/package_m7_qualification.py` creates a new package from that build plus the
pinned M1 reference bundle and existing video/UI/M7/ProRes fixture directories.
Its required arguments are documented by `--help`; use the preserved candidate-03
`identity.json` hashes to verify fixture equivalence, never substitute downloaded
files merely because their names match. A rebuild receives its own binary hashes.

For the actual candidate-03, from its package root:

```sh
python3 -B Tests/run_m7_qualification.py --platform sequoia \
  --completion-regressions --output /absolute/new/sequoia-output
```

On actual High Sierra use `/usr/local/bin/python3` and `--platform high-sierra`.
`--completion-regressions` reruns all short/reference/runner checks and the M6
sparse long-audition test; it explicitly omits the real M7 three-hour test. That
omission is justified here only by `long-programme-reuse.json` and the exact
integration equivalence review. Omit the flag when that equivalence does not hold.
Keep generated captures in ignored output folders for cross-platform comparison.

`Tests/run_m7_completion_video.py` runs the fixed short or endurance programme from
the frozen package. It requires `--candidate`, `--fixtures`, `--actual`, `--vfr`,
`--output`, `--platform`, `--phase short|endurance`, and explicitly coordinated
`--allow-audio`. It pins the approved built-in device/rate by platform and mutes
all routes. The precise invocations and expected media/build identities are in
each physical block's `runs.json`. Use new output directories; there is no retry
loop or diagnostic tracing. A temporary `caffeinate -di` process surrounded each
block to prevent idle sleep without changing power settings.

After all four physical block directories are present, the post-run text audit is:

```sh
python3 -B Tests/Reference/results/m7-completion-2026-10-25/audit_fresh_physical.py
```

It verifies the tested source manifest, exact case inventory, raw-report hashes,
platform/device/media/build bindings, processing gain, duration, eight-channel
presence, original strict verdict and revised acceptance. It recomputes the
existing memory checks and writes `fresh-physical-audit.json`. This verifier is a
post-run evidence tool, not a compiled input to candidate-03.

# M2 — versioned profiles and measured programme rules

Implementation continuation, 18 October 2026, in the existing logic working copy
`/Users/fuzz/Developer/LevelScope` on `codex/logic-module`, based on main
`514b7bce25a34ccd02af28a6a24100a44f8189a4`. The video worktree and all checkpoint
tags remain fixed. This milestone was explicitly approved with the narrow rule
interpretation below.

## Qualification status

**M2 qualification passes on actual Intel High Sierra and arm64 Sequoia.**
Candidate-02 is the final tested implementation, saved in the annotated checkpoint
`m2-measured-compliance-2026-10-18` and integrated into main.
Both platforms passed all automated suites and the full silent functional UI
checklist. The full UI checklist ran on candidate-01; a subsequent correction to
three static import/role messages was checked in candidate-02 on both platforms,
and all automated suites were rerun. The evidence records this coverage explicitly.
Existing M1/video acceptance remains preserved.

Evidence: [M2 qualification directory](../../Tests/Reference/results/m2-qualification-2026-10-18).

## Published checkpoint and main handover

- [Annotated M2 checkpoint](https://github.com/carlozaiana/LevelScope/releases/tag/m2-measured-compliance-2026-10-18):
  `m2-measured-compliance-2026-10-18`.
- Final main revision **at this publication** is the same commit as
  `m2-measured-compliance-2026-10-18^{commit}`. This immutable reference identifies
  the complete source, tests, qualification evidence and this handover. The tag's
  annotation records the full checkpoint and final-main commit SHA; resolve it with
  `git rev-parse 'm2-measured-compliance-2026-10-18^{commit}'`.
- Main advances by fast-forward from `514b7bce25a34ccd02af28a6a24100a44f8189a4`.
  There are no integration fixes or implementation differences from candidate-02.
  The logic workspace remains `/Users/fuzz/Developer/LevelScope` on
  `codex/logic-module`; the video working copy is untouched.
- [Qualification report pinned to the checkpoint](https://github.com/carlozaiana/LevelScope/blob/m2-measured-compliance-2026-10-18/Tests/Reference/results/m2-qualification-2026-10-18/qualification-summary.json)
  and `publication-review.json` preserve the acceptance and source/artifact audit.
  The latter records every pre-existing tag object and peeled commit for the final
  remote-preservation check. The earlier `preservation-audit.json` describes the
  pre-publication state; its unchanged-main/no-commit fields are historical.

The source manifest was reverified against all 91 implementation/test/manifest
files, all 230 packaged artifacts and all 73 pre-publication evidence files before
saving. The two-platform comparator was rerun against the retained final reports
and passed. Existing actual-platform automated/UI evidence is reused because the
tested implementation is byte-identical; publication changes only documentation
and evidence metadata. The qualification identity retains its original working
candidate/base-commit description rather than pretending the binaries were rebuilt
after this Git commit. A future rebuild must record its own binary identity.

| Check | Sequoia arm64 | Actual High Sierra x86_64 |
| --- | --- | --- |
| M2 rule/identity/lifecycle checks | 128 passed | 128 passed |
| Measurement regressions | 184 passed | 184 passed |
| Scaffold regressions | 205 passed | 205 passed |
| Video timing/transport + seven native fixtures | 329,054 passed | 329,054 passed |
| EBU runner tests | 8 passed | 8 passed |
| Pinned ITU files | 39/39 passed | 39/39 passed |
| Pinned EBU normative files | 66/66; 72 rules, 1,052 readings | 66/66; 72 rules, 1,052 readings |
| EBU supplementary alignment | 4/4 passed | 4/4 passed |
| Silent interactive UI | passed; affected status check repeated on candidate-02 | passed; affected status check repeated on candidate-02 |

Both official suites also passed comparisons against the accepted M1 baseline.
The direct M2 cross-platform comparison has zero errors: EBU whole-file metrics
match exactly; ITU's largest difference is 1.7763568394002505e-15 dB in true peak.
EBU curves/checkpoints and all 128 M2 outcomes match; the profile digest is identical.

Exact runtimes: Sequoia 15.7.9 (24G830), Apple M1 Max/Mac13,1, arm64, Python 3.9.6;
High Sierra 10.13.6 (17G66), Intel i7-3770K/iMac14,2, x86_64, existing Python 3.8.10.
Platform guards, executable architecture reports, process headers, dependency
listings and SHA-256 artifact verification are saved. No runtime was installed.

An initial sandboxed Sequoia bundle attempt could not create JUCE temporary test
fixtures or read CPU identity. It is not qualifying evidence. The complete rerun
outside that sandbox passed with the same binaries and no implementation change.
The initial local CTest also passed all four targets. Postprocessing additionally
checked that the comparator rejects wrong-platform and out-of-tolerance reports.

The High Sierra UI run found obsolete “Not measured” wording in the import-status
description after a successful measurement. Candidate-02 removes measurement claims
from the three import/role descriptions; the separate measurement status remains
authoritative. Its WAV import, role confirmation and completed-measurement display
were verified on both platforms with devices closed, routes muted and playback
stopped. This changes no measurement, evaluation, worker or transport behavior.
Candidate-01 reports and full UI observations remain in the evidence directory's
`candidate-01/` subdirectory; final reports and the change review are at its root.

Screen Sharing initially ignored automated remote clicks. After confirming its
Control Screen arrow, opening and cancelling the local View menu and sending the
remote click within the same CUA call restored repeatable input delivery. Short
held clicks and observing text entry before committing edits were reliable. No
persistent Screen Sharing setting, audio configuration or installed app was changed.
The user performed the initial WAV import; the agent confirmed roles and completed
all subsequent High Sierra checklist actions. Both isolated test apps were closed
after qualification.

## Candidate identity

At build time, the tested working candidate was based on the main commit above,
with uncommitted M2 changes. Those exact source files are now saved by the M2 tag.
The base commit alone does **not** reproduce this implementation.
`identity.json` captures the exact source files (including new files), compiler,
configuration, reference manifests and delivered artifact hashes. Its SHA-256 is
`9060672020a66a238328beb429a9022bf1fc1ea9f0fb8b050722c8afdfb29ab6`.

- Candidate source manifest digest: `9efedef8a34cab28971a60d7723d16bb50dc763fc1b649a3c3e4ed8c7d9229ea`.
- Compiled source/build fingerprint: `321e6230d26129f41c62d9c0a22171be2e28ef4851dcbf3872d39fdc40a393ea`.
- Profile definition digest: `ef6f5cc071a8c9d5a347da0dd61fcc46e9f24da0aed787bc9f3ea48cb3db7e66`.
- App executable SHA-256: `af0b45770102e69aac195066d371086dc832c189f7be94938a5ae0455edba305`.

The same universal artifacts ran on both machines. x86_64 minimum deployment is
10.13; arm64 minimum is 11.0. The postprocessing comparator is included in the final
candidate's source manifest. Final evidence/docs do not change the tested
implementation. The saved checkpoint includes the matching source/tests and final
evidence; rebuilding requires recording new artifact hashes.

Local bundle: `build/m2-qualification-2026-10-18/candidate-02`.
Remote bundle: `~/LevelScope-M2-qualification-2026-10-18/candidate-02`.
Both are separate from the preserved M1/video candidates and installed applications.

## Profile definition and applicability

The catalogue contains one authoritative definition:
`ebu-r128-2023-programme-linear-qc`, definition version `1.0.0`, schema `1`.
Primary reference: [EBU R128, November 2023](https://tech.ebu.ch/docs/r/r128.pdf).
The downloaded reference SHA-256 is
`4292cd2396e4cb386c3a6c1412cb4370fc108ef50956e78a24fbd8657ce2ca28`.

| Rule | Comparison | Provenance |
| --- | --- | --- |
| Programme loudness | -23.2 <= measured LUFS <= -22.8; intended normalization target remains -23 LUFS | clauses h, i, k, l |
| True peak | measured dBTP <= -1 | clause m |
| LRA | Informative; no acceptance ceiling | clause n |

The QC allowance is exactly ±0.2 LU. Clause m separately gives a ±0.3 dB
measurement tolerance for signals band-limited to 20 kHz. This is documented in
the definition but is never added to the peak ceiling. This evaluator does not
certify that bandwidth condition. Engine conformance and cross-platform test
tolerances are separate; neither widens the profile's acceptance interval.

The user must confirm pre-produced, whole-programme linear-audio QC. Applicability
starts unconfirmed and is reset on source replacement/clear. Live exceptions,
deliberate lower targets, short-form supplements and distribution-specific rules
are unsupported. Context is an explicit user assertion, not inferred from file
duration or suffix. Metadata, distribution encodings and complete broadcaster
delivery requirements are not evaluated. ATSC and Custom remain placeholders;
there is no general profile editor.

## Architecture and identity

- `Core/Proposal/Model/TargetProfile`: copied value definitions, explicit method,
  applicability, edition, provenance, limits and validation. SHA-256 covers the
  canonical versioned definition, including documentary meter provisions.
- `Core/Proposal/VerificationEngine`: deterministic selected-rule evaluation,
  independent of transport, devices and processing. These reusable files are linked
  into standalone support; the plugin does not consume them in M2.
- `StandaloneCompliance`: adapter from accepted immutable `AnalysisSnapshot` to
  measured-rule input. Binds content, ordered roles, full range, domain, plan,
  engine/build/platform and session identity, plus the exact accepted job/snapshot
  and evaluated values. Predicted and rendered-output domains are refused.
- `StandaloneSessionModel`: stores immutable measurements; derives evaluations
  from the selected profile and current applicability. A profile change reuses
  compatible measurements and does not change the audio revision or start a job.
  Current gain edits clear Current State measurement/evaluation and retain Source.
  Role changes, source replacement and clear remove the corresponding truth.
  Worker ticket/revision/cancellation checks still reject stale completions.
- Video identity, stream/offset, routing and monitor gain remain separate from
  audio measurement and evaluation identity. No audio/video timing code changes.

The measurement kernel, sample traversal and numerical contract remain M1. Values
are compared as binary64 against inclusive bounds formed from the profile target
and allowance; there is no epsilon or rounding before comparison. Adjacent values
can display identically while producing different verdicts. UI values are labelled
with their Source or Current State domain and show rule bounds/deltas.

Invalid/placeholder profiles, unsupported context/method/domain and stale/missing
snapshot identities produce `notEvaluated`. An unavailable, non-finite or inconsistent
required metric produces an `indeterminate` rule. Any known rule violation makes the
aggregate `fail`, retaining the other indeterminate rule details. Otherwise any
indeterminate rule prevents a pass. With no LRA constraint, unavailable/unstable LRA
does not prevent passing the required rules. Future explicit LRA constraints must
have their own provenance and require stable available LRA.

Success reads “Passes selected programme loudness and true-peak rules.” It does not
approve a delivery or a future encoded export. The public proposal/export services
remain unavailable; the existing export identity primitive is not a verdict.

## Acceptance tests

`Tests/ComplianceTests.cpp` covers independently specified exact/adjacent boundaries,
the -0.7 dBTP tolerance trap, display-equivalent values, invalid definitions and
provenance, missing/non-finite metrics, optional LRA, deterministic bindings,
unsupported contexts/methods, prediction separation, real Source/Current State
measurements, compatible-profile reuse, gain invalidation, stale/cancelled results,
role/source replacement/clear and audio/video/routing independence. It emits a
per-check JSON report and actual executable architecture.

Qualification also requires the existing measurement/scaffold/video/native suites,
runner checks and both pinned official reference suites on actual High Sierra and
Sequoia. Compare reports using the existing 0.01 LU/dB and 0.1 LU LRA budgets.
These comparison budgets must not affect profile decisions. Silent UI checks must
keep devices closed, routing muted and playback stopped. No VST3 target or installed
application is built/updated by this qualification.

Fresh configuration uses `build/logic-m2`; do not reuse or overwrite preserved
measurement/video products. The historical pinned JUCE 7.0.12 source is consumed
read-only at `build/juce-7.0.12` (commit
`4f43011b96eb0636104cb3e433894cda98243626`); no dependency checkout is moved or cloned.

```sh
cmake -S . -B build/logic-m2 -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_OSX_ARCHITECTURES=x86_64;arm64' -DCMAKE_OSX_DEPLOYMENT_TARGET=10.13 \
  -DLEVELSCOPE_BUILD_COMPLIANCE_TESTS=ON -DLEVELSCOPE_BUILD_MEASUREMENT_TESTS=ON \
  -DLEVELSCOPE_BUILD_SCAFFOLD_TESTS=ON -DLEVELSCOPE_BUILD_VIDEO_TESTS=ON \
  -DFETCHCONTENT_SOURCE_DIR_JUCE=/Users/fuzz/Developer/LevelScope/build/juce-7.0.12
cmake --build build/logic-m2 --target LevelScopeStandalone LevelScopeComplianceTests LevelScopeMeasurementTests LevelScopeScaffoldTests LevelScopeVideoTests -j 6
ctest --test-dir build/logic-m2 --output-on-failure
```

Package into a fresh directory (audio and binaries stay ignored):

```sh
python3 -B Tests/package_m2_qualification.py \
  --output build/m2-qualification-2026-10-18/candidate-NEW \
  --build build/logic-m2 --m1-bundle build/m1-qualification-2026-10-17 \
  --video-fixtures build/video-implementation/probe \
  --ui-fixtures build/m1-ui-fixtures-2026-10-17-v2
```

From that bundle root, use the existing platform's Python runtime, with a fresh
results path every time. The Sequoia run needs ordinary user access to macOS
temporary storage and CPU identity; a restrictive automation sandbox is unsuitable.

```sh
python3 -B Tests/run_m2_qualification.py --platform sequoia --output ../sequoia-NEW
# On the actual older Mac, after transferring to a new isolated folder:
/usr/local/bin/python3 -B Tests/run_m2_qualification.py --platform high-sierra --output ../high-sierra-NEW
# Back in the repository, after collecting both result directories:
python3 -B Tests/compare_m2_reports.py \
  build/m2-qualification-2026-10-18/sequoia-03 \
  build/m2-qualification-2026-10-18/high-sierra-02 \
  --output build/m2-qualification-2026-10-18/cross-platform-comparison-candidate-02.json
```

## Silent interactive checklist

Use only the isolated candidate and `ui-fixtures/` supplied with its bundle.
Keep Device closed, all output routes Muted and playback stopped throughout.
An automated suite pass does not substitute for these observations.

1. Import `stereo.wav`; confirm Left/Right roles. Select EBU programme QC, leaving
   applicability unchecked. Measure Source: it must remain NOT EVALUATED until
   applicability is confirmed; then PASS with the narrow selected-rule wording.
2. Initialize Current State; measure at 0 dB: same I/TP and a separate Current State
   PASS. Change trim to -6 dB: Current State immediately loses its result; Source
   remains PASS. Remeasure: programme loudness FAIL, peak PASS.
3. Select None and ATSC: NOT EVALUATED. Return to EBU: the previous Source and
   Current State evaluations return without a measurement job. Applicability
   confirmation in step 1 must reuse the existing measurement.
4. Change monitor gain while muted. Attach `video-fixtures/cfr-exact.mov`, apply
   a +1 s offset, then close the picture window. Audio results/evaluations must
   remain intact. Embedded audio stays excluded.
5. Replace Source with `stereo.aiff`, then `stereo.flac`; confirm roles and measure
   each. Replacement resets applicability and Current State. Import and confirm
   `unsupported-12ch.wav`: measurement is disabled with the supported-layout reason.
6. Import `long-stereo.wav`, confirm roles and start measurement. Cancel while
   active: discard the completion and never show a pass. Clear Source: remove
   source/current results. Reimport the long source, start measurement and close
   the app while active. Verify process exit and clean device-closed relaunch.

Each platform's `ui-smoke.json` records the final qualification basis, the affected
candidate-02 recheck and local screenshot hashes. Full candidate-01 observations
are preserved under `candidate-01/`, including all 16 functional UI checks and the
High Sierra display finding resolved by candidate-02. Screenshots remain in ignored
local build folders; reports contain their hashes. No reference audio or binaries
belong in the source checkpoint.

## Remaining scope

No automatic gain proposal, scene/dynamics processing, freezes, general profile
editing, processed audition, public export/report, persistence or new video feature.
Embedded soundtrack playback remains deferred. VST3 measurement parity stays in
[its separate milestone](16_Vst3MeasurementParityMilestone.md). A later deterministic
global-trim proposal can use this shared evaluator after M2 qualification/acceptance.
No M2 qualification issue remains open. Subsequent implementation requires its own
bounded milestone approval; publishing M2 does not authorize proposal generation.

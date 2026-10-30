# M17 qualification and handover

30 October 2026. **Candidate-04: Sequoia tested; High Sierra, video and three-hour
checks pending.** The approved Sequoia phase is complete, with no outstanding
required check in that phase. **Published development checkpoint:**
[`m17-job-lra-presets-sequoia-2026-10-30`](https://github.com/carlozaiana/LevelScope/tree/m17-job-lra-presets-sequoia-2026-10-30),
commit `e422f85ccbf7ae11876e6a31be2409937cb8381c`.
[Publication receipt](../../Tests/Reference/results/m17-publication-2026-10-30/README.md)
records remote verification, preserved packages and all 20 prior tags. Main integration
is not authorized. Main remains qualified M13
`a9be070ab9eb5cb4f591ff9d991d8999c3a016d0`; work remains in the established
`codex/logic-module` working copy. Read the approved scope in [60](60_M17JobLraRequirementsAndPresets.md)
and the pre-implementation impact plan in [61](61_M17ImpactQualificationPlan.md).

## Implemented boundary

`JobRequirements` is a separate, versioned definition around the unchanged EBU
R128-2023 programme QC profile. Its default has no LRA ceiling. An explicit ceiling
is stored as integer tenths, 0–1000: **0.0 LU is a strict ceiling**, not unset.
Measured values retain their full precision, with inclusive comparison and no
added epsilon, engine-test tolerance, or delivery allowance. The base IL/TP target,
allowance, ceiling, applicability and catalogue digest remain unchanged:
`ef6f5cc071a8c9d5a347da0dd61fcc46e9f24da0aed787bc9f3ea48cb3db7e66`.

The added requirement carries its own ID, revision, name, method, applicability,
user provenance and definition digest. Evaluation records retain the base profile
digest, added-requirements digest and combined rules digest separately. No user
preset acquires standards authority. Whole-programme LRA uses the existing method;
short-programme instability, silence, absent/non-finite/negative metrics, stale
inputs and predictions cannot satisfy a required rule. With no LRA rule selected,
unavailable LRA does not prevent the existing IL/TP pass.

The inline panel edits a draft. **Use EBU requirements** explicitly confirms the
definition and selects the compatible base profile while stopped with no conflicting
job. Checking the option does not supply a suggested ceiling. Edits mark the session
unsaved and invalidate proposals and publication approvals, while unchanged audio,
processing identity and compatible measurements are retained and reevaluated.
An identical confirmation is a no-op. M16 preservation preferences remain separate.

Portable `.levelscope-target` files are bounded, flat, versioned UTF-8 JSON; they
contain requirements only. Loading is transactional. Saving the active definition
excludes unfinished drafts and uses an owned temporary file followed by replacement.
Cancellation, malformed/unsupported/tampered data and open/write/replace failures
preserve the active session and previous valid file. Sessions contain an independent
complete definition, without a live dependency on an external preset.

Session schema **7** preserves coherent revision and runtime job guards. Schemas
1–6 migrate in memory with LRA unset; opening does not rewrite the original. Explicit
Save requires an M17 reader; Save As retains an older copy. An independently damaged
requirements record recovers otherwise valid media/settings with requirements
unavailable and dependent approvals disabled until explicit replacement. Bad caches
remain replaceable and do not discard valid requirements. Structural core corruption
still fails. Evaluations are reconstructed from compatible validated measurements.

All selected rules passing produces a genuine automatic no-op, retaining gain,
envelope, measurements, evaluation and audio revision. Otherwise, **any selected
LRA requirement makes changing automatic proposals unavailable in M17**, even when
LRA itself passes. The explanation says that correction for these requirements is
not supported yet, not that no solution exists. No private candidate pass is started
for this unsupported path. Manual all-channel processing and independent Current
remeasurement remain available. Plain-EBU automatic correction remains unchanged.

Both Current and the finalized, reopened export must pass the complete selected
requirements before publication. Source need not pass. A requirements-only edit can
reuse the exact rendered artifact and its valid measurement, but must refresh the
evaluation/report/approval. An old approval cannot publish. Render/verify and Publish
remain separate; PCM24 quantization/dither and decoded artifact measurement are
unchanged. JSON report schema **v5** and self-contained HTML separate base rules,
additional user requirements, measured domains, unavailable statuses and provenance.
Measurement-only reports never invent an exported artifact hash or verdict.

## Exact tested identity

- Frozen app: `/Users/fuzz/Developer/LevelScope/build/m17-qualification-2026-10-30/candidate-04/LevelScope Standalone.app`.
- Source/test/build manifest SHA-256: `853a5711c162fbbbfa9daff2dbc1acee3fee18ee94445e2ac0ae3407efdeb904`.
- Measurement/build identity: `5d1c4dc69d6df19c0fcf3dd39b53cffdd34495d23b52bed716657cb6e2c507b8`.
- Universal app executable SHA-256: `b84db4eff7945d5e264aa7f67bab7be153aede1758b6d5818e6d2cbbab680126`.
- Base HEAD: `06f3f75bcf956787124079d24f8fbd60b46b2fa6`, following published M16
  `73a7b25e71f87299fa2d93c1b1106c118c55beec`. **Neither identifies the M17
  implementation.** The frozen manifest is its exact tested source authority.
  Publication binds those same bytes to the new annotated checkpoint without a
  rebuild: checkpoint `e422f85ccbf7ae11876e6a31be2409937cb8381c`, annotated tag
  object `d243e4a68503bf719fc01ab2eefb2156d982ab20`. The publication receipt
  verifies all 277 committed manifest files and 293 frozen artifacts.
  Documentation/qualification receipts are outside the source digest and finalized
  after freezing. Original package identities retain their historical freeze-time fields.
- Actual execution: arm64 Apple M1 Max, Sequoia **15.7.9 / 24G830**. Release
  arm64+x86_64; all ten standalone/test executables have Intel minimum macOS 10.13.
  Apple clang 17.0.0, existing JUCE 7.0.12 commit
  `4f43011b96eb0636104cb3e433894cda98243626`. No dependency/deployment changes.

All 277 current source/test/build files match the candidate manifest. All 293
candidate-04 artifacts were verified, alongside earlier M17 candidates and the 283
M16 candidate-03 artifacts. Freeze-time pending status is intentionally retained
inside each immutable package; [final-status.json](../../Tests/Reference/results/m17-qualification-2026-10-30/final-status.json)
records the completed Sequoia phase separately. Full [textual evidence index](../../Tests/Reference/results/m17-qualification-2026-10-30/README.md).

## Fresh Sequoia qualification

All **twelve explicitly selected commands** pass their expected outcomes in
**153.69 seconds** on candidate-04. No recursive historical runner, Golem access,
video execution, physical-device opening, full cinema or three-hour workflow ran.

| Required check | Result |
| --- | --- |
| M17 requirements/presets/numerical boundaries/identity/session/manual/export/UI/resources | PASS: eight required groups, 372 assertions |
| Short measurement service, including legacy stale-requirements approval guard | PASS: 185 checks |
| Existing measured compliance | PASS: 182 checks |
| M16 preferences/support/full inspection/identity | PASS: seven groups, 256 checks |
| Global proposals | PASS: 302 checks |
| Short session migration/recovery/coherence | PASS: 1,006 checks |
| M12 short processing/offline/pre-monitor/export/audio timing | PASS: 334 checks |
| M13 short proposal/Apply/independent verification | PASS: 366 checks |
| Marker import/mapping/session compatibility | PASS: 487 checks |
| Short PCM export/report/gates | PASS: 4,872 checks |
| Frozen M16 reader rejects valid schema 7 | PASS: expected exit 1; exact old build/reason, checked core checksum, schema, applied plan and caches |
| Selected runner self-tests | PASS: eight tests; missing/duplicate/failed/mismatched evidence cannot qualify |

New tests include real 65-second mono/stereo/5.1 Source/Current/finalized exports,
under/exactly/over-60-second stability boundaries, exact/adjacent full-precision LRA
thresholds, unavailable/non-finite/negative/stale/predicted inputs, strict zero versus
unset, malformed/oversized/deep/duplicate/tampered presets, transactional failures,
independent session definitions, schemas 1–6 migration, late completions, manual
processing and cancellation, stale proposals/approvals, and artifact-only LRA failure
through the production publication predicate. Synthetic injected boundary failures
are distinguished from actual rendered/decoded artifact measurements.

Native candidate-04 checks cover default unset and readable invalid-value feedback;
safe measured-session reopen; genuinely unchanged confirmation with the correct
active notice and saved/unchanged state; already-pass no-op; preset save, malformed
load, then identical valid reload clearing the old error; strict zero invalidation;
SaveAs/New/reopen with exact gain/measurement/provenance. Moving only the disposable
external preset before reopening proves session independence. Current visibly passes
the selected rules after reopening, without a new measurement. The app was closed
stopped/device-closed/muted.

Impact-selected native reuse is explicit: candidate-01 supplies invalid-requirements
and damaged-cache recovery/preset cancellation observations; candidate-03 supplies
manual edit → unavailable Current → independent measured failure, plain-EBU explicit
Apply/verification, actual rendered artifact reevaluation and local publication,
measurement-only reports, byte-hash checks and Safari visual review. Their defects
remain recorded. Candidate-03 to candidate-04 changes only confirmation/load **display
refresh** in `StandaloneJobRequirements.cpp` and its test. All backend paths for these
reused UI workflows have identical source; their current-build automated tests are
also fresh on candidate-04. Affected native confirmation/load paths are rerun on
candidate-04. Do not relabel earlier UI execution as candidate-04 execution.

### Resource observations

| Observation | Measured | Approved bound |
| --- | --- | --- |
| Cached Source/Current reevaluation p95 | 0.2654 ms | 10 ms |
| Inline panel feedback p95 | 0.0201 ms | 100 ms |
| Short session verify/reopen | 118.92 ms; 66,814 bytes | 5 seconds; existing session cap |
| Native SaveAs / reopen including UI observation | 528 / 555 ms | 5 seconds |
| Candidate-04 disposable automated tree | 390.41 MiB allocated | 2 GiB planning budget |
| Frozen universal package | 169.42 MiB allocated | 1–3 GiB planning estimate |

Preset input is capped at 64 KiB, names at 80 Unicode characters and provenance at
2 KiB UTF-8. Retained definition payload is bounded and checked below 1 MiB; one
active value plus bounded request/document/chooser snapshots replaces any need for
an unbounded preset catalogue. This is a bounded-data check, **not a measured process
RSS delta**. Existing session/report caps remain. No new audio copy or private pass
is required solely by a requirements change. These are short workloads, not
three-hour or High Sierra resource evidence.

## Reproduction and preservation

Build only the named standalone/test targets listed in `package_m17_qualification.py`:
`LevelScopeStandalone`, `LevelScopeJobRequirementsTests`, `LevelScopeComplianceTests`,
`LevelScopeMeasurementTests`, `LevelScopeProposalTests`, `LevelScopeSectionProposalTests`,
`LevelScopeSessionTests`, `LevelScopeSectionProcessingTests`, `LevelScopeMarkerTests`
and `LevelScopeExportTests`. The established universal Release configuration enables
the corresponding test options, including `LEVELSCOPE_BUILD_JOB_REQUIREMENTS_TESTS`.
Configuration/compiler/dependency receipts are indexed with the package. No VST3
target is built or installed. Freeze into a **new** directory, never over a candidate.

```sh
python3 -B build/m17-qualification-2026-10-30/candidate-04/Tests/run_m17_qualification.py --candidate build/m17-qualification-2026-10-30/candidate-04 --output build/m17-qualification-2026-10-30/sequoia-reproduction-new
```

The runner requires actual arm64 macOS and verifies artifact hashes before execution.
Its output directory must not exist. Tests need access to their owned temporary
verified-source copies; no physical-device approval is implied. Native UI, preservation
and cross-build receipts remain separate required evidence.

At completion of qualification, all 20 existing tag objects/targets matched preflight
locally and directly on GitHub. Main, logic HEAD and the registered video worktree
location/HEAD were unchanged; no publication actions had yet occurred. Subsequent
authorized development-checkpoint publication is recorded separately. No plugin/app
installs, Pro40 access or system audio changes were performed. Unrelated Finder metadata is preserved. Only synthetic disposable media
and test exports were generated; raw media/builds remain outside Git evidence.

## Historical failures and reuse boundaries

Candidate-01's explicit runner stopped at a marker migration assertion expecting
“M16 reader”; the implemented schema-7 warning correctly said “M17 reader”. Candidate-02
changes only that test assertion and its binary. Its completion record reuses eight
passing candidate-01 commands with exact application/test binary equivalence, runs the
remaining four commands, and validates all twelve required records. The original
failed run remains a failure. Native review then found an incorrectly encoded dash
in the invalid-LRA warning. A further candidate changes that warning to plain ASCII
“0.0 to 100.0”. Candidate-03 passes the selected short manifest. Further native
inspection finds that confirming an unchanged definition can retain an obsolete
draft notice. Candidate-04 refreshes the confirmed display (also on identical preset
load) and adds a no-op feedback regression; the selected short manifest and affected
native UI are rerun. The earlier observations remain attributed to their executed
binaries.

The source audit compares M16 candidate-03 (manifest
`3ac2e19d8dd976a8fdc41490245ad9961af91e4bfe1a1d5afdc1ad5c03897f38`)
with the implementation. The detailed audit, rather than this abbreviated narrative,
is authoritative for file identities. Measurement arithmetic, programme decoding,
composition, processing, audition buffering, PCM containers and video remain unchanged.
All 24 Core/Measurement and Core/Processing files also match qualified M13. The older
measurement service gains one approval-identity guard; a fresh short service test
covers it. Shared evaluation/session/proposal/report/export callers are tested fresh.

Earlier reference/decoder/sample-parity/long-media evidence remains attributed to
its original candidate, not rerun or claimed as fresh M17 evidence. In particular,
M13 candidate-02 (`50eba19a2d66096607ec878c8eaa05e610580c7f5e8854d7b375d846a89c5f27`),
M14 candidate-02 (`9dd22a598281e83431039f210202a8a24fe49517294cc6bb8d83032c4975e013`),
and M15 candidate-05 (`38676a45d47ddba6eb54f7ab299f52f622f6be26554256781b8960749c609564`)
retain their original reference, bounded RF64/mono and DAW provenance. These do not
qualify new M17 features on deferred platforms or full-duration workloads.

## Deferred coverage and next-work handover

Actual Intel High Sierra/Golem, cross-Mac exchange, video loading/decoding/playback
and three-hour workflows remain **pending**. Use one frozen combined candidate after
M18 and no later than M19, with impact-selected requirements. Local Intel 10.13 build
success is not runtime qualification. No current dependency or failure warrants
lifting these deferrals. Coordinate physical audio separately. Occasional picture
skips remain an accepted development limitation when video testing resumes; audio
interruptions, accumulating drift, sustained freezes and identity/lifecycle failures
remain blockers. Physical full-7.1 listening is separately pending.

M17 adds no LRA-directed correction, new analysis arithmetic, waveform, progressive
analysis, dynamics, region suggestions, dialogue/intelligibility claims, 7.1.4 or VST3
parity. Keep the next timeline milestone useful and bounded, followed promptly by
region suggestions and selective LRA correction under separately defined methods
and acceptance. Avoid requiring a comprehensive timeline editor or DSP suite first.

The future timeline/session design must restore **historical explanations behind
applied proposals**, their relevant measured baseline and requirement/policy/method
provenance alongside exact processing settings. Bind the explanation to its applied
plan and distinguish subsequently edited processing. Do not restore private candidate
measurements or old verdicts as current verification. This is documented preparation
only; M17 does not implement complete explanation persistence. Existing saved applied
policy/binding and exact processing remain; save a report to retain full transient
inspection until that later milestone is implemented.

M18 is now [proposed for review](63_M18InteractiveTimelineAndProgressiveAnalysis.md):
a bounded interactive audio timeline, progressive waveform/M/S previews and historical
applied explanations. No M18 implementation has started. Its publication follow-up
is documentation-only; main and the tested M17 source remain unchanged.

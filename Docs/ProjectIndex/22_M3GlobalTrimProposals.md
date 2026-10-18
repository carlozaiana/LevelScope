# M3 — deterministic global trim, explicit application and measured verification

Approved implementation continuation, 18 October 2026. Working copy:
`/Users/fuzz/Developer/LevelScope`, branch `codex/logic-module`, based on published
M2/main `01a34833ec06a0d237e098faadb23b398581e8f9`. The existing video worktree and
all checkpoint tags are preserved. This document supersedes the earlier handover's
“proposal unavailable” description for this bounded standalone workflow only.

## Qualification status

**M3 acceptance criteria pass on both actual platforms, 18 October 2026.** The
qualified implementation is saved as `m3-global-trim-2026-10-18` on
`codex/logic-module` and integrated into main by fast-forward from the published
M2 commit above. There are no integration fixes or changes to the tested code.

The [qualification report](../../Tests/Reference/results/m3-qualification-2026-10-18/README.md)
contains the complete source manifest, artifact hashes, OS/runtime/architecture
records, test logs, UI observations and preservation audit.

| Check, per platform | Actual arm64 Sequoia | Actual Intel High Sierra |
| --- | --- | --- |
| M3 proposal/session/worker assertions | 302 passed | 302 passed |
| M2 / measurement / scaffold assertions | 128 / 184 / 205 passed | 128 / 184 / 205 passed |
| Video timing/transport and native media | 329,054 checks + 7 fixtures passed | 329,054 checks + 7 fixtures passed |
| EBU / M3 runner tests | 8 / 7 passed | 8 / 7 passed |
| Official ITU reference suite | 39 files passed | 39 files passed |
| Official EBU v05 | 66 normative files, 72 rules, 1,052 readings + 4 supplementary files passed | Same, passed |
| Silent M3 interactive workflow | Passed | Passed |

Sequoia is macOS 15.7.9 (24G830), Apple M1 Max, executing arm64. Golem is
Mac OS X 10.13.6 (17G66), Intel Core i7-3770K, executing x86_64. Its existing
Python 3.8.10/runtime dependencies were inspected and reused without installation.
Direct comparison passes: M2/M3 results and all eight scalar proposal bindings
match exactly; EBU differences are zero; maximum ITU true-peak difference is
1.7763568394002505e-15 dB, with other metrics identical. Final CTest also passes 5/5.

Final candidate-04 source digest:
`5f66aafc81be73b8184a88f662f1494c9059dfc25cefe43e8e42aac7629f38fd`.
Build fingerprint:
`556349482888d4a86ebea83d1f15403045475a67d787fe422f54abaffdbb88ef`.
Universal standalone executable SHA-256:
`4a88b2b3c70e2887633a3202018f804754377657d0d2374b910df1239d8868f5`.
The source manifest matches the checkpoint implementation; publication changes
are documentation/evidence only.

UI review corrected three misleading status messages: no-change advice no longer
suggests a new measurement; Current State status no longer caches a profile verdict;
and proposal job history now uses neutral completion wording after context/gain
edits. Candidate-02 also exposed an obsolete test assertion against the removed
cached text; candidate-03 checks the actual measured verdict and snapshot domain.
The rejected run and exact candidate-03-to-04 diff are retained. Full functional
UI paths ran on candidate-03; the final wording and affected draft/apply/verify,
no-change and invalidation paths were rechecked on candidate-04 on both systems.
All automated suites ran again on candidate-04. Unchanged format/unsupported-layout
UI coverage is explicitly inherited from qualified M2, not reported as a new run.

The gate-crossing fixture demonstrates why predictions are insufficient: a
-23 LUFS prediction measures -23.854939 LUFS after Apply. Both UIs show actual
FAIL, retain the -3.975769027009 dB applied gain and stop without correction.
Cancelled verification likewise retains gain and displays UNVERIFIED. No M3
acceptance check remains outstanding. Devices stayed closed, routes Muted and
playback stopped; both candidate apps were closed after qualification.

## Published checkpoint and continuation

- [Annotated M3 checkpoint](https://github.com/carlozaiana/LevelScope/tree/m3-global-trim-2026-10-18):
  `m3-global-trim-2026-10-18`.
- The checkpoint and final main revision **at this publication** are the same
  commit, `m3-global-trim-2026-10-18^{commit}`. The tag annotation records its full
  SHA. Resolve it with `git rev-parse 'm3-global-trim-2026-10-18^{commit}'`.
- [Qualification report pinned to the checkpoint](https://github.com/carlozaiana/LevelScope/blob/m3-global-trim-2026-10-18/Tests/Reference/results/m3-qualification-2026-10-18/README.md)
  and [publication review](../../Tests/Reference/results/m3-qualification-2026-10-18/publication-review.json)
  record the evidence audit, acceptance clarifications and five prior annotated
  tag objects/targets. Earlier preservation audits describe the pre-publication
  state; their unchanged-main fields remain historical.
- The latest fetched main was `01a34833ec06a0d237e098faadb23b398581e8f9`, already
  containing accepted measurement, compliance and reference-video work. The
  fast-forward adds M3 without altering that history. Continue in the established
  `/Users/fuzz/Developer/LevelScope` logic workspace after inspecting/fetching the
  then-current main; the video working copy remains untouched.

Publication reverified all 100 source/test/manifest files, 245 packaged artifacts
and 52 pre-publication evidence files. The two-platform comparator was rerun and
passed. Existing actual-platform automated/UI results remain valid because the
implementation is byte-identical to final candidate-04. Its original base-commit
and working-candidate identity is retained: the saved binaries were not rebuilt
after the checkpoint commit. Any future rebuild must record its own binary hashes.

## Scope and decision policy

The only method is `global-trim-m3.1`, acting on the existing `static-gain-v1`
processing plan. It requires valid measured Source and Current State, explicit
programme applicability, the authoritative versioned profile and available finite
integrated-loudness/true-peak metrics. It rejects predictions and unsupported methods.
Source may fail the selected rules; its measurement must still be valid.

The profile supplies the intended target, QC allowance and true-peak ceiling. For
EBU R128-2023 programme QC v1.0.0 these remain -23 LUFS, inclusive -23.2 to -22.8
LUFS, and strict -1 dBTP. Measurement-engine test tolerances are never added. LRA
is informative and is not steered. Profiles with an LRA constraint are unsupported
by this method. The profile's existing whole-programme linear-audio applicability
and distribution/live/short-form exclusions remain unchanged.

With measured Current State loudness `I`, peak `TP` and existing total gain `g`:

```
desired correction = target - I
lower = max(-60 - g, target - allowance - I)
upper = min(24 - g, target + allowance - I, peak ceiling - TP)
```

An already passing Current State returns **NO CHANGE**, retaining its gain,
measurement, evaluation and audio revision, without another measurement. Otherwise
an empty intersection returns **NO PREDICTED SOLUTION**, with no applicable plan.
The chosen correction is the nearest feasible point to the intended target.
The resulting total gain is `g + correction`, within the existing -60 to +24 dB
processing range. It affects every file channel, including LFE.

Arithmetic and application retain binary64 precision; displayed decimals never
feed back into the plan. Up to 16 adjacent representable-gain adjustments resolve
floating-point boundary rounding using forward strict profile comparisons. They
do not widen a rule or implement a measured-audio correction/retry loop. If none
is valid, the method reports no representable candidate.

Integrated loudness and true peak are scalar **predictions**, never measured passes.
Absolute-gate membership and float sample rounding can change the actual outcome.
No LRA prediction is supplied. A draft explains this uncertainty and its constraint.

## Architecture and identity

- `Core/Proposal/Methods/RuleBasedMethod`: deterministic value-only calculation,
  explanation and SHA-256 proposal binding. It consumes the M2 profile/evaluator.
- `StandaloneProposal`: adapts immutable measured snapshots. Their content hash,
  ordered roles, sample rate, channel count and full programme range must match
  each other, even when both snapshots separately validate against their own
  contexts. Current State must also match the actual current processing plan.
- `StandaloneSessionModel`: validates the exact current binding again when a
  worker result is accepted and when Apply is invoked. It reconstructs the plan
  from those inputs rather than trusting editable proposal output fields. The
  binding includes exact measured snapshot/evaluation identities and values,
  profile definition/version, method and full-precision current gain.
- `StandaloneMainComponent`: drafts asynchronously from copied values; job tickets,
  cancellation and revision checks discard stale completion. Draft/Apply/Discard
  are explicit controls. No-change and no-solution results cannot enable Apply.

Draft and Discard leave audio and measurements intact. Apply changes only the
Current State static gain, invalidates its measurement/evaluation and marks it
**UNVERIFIED**. Source remains intact. The user then invokes Re-measure Current
State through the qualified offline worker. M2 evaluates the actual new snapshot.
No prediction can enable a measured verdict, another proposal from unmeasured
audio, or export.

Cancelled/unsuccessful verification leaves the applied gain in place and Current
State explicitly unverified. Completed measurement with a rule violation displays
the actual measured FAIL. Neither path rolls back, retries, or applies another
correction. Detected external source changes invalidate source/current truth and
require explicit reimport; they do not silently reset the applied gain.

Profile/applicability edits discard a draft and reevaluate compatible existing
measurements without changing audio identity. Processing/role/source edits reject
stale drafts and invalidate the affected measurements. Video attachment/offset,
monitor gain and routing preserve the audio measurement and proposal identities.
Selected external audio remains authoritative; embedded video audio stays muted.
Source audition continues to play original audio, not the proposed/applied trim.

## Tests and reproducibility

`Tests/ProposalTests.cpp` covers exact and adjacent thresholds, deterministic plans,
gain limits, constrained/no-solution decisions, already-pass no-op preservation,
profile limits, missing/non-finite/stale/unsupported inputs, prediction exclusion,
optional LRA, individually valid but mismatched pairs, current-plan mismatch,
profile-versus-processing invalidation and video/monitor/routing independence.
Actual worker tests cover explicit application/verification, cancellation,
missing/modified source failures, an absolute-gate change that defeats prediction,
and gain on mono/stereo/5.1 (including LFE) at 44.1/48/96/192 kHz.

Use a fresh build directory and the pinned read-only JUCE 7.0.12 source (commit
`4f43011b96eb0636104cb3e433894cda98243626`). Never build the default/all target:
the VST3 target installs its output. Only these targets are required:

```sh
cmake -S . -B build/logic-m3 -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_OSX_ARCHITECTURES=x86_64;arm64' -DCMAKE_OSX_DEPLOYMENT_TARGET=10.13 \
  -DLEVELSCOPE_BUILD_PROPOSAL_TESTS=ON -DLEVELSCOPE_BUILD_COMPLIANCE_TESTS=ON \
  -DLEVELSCOPE_BUILD_MEASUREMENT_TESTS=ON -DLEVELSCOPE_BUILD_SCAFFOLD_TESTS=ON \
  -DLEVELSCOPE_BUILD_VIDEO_TESTS=ON \
  -DFETCHCONTENT_SOURCE_DIR_JUCE=/Users/fuzz/Developer/LevelScope/build/juce-7.0.12
cmake --build build/logic-m3 --target LevelScopeStandalone LevelScopeProposalTests LevelScopeComplianceTests LevelScopeMeasurementTests LevelScopeScaffoldTests LevelScopeVideoTests -j 6
ctest --test-dir build/logic-m3 --output-on-failure
```

`package_m3_qualification.py` creates a fresh candidate containing source/artifact
hashes, universal binaries, pinned ITU/EBU vectors, accepted M1 baseline reports,
seven video fixtures and silent UI fixtures. Generated media/builds stay ignored.
`run_m3_qualification.py` verifies artifact hashes, actual OS and architecture,
then runs M3/M2/measurement/scaffold/video/runner and both official reference suites.
Run from the bundle with `--platform sequoia` or `--platform high-sierra` and a
fresh `--output ../PLATFORM-NEW` directory. High Sierra uses its existing Python
3.8.10; it requires no runtime or app installation.

`compare_m3_reports.py LEFT RIGHT --output REPORT.json` requires both actual
platforms, identical candidate identity, M2/M3 outcomes and eight scalar proposal
plans/identities. It also compares reference results with the existing 0.01 LU/dB
and 0.1 LU LRA budgets, without changing profile decisions. Seven runner unit tests
include a one-ULP plan mismatch, missing outcomes and failed/wrong-architecture
reports. A working candidate is identified by its complete source manifest, not
by pretending its base commit already contains the uncommitted implementation.

## Required silent UI checks

Run the isolated candidate on each actual OS, with Device closed, routes Muted and
playback stopped. Do not install the app or touch preserved candidates.

1. Import `trim-required.wav`, confirm roles, select EBU and confirm applicability.
   Measure Source and Current at 0 dB. Both fail loudness; Draft becomes available.
2. Draft: inspect explanation, exact proposed gain and prediction labels. Current
   gain/measurement stays unchanged. Discard and redraft produce the same plan.
3. Apply: gain changes only now; old Current verdict disappears and UNVERIFIED is
   visible. Remeasure: actual Current passes; Source remains unchanged.
4. Draft the passing result: NO CHANGE, disabled Apply, retained gain/measurements
   and revision, and no new measurement requested.
5. Remove/restore applicability: draft is invalidated, measurements reused, no
   stale secondary verdict. Change processing: old Current result/draft disappears.
6. With a valid draft, change monitor gain while muted; attach reference picture
   and adjust offset. Measurements and draft remain valid, audio selection intact.
7. Import `gate-crossing.wav`; measure Source and Current at -50 dB. Draft/Apply
   predicts passing scalar values. Actual remeasurement fails loudness; applied
   gain remains, actual FAIL is shown, and processing stops without correction.
8. Use `long-stereo.wav` for verification cancellation after a gain Apply. Keep the
   applied gain and show UNVERIFIED. Replace/clear Source and check stale results
   cannot return. Close during measurement; confirm exit and clean silent relaunch.

## Remaining limitations and handover

No dynamics, scene processing, LRA steering, freezes, processed audition, public
export, persistence, embedded-soundtrack playback or VST3 parity. The separate
[VST3 parity milestone](16_Vst3MeasurementParityMilestone.md) remains unchanged.
No real-time or audible qualification is implied by silent offline checks.
Reference-picture timing/transport and its previously documented physical-sync
limits remain unchanged. Filesystem-read cancellation remains cooperative.

An actual measured selected-rule PASS covers only the named snapshot/profile rules,
not complete broadcaster delivery, metadata or a future encoded export. Future
changes must retain source-to-binary-to-result linkage; these qualification
artifacts must not be relabelled as a rebuild of a later commit.

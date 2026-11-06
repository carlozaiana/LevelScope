# M20 Logic interface and UI integration handover

6 November 2026. **Logic-stage review; M20 combined qualification is not complete.**
This handover covers the approved core/service/session work in
[70](70_M20RegionFormationAndPredictiveGainPlanning.md) and the pre-edit
[impact plan](71_M20LogicContractAndImpactPlan.md). UI controls, overlays and native
interaction qualification follow in the UI chat. No M20 push, tag or main merge is
authorized. Published main remains `22961f03de7919e1de631f8b13eb4b852108ac0b`.

## Starting point and ownership

The implementation stays in `/Users/fuzz/Developer/LevelScope`, branch
`codex/m20-region-planning`, following design handover
`b4ccb66b70ea4dcad227d106925fc04db493881e` and the pre-edit contract commit
`37d6449`. The final local implementation/receipt commits are recorded below after
checks complete. Integrate those commits into the existing UI working copy using
an appropriate `codex/` branch, preserving its pending work. Do not duplicate or
move either checkout, use a folder inside `build` as a worktree, or touch the video
working copy. UI owns controls/overlays; Logic remains the reviewer for constraints,
measurement, processing, ranking and session contracts. Strategy owns substantive
rule changes. No competing UI controls were implemented here; the only existing UI
edit is the necessary schema-9 migration/Save As wording.

## Frozen local integration package

- Package: `/Users/fuzz/Developer/LevelScope/build/m20-logic-2026-11-06/logic-handover-03`.
- App: the package's `LevelScope Standalone.app`. It contains the current ordinary UI;
  the new experimental method is available through the Logic interface, with no new
  automatic opt-in or advanced UI controls yet.
- Source manifest SHA-256:
  `fd200ebdc2d1e1f43715bafcba4d2f8e8a166d28283176fcc7c408f990ca6249`.
- Product build identity:
  `8c2397c14e757cd9f3c8580cd89c0a2bb065bc2e3de5a9eb97e30f03eb8f5254`.
- `identity.json` records exact source/test files and binary hashes. The source tree
  was not yet committed when packaged; its recorded base is provenance, not a claim
  that the base commit alone contains the implementation. The final commit must
  match the manifest. Preserve this package rather than replacing it with a rebuild.

## UI-facing immutable contract

`Source/Core/Proposal/RegionPlanning.h` defines value types; no UI/transport owns
or mutates a worker's inputs. `Source/Standalone/StandaloneRegionProposal.h` exposes
read-only attempt summaries. Use the existing proposal job/cancellation/completion
path in `StandaloneSectionProposal` and `StandaloneSessionModel`.

```cpp
using namespace levelscope::proposal;
auto settings = RegionPlanningSettings::defaults(
    RegionMode::automatic,
    sectionInputIdentity(session.measurementContext(MeasurementDomain::source)));
// Edit a private UI draft. Explicit confirmation only:
juce::String error;
if (!session.setRegionPlanning(settings, error)) showError(error);
// nullopt restores the unchanged M19/LRA-unset dispatcher.
// Copy request on the UI thread; execute generateSectionProposal in the existing job.
auto request = session.sectionProposalRequest();
```

- `RegionPlanningSettings`: explicit mode/selected ID or temporary source-frame
  span, support-edge ranges/pins, internal nominal-cut ranges/pins, margins,
  protections, split/minimum/separation/need/persistence controls, global permission
  and comparison intent. Time controls use integer milliseconds; explicit positions
  use exact source frames. Serialized coordinates are canonical decimal strings.
  Defaults/ranges are document 70's. Automatic global permission defaults on;
  selected modes default off. Unknown mandatory versions/protection fail closed.
- `RegionConstraints`: frozen roots H, permitted areas U, protection P and effective
  edge ranges. `RegionGeometry.support` means requested full changed-support edges;
  `RegionAssignment.nominal` and `geometry.cuts` mean processor plateau boundaries.
  Draw these as different concepts. Off-grid exact pins are never rounded.
- `RegionPlan.activity`: actual old/new-ramp support, affected/longest/margin frames,
  splits per original root, connected support regions, engine intervals/transitions,
  added/cumulative attenuation, absolute local dB-frame integral and global delta.
  None is an audio-quality score. Use programme frames for duration percentages.
- `SectionProposalResult::regionReview()`: immutable constraints, discovery
  alternatives/omissions, up to three actual complete attempts, scalar forecast and
  measured outcomes, forecast errors and resource/pass counts. `retainedAttempt`
  is out of range when no complete measured attempt is eligible. Do not index it
  without checking. No private candidate window populations are persisted as Current.
- `previewPlan()` remains inspection only. Only the existing `canApplySectionProposal`
  and explicit stopped-only `applySectionProposal` authorize processing changes.
  Apply revalidates full geometry and measured rules, then invalidates Current;
  independent Current measurement is still mandatory. Export independently measures
  the actual finalized PCM and retains its existing strict publication gates.

Requirements/preferences/region settings reuse compatible measurements but stale
proposals and export approvals. Automatic-mode binding excludes unrelated annotation
names and splits. A selected section's ID/bounds participate; changing them stales
its proposal. Audio identity includes the exact ordered input and processing plan.
Source/roles/replacement clears incompatible selection; job/runtime identity rejects
late completion after New/Open. Routing, monitor volume and video metadata do not
change proposal audio identity. Settings changes are stopped-only and blocked while
a job/persistence load is active. No state mutation occurs during planning.

## Search and prediction boundaries

Explicit M20 selection is required. The method validates matching measured
Source/Current, stable LRA, stereo/conventional 5.1 and the selected IL/TP/LRA rules
plus any M/S limits. All selected rules already pass: genuine no-op, zero private
passes. A sufficient permitted global-only correction is assessed first.

The auxiliary collector performs one streaming K-filter energy pass on the same
verified source copy, retaining 100 ms bins, channel maxima including LFE and labelled
EOF drain. It does not change the qualified meter. The approximate model anchors
window energies to the baseline, integrates the actual squared piecewise gain ratio,
and recalculates complete gate populations. Unknown TP/sampled M/S cannot confer
acceptance or exclude all plans. EOF filter memory, uniform bin energy and float
staging are explicit limitations. Actual minus predicted values are empirical errors,
not bounds or measurement tolerances.

The finite method explores at most 8 frozen roots, 12 nonmandatory cut nominations,
256 topology/edge constructions per root, 16 retained geometries per root, beam width
4, 64 distinct complete proxy plans (32/16/16), and 3 full candidate measurements.
A protected or unrepresentable envelope is rejected, never widened. Short-gap filling
is separately enabled/assessed and cannot cross protection. See document 71 for the
failed early broad-root case and the documented paired-gap discovery refinement.

Only an actual all-selected-rule pass enables Apply. Failed attempts use an anchored
whole-set minimum-LRA-excess +0.01 LU diagnostic shortlist, then other failed rules,
normalized excess and preservation costs. The 0.01 LU band never changes acceptance.
Failed bounded search means no acceptable candidate was found, not no solution exists.
Keep all measured attempts and restrictions accessible in Review.

## Sessions and historical explanations

Schema 9 embeds an independent `regionPlanning` definition (or explicit null).
Schema 8 and older migrate without selecting M20 or changing audio. Explicit saving
requires an M20 reader; Save As retains the original older-readable copy. Compatible
caches reconstruct evaluation; incompatible/damaged independent caches remain
replaceable. Existing stopped/device-closed/muted reopening and runtime-job isolation
are retained.

Applied explanation attachment version 2 stores M20 settings and the complete
historical narrative/scalar outcomes with before/after identities and requirement/
policy provenance. Version 1 remains readable and is still written for other methods.
Do not label old private metrics as current verification. Existing history status
identifies subsequently edited processing and changed requirements/preferences.
Ordinary global Apply remains independent of inactive experimental region settings.

## Initial UI integration checklist

1. Explicit experimental method selection with existing M19 as the unchanged default.
2. Private controls draft, validating stopped-only confirmation and readable errors.
3. Separate annotation, measured-root, allowed/protected-area, preview envelope and
   applied envelope layers. Numeric exact support edges versus nominal cuts/pins;
   optional margins and whole-programme global permission clearly labelled.
4. Independent counts for root splits, connected regions, affected duration, nominal
   intervals, actual transitions and omitted search work. Show reserve/deficit and
   requested versus achieved attenuation without rounding verdicts.
5. Source/Current inspection remains separate from audition. Lane selection inspects;
   ruler clicks use the stopped seek service; never implicit playback/device opening.
6. Predicted values and uncertainties beside actual attempt outcomes; all actual
   failures inspectable; no Apply from a proxy, stale result or unavailable metric.
7. Explicit Apply -> UNVERIFIED Current -> independent measurement. Cancellation or
   failure leaves applied settings and an honest unverified/measured-failure state.
8. Save/reopen exact settings and historical reasons; no restored private authority.

## Build and test entry points

Use only explicitly named standalone/test targets; no VST3 build/install. The local
build uses JUCE 7.0.12 at `/Users/fuzz/Developer/LevelScope-dependencies/JUCE-7.0.12`:

```sh
cmake -S . -B build/m20-ui-integration/compile \
  -DCMAKE_BUILD_TYPE=Release '-DCMAKE_OSX_ARCHITECTURES=arm64;x86_64' \
  -DCMAKE_OSX_DEPLOYMENT_TARGET=10.13 \
  -DFETCHCONTENT_SOURCE_DIR_JUCE=/Users/fuzz/Developer/LevelScope-dependencies/JUCE-7.0.12 \
  -DLEVELSCOPE_BUILD_REGION_PLANNING_TESTS=ON \
  -DLEVELSCOPE_BUILD_LRA_PROPOSAL_TESTS=ON \
  -DLEVELSCOPE_BUILD_JOB_REQUIREMENTS_TESTS=ON \
  -DLEVELSCOPE_BUILD_SECTION_PROPOSAL_TESTS=ON
cmake --build build/m20-ui-integration/compile --target \
  LevelScopeStandalone LevelScopeRegionPlanningTests LevelScopeLraProposalTests \
  LevelScopeJobRequirementsTests LevelScopeSectionProposalTests -j 6
```

`LevelScopeRegionPlanningTests` requires an explicit `--group`: `core`, `stereo`,
`surround`, `heterogeneous`, `failures`, `audio`, `population`, `gating`, `cinema`.
Use a fresh `--root` and `--output` per run. The last two require `--input` to the
existing bounded fixture, not the full film. Population testing models three-hour-sized
arrays without decoding three hours of audio. `Tests/Qualification/M20/check_logic.py`
requires every named fresh receipt, matching build and source-manifest files; missing
checks cannot produce Logic readiness, and it never marks combined M20 qualified.
`package_logic.py` freezes a new uniquely named package without modifying earlier ones.

## Qualification status and remaining work

Fresh receipts, explicitly reused evidence, resource observations and preservation
verification are recorded below. Historical failures remain historical. The unchanged
M19/M18 meter, processor, decoder and video evidence is dependency evidence only, not
fresh M20 qualification. No physical device, Golem or video execution occurs in this
Logic stage. Intel 10.13 cross-build success is not High Sierra runtime qualification.

After UI integration, freeze one combined candidate. Complete actual Sequoia/High
Sierra numerical/native/cancellation checks, cross-platform/session tests, native p95
<=100 ms responsiveness, cached planning targets <=10/20 s (30/60 s watchdog) and
additional model/search storage <=64 MiB. Evaluate the selected control variants and
freeze defaults before independent dialogue/ambience and music/effects holdouts;
coordinate listening separately. Preserve accepted occasional picture skips while
retaining audio continuity, authoritative-clock/drift and lifecycle blockers. Broaden
video/long-media testing only for an affected dependency or observed failure.


## Completed Logic-stage evidence

All 17 explicitly required groups pass (2,144 checks), with missing-receipt rejection
also exercised. See [textual evidence](../../Tests/Reference/results/m20-logic-2026-11-06/README.md),
[summary](../../Tests/Reference/results/m20-logic-2026-11-06/logic-check-summary.json)
and per-receipt provenance. Host: Mac Studio Mac13,1, Apple M1 Max, 32 GiB,
arm64 macOS 15.7.9 (24G830). Universal app contains x86_64/arm64; Intel load command
specifies macOS 10.13. No Intel execution or native M20 UI qualification is claimed.

| Frozen test | Outcome |
| --- | --- |
| Stereo / 5.1 complete workflow | LRA 5.910633862481 -> 4.400000001495 LU / 5.910633850473 -> 4.290033958071 LU; one complete candidate assessment each; exact independent Current metrics/peaks, measured finalized PCM24 passes |
| Broad root, two high parts and satisfactory middle | LRA 6.020599904473 -> 4.760704005700 LU against 4.9 LU; two splits and two disjoint changed supports; independent Current/artifact pass. Same bounded split-off search finds no passing plan. |
| Gate-crossing example | Complete retained measured LRA 6.999999884852 LU against 6.7: FAIL, no Apply. Separate model-only comparison predicts 7.0 -> 25.989999165642 LU with 2,971 previously excluded windows entering; no fourth candidate assessment. |
| Eight-minute cinema excerpt | Original excerpt SHA-256 `cfd57bea9a108db850b4e371de0f299f3396fbc1730e58e357380042a3064dd5` retained; LRA 23.435332665255 -> 22.385113195488 LU against 15: FAIL, no Apply; retained IL -23.143781961340 LUFS; three complete attempts, 64 proxies, Current unchanged. Not a full-film or Nuendo comparison. |
| New envelope audio | Stereo/5.1 at 44.1/48/96/192 kHz plus separate-file stereo and grouped-component 5.1 composition: exact offline/injected pre-monitor parity, all channels including LFE, chunk/seek/transition/EOF checks, no underruns |
| Permissions/session/failures | Exact pins, full old/new support, prior ramps/short radii, protections, bounds, splits versus components, additional/cumulative caps; cancel/source replacement/stale runtime; schema 8/9, invalid mandatory definitions and historical explanation; ordinary global Apply preserved |
| Existing methods | Explicit M19 core/stereo/surround/inputs/failures/gating; requirements/export/report/session checks; bounded section-proposal smoke. No recursive historical runners. |
| Three-hour-sized synthetic populations | 64 model plans in 1.727 s; analysis+sidecar allocation 35,846,608 bytes; whole-process peak footprint 44,008,896 bytes and max RSS 98,926,592 bytes (different OS accounting). This is cached model testing, not a three-hour media workflow or native UI responsiveness result. |

Production sources/app are byte-identical across local packages 01, 02 and 03.
Package 02 corrected the gate-test expectation: search need not select a known bad
gate-changing plan, so the model counterexample is now tested directly. Package 03
adds prior-envelope/short-radius boundary coverage. Only `RegionPlanningTests.cpp`
and its executable differ; core is fresh from 03, gating fresh from 02, other 15
checks retain their original package-01 receipts. Failed development and original
gate-assertion receipts remain visible, not rewritten as passes. The aggregate guard
checks the accepted manifest and product build identity and reports combined M20
qualification as false.

Observed worker planning allocations for these short/real workflows are about
1.4–3.6 MiB; cached planning 0.007–0.071 s, separate from source preparation,
auxiliary PCM analysis and complete candidate measurements. The owned build,
packages and fixtures occupy approximately 2.2 GiB. No frozen earlier package,
installed app/plugin, original media, Pro40 or system audio configuration was changed.
All 25 pre-existing tag objects/targets and main remain unchanged.

**Ready for UI integration of the Logic interface. Not ready for M20 checkpoint
publication or main integration.** Native UI, actual High Sierra/cross-platform,
selected control variants and independent real-material/listening review remain
required after UI integration.

# M18 dedicated UI/Timeline handover

30 October 2026. **Documentation and handover only; M18 implementation has not been
authorized or started.** No build, qualification, device opening or Golem access was
performed for this handover. This document supplements [63 — M18 proposal](63_M18InteractiveTimelineAndProgressiveAnalysis.md)
and incorporates the user's ruler-seek clarification. The new UI chat must obtain
implementation approval; preparing or copying this handover is not that approval.

## Ownership, exact starting point and preservation

The dedicated UI/Timeline chat will own M18 implementation, its candidate, qualification
and integration preparation. The logic-module chat stops overlapping implementation
and retains subsequent region/LRA proposal work. It can review measurement, processing,
identity and session interface changes. Propose small interface diffs for review before
shared contracts change; do not have both chats editing the same files concurrently.
The separate Video Implementation copy remains untouched. Main integration/publication
needs subsequent explicit approval after the combined candidate is qualified.

| Identity / location | Required value |
| --- | --- |
| Tested M17 checkpoint | `e422f85ccbf7ae11876e6a31be2409937cb8381c` |
| Annotated tag | `m17-job-lra-presets-sequoia-2026-10-30` |
| Tag object | `d243e4a68503bf719fc01ab2eefb2156d982ab20` |
| Candidate-04 source/test/build manifest SHA-256 | `853a5711c162fbbbfa9daff2dbc1acee3fee18ee94445e2ac0ae3407efdeb904` |
| Measurement/build identity | `5d1c4dc69d6df19c0fcf3dd39b53cffdd34495d23b52bed716657cb6e2c507b8` |
| Universal standalone executable SHA-256 | `b84db4eff7945d5e264aa7f67bab7be153aede1758b6d5818e6d2cbbab680126` |
| Original frozen app; read-only reference | `/Users/fuzz/Developer/LevelScope/build/m17-qualification-2026-10-30/candidate-04/LevelScope Standalone.app` |
| Main, to remain unchanged | M13 `a9be070ab9eb5cb4f591ff9d991d8999c3a016d0` |
| Existing logic copy | `/Users/fuzz/Developer/LevelScope`, `codex/logic-module` |
| Existing video copy | `/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation`, `codex/video-implementation` |
| Recommended NEW M18 copy | `/Users/fuzz/Developer/LevelScope-worktrees/UITimeline` |
| Recommended NEW branch | `codex/ui-timeline-m18` |

The proposed new folder and branch were absent at handover inspection. They have
**not** been created. Use a registered Git worktree, not a copied checkout and not
anything inside `build`. Main is older than M14–M17 and is the wrong checkout base.
Preserve all **21** existing checkpoint tag objects/targets, frozen packages, original
media, installations, plugin parameters/state/processing, Pro40 and system audio.

[62 — M17 qualification/handover](62_M17QualificationAndHandover.md) and the
[publication receipt](../../Tests/Reference/results/m17-publication-2026-10-30/README.md)
separate fresh/reused Sequoia evidence from pending coverage. M17's twelve selected
commands, native UI and Intel 10.13 cross-build passed; this is not High Sierra
qualification of M14–M17. Preserve the immutable package's historical freeze fields.
The later `fc4e4b0d54456cb32ab68672ede0a772b663d749` commit contains publication
documentation and ProjectIndex 63 only. This handover's own documentation-only commit
is reported separately in the originating chat; neither documentation commit is a
new tested binary or source-manifest identity.

### Recommended setup, for the new chat after setup is authorized

Inspect existing worktrees/branches first. If either proposed name already exists,
inspect it and stop on a mismatch; do not remove, reset or reuse it blindly. Read-only
fetch is reasonable to discover newer work, but do not replace the pinned M17 base
with `origin/main` or a moving logic branch. Report unexpected main/tag differences.
The following setup intentionally starts at the exact M17 commit, then adds only the
two documentation commits. It does not move the original logic checkout.

```sh
M18_LOGIC_ROOT=/Users/fuzz/Developer/LevelScope
M18_UI_ROOT=/Users/fuzz/Developer/LevelScope-worktrees/UITimeline
M18_TESTED_BASE=e422f85ccbf7ae11876e6a31be2409937cb8381c
M18_PROPOSAL_DOCS=fc4e4b0d54456cb32ab68672ede0a772b663d749
# Resolve the handover commit from this file's history, not a moving branch tip.
M18_HANDOVER_DOCS=$(git -C "$M18_LOGIC_ROOT" log -1 --format=%H -- Docs/ProjectIndex/64_M18UITimelineHandover.md)

git -C "$M18_LOGIC_ROOT" worktree list --porcelain
git -C "$M18_LOGIC_ROOT" branch --list codex/ui-timeline-m18
git -C "$M18_LOGIC_ROOT" rev-parse main 'm17-job-lra-presets-sequoia-2026-10-30^{commit}'
git -C "$M18_LOGIC_ROOT" show --stat "$M18_HANDOVER_DOCS"
git -C "$M18_LOGIC_ROOT" diff --name-only "$M18_TESTED_BASE" "$M18_HANDOVER_DOCS"
# Review: only Docs/ and the M17 publication-receipt directory may differ.
# Verify the handover hash against the originating chat's final response.

git -C "$M18_LOGIC_ROOT" worktree add -b codex/ui-timeline-m18 "$M18_UI_ROOT" "$M18_TESTED_BASE"
git -C "$M18_UI_ROOT" cherry-pick --ff "$M18_PROPOSAL_DOCS"
git -C "$M18_UI_ROOT" cherry-pick --ff "$M18_HANDOVER_DOCS"
git -C "$M18_UI_ROOT" diff "$M18_TESTED_BASE" -- Source CMakeLists.txt
git -C "$M18_UI_ROOT" status --short
```

The source/CMake diff must be empty; verify all 277 candidate-manifest entries for
an exact starting-source check. Do not use blanket branch merges to import future
logic changes. Keep future builds, synthetic media and candidate packages in the
new copy's own `build/`; preserve old build trees without reusing their CMake caches.
During M18, integrate required shared fixes deliberately and requalify affected code.
Return the frozen candidate/manifest and reviewed changes for integration; do not move
main or publish a checkpoint merely because a development phase passes.

## Read-first document and code map

Read 63 for complete acceptance/budgets; 14 for the measurement/identity contract;
23–24 for diagnostics/sections; 28/30 for verified audition; 35/36 for sessions;
39/40 for export/report gates; 44/47 for section processing/proposals; 53/55 for marker
mapping; 56/59 for preferences; 58 for the LRA/timeline roadmap; 60–62 for requirements
and M17 evidence. Read 18/27 for reference-video timing/ProRes and the optional
`Docs/VideoPaintDiagnosticMainIntegration-2026-10-24.md` investigation-tool boundary.
`Docs/ProposalEngineSpec_v2.md` remains the broader proposal architecture, not authority
to implement excluded strategies in M18. Later accepted milestone contracts take
precedence over obsolete historical status paragraphs in older documents.

- **Measurement / display / audition:** [14 — OfflineMeasurementContract](14_OfflineMeasurementContract.md), [23 — M4MeasuredDiagnostics](23_M4MeasuredDiagnostics.md), [24 — M5ProgrammeSections](24_M5ProgrammeSections.md), [28 — M6VerifiedCurrentStateAudition](28_M6VerifiedCurrentStateAudition.md), [30 — M6QualificationAndHandover](30_M6QualificationAndHandover.md).
- **Processing / proposals / markers:** [44 — M12SectionTrimProcessing](44_M12SectionTrimProcessing.md), [47 — M13BoundedSectionTrimProposals](47_M13BoundedSectionTrimProposals.md), [53 — M15MarkerImportAndSourceTimeMapping](53_M15MarkerImportAndSourceTimeMapping.md), [55 — M15QualificationAndHandover](55_M15QualificationAndHandover.md), [56 — M16ProposalPreferencesAndLimits](56_M16ProposalPreferencesAndLimits.md), [59 — M16QualificationAndHandover](59_M16QualificationAndHandover.md).
- **Sessions / reports / requirements:** [35 — PostM7BroadcastWorkflowProposal](35_PostM7BroadcastWorkflowProposal.md), [36 — M8SessionsQualificationAndHandover](36_M8SessionsQualificationAndHandover.md), [39 — M10VerifiedPcmExportAndReporting](39_M10VerifiedPcmExportAndReporting.md), [40 — M10QualificationAndHandover](40_M10QualificationAndHandover.md), [60 — M17JobLraRequirementsAndPresets](60_M17JobLraRequirementsAndPresets.md), [61 — M17ImpactQualificationPlan](61_M17ImpactQualificationPlan.md), [62 — M17QualificationAndHandover](62_M17QualificationAndHandover.md).
- **Roadmap / video:** [18 — VideoQualificationResults](18_VideoQualificationResults.md), [27 — ProResMainIntegration](27_ProResMainIntegration.md), [58 — LraWorkflowAndInteractiveTimelineDesign](58_LraWorkflowAndInteractiveTimelineDesign.md).

Paths below are relative to the new repository root; keep the existing structure.

| Entry point | Contract / expected M18 use |
| --- | --- |
| `Source/Standalone/MeasurementTimeline.{h,cpp}` | Existing measured M/S view, range/highlight/inspection, gate-aware display buckets; foundation, not a second engine |
| `Source/Core/Proposal/AnalysisEngine.{h,cpp}`, `Model/ProgramAnalysis.h` | `ProgramAnalysis`, `DiagnosticWindow`, `PeakLocation`, `prepareTimeline`, `inspectWindow`; full-resolution authority and EOF/gate semantics |
| `Source/Standalone/StandaloneDiagnostics{,Panel}.{h,cpp}` | Immutable `MeasuredDiagnostics`, profile-independent descriptive analysis and separately bound annotations; cancellation/generation guards |
| `Source/Standalone/StandaloneSections{,Panel}.{h,cpp}` | Existing section analysis/list, timeline linking and numeric gain draft/Apply controls; selection is not processing |
| `Source/Core/Proposal/Model/ProgrammeSections.h`, `Source/Core/Proposal/MarkerImport.{h,cpp}`; standalone marker adapters/panel | Annotation identity, exact frame boundaries, marker mapping/preview/confirmation; do not import plugin time offsets |
| `Source/Standalone/StandaloneMainComponent.{h,cpp}` | `beginMeasurement`, async completion, UI windows and session controls; natural controller integration point |
| `Source/Standalone/StandaloneServices.h`, `StandaloneJob.h` | Immutable worker `ServiceContext`, `AnalysisSnapshot`, domains, `JobTicket`, cancellation/progress and result ownership |
| `Source/Standalone/StandaloneMeasurementService.{h,cpp}` | Qualified worker sample loop, verification and final snapshot; optional preview observation must not alter results |
| `Source/Standalone/ProgrammeSource.{h,cpp}` (`ProgrammeReader`), `VerifiedSource.{h,cpp}`, `ProgrammeComposition.*` | Actual interleaved/assembled/composed input and verified-copy lifecycle; group validation remains all-or-nothing |
| `Source/Standalone/StandaloneProcessing.h`, `Source/Core/Processing/SectionGainPlan.h` | Shared `ProgrammeGainProcessor`, frame-based existing transitions, global gain, pre-monitor sample observation |
| `Source/Standalone/StandaloneSessionModel.h` | Freshness, document/runtime IDs, draft/application/invalidation, requirement reuse, proposal adoption and historical fields |
| `Source/Standalone/StandaloneAudition.{h,cpp}`, `StandaloneAuditionRequest.h`, `StandaloneAudioPanel.*`, `AudioPresentationClock.h` | Audition domain/identity, `synchronize`, stopped `seek`, `positionFrames`, device/routing/monitor separation |
| `Source/Standalone/SessionDocument.{h,cpp}`, `StandaloneSessionUI.cpp` | Schema 7, coherent `capture/save/read/verify/adopt/saved`, safe cache recovery; proposed historical-record schema migration |
| `Source/Standalone/StandaloneExport.*`, `StandaloneReport.cpp`, `StandaloneExportUI.cpp` | Current/artifact measured gates, report provenance/escaping; historical explanations are not export authority |
| `Source/Standalone/StandaloneJobRequirements.*`, `StandaloneCompliance.*`; `Core/Proposal/Model/JobRequirements.*` | M17 independent job requirements/presets and immediate reevaluation; unchanged EBU profile |
| `Source/Standalone/StandaloneVideoPanel.cpp`, `StandaloneVideo.*`, `VideoTiming.h` | `video::Time`, `floorTicks`, half-open picture coverage; video time is audio presentation time minus signed picture start offset |
| `CMakeLists.txt`, `Tests/package_m17_qualification.py`, `Tests/run_m17_qualification.py` | Named-target build and frozen source/artifact conventions; M17 runner is a selection example, not an automatic M18 runner |

## Concrete plugin reuse assessment

Inspected at the exact M17 implementation: `Source/VolumeHistoryComponent.h`, `.cpp`,
`_ViewNav.cpp`, `_Render.cpp`, `_History.cpp`, `Core/LevelScopeHistoryModel.{h,cpp}` and
`PluginEditor.{h,cpp}` (`MissionControlComponent`). The class is **not drop-in reusable**:
its constructor requires `LevelScopeAudioProcessor`, restores processor UI state,
registers APVTS threshold listeners, drains the processor FIFO and starts a timer.
The history model stores a DAW-indexed 60 Hz ring, live-derived loudness/gates and
plugin state, rather than immutable full-programme offline measurements.

| Candidate and concrete symbols | Decision for M18 |
| --- | --- |
| `_Render.cpp`: JUCE Path accumulation/stroking, clipped fills, grid and text drawing | **Direct presentation fragments:** reuse small drawing operations with explicit geometry/data supplied. No complete member can simply be linked without the processor-bound class. Keep extraction local and attributed; do not refactor the whole plugin. |
| `_Render.cpp`: `getTickStepSecondsWithHysteresis`, integer tick loop in `paint`; `_ViewNav.cpp`: `getTimeRulerArea` | **Adapt:** reuse tick spacing/hysteresis and distinct ruler hit zone. Replace 60 Hz indices with exact audio frames/rational labels, use common plot width, handle zero/small sizes and long programmes. Never inherit host/user timecode offsets or the hardcoded “Cubase common preroll” menu. |
| `_ViewNav.cpp`: `applyHorizontalZoom`, `panTime`, fit/reset | **Adapt:** cursor-anchored zoom and clamped scrolling are useful. Clamp against `[0,programmeFrames]`, not the furthest-written live ring. Preview coverage is distinct from duration. Preserve cursor frame across zoom, including the last frame. |
| `setFollowRightEdge`, auto-refollow in `paint` | **Adapt only:** use authoritative standalone playhead, not newest analysis/FIFO data. Recommend explicit Follow with manual pan disabling it and no surprise auto-refollow. Keep state updates out of paint. Analysis progress never drags the playhead. |
| `_History.cpp`: `recomputeGroupAbsFromChildren`, `selectBestLevelForCurrentZoom`, `buildVisibleGroupsForLevel` | **Adapt the reduction idea:** multiresolution extrema, stable grid anchoring and viewport budgets. Retain standalone included/excluded/ungated buckets and missing-data gaps. The plugin caps at up to 8,192 points; M18 retains its proposed 4,096-column bound. Do not import live ring overwrites or its float LUFS values as measurement authority. |
| `_Render.cpp`: `computeRepresentativeCurves`, min/max bands and representative polylines | **Presentation reference, not an accuracy guarantee:** representatives choose a centre/trend value within the range; bands may be hidden or suppressed by a threshold, and adjacent emitted points can bridge missing data. Keep both extrema visible, break gaps, preserve gate categories and use full-resolution inspection. Loudness bands are not PCM waveforms. |
| `_Render.cpp`: `markStaticBackgroundDirty`, `rebuildStaticBackgroundIfNeeded` | **Adapt:** JUCE cached backgrounds/scratch-path reuse. Cache keys need size, pixel scale, lane geometry/visibility, units, view transform and theme. Keep playhead/selection/live preview out of static caches; invalidate on every relevant change. |
| `_ViewNav.cpp`: rolling-lane divider/right-strip resizing and `setShowMomentaryCurve` / `setShowShortTermCurve`; `PluginEditor` visibility controls | **Adapt UI mechanics:** bounded splitters, cursor feedback and explicit lane toggles. Store only standalone presentation state; no processor persistence calls. Do not reuse plain lane-click curve toggling: lane clicks inspect/select. |
| `applyPersistedUiStateFromProcessor`, `syncPersistedUiStateToProcessor`, `LevelScopeHistoryModel::saveState/loadState` | **Do not reuse:** APVTS/DAW/plugin chunk persistence is distinct from `.levelscope` documents and runtime job IDs. M18 history is a bounded standalone session record. |
| `_History.cpp`: FIFO/bootstrap, running meters, `pushSecondSample`, `recomputeRollingLraForSecond`, `rebuildRollingLraStep` | **Deferred:** live 60 Hz/1 Hz sampling, rolling LRA and decay rules do not qualify standalone programme analysis. No rolling-LRA or cumulative-IL estimator in M18. |
| Threshold overlays, `handleThresholdMouse*`, `applyThresholdValuesDuringDrag`, `MissionControlComponent` processor/parameter controls | **Deferred:** automation gestures, DSP/latency/routing and plugin targets are not the M12 section draft model. Presentation ideas may inform later work only. |

Use `MeasurementTimeline` and qualified offline populations first. Extract a few pure
helpers only when they simplify M18; leave the plugin caller unchanged where feasible.
Do not add `PluginProcessor.cpp`, plugin `JuceHeader.h`, APVTS or `LevelScope` plugin
linkage to the standalone UI target. A comprehensive shared editor and VST3 measurement
parity are explicitly not prerequisites. Any genuinely shared helper edit needs impact
review, but must not become an excuse to build/install the VST3 target.

## Data, callbacks and authority boundaries

Use a thin presenter/controller between the passive view and session/transport.
Proposed callback names are design sketches, not existing APIs:

| View input / action | Required behavior |
| --- | --- |
| Immutable timeline binding | Runtime session ID + job/generation + source/member hashes + ordered roles + range/rate/frame count + domain + applied-plan key + data version |
| Final `ProgramAnalysis` / `AnalysisSnapshot` | Read-only validated final data; full-resolution inspection/rankings, final gates and metrics |
| Progressive batch | Bounded immutable coverage/sample/window batch; `ANALYSING`, gate membership pending, not an `AnalysisSnapshot` accepted for verification |
| `onInspect(frame)` / `onSelectInterval(begin,end)` | UI selection only; no audition switch, seek, revision, measurement, proposal or processing mutation |
| `onRequestSeek(frame)` from stopped ruler click or Seek button | Controller validates stopped state and current audition request, synchronizes identity and calls existing seek service; reports errors without fabricating movement |
| Section selection / numeric draft actions | Delegate to `beginSectionTrimDraft`, `setSectionDraftGain`, freshness checks, `applySectionTrimDraft` and discard; keep actual model method names/signatures authoritative |
| Playhead input | `StandaloneAudition::positionFrames()` and transport generation; a display timer observes it, never advances audio time |
| Follow / zoom / lane layout | Presentation-only state, separate from document audio/processing/requirements identities |

**Ruler behavior is an explicit requirement:** a primary click on the time ruler
positions the playhead while stopped. Measurement/waveform/gain-lane clicks only
inspect/select. Keep separate hit zones, and distinguish a click from a pan/resize
drag; do not transplant the plugin's double-click reset or plain-click visibility
shortcuts. Fit/visibility controls are explicit. During playback/buffering/preparation,
ruler seeking is disabled with a clear reason; it does not stop playback automatically.
Embedded video audio remains undecoded and unrouted. The existing
`StandaloneAudition::seek` resumes when `active()` was true: the new
controller must guard before calling it, not rely on `seek` to enforce stopped-only.
Also block conflicting jobs/loading and reject stale/unavailable audition identities.

Seek against `session.auditionRequest()` and existing routing, **not** the inspected
domain. Inspecting Current while auditioning Source does not switch processing.
Convert ruler pixels to the common viewport frame range, floor once to an integer
source frame, constrain to `[0,frameCount]`, and preserve EOF as an endpoint rather
than an extra sample. Numerical seconds/fractions use the existing rational parser.
A device-closed stopped seek updates the held position without opening a device.
Video continues to derive its frame from `audioFrame/sampleRate - startOnAudio` and
the same transport generation; signed offsets never move diagnostic/marker positions.

Waveform overviews must observe actual decoded programme samples: original Source,
or Current through `ProgrammeGainProcessor`, after composition/section/global processing
and before monitor routing/volume. Keep channels/roles and LFE; never average/downmix
for an overview without an explicit separately labelled display choice. Do not infer
processed extrema by multiplying Source bins across a varying envelope. The existing
`ProcessedAudioObserver::observe` is `noexcept` and also usable by real-time tests:
use bounded worker-side publication, not allocations/UI callbacks in the audio callback.
A new preview sink should have an explicit worker-only lifetime/error contract.

M/S points use current window-end coordinates and 100 ms plot cadence; detailed windows
retain half-open bounds. Sample-resolution maxima may lie between plotted windows.
EOF drain/padding stays identified and never extends programme duration. True-peak
locations retain eighth-frame coordinates. Decimated curves preserve visible extrema,
missing gaps and mixed gate membership; full-resolution values drive calculations.
Partial relative gates are unknown until final analysis, not estimates presented as
final inclusion. No rolling LRA/cumulative IL; final IL/LRA and under-60-second stability
remain separate from availability. Missing/silent/unavailable does not mean zero LUFS.

Automatic Source analysis starts once per established valid input identity, after
verification and explicit role/alignment/composition confirmation, when stopped and
no conflicting job runs. The current `beginMeasurement` explicitly stops audition;
auto-scheduling must first check idle/stopped and never stop playback opportunistically.
Reuse validated compatible cached measurements on reopening; a missing waveform offers
Build waveform without remeasurement. Cancel/failure suppresses auto retry for that
identity until explicit Analyze again or changed inputs. Do not initialize Current,
change applicability, infer a requirement or open devices automatically.

Use `SessionJob::finish` plus `acceptMeasurement` checks, not ticket ID alone. Runtime
`sessionId` is new on each opening; document/programme IDs and saved revision are
persistent and different. Preview jobs need equivalent generation guards on New/Open,
close, replacement, role/range/processing changes and cancellation. Source data remains
valid after Current changes; Current becomes stale/unverified. Requirements/policy
edits retain compatible measurements but stale dependent proposals/export approvals.
Video/monitor/navigation changes preserve audio identity. Publication of final data is
atomic; late previews cannot overwrite a final result or a new session.

## Draft, applied, verified and historical explanations

| State | Audio/result meaning |
| --- | --- |
| Annotation edit | Section/marker labels or boundaries changed; applied audio envelope unchanged |
| Draft gain or proposal | Preview/intent only; no processed audition or measured Current claim |
| Applied, UNVERIFIED | Existing stopped-only Apply changed the exact plan; old Current measurement/evaluation cannot remain current |
| Independently measured Current | Actual fresh result; show pass/fail/indeterminate against all selected rules, never substitute prediction/private assessment |
| Cancelled/failed remeasurement | Applied settings remain; Current stays unverified, no automatic rollback/retry |
| Already passes / identical apply | Genuine no-op: exact gains, measurements, evaluation, revision and prior historical explanation preserved |

Keep M17 LRA behavior: unset is default, zero is a strict ceiling, profile authority
is separate from user provenance. All selected rules passing gives a no-op; otherwise
changing automatic proposals with any selected LRA requirement remain unsupported.
Manual processing remains available. Navigation must not accidentally enable such
proposals, weaken artifact gates or add a hidden LRA ceiling.

At M17, schema 7 saves `lastAppliedProposalPolicy` and `lastAppliedProposalBinding`.
`lastAppliedProposalText` is transient. `applySectionProposal` records it; the direct
`applyTrim` path has no equivalent complete persisted explanation. M18 must cover both
changing automatic Apply paths. Capture baseline/provenance **before** processing
invalidation clears the measured Current/draft; commit history only with successful Apply.

Propose one independent versioned historical record (schema 8 session extension):
source/programme identity, before/after exact plan identities, relevant measured baseline
scalars/statuses/provenance, method, complete requirements/policy definitions and digests,
structured changed regions and reasons/limits/reserve disclosure. Bind reasons to their
original annotations/intervals, not renamed/moved current annotations. No private
candidate measurements/verdicts or old approvals may be installed as current verification.
Reopened labels distinguish “historical explanation for this plan”, “processing edited
since application”, and separately changed requirements/preferences. Restoring the same
plan does not restore verification. A new changing automatic Apply replaces the record;
manual edits retain its historical association. Source replacement must not attach it
to unrelated audio. No full undo/history browser.

Schemas 1–7 migrate without inventing missing reasons. Independently damaged history
can be discarded with a notice while valid settings/media/caches recover. Core corruption
still fails safely. Explicit Save explains newer-reader requirements; Save As preserves
older copies. Preserve coherent save revision and late-job guards: edits during save
remain unsaved. Reopen stopped, device-closed, routes-muted; saved verdicts are never
authority. Historical report text must be escaped and clearly separated from actual
Current/artifact evaluations. Existing Render/verify then Publish gates remain intact.

## Proposed initial layout and interaction checklist

The UI chat may refine spacing and grouping without changing contracts. Start with
one window/panel usable around 1200×800 and on Golem's smaller shared display. Use
scrollable/virtualized content and minimum lane sizes; do not require a new screen size.

```text
Programme / exact channel roles / requirements / saved-or-unsaved / job status
Inspect: Source | Current     Measured / UNVERIFIED / ANALYSING / unavailable
Audition: Source | Current    processing gain     device closed/open     Stop / Play
Zoom −/+  Fit programme/selection  Follow  Lane visibility  exact cursor / Seek
TIME RULER — stopped click seeks the auditioned domain
Markers / section boundaries                         | Section / window inspector
Actual per-channel PCM waveform                      | full-precision values
Measured M / S (or clearly unfinished preview)        | numeric section-gain draft
Applied section envelope + distinct draft overlay    | Apply / discard / remeasure
Horizontal scroll / coverage / exact selected interval| historical applied reasons
Separate LRA distribution/tails / final IL-LRA-TP and narrow selected-rule verdict
```

- Ruler click seeks once while stopped; verify first frame, EOF, extreme zoom and
  device-closed state. A drag or resize does not masquerade as a click/seek.
- Lane click selects/inspects; ranked-window/peak selection reveals exact locations.
  No hidden play, domain switch, gain change or measurement invalidation.
- Inspect and Audition labels can differ and remain obvious. Follow observes audition
  position; manual navigation disables follow without moving the playhead.
- Visibility toggles and lane resizing do not edit data; labels/units/status remain
  readable. Narrow extrema and gaps survive zoom; unsupported layouts stay explicit.
- Draft/applied overlays remain distinguishable without relying on colour alone.
  Apply is disabled while active; old passing results disappear after audio changes.
- Auto Source progress/cancel/retry, failed input and stale completion are visible;
  jobs never start playback or repeatedly restart after cancellation.
- Historical reasons survive Save/reopen; edited-plan/requirements notices are correct.
  Missing cache/waveform/explanation does not fabricate data or discard valid settings.
- No device opening during silent UI qualification. Coordinate later physical tests.

## Build and test entry points — commands for later approved work, not run here

Use existing JUCE 7.0.12 at `/Users/fuzz/Developer/LevelScope-dependencies/JUCE-7.0.12`,
commit `4f43011b96eb0636104cb3e433894cda98243626`; do not modify it. Keep C++17, Release,
universal arm64/x86_64, macOS 10.13 deployment, no-fast-math and fp-contract-off. No new
OS/dependency requirement is implied. CMake hashes all `Source` files into build identity,
so even UI/helper changes can invalidate conservative cache compatibility; do not falsify
the build ID to reuse caches. Migrate recoverable settings and remeasure where required.

After implementation approval, a fresh initial configuration can use:

```sh
M18_UI_ROOT=/Users/fuzz/Developer/LevelScope-worktrees/UITimeline
M18_BUILD="$M18_UI_ROOT/build/m18-dev-01"
cmake -S "$M18_UI_ROOT" -B "$M18_BUILD" -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_OSX_ARCHITECTURES=arm64;x86_64' -DCMAKE_OSX_DEPLOYMENT_TARGET=10.13 \
  -DFETCHCONTENT_SOURCE_DIR_JUCE=/Users/fuzz/Developer/LevelScope-dependencies/JUCE-7.0.12 \
  -DLEVELSCOPE_BUILD_VIDEO_TESTS=OFF -DLEVELSCOPE_BUILD_VIDEO_PAINT_DIAGNOSTIC=OFF \
  -DLEVELSCOPE_BUILD_MEASUREMENT_TESTS=ON -DLEVELSCOPE_BUILD_COMPLIANCE_TESTS=ON \
  -DLEVELSCOPE_BUILD_PROPOSAL_TESTS=ON -DLEVELSCOPE_BUILD_DIAGNOSTICS_TESTS=ON \
  -DLEVELSCOPE_BUILD_SECTIONS_TESTS=ON -DLEVELSCOPE_BUILD_AUDITION_TESTS=ON \
  -DLEVELSCOPE_BUILD_SESSION_TESTS=ON -DLEVELSCOPE_BUILD_EXPORT_TESTS=ON \
  -DLEVELSCOPE_BUILD_SECTION_PROCESSING_TESTS=ON -DLEVELSCOPE_BUILD_SECTION_PROPOSAL_TESTS=ON \
  -DLEVELSCOPE_BUILD_MONO_TESTS=ON -DLEVELSCOPE_BUILD_MARKER_TESTS=ON \
  -DLEVELSCOPE_BUILD_JOB_REQUIREMENTS_TESTS=ON
cmake --build "$M18_BUILD" --config Release --parallel 4 --target \
  LevelScopeStandalone LevelScopeMeasurementTests LevelScopeComplianceTests \
  LevelScopeProposalTests LevelScopeDiagnosticsTests LevelScopeSectionsTests \
  LevelScopeAuditionTests LevelScopeSessionTests LevelScopeExportTests \
  LevelScopeSectionProcessingTests LevelScopeSectionProposalTests LevelScopeMonoTests \
  LevelScopeMarkerTests LevelScopeJobRequirementsTests
```

For short edit/build iterations select only affected targets. `LevelScopeStandalone`
produces `LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app`; console
binaries are `<Target>_artefacts/Release/<Target>`. Add an explicit M18 timeline/progress
test target during implementation; it does not exist at this checkpoint. Do not build
`all`, the default target, `LevelScope` or `LevelScope_VST3`: CMake's plugin target sets
`COPY_PLUGIN_AFTER_BUILD TRUE`. Standalone support links video frameworks even in the
audio-only configuration; that is not authorization to load/decode/play video.

Existing short commands, run individually against a new owned results directory after
approval (relative executable names below refer to the build paths above):

| Test target / source | Selected entry point |
| --- | --- |
| `LevelScopeMeasurementTests` / `Tests/OfflineMeasurementTests.cpp` | `--report <new.json>`; sample and lifecycle primitives |
| `LevelScopeComplianceTests`, `LevelScopeProposalTests` | `--report <new.json>`; rule and global no-op/stale guards |
| `LevelScopeDiagnosticsTests`, `LevelScopeSectionsTests` | `--report <new.json>`; existing measured display/section contracts, including synthetic long analysis data |
| `LevelScopeAuditionTests` / `Tests/AuditionTests.cpp` | `--report <new.json>`; injected device, no `--long-source` during short phase |
| `LevelScopeSessionTests` / `Tests/SessionTests.cpp` | `--short --output <new.json>`; coherent migration/recovery |
| `LevelScopeSectionProcessingTests` | `--smoke --root <new-dir> --output <new.json>`; processing/seek/measurement parity |
| `LevelScopeSectionProposalTests` | separate `--smoke` and `--preferences`, each with `--root <new-dir> --output <new.json>` |
| `LevelScopeJobRequirementsTests` | `--root <new-dir> --output <new.json>`; M17 requirements/manual/export/UI guards |
| `LevelScopeMarkerTests` | `<new-dir> <repo>/Tests/Markers/nuendo-15.0.30-points.csv` |
| `LevelScopeExportTests` | `--artifacts <new-dir> --output <new.json>`; short real output and report gates |
| `LevelScopeMonoTests` / `Tests/MonoTests.cpp` | `--smoke --root <new-dir> --output <new.json>` includes bounded sparse container checks; omit `--long` / `--cinema` in the short phase |
| Programme/composition/multichannel/session-playback tests | Inspect/select the explicit audio-only phase arguments in their sources before invocation; do not start default long/video modes by assumption |

New M18 groups must test ruler versus lane actions, cached backgrounds, observation
on/off parity, preview ordering/cancellation, actual waveform extrema, historical
records and resource bounds. Expand affected existing caller cases, not just visual
snapshots. Use `Tests/EbuReferenceRunnerTests.py` and the explicit reference runners
below when shared measurement changes. `Tests/M17QualificationRunnerTests.py` is a
pattern for required-group accounting, not proof that an M18 manifest is complete.
Never invoke recursive historical runners or the frozen M17 runner against a new build;
it checks its own source/artifact identities and fixed selection.

Pinned reference entry points, only in the approved test phase:

```sh
M18_REFERENCE_RESULTS="$M18_UI_ROOT/build/m18-reference-01"
mkdir "$M18_REFERENCE_RESULTS"  # must be new; never overwrite an old run
python3 "$M18_UI_ROOT/Tests/run_ebu_measurements.py" --binary "$M18_BUILD/LevelScopeMeasurementTests_artefacts/Release/LevelScopeMeasurementTests" \
  --vectors /Users/fuzz/Developer/LevelScope/Tests/Ebu_Loudness_Testset_V05/ebu-loudness-test-setv05 \
  --report "$M18_REFERENCE_RESULTS/ebu.json" --arch arm64
python3 "$M18_UI_ROOT/Tests/run_reference_measurements.py" --binary "$M18_BUILD/LevelScopeMeasurementTests_artefacts/Release/LevelScopeMeasurementTests" \
  --vectors /Users/fuzz/Developer/LevelScope/build/measurement-reference/itu \
  --report "$M18_REFERENCE_RESULTS/itu.json" --arch arm64
```

Both vector directories above existed at handover inspection; use them read-only
and recheck paths/hashes before running. Raw vectors are ignored and will not appear
in a new worktree. Reference runners verify hashes
against `Tests/Reference/ebu-v05-manifest.json` and `itu-bs2217-manifest.json`.
Angle-bracket arguments in the entry-point table are placeholders to replace, not
literal shell redirection. Stop if a fresh output directory cannot be created.
Never change expectations to obtain a pass, copy media into Git or use a later rebuild
as the identity for old results.

## Qualification, resources and completion gate

Before coding, the UI chat records an impact manifest covering changed code, shared
callers/build dependencies, explicit fresh tests, justified reused evidence and pending
checks. Short Sequoia audio/UI work comes first; no Golem, video execution or three-hour
run during that phase absent a concrete explained dependency/failure. Build Intel
locally; do not call it actual High Sierra testing. Freeze a new combined M14–M18
candidate without replacing any old package, then plan actual High Sierra/cross-Mac,
selected H.264/ProRes and the consolidated long workflow **immediately after M18,
and no later than M19**. This handover does not lift deferrals or authorize execution.

Fresh coverage includes sample/result/window/gate parity with observation on/off,
all supported input forms/layouts/rates, exact frame navigation and audio-only seek/EOF,
stopped-only edits, Current invalidation/independent verification, preset/profile reuse,
strict actual-artifact export gates, new session history, failure recovery, silent native
UI and universal build/deployment inspection. Rerun pinned EBU/ITU on both actual Macs
if shared measurement observation is added, as anticipated. Retain missing/unstable LRA
semantics and M17's genuine no-op versus unsupported-correction distinction.

The combined long selection is one actual three-hour mono/192 kHz programme per Mac:
processing/proposal, independent verification, markers/requirements, timeline/progress,
explanation/session reopening, resource checks and verified PCM24 RF64 export beyond
4 GiB with independent decoder checks. Reuse that workload across requirements. Short
real stereo/5.1/7.1 plus three-hour-sized multichannel analysis/overview stress data
cover channel-density UI; the latter is not full-duration multichannel PCM evidence.
No automatic repeat of the full cinema proposal/export or older long workflows.

Selected video tests cover Source/Current, rational offsets, seeking, stop/resume/restart,
EOF, timeline open and analysis active/idle. Coordinate short muted physical runs
separately; prior milestone approval is not a new blanket device authorization.
Occasional picture skips are accepted and recorded. Audio interruptions, accumulating
drift, sustained freezes and identity/lifecycle/memory failures remain blockers. No
repeat-until-zero-skip loop. Full matrices/endurance are conditional on affected code
or concrete failures. Optional diagnostic tracing stays a separate executable; full
physical 7.1 listening remains separately pending.

| Budget inherited from 63 | Proposed acceptance |
| --- | --- |
| Preview delivery | <=5 updates/s, bounded/coalesced queue |
| View preparation / paint p95 | <=50 / <=100 ms; <=4,096 columns per curve/lane |
| Cancellation acknowledgement | <=1 s; stale completion never adopted |
| Both waveform domains | <=128 MiB; adaptive min/max levels, no complete PCM in RAM |
| Extra peak process RSS | <=256 MiB above equivalent final-analysis baseline |
| Observation overhead | <=15% wall time vs off, same host/fixture |
| Last-applied explanation | <=8 MiB, complete supported structured reasons, no silent truncation |
| Session / reports | Existing 128 MiB / 16 MiB; typical core save/parse <=5 s |
| Consolidated workflow | <=240 min Sequoia / <=300 min High Sierra |
| Disk planning | 2–4 GiB short fixtures, 1–3 GiB per package; >=40 GiB free for long workload |

Estimates, not measured M18 results: 15–30 minutes selected short automation plus
30–60 minutes UI per Mac; roughly 90–150/150–270 minutes for the one long workflow on
Sequoia/High Sierra, plus builds. Confirm actual storage/copy lifetimes and disk speed.
Owned partial caches/copies clean up on cancellation/write/disk-space failure; preserve
previous valid sessions. Large waveform construction cannot block UI/audio callbacks.

Until combined checks pass, report **“Sequoia tested; High Sierra, video and three-hour
checks pending.”** Preserve fresh/reused/deferred distinctions and historical failures.
Return exact source/build/app hashes, actual OS/CPU/architecture, explicit required-test
results, UI receipts, cross-platform comparison and observed resources. Main stays M13
until qualification and subsequent integration approval. No tag movement, history rewrite,
VST build/install, installation change or audio-system reconfiguration.

## Boundary and starter message

No freehand automation, rolling LRA/cumulative IL, new DSP, automatic segmentation,
selective LRA correction, coordinated neighbour strategy, dialogue/intelligibility
analysis, new marker adapters, 7.1.4/Atmos, extra export formats or VST3 parity. The
logic chat retains the next bounded region/LRA proposal work; it should not wait for
a comprehensive shared editor. Prepare its interface handback with immutable selected
regions, exact plan/state identities and current measured diagnostics, without secretly
implementing suggestions or correction here.

Suggested message to the new UI chat:

> Please take ownership of M18 planning using
> `/Users/fuzz/Developer/LevelScope/Docs/ProjectIndex/64_M18UITimelineHandover.md`
> and ProjectIndex 63. Start from tested M17 commit
> `e422f85ccbf7ae11876e6a31be2409937cb8381c`, with the documentation-only follow-ups;
> recommend `/Users/fuzz/Developer/LevelScope-worktrees/UITimeline` on
> `codex/ui-timeline-m18`. Review the handover and initial layout first and wait for
> my implementation approval. A stopped ruler click explicitly seeks; measurement
> lanes inspect only. Reuse standalone offline data and selectively adapt plugin
> presentation/navigation code. Preserve main at M13, all checkpoints/packages and
> the logic/video copies. Keep the Sequoia-first and combined post-M18 qualification
> plan. The logic chat owns subsequent LRA proposal work and shared-contract review.

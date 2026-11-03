# M18 dedicated UI/Timeline handover

**Current coordination summary — 3 November 2026**

M18 implementation has reached the short end-to-end acceptance gate on the
development M1 Max running Sequoia. The combined two-Mac qualification remains
pending. This summary gives the current position. Dated sections below retain
the earlier plans, build paths and per-slice evidence.

The active worktree is `/Users/fuzz/.codex/worktrees/a411/LevelScope`, branch
`codex/ui-timeline-m18`. The M18 implementation and development evidence are committed
as `4f02303a748b660c72471058111df8a1271758fc`, following
`cdac0c36b4259386417ecad507d107761eced864`, which integrates the logic module's
tested automatic Current initialisation. Preparation of the combined qualification
package has begun in `build/m18-qualification-2026-11-03`; its separate universal
build has completed. There is no completed frozen qualification package or published
M18 checkpoint yet. This handover update is subsequent documentation.
Main remains the qualified M13 commit `a9be070ab9eb5cb4f591ff9d991d8999c3a016d0`.

The advisory chat has since delegated the human-authorised combined qualification
to this UI/Timeline chat. The short driver pilots passed; long execution has not
started. [The explicit qualification selection and driver notes](../../Tests/Qualification/M18/README.md)
record the fresh groups, budgets, retained pilot failures and separate physical-device
coordination requirement. Actual test receipts determine completion.

Implemented behaviour includes:

- The approved native timeline layout, prominent Original / Current / Target
  values for IL, true peak, LRA, maximum momentary and maximum short-term loudness,
  direct target editing, complete presets including No Limit, immediate Unsaved
  preset naming, and resizable comparison, waveform and inspector areas.
- Automatic initial Current readiness from the identical completed Source
  measurement, with its origin inspectable. Applying later processing invalidates
  Current until an independent remeasurement succeeds.
- Progressive per-channel PCM waveform overviews from actual programme samples,
  with Source and Current sharing a combined 128 MiB budget. Reduction retains
  minima, maxima and their exact earliest sample positions. Deepest overview zoom
  exposes the represented interval. Measured M/S curves retain their own data,
  gates and coverage, with incomplete previews clearly identified. Build waveform
  can restore a missing overview while retaining valid measurements.
- Zoom, pan, Follow, an audio-timecode ruler above and seconds below, stopped ruler
  seeking through the existing transport, and Space Play/Stop. Inspection and
  audition remain independently selected. Main and Sections windows support
  click-drag-release section creation, draggable boundaries and numerical gain
  editing with explicit Apply and Discard.
- Separate remembered explicitly confirmed stereo and 5.1 role orders. Imports
  show suggestion provenance and label conflicts before Confirm roles and analyse.
  Reliable file labels take precedence. Reopened projects retain their own roles,
  physical output routing stays separate, and separate-file alignment confirmation
  remains required.
- Main-window Propose, Review, Apply and Discard, with disabled-action reasons,
  progress and specific outcomes. The gain lane distinguishes applied gains,
  manual drafts and automatic previews. Every proposed correction and its reason
  remains inspectable. Stale and rejected candidates cannot be applied.
- Complete last-applied global or section proposal explanations in schema-8
  sessions, using an independent checksummed attachment capped at 8 MiB within
  the existing 128 MiB session budget. Requirements, preferences, baseline,
  processing, reasons and available candidate results survive reopening. Subsequent
  edits identify historical association separately from Current verification.
- Fixes for preset-dialog dismissal, stale Target control availability, Space
  after Seek, Proposal limits window focus and repeated full measurement digests
  on the UI timer.

The existing proposal method still evaluates global gain first. The newly affected
duration allowance applies to local envelope changes, including transitions.
A section candidate uses desired global normalisation to identify local peak
deficits, measures the local plan, recalculates global gain, then measures the
combined candidate as required, with at most three private passes. Apply is explicit
and independent Current remeasurement follows. Automatic changing corrections with
selected LRA or M/S limits remain unsupported, although those limits are measured,
evaluated and saved, and an already-passing state can return no change.

The enabled Apple-silicon service now computes true peaks with bounded channel
parallelism for stereo and 5.1, preserving serial arithmetic, peak positions,
loudness order and proposal decisions. The Intel service retains serial execution.
Exact core/section/proposal equivalence, adversarial rounding and exception checks,
ThreadSanitizer and paired timing gates passed before enablement. The enabled build
also passed 946 focused measurement and proposal checks. [Prototype evidence](../../Tests/Reference/results/m18-bounded-peaks/README.md)

The latest acceptance used one 60-second stereo global-proposal job and one
60-second 5.1 section-proposal job, each repeated three times through the actual
native main component. Both passed fresh import, preview, roles, Original/Current
audition, section editing, explicit Apply, independent remeasurement, save/reopen
and verified PCM24 export. First partial waveform appeared at median 0.61 s / 1.23 s
from import. Play reached the first non-zero injected-output sample in 0.30 s /
0.71 s. The user selected injected output timing for this pass. Each stereo run
passed 61 checks and each 5.1 run passed 77. Independent bundle hashes and selected
PCM24 samples, including gain transitions, matched. [Acceptance evidence](../../Tests/Reference/results/m18-short-acceptance/README.md)

Current development app is
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m18-peak-dev/LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app`.
Its build identity is `93460e4031c2ec51765b4ad83dd7bc4e67551be1f94129b107ebbeea494a32ce`.
The universal app includes arm64 and x86_64, with Intel minimum macOS 10.13.
[Executable hashes](../../Tests/Reference/results/m18-bounded-peaks/candidate-03/development-build.json)

The next qualification block freezes the combined M14–M18 candidate, then runs
actual Intel High Sierra, cross-Mac comparisons and session exchange, selected
video regressions, and one consolidated three-hour mono/192 kHz workflow per Mac
including RF64 beyond 4 GiB and resource checks. Plan 63 defines the exact selection,
budgets and separately coordinated physical subset. The long tests have not started.
Integration and publication follow successful qualification and explicit approval.

Later work remains the visible-passage sample inspection reader, output selection
and opt-in reconnection, explicitly selected 5.1-to-stereo monitoring, and further
DSP/automation/meters. Faster initial waveform availability during import was
identified as a possible follow-up. Logic owns future region/LRA proposal work,
and the separate video worktree remains the video coordination point. Shared-file
changes should be coordinated with this UI/Timeline worktree before qualification.

Original handover 30 October 2026. **Implementation and the native layout approved
31 October 2026; M18 work is in progress.** The original handover performed no
build, qualification, device opening or Golem access. Current implementation
status and fresh validation are recorded at the end of this document.
This document supplements [63 — M18 plan](63_M18InteractiveTimelineAndProgressiveAnalysis.md).

**31 October planning update.** The managed UI worktree is
`/Users/fuzz/.codex/worktrees/a411/LevelScope` on `codex/ui-timeline-m18`.
The setup recipe below records the original recommendation. The revised M18 plan
specifies progressive waveform overviews, a later visible-passage sample reader,
and first-import readiness after integration of the logic module's tested change.
See [the revised interaction and implementation plan](63_M18InteractiveTimelineAndProgressiveAnalysis.md).

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
M18_HANDOVER_DOCS=$(git -C "$M18_LOGIC_ROOT" log --diff-filter=A -1 --format=%H -- Docs/ProjectIndex/64_M18UITimelineHandover.md)

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
Carry the frame of each minimum and maximum through every reduction. Label the
finest bin's actual frame interval and duration. Source and Current share a combined
128 MiB overview budget. The future sample reader uses the visible frame interval
and the inspected domain's exact processing identity.

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
identity until explicit Analyze again or changed inputs. On an eligible first import,
the integrated logic module initialises untouched, identical Current from the accepted
final Source measurement. Display both ready states and measurement origin from the model. After an applied processing change, Current needs
independent remeasurement. Applicability, requirements and device opening stay explicit.

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

## Build and test entry points — approved short implementation checks

Use existing JUCE 7.0.12 at `/Users/fuzz/Developer/LevelScope-dependencies/JUCE-7.0.12`,
commit `4f43011b96eb0636104cb3e433894cda98243626`; do not modify it. Keep C++17, Release,
universal arm64/x86_64, macOS 10.13 deployment, no-fast-math and fp-contract-off. No new
OS/dependency requirement is implied. CMake hashes all `Source` files into build identity,
so even UI/helper changes can invalidate conservative cache compatibility; do not falsify
the build ID to reuse caches. Migrate recoverable settings and remeasure where required.

The approved configuration uses the managed worktree path; the original recipe was:

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

## Approved native layout and initial implementation — 31 October 2026

The user approved the layout and implementation. The tested logic change
`9935ad725670958971c1aaed28d8e84f9b3d3e9a` is integrated as `cdac0c3`.
The UI consumes `acceptMeasurement` and never reinitialises Current afterward.
The approved waveform, 128 MiB combined budget and later visible-sample-read
boundary above remain part of M18's implementation contract.

The pinned comparison has Original, Current and Target rows and IL, maximum true
peak, LRA, maximum momentary and maximum short-term columns. Target values are
edited directly. Each optional value includes an explicit No Limit state. IL
retains an explicit tolerance. Editing a named profile immediately displays
Unsaved preset. Save as preset stores every field, including No Limit and M/S,
and adds the saved definition to the profile chooser. There is no separate
Requirements action and no LRA checkbox. User presets are user-specified
constraints, never presented as an authoritative broadcaster definition merely
because the user names one. A bounded vertical divider compresses comparison
row spacing without changing fonts or hiding values.

### First implementation slice and impact manifest

Integrate the logic handover, implement complete target definitions and their
measured rule/export guards, and replace the upper scaffold with the native
comparison/editor and resizable layout. Wire automatic initial Source scheduling
and measured M/S inspection using the existing services. Progressive PCM, the
combined gain lane and persisted applied explanations follow as subsequent M18
slices; this first slice is not M18 completion or qualification.

Fresh commands use `build/m18-dev-01`, Release, universal arm64/x86_64, Intel
10.13 and the existing JUCE 7.0.12 dependency. Build only named
`LevelScopeStandalone`, `LevelScopeJobRequirementsTests`, `LevelScopeSessionTests`,
`LevelScopeMeasurementTests`, `LevelScopeComplianceTests`, `LevelScopeProposalTests`,
`LevelScopeExportTests` and `LevelScopeSectionProposalTests` targets. Execute the
short argument forms documented above in a new owned M18 results directory.
Target changes affect evaluation, proposals, session/preset serialization and
export/report gates. Measurement arithmetic and device callbacks are unchanged
in this slice. Reuse the logic handover's provenance evidence, and rerun its
short integration checks on this build. No historical recursive runners, physical
audio, video decode, full-duration workflow or main integration are included.
Record fresh results and remaining M18 work here before handing back.

### Current first-slice behavior and format notes

- The native title bar reads LevelScope-Timeline. Open Mix offers the existing
  interleaved, separate-channel and confirmed composition paths. Audio setup
  remains a separate, explicit window. Space, Command-O and Command-S are wired
  outside text entry/modal contexts.
- Original and Current readiness are separate from rule verdicts. Complete
  measurement and initial-Source origin are available in comparison tooltips.
  The central surface currently displays the existing final measured M/S data.
  Progressive PCM and the consolidated gain/inspector lanes are still outstanding.
- All five Target fields are editable, including No Limit. IL includes its
  tolerance. New presets use definition schema 2 with every optional field stored
  explicitly. Legacy schema-1 EBU/LRA presets retain their original digest.
  A preset changes measurements' evaluation, never the measured audio identity.
- User constraints are labelled as user-specified. Maximum M/S use measured maxima
  and full-precision inclusive bounds. Missing required data is indeterminate.
  Zero selected constraints cannot produce a passing export approval. Changing
  automatic proposals with M/S or LRA constraints, or without IL/TP constraints,
  remain explicitly unsupported; all-selected-rules-pass is still a no-op.
- Session/core and programme schemas now use version 8 so older readers reject
  complete-target saves rather than silently losing new limits. Existing schemas
  migrate in memory; explicit Save retains the newer-reader notice. The optional
  historical explanation record is still to be added before M18 completion.
- Saving a session commits valid pending Target text and refuses invalid fields.
  The comparison divider varies spacing within 159–207 pixels, with fixed 16 px
  values and a minimum usable field height.

The first native launch was inspected with devices closed. Text encoding was
corrected after that inspection. Automated native input did not reliably change
JUCE controls through the computer-use tool, so interactive native verification
is not claimed from those attempts. The component/model regression tests cover
target edits, invalid input, no-op edits and compressed comparison geometry.

### First-slice validation — 31 October 2026

All nine selected short Sequoia/arm64 runs passed against the same measurement
build `d2b478709ba907c2eede7d9d4b6360205e2fa4f84b1c4f3d71d3104b6ea68b4b`.

| Group | Passed checks |
| --- | ---: |
| Measurement/services | 185 |
| Measured compliance | 182 |
| Global proposals | 303 |
| Short export/report | 4,872 |
| Section-proposal smoke | 366 |
| Proposal preferences | 256 |
| Complete targets/presets/UI components/scheduling | 489 |
| Short sessions, including initial Current provenance | 1,071 |
| Measured diagnostics | 369 |

The standalone and eight test executables built universal arm64/x86_64. Binary
inspection confirms Intel minimum macOS 10.13. Source and binary hashes, commands,
reports and retained failed development attempts are under
[`Tests/Reference/results/m18-development-2026-10-31`](../../Tests/Reference/results/m18-development-2026-10-31/README.md).
No physical device, video playback or full-duration test was run. Native launch
and layout inspection succeeded; interactive native UI qualification remains
pending. High Sierra execution, cross-Mac exchange and combined long/video checks
remain deferred as agreed. This is the first implementation slice, not completed
M18 qualification. Main remains qualified M13.

Next implementation work is the progressive per-channel PCM overview with exact
extrema frames and the shared 128 MiB budget, preview lifecycle/cancellation,
stopped ruler seeks and consolidated numeric gain/inspector lanes, followed by
historical applied explanations and the remaining M18 acceptance checks. The later
visible-passage sample reader remains out of M18 scope.

### PCM waveform slice — implementation and impact manifest, 31 October 2026

The next authorized slice adds actual-sample per-channel overview accumulation in
`OfflineAnalysisService` after the shared applied gain processor, a runtime
Source/Current cache, progressive publication and main timeline drawing/inspection.
The overview starts at 64 frames/bin and coarsens adaptively. A single 128 MiB pool
charges all live/retired worker and UI references: at most 60 MiB per pyramid plus
an 8 MiB display/metadata reserve. Immutable completed prefixes publish at most
five times/second, with immediate final publication. All pyramid/display merges
preserve minimum/maximum values and their earliest exact sample positions.
The inspector reports the represented half-open frame interval, duration and both
positions; unread passages remain pending. No sample-level detail read is added.

Source uses the verified confirmed programme. Current uses the shared actual
processor, including transitions. Accepted initial unity Current shares Source's
allocation only after verifying its measurement origin and audio identity.
Runtime/job/revision/input/roles/processing changes reject stale publications.
A cancellable Build waveform action replaces missing runtime data while preserving
accepted measurement snapshots. No waveform data is serialized in this slice.
M/S remain independently measured windows; no loudness or compliance is derived
from PCM extrema. Core measurement arithmetic and shared audio processing are unchanged.

Fresh required commands use `build/m18-dev-01`, Release universal arm64/x86_64,
Intel minimum 10.13, existing JUCE 7.0.12. Configure with the existing options plus
`-DLEVELSCOPE_BUILD_WAVEFORM_TESTS=ON -DLEVELSCOPE_BUILD_SECTION_PROCESSING_TESTS=ON`.
Build only `LevelScopeStandalone`, `LevelScopeWaveformCoreTests`,
`LevelScopeWaveformTests`, `LevelScopeMeasurementTests`, `LevelScopeDiagnosticsTests`,
`LevelScopeSessionTests`, `LevelScopeSectionProcessingTests`, `LevelScopeExportTests`.
Executables under `<Target>_artefacts/Release/<Target>` run as follows (the pure
core executable is directly under the build directory):

- WaveformCoreTests: no arguments; exhaustive PCM oracle, concurrency and pool bounds.
- WaveformTests: `--root <new-owned-root> --output <waveforms.json>`; Source/Current
  sample oracle, transitions, interleaved/assembly/composition, all four supported
  rates, mono/stereo/5.0/5.1/System-I 7.1, exact observer on/off measurement parity,
  cancellation/stale identity, allocation refusal, sharing, rebuild and UI inspection.
- MeasurementTests: `--report <measurement.json>`.
- DiagnosticsTests: `--report <diagnostics.json>`.
- SessionTests: `--short --output <sessions.json>`.
- SectionProcessingTests: `--smoke --root <new-owned-root> --output <processing.json>`.
- ExportTests: `--artifacts <new-owned-root> --output <export.json>`.

Reuse first-slice complete target/preset, compliance and proposal evidence because
those implementations are unchanged by this slice. Reuse qualified reference
measurement evidence because measurement arithmetic is unchanged; the fresh
measurement run and exact on/off parity guard the new service observation boundary.
Short automation budget 15–30 minutes plus build time. Preserve superseded evidence.
No default/all/plugin builds, physical devices, video, actual High Sierra or
three-hour end-to-end runs. This is development validation of the waveform slice,
not completion or full qualification of M18.

### Waveform slice completed and development-tested — 31 October 2026

The progressive PCM overview slice is implemented. Source and Current use actual
programme samples; Current includes the exact applied gain transitions. The initial
identical Current shares Source's allocation after provenance checks. Both domains,
obsolete worker references and display reserve remain within the combined 128 MiB
budget. The main timeline shows per-channel coverage, keeps both extrema and their
sample positions through reduction, and identifies the finest represented interval.
A missing overview can be rebuilt without changing valid measured results. M/S still
use their independently measured windows. The later visible-passage sample reader
remains deferred.

All seven required short development groups passed against build
`d6a5d2e5bcb87f26560dadbe17deb2fdd1fca7f00de08e533df2c8ee474fb80f`: waveform core 240,560,
waveform integration/UI 1,173, measurement 185, diagnostics 369, short sessions 1,071,
section processing smoke 334 and export 4,872 checks. The pure waveform core also
passed ThreadSanitizer and AddressSanitizer/UndefinedBehaviorSanitizer with no
findings. Final short-fixture observation overhead was 0.36%; the largest
sampled preparation + software paint was 4.479 ms. Every final release binary
is universal, with Intel minimum 10.13. Source/test/binary hashes and commands are
in [the waveform evidence](../../Tests/Reference/results/m18-waveform-2026-10-31/README.md).
`candidate-02` is final; the passing earlier candidate is retained and superseded
by a clarification that waveform rebuilding does not imply measurement invalidation.
Component render captures were inspected; no native interactive qualification is
claimed. No device, video, default/all/plugin target, main integration or other
worktree was changed by this slice.

Remaining M18 implementation includes progressive measured M/S publication,
time-ruler/transport and gain-lane/inspector integration, further resizable zones
and the bounded saved explanation history. Actual High Sierra, cross-Mac,
physical/video and full-duration/RSS qualification remain in the agreed combined
block after M18, no later than M19. This slice does not mark M18 complete.

### Progressive M/S and navigation slice — authorized follow-up

Scope: optional worker-only observation of the existing measured 400 ms M / 3 s S
windows at the existing 100 ms cadence; bounded sealed-prefix publication with
pending gate membership; final snapshot replacement and stale/cancel rejection.
Add an audio-timecode ruler (HH:MM:SS.mmm, with no invented SMPTE frame rate), retain
seconds below, pan/zoom/follow and a distinct stopped-only Seek action. Route seek
through the existing audition service using the audition domain, independently of
inspection. Reject ruler drags, playback/preparation/jobs/loading and stale requests.

Fix stale Target availability after transport changes. Draft presence itself does
not lock Target selection. Extend upper divider travel from 31 to 78 pixels at its
default height, retaining 16-point metrics; keep its information line below the
divider at a fixed height. Add a divider below the timeline scrollbar, giving the
released height to the existing lower inspector rather than changing fonts.

Focused validation selection: new LevelScopeTimelineProgressTests (exact observer
on/off results/windows and cadence, silence/gaps/drain, publication/lifecycle,
transport guards using injected devices, independent inspection, navigation and
layout/Target availability); LevelScopeMeasurementTests; LevelScopeDiagnosticsTests;
LevelScopeWaveformTests; LevelScopeJobRequirementsTests. Build only these named
targets plus LevelScopeStandalone, universal with Intel minimum 10.13. Reuse prior
unchanged export/session/proposal and waveform core/sanitizer evidence. Measurement
arithmetic/order/drain are unchanged; the new observation boundary gets fresh exact
on/off parity across supported rates/layouts and the existing numerical suite.
No historical full runners, physical device or video runs, actual High Sierra or
three-hour workflow. Those qualification deferrals remain as agreed. Record the
exact development-app path and a short interaction checklist at completion.

### Progressive loudness and navigation — development handover

The standalone now overlays the actual 400 ms M and 3 s S measurement windows
as analysis runs. Partial coverage, pending M gate membership and ungated S are
identified in the plot and inspector. Cancellation clears the partial data.
Final programme values, gate membership and identical-Source Current adoption
still require an accepted measurement. Source and Current waveform storage
remains subject to the existing combined 128 MiB budget.

Development app: `/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m18-dev-01/LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app`.
Use the [short manual checklist and rendered examples](../../Tests/Reference/results/m18-progress-navigation/README.md#manual-interaction-checklist)
to try graph inspection, stopped ruler/numeric seek, zoom/fit/pan/follow and the
two dividers. The Target controls now follow transport, loading and job state;
neither discarding a draft nor toggling programme applicability is required to
unlock them. At maximum compression the metrics keep their font size and the
Target information line keeps its own height beneath the upper divider.

[Final evidence](../../Tests/Reference/results/m18-progress-navigation/final-summary.json):
2,692 focused checks and four pinned reference curve cases passed. The app and
five test executables contain arm64/x86_64 and declare Intel minimum 10.13.
The app was not qualified on actual High Sierra. The combined physical/video,
cross-Mac and full-duration/RSS qualification remains deferred as agreed.
Gain-lane/section inspection and saved applied explanations remain later M18
work; the visible-passage sample reader belongs after M18.

### Section gain lane and transport focus — impact selection

Add a Current section-gain lane on the shared frame viewport and expose numerical
staging, explicit Apply and Discard in the main timeline. Use existing M12
validation, canonical plans, transition endpoints and invalidation. Selection is
inspection only. Test lane alignment, boundary/EOF selection, exact values,
stale/unsupported/active rejection, no-op, discard and independent remeasurement.
Test the stopped-Seek → Space → Space path across child focus and retain normal
text editing. Clarify playback preparation wording without changing verification,
device, decoder or DSP behavior. Build the named standalone, timeline interaction
tests and section-processing smoke target only. Reuse earlier measurement,
waveform, session/export and proposal evidence for the unchanged paths. Actual
High Sierra, cross-Mac, physical device/video and full-duration/RSS qualification
remain deferred; saved proposal history remains separate M18 work.

The gain/focus slice is implemented and passed **1,049 focused checks**: 715 for
timeline/progress/gain/focus and 334 for the existing section-processing smoke
suite. The main lane and Sections window share the same full-precision draft.
Applied gain is solid, draft is dashed; click either the lane or its section
navigation arrows to inspect without changing audition. Use Return or **Set
draft**, then explicit **Apply gains** or **Discard gains**. A changing Apply
marks Current unverified until remeasurement; a no-op and Discard preserve it.

The dense 10,000-section fixture prepared and painted in 6.041 ms on the
development Mac. The native JUCE focus regression executed both Space presses
after Seek, including the Seek button becoming disabled during playback. The
audio service continues to verify a new private playback copy; only its message
changed. For app path, captures, commands, hashes and the manual checklist, see
[the gain-lane evidence](../../Tests/Reference/results/m18-gain-lane/README.md).
Development only: this does not mark M18 or the combined qualification complete.

### Manual sections, preset dismissal and remembered roles

M18 adds one-gesture waveform annotation and shared stopped ruler seeking and
Play/Stop in the Sections window. Fresh stereo and 5.1 imports propose file
roles from explicit labels, then remembered confirmed order, then a suggested
layout for review. Both import paths require confirmation. A separate-file
group also requires explicit alignment. Reopening restores the project's roles.
Focused checks cover modal lifetime, section geometry, shared transport, role
provenance, conflict, alternating formats and project reopening. Reuse M18's
existing gain, measurement and audition evidence for the unchanged processing
paths. The combined hardware, cross-Mac, High Sierra and full-duration checks
stay scheduled after M18, no later than M19.

Development build `5811184d944d23b6b660dd6a4da55cfa24f81ffaa318414805d6c473fed92067`
passed 1,378 focused checks. These include both preset dialog exits, drawing
and moving shared annotations, stale gesture rejection, Sections ruler seek
and Play/Stop with injected audio, alternating stereo and 5.1, WAVE/FLAC/AIFF
label precedence, separate-file alignment and a real saved-project reopen.
The app and test binaries are universal with Intel minimum 10.13.
[Development evidence and manual checks](../../Tests/Reference/results/m18-sections-roles/README.md).
[Later playback roadmap](63_M18InteractiveTimelineAndProgressiveAnalysis.md#playback-roadmap-following-m18).

### Main proposal workflow and saved explanation — impact selection

Expose the existing global-first, section-led service through a main-window
Propose / Review proposal / Apply proposal / Discard proposal panel. Keep the
private assessment distinct from measured Current and require independent
Re-measure Current State after changing Apply. Outcomes identify readiness,
no-op, unsupported requirements, preference rejection, stale inputs, cancellation
and failure. Rejected/stale selections remain labelled and inspectable. The gain
lane shows applied, manual and automatic gains separately, including proposed
global gain and exact local-plus-global inspection through transitions. Original
proposal annotation intervals remain the review binding.

Both changing automatic Apply paths capture the immutable baseline before
invalidation. Session schema 8 gains an independent version-1, checksummed binary
attachment for the last successfully applied explanation. It has an 8 MiB limit,
inside the unchanged 128 MiB session budget; the 4 MiB core limit and independent
measurement-cache recovery remain unchanged. The record holds exact before/after
plans and identities, original annotations, structured corrections, full target
and preference definitions/digests, baseline and available candidate scalar
results/statuses/provenance, predictions and narrative. It contains no measured
window population or verification authority. Legacy sessions gain no invented
reasons. Bad history can be discarded independently with a notice. New changing
Apply replaces the record, no-op/Discard preserve it, and source replacement
clears it. Later edits identify historical processing/requirements/preferences
separately from Current's actual verification state.

Focused selection: new `--m18-workflow` checks in LevelScopeSectionProposalTests,
its existing `--smoke` and `--preferences` modes, global proposal tests, session
`--short`, timeline interaction tests and target requirement tests. These cover
actual numerical proposal generation/application, explicit main-panel actions,
automatic lane/rejected/stale views, 76.13% vs 20% diagnostic selection, 18 LU LRA
capability rejection, no-op, global and section history, independent remeasurement,
save/reopen, late-save edits, corrupt/missing/unrelated history and complete
10,000-row access within the explanation budget. Build only the named standalone
and selected test targets. Measurement/DSP arithmetic, device/routing, video,
waveform storage and actual-artifact export gates are unchanged. The High Sierra,
physical device/video, cross-Mac and full-duration/RSS qualification remains in
the agreed combined block after M18, no later than M19.

The new dense automatic-preview test identified quadratic geometry-key string
construction during freshness checks. Use a linear stream for the identical
established bytes and reuse that helper for manual geometry keys. A 10,000-row
integer-format oracle guards identity compatibility; no revision, bound input,
geometry or source-stamp check is removed from Apply.

### Proposal workflow and historical explanation — development handover

This slice is implemented. Final candidate 04 passed **3,737 focused checks**
against build `875fe8b80d1813513fc09c6ee0a4e5bd8e246bbd56d83128ac97a2b7794dbd8a`. The development app remains
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m18-dev-01/LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app`.
The final app and five selected test binaries are universal with Intel minimum
10.13. The dense automatic controls/lane test measured 14.042 ms p95;
10,000 complete corrections and maximum UTF-8 names fit the 8 MiB historical
attachment budget. Complete numerical/state, corruption recovery, save/reopen,
no-op, stale and main-control checks are recorded with the commands/hashes in
[the final evidence and manual checklist](../../Tests/Reference/results/m18-proposal-workflow/README.md).
Native main-window and rendered component layout were inspected. This is focused
development validation, not completion of M18 qualification. The High Sierra,
physical/video, cross-Mac and full-duration/RSS deferrals remain unchanged.

### Live proposal UI hang and limits window

Sampling the unresponsive M18 app showed repeated full measurement digests on the
UI timer. Proposal queries now reuse the digest of each immutable Source/Current
snapshot without retaining the windows or caching an Apply decision. Raising
Proposal limits now reuses its open editor and preserves unconfirmed values.
The original app recovered on Discard. [Focused evidence](../../Tests/Reference/results/m18-proposal-refresh/README.md)
records 1,136 passing checks, including live context/stale payload checks, large-window
refresh (1.02 ms p95), and native window focus.
The corrected development app is `/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m18-refresh-fix/LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app`. Qualification deferrals stay unchanged.

### Offline verification, measurement and Propose profiling

On the development M1 Max, a 60-second stereo programme took 1.40 s to analyse and 4.06 s for a fresh three-pass section proposal. The matching 5.1 jobs took 4.19 s and 12.34 s. The meter took 84% of Source analysis in both layouts. The peak routine appeared in 79–81% of whole Source worker samples and 87–89% of fresh proposal worker samples. The recommended first experiment is a bounded per-channel true-peak prototype with exact serial-equivalence checks. A duplicate hash accounts for about 6.5% of Source analysis and 2.2% of fresh section Propose on these warm local PCM files. This profiling used diagnostic executables and synthetic files on the development Mac. [Stage timings, method and recommendation](../../Tests/Reference/results/m18-offline-profile/README.md)

### Bounded per-channel true-peak development build

The Apple-silicon standalone now uses a bounded per-channel true-peak path for stereo
and 5.1. It preserves the serial filter arithmetic, frame/eighth-frame positions,
observer semantics and ordered loudness reduction. The development Mac passed exact
core, section-peak and proposal comparisons and matched stereo/5.1 timing gates
before enablement. The compact lane count preserves the serial engine layout.
The enabled app is `/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m18-peak-dev/LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app`.
[Timings, qualification boundaries and binary identity](../../Tests/Reference/results/m18-bounded-peaks/README.md).

### Short end-to-end acceptance on the development Mac

On 2 November 2026 the existing bounded-peak build passed one stereo global proposal and one 5.1 section proposal from fresh import through independent Current remeasurement, session reopen and verified PCM24 export. Each 60-second job ran three times. Median import-to-first-waveform was 0.61 s for stereo and 1.23 s for 5.1. Play-to-first-non-zero sample at the injected output was 0.30 s and 0.71 s respectively. The app restored the applied explanation and verified Current independently. Selected exported PCM samples, including the 5.1 gain transitions, matched the independent checker exactly. The combined two-Mac long, video and exchange qualification remains pending. [Short acceptance timings, workflows and evidence](../../Tests/Reference/results/m18-short-acceptance/README.md)

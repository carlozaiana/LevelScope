# M18 proposal — interactive audio timeline and progressive analysis

Revised 31 October 2026. **Implementation approved 31 October 2026; implementation in progress.**
The dedicated UI/Timeline handover in [64](64_M18UITimelineHandover.md) records
ownership, the separate M17-based checkout, plugin reuse assessment and the user
clarification that a stopped time-ruler click explicitly seeks. Preparing that
handover authorizes documentation only.
Continue from published M17 candidate-04, checkpoint
[`e422f85ccbf7ae11876e6a31be2409937cb8381c`](https://github.com/carlozaiana/LevelScope/commit/e422f85ccbf7ae11876e6a31be2409937cb8381c),
tag [`m17-job-lra-presets-sequoia-2026-10-30`](https://github.com/carlozaiana/LevelScope/tree/m17-job-lra-presets-sequoia-2026-10-30).
The publication/proposal follow-up changes documentation only. The logic checkout
remains on `codex/logic-module`. The dedicated UI/Timeline worktree is
`/Users/fuzz/.codex/worktrees/a411/LevelScope`, branch `codex/ui-timeline-m18`. Main remains qualified M13
`a9be070ab9eb5cb4f591ff9d991d8999c3a016d0` until combined qualification and subsequent
integration approval. Preserve all checkpoint tags, packages and the video worktree.

## Intended outcome and bounded scope

Make the existing measured results and section processing practical to inspect:
**Confirm programme → automatic Source analysis with visible progress → Source and
untouched, identical Current ready → select delivery requirements → Propose and
inspect → explicit Apply → independently measure edited Current → audition or
render and verify.** The tested logic-module handover is integrated and supplies the initial Current step. Playback is always an explicit
separate action. A user can cancel automatic analysis without a restart loop.

Deliver one reusable audio timeline with:

- Source/Current inspection selector; waveform overview, measured momentary and
  short-term curves, marker/section boundaries, selected interval and playhead on
  the same audio-source axis. Display the inspected domain separately from the
  auditioned domain; inspecting Current must not switch playback implicitly.
- Horizontal zoom, pan/scroll, fit programme/selection, exact numerical inspection
  and explicit **Seek to cursor**. **Clicking the time ruler while stopped explicitly
  positions the playhead through the existing seek service.** Clicking measurement
  lanes selects/inspects without seeking, playing, changing processing or opening a
  device. While playing, ruler seek is disabled; it must not stop/restart playback
  implicitly. The playhead observes the authoritative audio presentation position.
  Existing transport capabilities remain available.
- An applied section-envelope lane plus a visually distinct draft overlay. Selecting
  a section exposes its existing numeric gain control and full-precision stored
  value. Use M12's existing draft/Apply/discard operations and unchanged 100 ms
  transitions, ±12 dB manual section range and all-channel processing. No freehand
  breakpoint automation or new DSP. Annotation split/move/rename/import remains
  separate; changing an annotation never moves the applied audio envelope.
- Existing ranked windows and per-channel peak-location selection reveal their exact
  intervals/positions. Keep the LRA distribution/tails in a distinct panel linked
  by selection, not overlaid as a purported time-local LRA curve.
- A bounded saved explanation of the last explicitly applied automatic proposal,
  accessible on reopening alongside its actual applied-plan identity.

The first editable workflow remains mono with confirmed Left/Right/Centre, stereo
and conventional 5.1, matching existing section-processing/export capability.
Existing explicit System I 7.1 measurement/global-trim/audition inspection remains
supported; the timeline must not imply 7.1 section-processing/export support.
Per-channel data stays in confirmed file/programme order, including LFE. Interfaces
remain count/role driven with capacity for twelve channels, without a new 7.1.4 claim.
Interleaved, aligned separate-file and confirmed M&E/dialogue programme definitions
use the shared programme reader; there is no new importer or stem mixer.

## Existing implementation and proposed architecture

The M17 repository already contains `MeasurementTimeline` with measured M/S curves,
zoom, pan, numerical inspection, gate-aware extrema reduction and highlighted
windows. It deliberately has no session/transport reference. `StandaloneDiagnostics`
binds immutable final snapshots; `StandaloneSectionsPanel` has a measured timeline,
section list and gain draft controls. `ProgrammeGainProcessor` supplies the common
offline/audition/export processing, and `StandaloneAudioPanel` exposes the existing
seek service. Reuse these boundaries rather than create another measurement engine.

Add a thin timeline presenter/controller that supplies immutable view data and
explicit callbacks for inspection, section selection, draft actions and seek.
Keep sample-frame transforms, waveform data and view preparation testable apart
from native UI. Leave playback callbacks, buffering and the presentation clock
unchanged; the displayed playhead observes the authoritative audio position.
A display timer must never become a second transport clock.

### M18 waveform overview

Build the progressive overview from actual samples on the analysis worker, per channel
in confirmed programme order. Source represents the original confirmed
programme, including a composition's unity sum. Current represents the exact applied
plan through the shared processor, before monitor controls. Do not manufacture a
Current waveform by scaling Source extrema across a time-varying transition.
Draft envelopes are graphics only; they do not acquire a measured waveform or verdict.
Reuse the existing worker-side sample observer boundary where appropriate. Any
production observation must be optional, non-mutating and bounded; no test-capture
buffer may grow with the full PCM duration.

Each bin covers a half-open interval of absolute source frames. It stores the minimum
sample with its frame and the maximum sample with its frame. Merging bins
or reducing them into display columns carries the winning pairs together. Use the
earliest frame on equal extrema so positions agree across block sizes and zooms.
The finest stored bin width adapts to duration and channel count. All levels, extrema
positions and retained overview batches fit within **128 MiB for Source and Current
combined**. Identical domains may share immutable bins when the backend confirms
sample identity, with separate bindings for each domain.

At the deepest overview zoom, label the inspected channel and actual represented
interval `[begin, end)`, its frame count and duration, and each extremum at its own
frame. For example, 480 frames at 48 kHz span 10 ms. The last bin ends at real EOF,
and unfinished coverage is marked pending. The cursor can identify an exact frame
even when only bin extrema are available there. M/S continue to use their own measured
windows, plotted at their ends with their established gate and gap semantics.

Generate overviews during normal measurement, avoiding a second full read solely
to draw a completed measurement. If a session has valid measurements but
no overview, offer a cancellable **Build waveform** job using the same verified input
and processing; do not invalidate or repeat valid measurements. The first version
keeps overviews as replaceable runtime data, not a mandatory persistent cache or a
new media copy. Reuse an existing verified source only when its identity/lifetime
contract permits it. Missing/changed/unreadable members fail the complete programme;
there are no silent substitute channels.

### Sample inspection in a later milestone

Design the viewport around a domain, channel roles, exact visible frame interval and
requested resolution. Beyond the finest overview resolution, a future presenter can
request individual samples for that interval. Source reads the verified, confirmed
programme through the existing reader. Current reads the same programme and applies
the exact section and global processing at the passage's absolute frame positions,
including ramps crossing either edge, before monitor gain or routing.

The future read runs asynchronously with bounded blocks and a byte cap on outstanding
and retained samples. Keep the labelled overview visible while detail loads. Pan, zoom,
domain or plan changes supersede older requests. Accept replies only for the matching
runtime session, programme, applied plan and viewport generation. The returned samples
serve inspection only. Implementing and qualifying this reader belongs after M18.

## Progressive analysis and first-import readiness

After successful programme verification and explicit role/alignment/composition
confirmation, schedule Source measurement once for the current input identity when
stopped and no conflicting job is active. On session reopen, reuse compatible validated
measurements; start only when the required Source measurement is absent. Never invent
roles, applicability or job requirements. Expose Cancel and explicit Analyze again;
cancellation/failure suppresses automatic retries for that unchanged identity.
On first import, show Source as analysing and Current as awaiting the completed
Source measurement. The logic-module chat owns automatic initialisation of untouched,
sample-identical Current. The tested change is integrated. Bind
the timeline to the session model's accepted result. Show **Source ready** and
**Current ready** after that single completed Source analysis. In the inspector,
identify Source as measured directly and Current as initialised from that Source
measurement, with the originating measurement identity available. Readiness describes
measurement availability. Selected-rule verdicts are displayed separately.

The timeline consumes the backend's eligibility and provenance. A partial, failed,
cancelled or stale Source result cannot initialise Current. A later edit that changes
applied audio makes Current **UNVERIFIED** until an independent Current measurement
succeeds. Returning gain to unity after an edit still requires that measurement.
Preserve Source and show the new measurement origin when Current is remeasured.

Publish bounded incremental previews through immutable, generation-bound batches:
waveform coverage and the existing 400 ms M / 3 s S window values at their established
cadence. The unprocessed region stays visibly empty, not silent. Progress distinguishes
media verification from analysis; I/O or hashing may legitimately precede any curve.
Throttle UI updates; never copy the entire growing history on every processed block.
A small optional measurement-window observer may be needed because the current engine
only exposes completed `Result` vectors. It must leave arithmetic, operation order,
finish/drain semantics and final outputs identical with observation on or off.

Label all partial data **ANALYSING — incomplete, not a programme result**. During
analysis, relative-gate membership is not final: show it as pending, not included or
excluded under a gate inferred from the partial programme. Final acceptance atomically
replaces previews with the ordinary validated snapshot/diagnostics. Cancellation,
replacement, role/processing edits, New/Open and closing reject late batches and finals.
Partial data never satisfies compliance, proposals, Current verification or export.

**Recommended first-version boundary:** no rolling LRA or cumulative IL estimator.
Display final programme IL/LRA and stability only when complete; show “pending” before
that. This explicitly stages those optional displays from the wider design in [58](58_LraWorkflowAndInteractiveTimelineDesign.md).
Adding a new rolling statistic would require its own method/window definition and
validation, and is unnecessary for useful M/S inspection or subsequent LRA proposals.
Unavailable or unstable values, finite excluded windows and genuine zero PCM remain
distinguishable. The approved Target editor adds explicit user-specified IL, true-peak, LRA, maximum M and maximum S constraints. No limit is inferred from programme measurements.

## Identity, editing, transport and reopening

Bind previews/overviews to runtime generation, exact source/member identities, ordered
roles, programme range/rate, domain, applied processing identity and format version.
A new job ticket is never restored from disk. Source remains reusable after Current
gain changes; stale Current curves/overviews must not be presented as current. A job
requirement/profile edit only refreshes annotations/evaluation for compatible measured
audio. Monitor gain, routing and video attachment/offset cannot change audio identity.

Retain M4's exact half-open window intervals, plot-at-window-end convention, true-peak
subframe coordinates and EOF-padding labels. Clip highlights to real programme duration.
Keep audio frames authoritative for section edits and seeks; a timecode display uses
only an explicitly defined rational mapping. Video offsets do not move audio markers.
Use the existing signed-offset, held-picture and presentation-clock contract.

Draft gain changes, Apply, source/role/routing changes and audition-domain switches
remain stopped-only. Numeric editing must not round existing full-precision gains.
Inspection can continue while playing; it never changes processing or transport.
Monitor volume remains adjustable. Devices open only by explicit action, with the
existing overload stop/recovery policy and no hidden clipping or attenuation.

For explanation persistence, introduce one bounded, versioned historical record of
the most recent changing automatic Apply, covering both the direct global path and
section/global-only proposal path. Store reasons/constraints, relevant baseline scalar
metrics and their availability/provenance, source and before/after plan identities,
method version, complete selected requirement and policy definitions/digests, and
structured changed-region information. Do not save private candidate sample histories,
measurements or verdicts as current verification. “Already passes” remains a genuine
no-op and does not replace the previous record or dirty the session.

On reopen label it **Historical explanation — applies to this processing plan**, or
**Historical explanation — processing edited since application**. Also identify changed
requirements/preferences independently; they must not rewrite past reasons. Matching
a plan again does not make the old result a present verification. Manual edits retain
the record with its association; a new automatic Apply replaces it. Source replacement
must not attach old reasons to unrelated audio. An undo/history browser is excluded.

Propose session schema 8 for this independent optional record. Schemas 1–7 retain exact
settings with “explanation unavailable in this older session”; do not fabricate it from
current values. Opening never rewrites originals; explicit Save explains the newer-reader
requirement and Save As preserves an older copy. An independently damaged history record
is discarded with a notice while valid core settings/caches remain recoverable. Keep
coherent-save revision guards, safe reopened stopped/device-closed/muted state and
measurement-based reconstruction of all verdicts. Export/report explanations remain
historical and escaped, and cannot authorize publication.

## Implementation order

1. Establish the passive frame viewport and its overview resolution/coverage contract.
   Reserve the visible-passage request shape for the later inspection reader.
2. Accumulate per-channel extrema and their frames during analysis. Connect bounded
   preview updates and the finest-bin inspector to that data.
3. The tested logic-module change and evidence have been reviewed and the change
   integrated in the UI worktree. Bind first-
   import readiness and measurement origin to its session state. This dependency
   gates the two-ready workflow. Eligibility and adoption belong to the logic module.
4. Connect existing stopped seek, numeric drafts, explicit Apply and remeasurement.
   Complete historical explanation persistence, then run the impact-selected checks.

## Acceptance tests

1. **Final-result parity:** observer on/off and differing batch/block sizes produce
   identical final samples, windows, gates, scalar metrics and identities. Run pinned
   EBU/ITU references if the window observer touches shared measurement code, as
   anticipated. No preview is accepted by any measured-rule/proposal/export gate.
2. **Waveform truth:** actual Source/Current min/max compared with processed samples
   for mono/stereo/5.1 and supported 7.1 global trim, composition, assembled members,
   all rates, LFE, positive/negative gain, transitions, silence, over-full-scale finite
   values, odd EOF and narrow impulses. Include anti-phase stereo to prove there is
   no implicit downmix. Compare both extrema and their exact positions after each
   reduction, including ties, impulses on bin edges and a shortened final bin.
   Check the finest-zoom interval label and combined storage with both domains loaded.
3. **Progress/lifecycle:** deterministic cancellation mid-verification and mid-analysis,
   missing/changed members, rapid reopen/source changes, stale batches, no retry loop,
   closing during work, bounded queues and allocation failures. Cached reopen performs
   no unnecessary measurement. Missing waveforms cannot invalidate sound measurements.
   With the tested logic change integrated, confirm one completed Source pass gives
   both ready states and the correct origin. Reject premature or repeated adoption,
   including after an audio edit. Verify independent Current remeasurement and origin.
4. **Navigation:** frame-zero/end, fractional-time floor convention, long coordinates,
   exact window selection, no phantom edge markers, EOF padding and viewport clipping.
   Explicit seeks through the existing service preserve audio-only continuity, duration,
   stop/restart/EOF and pre-monitor sample parity. Distinguish the explicit stopped
   ruler click (one seek) from lane clicks/zoom/selection (no seek). Test ruler clicks
   during playback and jobs, exact first/last frames, and device-closed seeking.
   No navigation changes measurements, proposal validity or processing.
5. **Editing/requirements:** annotation edits preserve audio; draft/apply/discard use
   existing M12 semantics; no-op preserves revisions; Apply makes Current UNVERIFIED.
   Cancellation leaves the applied gain in place; completed failing measurement shows
   the real failure. Profile/preset changes reuse measurement but invalidate outdated
   approvals. M17 selected-LRA changing proposals remain explicitly unsupported.
6. **History/sessions/reports:** both proposal paths round-trip exact provenance, no-op
   does not mutate it, subsequent audio/requirements edits show the correct association,
   damaged history recovers safely, old schemas migrate without rewrites, old readers
   reject new saves, coherent saves survive late jobs. No old private result is restored
   as verification; report text escapes names/labels and retains narrow result wording.
7. **Native UI and resources:** actual zoom/selection/seek/draft/Apply/remeasure/reopen,
   progressive/cancel/error states, dense marker/gain lanes and channel labels. Lists
   remain virtualized. Device-closed silent checks and injected devices/captures first;
   physical devices or listening only after separate coordination.

## Qualification sequence and budgets

**Execution update — 3 November 2026.** The human-authorised combined selection
passed on actual arm64 Sequoia and Intel High Sierra: fresh short/reference gates,
bidirectional sessions, selected injected/native/muted-physical video, actual
three-hour mono/192 kHz RF64 workflows, independent export samples beyond 4 GiB,
complete plain/observed numerical comparisons and synthetic six/eight-channel
display stress. [The frozen candidate and results](../../Tests/Reference/results/m18-qualification-2026-11-03/README.md)
record the exact scope, resource measurements and test-only corrections. The user
subsequently approved publication and main integration, now completed through
[`m18-interactive-timeline-2026-11-03`](https://github.com/carlozaiana/LevelScope/tree/m18-interactive-timeline-2026-11-03).
[The publication receipt](../../Tests/Reference/results/m18-publication-2026-11-03/README.md)
records the unchanged tested product and preserved checkpoints. The exclusions
and conditional checks below still apply.

Before implementation, record the final impact manifest: explicit required commands,
shared callers, fresh/reused evidence and budget. Start with short targeted Sequoia
audio tests. Build only the standalone and selected test targets, universal arm64/x86_64
with established Intel 10.13 dependencies. Proposed fresh selection: new timeline/
progress tests, measurement/reference checks if the observer changes, diagnostics,
sections/markers, short sessions, M12/M13/M16/M17 selected cases, short export/report
and injected audition. Do not recursively invoke historical milestone runners. Every
required group must be present and pass; preserve failed runs and superseded candidates.

Freeze a combined M14–M18 candidate after Sequoia stabilization, then complete actual
Intel High Sierra execution, cross-platform comparisons/session exchange and selected
deferred video/full-duration checks. **Plan this immediately after M18, no later than
M19**, rather than perpetuating the deferral. Before that block retain the status
**Sequoia tested; High Sierra, video and three-hour checks pending.** This proposal
itself authorizes neither execution nor main integration.

Use one consolidated actual three-hour mono/192 kHz workflow per Mac to close M14's
full-duration mono/RF64 gap while exercising processing/proposals, independent Current
verification, markers, requirements, timeline/progress, saved explanations, sessions,
resources and verified PCM24 RF64 export. Include independent decoder/sample-position
checks beyond 4 GiB. Share fixtures across these checks; do not run separate equivalent
historical long workflows. Short real multichannel cases plus three-hour-sized 6/8-channel
analysis/overview stress data cover channel-density UI bounds; explicitly do not claim
these are full-duration multichannel PCM end-to-end tests. Record original identities
and equivalence rationale for older long multichannel evidence.

Select H.264/ProRes Source/Current seek, offset, stop/resume/restart/EOF and timeline-open
regressions on both Macs, including injected timing and a separately coordinated short
muted physical subset. Since the timeline observes transport and adds UI work, measure
UI/paint pressure with analysis both active and idle. Full codec/colour matrices and
30-minute endurance are conditional on timing/buffering/API changes or concrete failures.
Occasional picture skips remain accepted and recorded; audio discontinuity, accumulated
drift, sustained freezes and identity/lifecycle failures still block completion. Keep
optional tracing confined to its separate executable and full physical 7.1 listening
separately pending.

Proposed budgets (acceptance targets, not measured M18 results):

| Resource / operation | Proposed limit |
| --- | --- |
| Incremental preview publication | At most 5 updates/s; bounded coalesced batches, never a growing UI queue |
| View preparation / paint | p95 <=50 / <=100 ms; at most 4,096 display columns per curve/lane |
| Cancellation acknowledgement | <=1 s; no late adoption |
| Waveform storage, both domains combined | <=128 MiB, adaptive multiresolution bins; no full PCM in RAM |
| Additional peak process RSS | <=256 MiB above the equivalent final-analysis baseline |
| Source analysis with observation vs observation off | <=15% wall-time overhead on the same host/fixture |
| One historical applied record | <=8 MiB; complete structured reasons for supported plan sizes, no silent truncation |
| Session / reports | Existing 128 MiB / 16 MiB bounds; typical core save/parse <=5 s |
| Three-hour combined workflow | <=240 min Sequoia / <=300 min High Sierra; retain measured substep times |

Suggested workload limits: short fixtures 65–120 seconds, extremes/10,000 annotations
through bounded synthetic data, no automatic full cinema rerun. Budget 15–30 minutes
selected automation and 30–60 minutes native UI per Mac, plus build time. Plan roughly
90–150 minutes Sequoia / 150–270 minutes High Sierra for the single long workflow,
subject to measured disk/CPU speed. These estimates use M14's 35.16-minute 90-minute
mono run and M13's 117.3/202.1-minute three-hour composition runs as context, not proof.

Reserve approximately 2–4 GiB for short fixtures, 1–3 GiB per frozen package, and
**40 GiB free per test volume** before the consolidated workflow: three hours of
mono/192 kHz PCM24 is about 5.8 GiB per full audio copy. Check actual staging/verified
copy lifetimes before running. Delete only newly owned disposable intermediates;
never overwrite packages or relocate existing evidence as part of M18. Missing space
or failed I/O must cancel cleanly and preserve the previous valid session.

## Exclusions and decisions for approval

No automatic segmentation, region suggestions, LRA-directed correction, coordinated
neighbour adjustments, new proposal strategy, limiter/compressor, per-component DSP,
dialogue/intelligibility analysis, freeze controls, marker-format expansion, 7.1.4,
Atmos, new export formats or VST3 parity. The next proposal should add bounded region
suggestions and selective LRA correction using this inspection surface; it should not
wait for a comprehensive editor or DSP suite.

M18 combines one automatic initial Source analysis with the logic module's tested
initialisation of identical Current, progressive waveform and M/S inspection, numeric
section drafts and one historical applied explanation. Programme IL/LRA require final
measurements. The visible-passage sample reader is a later milestone. Combined deferred
qualification follows M18.

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

### Implemented interaction

During an ordinary Original or Current analysis the worker publishes only complete
measured M/S windows from the verified programme being processed. A bounded preview
keeps the existing 100 ms cadence and never becomes a measurement snapshot. The
inspector reports the actual window and its frame interval. M has no final gate
membership during analysis; S remains ungated. Silence and missing windows break
the connecting line. Cancellation, stale inputs and completion retire the preview.
Only the accepted result supplies final gate membership, programme values and
first-import Current initialisation. PCM overview storage still has its existing
128 MiB combined Source/Current budget.

The upper audio timecode and lower seconds scale describe the same viewport.
Clicking the graph or pressing Return in the seconds field inspects; a stopped
ruler click or explicit **Seek** moves the existing audition transport without
starting it or switching domains. Fractions are accepted and the displayed
decimal retains the exact inspected frame. A ruler drag pans. **Fit selection**
fits the represented bin or measured interval; **Whole programme** restores the
full view. **Follow** tracks audition position and releases on manual pan.

The comparison divider sits above the fixed Target information line. At minimum
height the profile and save controls share the spare header cell, leaving all
five Target fields usable. The second divider sits under the timeline scrollbar.
Both accept dragging or arrow keys; **Session / tools → Reset layout** restores
their defaults. Preset editing follows actual playback/preparation, load and job
state: Stop or job completion releases it without discarding an existing draft
or changing programme applicability.

All five focused suites passed (2,692 checks), plus the four pinned EBU v05
M/S curve files. See the [development evidence and manual checklist](../../Tests/Reference/results/m18-progress-navigation/README.md).
Gain-lane/section inspection, historical applied explanations and remaining M18
acceptance are still open. The later visible-passage sample reader and agreed
combined qualification deferrals are unchanged.

### Timeline section gain and focus slice

The main timeline's Current section-gain lane uses the same viewport as the
per-channel waveforms and measured M/S. The solid applied curve and dashed manual
draft show the section processing stage before global trim. The section name,
half-open frame range, actual applied dB at the cursor and full-precision stored
setting are inspectable without seeking or switching domains. The Section window
and main lane edit a single runtime draft. **Set draft** or Return stages a
number; only explicit **Apply gains** can change Current. Discard or a genuine
no-op retains accepted measurements. Changing Apply clears Current verification
until independent remeasurement. Annotation edits do not rewrite the envelope;
a draft based on stale geometry cannot be applied.

The lane retains the exact canonical shortened transition endpoints. Visible
geometry keeps at most four time-ordered knots (first, last and both extremes)
per column, with up to 4,096 columns; full-precision inspection is independent of
that display reduction. A click selects a section, a drag pans, and **Fit selection**
fits its exact interval. The lower divider maintains room for the new controls.

After explicit Seek or playback toggle, keyboard focus is held by the timeline
rather than the Seek button that disables during playback. Space toggles audition
there; focused text inputs and modal controls retain their own keyboard behavior.
Playback's protected-copy/hash/decode verification is unchanged; the preparation
message now reads “Verifying audio for playback...”.

See the [focused evidence and manual checklist](../../Tests/Reference/results/m18-gain-lane/README.md).
Saved automatic-proposal explanations and wider M18 acceptance remain outstanding;
the later sample reader and qualification deferrals are unchanged.

### M18 manual section editing and remembered speaker roles

A horizontal drag on either waveform proposes the exact source-frame interval.
Release commits one annotation edit. Drag an internal boundary to move it between
its two neighbours. Existing section names and outer boundaries are preserved.
The numeric frame editor stays available in Sections. Both windows coordinate their audition through
the existing transport, and allow a stopped click on the upper timecode ruler.
The ruler, scrollbar, wheel and Alt-drag provide panning.

Fresh stereo and 5.1 imports keep separate last-confirmed orders on each Mac.
Show the proposed assignment for every file channel. Prefer reliable file
labels, show any disagreement with memory, and offer **Confirm roles and
analyse** and **Edit order**. When neither labels nor memory provide an order,
channel count selects a reviewable standard suggestion. Reopened projects use
their own stored assignments. Separate mono files require the existing explicit
sample-zero alignment confirmation. Apply this memory only to speaker roles.
Physical output assignments continue through audio setup.

### Playback roadmap following M18

**Output selection and reconnection** is a separate playback change. Persist
the selected interface on each Mac using a stable device identity. Selecting
an output connects it when the programme sample rate is known. Provide an
opt-in **Automatically reconnect my last output device** preference. Opening
or reconnecting leaves transport stopped. Explicit Disconnect suppresses
automatic attempts for the rest of the session. Report missing, busy and
rate-incompatible devices distinctly. Keep offline analysis and the output
chooser available in each case. Test rate changes, missing interfaces, busy
devices, explicit disconnect and a new session with the preference on and off.

**Stereo monitoring of 5.1** is a later, explicitly selected listening mode
for either Source or Current. Label the fold-down at the transport and device
panel. Proposed monitor matrix before the common headroom trim is L and Ls
at unity to left, R and Rs at unity to right, and C at −3.01 dB to both.
Default LFE off. **Include LFE in stereo monitoring** adds LFE at −6.02 dB
to each side. Reserve a common 11 dB trim in either setting so coherent
full-scale inputs fit even with LFE enabled. Verify impulse mapping, coherent
sums and both Source and Current paths before enabling it. Keep the matrix
in the listening path after programme processing. Measurement, proposal
verification and 5.1 export consume the full programme.
### M18 main proposal workflow and applied explanations

The main timeline now exposes **Propose**, **Review proposal**, **Apply proposal**,
**Discard proposal**, **Proposal limits…** and **Last applied explanation**. The
existing section-led service still considers global gain first. Apply remains
explicit and is followed by separate **Re-measure Current State**. The outcome
panel leads with readiness, no-op, the limiting requirement/preference, stale
inputs, cancellation or failure, and explains unavailable actions. A rejected
local selection does not imply that its combined candidate was measured.

The existing Current gain lane uses green solid applied gains, amber dashed
manual drafts and cyan dotted automatic previews. The proposed global gain is
shown separately. Section inspection includes exact local plus global gain at
the inspected frame, original intervals, peak attribution, requested/planned
reduction, reserve and uncovered deficit. Its arrows traverse proposed changes.
Rejected and stale previews are labelled, and their Apply action is unavailable.
Complete paginated Review remains available for every correction.

A changing global or section Apply captures its complete historical explanation
before clearing the baseline Current measurement. The session attachment keeps
exact plans/identities, complete requirements/preferences, scalar measurement
statuses and provenance, predictions, original annotations and structured reasons.
A later changing proposal replaces it. Manual edits and requirement/preference
changes retain it with a historical association label. Reopening exposes it in
the main window while Current verification is assessed independently. No private
candidate result is installed as measured Current. The attachment has its own
8 MiB bound and checksum inside the existing 128 MiB session budget. Missing or
bad historical data is reported independently of valid settings and caches.

The extra main-window controls have a larger default window and inspector minimum;
the comparison and waveform dividers remain available. See the
[focused development evidence and manual checklist](../../Tests/Reference/results/m18-proposal-workflow/README.md).
The visible-passage sample reader stays after M18. High Sierra, physical device/video,
cross-Mac and full-duration/RSS qualification retain their agreed deferrals.

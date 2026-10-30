# M18 proposal — interactive audio timeline and progressive analysis

30 October 2026. **Proposed for review; implementation is not authorized.**
Continue from published M17 candidate-04, checkpoint
[`e422f85ccbf7ae11876e6a31be2409937cb8381c`](https://github.com/carlozaiana/LevelScope/commit/e422f85ccbf7ae11876e6a31be2409937cb8381c),
tag [`m17-job-lra-presets-sequoia-2026-10-30`](https://github.com/carlozaiana/LevelScope/tree/m17-job-lra-presets-sequoia-2026-10-30).
The publication/proposal follow-up changes documentation only. Work remains in the
established `codex/logic-module` copy. Main remains qualified M13
`a9be070ab9eb5cb4f591ff9d991d8999c3a016d0` until combined qualification and subsequent
integration approval. Preserve all checkpoint tags, packages and the video worktree.

## Intended outcome and bounded scope

Make the existing measured results and section processing practical to inspect:
**confirm programme → automatic Source analysis with visible progress → inspect the
audio timeline → draft existing section gains → explicit Apply → independently
measure Current → inspect/audition the result**. Playback is always an explicit
separate action. A user can cancel automatic analysis without a restart loop.

Deliver one reusable audio timeline with:

- Source/Current inspection selector; waveform overview, measured momentary and
  short-term curves, marker/section boundaries, selected interval and playhead on
  the same audio-source axis. Display the inspected domain separately from the
  auditioned domain; inspecting Current must not switch playback implicitly.
- Horizontal zoom, pan/scroll, fit programme/selection, exact numerical inspection
  and explicit **Seek to cursor**. Clicking a plot selects/inspects; it does not seek,
  play or open a device. First-version timeline seeking is stopped-only and uses
  the existing seek service. Existing transport capabilities remain available.
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

Build a versioned, bounded multiresolution min/max waveform overview from the actual
programme samples, off the UI/audio callback. Source represents the original confirmed
programme, including a composition's unity sum. Current represents the exact applied
plan through the shared processor, before monitor controls. Do not manufacture a
Current waveform by scaling Source extrema across a time-varying transition.
Draft envelopes are graphics only; they do not acquire a measured waveform or verdict.
Reuse the existing worker-side sample observer boundary where appropriate. Any
production observation must be optional, non-mutating and bounded; no test-capture
buffer may grow with the full PCM duration.

Generate overviews during normal measurement where possible, avoiding a second full
read solely to draw a completed measurement. If a session has valid measurements but
no overview, offer a cancellable **Build waveform** job using the same verified input
and processing; do not invalidate or repeat valid measurements. The first version
keeps overviews as replaceable runtime data, not a mandatory persistent cache or a
new media copy. Reuse an existing verified source only when its identity/lifetime
contract permits it. Missing/changed/unreadable members fail the complete programme;
there are no silent substitute channels.

## Progressive analysis and authority

After successful programme verification and explicit role/alignment/composition
confirmation, schedule Source measurement once for the current input identity when
stopped and no conflicting job is active. On session reopen, reuse compatible validated
measurements; start only when the required Source measurement is absent. Never invent
roles, applicability or job requirements. Expose Cancel and explicit Analyze again;
cancellation/failure suppresses automatic retries for that unchanged identity.
Do not initialize or modify Current processing automatically. Applying a draft keeps
Current UNVERIFIED until the existing independent Current measurement succeeds.

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
distinguishable. No new compliance rule or inferred LRA target is introduced.

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

## Acceptance tests

1. **Final-result parity:** observer on/off and differing batch/block sizes produce
   identical final samples, windows, gates, scalar metrics and identities. Run pinned
   EBU/ITU references if the window observer touches shared measurement code, as
   anticipated. No preview is accepted by any measured-rule/proposal/export gate.
2. **Waveform truth:** actual Source/Current min/max compared with processed samples
   for mono/stereo/5.1 and supported 7.1 global trim, composition, assembled members,
   all rates, LFE, positive/negative gain, transitions, silence, over-full-scale finite
   values, odd EOF and narrow impulses. Include anti-phase stereo to prove there is
   no implicit downmix. Zoom reduction preserves extrema and source-frame locations.
3. **Progress/lifecycle:** deterministic cancellation mid-verification and mid-analysis,
   missing/changed members, rapid reopen/source changes, stale batches, no retry loop,
   closing during work, bounded queues and allocation failures. Cached reopen performs
   no unnecessary measurement. Missing waveforms cannot invalidate sound measurements.
4. **Navigation:** frame-zero/end, fractional-time floor convention, long coordinates,
   exact window selection, no phantom edge markers, EOF padding and viewport clipping.
   Explicit seeks through the existing service preserve audio-only continuity, duration,
   stop/restart/EOF and pre-monitor sample parity. Mere clicks/zoom/selection do nothing
   to transport, measurements, proposal validity or processing.
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

Recommended decisions for review are: automatic **Source-only** analysis with explicit
retry after cancellation; numeric existing-section edits rather than freehand automation;
progressive waveform/M/S but final-only programme IL/LRA; a single last-applied historical
explanation record; and combined deferred qualification immediately after M18. These
are concrete proposed boundaries, not yet implemented features.

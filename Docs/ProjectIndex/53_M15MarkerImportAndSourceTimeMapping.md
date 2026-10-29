# M15 proposal — marker import and explicit audio-source time mapping

Status, 29 October 2026: **Approved for implementation and impact-selected Sequoia qualification.**
The user approved both explicit DAW readers and the clarifications below. Nuendo 15.0.30’s six-point CSV and 25 fps setting are confirmed. A second actual export confirms
its two-second export-coordinate origin. A Pro Tools 10.3.10 export is now supplied;
its parser and focused UI checks passed on candidate-05. The user confirmed 24 fps,
project/audio/trailer start at zero, and the first two visible marker positions. Actual interoperability is recorded separately.

## Starting point and qualification boundary

Continue in `/Users/fuzz/Developer/LevelScope` on `codex/logic-module` from the
published M14 development checkpoint
[`m14-mono-rf64-sequoia-2026-10-29`](https://github.com/carlozaiana/LevelScope/tree/m14-mono-rf64-sequoia-2026-10-29),
implementation commit `3d98284c2c0b7de92cfee70cb0007d43d1f04d1f`, plus its
publication-only receipt `be4bf23c8bb9f661ae1ede30f3b6e79a067aea37`.
M14 candidate-02 source manifest:
`9dd22a598281e83431039f210202a8a24fe49517294cc6bb8d83032c4975e013`.
Read [50](50_M14MonoProcessingAndVerifiedExport.md),
[51](51_M14ImpactQualificationPlan.md) and [52](52_M14QualificationAndHandover.md).

Main remains qualified M13 `a9be070ab9eb5cb4f591ff9d991d8999c3a016d0`.
Preserve every checkpoint, frozen package, historical receipt, original media,
installation, the separate video worktree, Pro40 and system audio configuration.
No publication or main integration is part of this proposal approval.

Retain the agreed deferrals: actual Golem/High Sierra execution and UI, cross-Mac
exchange, all video loading/decoding/playback, and full three-hour qualification.
Continue the established universal/Intel cross-build and deployment inspection;
this does not establish High Sierra runtime qualification. Combined qualification
is planned after M18, no later than M19, on one frozen candidate with an explicit
impact-selected manifest. Each intervening milestone still needs separate approval.
Earlier testing needs a concrete API/dependency, arithmetic, buffering/clock/latency
change or observed failure, explained before revising these deferrals.

## Recommended outcome

Import editorial point markers, or paste a list of positions, preview their exact
mapping to audio-source samples, and explicitly replace the annotation section map.
This removes repetitive manual boundary entry while retaining separate annotation,
applied-processing and measured-result identities. Markers never certify scenes,
select gains, start a proposal, remeasure, seek, play, open a device or move video.

The cinema test illustrates the usefulness of explicit mapping: the user supplied
video frames at 25 fps with frame zero at audio sample zero. Those six boundaries
created seven sections. Its proposal used four changed sections and was rejected
because affected envelope duration was 89.529560%, above 20%; the eight-changed-
section limit was not reached. The user wants parameter review later. **M15 does
not change the 3 dB cap, eight changed sections, 20% duration, reserve, pass budget,
manual +/-12 dB range, transitions, profile limits or already-pass behavior.**
The existing 10,000-annotation capacity is a different limit from eight changed
sections in a single automatic proposal.

## Approved DAW routes and export instructions

- **Nuendo CSV:** File > Export > CSV Marker, select one point-marker track,
  UTF-8 encoding (optional BOM), comma/semicolon/tab separator, headings enabled,
  and position/name attributes. The importer requires explicit position/name column
  selection and shows every unused column. Select a marker-type column when present;
  cycle/range/duration/end fields cannot be silently treated as point markers.
- **Pro Tools:** File > Export > Session Info as Text, with **Include Markers**
  selected. A dedicated reader locates the marker table and validates its structure;
  this is not generic CSV import. Unrelated session/track text is not imported.
  Supported marker fields and encoding are documented with the reader fixtures.
  Unsupported marker structures produce a clear error.

For each actual DAW export, record exact version, export choices, encoding, separator,
column names, time format, programme origin and positions visible in the DAW. Verify
the marker field’s actual units independently of a general session-export setting.
The user supplied `/Users/fuzz/N15/6NuendoMarkers.csv` from Nuendo 15.0.30:
UTF-8 BOM/comma/headings, six `Timecode In` positions, empty `Timecode Out`,
Description and Length. The user confirmed default UTF-8/comma export at 25 fps,
project start `-1 23:59:58`, programme start `00:00:00:000`, and first two visible
markers at `00:00:23:05` / `00:00:29:04`. CSV labels are exactly two seconds later.
An explicit export-coordinate programme origin `00:00:02:00` matches those two
observations. The second user export `TimeOffsetMarkers.csv` adds a named marker
at displayed programme zero; it exports at `00:00:02:00`, confirming this mapping
with the other six rows unchanged. The exact Time Offset control value is not
separately reported. Do not call this an established Nuendo bug or silently compensate.
The user supplied `/Users/fuzz/N15/TUTTI_UCCIDONO_V2.txt` from Pro Tools 10.3.10.
It contains 45 point rows under `M A R K E R S  L I S T I N G`, after extensive
track/file listings. Headers declare 48 kHz, 24 Frame and session start zero.
`LOCATION` is NDF display timecode; `TIME REFERENCE` contains more precise sample
positions and every row declares `UNITS=Samples`. Explicit audio-sample units
select that reference, requiring one matching sample-rate header and consistent
row units. Other explicitly selected units use LOCATION, with a visible subframe
precision warning. No automatic field/origin inference or resampling. The first
LOCATION labels are 12:05 and 33:13; sample references are 587,888 and 1,610,400.
The user confirmed project and actual audio start at 00:00:00:00, with a trailer
there, and first visible markers at 00:00:12:05 / 00:00:33:13. Origin zero is correct;
“Start Audio” is an internal marker and must not shift the import.
The minimal verbatim fixture excludes unrelated project/track paths; the full
original is checked locally. See the fixture record for exact hashes and scope.
Synthetic fixtures alone qualify parser behavior only; missing independent mapping
confirmation leaves actual programme interoperability pending without stopping work. No native project or complete marker-track reconstruction
is included. [Avid reference guide](https://resources.avid.com/SupportFiles/attach/Pro_Tools/12.0/Pro_Tools_12_0_Reference_Guide.pdf)
documents Session Info as Text / Include Markers; it does not establish the user’s
current version or its concrete exported marker-field units.

## Bounded inputs and UI

1. Add **Import markers…** and **Paste positions…** to the existing section editor.
   A confirmed valid audio programme is required; an attached picture is optional.
   Use the existing supported annotation layouts/rates without expanding processing
   capability: mono, stereo, 5.1 and existing 7.1 annotation behavior. M15 does not
   enable section processing or export for otherwise unsupported layouts.
2. Accept UTF-8 text, optional BOM, and CSV/TSV with an explicitly chosen comma,
   semicolon or tab delimiter. Allow a header row, a required position column and
   an optional label column; support quoted delimiters/escaped quotes and LF/CRLF.
   Plain paste mode accepts one position per line or a comma-separated integer
   list, including the six supplied cinema numbers. Do not infer units from size.
3. Support **point markers only**. Range/cycle-marker interpretation, multiple
   marker tracks, automatic column guessing, general spreadsheet import, cue sheets,
   embedded WAV markers, EDL/XML/AAF/MIDI/DAW-project parsing remain outside M15.
   Additional CSV columns are explicitly listed as unused in preview; selecting
   point positions is not a claim to preserve ADR, duration or track semantics.
4. Show original row/position/label, selected units and origin, exact rational time,
   resulting audio frame, any quantization error, and the resulting complete section
   table in preview; the committed map updates the existing source-time timeline. Show added/removed boundaries and
   changed names before confirmation. These are proposed annotations until commit.
5. First version offers **Replace annotations** only, with Preview, Confirm and
   Cancel. No implicit merge or gain remapping. Require stopped transport and no
   active session/measurement/proposal/export job at confirmation; never stop or
   cancel another operation implicitly. Confirmation is one atomic document edit.

A practical CSV example is `position,label`, with rows such as
`42397,Peak region A`; the user explicitly selects video-frame counts, 25/1 fps
and origin zero. A label names the section starting at that marker. With no label,
use deterministic section names. The initial interval has its own explicit name.

Steinberg documents CSV marker export with selectable attributes, encoding and
separator: [Nuendo CSV export](https://www.steinberg.help/r/nuendo/15.0/en/cubase_nuendo/topics/markers/markers_exporting_markers_as_csv_file_t.html),
[export options](https://www.steinberg.help/r/nuendo/15.0/en/cubase_nuendo/topics/markers/markers_export_options_r.html).
Use column mapping rather than a hard-coded localized header. One representative
user-exported point-marker file is needed before claiming compatibility with that
particular Nuendo export setup. Its absence must not turn synthetic CSV evidence
into a Nuendo compatibility claim; canonical CSV/paste support can be qualified
independently. No exact version of the earlier Nuendo cinema measurement is inferred.
CSV quoting follows the bounded rules described above, informed by
[RFC 4180](https://www.rfc-editor.org/rfc/rfc4180); semicolon/tab are explicit extensions.

## Nuendo 15 manual review and additional-format scope

The 29 October manual review identifies these interchange routes. A general list
of formats containing markers must not be read as a claim that Nuendo exports each
one through a marker-export command.

| Route | Documented export capability | M15 decision / later work |
| --- | --- | --- |
| [CSV Marker](https://www.steinberg.help/r/nuendo/15.0/en/cubase_nuendo/topics/markers/markers_exporting_markers_as_csv_file_t.html) | Selected marker tracks, attributes, encoding and delimiter | Current explicit point-reader; actual origin validation in progress |
| [MIDI](https://www.steinberg.help/r/nuendo/15.0/en/cubase_nuendo/topics/markers/markers_exporting_markers_via_midi_c.html) | Point/cycle starts and track assignment | Separate adapter needs tick/tempo or SMPTE timing validation; cycle ends cannot be assumed preserved |
| [XML track archive](https://www.steinberg.help/r/nuendo/15.0/en/cubase_nuendo/topics/markers/markers_exporting_markers_as_part_of_a_track_archive_c.html) | Marker IDs, attributes, track assignment, cycle start/end | Candidate later adapter; retain point-only boundary until range semantics approved |
| [WAV marker chunk](https://www.steinberg.help/r/nuendo/15.0/en/cubase_nuendo/topics/export_audio_mixdown/export_audio_mixdown_options_r.html) | Audio Mixdown can embed markers from selected tracks | Separate metadata reader with export-range/sample-origin and container checks; no need to decode media to import annotations |
| [Cue sheet CSV/TXT](https://www.steinberg.help/r/nuendo/15.0/en/cubase_nuendo/topics/exchanging_files_with_other_applications/exchanging_files_with_other_applications_cue_sheet_c.html) | Project information and marker-track list among other sections | A different sectioned document from CSV Marker; require its own tested reader |
| [CMX3600 EDL](https://www.steinberg.help/r/nuendo/15.0/en/cubase_nuendo/topics/markers/markers_scene_detection_with_edl_files_c.html) | This page documents import from an editing system | Remains a later bounded EDL proposal; not evidence of Nuendo marker-track export |

[Track archives use XML](https://www.steinberg.help/r/nuendo/15.0/en/cubase_nuendo/topics/track_handling/track_handling_exporting_tracks_as_track_archive_c.html),
but project-specific settings such as tempo are not generally included.
[Project Setup](https://www.steinberg.help/r/nuendo/15.0/en/cubase_nuendo/topics/project_handling/project_handling_project_setup_dialog_r.html)
distinguishes Project Start Time and the Time Offset applied to ruler/position
displays. The consulted CSV pages do not define the export field's origin relative
to those controls. The two-second interpretation is verified for this setup by the supplied
programme-zero marker and visible-position observations, not asserted as a vendor
guarantee for other settings or versions.

The user confirmed: finish M15 with CSV/Pro Tools text; keep the other formats
pending while assessing which are important. No additional format implementation
is approved. Any later adapter needs real exported fixtures, explicit time mapping
and the same preview / annotation-only transaction, under a separate bounded plan.

## Exact mapping and edge behavior

Recommended supported position representations:

- Whole audio sample-frame indices at the verified source rate.
- Finite exact seconds (up to nine decimal fractional digits, or numerator/denominator), parsed as integers/rationals rather than
  binary floating-point approximations.
- Whole video-frame counts at an explicitly selected rational rate: 24/1, 25/1,
  30/1, 24000/1001 or 30000/1001 fps.
- Non-drop `HH:MM:SS:FF` labels at those rates, with an explicit nominal label
  base (24, 25 or 30) and a programme-start label. No rate inference from labels.

For frame labels/counts, let F be the parsed frame index and F0 the explicitly
confirmed index corresponding to audio sample zero. Convert using checked integer
arithmetic: `sample = floor((F - F0) * rate_denominator * sample_rate / rate_numerator)`.
Seconds use the same rational subtraction and floor rule. Show when the rational
position is not sample-integral: the boundary is less than one sample earlier,
never silently described as sample-exact. For exact mappings, show zero error.
Audio-frame mode subtracts the confirmed integer origin without quantization.
Negative intermediate times and positions beyond real EOF are rejected before
rounding; do not clamp, wrap, trim or silently repair them.

This is editorial time mapping, not a decoded video-frame lookup. It does not
infer VFR frame positions, inspect a video file, copy its offset, or change the
existing audio/video clock. Reference attachment and signed offset must leave
all imported audio positions unchanged. No changes to presentation/transport
code are planned; reuse checked rational primitives without changing their behavior.

Drop-frame labels, 24-hour wrap inference, pull-up/down, retiming and mixed-rate
lists are deferred. Reject unsupported/ambiguous forms, zero denominators,
non-finite values, overflow and invalid frame fields with row-specific errors.

Sort valid points by mapped audio frame and disclose reordered rows. Identical
positions with identical labels may coalesce, visibly; conflicting labels at the
same mapped frame block confirmation. The same rule applies when two different
fractional positions quantize to one sample. Preserve all row provenance for any
coalesced points. Positions at zero and exact EOF are existing endpoints, not
zero-length sections; preview labels them accordingly. An EOF label is retained
as endpoint provenance, not assigned to a nonexistent following section.

Construct half-open `[begin,end)` sections covering `[0,source_frames)` exactly,
with no gaps/overlap/empty intervals. Reuse IDs for unchanged exact intervals;
allocate new IDs monotonically and never reuse removed IDs. An identical import
(map, labels and provenance) is a true no-op. A provenance-only update changes
only document metadata; a name-only change does not invalidate audio or geometry.

## Architecture and identity

Add a small pure, versioned importer/converter beside the existing section model:
immutable parsed input -> resolved mapping -> preview/diff -> transactional commit.
It must be testable without audio decoding, a window, an audio device or a video
reader. The UI supplies confirmed source identity, rate, frame count, explicit
mapping and current annotation map. Large text parsing stays off the UI thread.

A preview binds to its input bytes/hash, mapping version/options, source content,
ordered roles, programme range, annotation revision and runtime session/job token.
Changed source/roles/range, competing annotation edits, changed file bytes or a new
session invalidate pending confirmation. Late completion cannot affect another
session. Profile/video/monitor changes never alter the mapping or audio identity.

Committing annotations increments the persistent document revision, not the audio
revision. Source/Current analyses, compliance verdicts, processing-plan identity,
audition samples and rendered-audio identity remain valid when audio is unchanged.
Geometry-dependent section diagnostics and proposal attribution are recomputed;
manual drafts and section-led proposals bound to older geometry become visibly
stale and cannot Apply. Retain their explanation until explicit discard/recreation.
Global-only proposal validity follows its existing audio/profile binding.

Existing applied section gains and their captured annotation provenance keep their
original absolute intervals and transitions. Import must never remap or remove
these applied gains. The existing warning that current annotations differ from
applied processing remains visible. This includes import while a non-unity section
plan already exists. No proposal generation or processing Apply is part of import.

A report request captured before an annotation change is stale as a document/report
request; it must not publish obsolete labels unnoticed. An unchanged rendered audio
artifact remains reusable under existing audio identity checks. Refresh its report
and require existing explicit publication checks, without rerendering unchanged audio.

## Session and provenance design

Persist one latest import receipt with importer/schema version, marker-file basename
and SHA-256 (or pasted-text hash), original positions/row references, explicit unit,
rate/origin/rounding convention, normalized sample boundaries, labels/section IDs
and confirmation binding. Retain enough data to explain the mapping without the
original marker file. Do not embed unselected ADR/dialogue columns or the whole raw
CSV. Later manual edits retain the import receipt as historical provenance and
clearly show that the current annotations differ. An import-history browser is deferred.

Use session schema 5 for this explicit provenance, preserving readers for
schemas 1–4 and the existing replaceable-cache policy. Saving with M15 requires an
M15 reader; warn before overwriting an older session and offer Save As. Opening
never rewrites the original. Prove that the frozen M14 reader rejects schema 5
clearly rather than silently losing provenance on resave. Do not change the audio
programme/assembly meaning or measurement algorithm/cache format unnecessarily.

The current build identity hashes all Source files: an M15 source change normally
changes that identity even if measurement arithmetic is untouched. Preserve the
conservative cache policy. Old-build sessions recover settings/media/annotations,
discard incompatible caches with explanation and require remeasurement. Within
the same compatible build, importing markers does not force audio remeasurement.
Reopen remains stopped/device-closed/routes-muted; no marker file is needed to
recover committed annotations. Reports include versioned mapping provenance with
labels escaped as text, separate from measured-rule authority.

## Acceptance tests and impact-selected Sequoia phase

1. **Exact cinema boundary fixture, no film rerun:** inputs
   `42397,42575,45099,47283,130243,130547`, 25/1 fps and origin zero at 48 kHz map
   to `81402240,81744000,86590080,90783360,250066560,250650240`. They produce
   seven sections ending at `269875200`, exactly matching the completed local test.
   This is a mapping fixture, not another full-film measurement or video test.
2. **Mapping oracle:** independently calculated audio-frame, seconds, rational-rate
   and non-drop-label cases, including nonzero programme origins, exact samples,
   fractional-sample floor/error reporting, endpoints, last sample, overflow,
   malformed/negative/out-of-range input and explicit rejection of drop-frame/VFR
   assumptions. Test 44.1/48/96/192 kHz; no floating-point tolerance can move a split.
3. **Parsing/preview:** headers, selected delimiters, quoting, escaped quotes,
   LF/CRLF/BOM, Unicode, empty/malformed rows, unsupported encoding, duplicate and
   conflicting labels, sort disclosure, limits and pasted lists. No silent row
   deletion/truncation or partial commit. Escaped labels cannot inject report HTML.
4. **Transaction/lifecycle:** cancel/error leaves the prior map and session intact;
   stale source, annotation revision, input bytes or runtime rejects confirmation;
   a failed commit changes nothing; closing/new/open during parsing cannot adopt a
   late result. Repeated identical import is a no-op. Save captures one revision;
   edits made while saving stay unsaved afterward.
5. **Audio invariance:** short mono/stereo/5.1 cases before/after import retain
   exact Source/Current metrics, evaluations and processing identity. Include a
   non-unity applied section plan, all-channel/LFE sample captures, audio-only
   seek/EOF/continuity and unchanged M13 limits/no-change/stale-Apply safeguards.
   Check existing 7.1 annotation behavior without enabling unsupported processing.
6. **Sessions/reports/export:** schema migration, Save As, exact annotation/provenance
   roundtrip, same-build cache reuse, old-build cache discard, bad cache recovery,
   missing original marker file, write/space/cancel failures and old-reader rejection.
   Short verified export demonstrates unchanged PCM/TPDF samples and audio identity;
   report refresh uses the new annotations and strict existing Current/artifact gates.
7. **Silent Sequoia UI:** file/paste, units/rate/origin, readable row errors and
   quantization, preview/cancel/confirm, stale drafts, imported intervals on timeline,
   save/reopen and no implicit transport/device action. No video media is loaded.
8. **Bounded resources:** up to 8 MiB text, 10,000 input rows and 10,000 resulting
   sections, with the existing 128-character section-name bound. No silent truncation.
   Keep the existing 4 MiB session-core and 128 MiB total-session limits; preflight
   serialized annotations/provenance before commit and reject an over-budget import.
   Target <=64 MiB added importer RSS, <=2 s parse/map for the maximum representative
   fixture, <=100 ms UI heartbeat p95 and <=1 s cancellation. Use virtualized preview;
   avoid repeated per-marker whole-map copying and hashing. Record actual numbers.

Before implementation, enumerate the required fresh runner cases and shared callers.
Run one explicit marker-import test target plus selected section/diagnostic, session,
M12/M13 identity, short audition and report/export regressions. Build only named
standalone/test targets, including the Intel slice; never default-all/VST3/install.
A missing required receipt blocks completion. Estimate 15–30 minutes automated
short checks and 15–25 minutes silent UI, excluding builds/fixes/operator delays;
allow 2 GiB disposable fixture space plus separately preflighted build headroom.
No new sustained audio workflow is required for a metadata-only import feature.
Large frame-index/marker stress fixtures do not qualify three-hour audio execution.

Reuse M14 candidate-02 measurement/reference, mono/RF64 and 90-minute evidence only
with its exact identity and unchanged dependency rationale; it is not fresh M15
qualification. Pinned EBU/ITU reruns are conditional on actual changes to shared
measurement/input/sample processing or uncertain dependencies, not automatic for
annotation-only code. Inspect session/build callers as well as directly edited files.
After fixes, rerun affected checks and widen coverage only for concrete findings.

Actual High Sierra/cross-Mac exchange, all video execution and full three-hour checks
remain pending for the combined window. Accepted occasional picture skips retain
recorded counts/timing and the existing audio-continuity, drift, freeze, identity,
memory and lifecycle blockers. No physical-device opening or audible test is needed
for M15; any newly justified physical test requires separate coordination.

## Exclusions and approval decisions

No proposal-parameter editor or changed limits, automatic segmentation, peak-event
inference, gain selection/application, neighbouring-section optimization, limiter,
compression, dialogue analysis, component DSP, new audio import/export formats,
7.1.4/Atmos, VST3 parity or new video feature. Bounded EDL remains a later proposal;
marker import does not become a prerequisite for useful future dynamics-led methods.

Approved defaults: point-marker Nuendo CSV / Pro Tools text / paste; explicit units/rate/origin;
non-drop forms only; visible floor-to-sample conversion; replace-with-preview;
schema-5 provenance and Save As; impact-selected Sequoia checks with existing deferrals.
A real Nuendo point-marker export is useful for a specific interoperability claim,
but is not required to review or approve this canonical importer design.

On later Sequoia-phase completion, report candidate identity, fresh/reused evidence,
resource results and pending checks as **“Sequoia tested; High Sierra, video and
three-hour checks pending.”** Stop for review before checkpoint publication, main
integration or starting the following milestone. M15 implementation is approved;
additional formats remain pending a separate usefulness/scope review.

## Subsequent roadmap — separately approved milestones

1. Bounded proposal-preference review and configurable controls. Consider section
   count, affected envelope duration, additional and cumulative gains, neighbouring
   corrections, transitions and resources together. The cinema test changed four
   sections covering 89.53%; its 20% duration gate blocked the attempt, not the
   eight-section gate. Include marker-density sensitivity: splitting the same passage
   more finely must not automatically count as proportionally greater audible
   alteration. Preserve the current policy version; configurable policies need their
   own explicit version, identity, tests and approval. No limits change in M15.
2. Standalone true-peak limiter qualification, before depending on it in proposals.
3. Basic measured-audio region suggestions before the broader enhanced section-led
   strategy. Suggest editable regions, not narrative scenes. Protect user-confirmed
   or protected boundaries and require review/adoption. Dynamics-led proposals must
   remain usable without markers or sections. Further DSP/strategies and bounded EDL
   import remain separate approvals.

Implementation uses [54 — impact plan](54_M15ImpactQualificationPlan.md). Completion
requires reporting fresh/reused/pending evidence, exact candidate identity and real
DAW versions/settings. Stop for review before publication, main integration or the
next milestone.

Concrete reader dialects, export instructions and the pending actual-version matrix:
[Tests/Markers](../../Tests/Markers/README.md).

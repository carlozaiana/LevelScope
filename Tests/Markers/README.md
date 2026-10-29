# M15 point-marker import fixtures and interoperability record

The checked-in examples include synthetic parser fixtures and two unchanged actual
Nuendo 15.0.30 CSV exports plus verbatim timing headers and the marker table from
the user’s Pro Tools 10.3.10 export. Record exact DAW
version/settings and cross-check visible positions and actual exported field units
before claiming interoperability; a parser pass is not a general DAW-version claim.

## Nuendo CSV subset

File > Export > CSV Marker: export one point-marker track, choose UTF-8 (optional
BOM), comma/semicolon/tab delimiter, headings enabled, position and name attributes.
Select the exported position/name columns explicitly in LevelScope (one-based).
A missing name can use column 0. If a Type/Marker Type column exists, select it too.
Supported type values: Point, Marker, Position (case insensitive). End/End Position/Timecode Out
must be empty or equal Position; Duration/Length must be empty or an explicit zero.
Other types or nonzero/ambiguous ranges block the import. Unknown columns are listed
as unused; the user must confirm the point-only subset. Localized column headings
can be selected by index; unknown marker-type dialects require an adapter update,
not silent reinterpretation. CSV supports quoted separators and doubled quotes.
Marker names are single-line text of at most 128 characters. Empty rows within the
marker table and malformed quoting/field counts are rejected. LF and CRLF supported.

## Pro Tools text subset

File > Export > Session Info as Text, **Include Markers** selected. The dedicated
reader recognizes MARKERS / M A R K E R S and MARKERS LISTING (including its spaced
form), with a tab-separated header containing `#` (or `Number`), `LOCATION` (or
`Position`) and `NAME`. Column order is not fixed. Fixed-width padding within the
tab fields is retained as original text. Empty COMMENTS fields are supported.

Choose the field through the explicit position-units selection:

- **Audio sample frames:** uses `TIME REFERENCE`, requiring `UNITS=Samples` on
  every row and one `SAMPLE RATE` header exactly matching the programme rate.
  Enter programme origin in that sample coordinate system. No resampling and no
  fallback to LOCATION. Missing/duplicate headers, mixed units, fractional or
  overflowing sample counts reject the complete import.
- **Other units:** uses `LOCATION` only, with the explicitly selected rate/origin.
  For this actual export choose NDF timecode and 24/1. The displayed timecode
  omits subframe precision; the preview discloses this and the unused reference
  field. `UNITS` describes TIME REFERENCE, not LOCATION. No session-header inference.

Other marker columns are disclosed as unused. A type column, if present, must
describe point markers. End/length/duration tables and ambiguous section/row
structures are rejected. Surrounding track/file listings are not imported as
markers. Export UTF-8 or explicitly convert to UTF-8 while retaining the original;
UTF-16/legacy encodings remain unsupported. Musical bars/beats, mixed reference
units, ranges and project/memory-location reconstruction are outside this reader.

Confirm the selected field against DAW-visible positions and the actual programme
origin. A marker named “Start Audio” does not automatically establish sample zero.

## Mapping and migration

Select audio-frame integers, exact seconds (up to nine decimal fractional digits,
or a numerator/denominator), video-frame integers or NDF HH:MM:SS:FF. Select exact
24, 25, 30, 24000/1001 or 30000/1001 rate and the origin representing audio sample
zero. For example, timecode origin `01:00:00:00` is subtracted before conversion.
Mapping floors to an audio sample and displays the exact earlier adjustment.
Drop-frame, inferred wrap, VFR lookup and mixed rate are rejected/deferred.

Preview original rows, complete section map and complete changes/disclosures.
Confirm replaces annotations only. Applied processing retains its original absolute
intervals; compatible measurements stay valid. Schema-5 sessions retain mapping
provenance without requiring the marker file. Saving requires an M15 reader; Save As
retains an original M14-readable session. Cached results from a different build are
conservatively discarded; media/settings/provenance remain recoverable.

## Actual interoperability status

| DAW | Exact version / export settings | Actual units / positions verified | Status |
| --- | --- | --- | --- |
| Nuendo | 15.0.30, default UTF-8/comma export, 25 fps (user reported); BOM/headings observed | First two ruler positions are two seconds before exported labels; explicit export-coordinate origin `00:00:02:00` reproduces both | Actual CSV structure and programme-zero mapping confirmed for this export setup |
| Pro Tools | 10.3.10 (user reported), Session Info with most other information; UTF-8/tab observed | Header declares 48 kHz/24 Frame/start 00:00:00:00; 45 LOCATION labels and Samples references | Candidate-05 actual-file parser/UI passed; user confirmed 24 fps, zero programme origin and first two visible positions |

The actual six-point Nuendo fixture has SHA-256
`24941511d781d103f801a348b022001f3eda2b33887d694b37a87393f1b835b7`.
Position is column 2 (Timecode In); name is column 4 (Description, empty in this
example). Timecode Out and Length are empty. This confirms the observed export
structure; it does not by itself establish frame rate or programme origin.

On 29 October the user confirmed 25 fps and a project start of `-1 23:59:58`,
with programme audio/video zero at `00:00:00:000`. The first two ruler positions
are `00:00:23:05` and `00:00:29:04`, versus exported `00:00:25:05` and
`00:00:31:04`. Subtracting an explicit export-coordinate origin of `00:00:02:00`
(50 frames) reproduces both exactly. At 48 kHz those positions are audio frames
1,113,600 and 1,399,680. The same offset is applied to the other four rows. A second export now confirms
programme zero directly; no universal Nuendo offset rule is established.

The test retains a zero-origin raw-label case and separately tests the explicit
two-second mapping, unchanged original labels and provenance round trip. Earlier
candidate records using origin zero are raw-coordinate tests, not evidence that
the markers align with the displayed programme. No input file is rewritten and
no automatic vendor-specific correction is applied.

The second original `/Users/fuzz/N15/TimeOffsetMarkers.csv` is preserved verbatim
as `nuendo-15.0.30-programme-zero.csv` (368 bytes, UTF-8 BOM/comma), SHA-256
`f4709037d4cbe30436be995b4fa3cf5cd9ce1795575ee97594fca60216e09cfd`.
The user reports the same project/settings plus a marker at displayed programme
`00:00:00:00`. Its name is `AudioVideoStart at 00:00:00`; its exported position is
`00:00:02:00`. All six earlier rows are unchanged. This independently confirms
that explicit export-coordinate origin for this setup: the anchor maps to audio
sample zero and introduces no empty section. Both original files remain unchanged.
The exact numeric Project Setup Time Offset control was not supplied separately;
its general relationship to CSV export is not claimed. The observed anchor and
first two visible positions establish the mapping without inferring that control.

The actual-file regression checks named anchor, original six rows, exact mapped
samples, endpoint behavior and receipt round trip. Candidate-01/02 records remain
historical; the added anchor tests are frozen separately. Pro Tools evidence follows. See ProjectIndex 53 for the manual review:
additional formats remain pending a usefulness review, with no new adapters in M15.

## Actual Pro Tools 10.3.10 export

Original `/Users/fuzz/N15/TUTTI_UCCIDONO_V2.txt`: 672,421 bytes, 6,391 lines,
UTF-8 without BOM, SHA-256
`10aecddf80f7dcb668315f68cfd4782a10fe0b218a3c2d487529267fbc60dc25`.
The user reports version 10.3.10 and exporting most other information; the marker
list is at the end. The exact export-dialog checkbox inventory was not supplied.
The file itself proves the marker table and surrounding listings are present.

The minimal fixture `pro-tools-10.3.10-markers.txt` contains verbatim original lines
2–5 (sample rate, bit depth, session-start TC, TC format), a separating LF, and
lines 6,345–6,391 (heading, header, all 45 markers). It excludes the unrelated
project name, media paths, track and clip listings. It is 3,962 bytes, SHA-256
`f9f5bf49f5ddda6707e9031b587d78a383797fb84eb2af5b64b29c2ba8aaa460`.
The full original is tested locally against the same selected marker content; it
is not copied to the repository or frozen package. No original bytes are modified.

First two LOCATION labels: `00:00:12:05`, `00:00:33:13`. At explicit 24 fps,
48 kHz and origin zero these map to samples 586,000 and 1,610,000. Their actual
TIME REFERENCE / Samples values are 587,888 and 1,610,400, retaining 1,888 and
400 subframe samples respectively. Last sample reference is 41,343,984. All 45
references lie within their corresponding displayed frame. The user independently confirmed project/audio/trailer start at zero, 24 fps and
those first two displayed positions. Thus origin zero reproduces programme
coordinates for this supplied setup; “Start Audio” is an internal marker.

Regression coverage exercises exact references, explicit origin subtraction,
timecode-only precision disclosure, full/trimmed-table equality, malformed
reference/rate/units/ambiguous tables, original text and receipt round trips.
The independent user confirmation is recorded separately from parser results.
No offset is inferred from the “Start Audio” name.

Candidate-05 ran the complete original and minimal fixture in its fresh required
manifest: 535 marker checks in ten groups, with all eight short regression commands
passing. Silent native UI also tested both fields, cancellation, explicit commit,
Save As and safe reopening. Together with the user confirmation this establishes the bounded reader and
programme mapping for this supplied setup, not every Pro Tools export dialect.

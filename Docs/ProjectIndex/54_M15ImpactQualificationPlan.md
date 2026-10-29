# M15 impact-selected qualification plan

Approved implementation begins 29 October 2026, from logic receipt
`be4bf23c8bb9f661ae1ede30f3b6e79a067aea37` / published M14 candidate-02.
Main remains M13 `a9be070ab9eb5cb4f591ff9d991d8999c3a016d0`.

## Affected behavior and required fresh evidence

| Area / shared callers | Required checks |
| --- | --- |
| Pure point readers / exact mapping | Nuendo CSV, dedicated Pro Tools marker table, paste; encoding/quoting/unsupported structures; independently calculated rate/origin/quantization/overflow cases; six cinema boundaries |
| Annotation replacement / identity | Atomic preview/cancel/commit/no-op; stale source/map/runtime/bytes; ID preservation; non-unity applied gains unchanged; geometry drafts/proposals stale; measured snapshots/compliance/audio revisions unchanged |
| Session schema 5 / reports | Legacy 1–4 recovery, original untouched, Save As, frozen M14 rejection, receipt validation, same-build cache reuse/conservative old-build rejection, late completion and failed saves; escaped report provenance, refreshed report with unchanged rendered artifact |
| Shared section/proposal/export callers | Named short section, session, M12/M13 identity, audition sample capture and PCM export regressions; no recursive historical runner |
| UI / resource lifecycle | Silent import/paste/mapping/preview/diff/cancel/confirm/stale/save/reopen; device closed; virtualized 10,000-row view; <=8 MiB input, <=10,000 rows/sections, <=128-character names; <=4 MiB core/128 MiB session; <=64 MiB incremental importer RSS, <=2 s representative maximum parse/map, <=100 ms heartbeat p95, <=1 s cancellation |
| Actual DAW interoperability | User Nuendo CSV and Pro Tools Include Markers text exports, exact version/settings and actual marker-field units verified against visible DAW positions; Nuendo 15.0.30 actual exports and 25 fps/programme-zero mapping now verified; Pro Tools 10.3.10 actual parser/UI checks passed with user-confirmed 24 fps, zero programme origin and first two visible positions |
| Build compatibility | Explicit standalone/test targets, arm64+x86_64 slices; established Intel 10.13 deployment and dependencies. No VST3/default-all/install |

The fresh runner will enumerate required cases and reject missing receipts. Run
short targeted checks during development, then freeze one candidate and run the
selected manifest plus UI against that identity. Fixes require affected reruns.
No video media is opened by these tests. Tests with hidden historical recursion
or optional video paths must not be used without explicit selection.

## Reuse and deferrals

M14 candidate-02 source manifest
`9dd22a598281e83431039f210202a8a24fe49517294cc6bb8d83032c4975e013`
supplies historical mono/RF64, EBU/ITU and 90-minute evidence. Reuse only after
confirming unchanged measurement, input, processing, quantization and transport
dependencies; this is not fresh qualification of M15. Annotation/session/report
changes get fresh targeted evidence. Build identity changes conservatively invalidate
old-build caches even when measurement arithmetic is unchanged.

No new 90-minute workload. Actual Golem/High Sierra, cross-Mac exchange, all video
execution and three-hour qualification remain deferred after M18, no later than
M19. An API/dependency, arithmetic, buffering/clock/latency change or observed
failure must be explained before revising this boundary. Intel build success is
not High Sierra execution. Occasional picture skips remain an accepted limitation;
audio continuity, drift, sustained freezes, identity and lifecycle blockers remain.

## Budget and preservation

Initial free space: approximately 171 GiB. Allow 2 GiB disposable short fixtures
plus separate universal build headroom. Expected short tests 15–30 minutes and
silent UI 15–25 minutes, excluding builds, fixes and real-export availability.
Record observed wall time/RSS, not estimates as results. Baseline refs/status are
saved under `build/m15-qualification-2026-10-29/`. Preserve frozen M14 packages,
all 18 existing annotated tags, original media, video worktree and installations.
No physical device tests are planned. Stop for readiness review before publication.

## Actual Pro Tools adapter follow-up, 29 October

The supplied 10.3.10 export uses `M A R K E R S  L I S T I N G`, fixed-width
tab-separated columns, `LOCATION` timecode and numeric `TIME REFERENCE` with
`UNITS=Samples`. Add this heading and an explicit unit-driven field policy:
audio-sample selection uses the sample reference only with one matching session
sample-rate header and Samples on every row; other selected units use LOCATION.
No origin inference from the marker name or general session headers. Disclose
timecode subframe loss. Keep minimal verbatim headers/table as a regression fixture
and verify the complete original locally without checking its unrelated listings in.

Affected core reader and UI disclosure require a new frozen candidate, all eight
short commands fresh (including the required actual-PT acceptance group), universal
cross-build inspection and focused silent reader/UI checks. Existing processing,
measurement, transport and video algorithms remain unchanged; earlier unaffected
UI sequences may be reused with their original identities. Estimated additional
short tests under three minutes, build approximately 5–10 minutes, <100 MiB extra
disposable UI media; frozen package/build storage separate. No deferred tests.

# Timeline inspection convenience — 7 November 2026

**Follow-up:** [gate-membership correction and plain A shortcut](TimelineGateMembershipFollowup-2026-11-07.md) supersede the original red-marker rule and chord described in this frozen-candidate handover.

Continuation of the preserved M20 UI candidate on `codex/m20-ui-integration`.
This is a display/interaction slice before the user's independent listening review.
M20 publication and main integration remain pending. Candidate 04 and all session files found at the start of this slice remain preserved.

## Behavior

- A thin vertical blend control uses the existing left graph margin. Centre keeps
  both layers at their previous opacity; upwards fades only PCM, downwards fades
  only loudness. Double-click restores centre. Inspection data and processing do
  not depend on visibility. The reusable Sections timeline has the same control.
- The lower divider travels down as well as up. Its minimum inspector height is
  150 px, with the previous 290 px default. When compressed, a viewport retains
  access to every proposal, measurement and detail control. Reset layout resets
  the dividers; it does not change processing.
- Applied and draft/proposed plateau gains appear below the gain graph. Narrow
  intervals keep boundary marks and remain accessible by zoom and previous/next.
  The lane selector distinguishes editable annotation sections from applied and
  proposed processing regions. Clicking a readout selects that interval directly.
  Applied regions remain available after reopening without a live proposal.
  Selected bounds are nominal half-open processing intervals; ramps can cross
  their edges. The detailed readout retains cursor gain and separate global gain.
- Range inspection reports per-channel PCM sample-peak magnitudes and original
  sample positions. A covering overview maximum is exact for the selected range
  only if coverage is accepted/complete and its extremum lies inside the range.
  Otherwise the covering bin interval and limitation are explicit. Programme
  true-peak maxima are shown when their measured position lies in the selection;
  missing local true-peak information is never inferred from PCM. This adds no
  detailed audio reader or independent regional true-peak measurement.
- Cmd+Shift+A switches inspection and audition together between Original/Current
  while stopped and idle. Text entry, active playback/preparation and jobs retain
  their guards. The independent dropdowns remain available. The shortcut never
  presses Play, changes a gain or verifies Current.
- Programme filename is above session status, immediately left of Audio setup;
  its tooltip retains the path. Assembled/composed programmes retain their names.
- Comparison rows are Target, Original, Current. Preset selection/save remain in
  the metric-header lane at every height. Measured/target fonts remain 16 px.
- Command-wheel zoom preserves the source frame under the pointer, subject to
  programme edges and the existing overview minimum resolution.
- Backgrounds reuse the plugin palette: graph RGB (16,30,50), panels (12,22,36),
  cards (10,18,28). Semantic curve/proposal colours remain distinct.
- Short-term loudness remains ungated data. Blue/blue-grey describes measured
  LRA-population membership at matching full-resolution window ends, including
  absolute and relative exclusion. Pending previews have no final gate colours.
- The relative IL gate remains visible in both inspection domains. Only Current
  adds bright-red, double-width segments where completed measured momentary
  energies differ from Original at matching sample-window ends. Display cells
  cover the preceding 100 ms cadence interval. A global trim can colour nearly
  the entire programme; this is not a local-gain map or an audibility judgment.
  Incomplete, missing or unverified Current measurements cannot create red parts.

The relative IL gate is 10 LU below its absolute-gated population, not a fixed
-10 LUFS line. The LRA relative gate is 20 LU below its own absolute-gated
population. The -70 LUFS absolute gate remains separate. See
[EBU Tech 3343](https://tech.ebu.ch/docs/tech/tech3343.pdf) and
[EBU Tech 3342](https://tech.ebu.ch/files/live/sites/tech/files/shared/tech/tech3342.pdf).

## Preservation and qualification scope

No measurement, processing, proposal/search/ranking, export, session-schema or
monitoring-policy change. The 128 MiB combined overview budget is unchanged.
Existing M20 decisions and the independent comparison protocol are retained.
The new build identity requires independently verified listening copies. The old
M20 comparison pathname now holds a later DialogAmbiance session with manual gains
(saved before this slice); it is not overwritten. A test-only refresh utility first
assumed the original history was present and crashed on that missing pointer, after
successfully refreshing M19. The missing-history guard is corrected. To recover a
valid paired comparison, the fixed music-only protocol is repeated with unchanged
settings and algorithms, and both applied plans, all six scalar results, requirements,
preferences, method settings and decisions are compared exactly with the committed
original receipts. New explanation provenance identifies the new assessment; it is
not presented as recovered bytes of the unavailable original M20 attachment. No
policy is tuned against these files. Original numerical evidence remains the
predeclared run; recovery is a reproducibility check, not a new independent holdout.

Focused checks cover the new interactions, min/max positions and incomplete data,
actual native UI, unchanged explicit Apply/remeasure, stale/rejected candidates,
save/reopen, dense layouts and target editing. Native/device tests use injected
outputs and isolated preferences; no physical device opening or playback is part
of this slice. Unchanged long-media, video, cross-Mac exchange and numerical
qualification evidence is reused. Actual High Sierra checks are limited to the
affected UI/target/timeline behavior.

## Manual checklist

1. Open the preserved new app, then the matching reverified M19/M20 session copy.
   Use Cmd+Shift+G in Open session to enter its full path.
2. Build waveform for each inspection domain if its overview is absent. Drag the
   left slider up/down and double-click centre. Check stereo/5.1 labels stay clear.
3. Drag both dividers; scroll the compressed lower inspector to reach its controls.
4. Click Applied/Proposed gain readouts, then use previous/next or Fit selection.
   Select sections to edit a manual draft. Check the peak readout's domain and
   represented interval; Apply and independent Current remeasurement stay explicit.
5. While stopped, use Cmd+Shift+A; confirm both dropdowns change. Command-wheel
   around a point under the pointer. Inspect Current for red measured differences;
   Original retains its ordinary relative-gate line. Review blue-grey S exclusions.
6. Audition Current using the same interface/routing/monitor volume in both music
   sessions. Compare 0:48–1:03 and 1:17–1:40. Listening judgment remains the user's.

## Preserved package and exact evidence

- App implementation commit: `a37938e33ba764ebffff901b7b4cf3864a5b95cc`.
- Final package/test-tool commit: `086f316d8f7e1a6387697487f140ececf3dd60fc`.
- Development app: `/Users/fuzz/.codex/worktrees/a411/LevelScope/build/ui-inspection-2026-11-07/candidate-02/LevelScope Standalone.app`.
- Measurement build: `7910deeff7be3718604a7cd02992894f561266157b72d6d5d38f8dc8068705e5`.
- App SHA-256: `8b4bcc28e1a638ed958daa9adbf52cd06aba16adf9a6085ae20af693cd7599a4`.
- Universal arm64/x86_64; Intel minimum macOS 10.13.

Candidate 01 was tested on both Macs. Candidate 02 has byte-identical application
code/binary and adds only the corrected test utility and music recovery checks.
Both packages remain frozen. The receipt proves this equivalence explicitly.

| Focused group | Sequoia 15.7.9 arm64 | High Sierra 10.13.6 Intel |
| --- | ---: | ---: |
| Timeline/appearance/navigation | 913 | 913 |
| Native proposal/inspection/persistence UI | 103 | 103 |
| Target requirements/presets | 489 | 489 |
| Total | **1,505** | **1,505** |
| Native review heartbeat p95 | 28.6 ms | 50.8 ms |
| Dense-history paging heartbeat p95 | 30.4 ms | 70.0 ms |
| Dense 10,000-annotation lane paint p95 | 12.4 ms | 34.4 ms |
| Cancellation | 29.5 ms | 65.9 ms |

As in the earlier qualification, maximum heartbeat gaps include cold setup,
snapshots and test-only full-text reconstruction (up to 1.895 s on High Sierra).
The p95 results do not claim all cold operations complete within 100 ms. No physical
output was opened. High Sierra tar emitted extended-attribute warnings; the test
runner independently verified every frozen artifact hash before and after the run.

Fixed music recovery passed **65 checks** on Sequoia. Both newly generated outcomes
exactly match the committed original plans, six measured scalar values, requirements,
preferences, method settings, pass decisions and private-pass counts. Each applied
Current also exactly reproduces its new private assessment's sample/true-peak
positions. Both new sessions reopen verified. The old M20 pathname's later manual
DialogAmbiance work remains byte-exact, as does the old M19 file. The corrected
refresh utility is separately checked to reject the changed input cleanly before
measurement or writing. The original test-only crash is retained as a failure.

New listening sessions:

- M19: `/Users/fuzz/.codex/worktrees/a411/LevelScope/build/ui-inspection-2026-11-07/listening-recovery-02/MusicInBackAndForegroundSomeDialog/m19/inspection.levelscope`
- M20: `/Users/fuzz/.codex/worktrees/a411/LevelScope/build/ui-inspection-2026-11-07/listening-recovery-02/MusicInBackAndForegroundSomeDialog/m20/inspection.levelscope`

[Acceptance and preservation receipt](../Tests/Reference/results/ui-inspection-2026-11-07/acceptance.json).
The evidence folder contains exact platform receipts, rendered screenshots, recovery
results/reasons, frozen identities and a SHA-256 inventory. Original audio is not
committed. Main/origin remain `22961f03de7919e1de631f8b13eb4b852108ac0b`.

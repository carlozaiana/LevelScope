# Timeline inspection convenience — 7 November 2026

Continuation of the preserved M20 UI candidate on `codex/m20-ui-integration`.
This is a display/interaction slice before the user's independent listening review.
M20 publication and main integration remain pending. Candidate 04 and the original
M19/M20 comparison sessions remain preserved.

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
The new build identity requires independent remeasurement of the already-applied
music plans; new listening copies preserve the original plans/history/requirements
and compare all six measured scalars and sample/true-peak positions exactly.
No proposal is regenerated or tuned against the held-out files.

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

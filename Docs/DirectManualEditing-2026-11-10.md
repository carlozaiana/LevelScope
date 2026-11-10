# Direct manual editing — 10 November 2026

**Build update:** candidate-02 below is superseded for use by the
[Re-measure preview crash correction](RemeasurePreviewCrash-2026-11-10.md).
The original package and its qualification receipts remain preserved.

Continuation of the reviewed manual-workflow package on `codex/direct-manual-editing`.
The previously preserved package, running user app, listening sessions and main are
not replaced. Publication/main integration remains for review.

## Agreed behaviour

Manual section gains commit on Return or focus loss. Their editable readout appears
under the selected section when there is room; a bottom-row field remains available
for narrow sections and multi-selection. Values replace the absolute local setting,
not an additional gain delta. Other passages and global gain remain intact. Manual
Apply/Discard controls are removed from both windows; proposals retain their existing
explicit Review/Apply/Discard and independent Current measurement workflow.

Waveform left-drag draws a section. Gain-lane drag and waveform right-drag select only
fully enclosed sections. Partially enclosed sections remain untouched. One numeric
entry can set all selected sections, with Mixed shown for differing settings. Existing
pointer-centred Command-wheel zoom, Shift-Command-wheel horizontal scrolling and
ruler pan/seek remain available. Text editing never receives destructive section keys.

Backspace or Delete resets/merges the selected passage immediately. Zero neighbours
need no extra confirmation. Different nonzero neighbours offer Set to 0 dB, Merge left
and Merge right (with their values); equal nonzero neighbours offer one merge choice.
Cancel changes nothing. Each successful action combines geometry and gain into one
undoable edit. Global gain stays unchanged. Necessary gain boundaries remain after a
reset to zero. Redundant adjacent boundaries are removed when their gains match.
The secondary Edit-menu command Remove section boundaries only preserves processing
and can intentionally leave a Mixed annotation over different applied gains.

Command-Z/Shift-Command-Z undo/redo committed manual edits. Return moves focus out of
the gain editor so the next Undo acts on the committed edit. While typing, text-editor
Undo remains local to that field. The session-local journal has no step-count cap;
it stores bounded edit records and an on-disk index, never audio or waveform copies.
History ends at New/Open/restart. Editing after Undo creates a new redo branch.
Unrelated processing/annotation changes cannot be overwritten by an old history.
Disk-write failure refuses the edit. History is not portable programme content and
does not expand `.levelscope` storage or the combined waveform budget.

## Waveform and measurement authority

Changed processing immediately marks Current **Needs re-measurement** and clears
all five measured values, including true peak. No background peak or loudness pass
is triggered. Playback uses the committed processor; measuring remains explicit.
Editing and history operations remain stopped-only and blocked by active jobs.

The Current view reuses a compatible completed Source or prior Current overview.
Its retained extrema are adjusted by the ratio between the old and new gain envelope
at their retained sample positions, including the unchanged canonical transition rule.
It is visibly a CURRENT GAIN PREVIEW, never an accepted measurement. Unchanged sample
positions retain their data. A time-varying interval can contain an unretained extremum;
therefore this preview is not an exact regenerated waveform or an exhaustive peak map.
Inspection cannot report baseline range extrema as newly measured Current peaks.
No new full-programme sample read or overview slab is introduced; the existing shared
128 MiB budget remains. Incomplete analysis drawing retains its separate pending status.

Peak warnings default on and remain selectable without a target. A finite Max Peak
supplies the numeric sample ceiling; No Limit uses strict sample magnitude > 0 dBFS.
The tooltip retains HH:MM:SS.mmm, channel, state, peak, threshold and reduction guidance,
and removes the redundant seconds/sample-frame line. TP! uses existing measured
per-channel maxima. Original maxima above 0 dBTP can be exposed without remeasurement,
but this does not locate every inter-sample over. Unverified Current has no measured
TP markers. No overs in this view is display feedback; previews are explicitly labelled
and unavailable peak data never supplies a clean measurement verdict.

## Commands and transitions

Native macOS File, Edit, View, Transport and Analysis menus reuse existing actions.
Command-R explicitly remeasures Current. Preferences opens the existing audio setup,
including Mac-local monitoring/device preferences; proposal preferences remain under
Analysis. This is the command foundation, not a custom keyboard-assignment editor.

Existing nominal 100 ms linear-in-dB transitions (shortened by adjacent interval
lengths) are shaded in the gain lane. Section-boundary dragging remains available.
Adjustable transition lengths and dedicated ramp handles remain the next separate
processing slice: they need versioned session/processing identity, exact export and
measurement agreement, and must not silently alter existing sessions when a default
preference changes. No core DSP/proposal strategy or processing schema changes here.

## Short manual check

1. Open a programme and select **Inspect Current**. Left-drag a section, click its
   gain readout, enter `-2`, and press Return. Current should immediately show
   **Needs re-measurement**, with five dashes and a labelled waveform preview.
   Draw the next section immediately and enter another value by clicking away.
2. Use Command-Z and Shift-Command-Z to undo/redo gain and section changes. After
   Return, Undo acts on the committed gain; while typing it acts on the field text.
   Undo is session-local and starts fresh after New/Open/restart.
3. Drag across the gain lane, or right-drag across the waveform. Only completely
   enclosed sections should be selected. Enter a gain for the selection, then try
   Backspace/Delete. Zero neighbours merge directly; differing nonzero neighbours
   offer explicit reset/merge choices. One Undo restores the whole operation.
4. Check **Peak warnings** with a finite Max Peak and with No Limit. Hover a circle:
   TC milliseconds remain, and the raw seconds/sample-frame line is gone. No Limit
   uses sample peaks above 0 dBFS; Original can additionally show its already
   measured per-channel true-peak maxima above 0 dBTP.
5. Use Analysis > Re-measure Current State or Command-R when ready. Only a completed
   independent measurement restores the measured values and verdict. Save/reopen
   to check the applied gains. Inspect the shaded transitions and native menus.

## Qualification and preserved build

Preserved test app:
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/ui-direct-editing-2026-11-10/candidate-02/LevelScope Standalone.app`

Exact tested source commit: `082d45ff7d88768f44bd9b557f16da561f2b7d93`.
Measurement build: `b5a537b6aefe2d5fbe4a7d3b6d0d033eb87c0a5dd0e9fce37c03cc59a57079cf`.
Universal arm64/x86_64 package, Intel deployment minimum 10.13. Golem has the same
package at `/Users/fuzz/LevelScope-UI-direct-editing-2026-11-10/candidate-02/LevelScope Standalone.app`.

| Focused group | Sequoia 15.7.9, arm64 | High Sierra 10.13.6, x86_64 |
|---|---:|---:|
| Timeline, manual model, history and persistence | 1,992 pass | 1,992 pass |
| Main UI, proposal integration and commands | 144 pass | 144 pass |
| Target requirements and presets | 489 pass | 489 pass |
| Total | **2,625 pass** | **2,625 pass** |

The timeline group includes native Return-to-Undo dispatch, 260 successive journalled
edits and undos, redo branching, multi-section gain replacement/deletion, neighbour
choices and cancellation, stale bindings, save/reopen, preview budget and retained
positions, progressive views, strict No Limit sample warnings and measured TP-only
overs. The UI group checks all five unavailable Current cells, Unicode menu labels,
shared commands, retained preview after editing, explicit remeasurement and existing
proposal review/apply/history behaviour.

In the synthetic 10,000-section case, direct gain commit with journal writing took
16.0 ms on Sequoia and 34.5 ms on High Sierra. Three-hour synthetic loudness-preview
preparation p95 was 1.36/3.92 ms respectively, retaining 3,455,488 bytes. These are
focused UI fixtures, not renewed long-programme measurement or listening qualification.

Screenshots were reviewed on both actual hosts, including the direct gain field,
selection, reset/merge and boundary dialogs, abbreviated peak tooltip, gain preview
and unavailable Current values. Candidate-01 remains preserved; visual review caught
a Unicode encoding defect in unavailable values and command labels. Candidate-02
corrects it, and all three focused groups were rerun on each host.

[Acceptance receipt](../Tests/Reference/results/ui-direct-editing-2026-11-10/acceptance.json)
records the exact package hashes, host results and preservation checks. All 427
artifacts in the previously delivered package remain byte-identical. A source audit
confirms 168 unchanged product files, including all core measurement, proposal and
DSP files. Both `main` and `origin/main` remain at
`22961f03de7919e1de631f8b13eb4b852108ac0b`.

Checks used injected outputs and isolated settings; no physical interface was opened
and the user's running app was not restarted. Unchanged engine, monitor, device,
video, long-media and independent listening evidence is reused. Publication and main
integration remain for review; this is not a fresh full M20 qualification.

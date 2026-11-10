# Direct manual editing — 10 November 2026

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

## Qualification and preserved build

Pending final focused receipts and immutable package identity. Actual Sequoia and
High Sierra checks use injected outputs/isolated settings, with no physical interface
opening. Unchanged engine, monitor, device, video, long-media and independent listening
evidence is reused; this does not claim fresh full M20 qualification.

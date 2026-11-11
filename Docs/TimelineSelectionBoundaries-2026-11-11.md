# Timeline selection and boundary editing — 11 November 2026

This continuation preserves the direct manual editing workflow and the corrected
remeasurement preview lifecycle. It adds four focused interaction changes:

- Starting measurement transfers keyboard focus from its initiating control to the
  timeline before that control is disabled. Repeated Command-Z and Shift-Command-Z
  work after measurement without another click. Completion does not steal focus
  from a field the user subsequently chose to edit. Undo history remains
  session-local; changing processing invalidates Current verification independently.
- Command-click in the waveform or gain lane toggles individual sections. An
  ordinary click selects one section. Entering one absolute local dB value applies
  it to all selected sections as one undoable action. Unselected gaps and global
  gain remain intact. Clicking a selected gain readout retains the multiple
  selection; its common editor shows Mixed if needed.
- View > Gain transition shading (Current) toggles subtle waveform bands. The
  default is on; Original has none. The displayed support is taken from the existing
  applied gain envelope, including shortened transitions near close boundaries.
  This neither changes DSP nor reads additional audio.
- Command-drag section boundaries to accumulate a shared pending edit. Releasing
  Command allows the existing zoom/pan controls; the batch remains pending. Enter
  opens exactly Move gain with sections and Cancel; Escape cancels directly.
  Dashed amber boundaries and an amber gain preview distinguish pending geometry
  from applied audio. Confirmation commits the entire batch as one Undo step.

A boundary batch is derived from the original sections and the final positions,
not accumulated intermediate mouse positions. Expanding a section extends its
edge gain; shrinking it removes that edge's span. Surviving interior gain changes
stay at their original samples. Boundaries cannot cross or collapse sections.
Main and Sections windows share one batch. Playback, measurement, proposal and
session saving/replacement wait until Enter or Escape resolves it. Gain-changing
confirmation invalidates Current; an annotation-only edit does not manufacture a
measurement change. Stale programme/revision bindings cannot apply.

For separate selections, Delete offers a scoped reset to 0 dB and preserves the
outer boundaries protecting unselected passages. The secondary Remove section
boundaries only action requires a contiguous selection. Neither operation treats
the bounding span of separate selections as selected audio.

## Validation and preserved build

Source and focused tests: `c4b4620de88fb74f44016832fde8d9cba61160b6`.
Preserved app:
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/ui-selection-boundaries-2026-11-11/candidate-01/LevelScope Standalone.app`.

Measurement build: `cc96acef55316d0ad8d0e1ced47587b1c5ec661b188ee9da1c00694d962284e9`.

| Focused checks | Sequoia 15.7.9 arm64 | High Sierra 10.13.6 x86_64 |
|---|---:|---:|
| Timeline, selection, transitions and batch editing | 2,299 pass | 2,299 pass |
| Main UI, native focus/Undo and workflow guards | 163 pass | 163 pass |
| Total | **2,462 pass** | **2,462 pass** |

[Acceptance receipt](../Tests/Reference/results/ui-selection-boundaries-2026-11-11/acceptance.json)
records matching frozen source/build/package identities and unchanged hashes before
and after both host runs. All 425 artifacts in the previous crash-fix candidate
remain byte-identical; earlier packages are also left in place.

The High Sierra copy is preserved at
`/Users/fuzz/LevelScope-UI-selection-boundaries-2026-11-11/candidate-01/LevelScope Standalone.app`.
Both hosts used synthetic fixtures, isolated settings and injected outputs.
The check scope covers
native focus and keyboard dispatch, sparse selection and absolute gains, guarded
Delete, exact envelope support and visibility, batch gestures, zoom persistence,
confirmation/cancellation, shared windows, stale inputs and atomic Undo/Redo.

The existing core measurement/DSP, proposal decisions, export, channel roles,
device/routing/monitor-volume, video, long-duration and listening evidence is reused.
No physical audio device opening or audible test is part of this slice. The
combined 128 MiB waveform allocation is unchanged. Earlier packages remain intact.

## Manual checklist

1. Apply several gain edits, Re-measure, then immediately press Command-Z twice.
   Shift-Command-Z should redo them; processing changes require remeasurement.
2. Click section 3, Command-click 6 and 8, enter a local gain, then press Enter.
   Only those sections should change. Undo should restore all three together.
3. Inspect Current and zoom into a gain boundary. Toggle transition shading in
   View; switch to Original to confirm that Current's processing bands disappear.
4. Command-drag two boundaries, release Command and zoom, then adjust another.
   Enter should offer only Move gain with sections and Cancel. Confirm, then Undo;
   repeat once using Escape to cancel. Check the same workflow in Sections.

The complete continuation history is published for review on
[`codex/timeline-selection-boundaries-2026-11-11`](https://github.com/carlozaiana/LevelScope/tree/codex/timeline-selection-boundaries-2026-11-11),
including the preceding manual-workflow and crash-fix commits. The previous local
branch is retained. No merge or main update is part of this publication;
`main` remains `22961f03de7919e1de631f8b13eb4b852108ac0b`.

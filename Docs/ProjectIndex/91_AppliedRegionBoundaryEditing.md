# Applied processing region boundary editing — 22 November 2026

This is the first workflow-simplification slice after published integrated M20
(`1056ad4ace536abfd68904f083c7fc56a1a8999f`). Development stays on
`codex/applied-region-boundaries-2026-11-22`; main, existing tags, frozen packages,
original media/sessions, listening material and both manuals remain untouched.
Publication and main integration require a separate later instruction.

## Interaction and processing contract

Choose **Select applied regions** in the gain lane. Small tabs at the top of the
gain plot identify the internal nominal boundaries. Drag a tab horizontally;
the vertical stem identifies the cut and the shaded band identifies the actual
centred transition. The amber preview uses the complete edited envelope,
including any neighbouring transitions whose radius changes. Current inspection
also previews those transition bands over the waveform. Original inspection
continues to show Original.

Release commits one processing edit. Escape, loss of keyboard focus, hiding the
view, or changing the viewport cancels the unfinished drag. A click or an
out-and-back gesture preserves the exact original sample cut and has no effect
on processing, measurement authority or Undo. The gesture is anchored to the
original sample position rather than quantising a click to a display pixel.

Applied intervals form a contiguous partition, including unity intervals.
One shared cut changes the preceding interval's end and following interval's
start; both gain values, global gain, all other cuts and section annotations
stay fixed. Programme endpoints are fixed. Cuts cannot cross their neighbours:
the adjoining plateaus retain at least two audio samples, allowing the existing
processor's transition radius to remain valid. Its unchanged 50 ms maximum
half-width shortens automatically where necessary. No smoothing or measurement
mathematics changes.

The live preview is transient and owns no audio. It does not launch a measurement
or change the applied processing. Playback/preparation, measurement, pending
section edits and unapplied proposals retain their existing editing constraints.
While dragging, playback, competing edits, Save, New/Open and measurement wait
for release or cancellation. Runtime/revision/state bindings reject stale gestures.

A real commit uses the existing manual-edit journal, processing invalidation,
gain provenance and session writer. Current immediately needs remeasurement;
previous Current measurements and old proposal authority cannot verify or apply
to the changed processing. Undo/Redo remain available after remeasurement and
do not resurrect measurement or proposal authority. The full last applied
explanation remains available, labelled historical after processing changes.
Saving/reopening preserves the edited envelope; explicit **Re-measure Current
State** verifies it independently.

The Logic ownership review confirms that `protectedLocal`, pins, margins, root
ownership and proposal duration/gain allowances constrain the planner, not a
deliberate manual edit. Those stored source-frame settings remain unchanged.
The next Propose enforces them against independently remeasured edited Current.
Core processing, proposal dispatch, search limits, defaults and strategy are
unchanged. Section drawing, section-marker dragging and Command boundary batches
retain their own interaction and geometry; applied handles are available only
in the explicitly selected applied-region mode.

## Focused qualification

Affected checks cover native timeline gestures, real sample processing, shared
cuts and shortened neighbouring ramps, cancellation/no-ops, stopped editing,
stale bindings, Undo/Redo after remeasurement, save/reopen, main-window guards,
complete historical explanation preservation, and selection/seek/zoom/section
editing conflicts. The supplied schema-3 passing session and audio are copied
before native checks. Hardware outputs are not opened.

Exact candidate identity and actual Sequoia/High Sierra results will be recorded
with the preserved build. Unaffected integrated-M20 evidence from documents
89–90 is inherited rather than rerun. This does not add an independent-programme
correction or listening acceptance claim. Broad video, full-film and three-hour
endurance remain deferred, without an automatic 90-minute substitute. The
inherited completion pause remains a separate limitation. Existing limited
dropped-video-frame acceptance does not waive audio interruption, accumulating
drift or sustained freezes.

## Manual check

1. Open a copy of a session with applied local gains; choose Current inspection
   and **Select applied regions**. Zoom around a transition.
2. Drag its top-edge tab and inspect the moving cut and ramp. Escape once, then
   repeat and release. Check unchanged gain values and **Needs remeasurement**.
3. Undo/Redo, save/reopen, then explicitly remeasure. Undo once more after
   measurement. Inspect the last applied explanation's historical status.
4. With your normal monitoring setup and a comfortable level, explicitly Play
   before the edited entrance/return. Listen across both ramp ends for clicks,
   abrupt changes or an unintended neighbouring gain change. Playback is never
   started by a drag.

Next: simplify the Propose/apply/remeasure workflow in a separate slice. Automatic
application, alternative switching, extension-default evaluation and planning
strategy changes are outside this implementation.

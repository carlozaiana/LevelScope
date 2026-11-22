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

Preserved candidate 01 passed on actual Sequoia 15.7.9 arm64 and actual High
Sierra 10.13.6 (17G66) x86_64. The universal app and both test executables retain
an Intel minimum of macOS 10.13. Both platform snapshots were inspected.

| Affected group | Sequoia | High Sierra |
| --- | ---: | ---: |
| Timeline/manual processing, navigation and existing section editing | 2,580 passed | 2,580 passed |
| Native schema-3 applied boundary / history / save-reopen | 24 passed | 24 passed |
| Main UI, pending proposal guards, independent verification and Undo | 167 passed | 167 passed |
| Total assertions in these affected groups | 2,771 passed | 2,771 passed |

The two independently edited native sessions also have exactly equal serialized
Current processing, annotation geometry, region settings, requirements and proposal
preferences. This is a focused edited-state comparison, not a new full interchange
or numerical qualification. The unchanged writer and measurement/processing Core
retain integrated-M20 evidence. No physical output was opened and no audible test
ran. Unaffected integrated-M20 evidence from documents 89–90 is inherited rather
than rerun. This does not add an independent-programme
correction or listening acceptance claim. Broad video, full-film and three-hour
endurance remain deferred, without an automatic 90-minute substitute. The
inherited completion pause remains a separate limitation. Existing limited
dropped-video-frame acceptance does not waive audio interruption, accumulating
drift or sustained freezes.

## Exact build and preservation

- Implementation commit: `240da530302af8f539d0351b7332798e2373314a`.
- Measurement build: `6a77c53eed0d1e216ecf78750bcc7433a6c5c34455dfd5996c85d6f8334ab819`.
- Source manifest SHA-256: `5dc388feeaab2d5bc97270b9ad6cc9503bc39140df0b123545efacd778b15c04`.
- Package identity SHA-256: `64d5f870f5382f332f58564b945f0ebf61012218570dc1bb970bb06d65171942`.
- App: `/Users/fuzz/.codex/worktrees/a411/LevelScope/build/applied-region-boundaries-2026-11-22/candidate-01/LevelScope Standalone.app`.
- [Acceptance receipt](../../Tests/Qualification/AppliedBoundaries/results/2026-11-22/acceptance.json) and adjacent per-platform receipts include exact commands, hashes and test results.

Final checks confirm all 27 pre-existing tags, nine prior packages and 2,572
protected files unchanged. Local and remote main remain `1056ad4ace536abfd68904f083c7fc56a1a8999f`;
the integrated M20 tag object and target also read back unchanged. Original
session/audio hashes match before and after both native runs. The branch and
this candidate have not been published. This subsequent handover/evidence commit
does not alter the tested Source tree or rebuild the frozen app.

### Development and transfer observations

The first development run in the default sandbox stopped at the pre-existing
worker-prefix readiness check, before any new boundary checks. Running with
normal verified-source cache access reached the new checks. That second run
passed the envelope, cancellation, Undo, persistence and neighbour checks but
exposed a test fixture with no waveform ruler; the fixture now supplies its
actual measured waveform. Both original receipts are preserved. A test-only
mixed-type `auto` declaration was also corrected during compilation. These are
not suppressed failures in the final candidate: the complete frozen groups
passed on both platforms.

High Sierra's older tar rejected modern macOS extended-attribute metadata in
the first transfer archive, before tests started. The initial extraction remains
preserved separately. A GNU-format archive without extended attributes supplied
the identical candidate bytes; the runner verified every packaged artifact
before and after the successful High Sierra run. No product or test criterion
was changed for that transfer.

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

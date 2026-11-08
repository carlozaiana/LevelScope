# Manual timeline workflow — 8 November 2026

Continuation on `codex/m20-ui-integration`. Publication and main integration remain
pending. Previous packages and the user's running app/session are preserved.

## File-limit evaluation

The main-window programme applicability checkbox is removed. New imports use the
existing supported whole-file evaluator automatically, with the visible scope
statement: PASS/FAIL applies to this file and the selected measured limits. This
adds no audio read, measurement pass or processing. An excerpt's result describes
that file, not its parent programme or unselected broadcaster requirements.

Reopening a legacy session with applicability `unconfirmed` adopts this automatic
scope. Explicit unsupported applicability values remain unsupported. Target
methods, tolerances, missing/unstable metrics, processing identities and export
verification are unchanged. Historical proposal explanations are not rewritten.

## Peak hover

Hover an orange circle for Original/Current, channel role, timecode, exact sample
frame, sample peak in dBFS, Target in dBTP and approximate additional reduction.
There is no redundant "Above numeric ceiling" line. TP! hover uses the available
measured programme true-peak maximum and its eighth-frame position instead.

The subtraction refers to the inspected signal, already including applied global
and local gains for Current; monitoring gain is excluded. Sample-based guidance
is not a true-peak compliance promise. Both hover and markers disappear with the
waveform layer and are unavailable for incomplete/unverified states or No Limit.
Hit testing uses existing bounded display columns and reads no audio.

## Manual gain entry and annotation edits

- Remove Set draft in the main lane and both separate draft-creation/set buttons
  in Sections. Return or leaving a valid gain field stages its value. Clicking
  Apply also commits a valid pending field before applying. Invalid input cannot
  accidentally apply the previous value. Escape restores the displayed value;
  Discard retains applied audio. Values are absolute local dB before global gain.
- The explicitly chosen section/applied/proposed selection mode is retained when
  clicking graph or numeric readouts. An unavailable mode falls back to an
  available one; explicit candidate review can still select that candidate.
- Every new manual draft starts from the existing gain plan. Editing one section
  replaces its nominal local gain while preserving all other plateau values,
  including after annotations are drawn or moved. A section spanning different
  applied values displays Mixed. Drawing a section never alters processing.
- Deleting an annotation with applied local processing offers Keep applied gain,
  Remove applied gain, or Cancel. Keep merges the annotation into its neighbour
  without changing processing; annotation Undo remains available. Remove drafts
  local 0 dB in its span, preserving global gain.
- Moving a boundary in either waveform or the numerical Sections editor offers
  Keep gain in place, Move gain with section, or Cancel. With-gain movement extends
  the plateau on the expanding side across the swept interval. Other nominal
  plateaus stay fixed. Existing canonical smoothing recomputes boundary ramps,
  whose support can extend beyond the nominal interval. This does not introduce
  a new free-floating clip-drag interaction.
- With-gain edits display provisional annotations and a manual gain draft as one
  pending operation. Apply commits both; Discard restores both. Saving is blocked
  until resolved. Further geometry/global-gain edits cannot modify this paired
  draft. Actual changed processing leaves Current unverified until independently
  remeasured; annotation-only and audio no-op edits preserve valid measurements.
- A delayed dialog answer is checked against the opening, audio revision and full
  annotation map; stale answers cannot mutate newer state. Proposal candidates
  remain immutable and require their existing explicit Apply/Discard workflow.

The same session schema, core storage and waveform budgets remain in use. The
additional rollback map is transient; it is never saved as applied processing.

## Focused qualification

Exercise automatic evaluation and legacy reopening, exact sample/TP hover and
visibility, numeric Return/focus-loss/invalid input, fixed selection, nested sections,
preserved gains after geometry changes, both dialog choices/Cancel/stale answers,
paired Discard/Apply, independent Current measurement and save/reopen. Reuse
unchanged numerical, monitoring, device, video and long-media evidence. Actual
Sequoia/High Sierra checks use isolated settings and injected outputs only.

Final package identity, receipts and manual checklist follow after qualification.

# Timeline monitoring, deletion and peak overlay — 7 November 2026

Local continuation of the tested UI branch `codex/m20-ui-integration`. Earlier
packages and sessions remain preserved; publication/main integration is pending.

## Monitor overload

The main timeline shows a persistent, high-contrast **PLAYBACK STOPPED — Monitor
overload** banner with the confirmed speaker role, file channel, physical output,
source seconds and exact sample frame. A Monitor volume button opens Audio setup.
The instruction is to lower monitor volume, then explicitly press Play again.
The banner survives Stop/volume adjustment; a new source/audition identity or an
explicit retry clears the old event. A repeated overload produces a fresh warning.

The existing callback protection is unchanged: the whole overloaded callback is
silenced without consuming its source frames. The warning is assembled on the
message thread, with no new callback allocation, locks or device actions. It does
not change processing, monitoring preferences or verification, and never reconnects
or starts playback automatically.

## Backspace in the main gain lane

Click a section/region in the lane, then press plain Backspace. Text editors keep
Backspace for typing. Every deletion consumes its explicit selection, preventing
key repeat from deleting successive sections. Busy/stale inputs are guarded.

- Annotation section: delete its identity/name and merge its span into its preceding
  neighbour; the first merges into the next, retaining the next section's identity.
  Full programme coverage remains; the final whole-programme annotation cannot be
  deleted. Audio/gains/measurements are unchanged. One-level **Undo delete** restores
  the prior annotations; subsequent unrelated annotation edits invalidate that undo.
- Applied gain region: stage 0 dB **local** gain on its existing nominal interval.
  Keep global gain and other plateau values. Show the removal as a manual draft.
  **Apply gains** changes Current and marks it unverified; **Discard gains** retains
  the applied plan. Existing smoothing recomputes boundary transitions on Apply,
  including their support outside nominal boundaries. No gaps/audio cuts are made.
  The draft uses the applied intervals even when annotation geometry differs.
- Proposed/rejected/stale candidate: read-only as before. Propose draws a preview;
  it does not apply it. Explicit Apply proposal is retained. Backspace cannot edit
  a privately verified candidate and still present that candidate as verified.

The new removal draft is transient and uses the existing plan, application and
session-persistence paths. It cannot accidentally enter the annotation-index gain
editor. Save/reopen keeps an applied removal unverified until independently measured.
The last applied proposal's complete explanation remains intact and becomes historical
when the manual removal changes its processing. No session schema or DSP rule change.

## Optional over-target peak display

**Peaks over target**, initially off, uses the selected state's completed measurement
and matching accepted waveform. It also appears in the reusable Sections timeline.

- Orange circles mark retained PCM sample extrema whose absolute amplitude exceeds
  the numeric Max Peak ceiling. They retain exact source sample positions, including
  negative extrema. Dashed amplitude lines locate the ceiling on each channel.
- **TP!** labels flag measured per-channel programme true-peak maxima above the target,
  at their measured eighth-frame positions, on the corresponding waveform channel.
- The header/tooltip distinguish sample extrema from true-peak maxima. Overview
  reduction does not retain every event and cannot establish every inter-sample
  exceedance. Zooming improves the retained sample detail; this is not a complete
  true-peak event map or a compliance verdict.

No Limit has no ceiling to display. Incomplete or unverified states do not acquire
highlighted results. Target edits update the display only; monitor volume is not
part of programme peaks. Visibility follows the waveform layer's blend. Painting
uses the existing bounded columns, with no extra audio reads or peak buffer; the
combined Source/Current waveform budget remains 128 MiB.

## Qualification scope

Focused checks cover overload/recovery using injected outputs; native warning
visibility; Backspace selection/text-entry/repeat/blocked states; first/interior/last
annotation removal and Undo; applied interval removal, discard/apply and unchanged
neighbours/global gain; stale drafts; immutable proposed candidates; historical
explanations/save/reopen; and target/no-limit/stale/boundary sample markers.

Actual Sequoia and High Sierra checks use isolated settings and injected outputs.
No physical opening or audible test is part of this slice. Unchanged M19/M20 engine,
proposal and export evidence and the agreed qualification deferrals remain in force.

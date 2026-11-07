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

## Preserved development build

- Product/test commit: `3c853e4c855b047ed990a899fde586fa2d952aaf`.
- Main implementation: `6d58347bb29c0f121bc1adea1db4c61a9c74123a`; the final
  follow-up corrects the disabled removal field and annotation-reset message.
- App: `/Users/fuzz/.codex/worktrees/a411/LevelScope/build/ui-safety-editing-2026-11-07/candidate-02/LevelScope Standalone.app`.
- High Sierra copy: `/Users/fuzz/LevelScope-UI-safety-editing-2026-11-07/candidate-02/LevelScope Standalone.app`.
- Measurement build: `160f9c55675c802e775cd1fac308cc0caddbbf0a858547dde0efcdcbff412141`.
- Universal arm64/x86_64 app and tests; Intel deployment minimum remains 10.13.

Candidate 01 and prior inspection packages remain preserved. Candidate 02 is the
review build. The running user app was not restarted, no physical device was
opened, and the existing listening session files were not rewritten. `main` and
`origin/main` remain at `22961f03de7919e1de631f8b13eb4b852108ac0b`.

Final candidate 02 passed the same affected groups on both actual hosts:

| Focused group | Sequoia 15.7.9 arm64 | High Sierra 10.13.6 Intel |
| --- | ---: | ---: |
| Timeline, deletion and peak display | 1,199 | 1,199 |
| Native UI, proposals and save/reopen history | 128 | 128 |
| Injected monitor guard and recovery | 111 | 111 |
| Total checks | **1,438** | **1,438** |

The native overload banner and removal draft were visually reviewed on both
hosts. The peak evidence includes an independently measured inter-sample-only
exceedance: TP! appears without inventing a PCM sample exceedance. A source audit
verifies 170 unchanged files against the task baseline; the ten changed product
files and exact scope are recorded. Prior candidate 04's 420 artifact hashes and
both existing listening-session hashes match their state at the start of this task.

[Receipts, manifests and rendered evidence](../Tests/Reference/results/ui-safety-editing-2026-11-07/README.md)
include the intermediate candidate and the final readout correction. Full-duration,
video, cross-Mac session exchange and physical/listening qualification were not
expanded. Existing evidence and agreed deferrals remain in force.

## Manual checklist

1. Save your current session before quitting the older running app and opening the
   new preserved build. Existing sessions and comparison copies have not been
   rewritten. This build has a new measurement identity: if saved measurements are
   rejected as stale, measure Source and Current again before inspecting peak results.
2. In the main gain lane choose **Select sections**, click an annotation, then
   Backspace. Check the neighbour takes over its span and **Undo delete** restores
   it. Applied gains and valid measurements should stay unchanged. Backspace inside
   a text field continues to edit text; it does not delete a lane item.
3. Choose **Select applied regions**, click a nonzero region, then Backspace. Check
   the amber draft and disabled draft field show 0 dB local gain. Try **Discard gains**;
   repeat and choose **Apply gains**. Current must become unverified. Removing
   attenuation can raise the level; global gain and neighbouring plateau values
   remain unchanged, while existing smoothing recomputes the boundary ramps.
4. Set a finite Target Max Peak, select a completed measured state and enable
   **Peaks over target**. Zoom around an orange sample marker or **TP!** marker;
   toggle the layer off again. Check No Limit or unverified Current disables it.
   The tooltip explains the difference between retained PCM extrema and measured
   programme TP maxima; it does not promise every inter-sample event.
5. If the monitor guard stops playback during normal testing, check the prominent
   banner names the role/channel/output and position. Use **Monitor volume...** to
   lower the monitor level, then explicitly Play. Lowering the level alone must
   not resume playback. There is no need to deliberately provoke a physical overload.

Propose still creates a preview and **Apply proposal** remains explicit. The
earlier supplied comparison sessions already contained their applied proposals.
Deleting an annotation or drafting a region removal does not turn a private
candidate measurement into Current verification.

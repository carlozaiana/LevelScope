# Applied-region boundaries in the main graph — 24 November 2026

This focused UI continuation starts from the exact accepted shared-source
candidate-03 handover `f9eeae554421b4ef4f8b48d53a69785aa3f6af46`, including
implementation `614377f4e573bdd0fc34726584418b9333d092c2` and its HFS+ integrity
correction. Work remains on `codex/graph-region-boundaries-2026-11-24`. Main,
existing tags, earlier packages, manuals and original media/sessions are unchanged.
Publication and main integration are not authorized by this slice.

## Interaction

While stopped, choose **Inspect Current** and **Select applied regions**. The
mint vertical lines and small top tabs in the waveform/curves graph mark the
actual internal processing cuts. Hover gives a horizontal resize cursor; drag
anywhere on that line. The amber cut and shaded transition follow the pointer,
and the gain lane previews the same edited envelope. Active preview bands remain
visible even when normal transition shading has been turned off.

Both entry points use the existing gain-lane boundary transaction, not a second
editing implementation. Release changes the actual processing once and creates
one Undo entry. The shared cut moves both adjoining interval edges; their gains,
global gain, other nominal cuts, annotations and planner settings remain fixed.
Programme ends cannot move and neighbouring cuts cannot cross. The processor
continues to resize shortened ramps with the existing two-sample plateau minimum.
Document 91's deliberate manual-edit versus planner-protection semantics remain
unchanged.

**Escape while the mouse button is still held** cancels the preview. After
release, use Undo. No-motion and out-and-back gestures preserve the exact sample
position and verification. Focus loss, hiding, viewport/size changes, changing
away from applied-region selection or leaving Current inspection cancel an
unfinished graph edit. Existing job, playback, pending-edit and stale-session
guards apply to both entry points.

Selection mode decides which overlapping line is editable: section mode edits
section markers, applied-region mode edits Current processing. Original's graph
does not expose Current-processing handles; the gain lane retains its existing
explicit Current-editing semantics. Alt-drag pans even over a handle. Ruler seek,
pointer-centred Command-wheel zoom, region inspection, section ranges and Command
section-boundary batches retain their separate gestures.

A real release immediately removes Current verification and invalidates old
proposal authority. The previous measured explanation remains historical. The
next explicit Play uses the edited Current before remeasurement, through the
accepted shared-source preparation path. A drag does not start playback, open an
output, change routes, run the planner or launch measurement on mouse movement.
Independent **Re-measure Current State** remains the authority for the edited
Current. Save/reopen retains the edited envelope and complete historical explanation.

## Focused results and limits

The exact preserved candidate passed on actual Sequoia 15.7.9 arm64 and
High Sierra 10.13.6 (17G66) x86_64:

| Group | Sequoia | High Sierra |
| --- | ---: | ---: |
| Graph/gain-lane equivalence, exact samples, cancellation, Undo/Redo, section editing and navigation | 943 | 943 |
| Native main-window graph, supplied schema-3 history/save/reopen, stereo and 5.1 edited playback | 563 | 563 |
| Total | 1,506 | 1,506 |

The playback checks explicitly Play unverified edited Current and compare every
captured pre-monitor sample with the committed gain processor, including all
transitions. They pass in stereo and 5.1 without another full-source copy or
verification. A subsequent normal Current measurement succeeds. Outputs are
injected: no physical audio device was opened, no sound or latency benchmark
was run. Native preview and committed-state screenshots were visually inspected.
Commands, identities and individual assertions are in the
[both-platform acceptance receipt](../../Tests/Qualification/GraphBoundaries/results/combined-2026-11-24/acceptance.json).
The initial Sequoia-only receipt and all earlier evidence are preserved unchanged.

The universal app and test binaries contain arm64 and x86_64, with Intel minimum
macOS 10.13. After the user restored Golem, both new graph groups passed there
using the identical package, without rebuilding or changing any test. All package
artifacts were verified before and after each platform run; input session and
audio hashes remained unchanged. High Sierra native screenshots were inspected.
Candidate 03's 227,322 checks per actual platform remain inherited evidence for
unchanged source integrity, playback, processing and proposal paths. No
historical full-suite repeat was necessary for these four UI
production files. Measurement mathematics, planner/search/defaults, dispatch,
session writer, source service, routing and transport implementations are unchanged.

Broad video, full-film and three-hour endurance remain deferred, with no automatic
90-minute substitute. The known completion pause remains separately documented.
Limited acceptance of occasional dropped video frames does not waive audio
interruptions, accumulating drift or sustained freezes.

Before this slice the user reported immediate playback after import/drawing and
seek/replay on a 93-minute 5.1 file in the accepted shared-source build. This is
user interaction feedback, not an endurance run, physical latency measurement or
acceptance of this new graph build. New listening acceptance remains the user's
separate check.

## Identity and preservation

- Implementation: `f801ac06307f25f423a953a7e35b59c28013ccc7`.
- Source tree: `31be38bd890fee3667a6a976fbc3524a2b26e609`.
- Measurement build: `fcb60b72c950d8b21e3c8ae4be8ce57bf95f6991c27aee84073d4d092bd78ea7`.
- Source manifest SHA-256: `460359c1b4d24cc61c1a9be3176933b8ff725b96af1925718a0341ab5d42556c`.
- Package identity SHA-256: `cecac62f3f98a5c92e1b87d938ab94f813b32f7cd2984cf8e7156a7bf7af43cf`.
- App: `/Users/fuzz/.codex/worktrees/a411/LevelScope/build/graph-region-boundaries-2026-11-24/candidate-01/LevelScope Standalone.app`.

The new package contains the exact built app and two companion executables,
frozen source, build commands and hashes. Prior packages are not rebuilt.
Preservation checks cover the preceding eleven packages/3,638 protected files,
all 477 shared-source candidate-03 artifacts, its identity, all 28 tags and main
at `18adeb520006b7fe7e0a22300cc26abd28eec51c`. No discrepancies were found.
The handover is documentation/evidence only and retains the tested Source tree.

One initial development native run stopped at an inherited fixture assertion
that expected no preparation worker after Open. The accepted background source
preparation was still finishing. The graph-mode fixture now waits for that
existing preparation before testing direct disabled transport callbacks during
a drag. Product logic and guard criteria were unchanged. That failed receipt,
the successful development runs and both frozen-package groups on both Macs are preserved.

## User check

1. Open a copy with applied gains. Select **Inspect Current** and **Select applied
   regions**, then zoom around a transition.
2. Drag its line in the main graph. Check that the cut and transition move in both
   graph and gain lane. While still holding the mouse, press Escape; then try
   another drag and release. Release should show **Needs remeasurement**.
3. Undo/Redo, switch to section mode and check the section markers remain distinct.
   Try ruler seek, Alt-drag pan and Command-wheel zoom.
4. Using your normal monitoring setup and a comfortable level, explicitly Play
   across the edited transition before remeasurement. Listen for the intended
   change and for clicks or unintended neighbouring changes. Save/reopen and
   independently remeasure when ready.

Remaining: separate user interaction/listening feedback. The focused Golem
graph checks are complete. Publication/main integration remain separate.

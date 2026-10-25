# Additional candidate-01 silent UI observations

Recorded 2026-10-24. **Partial evidence only; M7 and candidate-02 qualification
remain pending.** The earlier in-progress record is preserved. Candidate-01's
three-hour audition failure remains a failure. These observations cover unchanged
UI/measurement behavior; candidate-02 requires its own affected audition checks,
automated suites and muted physical video evidence.

All interactions below used the isolated candidate-01 app, with devices closed,
routes muted except for the earlier explicitly selected partial-monitoring map,
and playback stopped. No audio device was opened by these UI checks.

On actual High Sierra:

- The malformed eight-channel mask reported “channel mask does not describe 8
  file channels”; the previous Source/Current measurements and gain remained.
- The unsupported SDDS mask imported unassigned and unconfirmed, clearing Current.
  The WAVE/FLAC preset still required Confirm roles before measurement. Original
  `0xff` metadata and the explicit manual interpretation remained visible.
- Editing Channel 1 to ambiguous `Ls` invalidated Source and disabled confirmation.
  Explicit `M+030` enabled confirmation and normalized to Left.
- Unlabelled AIFF required manual roles. Native FLAC suggested the supported
  original file order, required confirmation, and measured correctly.
- Cancelling long FLAC import preserved the previous valid FLAC and measurements.
- Cancelling Source measurement discarded its result. Cancelling Current
  measurement left it UNVERIFIED, retaining gain without automatic retry.
- Clear Source removed Source and Current. Closing during active measurement
  removed the app window; Finder remained responsive and a subsequent process
  check found no Standalone process.

On arm64 Sequoia:

- The user selected the first explicit System I preset for unlabelled WAV and
  confirmed roles. The UI recorded manual interpretation and completed Source
  measurement with the device closed.
- Unlabelled AIFF required manual roles. Native FLAC suggested the supported
  original file order, required confirmation and completed measurement.
- Cancelling long FLAC import preserved the prior measured FLAC source.
- Source measurement cancellation discarded the result. Current cancellation
  displayed UNVERIFIED with the applied gain retained and no automatic retry.
- Clear Source removed both states. Closing during active measurement terminated
  the app; the subsequent accessibility call returned `procNotFound`.

The detailed local accessibility observations and screen captures remain in the
ignored `build/m7-qualification-2026-10-24/ui-evidence` folder. Remote desktop
screenshots include unrelated desktop material and are not published here.

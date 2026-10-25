# Candidate-02 affected silent UI observations

Recorded 2026-10-24 on arm64 Sequoia and actual Intel High Sierra.
**These observations do not establish complete M7 qualification.** Long-programme,
cross-platform and muted physical-video results must be assessed separately.

Candidate source digest:
`c3df888f567a69adb83cbb6c7d2c7a1dcdf1da9498d43d688ffe890786c0e43e`.
Standalone binary SHA-256:
`c1686278deb2ea36b926ccb5cb9df10b169bec903594e38680fc5cbfa4817248`.
Both apps ran from their isolated `candidate-02` packages. Devices remained
closed, no output was selected in these fresh sessions, all routes were Muted,
and no playback was started.

Observed on both Macs:

- Native eight-channel FLAC import suggested L/R/C/LFE/rears/sides in original
  file order and still required explicit role confirmation.
- Source preparation completed without a selected device or a measurement.
  It displayed “Verified source ready; playback stopped”.
- An exact stopped seek to `383999/48000` was accepted without playback or
  device opening. The automated suite separately exercises actual decoded
  samples at seek positions, including the three-hour final sample.
- Current was initialized and its gain set to `-3.305074118944` dB. Explicit
  Current selection displayed that full gain and UNVERIFIED. Preparing again
  reused the verified source and displayed “device and playback unchanged”.
- Source and Current measurements completed with no device selected. Current
  became “Measurements current”; with no profile selected, its evaluation
  remained NOT EVALUATED. Measured availability did not imply compliance.
- Replacing the source with the long eight-channel clock FLAC cleared Current.
  After confirmation and initialization, preparing Current showed progress.
  Stop during preparation displayed “Source preparation cancelled; stopped”.
- Starting preparation again and closing while progress was active terminated
  the app normally. Subsequent process inspection found no Standalone process.

Additional Sequoia identity observations:

- Adjusting monitor volume to -15 dB retained the Current processing gain and
  unverified state.
- A measured Source was invalidated by changing Channel 1 from Left to a
  duplicate Right. Confirmation and measurement became unavailable.
- Restoring Left required explicit confirmation; confirming did not restore
  the discarded measurement.

The [source comparison](source-differences.json) supports reuse of candidate-01's
unchanged non-audition UI observations, including applicability, M3 application,
actual measured pass/failure, all eight diagnostic/section peaks, partial monitor
routing, video-offset identity, malformed import preservation, manual presets,
and Source/Current cancellation. Candidate-01's failed long seek is not reused
as passing evidence.

Detailed accessibility records and remote screenshot hashes are retained locally
under `build/m7-qualification-2026-10-24/ui-evidence`. Screenshots containing
unrelated remote desktop content are excluded from repository evidence.

# Main integration and accepted milestone handover

Integration continuation, 17 October 2026. Historical M1 and video reports remain
unchanged; this document records the additional gate and subsequent integration.

## Measurement stage completed

`main` was fast-forwarded from `556d387ad772dc5e604ead95f1525f60e66fd5f9` to qualified
measurement checkpoint `b29bc7c183743f35c6271d2d7ec8fb5ba1464591`, then pushed and
verified directly on GitHub. Fetched remote main had no newer changes. Neither
working copy was switched, reset, relocated or removed. No redundant merge commit
was created. Existing M1 evidence was verified: all 25 evidence hashes and all 76
identity-file hashes match their tested implementation. The only manifest-file
change from that implementation to the measurement checkpoint is its README.

Standalone measurement acceptance does not depend on VST3 measurement parity.
[Plugin parity remains a separate milestone](16_Vst3MeasurementParityMilestone.md).
No plugin target, installed application, Pro40 component or system audio setting
was changed.

## Final-candidate video acceptance addendum

The remaining qualification gap was a real 30-minute run of final candidate 06 on
each actual platform. The original endurance/listening evidence used harness 05,
before the MP4 compressed output-timestamp correction. Both new runs used the exact
preserved final-candidate harness, SHA-256
`291bc99b71815056ce9dcb8cf5e20974a48b79bf66ea042b56e33533e926b231`, with implementation
`54e2aa9d2773229dd666bb36df410ff932e9f9bf`. No source change or rebuild was needed.

| Runtime | Source seconds | Paints | Max CPU paint ms | Late polls | Skipped frames | Max selected PTS skew seconds |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| sequoia-arm64 | 1800.007771 | 107,766 | 12.820625 | 17 | 0 | 0.000000000 |
| high-sierra-x86_64 | 1800.008481 | 107,880 | 0.386965 | 41 | 0 | 0.000000000 |

Both reports pass with zero blank paints, one-frame-bound failures, paint-deadline
failures and audio underruns. Both use the same 54,000-frame, 30000/1001 reference
movie; HAL output timestamps are available. Routes were muted throughout: Mac Studio
Speakers at 48 kHz and Golem Built-in Output/Headphones at its existing 44.1 kHz.
Before/after audio-device inventories match. Playback windows were observed advancing;
this is not continuous optical observation. Temporary caffeinate assertions kept
the runs awake and ended with their processes; no persistent power setting changed.

These are internal presentation-clock and UI paint-selection measurements. The
expected-frame calculation shares the rational timing mapping; independent decoded
burn-in checks validate frame/index identity. Zero selected-PTS skew does not measure
display scanout, display processing, DAC/acoustic latency or physical lip sync.
The user's earlier listening observation, “Heard clicks; flashes seemed aligned,”
remains attributed to harness 05. No new audible test or physical measurement is
claimed by these muted runs.

## Acceptance criteria and reused evidence

| Requirement | Evidence and disposition |
| --- | --- |
| Six frame rates, hour-long rational simulation, VFR, signed offsets and drop-frame boundaries | Existing final-candidate tests on both actual platforms; unchanged source |
| Indexed PTS, coverage, seek and independent burn-in agreement | Existing final-candidate native tests, including all 54,000 long-movie PTS/coverage pairs |
| Real 30-minute internal timing bound on both OS versions | New final-candidate runs above close the harness-version gap |
| VFR and MP4 edit-list timing on real outputs | Existing final-candidate 30-second VFR and 4-second MP4 runs on both systems remain valid |
| Stop/resume/restart, seeks, EOF, replacement, stale work, decode failure and bounded queues | Existing injected-device and native lifecycle checks on both systems |
| Embedded audio excluded; selected audio, routing and measurement identity preserved | Existing native two-stream/embedded-audio fixture, transport/independence checks and M1 suites |
| Product UI and listening observations | Existing recorded observations only; High Sierra product UI coverage is narrower than Sequoia |
| Physical synchronization measurement | Not part of this internal one-source-frame acceptance gate; remains unmeasured |

The final candidate satisfies the defined reference-picture qualification gate.
The integration regression pass and final main handover are recorded below once
the accepted branch revision has been integrated locally and verified.

## Integration status

Measurement is on remote main. The final-candidate video gate is complete; publish
this acceptance addendum on `codex/video-implementation`, then fast-forward local
main if remote ancestry still permits it. Run the named video and M1 regression
executables against the exact preserved binary hashes before publishing video on
main. Do not force a branch update or move an existing checkpoint tag.

## Evidence, preserved locations and limits

New evidence is in
[`Tests/Video/results/main-integration-2026-10-17`](../../Tests/Video/results/main-integration-2026-10-17),
including executable/source bindings, before/after device inventories and retained
tag targets. Original records remain in
[`video-qualification-2026-10-17`](../../Tests/Video/results/video-qualification-2026-10-17)
and the M1 result directory, with their original tested identities and paths.

The video working copy remains
`/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation`. Historical tested
products remain under `/Users/fuzz/Developer/LevelScope/build/video-implementation`,
including its `preserved/video-checkpoint-2026-10-17` snapshot. New raw integration
results use `/Users/fuzz/Developer/LevelScope/build/video-integration-2026-10-17`.
Golem's preserved copies are read-only inputs; new reports use only
`/Users/fuzz/LevelScope-video-integration-2026-10-17`. No original binary, fixture,
qualification report or CMake cache was overwritten. No build was performed.
Committed text copies normalize line endings and trailing whitespace; raw outputs
remain intact in the new integration-results folders and report values are unchanged.
Future builds still require a fresh directory inside the relevant working copy;
the old absolute-path CMake cache is historical and not reusable.

Supported picture remains bounded local MOV/MP4/M4V, SDR H.264, square pixels/full
aperture, at most 3840x2160, two million indexed frames and eight queued pictures.
Embedded-soundtrack playback/selection, HDR/general codecs, video editing/export,
session persistence and VST3 parity remain outside this milestone. Filesystem-read
cancellation is cooperative with no hard OS-read timeout. Measurement retains the
M1 scope and limitations in [document 17](17_M1QualificationResults.md); this is not
a delivery-compliance, processed-audition or public-export qualification.

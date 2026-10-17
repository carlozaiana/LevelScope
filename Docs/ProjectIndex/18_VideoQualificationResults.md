# Reference-picture implementation and qualification

Prepared 17 October 2026 for `codex/video-implementation`, based on annotated
`measurement-milestone-2026-10-17` (peeled commit
`b29bc7c183743f35c6271d2d7ec8fb5ba1464591`). See
[the contract](15_VideoImplementationHandover.md) and
[reproduction instructions](../../Tests/Video/README.md).

## Implemented behavior

The standalone application attaches local MOV/MP4/M4V SDR H.264 reference picture,
selects a video track, reports content identity, frame coverage, rational timing,
VFR and available timecode, and accepts exact signed offsets. Positive offsets
start picture later: `videoTime = audioTime - offset`. Selected audio remains
authoritative; embedded audio is counted as metadata and never decoded or routed.

The audio presentation clock maps consumed source spans to HAL output timestamps,
adding device/stream hardware latency once. Its explicitly labelled fallback uses
callback time plus total reported latency. Stop retains position, Play resumes,
Restart begins at zero, and seek flushes and changes generation before preroll.
Video work runs off the audio callback, queues at most eight pictures, clears stale
generations and displays explicit buffering, missing-file and unsupported states.
Attaching video and changing offsets leave measurement/routing identity unchanged.

The supported raster is square-pixel/full-aperture H.264 up to 3840x2160, with a
two-million-frame index limit. Codec support is bounded; this is not a general media
player. Session persistence, video export/editing, HDR, embedded-soundtrack selection
or mixing, and VST3 measurement parity remain outside this milestone. File reads
use cooperative cancellation; a blocked OS filesystem read has no hard timeout.

## Builds and retained evidence

Evidence is under
[`Tests/Video/results/video-qualification-2026-10-17`](../../Tests/Video/results/video-qualification-2026-10-17).
The `identity` directory binds reports to sources, executable hashes and generated
fixture hashes. Binaries/media are kept in the isolated build/test directories,
not Git. Committed text copies use LF endings and remove trailing whitespace;
report values are unchanged, and raw outputs remain in the working evidence folders.
All binaries are universal arm64/x86_64 with macOS 10.13 deployment minimum;
actual execution used Sequoia 15.7.9 arm64 and High Sierra 10.13.6 x86_64.

Final candidate 06 uses measurement build identity
`f0affb955695cc2b104e89fdac22eeb968813a4e9ca45093b098ccc3c84fcf42`.
Its standalone executable SHA-256 is
`5d6253740a7b78266973e4db4a294d64978fb0a7f83032cd6c9809541583bed1`;
its playback harness is
`291bc99b71815056ce9dcb8cf5e20974a48b79bf66ea042b56e33533e926b231`.
The 30-minute and listening runs use the preceding harness 05,
`d7bef4f1f233caf2f75e992ed9ace55c962057d50c5f7dac09e3f73b1ada5c30`.
This distinction is deliberate: while those runs continued, compressed indexing
was corrected to use output timestamps for MP4 edit lists. Final-candidate tests
verify all 54,000 PTS/coverage pairs of the same long MOV and decoded burn-ins across
its extent. Final-candidate short VFR/MP4 hardware runs also pass on both systems.
The final executable itself is not claimed to have run for 30 minutes.

## Numeric and lifecycle results

Both actual systems pass:

| Check | Result |
| --- | --- |
| Rational timing, transport, native fixtures, long index | 329,054 checks |
| Measurement unit checks | 184 checks |
| Hardware-independent scaffold | 205 checks |
| Official ITU vectors | 39/39; zero baseline comparison errors |
| EBU v05 minimum requirements | 66/66 files, 72 rules, 1,052 readings |
| EBU supplementary alignment | 4/4 project sanity checks |
| Reference runner validation | 8 tests |

Sequoia CTest additionally passes all three registered suites. Native tests verify
frame-number burn-in against indexed PTS, random and sequential decoding, MOV and
actual MPEG-4 containers (both .mp4 and .m4v extensions), VFR gaps, timecode origin,
multiple streams with embedded audio, cancelled/stale work, queue bounds, clear,
close during work, corrupt media and deletion before a new seek. Injected-device
tests exercise resume/restart, repeated seek, EOF draining, device failure,
underrun-source holds and all-muted harness output. Six frame rates are checked
over simulated hour-long timelines, including signed rational offsets and overflow.

M1 reports are compared with their pinned baseline reports and manifests. No
measurement-kernel or plugin source changed. Existing M1 files are retained intact.

## Playback observations

Both real-time 30-minute fractional-rate runs completed successfully, using the
54,000-frame 30000/1001 movie and source-sample-timed clicks. These used harness 05
as identified above, under concurrent build/test/UI activity.

| Platform | Source elapsed | Paints | Max CPU draw | Late polls | Blank / bound / deadline / underrun failures |
| --- | ---: | ---: | ---: | ---: | --- |
| Sequoia | 1800.000333 s | 107,790 | 7.506417 ms | 17 | 0 / 0 / 0 / 0 |
| High Sierra | 1800.005941 s | 107,879 | 0.406 ms | 41 | 0 / 0 / 0 / 0 |

Both runs report zero selected-frame PTS skew, zero first/last-minute mean-skew
difference and zero skipped indexed frames. Late polls include preroll; no blank
paints were recorded after the two-second preroll exclusion. The maximum draw
duration stayed below the 33.3667 ms CFR source-frame duration on each machine.

Final-candidate short runs, with all routes muted:

| Platform | Material / duration | Paints | Max CPU draw | Blank / bound / deadline / underrun failures |
| --- | --- | ---: | ---: | --- |
| Sequoia | VFR / 30 s | 1,678 | 1.509125 ms | 0 / 0 / 0 / 0 |
| High Sierra | VFR / 30 s | 1,680 | 0.269152 ms | 0 / 0 / 0 / 0 |
| Sequoia | MP4 / 4 s | 120 | 0.454958 ms | 0 / 0 / 0 / 0 |
| High Sierra | MP4 / 4 s | 120 | 0.246269 ms | 0 / 0 / 0 / 0 |

Each recorded maximum selected-frame PTS skew is zero. Selection samples the clock
at paint time. Its expected-frame calculation shares the rational mapping verified
by independent burn-in tests; zero internal skew does **not** establish zero physical
latency. CPU draw completion is checked against the expected frame's actual duration.
Neither UI paint nor native-layer submission measures display scanout, processing,
DAC/acoustic delay or optical/acoustic lip sync. Window occlusion/focus changed
during concurrent UI checks; this is not a continuous optical observation.

Mac Studio used its named speakers at 48 kHz, no inputs, with -30 dB monitor gain.
Golem used Built-in Output/Headphones at 44.1 kHz, no inputs, with endurance routes
muted. Reported total output latencies were 829 and 621 samples respectively, with
HAL timestamps available. The user connected/tested headphones and selected that
44.1 kHz output; the tests matched it rather than changing system configuration.

The user observed a 20-second headphone run on Golem's own Screen 1 and reported:
“Heard clicks; flashes seemed aligned.” Its source peak amplitude was 0.4 with
-30 dB monitor gain. This is a useful subjective observation, not a measured
millisecond bound. Earlier attempts were inconclusive (obscured picture, then
inaudibly quiet clicks) and are not counted as listening passes.

Golem has a GTX 670 and three 60 Hz displays (SyncMaster 1920x1080 main,
PL2773HD 1920x1080, PB248 1920x1200). Its own Screen 1 was used for listening.
Sequoia's system profiler reported M1 Max graphics but no display model; its physical
panel and scanout latency are therefore unqualified. Hardware reports are retained.

## UI observations and discovered defects

On Sequoia the candidate UI was observed importing and measuring Source/Current
audio with device closed, preserving those results through reference attachment,
choosing the second video stream, accepting decimal/fraction offsets, seeking while
stopped and blanking outside coverage. The final candidate displays 01:00:00;00 at
frame zero and 01:00:01;00 at a -1.001 s offset; +1.001 s blanks picture. Clear resets
offset/metadata, and the remuxed MP4 now starts at frame zero/PTS zero. High Sierra's
candidate launch and the shared decoder/transport/paint harness were exercised;
only directly observed GUI actions are claimed here.

A VFR probe exposed stale timer selection across a long-to-short frame boundary:
nine one-frame-bound failures occurred in four seconds. Selection was moved to the
start of paint in both product and harness. The failure log is retained under
`superseded`; subsequent VFR runs pass. Earlier endurance attempts stopped for this
fix have no completion reports and are not counted as completed runs.

Another regression fixture exposed MP4 raw compressed timestamps shifted two
frames by the container edit list. Indexing now uses Core Media output timing,
matching decoded presentation timestamps. The VFR fixture also has a genuine final
sample gap; a test expectation was corrected after inspecting its sample durations,
without extending the product's last frame into that gap.

## Preservation and delivery

The original checkout remains on `codex/logic-module` at `b29bc7c...`; main remains
`556d387ad772dc5e604ead95f1525f60e66fd5f9`. Existing annotated baseline tags and M1
evidence are unchanged. Only named standalone/test targets were built. Installed
apps/plugins, Pro40 code/drivers/services, and system audio defaults were not changed.
The existing shared Golem SSH handle was reused without closing or modifying it.

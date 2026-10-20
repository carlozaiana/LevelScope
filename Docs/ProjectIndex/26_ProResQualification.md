# Bounded SDR ProRes qualification — 20–21 October 2026

**Result: the bounded extension meets its required acceptance criteria on both Macs.**
The final candidate passes all 30 required profile/rate combinations, both actual-file
playback runs, native/isolation/UI/regression checks, and all six fresh 1800-second
endurance runs. This verdict uses the explicitly user-approved memory-method
revision in [25_ProResAcceptance.md](25_ProResAcceptance.md). Earlier failures retain
their original verdicts; the zero-skipped-frame gate was never relaxed.

The source and evidence checkpoint is [video-prores-checkpoint-2026-10-21](https://github.com/carlozaiana/LevelScope/tree/video-prores-checkpoint-2026-10-21)
on [codex/video-implementation](https://github.com/carlozaiana/LevelScope/tree/codex/video-implementation).
This extension is not integrated into main by this task.

## Scope and exact candidate

Development incorporated published main `6d304dbecd8fa47f56b0260105aa9d07f96ba1bc`,
including qualified M5, before implementation. The frozen tested source is
`6662091f1dfaff948a4a3eb3f000f47e37cfdc21` (candidate 04). Production source last
changed at `3920d08f849d2420509e1d789c5b8208e50d79af`; the later commit changes the
approved memory harness and its runner/verifier. Subsequent checkpoint changes
contain documentation and text evidence only.

Support is local MOV containing ProRes 422 Proxy, LT or standard 422; progressive,
square pixels, full centred aperture and explicit Rec.709 SDR. Required qualified
combinations are 1920×1080 at 24000/1001, 24, 25, 30000/1001 and 30 fps.
AVFoundation decodes into the existing 8-bit reference preview. Its existing raster
ceiling is not a qualification claim for 4K or higher frame rates.

Selected audio remains the authoritative clock. Exact PTS, signed offsets, timecode,
seeking, stream selection, transport generations and the eight-picture queue are
preserved. Embedded audio is counted as metadata, **never decoded or routed**.
The ProRes path uses separate immutable native image backing for each picture,
avoiding repeated Core Graphics copies during paint. H.264 image backing and
measurement/audio algorithms are unchanged.

The same universal executable bytes run natively on arm64 Sequoia 15.7.9 (24G830)
and x86_64 High Sierra 10.13.6 (17G66). JUCE 7.0.12 is pinned at
`4f43011b96eb0636104cb3e433894cda98243626`; deployment minimum is macOS 10.13.

| Identity | SHA-256 |
| --- | --- |
| Candidate 04 identity JSON | `45219de9006b845230c5e1296859bff2c483beb49c4b827b8c33c552f8381285` |
| Candidate source manifest digest | `de0a4e2a18bd9e250b36aafdb1f474be447316e7283d8698ed7a7672cfae2db6` |
| Standalone executable | `5c26cb8faacb1e7b26b117f1b89b7d53049012a3da947c422f760e32fe257208` |
| VideoTests executable | `5aa80e7800e2bb8207563627e31f33cb924616e1d923e45efe915bc1c393c8e1` |
| Video Qualification executable | `3c881134238db1b2335b19ff5b2e57cb27879d76a32d8dfb4259e52be47af82f` |
| Fixture manifest | `398bb706268afea61b33529de55f6eda30bb61499495a47bcae3c2b7c24814ce` |
| Original 51_SyncTest.mov | `0870d5ab93d1492f23524ed1fff9319af1a6a9cf71882543f2190929afeb8e24` |

Full executable/source hashes are in the committed
[candidate identity](../../Tests/Video/results/prores-2026-10-20/identity/candidate-04.json).
Golem's transfer verification checks all 286 package artifacts and fixture identities.

## Colour and unsupported media

Exactly one explicit `nclc` or limited-range `nclx` colour declaration must specify
primaries/transfer/matrix 1/1/1, and effective AVFoundation tags must agree. Missing,
duplicate, contradictory, HDR or full-range declarations are rejected with an
explanation. Decoded pixel-buffer tags, when present, must also agree. One-field
progressive metadata, square pixels and full centred aperture are required.
No silent codec-based Rec.709 assumption is made; H.264 policy remains unchanged.

## Native, isolation and product UI results

All 53 final native cases pass per platform. Coverage includes the fifteen required
combinations; exact timestamps/coverage; sequential independent frame identity;
random seeks; transport generations; queue bounds; VFR; multiple streams;
RGB/grey conversion; replacement/clear lifecycle; and colour/field/aperture/aspect/
codec rejection. The actual movie passes 9709 native checks, including all-frame
decoding, random/sequential pixel agreement and timecode.

The original movie remains read-only at
`/Users/fuzz/Developer/LevelScope/Tests/51_SyncTest.mov`, outside Git:
375,714,116 bytes, standard 422, 1080p25, 2772 frames/110.88 seconds, PCM audio and
timecode origin 00:00:00:00. Both platforms also pass 100 seconds of actual playback.

Independent isolation uses nonzero selected 337/613 Hz audio captured through an
injected sink before physical output. The captured peak is 0.00312296 and all 2048
compared samples remain bit-identical after movie attachment, seek, stream switch
and clear. Generated movies contain distinctive 1973 Hz embedded audio; it does
not enter the selected output. The actual movie is also checked. Runtime observation
records video reader outputs and zero audio reader outputs. Physical output count
is zero in this test; source identity, routing and independence checks pass.

Product UI observations on both Macs cover the actual picture and exact hash,
frame/timecode display, offset −1 selecting frame 25/PTS 1/TC 1, offset +1 blanking
outside coverage, clear, readable missing-colour rejection, supported-file recovery
and LT stream 2. Devices remained closed and monitor gain stayed −30 dB. Sequoia
also checks invalid offset text. These are agent observations, separate from timing
statistics and from subjective listening. Private desktop captures are not in Git.

## Final playback matrix

Each platform completes 15 generated 18-second runs plus the actual movie for 100
seconds. Counts start after the predeclared two-second preroll. All reports show
zero actual-paint skipped frames, blanks, bound failures, paint-deadline failures
and audio underruns.

| Platform | Runs | Actual skipped frames | Max PTS skew (ms) | Max CPU paint (ms) | Max sampled RSS (MiB) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Sequoia | 16/16 | 0 | 0.000 | 0.112 | 541.69 |
| High Sierra | 16/16 | 0 | 0.000 | 0.058 | 143.21 |

The independent verifier checks fixture/audio/executable hashes, native architecture,
named device/sample rate, hardware presentation-clock use, full duration and unique
painted-frame coverage. A one-frame endpoint allowance accounts only for termination
phase; the allowance for skipped interior frames is zero.

## Six final endurance runs

Each profile runs 1800 seconds at 1080p30000/1001 on each Mac, using final candidate
04 and the same fixture identities. All six have zero blanks, timing-bound failures,
paint-deadline failures and audio underruns, and pass attachment and memory limits.
Growth is last-minute mean minus first-minute mean, independently recomputed from
raw coherent VM snapshots. All VM queries and required sample counts pass.

| Platform | Profile | Audio-clock seconds | Unique painted frames | Skips | Peak RSS MiB | Non-reusable growth MiB | Physical-footprint growth MiB |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Sequoia | PROXY | 1800.030 | 53887 | 0 | 684.73 | -7.522 | +2.572 |
| Sequoia | LT | 1800.041 | 53887 | 0 | 676.34 | -10.519 | +1.433 |
| Sequoia | 422 | 1800.039 | 53887 | 0 | 612.36 | -5.979 | -2.396 |
| High Sierra | PROXY | 1800.032 | 53887 | 0 | 136.91 | +0.967 | +0.539 |
| High Sierra | LT | 1800.032 | 53887 | 0 | 144.22 | +1.950 | +1.652 |
| High Sierra | 422 | 1800.069 | 53887 | 0 | 149.73 | +1.354 | +1.099 |

| Platform | Profile | Max PTS skew ms | First-minute mean skew ms | Last-minute mean skew ms | Max CPU paint ms | First picture s |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Sequoia | PROXY | 0.000 | 0.000 | 0.000 | 2.626 | 54.721 |
| Sequoia | LT | 0.000 | 0.000 | 0.000 | 1.722 | 126.282 |
| Sequoia | 422 | 0.000 | 0.000 | 0.000 | 0.228 | 186.593 |
| High Sierra | PROXY | 0.000 | 0.000 | 0.000 | 0.254 | 87.660 |
| High Sierra | LT | 0.000 | 0.000 | 0.000 | 0.309 | 194.884 |
| High Sierra | 422 | 0.000 | 0.000 | 0.000 | 0.594 | 283.310 |

The memory gate retains total sampled RSS ≤768 MiB and requires both non-reusable
resident and physical-footprint mean growth ≤32 MiB. First-window seconds are 2..61;
last-window seconds are strictly later than duration minus 60, with at least 55
valid samples in each. No allocator relief or artificial memory pressure is used.
Eleven independent verifier checks confirm rejection of missing/invalid
samples, falsified rollups, excessive retained/footprint growth and excessive raw
RSS even when it is reusable. An additional audit checks every final report’s
OS-reported process RSS high-water mark, including loading, against the same
768 MiB ceiling; all pass. The endurance table reports this conservative peak.

## Large-file loading and cancellation

The fully materialized standard-422 file is 33,008,769,739 bytes. The predeclared
size-scaled first-picture budget is 527.269 seconds; both Macs pass.

| Measurement | Sequoia | High Sierra |
| --- | ---: | ---: |
| Whole-file hash | 171.332 s | 264.735 s |
| Index | 13.154 s | 13.017 s |
| First picture, including hash/index | 184.637 s | 277.960 s |
| Hash cancellation join | 0.000190 s | 0.000227 s |
| Index cancellation | 0.000502 s | 0.001718 s |
| Clear/destruction join | 0.001364 s | 0.000965 s |
| Sampled native loading peak RSS | 40.219 MiB | 29.926 MiB |

Cancellation is cooperative on responsive local storage; no hard deadline for a
blocked operating-system read is claimed. OS cache state was uncontrolled.
Generated endurance files contain 54,000 exact-PTS frames over 1801.8 seconds,
with a repeating 256-frame compressed picture cycle, detailed noise/ramp content,
colour patches and identity bits. Proxy is 10,329,269,771 bytes and LT is
23,155,290,283 bytes. Repeating content is a disclosed limitation, not natural footage.

## Regression evidence and exact reuse

M1–M5, scaffold, H.264 and official ITU/EBU regressions pass on both actual systems;
cross-platform numeric and performance comparisons pass. Sequoia candidate-03 native,
regression and UI evidence is reused only because all production source and all
eight relevant executables are byte-identical in candidate 04. Final playback and
all six endurance runs are fresh candidate-04 results; Golem's full checks are fresh.

The ordinary M5 comparator reports two package identity/source-manifest mismatches.
Its raw output remains preserved. A separate proof validates every artifact, limits
source differences to the unused playback harness/runner/verifier, and accepts
only those two identity differences while retaining every numeric/performance check.
No historical report is rewritten to pretend its package identity changed.

## Earlier failures and approved revision

- Initial cropped-negative fixtures lacked a crop declaration and duplicated valid
  seed bytes. Their three failures remain recorded. Replacement `*-cropped-v2.mov`
  files explicitly declare a 1918×1080 aperture and are correctly rejected on both
  Macs. Original files/manifests remain preserved as invalid test inputs.
- Candidate 01's Sequoia Proxy 29.97 matrix run skipped one painted frame. Its
  cause was not proven. Candidate 02 moved memory sampling off the UI thread and
  added paint-gap diagnostics; both later full matrices passed.
- Candidate 02's Sequoia Proxy 1800-second run had zero timing failures but failed
  the original RSS-growth gate: approximately 120.3 MiB growth against a 32 MiB limit. Its LT run was stopped
  for investigation and 422 was not run. Golem candidate 02 Proxy passed but does
  not substitute for a final endurance result.
- Candidate 03's 300-second stabilization also failed the original RSS-growth
  gate: first mean 529,670,690, last 677,251,801, peak 686,243,840 bytes. Native image
  backing reduced paint copy cost but alone did not resolve that original verdict.
- Diagnostic-only runs found stable live allocations/physical footprint. In the
  180-second diagnostic RSS grew 48.964 MiB while reusable pages grew 50.047 MiB;
  footprint fell 1.266 MiB and live heap fell 4.900 MiB. A diagnostic allocator
  relief request released zero bytes. These experiments are ineligible evidence.
- The user explicitly approved the revised memory gate on 20 October. All earlier
  failures remain failures under their original method. Final candidate 04 contains
  no relief calls and passes fresh stabilization and all six fresh endurance runs.

## Preserved locations and limitations

Source: `/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation`,
branch `codex/video-implementation`.
Local build/candidates/fixtures/raw evidence:
`/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/prores-20261020`.
Golem: `/Users/fuzz/LevelScope-ProRes-qualification-2026-10-20`.
Candidates 01–04 and diagnostic packages remain local; Golem preserves 01, 02 and 04.
Candidate 03 was local only. Historical caches are evidence, not reusable build
destinations. Subsequent work needs a fresh build directory and new identities.

All current hardware runs use Mac Studio Speakers at 48 kHz or Golem Built-in
Output at 44.1 kHz, with every route muted. Internal audio-clock/paint results do
not measure photons, acoustic latency or physical optical/acoustic synchronization.
No new subjective listening claim is made. The preview is 8-bit reference picture,
not a colour-grading or display-calibration claim.

4K/higher rates remain outside the required support claim. Optional soundtrack
playback, VLC/LibVLC, additional codecs, ProRes HQ/4444/RAW, HDR/log conversion,
interlacing, MXF, export/transcoding and VST3 parity remain outside this milestone.
Existing checkpoint tags, installations, historical Golem copies, Pro40 and system
audio configuration are preserved. Media, binaries, private captures and credentials
are excluded from Git. The original logic branch and tracked files are preserved; unrelated untracked files remain in place.

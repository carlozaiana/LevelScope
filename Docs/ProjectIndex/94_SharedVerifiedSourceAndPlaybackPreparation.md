# Session-wide verified source and background playback preparation

This development starts at qualified main `18adeb520006b7fe7e0a22300cc26abd28eec51c`. The reference production implementation is `79b55eb177a01e56fdbe86f843e03140df356bcd`. The shared-source service and cross-session export guard are at `e54ee83be3d0fb9d697e0b4a40e9a4a1d15d6a91`.

## Candidate 03 qualified on both Macs, 24 November 2026

Use [LevelScope Standalone.app](</Users/fuzz/Developer/LevelScope/build/session-source-2026-11-23/candidate-03/LevelScope Standalone.app>) for personal testing. The same universal package passed all 11 focused groups and 227,322 assertions on each actual Mac. Sequoia 15.7.9 (24G830) ran arm64 and Golem ran High Sierra 10.13.6 (17G66) x86_64. Each host verified all 477 manifest entries before and after execution. The rational clock cases account for most of the assertions. [Sequoia receipt](../../Tests/Qualification/SharedSource/results/combined-2026-11-24/sequoia/receipt.json) · [Golem receipt](../../Tests/Qualification/SharedSource/results/combined-2026-11-24/high-sierra/receipt.json)

| Identity | Value |
| --- | --- |
| Tested implementation | `614377f4e573bdd0fc34726584418b9333d092c2` |
| Source tree | `faec0b86d24f487ec10eb127d17ee05281688333` |
| Measurement build | `a4ff98da4b6d20c2ad8e52492db006b7bb7c529b224515f9b48f5ebc220ac20a` |
| Package identity SHA-256 | `94c42f1eca24de3adc671723b4b01cd4bcc986fe602b4c6c6d25484379079843` |
| App executable SHA-256 | `a6f13de9d818b31bd3a72ac19d1a6500241b772b488219c53c6343d9b9e5a91d` |

Candidate 02 missed a same-size mutation of its private source on Golem’s Journaled HFS+ volume. Two writes within one filesystem timestamp interval can have identical timestamps. Candidate 03 watches the original and its protected copy for kernel file-change notifications. The added test overwrites a file through an already-open descriptor, matches all the timestamp fields and checks that both consumers reject the old version. All 289 source-service assertions passed on each host. [Candidate 02 failure](../../Tests/Qualification/SharedSource/results/combined-2026-11-24/candidate-02-high-sierra/service.json) · [Candidate 03 integrity checks](../../Tests/Qualification/SharedSource/results/combined-2026-11-24/high-sierra/service.json)

| Host and programme | First Play, reference → candidate 03 | Subsequent Play, reference → candidate 03 | Background preparation |
| --- | ---: | ---: | ---: |
| Sequoia stereo | 912.45 → 1.31 ms | 2.61 → 3.49 ms | 2.62 ms |
| Sequoia 5.1 | 1,844.28 → 1.51 ms | 2.63 → 4.00 ms | 2.65 ms |
| Golem stereo | 1,613.22 → 1.76 ms | 3.15 → 4.85 ms | 3.22 ms |
| Golem 5.1 | 3,118.76 → 1.80 ms | 3.14 → 4.87 ms | 3.22 ms |

These are single elapsed observations from an explicit Play request to the first nonzero sample at injected output with 1 ms engine polling. Physical output and native UI scheduling fall outside the timing measurement. Each session made one verified input copy during initial analysis. Background preparation, first and subsequent Play, and two subsequent proposal searches each added zero full source copies. Repeated Propose on Golem spent 0.079 ms per search preparing its input, previously 2.31 s and 1.56 s. The complete searches took 17.93 s and 17.94 s. [Timing and proposal comparisons](../../Tests/Qualification/SharedSource/results/combined-2026-11-24/summary.json)

On each host, Source, every candidate receipt and selection, independently measured Current, rendered-output metrics and exported WAV bytes exactly match the qualified reference. Across hosts, Source, candidate plans and receipts, Current and the exported file also match exactly. Of 16,480 fields in the rendered measurement, only channel 2 true peak differs. Its linear amplitudes were 0.8911029454659865 on Sequoia and 0.8911029454659866 on Golem, one binary64 step apart. All peak positions and programme metrics match. The same-platform reference digests reproduce this difference exactly. The exported file SHA-256 on both Macs is `a06057a7bbec3a36889f814763a19f5d33d52963716de4ba4e27186687cbbb11`. [Full-precision comparison](../../Tests/Qualification/SharedSource/results/combined-2026-11-24/cross-platform.json)

The development branch is `codex/session-verified-source-2026-11-23`. Both manuals, all 28 pre-existing tag objects and targets, and all 3,638 protected files across 11 earlier packages were preserved. Main and tracking main were `18adeb520006b7fe7e0a22300cc26abd28eec51c` at handover. The candidate 02 failure and both frozen packages are retained. [Package identity](../../Tests/Qualification/SharedSource/results/combined-2026-11-24/package/identity.json) · [Preservation audit](../../Tests/Qualification/SharedSource/results/combined-2026-11-24/preservation.json)

## Source lifetime

The live session holds one protected file for each verified physical audio asset. The first consumer streams a copy to private disk, validates the audio and compares the complete copied-file SHA-256 with the imported content identity. An extra whole-file hash previously performed immediately before that validation was removed. Import and session opening still verify the referenced media.

Source measurement, Current measurement, waveform rebuilding, proposals, playback and export rendering each open independent readers. Separate-file ordering and composed-programme maps are reconstructed from the current interpretation. Requirements, gains, roles and region boundaries retain their existing result and processing bindings. The existing initial-unity Current adoption remains available. Subsequent candidate measurements, independent applied Current verification and the encoded export measurement each run on their own output. Final delivery retains its original-file and output-hash checks.

The content identity comes from the verified copy. On macOS, the service opens a kernel vnode watch before validating each original and protected copy. The runtime version combines pathname identity, size and timestamps with this watch. A write or metadata event latches invalidation for every consumer of the version. Releasing the last reference closes its descriptors. A changed or missing original blocks use and requires media resolution. A missing or altered private copy is verified anew from the original before reuse. Original media stays read-only to the application. An explicit older export lease can only be reused when the receiving store owns the corresponding protected files.

Concurrent acquisitions of the same asset have one producer. A waiting consumer can cancel independently. If the copying consumer cancels, its partial file is discarded and a remaining consumer can establish trust. New, Open and replacement get a new store. A live recipe edit can retain already verified files for unchanged members. A worker or reader can finish against its own lease without granting a newer session authority. The last owner removes its private file and leaves unrelated replacement files alone. Opening a saved project establishes fresh trust in its referenced media.

## Playback

After initial Original analysis and its progressive overview finish, the main window automatically prepares the selected Source or Current reader, processor and a prefix of up to 4,096 frames. Readiness is visible with the outputs closed. At a stopped seek the prefix belongs to that exact frame. A gain, boundary, domain, runtime or source-owner change expires an incompatible prefix. Existing explicit output connection, routing, mute and volume choices retain their own controls.

A Play request can join preparation. Stop or Cancel revokes delayed playback. Starting output checks the live session, routing and device. The reference-picture gate still controls callback release. A prepared reader does not grant Current a numerical PASS. Physical playback and listening are user actions.

## Qualification boundary

The focused clock and transport group uses synthetic audio, injected output and no movie. Broad video, full-film and three-hour endurance remain deferred without a 90-minute substitute. The previously recorded completion pause remains a separate limitation. These results make no new listening or real-device latency claim. The proposal engine, rules, search order, supported layouts and margin defaults retain their existing behavior.

## Earlier candidate 02 on Sequoia

The preserved candidate 02 build is [LevelScope Standalone.app](</Users/fuzz/Developer/LevelScope/build/session-source-2026-11-23/candidate-02/LevelScope Standalone.app>). Its preserved directory contains the app, universal test executables, source snapshot, build commands and the [package identity](../../Tests/Qualification/SharedSource/results/sequoia-02-2026-11-23/package/identity.json).

| Identity | Value |
| --- | --- |
| Tested implementation | `958871a9f74f3efc7dcac80249490ad127c6a672` |
| Git Source tree | `1067f6ee334361493ba7c627064476563a756d37` |
| Measurement build | `3fc9ebd577afa62dc62fc04e12ef52d1abc03042379852d0ebe6d17120e25e86` |
| Package identity SHA-256 | `bf6ce039c0c8ccde854e359e121fe30bb1a959acd959bb8044247f6ffcafd101` |
| App executable SHA-256 | `18093833439afb5159fb197f7c14da20e00a406e34e836894bc8b34aa7862678` |
| Source and tests manifest SHA-256 | `12e4c1ce783d142840f2d8bf7fa663bfa1b947f506e11e095b10ae4605356cf6` |

All 11 groups passed on macOS 15.7.9 (24G830), arm64. The runner verified all 476 manifest entries before and after execution. The Intel slice targets 10.13. The independent host record is in the [Sequoia receipt](../../Tests/Qualification/SharedSource/results/sequoia-02-2026-11-23/receipt.json).

| Group | Passing assertions |
| --- | ---: |
| service | 283 |
| feedback | 940 |
| audition | 4,069 |
| devices | 450 |
| timeline | 2,580 |
| clock | 218,104 |
| ordinary-ui | 167 |
| playback-ui | 668 |
| propose | 17 |
| timing-stereo | 19 |
| timing-surround | 19 |
| **Total** | **227,316** |

Most clock assertions exercise rational time arithmetic. The synthetic transport cases check exact Source and Current samples and the played frame, explicit Stop, live and stopped seeks, device loss, routing or source replacement, and the reference-picture preroll gate. Candidate 01 failed on live seek. The replacement keeps one active worker and remembers the latest requested frame. A newer request supersedes the earlier one, and Stop revokes any delayed resume. Both the first package and its failed receipt are preserved.

## Earlier candidate 02 timings on Sequoia

Each Play observation measures elapsed time from an explicit request until the first nonzero sample at the injected output. The same input and development Mac were used with one-millisecond engine polling. Main-window timer scheduling, device buffers and physical audio latency lie outside this measurement. Durations are single observations for each case.

| Input | First Play before | First Play after background preparation | Background preparation | Subsequent Play before | Subsequent Play after | Extra bytes copied for Play, before → after |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Music/dialogue stereo | 912.45 ms | 1.35 ms | 2.58 ms | 2.61 ms | 3.73 ms | 79,375,076 → 0 |
| Busy 5.1 | 1844.28 ms | 1.37 ms | 2.58 ms | 2.63 ms | 3.71 ms | 160,617,696 → 0 |

Initial Source analysis established one verified file per interleaved input. Background preparation, first Play and subsequent Play shared that file. The four-asset composed fixture also reused its four files for waveform generation, Current measurement, render and finalization. Reopening established independent trust.

| Stereo Propose run | Source preparation before | Source preparation after | Full search before | Full search after | Additional source copies before → after |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 918.42 ms | 0.052 ms | 10.79 s | 9.75 s | 1 → 0 |
| 2 | 909.18 ms | 0.057 ms | 10.88 s | 9.86 s | 1 → 0 |

Each qualified baseline proposal copied and validated 79,375,076 bytes. Both new runs recorded zero full source copies, validations or validation hash passes. Each still assessed three candidate plans. Attempt order, plan boundaries and gains, full numerical payloads, selected result, independently measured applied Current and final encoded export matched the qualified implementation exactly. The rendered file's SHA-256 was `a06057a7bbec3a36889f814763a19f5d33d52963716de4ba4e27186687cbbb11`. See the [proposal comparison](../../Tests/Qualification/SharedSource/results/sequoia-02-2026-11-23/proposal-parity.json) and [playback comparison](../../Tests/Qualification/SharedSource/results/sequoia-02-2026-11-23/playback-comparison.json).

The [preservation audit](../../Tests/Qualification/SharedSource/results/sequoia-02-2026-11-23/preservation.json) checked 3,638 protected files across 11 earlier packages and the objects and targets of all 28 existing checkpoint tags. Main remains at `18adeb520006b7fe7e0a22300cc26abd28eec51c`. Development and handover are local.

## Try the build

1. Open candidate 03 from the link at the top and import the existing music/dialogue stereo WAV, confirming its channel roles. Let Original analysis and the progressive waveform finish. Look for **Playback ready** beside the transport position.
2. Connect your usual output if needed and press Play. Stop, seek to a different passage, then press Play again. The stopped cursor can move without starting audio.
3. Switch to Current, apply a gain or boundary edit and play the passage again. Run Propose with your usual settings. Current still requires its independent verification.
4. To check cancellation, press Play while preparation is visible and press Stop before it completes. Later completion must leave transport stopped.

## Golem run provenance

Golem received candidate 03 in `/Users/fuzz/LevelScope-shared-source-20261123`. The input bundle, qualified reference drivers and two WAVs matched their expected checksums. The corrected package archive SHA-256 was `6b390d1a3547402c4bc7573abc953f13b2af941239ec1db2305d180d44c37d0a`. The returned receipt archive SHA-256 was `d32e335f1a781fe55e2a1c7df69440c31e7c3603f31ffcf0a45a71ec5da1bb89`. The same finite runner performed the earlier candidate 02 run (`high-sierra-01`) and the final candidate 03 run (`high-sierra-02`). [Host identity](../../Tests/Qualification/SharedSource/results/combined-2026-11-24/host.txt)

```sh
/usr/local/bin/python3 qualify_host.py candidate-03 94c42f1eca24de3adc671723b4b01cd4bcc986fe602b4c6c6d25484379079843 high-sierra-02 \
  --benchmark "$PWD/inputs/stereo.wav" --timing "$PWD/inputs/stereo.wav" "$PWD/inputs/surround.wav"
```

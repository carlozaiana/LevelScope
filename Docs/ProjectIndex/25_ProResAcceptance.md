# Bounded SDR ProRes acceptance — 20 October 2026

Status: implementation and qualification in progress; no acceptance claim yet.
Base: published main 6d304dbecd8fa47f56b0260105aa9d07f96ba1bc (includes qualified M5).
The prior video checkpoint and all other tags, installed products and historical evidence remain immutable.

## Required scope

MOV; ProRes 422 Proxy, LT and standard 422; progressive, square pixels, full aperture;
Rec.709 SDR; existing 8-bit reference preview. Required real-time combinations:
1920x1080 at 24000/1001, 24, 25, 30000/1001 and 30 fps on actual Sequoia arm64
and High Sierra x86_64. 3840x2160 and 60000/1001 probes are exploratory only.
Selected audio owns the clock. Embedded audio remains undecoded and unrouted.

## Original acceptance limits (defined before qualification)

- After the existing two-second preroll, **zero skipped source frames at actual UI
  paints** for required combinations. At <=30 fps on the existing 60 Hz displays
  every source frame has at least two refresh opportunities; no frame dropping is
  an intended operating mode. Timer/worker selection skips are recorded separately
  and cannot substitute for the actual-paint count.
- Zero blank paints, one-source-frame-bound violations, paint-duration deadline
  failures and audio underruns. Record first/last-minute and maximum timing skew.
  Complete 1800 seconds per profile/platform using the same final candidate.
- Healthy local-storage attachment budget: 30 seconds + file bytes / 64 MiB/s +
  10 seconds per 100,000 indexed frames. This explicitly includes whole-file hashing
  and indexing and permits conservative sustained hard-disk throughput. Report
  hash, index and first-picture time separately. Large fixture must exceed 8 GiB.
- Cooperative cancellation completes within two seconds on responsive local storage;
  clearing the reference removes its picture immediately and stale work never returns.
  No hard OS-read timeout is claimed for unavailable/failing storage.
- Peak RSS <=768 MiB for required 1080p playback; after warmup, last-minute mean RSS
  must be within 32 MiB of first-minute mean. The bounded queue remains <=8 images.
  Report attachment peak separately, including any temporary double indexing in tests.
  **The RSS growth subcriterion is superseded by the explicit approved revision below.**
  Earlier runs retain their original verdicts; the peak ceiling is unchanged.
- Picture tests: dominant RGB patch component >210, other components <40; neutral
  128 patch within 12 eight-bit levels in every channel. These allow lossy Proxy
  and decoder rounding, while detecting channel/range errors. Independent burn-in
  identity must match exactly; actual-file sequential/seek pixel samples must match.

## Explicit memory-method revision — approved 20 October 2026

The user approved this correction after the Sequoia diagnostic demonstrated that
RSS growth was accounted for by OS-reusable pages. In the 180-second diagnostic,
first/last-minute RSS increased 48.964 MiB and reusable pages increased 50.047 MiB;
physical footprint decreased 1.266 MiB and live heap decreased 4.900 MiB.
The original failed reports are preserved and are **not** retroactively accepted.
See [26_ProResQualification.md](26_ProResQualification.md) and the raw memory-investigation evidence.

For fresh final-candidate qualification, retain total sampled playback RSS <=768 MiB.
Require last-minute mean minus first-minute mean <=32 MiB **independently for both**:

- non-reusable resident memory: `TASK_VM_INFO.resident_size - reusable`;
- physical footprint: `TASK_VM_INFO.phys_footprint`, including the OS accounting
  for compressed/private memory and relevant mappings.

One coherent VM snapshot per playback second is collected on a background thread.
Record every raw sample, total RSS and reusable pages as well as both growth metrics.
The growth windows remain seconds 2..61 and strictly later than duration minus 60;
require at least 55 valid samples in each window for runs >=120 seconds. Any VM
query failure, incomplete snapshot, invalid counters or insufficient samples fails
the gate. Short matrix runs retain the total-RSS ceiling; the 1800-second runs
assess both growth limits. The queue bound and all timing/skip gates are unchanged.
No allocator pressure relief, forced cache purge or artificial memory pressure is
requested in qualification or production. The diagnostic-only allocator experiments
remain separately labelled and ineligible. All six endurance runs must be fresh.

This distinguishes reclaimable resident pages from memory still held by the process;
it does not remove raw RSS reporting or the original total-memory ceiling. Apple's
[task ledger definitions](https://github.com/apple-oss-distributions/xnu/blob/xnu-11417.140.69/osfmk/kern/task.c#L1175-L1213)
and [reusable-page reclamation](https://github.com/apple-oss-distributions/xnu/blob/xnu-11417.140.69/osfmk/vm/vm_resident.c#L6762-L6780)
provide the accounting context; the actual two-platform runtime samples decide acceptance.

## Colour interpretation

ProRes requires one explicit MOV sample-description `colr` atom: nclc/nclx with
Rec.709 primaries/transfer/matrix (1/1/1); nclx full-range declarations are outside
scope. All three effective AVFoundation format-description tags must also agree
with Rec.709. Missing, duplicate, conflicting, HDR or other tags are rejected with
an explanatory error. This prevents OS-supplied defaults from silently admitting
untagged material. Any contradictory decoded pixel-buffer tags also fail explicitly.
Explicit one-field/progressive metadata is required. H.264's existing policy is unchanged.

## Evidence and regression requirements

Include the original /Users/fuzz/Developer/LevelScope/Tests/51_SyncTest.mov read-only:
SHA256 0870d5ab93d1492f23524ed1fff9319af1a6a9cf71882543f2190929afeb8e24,
375714116 bytes, ProRes 422, 1920x1080, 25 fps, 2772 frames/110.88 s.
Keep it out of Git. Generated detailed fixtures, their manifests and report hashes
bind all observed results to exact source, binaries, OS, architecture and media.
Use independently nonzero selected audio captured through the injected test sink;
compare all captured samples with the baseline, with a distinct embedded tone present.
A test-only Objective-C observer counts actual AVAssetReaderTrackOutput constructions.
Do not treat silence or a hard-coded output counter as isolation evidence.

Run native identity/seeking/coverage/timecode/VFR/lifecycle checks, real product UI
checks, H.264 regressions and current M1–M5 tests plus official ITU/EBU comparisons.
Keep internal scheduler results, user listening and physical synchronization separate.
Physical optical/acoustic latency remains unmeasured.

Use only fresh build/prores-20261020 products in the VideoImplementation worktree.
No VST3 build/install, VLC, additional codecs, HDR/log conversion, interlacing, MXF,
export/transcoding or embedded-soundtrack playback is part of this milestone.

# M18 development Mac offline job profile

The first implementation I recommend is a bounded per-channel true-peak prototype. The same peak calculation dominates Source measurement, independent Current remeasurement and all three private section passes. Keep each channel's filter history continuous across blocks and reduce the loudness contributions in the existing channel order. Compare all scalar results, windows, peak positions, tie handling and proposal decisions with the serial engine before adopting the prototype.

## Conditions

1 November 2026. Apple M1 Max, 10 physical cores, 32 GiB memory, macOS 15.7.9 (24G830), native arm64 Release. Both fixtures are 60 seconds of local, interleaved 48 kHz 24-bit PCM with 20 sections. The 5.1 order is L, R, C, LFE, Ls, Rs. Each run starts in a fresh process. The OS file cache is warm. The fresh section proposal finds one local correction and measures attribution, local and combined plans. A second proposal reuses its valid attribution. The harness applies that result and independently remeasures Current. Initial untouched Current uses Source's completed measurement.

One warm-up per layout precedes three timed repetitions of each executable per layout. Control runs use the unmodified Release service objects. A companion executable times stage boundaries. Each corresponding Source, Current, candidate, proposal and plan identity matched across both builds and all repetitions. For Source, Current and both section proposal jobs, instrumented median times were within 1% of the corresponding control medians. The numeric build identity, binary hashes and individual measurements are in [summary.json](summary.json). The [harness and reproduction instructions](../../../Profiling/OfflineStages/README.md) describe the fixture, call boundaries and timer.

## Time per operation

Seconds, median with the observed minimum and maximum across the three control runs. The import is separate. Source and Current include verification, PCM decoding, progressive waveforms and measured loudness. Proposal entries include job execution, synchronous global preflight and acceptance. Window polling and subsequent UI updates are outside these totals.

| Operation | Stereo | 5.1 |
|---|---:|---:|
| Initial import validation | 0.101 (0.101–0.103) | 0.304 (0.303–0.306) |
| Source analysis | 1.397 (1.389–1.398) | 4.192 (4.158–4.192) |
| Global-only Propose, 0 passes | 0.003 (0.003–0.003) | 0.003 (0.003–0.003) |
| Fresh section Propose, 3 passes | 4.062 (4.054–4.095) | 12.340 (12.312–12.549) |
| Repeated section Propose, 2 passes | 2.605 (2.601–2.630) | 7.799 (7.721–8.033) |
| Re-measure Current after Apply | 1.416 (1.402–1.433) | 4.200 (4.169–4.245) |


## Source analysis and verification

Seconds, medians from the instrumented runs. The copy, first hash, payload scan and validation hash are nested inside verified-copy preparation. Reading blocks and applying gain occur after preparation. The meter includes all channels' filters, peaks and loudness windows. Setup, small progress updates and cleanup account for the rest.

| Operation | Stereo | 5.1 |
|---|---:|---:|
| Verified-copy preparation, inclusive | 0.198 | 0.595 |
| Copy and flush | 0.006 | 0.015 |
| First SHA-256 scan of the copy | 0.091 | 0.276 |
| Validation decode and finite-sample scan | 0.009 | 0.027 |
| Validation SHA-256 scan of the same copy | 0.092 | 0.276 |
| Read programme blocks for analysis | 0.009 | 0.024 |
| Apply processing | 0.003 | 0.009 |
| Build per-channel waveform | 0.006 | 0.018 |
| Feed the meter | 1.178 | 3.529 |
| Finalize meter and gates | 0.001 | 0.002 |

Feeding the meter takes 84.4% of the stereo Source job and 84.4% of the 5.1 job. The Source worker uses 1.00 and 1.00 core equivalents respectively, from its own CPU time divided by elapsed time. Hashing the same copy twice accounts for about 92–93% of verification. One of those scans accounts for 6.5% of stereo Source analysis, 6.6% of 5.1 Source analysis and 2.2% of fresh section Propose in both layouts. These fractions describe the measured cost of the existing scan.

## Where fresh section Propose spends its time

Seconds, medians from the instrumented runs. Each pass includes its own block reads, processing, measurement and finalization. A single verified copy is shared across the passes of one proposal.

| Operation | Stereo | 5.1 |
|---|---:|---:|
| Prepare verified copy | 0.199 | 0.597 |
| Measure baseline and collect section peaks | 1.462 | 4.601 |
| Measure proposed local gains | 1.209 | 3.581 |
| Measure local and global gains together | 1.207 | 3.614 |

The global-only case uses a separately selected, attainable IL target (−23.7 LUFS stereo and −24.0 LUFS 5.1) with the same −1 dBTP ceiling. The section cases use −23 LUFS and −1 dBTP. Since local candidates and combined gain depend on prior results, a second section proposal saves a pass only when its measured Current snapshot and annotation geometry still match.

## Where the worker executes

Independent 3 ms macOS stack samples were collected in separate uninstrumented runs. An initial sampler missed the short Source worker while loading symbols. The reported samples use a separate uninstrumented launcher with a startup delay outside every timed job. Its measured results and plans matched the control. Counts are inclusive shares of stacks for the named worker throughout that job, including verification. The true-peak routine handles sample peaks, interpolation and observer callbacks. The collector percentage is nested in the true-peak routine. Stack sampling identifies hot code but does not distinguish individual instructions within that routine.

| Operation | Stereo | 5.1 |
|---|---:|---:|
| Source worker in addPeak | 78.7% | 80.5% |
| Fresh Propose worker in addPeak | 87.4% | 89.0% |
| Fresh Propose worker in the section true-peak collector | 6.8% | 7.5% |
| Re-measure Current worker in addPeak | 79.8% | 80.0% |

## Timeline responsiveness and scope

The main JUCE thread refreshed proposal controls and the gain lane and rendered an offscreen timeline at a requested 100 ms interval during the background jobs. Each figure below pools the three uninstrumented runs for that layout.

| Operation | Stereo | 5.1 |
|---|---:|---:|
| Refresh interval, p95 | 110.2 ms | 109.6 ms |
| Longest interval including end of job | 110.8 ms | 110.8 ms |
| Render time, p95 | 4.9 ms | 8.2 ms |

These are short development timings with local PCM and a warm filesystem cache. Cold or external storage, compressed audio, live playback and full-window input latency are different workloads. Actual High Sierra, video, cross-Mac and full-duration qualification continue on the existing schedule.

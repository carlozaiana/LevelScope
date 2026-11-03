# Bounded true-peak development gate

The prototype uses the existing peak calculation and interpolation coefficients.
A synchronous dispatch loop divides channels across at most six lanes, capped
by the channel count and the reported logical core count minus one. It borrows
one decoded input block and joins before that block can be reused. Each channel
owns its peak history and maxima. The original frame-major loudness reduction,
short tails and EOF drain stay on the calling worker. The service gate targets
Apple-silicon stereo and 5.1. Other channel layouts, unknown core counts, and
observers that require serial callbacks use the serial path.

The first candidate passed ordinary measurements and seven proposal comparisons
but failed an adversarial comparison with downward rounding after pool threads
had been used. The replacement captures the caller's floating-point environment
per block and restores each system worker's prior environment even on an observer
exception. The core oracle exercises all four rounding modes.

Before enabling the development build, require bit equality of all metric values
and statuses, gate counts and thresholds, window energies and end frames, sample
peaks and true peaks and their positions. Compare per-channel observer events
including ties, interpolation phases and padding. Cover multiple block sizes,
shuffled stereo/5.1 orders at 44.1, 48, 96 and 192 kHz, empty and very short input,
silence, fallback layouts, legacy observer ordering, non-finite samples and an
observer error on a worker. Run a focused ThreadSanitizer audit. Compare complete
proposal explanation and section peak map digests, statuses, Apply availability,
private candidate measurements, processing plans and independently remeasured
Current on the same PCM. Exercise stale inputs, preference rejection, unsupported
correction, measured rejection, global-only, section-ready and already-meets.

Use the previous 60-second PCM24 fixtures and 10 Hz offscreen timeline render
on the M1 Max. Alternate uninstrumented native serial and bounded service runs
with one warm-up per layout/mode and three timed repetitions. Require at least a
10% median elapsed-time improvement on Source, independently measured Current,
fresh section Propose and repeated section Propose in each layout. Global-only
Propose must stay below 10 ms and the p95 heartbeat gap must stay within 20 ms
of its matched serial run, with the longest gap below 200 ms. The timed jobs use
warm local files. Raw reports live under the candidate directories.

Candidate 02 passed the numerical, proposal and timing gates. An extra run of
the frozen serial binary exposed a timing regression in the rebuilt fallback.
Keeping the lane count in one byte of existing padding preserved the original
object size and observer offsets. An isolated stereo comparison restored fresh
serial Propose from about 4.7 s to 4.07 s. Candidate 03 includes this layout fix.
Its additional gate compares serial and bounded medians with a fresh spot
check of the previous hash-verified binary. The fallback may be at most 10%
slower and the bounded path must be at least 10% faster for all four substantial
jobs. The earlier executable matched its previously recorded measurement and
plan identities. Historical checks are in [the frozen controls](candidate-02/prior-control).

Actual High Sierra, cross-Mac, hardware/video and full-duration/RSS qualification
remain in the combined M18/M19 block.

## Final development results

Candidate 03 is the enabled prototype on the development M1 Max (10 cores, 32 GiB,
macOS 15.7.9). The probe used 60-second interleaved 48 kHz PCM24, the previous fixture
hashes, the same 20 sections, one warm-up and three alternating timed repetitions per
layout and execution choice. Values below are median seconds, except for the global
case in milliseconds. The frozen earlier binary was run once per layout as an
independent comparison. Its input, measurement and plan identities matched the
earlier profiling report. Source and Current include verified-copy preparation. Section proposal times measure service execution. The global-only case also includes synchronous preflight and acceptance.

| Job | Serial stereo | Bounded stereo | Serial 5.1 | Bounded 5.1 |
|---|---:|---:|---:|---:|
| Source analysis | 1.39 | 0.85 | 4.18 | 1.34 |
| Fresh section proposal, three passes | 4.07 | 2.33 | 12.19 | 3.97 |
| Repeat proposal, two passes | 2.61 | 1.53 | 7.89 | 2.08 |
| Independent Current remeasurement | 1.41 | 0.87 | 4.18 | 1.35 |
| Global-only proposal (ms) | 2.78 | 2.75 | 2.76 | 2.79 |

The 100 ms JUCE heartbeat had bounded-mode p95 gaps of 105.0 ms and 105.3 ms for stereo and 5.1 respectively. These measurements cover offscreen timeline rendering during offline jobs.

The [paired timings and ranges](candidate-03/timings/summary.json) and
[exact decision comparisons](candidate-03/decisions/passed.json) passed while the
application setting was OFF. The normal development build was then configured
with `LEVELSCOPE_ENABLE_BOUNDED_TRUE_PEAK=ON`. The [enabled-service comparison](candidate-03/enabled/passed.json)
exercised the ordinary service objects at this build identity and checked concurrent
execution, faster elapsed times and complete serial parity in each layout.
The independent section-peak event oracle, cancellation, explicit Apply,
remeasurement and persistence are covered by the enabled measurement and proposal
suites (946 checks), in addition to the core parity and
ThreadSanitizer checks.

Development app: `/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m18-peak-dev/LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app`.
Enabled identity: `93460e4031c2ec51765b4ad83dd7bc4e67551be1f94129b107ebbeea494a32ce`. The app and selected tests are universal with
Intel minimum 10.13. The Intel application uses the serial service. Previous-build
session caches are validated against the new build binding and can require
remeasurement. [Binary manifest](candidate-03/development-build.json).

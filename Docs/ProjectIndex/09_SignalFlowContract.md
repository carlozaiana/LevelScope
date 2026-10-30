# Signal Flow Contract

Purpose: freeze operational order-of-operations and path separation so DSP, UI, and future compliance/proposal systems reference one signal contract.

## Standalone M6 path (qualified candidate 01; publication pending)

The numbered plugin contract below remains separate. Standalone M6 uses an owned,
hash-verified source copy, float32 decoding, and the shared `StaticGainProcessor`
operation `float(double(sample) * pow(10.0, gainDb / 20.0))`. Every ordered source
channel, including LFE, receives the same trim. Source uses zero dB; Current uses
the exact typed M1 plan. There is no plugin graph, limiter, downmix or resampling.

The worker feeds a bounded 12-channel-capacity FIFO. Consumed pre-monitor samples
are the qualification capture point. Monitor smoothing and explicit one-to-one
output routing follow that point and stay outside processing/measurement identity.
An actual routed monitor sample outside [-1, 1] suppresses its entire output
callback and stops playback. Explicit recovery resets FIFO and ramp without
altering processing gain or measurements. Suppressed spans never advance the
audio presentation clock. Signed video offsets and rational source coordinates
retain the accepted video contract; embedded video audio remains excluded.

See [M6 scope, lifecycle and acceptance tests](28_M6VerifiedCurrentStateAudition.md).
Actual two-platform results and exact candidate identity are recorded in
[the qualification handover](30_M6QualificationAndHandover.md).

## 1) End-to-end order of operations

1. Host calls `processBlock()` with input audio (+ optional MIDI).
2. Plugin resolves transport/timebase state and discontinuity flags.
3. Plugin creates `ProcessContext` (sampleRate, numSamples, channel set, transport/discontinuity, timing indices).
4. **Input-side analysis path** updates loudness/history truth from unprocessed input (current behavior).
5. Plugin runs active `ProcessorCore` module graph on audio buffer in order.
6. Output-side snapshots/meters are captured for UI-safe polling.
7. Plugin returns processed audio to host.

## 2) Analysis path vs processing path

### Analysis path (current)
- Source: pre-module-chain input.
- Outputs: timeline history, loudness metrics (M/S/I/LRA variants), analysis-driven UI overlays.
- Contract note: this is currently **input truth**, not post-processed compliance truth.

### Processing path (current)
- Source: same input block, then module graph mutates buffer in-place.
- Composition: module graph owned by `ProcessorCore` (RCU-style swap model, audio thread reads active snapshot).
- Outputs: post-chain audio + stage/output meter snapshots.

### Planned path split
- Future IN snapshot vs OUT live analysis track split remains planned; this contract will expand when OUT analysis is formalized.

## 3) Latency introduction and compensation

### Latency-introducing points
- Limiter lookahead and any stage oversampling/lookahead design options.
- Any future linear-phase or analysis-window buffering modules.

### Compensation contract
- Plugin reports host latency and updates when latency-affecting params change.
- Update trigger mechanism: APVTS listener + message-thread timer policy already integrated.
- UI playback lock: controls that can change structure/latency are blocked during effective playback.

## 4) Transport/discontinuity handling notes

Implemented safety behaviors:
- Transport start ramp guard to avoid startup spikes/artifacts.
- Stop-time silence freeze to avoid corrupting analysis history.
- Discontinuity detection for seek/loop conditions.
- Effective playback state combines transport intent + callback freshness for UI lock decisions.

Operational implication:
- Any module or analysis extension must treat discontinuity boundaries explicitly (reset/slew/freeze behavior) and must not smear pre/post-seek measurement continuity unless intentional.

## 5) Change-control checklist (for future edits)

When changing flow order, latency, or transport semantics, update:
1. this contract,
2. `Docs/ComplianceMatrix.md` status rows tied to timing/transport,
3. validation notes in `08_TestingBuildAndRiskNotes.md` (or successor automated suite docs).

## 6) Historical standalone scaffold source audition (15 October 2026)

The original scaffold decoded Source directly into a bounded SPSC buffer, applied
separate monitor gain/routing, and restarted from frame zero on Play. This is
historical behavior, superseded by the accepted transport and M6 contracts above.
Current standalone audition reads a verified private copy and applies the selected
Source/Current operation before monitoring. Play resumes the held source position;
Restart explicitly begins at zero. Underruns do not skip source frames, and no
automatic downmix or sample-rate conversion is introduced. Current qualification,
including completed listening checks, is recorded in the M6 handover.

## 7) Standalone offline measurement (15 October 2026)

Worker: private byte snapshot → content hash/payload validation → decode in original
file order → identity or explicit all-channel static gain → OfflineMeasurement.
Source measurement requires identity processing. Current State measures actual gain-
processed samples. The internal render harness encodes PCM24, closes/reopens the file
and measures the encoded samples with the same kernel. Full-program measurement uses
the source sample clock and its own reset/window/finalization rules, not plugin host
transport guards or the audition ring. M1 added this path without changing the
plugin callback. M6 shares snapshot/gain primitives with audition; measurement
windowing and finalization remain offline.
See [14_OfflineMeasurementContract.md](14_OfflineMeasurementContract.md).

## M16 proposal preferences remain outside sample processing

The standalone signal remains confirmed programme/component sum → applied canonical
section envelope → global trim → measurement or pre-monitor audition/export path.
M16 selects and describes a candidate using versioned proposal preferences; the
preferences themselves never enter the sample pipeline or measured-audio identity.
Policy edits preserve audio, while pending proposal/report approvals become stale.
Sessions distinguish active policy, historical application policy and exact applied
audio. Independent measured verification and one final export dither remain unchanged.
See ProjectIndex 59 for phase status; no new video/High Sierra qualification is claimed.

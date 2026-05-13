# Signal Flow Contract

Purpose: freeze operational order-of-operations and path separation so DSP, UI, and future compliance/proposal systems reference one signal contract.

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

# DSP / Audio Module Inventory

## Core module host
- `ProcessorCore`: owns/swaps active module graph and drives per-block module processing.

## Modules (IAudioModule)
1. **LevelerModule**
   - Purpose: channel/group level stabilization toward target LUFS with rate limits and optional host-gain lane capture/control.
   - Files: `LevelerModule.h/.cpp`, `LevelerParamIDs.h`.

2. **MultiThresholdDynamicsModule (MTDM)**
   - Purpose: threshold-zone based dynamics structure around T0..T3 with upward/downward/limiter sections, channel-mask policies, zone audition and bypass controls.
   - Files: `MultiThresholdDynamicsModule.h/.cpp`, `MultiThresholdDynamicsParamIDs.h`, `MultiThresholdDynamicsParamGroups.h`.

## DSP building blocks used by modules
- `SpectralUpwardCompressor` (FFT-band upward processing)
- `BroadbandUpwardCompressor` (time-domain upward path)
- `BroadbandDownwardCompressor` (downward compression)
- `LookaheadLimiter` (peak limiting, detector oversampling options)
- `UpwardGainLaw` (shared upward curve math)
- `BS1770MomentaryLufsDetector` + `BS1770KWeighting` (loudness control-signal measurement)

## Supporting analysis/state
- `LevelScopeHistoryModel` for timeline history.
- `RunningLoudnessStats` for running integrated/LRA summaries.
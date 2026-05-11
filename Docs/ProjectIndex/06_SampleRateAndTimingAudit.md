# Sample Rate & Timing Audit

## Conclusion
The codebase is largely **sample-rate aware/independent**, but includes multiple **fallback defaults** to 48 kHz (and one 44.1 kHz fallback in `prepareToPlay`) when host/sample-rate inputs are invalid or unavailable.

It is **not hard-locked** to 48 kHz for normal operation.

## Evidence summary
- Processing context carries dynamic `sampleRate` per block (`ProcessContext`).
- Modules/DSP `prepare(...)` methods receive host sample rate and compute time constants from it.
- Many DSP classes initialize internal `fs` to `48000.0` as default seed/fallback.
- `PluginProcessor::prepareToPlay` fallback uses `44100.0` if invalid sampleRate arrives.

## Files with explicit 48k fallback/defaults
- `Source/Core/Processing/Modules/LevelerModule.cpp/.h`
- `Source/Core/Processing/DSP/LookaheadLimiter.cpp/.h`
- `Source/Core/Processing/DSP/BroadbandUpwardCompressor.cpp/.h`
- `Source/Core/Processing/DSP/BroadbandDownwardCompressor.cpp/.h`
- `Source/Core/Processing/DSP/SpectralUpwardCompressor.cpp/.h`
- `Source/Core/BS1770KWeighting.h`
- `Source/Core/Processing/DSP/BS1770MomentaryLufsDetector.h`
- `Source/Core/Processing/IAudioModule.h`
- `Source/Core/Processing/ProcessContext.h`

## Risk assessment
- Low risk for standard hosts (valid sample rate provided).
- Medium risk for edge/offline/tooling contexts that may create invalid or delayed sample-rate initialization.

## Recommendations
1. Standardize invalid-rate fallback policy (choose one value globally).
2. Add assertions/log counters in debug for invalid sampleRate path usage.
3. Add tests across 44.1 / 48 / 96 / 192 kHz for detector windows, attack/release calibration, and limiter lookahead timing.
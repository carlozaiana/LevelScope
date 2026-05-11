# Parameter Registry Map

Primary registration file: `Source/PluginProcessor.cpp::createParameterLayout()`.
Primary ID/range sources: `LevelerParamIDs.h`, `MultiThresholdDynamicsParamIDs.h`.

## Leveler (`lvlr.*`)
- `enabled` (bool, default 0)
- `targetLufs` (-48..-12, default -27)
- `maxBoostDb` (0..24, default 12)
- `maxCutDb` (0..24, default 12)
- `measChoice` (Auto/Momentary/Short-term; non-automatable)
- `modeChoice` (Adaptive/Learn-Hold; non-automatable)
- `learn01` (bool; non-automatable)
- `rateUpDbPerSec` (0.1..24, default 1)
- `rateDownDbPerSec` (0.1..24, default 3)
- `controlModeChoice` (Internal/Host Gain; non-automatable)
- `hostGainDb` (-24..24, automatable lane)
- `captureToHost01` (bool; non-automatable)

## MTDM high-level + zone thresholds
- `enabled` (bool)
- legacy placeholders: `thresholdDb`, `ratio`
- thresholds: `t0Lufs`, `t1Lufs`, `t2Lufs`, `t3Lufs` (all -80..0 with defaults in ParamIDs)

## Upward stage
- SUC controls: amount, maxBoost, curve, low/high knee, attack/release
- SUC advanced: FFT size choice, bands/oct choice, min/max freq, calibration trim, curve type
- Mode/structure: `upwardModeChoice`, `upEnabled01`, `upBypass01`

## Downward stage
- `downEnabled01`, `downRatio`, `downKneeDb`, `downAttackMs`, `downReleaseMs`, `downMakeupDb`, `downBypass01`

## Limiter stage
- `limEnabled01`, `limCeilingDb`, `limLookaheadMs`, `limReleaseMs`, `limAttackMs`, `limDriveDb`, `limOversamplingChoice`, `limBypass01`

## Routing / policy / audition
- MC policy: `mcPolicyChoice`, `dialogDetectorChoice`, `dialogApplyChoice`
- LFE policy: `lfeInDetector`, `lfeInApply`
- Zone solo/mutes + combinable zone audition bools (`zoneAud.*`)

## Parameter use sites
- Registration: `PluginProcessor.cpp`
- DSP bindings: `PluginProcessor.cpp` via `getRawParameterValue(...)` and module `bindParameters(...)`
- UI attachments: `PluginEditor.cpp` (SliderAttachment, ButtonAttachment, ComboAttachment)
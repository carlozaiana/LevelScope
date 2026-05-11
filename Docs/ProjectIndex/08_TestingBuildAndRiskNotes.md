# Build, Testing, and Risk Notes

## Build system snapshot
- Build entry: `CMakeLists.txt`
- Current primary artifact: JUCE-based VST3 plugin

## Existing quality anchors
- Architecture contract file: `Docs/ArchitectureSnapshot.md`
- Parameter IDs/defaults/ranges centralized for Leveler and MTDM
- Additive state-chunk compatibility strategy

## Current gaps
- No dedicated automated test suite files found in repository snapshot.
- No explicit MIDI/OSC mapping subsystem yet.
- No localization framework yet.

## Technical risk hotspots
1. **RT safety drift risk** as features grow (must preserve no-lock/no-allocation audio-thread rules).
2. **State compatibility risk** if parameter IDs/chunk schemas change without migration.
3. **Sample-rate fallback consistency** (mixed default values).
4. **Complex UX coupling** between thresholds, masks, zone audition, and module bypass logic.

## Suggested next documentation increments
- `Docs/ProjectIndex/09_SignalFlowContract.md` (block diagram + order-of-operations)
- `Docs/ProjectIndex/10_ParameterCapabilitiesMatrix.md` (automatable, realtime-safe, preset-scoped, snapshot-scoped)
- `Docs/ProjectIndex/11_StandaloneRequirements.md` (I/O, transport, batch, render determinism)
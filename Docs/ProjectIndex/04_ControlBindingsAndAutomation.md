# Control Bindings, Automation, Serialization, Presets, External Control

## Parameter registration
- Centralized in `PluginProcessor.cpp` (`createParameterLayout`).
- IDs/ranges/defaults come from `LevelerParamIDs.h` and `MultiThresholdDynamicsParamIDs.h`.

## UI bindings
- APVTS attachments in `PluginEditor.cpp` bind controls to parameter IDs.
- History/threshold interactions also route through parameter objects (threshold and zone controls).

## Automation hooks
- Automatable status explicitly set during parameter creation (`withAutomatable(false)` for selected params).
- Runtime control reads raw atomics (`getRawParameterValue`) in processor/module bindings.

## Config/state serialization
- State chunks documented in `Docs/ArchitectureSnapshot.md`:
  - `APVS` (APVTS), `MODG` (module graph), `UIST` (UI state), plus baseline history chunks.
- Backward compatibility strategy is additive, unknown chunks ignored.

## Preset systems
- Settings presets (`.lscsettings`) contain APVS+MODG.
- Snapshot/project semantics include broader state (per architecture doc).

## MIDI mapping status
- No dedicated MIDI mapping subsystem found in current source audit.
- `ProcessContext` can carry `MidiBuffer*`, but no explicit CC/Note-to-parameter mapping layer detected.

## OSC mapping status
- No OSC subsystem/files or endpoint mapping found in current source audit.

## Suggested next metadata docs (future)
- Explicit parameter capability matrix (automatable/serializable/realtime-safe/edit-lock during playback).
- Formal external-control layer spec (MIDI/OSC) before implementation.
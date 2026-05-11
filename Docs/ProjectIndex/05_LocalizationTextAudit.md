# Localization / User-Facing Text Audit

This audit flags files likely containing UI-visible text (labels, buttons, menu choices, tooltips, status text, alerts, readouts).

## High-priority text-bearing files
- `Source/PluginEditor.cpp`
- `Source/PluginEditor.h`
- `Source/VolumeHistoryComponent.cpp`
- `Source/VolumeHistoryComponent_Render.cpp`
- `Source/VolumeHistoryComponent_ViewNav.cpp`
- `Source/VolumeHistoryComponent_History.cpp`
- `Source/LoudnessStatsComponent.cpp`
- `Source/LoudnessStatsComponent.h`
- `Source/DynamicsCurveComponent.cpp`
- `Source/LevelerCurveComponent.cpp`
- `Source/PluginProcessor.cpp` (parameter display names/choices that hosts may display)

## Medium-priority text-bearing files
- `Source/Core/Processing/Modules/*ParamIDs.h` (parameter semantics/comments; can influence displayed naming)
- `Docs/*.md` (documentation text, not runtime UI)

## Recommended extraction strategy
1. Introduce string IDs (e.g., `ui.mission.target_lufs`, `ui.button.follow`).
2. Replace literals in UI components with lookup calls.
3. Add language packs (JSON/TOML/ValueTree) loaded at startup.
4. Keep parameter IDs stable; only localize display names/tooltips.
5. Add pseudo-locale test mode to detect hardcoded leftovers.

## Current status
- Localization framework not yet present.
- English text is currently embedded in source files.
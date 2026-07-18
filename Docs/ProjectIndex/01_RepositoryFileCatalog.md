# Repository File Catalog

Last reviewed: **2026-07-18**

This file maps the current repository structure to responsibilities. It is intended as a quick orientation guide before editing code or documentation.

Generated build folders, downloaded dependencies, IDE files, and packaged artifacts are intentionally excluded.

## Top-level repository files and folders

- `.circleci/config.yml` — Legacy/alternate CircleCI configuration.
- `.github/workflows/build.yml` — GitHub Actions macOS build workflow. Builds the VST3 plugin and standalone app, then uploads collected artifacts.
- `CMakeLists.txt` — Main CMake project file. Fetches JUCE, defines `LevelScopeCore`, defines the `LevelScope` VST3 plugin target, and defines the `LevelScopeStandalone` JUCE GUI app target.
- `gitignore` — Repository ignore rules.
- `Docs/` — Architecture, compliance, proposal-engine, and ProjectIndex documentation.
- `Source/` — C++ source for shared Core, plugin wrapper/UI, and standalone app shell.

## Build targets

### `LevelScopeCore`

Static library with reusable non-GUI code.

Owned by:

- `Source/Core/**`

Used by:

- `LevelScope` VST3 plugin
- `LevelScopeStandalone` standalone app

Guardrail:

- Core should remain reusable by plugin, standalone, and future headless tests.
- Core should not depend on plugin wrapper classes, APVTS, plugin editor UI, or DAW transport.

### `LevelScope`

Existing JUCE VST3 plugin target.

Owned by:

- `Source/PluginProcessor.*`
- `Source/PluginEditor.*`
- plugin-specific UI components in `Source/`

Guardrail:

- This target must remain working and unchanged unless a task explicitly requires plugin changes.
- Do not introduce standalone workflow dependencies into the plugin target.

### `LevelScopeStandalone`

Standalone JUCE GUI app target.

Owned by:

- `Source/Standalone/**`

Current status:

- Opens as a native app window.
- Provides standalone workflow scaffolding.
- Supports source-file selection metadata.
- Provides Source / Target / Current State / Export placeholder pages.
- Provides target-profile family placeholders.
- Provides current-state initialization scaffold.
- Provides workflow readiness / blocked-next-action display.

Not yet implemented:

- audio decoding
- source measurement
- current-state re-measurement
- authoritative target compliance values
- proposal engine
- render/export/report generation

## Documentation files

### Root docs

- `Docs/ArchitectureSnapshot.md` — Implemented architecture snapshot and guardrails.
- `Docs/ComplianceMatrix.md` — Human-readable compliance matrix documentation.
- `Docs/ComplianceMatrix.schema.json` — Schema for compliance matrix data.
- `Docs/ComplianceMatrix.yaml` — Machine-readable compliance matrix data.
- `Docs/ConnectionTest.md` — Connection/test note document.
- `Docs/ProposalEngineSpec_v2.md` — Standalone-first proposal-engine specification. Covers Source / Target / Current State, proposal workflow states, measurement truth, profile ownership, re-measurement, verification, and export readiness.

### ProjectIndex docs

- `Docs/ProjectIndex/README.md` — ProjectIndex navigation and recommended reading order.
- `Docs/ProjectIndex/00_ProjectOverview.md` — High-level project overview and roadmap context.
- `Docs/ProjectIndex/01_RepositoryFileCatalog.md` — This file.
- `Docs/ProjectIndex/02_DSPModuleInventory.md` — DSP and module inventory.
- `Docs/ProjectIndex/03_ParameterRegistryMap.md` — Parameter IDs, ranges, defaults, automation/state mapping.
- `Docs/ProjectIndex/04_ControlBindingsAndAutomation.md` — UI controls, APVTS bindings, and automation notes.
- `Docs/ProjectIndex/05_LocalizationTextAudit.md` — User-visible text and localization audit.
- `Docs/ProjectIndex/06_SampleRateAndTimingAudit.md` — Sample-rate and timing behavior audit.
- `Docs/ProjectIndex/07_StandardsCoverageAndRoadmap.md` — Standards coverage, compliance status, and roadmap.
- `Docs/ProjectIndex/08_TestingBuildAndRiskNotes.md` — Testing, build, and risk notes.
- `Docs/ProjectIndex/09_SignalFlowContract.md` — Signal-flow contract and routing assumptions.
- `Docs/ProjectIndex/10_ParameterCapabilitiesMatrix.md` — Parameter capability and automation matrix.
- `Docs/ProjectIndex/11_StandaloneRequirements.md` — Standalone product requirements and workflow expectations.
- `Docs/ProjectIndex/12_CodexTaskBriefTemplate.md` — Template for future Codex/web-only implementation tasks.

## Shared Core source

### `Source/Core`

- `Source/Core/PublicAPI.h` — Core public-facing umbrella/context header.
- `Source/Core/BS1770KWeighting.h` — BS.1770 K-weighting filter implementation.
- `Source/Core/LevelScopeHistoryModel.h` — Loudness/history timeline storage model.
- `Source/Core/LevelScopeHistoryModel.cpp` — Loudness/history timeline storage implementation.
- `Source/Core/RunningLoudnessStats.h` — Running integrated loudness and loudness-range statistics utilities.
- `Source/Core/RunningLoudnessStats.cpp` — Running loudness statistics implementation.

### `Source/Core/Processing`

- `Source/Core/Processing/ProcessContext.h` — Per-block processing context passed through the module chain.
- `Source/Core/Processing/IAudioModule.h` — Common module interface for prepare/process/state/bypass behavior.
- `Source/Core/Processing/ProcessorCore.h` — Core module graph host and processing-chain coordinator.

### `Source/Core/Processing/Modules`

- `Source/Core/Processing/Modules/LevelerModule.h` — Leveler module interface and state.
- `Source/Core/Processing/Modules/LevelerModule.cpp` — Leveler module implementation.
- `Source/Core/Processing/Modules/LevelerParamIDs.h` — Leveler parameter IDs, defaults, and ranges.
- `Source/Core/Processing/Modules/MultiThresholdDynamicsModule.h` — Multi-threshold dynamics module interface and state.
- `Source/Core/Processing/Modules/MultiThresholdDynamicsModule.cpp` — Multi-threshold dynamics module implementation.
- `Source/Core/Processing/Modules/MultiThresholdDynamicsParamIDs.h` — Multi-threshold dynamics parameter IDs, defaults, and ranges.
- `Source/Core/Processing/Modules/MultiThresholdDynamicsParamGroups.h` — Parameter grouping metadata for module/UI organization.

### `Source/Core/Processing/DSP`

- `Source/Core/Processing/DSP/BS1770MomentaryLufsDetector.h` — Momentary LUFS detector helper.
- `Source/Core/Processing/DSP/BroadbandDownwardCompressor.h` — Broadband downward compressor DSP interface/state.
- `Source/Core/Processing/DSP/BroadbandDownwardCompressor.cpp` — Broadband downward compressor implementation.
- `Source/Core/Processing/DSP/BroadbandUpwardCompressor.h` — Broadband upward compressor DSP interface/state.
- `Source/Core/Processing/DSP/BroadbandUpwardCompressor.cpp` — Broadband upward compressor implementation.
- `Source/Core/Processing/DSP/LookaheadLimiter.h` — Lookahead limiter DSP interface/state.
- `Source/Core/Processing/DSP/LookaheadLimiter.cpp` — Lookahead limiter implementation.
- `Source/Core/Processing/DSP/SpectralUpwardCompressor.h` — Spectral upward compressor DSP interface/state.
- `Source/Core/Processing/DSP/SpectralUpwardCompressor.cpp` — Spectral upward compressor implementation.
- `Source/Core/Processing/DSP/UpwardGainLaw.h` — Upward gain-law helper.

## Plugin wrapper and plugin UI source

These files are currently plugin-owned and should not be modified for standalone workflow work unless explicitly required.

- `Source/PluginProcessor.h` — Main JUCE `AudioProcessor` public interface, APVTS ownership, processing helpers, state/preset helper declarations.
- `Source/PluginProcessor.cpp` — Main plugin implementation: APVTS parameter layout, module binding, processing orchestration, host/latency behavior, state and preset I/O.
- `Source/PluginEditor.h` — Main plugin editor public interface and component member declarations.
- `Source/PluginEditor.cpp` — Main plugin UI implementation: controls, attachments, layout, interaction logic.
- `Source/DynamicsCurveComponent.h` — Dynamics curve visualization/control component interface.
- `Source/DynamicsCurveComponent.cpp` — Dynamics curve visualization/control implementation.
- `Source/LevelerCurveComponent.h` — Leveler curve visualization/control component interface.
- `Source/LevelerCurveComponent.cpp` — Leveler curve visualization/control implementation.
- `Source/LoudnessStatsComponent.h` — Target/current loudness stats readout component interface.
- `Source/LoudnessStatsComponent.cpp` — Target/current loudness stats readout implementation.
- `Source/VolumeHistoryComponent.h` — Volume/loudness history display component interface.
- `Source/VolumeHistoryComponent.cpp` — Volume/loudness history component core implementation.
- `Source/VolumeHistoryComponent_History.cpp` — History-specific implementation split for the volume history component.
- `Source/VolumeHistoryComponent_Render.cpp` — Rendering-specific implementation split for the volume history component.
- `Source/VolumeHistoryComponent_ViewNav.cpp` — View/navigation-specific implementation split for the volume history component.

## Standalone app source

These files are standalone-owned. They should not include or depend on `PluginProcessor`, `PluginEditor`, APVTS, DAW transport, or plugin state chunks.

- `Source/Standalone/StandaloneApp.cpp` — JUCE standalone application entry point and main window setup.
- `Source/Standalone/StandaloneMainComponent.h` — Standalone main component interface and UI member declarations.
- `Source/Standalone/StandaloneMainComponent.cpp` — Standalone workflow-shell UI implementation.
- `Source/Standalone/StandaloneSessionModel.h` — Standalone session state model for Source, Target, Current State, and selected workflow page.
- `Source/Standalone/StandaloneTargetProfiles.h` — Standalone target-profile family scaffold and placeholder profile metadata.
- `Source/Standalone/StandaloneWorkflowReadiness.h` — Standalone readiness checklist and blocked-next-action logic.

## Ownership summary

### Shared / reusable

- `Source/Core/**`

These files should remain usable by plugin, standalone, and future tests.

### Plugin-only for now

- `Source/PluginProcessor.*`
- `Source/PluginEditor.*`
- plugin UI components directly under `Source/`

These files may eventually share extracted UI/model helpers, but they should remain plugin-owned until a deliberate refactor is planned.

### Standalone-only for now

- `Source/Standalone/**`

These files hold the current standalone workflow scaffold. They may later call shared Core measurement/proposal/export services, but the scaffold itself should not be moved into Core prematurely.

## Current standalone development sequence

Recommended conservative order:

1. Standalone workflow shell.
2. Source import + measurement view scaffold.
3. Target profile system.
4. Current-state re-measure flow scaffold.
5. Workflow readiness gates.
6. Source measurement implementation.
7. Authoritative target profile loading.
8. Current-state re-measurement.
9. Proposal-engine v1.
10. Render/export/report generation.

## Guardrails

- Keep the VST3 plugin working.
- Keep the standalone app independent from plugin wrapper classes.
- Keep Core reusable and non-GUI.
- Do not add proposal-engine logic inside UI files.
- Do not claim compliance support until authoritative profile values and measurement verification are implemented.
- Treat placeholder target profiles and picker file patterns as scaffolding, not final product promises.
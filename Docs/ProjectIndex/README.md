# ProjectIndex

This folder is the persistent orientation layer for future development chats. Read this before touching code.

## Recommended reading order
1. `00_ProjectOverview.md`
2. `01_RepositoryFileCatalog.md`
3. `02_DSPModuleInventory.md`
4. `03_ParameterRegistryMap.md`
5. `04_ControlBindingsAndAutomation.md`
6. `05_LocalizationTextAudit.md`
7. `06_SampleRateAndTimingAudit.md`
8. `07_StandardsCoverageAndRoadmap.md`
9. `08_TestingBuildAndRiskNotes.md`

## Task-based quick start
- **Implement DSP/module behavior:** 02 + 03 + 06
- **UI overhaul/localization:** 01 + 04 + 05
- **Standalone implementation / next logic-module chat:** [13_StandaloneScaffoldHandover.md](13_StandaloneScaffoldHandover.md), then 00 + 09 + 11 and the proposal spec
- **Video implementation:** [15_VideoImplementationHandover.md](15_VideoImplementationHandover.md), after the measurement qualification gate
- **VST3 measurement parity:** [16_Vst3MeasurementParityMilestone.md](16_Vst3MeasurementParityMilestone.md), a separate future milestone
- **Automation/preset/state compatibility:** 03 + 04 + 08

## Additional design specs

* [14_OfflineMeasurementContract.md](14_OfflineMeasurementContract.md) — first standalone measurement milestone, numerical boundaries, source/current/render identity and video timing design; read with the continuation evidence in the handover.

* [`../ProposalEngineSpec_v2.md`](../ProposalEngineSpec_v2.md) — detailed standalone-first proposal-engine design, including workflow states, shared engine modules, separate proposal methods, Source/Target/Current State model, freeze-aware re-measurement, and verification flow.

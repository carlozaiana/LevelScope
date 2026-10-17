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
- **Working-copy locations, GitHub video checkpoint and preserved artifacts:** [19_DevelopmentWorkflow.md](19_DevelopmentWorkflow.md)
- **Implement DSP/module behavior:** 02 + 03 + 06
- **UI overhaul/localization:** 01 + 04 + 05
- **Standalone implementation / next logic-module chat:** [13_StandaloneScaffoldHandover.md](13_StandaloneScaffoldHandover.md), then 00 + 09 + 11 and the proposal spec
- **Video implementation:** [15_VideoImplementationHandover.md](15_VideoImplementationHandover.md), continuing from `video-checkpoint-2026-10-17` in the video worktree
- **Implemented reference picture and qualification:** [18_VideoQualificationResults.md](18_VideoQualificationResults.md), including actual two-platform tests and physical-sync limits
- **M1 qualification evidence:** [17_M1QualificationResults.md](17_M1QualificationResults.md)
- **VST3 measurement parity:** [16_Vst3MeasurementParityMilestone.md](16_Vst3MeasurementParityMilestone.md), a separate future milestone
- **Automation/preset/state compatibility:** 03 + 04 + 08

## Additional design specs

* [14_OfflineMeasurementContract.md](14_OfflineMeasurementContract.md) — first standalone measurement milestone, numerical boundaries, source/current/render identity and video timing design; read with the continuation evidence in the handover.

* [`../ProposalEngineSpec_v2.md`](../ProposalEngineSpec_v2.md) — detailed standalone-first proposal-engine design, including workflow states, shared engine modules, separate proposal methods, Source/Target/Current State model, freeze-aware re-measurement, and verification flow.

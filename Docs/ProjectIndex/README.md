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

- **Video paint diagnostics / integrated tooling; intermittent skips unresolved; M7 frozen:** [Diagnostic integration handover](../VideoPaintDiagnosticMainIntegration-2026-10-24.md)
- **M6 verified Current State audition / qualified, published and integrated into main:** [28_M6VerifiedCurrentStateAudition.md](28_M6VerifiedCurrentStateAudition.md)
- **M6 candidate identity, completed qualification and publication handover:** [30_M6QualificationAndHandover.md](30_M6QualificationAndHandover.md)
- **Future channel-based multichannel support / deferred:** [29_MultichannelPcmMilestone.md](29_MultichannelPcmMilestone.md)
- **Bounded SDR ProRes / qualified and integrated into main:** [27_ProResMainIntegration.md](27_ProResMainIntegration.md) and [26_ProResQualification.md](26_ProResQualification.md), with [frozen acceptance criteria](25_ProResAcceptance.md)
- **M5 programme sections / qualified, published and integrated into main:** [24_M5ProgrammeSections.md](24_M5ProgrammeSections.md)
- **M4 measured diagnostics / published checkpoint:** [23_M4MeasuredDiagnostics.md](23_M4MeasuredDiagnostics.md)
- **M3 global trim checkpoint:** [22_M3GlobalTrimProposals.md](22_M3GlobalTrimProposals.md)
- **M2 profiles and measured rules:** [21_M2MeasuredCompliance.md](21_M2MeasuredCompliance.md)
- **Accepted milestones on main / subsequent logic-module work:** [20_MainIntegrationAcceptance.md](20_MainIntegrationAcceptance.md)
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

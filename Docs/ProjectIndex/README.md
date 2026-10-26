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

- **M7 channel-based PCM 7.1 / qualified, published and integrated into main:** [checkpoint](https://github.com/carlozaiana/LevelScope/tree/m7-pcm-7.1-2026-10-25); [proposal](31_M7ChannelBasedPcm71Proposal.md); [qualification record](32_M7QualificationAndHandover.md); [completion decision](34_M7CompletionAcceptance.md)
- **Optional video paint diagnostics / integrated investigation tooling:** [diagnostic integration handover](../VideoPaintDiagnosticMainIntegration-2026-10-24.md); [historical M7 video handover](33_M7VideoPaintTimingHandover.md)
- **M6 verified Current State audition / qualified, published and integrated into main:** [28_M6VerifiedCurrentStateAudition.md](28_M6VerifiedCurrentStateAudition.md)
- **M6 candidate identity, completed qualification and publication handover:** [30_M6QualificationAndHandover.md](30_M6QualificationAndHandover.md)
- **M8 sessions / qualified, published and integrated into main:** [checkpoint](https://github.com/carlozaiana/LevelScope/tree/m8-reproducible-sessions-2026-10-25); [36_M8SessionsQualificationAndHandover.md](36_M8SessionsQualificationAndHandover.md)
- **M9 separate-file stereo/5.1 import — candidate-04 qualified, published and integrated into main:** [checkpoint](https://github.com/carlozaiana/LevelScope/tree/m9-separate-file-stereo-5.1-2026-10-25); [37_M9SeparateFileStereo51Import.md](37_M9SeparateFileStereo51Import.md); [38_M9QualificationAndHandover.md](38_M9QualificationAndHandover.md)
- **M10 verified PCM export and reporting — candidate-03 qualified and published:** [39_M10VerifiedPcmExportAndReporting.md](39_M10VerifiedPcmExportAndReporting.md)
- [40 - M10 qualification and handover](40_M10QualificationAndHandover.md) - actual Sequoia/High Sierra evidence, tested identities, published checkpoint/main mapping and limitations.
- **Approved M8 contract and following stereo/5.1 roadmap:** [35_PostM7BroadcastWorkflowProposal.md](35_PostM7BroadcastWorkflowProposal.md)
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

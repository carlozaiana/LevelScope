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

- **M18 proposed — awaiting approval, no implementation:** [63 — interactive audio timeline and progressive analysis](63_M18InteractiveTimelineAndProgressiveAnalysis.md). Reuse existing processing; waveform/M/S previews, exact navigation, section-gain drafts and restored historical proposal explanations. Combined actual High Sierra/video/three-hour qualification planned after M18, no later than M19.

- **M17 candidate-04 published development checkpoint — Sequoia tested; High Sierra, video and three-hour checks pending:** [62 — qualification and handover](62_M17QualificationAndHandover.md). Optional explicit job LRA ceiling, independent named presets, reused measured evaluation and strict all-selected-rule export gates. [60 — approved scope](60_M17JobLraRequirementsAndPresets.md); [61 — impact plan](61_M17ImpactQualificationPlan.md). [Checkpoint](https://github.com/carlozaiana/LevelScope/tree/m17-job-lra-presets-sequoia-2026-10-30), commit `e422f85ccbf7ae11876e6a31be2409937cb8381c`; [publication receipt](../../Tests/Reference/results/m17-publication-2026-10-30/README.md). Main stays M13. Combined qualification follows M18, no later than M19.

- **M16 candidate-03 published development checkpoint — Sequoia tested; High Sierra, video and three-hour checks pending:** [checkpoint](https://github.com/carlozaiana/LevelScope/tree/m16-proposal-preferences-sequoia-2026-10-30); [56 — bounded proposal preferences and configurable limits](56_M16ProposalPreferencesAndLimits.md). Explicit additional/cumulative gain, reserve, changed-region and duration budgets; density-aware accounting, policy identity and session/report provenance. [59 — exact candidate identity and handover](59_M16QualificationAndHandover.md). Main remains M13; High Sierra/video/three-hour qualification deferred.

- **M15 candidate-05 published development checkpoint — Sequoia tested; High Sierra, video and three-hour checks pending:** [checkpoint](https://github.com/carlozaiana/LevelScope/tree/m15-marker-import-sequoia-2026-10-29); [53 — marker import and explicit audio-source time mapping](53_M15MarkerImportAndSourceTimeMapping.md) — explicit Nuendo CSV / Pro Tools Session Info as Text point-marker readers and paste, exact preview/confirmation, annotation-only identity and session provenance. M16 is now published as a development checkpoint; see documents 56–59. Main remains M13 and High Sierra/video/three-hour qualification remains pending. See [54 — impact plan](54_M15ImpactQualificationPlan.md) and [55 — qualification/handover](55_M15QualificationAndHandover.md).
- **M14 candidate-02 published development checkpoint — Sequoia tested; High Sierra, video and three-hour checks pending:** [checkpoint](https://github.com/carlozaiana/LevelScope/tree/m14-mono-rf64-sequoia-2026-10-29); [50 — complete mono processing and verified PCM export](50_M14MonoProcessingAndVerifiedExport.md) — existing manual/bounded section workflow extended to mono, explicit role-preserving PCM24 output, odd-byte padding, sessions and impact-selected Sequoia qualification. Bounded direct cinema-layout RF64 import is included; High Sierra/video/three-hour qualification deferred. [51 — impact plan](51_M14ImpactQualificationPlan.md); [52 — handover/status](52_M14QualificationAndHandover.md).
- **M7 channel-based PCM 7.1 / qualified, published and integrated into main:** [checkpoint](https://github.com/carlozaiana/LevelScope/tree/m7-pcm-7.1-2026-10-25); [proposal](31_M7ChannelBasedPcm71Proposal.md); [qualification record](32_M7QualificationAndHandover.md); [completion decision](34_M7CompletionAcceptance.md)
- **Optional video paint diagnostics / integrated investigation tooling:** [diagnostic integration handover](../VideoPaintDiagnosticMainIntegration-2026-10-24.md); [historical M7 video handover](33_M7VideoPaintTimingHandover.md)
- **M6 verified Current State audition / qualified, published and integrated into main:** [28_M6VerifiedCurrentStateAudition.md](28_M6VerifiedCurrentStateAudition.md)
- **M6 candidate identity, completed qualification and publication handover:** [30_M6QualificationAndHandover.md](30_M6QualificationAndHandover.md)
- **M8 sessions / qualified, published and integrated into main:** [checkpoint](https://github.com/carlozaiana/LevelScope/tree/m8-reproducible-sessions-2026-10-25); [36_M8SessionsQualificationAndHandover.md](36_M8SessionsQualificationAndHandover.md)
- **M9 separate-file stereo/5.1 import — candidate-04 qualified, published and integrated into main:** [checkpoint](https://github.com/carlozaiana/LevelScope/tree/m9-separate-file-stereo-5.1-2026-10-25); [37_M9SeparateFileStereo51Import.md](37_M9SeparateFileStereo51Import.md); [38_M9QualificationAndHandover.md](38_M9QualificationAndHandover.md)
- **M10 verified PCM export and reporting — candidate-03 qualified and published:** [39_M10VerifiedPcmExportAndReporting.md](39_M10VerifiedPcmExportAndReporting.md)
- [40 - M10 qualification and handover](40_M10QualificationAndHandover.md) - actual Sequoia/High Sierra evidence, tested identities, published checkpoint/main mapping and limitations.
- **M11 candidate-02 qualified, published and integrated into main:** [checkpoint](https://github.com/carlozaiana/LevelScope/tree/m11-aligned-stem-composition-2026-10-26); [41_M11AlignedStemComposition.md](41_M11AlignedStemComposition.md) — aligned M&E/dialogue composition and future processing boundaries; [42_M11ImpactQualificationPlan.md](42_M11ImpactQualificationPlan.md) — pre-implementation impact-selected fresh/reused evidence plan, budgets and completion gates.
- [43 — M11 implementation, qualification and handover](43_M11QualificationAndHandover.md) — implemented composition, completed actual-platform qualification, published checkpoint/main mapping and remaining product limits.
- **M12 candidate-03: Sequoia tested; High Sierra and video qualification pending.** [Development checkpoint](https://github.com/carlozaiana/LevelScope/tree/m12-section-trim-sequoia-2026-10-27); [44 — explicit section-trim processing and measured verification](44_M12SectionTrimProcessing.md) — manual whole-programme section gains, shared processing, M3 continuation, sessions/audition/export. Originally published on `codex/logic-module` with main held at M11; now included in the combined M13 checkpoint. Automatic section-gain selection is subsequently implemented in M13; [49](49_M13QualificationAndHandover.md) records the separate combined qualification.
- [45 — M12 change-impact qualification plan](45_M12ImpactQualificationPlan.md) — required Sequoia audio checks, bounded resources, historical reuse and explicit deferrals.
- [46 — M12 implementation and qualification handover](46_M12QualificationAndHandover.md) — current evidence and phase status; no main integration before combined post-M13 qualification.
- **M13 candidate-02 qualified and published:** [checkpoint](https://github.com/carlozaiana/LevelScope/tree/m13-bounded-section-proposals-2026-10-28); [47 — first bounded section-led method and measured verification](47_M13BoundedSectionTrimProposals.md). Three-strategy architecture, shared configurable-preference design and optional-marker roadmap; M13 implements fixed conservative defaults and no limiter/dynamics. The single cinema evaluation completed with an honest policy-limited failure; completed actual-platform evidence and limitations are in [49](49_M13QualificationAndHandover.md). Main receives the combined M12/M13 checkpoint by fast-forward.
- [48 — M13 change-impact qualification plan](48_M13ImpactQualificationPlan.md) — explicit affected tests, bounded workloads, resource budgets and combined actual-platform qualification gates.
- [49 — M13 implementation and qualification handover](49_M13QualificationAndHandover.md) — completed qualification, exact candidate identity, fresh/reused evidence, cinema-copy provenance and remaining limitations.
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

* [`../ProposalEngineSpec_v2.md`](../ProposalEngineSpec_v2.md) — revised three-strategy design: dynamics-led, section-led and mixed; separate selection algorithms and processing preferences; shared truth/identity, candidate/Current/artifact verification and export gates. “Combined-analysis” replaces the ambiguous design term “hybrid”; bounded M13 implementation/qualification is subsequently approved; subsequent publication and main integration are authorized and recorded by the M13 checkpoint.

- [57 — M16 impact qualification plan](57_M16ImpactQualificationPlan.md)
- [58 — LRA workflow and interactive timeline design](58_LraWorkflowAndInteractiveTimelineDesign.md)

- [59 — M16 qualification and handover](59_M16QualificationAndHandover.md)

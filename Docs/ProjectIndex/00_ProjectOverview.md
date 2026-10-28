# Project Overview

M13 qualification update, 28 October 2026: bounded candidate-02 is implemented
and qualified on actual Sequoia/High Sierra; [ProjectIndex 47](47_M13BoundedSectionTrimProposals.md) is authoritative for this scope. The broader
no-change rule covers all applicable explicitly selected, method-defined objectives;
M13 still evaluates only supported programme rules. The intended cinema comparison
is user-reported Nuendo statistics supplied 28 October, while the historical screenshot
is Cubase 10 supplementary evidence. No exact Nuendo version is claimed. The user authorized checkpoint publication
and combined M12/M13 main integration on 28 October. The preserved M12 qualification
status remains unchanged; [49](49_M13QualificationAndHandover.md) records the completed combined
M12/M13 qualification. The [published M13 checkpoint](https://github.com/carlozaiana/LevelScope/tree/m13-bounded-section-proposals-2026-10-28)
identifies the combined main revision at publication.

## Product goal
LevelScope is evolving from a VST3 plugin into a standalone loudness-conforming app for film/TV post-production. It aims to transform cinema/TV mixes to broadcaster specs with minimal sonic damage.

## Measurement and profile targets
- The current versioned EBU profile evaluates programme loudness and true peak.
- LRA is informative for this profile; there is no universal LRA ceiling.
- Future dialogue-specific methods and intelligibility objectives require separate
  definitions and validation; they are not current compliance rules.

Standards priority: EBU R128 and ATSC A/85.

## Current implemented state
- Working JUCE VST3 plugin
- Loudness history display (Momentary, Short-term, gate, running LRA)
- 4 draggable threshold lines (T0..T3)
- Processing modules present: Leveler, Upward compressor (spectral + broadband mode), Downward compressor, Lookahead limiter (optional oversampling)
- Preset/snapshot infrastructure and APVTS-based parameter system
- Standalone scaffold: validated background media import, confirmed file-channel roles, explicit output mapping/device selection and original-source audition; cancellable jobs and future service interfaces. See [scaffold handover](13_StandaloneScaffoldHandover.md) for actual validation and limitations.
- Qualified standalone M1 measurement, synchronized reference picture, M2 measured programme rules, [M3 global trim](22_M3GlobalTrimProposals.md), [M4 diagnostics](23_M4MeasuredDiagnostics.md), [M5 in-memory sections](24_M5ProgrammeSections.md) and [bounded SDR ProRes](27_ProResMainIntegration.md) are published on main.
- [M6 verified Current State audition](28_M6VerifiedCurrentStateAudition.md) candidate 01 is qualified on actual Intel High Sierra and arm64 Sequoia, published as `m6-verified-current-state-audition-2026-10-24` and integrated into main. A shared static-gain sample operation serves offline measurement and the bounded playback worker. Source/Current selection, ordered roles, processing gain and routing are stopped-only; monitor volume remains live. Playback never establishes measurement or compliance. [The handover](30_M6QualificationAndHandover.md) records exact source/binary identity and completed evidence. Section processing and public export were outside M6; VST3 parity remains deferred.
- [M8 reproducible sessions](36_M8SessionsQualificationAndHandover.md) is qualified,
  published and integrated into main. Documents preserve settings, section
  annotations, media identities and rational video alignment; replaceable caches
  are accepted conservatively. Reopening is stopped, device-closed and routes-muted.
- [M9 separate-file stereo/5.1 import](38_M9QualificationAndHandover.md) is implemented
  and qualified in frozen candidate-04 on actual Sequoia and High Sierra, published
  as `m9-separate-file-stereo-5.1-2026-10-25` and integrated into main. A versioned programme assembly binds explicit ordered
  roles and sample-zero alignment confirmation to every mono member's identity.
  An all-or-nothing verified group and bounded worker reader feed measurement and
  audition. Session schema 2 supports grouped recovery and explicit per-file relink;
  existing interleaved stereo/5.1/7.1 behavior remains covered by regressions.
- [M10 verified PCM export/reporting](40_M10QualificationAndHandover.md) is qualified
  and published. Current State and the finalized reopened PCM24 artifact must both
  pass the selected rules before explicit local publication. Final quantization
  applies deterministic TPDF once; measurement-only reports remain separate.
- [M11 aligned M&E/dialogue composition](43_M11QualificationAndHandover.md) candidate-02
  is qualified on actual Sequoia/High Sierra, published as
  `m11-aligned-stem-composition-2026-10-26` and integrated into main. Original component
  identities and exact mappings feed one shared summing path before global trim;
  schema 3 retains older-session migration. Future component DSP, whole-programme
  processing plans and non-audible analysis references have distinct identity
  boundaries. Unity component gain is only the M11 scope boundary.
- [M12 explicit section processing](46_M12QualificationAndHandover.md):
  **Sequoia tested; High Sierra and video qualification pending.** Candidate-03 is
  published as a development checkpoint on `codex/logic-module` while main was held
  at M11; the combined M13 checkpoint now includes this work. Typed manual section gains follow composition and precede
  global trim through the same measurement/audition/export processor. Drafts,
  applied settings and measured results remain distinct; schema 4 and report v3
  preserve the exact plan. The combined M13 candidate has now completed actual
  High Sierra and selected video qualification; its separate evidence is in [49](49_M13QualificationAndHandover.md).
- [M13 bounded section-trim proposal](47_M13BoundedSectionTrimProposals.md) is implemented
  and qualified as candidate-02. The 28 October design revision makes it the first
  bounded section-led method: global-first planning, bounded local peak correction,
  private candidate assessment, explicit Apply and measured Current verification.
  It has no limiter/dynamics; the 0.5 dB reserve is inside the 3 dB cap.
  Combined Sequoia/High Sierra and selected video qualification passed; the user authorized checkpoint
  publication and main integration on 28 October. The [M13 checkpoint](https://github.com/carlozaiana/LevelScope/tree/m13-bounded-section-proposals-2026-10-28)
  identifies the unchanged tested implementation now integrated by fast-forward.

## Medium-term direction
- Approved M7 implementation extends the standalone measurement/Current path to
  explicitly confirmed BS.2051 System I PCM 7.1. Combined candidate-03 is qualified
  on both actual Macs as of 25 October, including main's optional diagnostic tools;
  candidate-02 and earlier evidence remain preserved. The annotated checkpoint
  [`m7-pcm-7.1-2026-10-25`](https://github.com/carlozaiana/LevelScope/tree/m7-pcm-7.1-2026-10-25)
  is published and integrated into main.
  [The revised acceptance](34_M7CompletionAcceptance.md) permits occasional picture
  skips while retaining audio, synchronization, identity and resource gates.
  [The qualification record](32_M7QualificationAndHandover.md) distinguishes offline,
  injected routing, video/UI and separately pending full physical 7.1 listening.
- Defer 7.1.4/immersive expansion and prioritize the practical stereo/5.1 workflow:
  qualified M8 sessions, M9 channel assembly and M10 verified PCM export/reporting.
  M11 adds explicit confirmed composition. Later whole-programme processing and
  dialogue levelling relative to M&E need bounded proposals, explicit application
  and verification of both the intended improvement and the complete programme.
  Analysis-only references must remain non-audible; channel assembly, stem
  composition and dialogue-assisted analysis are separate capabilities.
- [Shared proposal-engine design](../ProposalEngineSpec_v2.md): dynamics-led,
  section-led and mixed processing strategies, separate from algorithms and user
  preferences. All share profiles, measurements, processing, identities and strict
  verification/export gates. Available processors need not be active.
- Following M13: separately review mono section/export support and optional basic
  marker import; qualify standalone limiter, downward and upward DSP in small steps;
  then dynamics-led proposals, enhanced section-led proposals and the mixed strategy.
  Avoid assigning future milestone numbers prematurely. Generate the selected
  strategy only; comparisons are explicit, bounded and use the same baseline.
- Markers are optional. Editorial edits, processing sections, loudness regions and
  peak events retain distinct meanings/provenance. Component-relative dialogue/M&E
  processing and validated intelligibility objectives follow the broader workflow.

## Platform direction
- Current format: VST3
- Standalone development now runs on an Apple Silicon Mac Studio with macOS Sequoia; earlier development was on High Sierra. The configured 10.13 deployment minimum is distinct from OS versions actually tested. Later: Windows qualification.

# Project Overview

## Product goal
LevelScope is evolving from a VST3 plugin into a standalone loudness-conforming app for film/TV post-production. It aims to transform cinema/TV mixes to broadcaster specs with minimal sonic damage.

## Primary compliance targets
- Integrated Loudness (LUFS)
- Max Peak / True Peak
- Loudness Range (LRA)
- Future: Dialog Intelligibility

Standards priority: EBU R128 and ATSC A/85.

## Current implemented state
- Working JUCE VST3 plugin
- Loudness history display (Momentary, Short-term, gate, running LRA)
- 4 draggable threshold lines (T0..T3)
- Processing modules present: Leveler, Upward compressor (spectral + broadband mode), Downward compressor, Lookahead limiter (optional oversampling)
- Preset/snapshot infrastructure and APVTS-based parameter system
- Standalone scaffold: validated background media import, confirmed file-channel roles, explicit output mapping/device selection and original-source audition; cancellable jobs and future service interfaces. See [scaffold handover](13_StandaloneScaffoldHandover.md) for actual validation and limitations.
- Qualified standalone M1 measurement, synchronized reference picture, M2 measured programme rules, [M3 global trim](22_M3GlobalTrimProposals.md), [M4 diagnostics](23_M4MeasuredDiagnostics.md), [M5 in-memory sections](24_M5ProgrammeSections.md) and [bounded SDR ProRes](27_ProResMainIntegration.md) are published on main.
- [M6 verified Current State audition](28_M6VerifiedCurrentStateAudition.md) candidate 01 is qualified on actual Intel High Sierra and arm64 Sequoia, published as `m6-verified-current-state-audition-2026-10-24` and integrated into main. A shared static-gain sample operation serves offline measurement and the bounded playback worker. Source/Current selection, ordered roles, processing gain and routing are stopped-only; monitor volume remains live. Playback never establishes measurement or compliance. [The handover](30_M6QualificationAndHandover.md) records exact source/binary identity and completed evidence. Section processing, public export and VST3 parity remain deferred.
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
  qualified M8 session saving/reopening and M9 separate-file import, then a
  bounded proposal for verified PCM export and artifact-bound reports. Stem
  composition and dialogue-assisted analysis follow as separate proposals;
  channel assembly does not implement either one.
- Standalone app with analyze/propose/apply/export workflow
- Intelligent logic/policy engine for minimal-intervention target matching
- Later: stem workflow + marker-track/scene-aware processing

## Platform direction
- Current format: VST3
- Standalone development now runs on an Apple Silicon Mac Studio with macOS Sequoia; earlier development was on High Sierra. The configured 10.13 deployment minimum is distinct from OS versions actually tested. Later: Windows qualification.

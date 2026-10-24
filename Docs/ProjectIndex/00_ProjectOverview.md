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
- [M6 verified Current State audition](28_M6VerifiedCurrentStateAudition.md) candidate 01 is qualified on actual Intel High Sierra and arm64 Sequoia, published as `m6-verified-current-state-audition-2026-10-24` and integrated into main. A shared static-gain sample operation serves offline measurement and the bounded playback worker. Source/Current selection, ordered roles, processing gain and routing are stopped-only; monitor volume remains live. Playback never establishes measurement or compliance. [The handover](30_M6QualificationAndHandover.md) records exact source/binary identity and completed evidence. Section processing, persistence, public export and VST3 parity remain deferred.

## Medium-term direction
- Standalone app with analyze/propose/apply/export workflow
- Intelligent logic/policy engine for minimal-intervention target matching
- Later: stem workflow + marker-track/scene-aware processing

## Platform direction
- Current format: VST3
- Standalone development now runs on an Apple Silicon Mac Studio with macOS Sequoia; earlier development was on High Sierra. The configured 10.13 deployment minimum is distinct from OS versions actually tested. Later: Windows qualification.

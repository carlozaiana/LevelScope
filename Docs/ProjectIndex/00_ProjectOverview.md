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
- Qualified standalone M1 measurement, synchronized reference picture, M2 measured programme rules and [M3 global trim](22_M3GlobalTrimProposals.md) are published on main. M3 adds deterministic proposals, explicit application and actual measured verification; its handover records the checkpoint and qualification evidence. Processed audition, public export and VST3 parity remain deferred.

## Medium-term direction
- Standalone app with analyze/propose/apply/export workflow
- Intelligent logic/policy engine for minimal-intervention target matching
- Later: stem workflow + marker-track/scene-aware processing

## Platform direction
- Current format: VST3
- Standalone development now runs on an Apple Silicon Mac Studio with macOS Sequoia; earlier development was on High Sierra. The configured 10.13 deployment minimum is distinct from OS versions actually tested. Later: Windows qualification.

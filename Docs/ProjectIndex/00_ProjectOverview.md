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

## Medium-term direction
- Standalone app with analyze/propose/apply/export workflow
- Intelligent logic/policy engine for minimal-intervention target matching
- Later: stem workflow + marker-track/scene-aware processing

## Platform direction
- Current format: VST3
- Planned: macOS standalone first (including older macOS support), later Windows
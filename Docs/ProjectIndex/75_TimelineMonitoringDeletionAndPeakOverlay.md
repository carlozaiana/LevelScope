# Timeline monitoring, deletion and peak overlay — 7 November 2026

Continue on `codex/m20-ui-integration`; publication/main integration remains pending.
See the [complete handover](../TimelineSafetyAndEditing-2026-11-07.md) for preserved
build, exact commit, qualification receipts and manual checks.

- Persistent main-window monitor-overload warning names speaker/file channel/output
  and source position, with recovery advice and a Monitor volume button. It never
  starts playback, changes processing or overrides connection intent.
- Click a main gain-lane item then Backspace. Annotation deletion merges its span
  into a neighbour, with one-level Undo delete and no audio/measurement change.
  Applied-region deletion stages 0 dB local gain on the applied interval, preserving
  global/other plateau values. Explicit Apply/Discard and independent remeasurement
  remain required. Candidate proposal regions remain immutable previews.
- Optional, initially hidden over-target peak layer: retained PCM sample extrema,
  numeric amplitude ceiling lines and measured per-channel programme TP maxima.
  Not an exhaustive inter-sample event map. No Limit, previews and unverified Current
  cannot acquire peak warnings. Existing overview budget remains 128 MiB.

No proposal-search, measurement/DSP, export or session-schema change. An applied
manual removal keeps the last proposal explanation intact and historical. Testing
uses injected outputs only, on actual Sequoia and High Sierra, with unchanged
qualification evidence/deferrals reused.

Final product/test commit: `3c853e4c855b047ed990a899fde586fa2d952aaf`.
Preserved package: `build/ui-safety-editing-2026-11-07/candidate-02/LevelScope Standalone.app`.
Actual Sequoia 15.7.9 arm64 and High Sierra 10.13.6 Intel each pass 1,438 focused
checks (1,199 timeline + 128 native UI/history + 111 injected monitor). The final
removal draft field displays 0 dB and never claims neighbouring gains are reset.
The handover links the full evidence and manual checklist. Existing running app,
session files, prior packages and main remain unchanged.

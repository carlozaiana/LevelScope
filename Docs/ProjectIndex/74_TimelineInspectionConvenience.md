# Timeline inspection convenience — 7 November 2026

Implemented locally on `codex/m20-ui-integration`, before the user's M19/M20 listening
review. Main publication/integration remains pending.

Read [the complete handover](../UIInspectionConvenience-2026-11-07.md) for exact app,
session paths, limits, source commits, test evidence and manual checklist.

The slice adds the left-margin waveform/loudness blend, downward divider travel
with scroll-accessible lower controls, aligned gain readouts, selection of applied
processing regions after reopen, qualified per-channel sample-peak inspection,
stopped Cmd+Shift+A for paired inspection/audition switching, programme filename,
Target/Original/Current row order, fixed header preset controls, pointer-anchored
wheel zoom and the existing plugin background palette. Short-term colours describe
measured LRA membership. The IL relative gate remains in both domains; only Current
adds red double-width parts for actual measured momentary differences.

No measurement, processing, proposal, export, device-policy or session-schema changes.
The combined waveform overview budget remains 128 MiB. Detailed sample reading and
new regional true-peak measurement remain separate work.

The app implementation is `a37938e33ba764ebffff901b7b4cf3864a5b95cc`; final preserved
candidate 02 is `086f316d8f7e1a6387697487f140ececf3dd60fc`, with byte-identical product
binary and only test-utility recovery changes. Each actual Mac passes 1,505 focused
checks. Music recovery passes 65 additional reproducibility checks.

The previous M20 music-session pathname had been saved with later DialogAmbiance
manual edits. It is preserved. Separate new music copies reproduce the archived
plans/results exactly and reopen independently verified; regenerated explanation
provenance is identified as new. Listening quality remains for the user to assess.

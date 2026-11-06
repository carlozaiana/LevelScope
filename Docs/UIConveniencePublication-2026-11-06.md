# UI convenience publication — 6 November 2026

The user approved publication and main integration of all three completed UI slices.
The continuation branch and main were published together with the annotated checkpoint,
using a normal fast-forward from M19. No product change, rebuild or repeated runtime
qualification was needed. This subsequent receipt changes documentation only.

| Identity | Value |
| --- | --- |
| Published checkpoint | [`ui-convenience-2026-11-06`](https://github.com/carlozaiana/LevelScope/tree/ui-convenience-2026-11-06) |
| Checkpoint / initially integrated main | `b7350128d7b65507a5db434e38ecd73b8f2c4cfd` |
| Annotated tag object | `8edcd9000d62556b7b253658f6ca7c1f399f1469` |
| Previous main | `5fc8c6ee95e705077cf996cecbc6543dbc20d257` — M19 and Intel optimisation |
| Published continuation | `codex/ui-device-preparation-m19` |
| Combined tested product | `96159ebd79f762bea147b41ad696fd80f449dfc9` |
| Final main revision | This subsequent documentation-only receipt, fast-forwarded onto main and the continuation. The final response and local `build/ui-convenience-publication-2026-11-06/final-verification.json` record its exact hash; no commit embeds its own hash. |

## Short handover for the other chats

1. **Output connection:** remember the last successfully opened interface per Mac,
   using stable identity where available. Restore the visible selection. Explicit
   selection opens once the programme rate is known; automatic reconnect is opt-in,
   initially off. Disconnect suppresses automatic attempts for the rest of the app
   session. Missing, busy, ambiguous or incompatible devices are explained without
   substitution. Offline analysis remains available; connection never starts playback.
2. **Initial import and analysis:** a main-window banner shows validation, preparation /
   verification, measurement and finalisation, with meaningful per-stage progress and
   cancellation. Progressive waveform/loudness previews retain incomplete labels and
   completed-result authority. Cancellation does not auto-retry or install partial results.
3. **Monitoring memory:** speaker-role routes are local to each Mac/interface and stereo
   or 5.1 layout. Restore only after role confirmation and an exact output-configuration
   match. Reordered files keep intended destinations; deliberate mutes, including
   all-muted configurations, survive. Changed outputs remain muted for review.
   Explicit user monitor volume is remembered per interface; internal resets do not
   overwrite it. Existing range, smoothing and overload safeguards remain. These
   preferences are separate from portable sessions, processing and measurements.

M19's proposal, measurement, history and export implementation is unchanged. Strategy /
Audio Quality remains design-only. First-playback preparation feedback, session-wide
verified-source reuse, optional project-owned audio, earlier waveform generation,
the inherited completion pause and 5.1 stereo monitoring remain separate work.

## Preserved builds and qualification

Use the combined tested app without rebuilding over it:

`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/ui-monitoring-memory-2026-11-05/candidate-01/LevelScope Standalone.app`

Its identical High Sierra package remains at:

`/Users/fuzz/LevelScope-UI-monitoring-2026-11-05/candidate-01/LevelScope Standalone.app`

Measurement build: `4a44fc52dd01be5947cc5627041391dac42ac5e2ee76044c60dc2d8f0faa3b9a`.
App executable SHA-256: `62c76abdbaa2ec0ed7a299d56eec4374de34910b2b7eeca9b8a4232bfde0e55a`.
The executable remains the original tested product revision, not a rebuilt publication revision.

The combined candidate passed **7,270 focused checks on each Mac**: Sequoia 15.7.9
arm64 and actual High Sierra 10.13.6 Intel. These comprise 445 device/routing/volume,
4,059 injected audition/smoothing, 946 preparation/workflow and 1,820 session-playback
checks. Earlier slice-specific waveform/timeline and unchanged M19/Intel, full-duration
and video evidence are reused. No physical opening or audible test occurred for these
slices or publication. Publication only rechecked identities and preservation.

The earlier `build/ui-device-settings-2026-11-05/candidate-01` and
`build/ui-preparation-feedback-2026-11-05/candidate-01` packages remain intact, together
with their Golem copies. Publication rehashed all **397 + 397 + 404** indexed local
artifacts, all three identity files, and the **68** original evidence files. All **383**
combined candidate source-manifest entries match the published checkpoint. No frozen
package, evidence file, original media, installed app or audio setting was changed.
Golem was not accessed for publication.

All 24 pre-existing annotated tag objects and targets and seven unrelated remote branch
refs were preserved. The UI branch remains separate; logic/video worktrees were not
switched or edited. Existing Finder files were left alone. Generated raw reports retain
their original line endings and trailing whitespace; source/document checks pass.

See [publication verification](../Tests/Reference/results/ui-convenience-publication-2026-11-06/checkpoint-verification.json),
[connection/preparation details and manual checks](UIDeviceAndPreparation-2026-11-05.md),
and [monitoring details and manual checks](UIMonitoringMemory-2026-11-05.md).

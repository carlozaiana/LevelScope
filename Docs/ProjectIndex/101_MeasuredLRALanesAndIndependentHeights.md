# Measured LRA lanes and independent timeline heights

Authorized implementation starts from accepted UI handover `583bb6114728483f28a11a62b5c77a669d0433c5`, including implementation `5537a28bff061b90e83c675d15e4c6cfc080fe76` and Source tree `9cc787a660e0f26b87a4ae2dd82365ab2b491ce7`. The user accepted that build as working as expected. This acceptance is separate from automated qualification. Main/origin main is `b5ff52b9c87422b28f27a089c66b746fd2447a05`; the accepted navigation fixes are still local. This publication gap does not change the implementation baseline. No publication or main merge is authorized for this slice.

UI owns timeline presentation, independent height controls, view preferences and native integration on `codex/lra-lanes-resizing-2026-11-24`. Logic owns `Core/Proposal/LraTimelineDiagnostics`, the `StandaloneDiagnostics` service contract and diagnostic calculations/tests on `codex/measured-lra-diagnostics-2026-11-24`, from the same accepted handover. The integrated Logic handover is `68a2cd1f50ed8d651957e484636421dcea301f4e`, retained as an ancestor. Classification stays in that service; proposal and measurement mathematics are unchanged.

Every displayed timeline lane gets a direct vertical divider, including the five existing proposal-context bands, local-gain lane and three measured LRA lanes. The shared waveform/loudness graph has its own divider. The outer inspector divider controls available timeline viewing space; vertical overflow is scrollable without shrinking other chosen lane heights. Reset layout restores defaults. Sizes belong to local view preferences, separate from programme/session data, gains and measurement identities. Source/Current time alignment, horizontal navigation, section editing, Enter-only gain confirmation and outside-click cancellation remain intact.

LRA lanes use completed, matching measured diagnostics only:

- Potential high-side passages associated with programme LRA above the selected ceiling are an inspection heuristic, never a section compliance verdict or a guarantee of correction.
- The upper percentile rank tail distinguishes values strictly above P95 from tied boundary values assigned to excluded ranks. It is not exactly 5% of programme duration and is unrelated to clipping.
- Quiet material below P10 is separate from absolute/relative-gate exclusions, zero energy and insufficient-duration windows. Low percentile rank is not proof of inaudibility. Global gain can change absolute LUFS without changing relative rank.

Logic defines exact window support, EOF padding, eligibility, ties, thresholds and original-frame mapping. The UI queries its bounded index when measurement, viewport, target or width changes; paint uses cached columns. Unverified Current, unresolved Source, provisional measurement and insufficient/unstable data must not display stale completed diagnostics as valid.

The running unsaved user session, original media/sessions, manuals, installations, existing tags and frozen test builds are protected. The boundary-behaviour preference and automatic selection-mode switching from document 100 remain separate later work, as do near-miss application, rule redesign and margin defaults.

Qualification is focused on new data categories and layout/interaction, plus directly affected regressions on the identical integrated package on Sequoia and actual High Sierra. Do not repeat unaffected historical suites in both chats. No physical device or audible testing without coordination. Broad video/full-film/three-hour endurance remain deferred with no automatic 90-minute substitute. The inherited completion pause remains separate. Limited occasional dropped-frame acceptance does not waive audio interruptions, accumulating drift or sustained freezes.

Implementation details: the main and Sections timelines share Logic’s read-only LRA index, with independent view heights. Pixel summaries are cached; exact window inspection lists clipped half-open audio sample support as well as timecode and LUFS. Absolute and relative gate exclusions occupy separate strips; tied upper ranks are hatched. The report remains readable across routine refreshes and expires when its measurement/domain/target changes. Source and Current diagnostics never use a proposal prediction as measured authority.

Heights are remembered separately for main and Sections in the application’s local `LevelScope/View Settings/timeline.json`. Drag a divider or use its arrow keys; double-click resets that height. Session / tools → Reset layout resets all view heights. The existing outer divider changes the visible timeline area, with a vertical scrollbar for overflow. The waveform and loudness curves remain one overlaid plotting area.

Focused development runs are preserved, including a renderer test sampling the playhead overlay instead of its intended gate strip, a test-only null Undo-journal dereference before the first edit, and report-retention failures. Those fixtures and product report handling were corrected before candidate qualification; original receipts/crash report remain available. No live user session was used.

## Preserved implementation and build

- UI implementation: `fbd650cf78e8dc8b36bfaee6d8626a9aec816e1f`.
- Candidate-02 checkpoint, including the coordinate-fixture correction: `61c93ae86d49e4d35bfec798f6568344f5cb8fc5`.
- Production Source tree: `97492cc70c3ef845dbfcbe65e946413380e7f4f6`.
- Measurement build: `86753c73c79a7096f240a0beda16d792517bfe8edcda7861d6f6b155b8e8cd9a`.
- Candidate-02 identity SHA-256: `843f513b6f0afed08b13eda0ca13692721b2506aa65df5dfaddd500a54b5cab6`.
- App executable SHA-256: `e4856a517a4e4d93655d544600005764794906b08b5ff794e4dc37b1781bd368`.
- Source manifest SHA-256: `3c63053f40e16f67f5188f3fed8309e2a52ea93931097a1e1c15f11a4be8f534`.
- Preserved app: `/Users/fuzz/.codex/worktrees/a411/LevelScope/build/lra-lanes-resizing-2026-11-24/candidate-02/LevelScope Standalone.app`.

The app and all three test executables contain arm64 and x86_64; the Intel minimum is macOS 10.13. Candidate 01 is preserved. Candidate 02 has identical production Source, measurement identity and app bytes; only `LevelScopeTimelineProgressTests` changed. Its section gesture fixture previously reused main-window pixel coordinates after the Sections parent established a different plot width. The corrected fixture maps the same exact source samples through each actual width, retaining all original boundary/audio assertions. Candidate 01’s failed receipt and the passing corrected development run are retained.

## Focused qualification

All six candidate-02 groups pass on **each actual platform**: Sequoia 15.7.9 arm64 and High Sierra 10.13.6 x86_64, **2,762 assertions per host**, with identical preserved package identity. The runner verifies every packaged artifact before and after execution.

| Group | Assertions per host |
| --- | ---: |
| Measured LRA data contract (`--lra-timeline`) | 62 |
| LRA native lanes/layout/validity (`--lra-timeline-ui`) | 225 |
| Gain-entry/navigation and popup viewport | 221 |
| Gain entry, selection and manual editing | 1,079 |
| Graph/applied/section boundary interactions | 943 |
| Current initialization and sections-from-regions native workflow | 232 |
| Total | 2,762 |

The new native checks include all ten independent height controls, separate local main/Sections preferences, Reset layout, vertical overflow, exact horizontal range, category strips, cached painting, exact shared click mapping including exclusive EOF, copied-report retention/invalidation, Original versus edited Current, unavailable/in-progress states, save/reopen and independent remeasurement. View changes preserve processing, measurement objects, Undo count, source range and stopped/disconnected transport. Existing gain-entry/popup/graph checks retain the accepted interaction semantics.

Only device-free/injected-output tests ran. The source-integrity/HFS+, proposal algorithm and numeric, export and unchanged playback qualification are inherited, not rerun or extended by this UI evidence. No physical-output, latency, independent-programme correction, listening or endurance claim is added. The new lane build awaits user acceptance. The previously accepted navigation fixes remain recorded separately.

Evidence is under `Tests/Qualification/LraTimeline/results/2026-11-24`, with frozen manifests, raw receipts, screenshots, development failures and preservation audit. Local and remote main remain `b5ff52b9c87422b28f27a089c66b746fd2447a05`; all 31 existing tags and 7,864 protected files across the prior packages/media/manuals/evidence match the preflight. This branch and candidate remain unpublished.

## Short manual check

1. Open a measured programme, set an LRA ceiling and inspect the three lanes below the gain lane. Scroll vertically if needed. Compare Original and Current; click a coloured passage to read its measured windows.
2. Drag the divider directly beneath the waveform/curves, then the gain lane, each visible context band and each LRA lane. Other chosen heights and horizontal zoom/scroll should stay fixed. Repeat in Sections / gain.
3. Change the target after inspecting a lane: the old report expires and the guidance updates. Change a Current gain deliberately with Enter: Current diagnostics disappear until independent remeasurement; Original stays available.
4. Check ruler seeking, section/applied-boundary dragging and Undo/Redo. An unconfirmed gain is cancelled by clicking elsewhere, and that click performs its normal action.
5. Reopen the app to check remembered heights. Double-click a divider to reset that height, or use Session / tools → Reset layout for all heights.

Broad video/full-film/three-hour endurance remain deferred without a 90-minute substitute. The completion pause is still a separate limitation. No publication/main integration is included in this slice.

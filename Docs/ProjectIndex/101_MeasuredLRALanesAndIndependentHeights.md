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

Frozen identities, final two-host results and manual checklist: pending qualification.

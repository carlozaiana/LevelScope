# Gain entry and section-boundary navigation

This UI follow-up starts from handover `b5ff52b9c87422b28f27a089c66b746fd2447a05` and its preserved candidate 02. The user's open, unsaved Tutti session is not a test fixture: do not close, reload, replace, save over, start sound in, or claim to repair that running document. Reproduce in owned test sessions. Main and existing tags remain unchanged; publication needs later approval.

## Confirmed interaction rule

The user's clarification supersedes the earlier click-away-to-apply behaviour. **Enter confirms a valid numeric local gain. Escape or clicking elsewhere cancels every unconfirmed entry**, whether valid or invalid, restoring the actual applied value. For a heterogeneous selection this restores `Mixed` and each section's own gain. The outside click must perform its normal action immediately, including seek, selection or starting a boundary drag. Cancelling changes no processing, measurement identity or Undo history. Invalid Enter must show concise validation and must never treat blank/whitespace as zero. Confirmed Enter remains a normal undoable manual edit requiring remeasurement. An existing explicit gain Apply control may also confirm; no additional main-window control is required.

Apply this to the main gain lane, Sections gain entry, generic pending-entry helpers, dropdown changes and Propose. Numeric typing remains numeric; Space while a numeric gain editor is focused must not accidentally replace selected text. Dialog completion should return focus to the graph for transport shortcuts.

## Boundary-popup correction

A section's left edge is stored as the preceding section's end. That implementation detail must not select the preceding section or reveal its distant beginning. Keep the intended adjacent section selected and retain the exact visible range for Move gain, Keep gain and Cancel, including delayed refresh. Preserve the existing underlying annotation/gain transaction, shared-edge constraints and Undo/Redo. Explicit Fit/Zoom/seek remains navigation.

## Deferred user requests — plan only

1. Add a **general settings-menu preference**, not a main-window control, for section-boundary behaviour: **Ask each time / Always move gain with boundaries / Keep gain in place** as appropriate. The popup should offer **Always use this behaviour**, so the first choice can establish it. Make the preference discoverable and resettable later; retain validation and Undo when confirmation is suppressed. No preference, default change or checkbox is implemented in this fix.
2. After Propose actually applies processing regions, automatically choose **Select applied regions**. Diagnostic, rejected and merely inspected proposals must not trigger it. After **Create sections from applied regions**, choose **Select sections**. Actual insertion already selects sections in the existing code; an already-matching no-op currently retains its previous mode and multi-selection. The future work includes defining and implementing the no-op mode transition without confusing selection intent. No new automatic mode switching is implemented in this fix.

LRA diagnostic lanes, near-miss starting points, acceptance-rule redesign and margin-default work remain separate. Both manuals stay unchanged. Broad video/full-film/three-hour endurance remain deferred with no automatic 90-minute substitute. The inherited completion pause remains a separate limitation; occasional dropped video frames do not waive audio interruptions, accumulating drift or sustained freezes.

## Implementation and qualification

Implementation `5537a28bff061b90e83c675d15e4c6cfc080fe76` on `codex/timeline-input-boundary-fixes-2026-11-24` has production Source tree `9cc787a660e0f26b87a4ae2dd82365ab2b491ce7`. It retains the complete preceding Current-initialization/section-copy slice. Document 99's requirement to finish a pending entry before copying sections is superseded by the explicit cancellation rule above.

The lockout was caused by generic navigation calling the gain-entry commit helper, which returned false for whitespace. Navigation now cancels unconfirmed text. Both numeric editors also receive outside mouse clicks even when the clicked control retains keyboard focus. Enter is the only normal text-confirmation path. The dialog is destroyed before focus returns to the graph, and Space in a deliberately focused gain field cannot erase its selected number. This addresses the reported focus route without claiming how the user's original whitespace was entered.

A left boundary still edits the preceding section's end in the existing model. The UI separately retains the intended selected section and captures/restores the exact visible range across completion. Cancel retains selection; the Sections panel no longer interprets the cancellation's zero ID as a request to select the first section. The model transaction and processing mathematics are unchanged.

Preserved app:
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/timeline-input-boundary-fixes-2026-11-24/candidate-01/LevelScope Standalone.app`

- Package identity SHA-256: `37a12292261e3c5dba77438ac649403cef9f876d829b4f4b20c7746df4a1d42e`.
- Measurement build: `182b0c4836d677c328d2e6457bbe95dc42a9f3977852e115976478e83c3429eb`.
- App executable SHA-256: `ebcf678ff684c134c2e49ebe60b40780758c38023c20b771fbf7e85093665b4e`.
- Universal arm64/x86_64; Intel minimum macOS 10.13 verified in Mach-O. Qualification checks all 473 frozen artifacts before and after each host run.

| Focused group | Sequoia 15.7.9 arm64 | High Sierra 10.13.6 x86_64 |
| --- | ---: | ---: |
| Native gain input, click-away, boundary dialogs and focus | 221 | 221 |
| Gain editor and manual workflow | 1,079 | 1,079 |
| Graph boundaries, sections, navigation and batch editing | 943 | 943 |
| Current initialization and section-copy native regression | 232 | 232 |
| **Passed** | **2,475** | **2,475** |

The same frozen candidate passed all four groups on each actual Mac. Coverage includes valid/invalid/blank/whitespace entries, Mixed selection, Enter/Escape/blur, a real native mouse click on a button that retains keyboard focus, first-click seek/select/pan/drag, dropdown and Propose cancellation, both edges with Move/Keep/Cancel in main and Sections views, exact viewport and selection after refresh, popup return focus, Space, Undo/Redo, gain/envelope and measurement retention or invalidation, save/reopen and the preceding mono/stereo/5.1 workflow. Screenshots from both hosts were inspected. The existing processing tests establish unaffected audio semantics; no new physical output or audible test ran.

[Combined receipt and artifact hashes](../../Tests/Qualification/GainNavigation/results/2026-11-24/acceptance.json) · [Reproduction instructions](../../Tests/Qualification/GainNavigation/results/2026-11-24/README.md).

Three development fixture failures are preserved. Native 01 incorrectly required nonempty gain text after Keep gain, although a section spanning different local gains correctly displays Mixed. Native 02/03 sent mouse events to an obscured owned window; macOS hit testing dispatched neither the button nor the outside-click listener. Positioning the owned fixture visibly fixed the test. Native 04 and all frozen runs passed; no frozen candidate failed. These were test-harness corrections, not weakened processing or navigation checks.

The preservation audit verified 7,312 protected files, all 30 existing tag objects/targets, and unchanged main `c29b19940253d615ba4c3d3909eafa9b32bcc2f2`. The preceding candidate 02, other packages, original media/sessions, evidence and manuals remain intact. The user's running unsaved Tutti session was never operated on. User acceptance of this new build remains pending and separate from automated results. No publication or main integration occurred.

Unaffected model, proposal, measurement, mono, source-integrity/HFS+, transport and export evidence from documents 91–99 is reused. There is no new independent-programme correction, listening, physical-output latency or endurance claim. All deferrals and the separate completion-pause limitation above remain.

## Short manual check

1. Type a different local gain, then click the ruler or start dragging a boundary without Enter. The old gain returns and the click/drag works immediately. Repeat with invalid text and a Mixed selection.
2. Type a valid gain and press Enter: it applies and Current needs remeasurement. Undo/Redo should remain normal. Escape cancels unconfirmed text.
3. Zoom closely around a selected section's left and right edges. Try Move gain, Keep gain and Cancel. Selection, zoom and scroll should remain steady in both timeline windows.
4. After the boundary dialog closes, Space should target Play/Stop rather than enter whitespace into the gain field. Test playback only with your normal chosen monitoring setup.

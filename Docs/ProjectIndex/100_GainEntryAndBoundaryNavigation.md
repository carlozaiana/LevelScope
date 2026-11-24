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

In progress. The final handover will identify the separate preserved app, source and package identities, focused native checks on each available actual platform, preserved failure receipts and remaining limitations. No physical output or audible testing is included.

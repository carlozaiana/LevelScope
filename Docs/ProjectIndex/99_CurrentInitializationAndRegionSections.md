# Current initialization and sections from applied regions

Propose can now start when Current has been cleared. Its explicit click creates unity Current from the confirmed Original, shows preparation, and continues the existing proposal and independent verification workflow. An existing Current always remains the starting point, including unmeasured manual gains. Valid, freshly accepted identical Source measurements may initialize unity Current through the existing model semantics; a restored Source cache alone cannot grant Current verification. Cancelling after initialization retains that unity Current and says so.

The gain lane has **Create sections from applied regions**. It copies actual differing-gain nominal cuts, including hand-adjusted cuts, into the section map once. Existing cuts, IDs and names are retained on the original left-hand pieces. It selects the first new section for the existing absolute local-gain editor. A matching map reports no change and preserves the selection. A dirty gain entry must be completed or cancelled first. No local regions, busy work, stale input or the section limit have brief blocker tooltips.

Creating sections changes annotations only. It retains exact processing and ramp support, historical explanations and valid measured Current, including PASS. One Undo/Redo changes only those annotations. The model and operation bind a monotonic runtime section-map generation, so an old proposal or replacement choice cannot regain authority after copy/Undo. Save/reopen preserves sections and processing; runtime replacement authority is not restored. The copy is not an ongoing link: moving an applied boundary later does not automatically move a section.

The gain editor displays absolute local dB, before global trim. Entering -1 dB in a section currently at -2 dB replaces the local setting; it does not add another -1 dB. A genuine gain or gain-boundary edit immediately requires independent Current remeasurement.

## Integration and evidence

Baseline main c29b19940253d615ba4c3d3909eafa9b32bcc2f2; Logic transaction handover 389e5154b0e294d1b9ab65ffd33d921533f9cbf7 is retained as an ancestor. UI branch codex/current-and-region-sections-2026-11-24. The UI implementation is `b766aaedb7694a5e23deb058fbe886af3e5c870f`. The tested checkpoint is `0fdb7d54d88d82dcab67bd3b90d46d3aa81ca891`; its final change updates the first-import regression to expect initialization before Source preparation. Both have production Source tree `9bde2e81e439c6a85f035320d0a26c5b33e20629`.

Preserved app:
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/current-and-region-sections-2026-11-24/candidate-02/LevelScope Standalone.app`

- Candidate identity SHA-256: `8203a98b4ef7481ae67557a77ccd851a0047b66f07edf3bd5a5b9cbe60303c12`.
- Measurement build: `aff7de74907a71958e3a8d11383a6e0a03908dd96c2839322741ef373bff4c5a`.
- App executable SHA-256: `485cdcaa27cb0310b727d097807e7284dc614250bbef4b6905e8146b610775cd`.
- Universal arm64/x86_64 app and tests; Intel minimum macOS 10.13 verified from Mach-O load commands. The exact package's 473 artifacts were checked before and after each host run.
- No session schema change. The source-derived measurement build advances; saved measurements retain the normal build/identity checks. Runtime proposal authority and the map generation are not restored from saved historical explanations.

| Affected group | Sequoia 15.7.9 arm64 | High Sierra 10.13.6 x86_64 |
| --- | ---: | ---: |
| Current initialization and applied-region section model | 77 | 77 |
| Existing section-processing smoke | 334 | 334 |
| New native main-window workflow, mono/stereo/5.1 | 232 | 232 |
| Graph boundary, section editing and navigation regression | 943 | 943 |
| Proposal orchestration, alternatives, cancellation and stale/retry guards | 187 | 187 |
| **Passed** | **1,773** | **1,773** |

The new native group covers cleared and manually edited Current; real blockers; cancellation during initialization and preparation; fresh Source versus a reopened cache; same-session/runtime guards; graph-adjusted nominal cuts; existing named annotations; selection and repeat no-op; retained PASS; annotation Undo/Redo; actual processed float samples; save/reopen; dirty gain entry; absolute local gain editing and independent remeasurement. The model checks additionally include an existing editorial cut six frames into a transition, map limits, stale/busy guards and proposal authority through immediate Undo.

[Combined receipt and hashes](../../Tests/Qualification/CurrentSections/results/2026-11-24/acceptance.json) · [Reproduction commands](../../Tests/Qualification/CurrentSections/results/2026-11-24/README.md).

Two test-fixture failures are preserved: the first new native run read the gain-entry state before delivering its queued text callback; candidate 01's existing operation test expected Source preparation before the newly visible initialization phase. Corrected fixtures passed, with unchanged production Source and an identical app executable. Candidate 01 and all development receipts remain available. Native screenshots on both hosts were inspected for the command, section gain editor and measured status.

The preservation audit checked all 30 existing tag objects/targets and 6,292 protected files. Main remains `c29b19940253d615ba4c3d3909eafa9b32bcc2f2`. No release refs or other checkouts were changed. New user acceptance is pending and is separate from these automated results.

The model review covered the initialization order, map generation binding, cancellation messaging and no-op selection preservation. The inherited operation test's old mono-rejection assertion was updated to the supported mono LRA pass qualified in document 97; method dispatch and mathematics are unchanged.

## Scope retained

Documents 91-98 retain authority for unaffected DSP, proposal methods, mono capability, source reuse and integrity (including HFS+), transport, output routing and export. This slice runs affected annotation/model/native checks; it adds no new physical-output, listening, independent-programme correction or endurance claim. Both manuals and all existing packages are untouched. Main and publication remain for separate approval.

Broad video, full-film and three-hour endurance remain deferred, with no automatic 90-minute substitute. The inherited completion pause remains separately documented. Limited acceptance of occasional dropped video frames does not waive audio interruptions, accumulating drift or sustained freezes.

## Short manual check

1. With a confirmed mix, choose Session / tools → Clear Current State, then click Propose. Watch preparation and the final measured result. Cancelling during preparation retains unity Current with an explicit message.
2. Make a deliberate manual gain edit and click Propose. It prepares the actual edited Current without reverting to Original.
3. With local applied regions, optionally drag an applied cut, then click Create sections from applied regions. New section cuts follow nominal applied cuts. If Current was already measured PASS, it stays PASS; sound is unchanged.
4. Repeat the command, Undo and Redo, then save/reopen. Repetition should make no further section or selection change.
5. Select a copied section, enter an absolute local dB value and press Return. Current needs remeasurement. Re-measure explicitly, and listen around the edited transitions using your normal monitoring setup.

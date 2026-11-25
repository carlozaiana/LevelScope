# Navigation, measured LRA lanes and independent resizing — 25 November 2026

The user accepted the integrated candidate: **“Tested and works as expected.”** After review of qualification, ancestry and the exact tested Source, the user explicitly authorized publication and main integration: **“Ok please tell the UI chat to do the publication.”** Coordinating chat `01a13c8f-90a9-7150-813f-cb057626ef7e` relayed that authorization to the UI owner. User acceptance is separate from automated results; no additional listening or individual checklist outcomes are inferred.

This combined checkpoint includes the previously unpublished Enter-only gain-entry / outside-click cancellation and section-boundary navigation fixes, Logic’s measured LRA diagnostics, and the native lanes with independent graph/lane resizing. It closes the publication gap recorded in [100](100_GainEntryAndBoundaryNavigation.md) and [101](101_MeasuredLRALanesAndIndependentHeights.md). Those documents and their original qualification receipts remain unchanged as historical evidence; this publication record supersedes their publication/acceptance status.

## Published identities

- Published [development branch](https://github.com/carlozaiana/LevelScope/tree/codex/lra-lanes-resizing-2026-11-24) and [main](https://github.com/carlozaiana/LevelScope/tree/main).
- Annotated [checkpoint `m20-navigation-lra-lanes-resizing-2026-11-25`](https://github.com/carlozaiana/LevelScope/tree/m20-navigation-lra-lanes-resizing-2026-11-25).
- Checkpoint target: `d3fd6b3ab81bc51272f8003c3c969a9a09a4c31b`.
- Annotated tag object: `2320cb775fa4ad81db8290e8b0b4000e04d588f4`.
- Included navigation handover: `583bb6114728483f28a11a62b5c77a669d0433c5`.
- Included Logic handover: `68a2cd1f50ed8d651957e484636421dcea301f4e`.
- UI implementation: `fbd650cf78e8dc8b36bfaee6d8626a9aec816e1f`.
- Frozen candidate-02 checkpoint: `61c93ae86d49e4d35bfec798f6568344f5cb8fc5`.
- Two-platform qualification handover: `670e987fc4bc84eaa4d5f222ea8b150d5cb9ed7b`.
- Tested Source tree: `97492cc70c3ef845dbfcbe65e946413380e7f4f6`.
- Candidate identity SHA-256: `843f513b6f0afed08b13eda0ca13692721b2506aa65df5dfaddd500a54b5cab6`.
- Measurement build: `86753c73c79a7096f240a0beda16d792517bfe8edcda7861d6f6b155b8e8cd9a`.

Preserved app:

`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/lra-lanes-resizing-2026-11-24/candidate-02/LevelScope Standalone.app`

The combined checkpoint is published and integrated into main, closing the earlier navigation publication gap. Main was still `b5ff52b9c87422b28f27a089c66b746fd2447a05`, with no intervening work. The existing dedicated integration checkout at `/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m18-publication-2026-11-03/integration` performed a normal fast-forward to `d3fd6b3ab81bc51272f8003c3c969a9a09a4c31b`. The branch and new tag were published atomically, read back, then main was pushed without force. No history was rewritten and no existing tag moved.

This following documentation/evidence-only receipt commit advances main and the development branch; the annotated checkpoint stays fixed on `d3fd6b3`. Final remote readback and the exact final-main identity are retained externally in `build/lra-navigation-publication-2026-11-25/remote-final.txt` and `final-verification.json`, avoiding a self-referential commit or retagging. No tested app was rebuilt, replaced or repackaged.

## Publication verification

The [preflight](../../Tests/Qualification/LraTimeline/results/publication-2026-11-25/preflight.json) verified the frozen candidate’s **471 committed source-manifest entries**, exact production Source/CMake match, qualification receipts, ancestry and clean integration checkout. The [checkpoint readback](../../Tests/Qualification/LraTimeline/results/publication-2026-11-25/remote-checkpoint.txt) and [publication audit](../../Tests/Qualification/LraTimeline/results/publication-2026-11-25/checkpoint-verification.json) passed **18,331 identity/integration/preservation checks**. These are publication checks, separate from the runtime assertions below.

All **31 previous tag objects/targets**, **35 local candidate packages** and **18,151 inventoried protected files** remain unchanged, including inherited external packages, manifests, qualification receipts, original media/sessions, listening materials and manuals. Other worktree heads, branches, existing edits and untracked files remain preserved. The user’s live session and installations were not operated on. The new tag is the only tag addition.

[Publication receipt](../../Tests/Qualification/LraTimeline/results/publication-2026-11-25/publication-receipt.json) · [User acceptance and authorization](../../Tests/Qualification/LraTimeline/results/publication-2026-11-25/user-acceptance.json).

## Qualification and limits

The identical candidate passed **2,762 focused assertions on each actual platform**: Sequoia 15.7.9 arm64 and High Sierra 10.13.6 x86_64. The six groups cover the measured LRA contract, independent heights and native inspection/expiry, gain input/navigation, manual gain entry, graph/section boundary editing, and Current initialization/sections from applied regions. Unchanged qualification is reused as documented in [101](101_MeasuredLRALanesAndIndependentHeights.md). Publication requires only focused identity, ancestry, integration and preservation checks; no source change requires repeated runtime testing.

Broad video, full-film and three-hour endurance remain deferred, with **no automatic 90-minute substitute**. The inherited completion pause remains separate. Limited acceptance of occasional dropped video frames does not waive audio interruptions, accumulating drift or sustained freezes. No new physical-output, listening, latency or endurance claim is added. Publication requires no Golem access, audible test or interaction with the user’s live session.

Boundary-behaviour preferences, automatic mode switching, measured near-miss application and proposal-rule redesign remain separate later work. This publication starts no new feature.

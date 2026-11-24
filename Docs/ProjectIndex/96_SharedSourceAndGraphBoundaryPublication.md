# Shared source and main-graph boundary publication — 24 November 2026

The qualified shared-source/background-playback work and main-graph applied-region
editing are published together and integrated into main. The user accepted the
graph interaction: **“Works as expected, please merge into main.”** That acceptance
is recorded separately from automated checks, without inferring additional
listening results or completion of individual checklist items. This publication
record supersedes the pre-publication status in documents 94–95; their original
qualification evidence remains unchanged.

The application shares a verified source asset across session consumers and
prepares playback in the background while stopped. Current inspection can now
drag actual applied-region boundaries directly in the waveform/curves graph,
using the gain lane's existing processing transaction, Undo and invalidation.
An edited Current can be auditioned before independent remeasurement; playback
does not confer a numerical PASS. Neither preparation nor editing starts sound.

## Published identities

- Annotated [checkpoint `m20-shared-source-graph-boundaries-2026-11-24`](https://github.com/carlozaiana/LevelScope/tree/m20-shared-source-graph-boundaries-2026-11-24).
- Checkpoint/both-platform handover: `2590d3fbdaaf437e313f8646021039d3467c4ef2`.
- Annotated tag object: `6cf0a897b673817c1d4fe419e17413493ac0a6ac`.
- Graph implementation: `f801ac06307f25f423a953a7e35b59c28013ccc7`.
- Included shared-source handover: `f9eeae554421b4ef4f8b48d53a69785aa3f6af46`.
- Included shared-source candidate-03 implementation and HFS+ integrity fix:
  `614377f4e573bdd0fc34726584418b9333d092c2`.
- Tested production Source tree: `31be38bd890fee3667a6a976fbc3524a2b26e609`.
- Published [development branch](https://github.com/carlozaiana/LevelScope/tree/codex/graph-region-boundaries-2026-11-24) and [main](https://github.com/carlozaiana/LevelScope/tree/main).

GitHub main was still `18adeb520006b7fe7e0a22300cc26abd28eec51c`, an ancestor of
the checkpoint; there was no intervening main work. The existing main integration
checkout at
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m18-publication-2026-11-03/integration`
performed a normal fast-forward. Main, the development branch and the new tag
were published atomically without force-push, history rewriting or moving any
existing tag. Other chats' checkouts and edits were preserved.

This following publication-receipt commit contains documentation/evidence only.
Main and the development branch advance to it, while the checkpoint tag stays
permanently on `2590d3f`. Final remote readback and its containing commit identity
are retained externally in
`build/shared-source-graph-publication-2026-11-24/remote-final.txt` and
`build/shared-source-graph-publication-2026-11-24/final-verification.json`, avoiding
a self-referential commit or retagging.

## Preserved applications

The current combined test app is unchanged:

`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/graph-region-boundaries-2026-11-24/candidate-01/LevelScope Standalone.app`

- Package identity SHA-256: `cecac62f3f98a5c92e1b87d938ab94f813b32f7cd2984cf8e7156a7bf7af43cf`.
- App executable SHA-256: `7e276e193f2f553f87f0f9ae191866f86d10699caf7612282b408afba0776ba5`.
- Source manifest SHA-256: `460359c1b4d24cc61c1a9be3176933b8ff725b96af1925718a0341ab5d42556c`.
- Measurement build: `fcb60b72c950d8b21e3c8ae4be8ce57bf95f6991c27aee84073d4d092bd78ea7`.

The accepted shared-source baseline remains at:

`/Users/fuzz/Developer/LevelScope/build/session-source-2026-11-23/candidate-03/LevelScope Standalone.app`

Its package identity is
`94c42f1eca24de3adc671723b4b01cd4bcc986fe602b4c6c6d25484379079843`.
No frozen app was rebuilt, replaced or repackaged. Shared-source candidates 01
and 02 and their earlier failure evidence also remain preserved.

## Verification and scope

The [publication preflight](../../Tests/Qualification/GraphBoundaries/results/publication-2026-11-24/preflight.json)
verified all **456 committed source-manifest entries** against graph candidate 01.
The checkpoint and following receipt contain the exact tested production Source
and CMake configuration. The
[remote checkpoint readback](../../Tests/Qualification/GraphBoundaries/results/publication-2026-11-24/remote-checkpoint.txt)
and [5,837-check publication audit](../../Tests/Qualification/GraphBoundaries/results/publication-2026-11-24/checkpoint-verification.json)
passed: all **28 existing tag objects/targets**, unrelated remote branches,
**15 packages** and **5,666 inventoried protected files** are unchanged. This
includes source manifests, qualification receipts, original media/sessions,
listening material and manuals. Other worktree heads, branches and local edits
are preserved; Finder metadata differences, if any, are reported rather than reset.

Graph candidate 01 passed **1,506 focused assertions on each actual platform**,
Sequoia 15.7.9 arm64 and High Sierra 10.13.6 Intel, using the identical universal
package. See [95](95_GraphAppliedRegionBoundaries.md) and the preserved
[both-platform receipt](../../Tests/Qualification/GraphBoundaries/results/combined-2026-11-24/acceptance.json).
Unchanged shared-source, processing, proposal and transport paths inherit
candidate 03's qualification from [94](94_SharedVerifiedSourceAndPlaybackPreparation.md):
11 groups and 227,322 assertions per actual platform, including the HFS+
same-timestamp write correction. These are separate evidence sets, not a new
combined whole-suite run. Publication performed only identity, ancestry,
integration and preservation checks; no source change required more runtime tests.

The [user acceptance](../../Tests/Qualification/GraphBoundaries/results/publication-2026-11-24/user-acceptance.json)
is separate from those automated results. Earlier immediate Play/seek/replay
feedback on a 93-minute 5.1 file remains interaction feedback, not an endurance
run or physical-output latency benchmark. Publication opened no physical output,
started no audible test and did not access Golem.

Broad video, full-film and three-hour endurance remain deferred, with **no
automatic 90-minute substitute**. The documented completion pause remains a
separate limitation. Limited acceptance of occasional dropped video frames does
not waive audio interruptions, accumulating drift or sustained freezes.

No new feature milestone starts in this publication task. Proposal strategy,
search limits, defaults and measurement mathematics are unchanged.

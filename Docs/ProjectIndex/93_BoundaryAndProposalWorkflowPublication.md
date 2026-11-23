# Combined boundary editing and Propose workflow publication — 23 November 2026

The user explicitly authorized publishing the accepted applied-region boundary
editing and streamlined Propose/apply/independently-verify workflow together and
integrating them into main. Publication is complete. This record supersedes the
pre-publication branch/main status in documents 91–92 without changing their
original qualification evidence.

## Checkpoint and implementation

- Annotated [checkpoint `m20-boundary-editing-propose-verify-2026-11-23`](https://github.com/carlozaiana/LevelScope/tree/m20-boundary-editing-propose-verify-2026-11-23).
- Checkpoint/handover commit: `d5f6e86801c675794bb0b37281afc44bd30aa961`.
- Annotated tag object: `7275a2554f96a33867573ccc370078b31569d0e9`.
- Workflow implementation: `79b55eb177a01e56fdbe86f843e03140df356bcd`.
- Included boundary implementation: `240da530302af8f539d0351b7332798e2373314a`;
  accepted boundary handover: `c2af1853c783d559a0414a49c7bb2f280ded7385`.
- Tested production Source tree: `0a96e9ab0b47847d301c7b273577b945a4c1921a`.
- Published [development branch](https://github.com/carlozaiana/LevelScope/tree/codex/propose-verify-workflow-2026-11-22) and [main](https://github.com/carlozaiana/LevelScope/tree/main).

Main was still `1056ad4ace536abfd68904f083c7fc56a1a8999f`, an ancestor of the
checkpoint. There was no intervening main work to resolve. The existing isolated
main checkout at
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m18-publication-2026-11-03/integration`
performed a normal fast-forward. Main, the new development branch and annotated
tag were published atomically, without force-push, history rewriting or moving
any existing tag. Other chats' checkouts and branches were not changed.

This following publication-receipt commit contains only documentation and
evidence. Main and the development branch advance to that receipt commit; the
checkpoint tag remains permanently on `d5f6e86`. The final direct remote readback
and verification are retained in
`build/boundaries-workflow-publication-2026-11-23/remote-final.txt` and
`build/boundaries-workflow-publication-2026-11-23/final-verification.json`.
Those external receipts identify their own containing Git revision without a
self-referential commit or a moved tag.

## Preserved tested applications

The current combined application is the original universal workflow candidate 01:

`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/propose-verify-workflow-2026-11-22/candidate-01/LevelScope Standalone.app`

- Candidate identity SHA-256: `4eb921611c8cec84e98b79ef7cb94aa9fd09da6b27bfbbeed958e35706e68ce3`.
- Source manifest SHA-256: `eb883cde76149331d7016dc8e3d6cc6dcc5aeb0f0bd4263e5014fafe482fd1e4`.
- Measurement build: `9b379f0479f69a6e2ded23dc778ad5c50c5e92a20e4076a5936e06b5a82698fd`.

The separately accepted boundary-editing app remains at:

`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/applied-region-boundaries-2026-11-22/candidate-01/LevelScope Standalone.app`

Its candidate identity remains
`64d5f870f5382f332f58564b945f0ebf61012218570dc1bb970bb06d65171942`.
Neither app was rebuilt, replaced or repackaged. Frozen identities retain their
original build-time status fields; subsequent qualification, user acceptance and
publication authority are recorded separately.

## Verification and acceptance

[Preflight](../../Tests/Qualification/ProposalWorkflow/results/publication-2026-11-23/preflight.json)
verified all **447 committed source-manifest entries** against the frozen
workflow candidate. The checkpoint's production Source and CMake configuration
match the tested implementation exactly. The
[remote checkpoint readback](../../Tests/Qualification/ProposalWorkflow/results/publication-2026-11-23/remote-checkpoint.txt)
and [3,813-check publication audit](../../Tests/Qualification/ProposalWorkflow/results/publication-2026-11-23/checkpoint-verification.json)
passed. All **27 existing tag objects and targets**, unrelated remote branches,
**11 packages** and **3,638 inventoried protected files** remain unchanged,
including qualification evidence, source manifests, original media/sessions,
listening materials and manuals. Unrelated worktree heads, branches and edits
are preserved; Finder metadata, if changed, is reported separately and not reset.

The identical workflow candidate passed **2,957 affected assertions on each
actual platform**, Sequoia 15.7.9 arm64 and High Sierra 10.13.6 Intel. See
[92](92_ProposeApplyVerifyWorkflow.md) and the unchanged
[both-platform acceptance receipt](../../Tests/Qualification/ProposalWorkflow/results/2026-11-23/acceptance.json).
Earlier unchanged M20 and boundary evidence is reused as documented. Publication
ran only identity, ancestry, preservation and integration checks: no new source
change warranted repeated runtime qualification. No physical audio device was
opened, no audible test was started, and Golem was not accessed for publication.

The user reports that the new workflow was tested and works as expected. This
[human acceptance](../../Tests/Qualification/ProposalWorkflow/results/publication-2026-11-23/user-acceptance.json)
is recorded separately from automated checks; no additional listening result or
individual checklist completion is inferred. Prior boundary editing, Undo/Redo
and listening acceptance remains separately recorded. Escape cancels while the
mouse is held during a drag; after release, Undo reverses the committed edit.

This does not add independent-programme correction success or physical-device
qualification. The earlier independent searches found no passing correction and
correctly blocked selection/Apply; those results remain failures to find an
eligible correction, not successful corrections or listening acceptance.

## Retained limits and follow-on work

Broad video, full-film and three-hour endurance remain deferred, with **no
automatic 90-minute substitute**. The known completion pause remains a separate
limitation. Existing limited acceptance of occasional dropped video frames does
not waive audio interruptions, accumulating drift or sustained freezes.

The accepted ordinary workflow is now available: Propose from actual Current,
apply an eligible measured passing plan while stopped, independently verify
Current, inspect or listen, and adjust applied boundaries if necessary. A manual
boundary edit immediately requires remeasurement. Complete explanations,
bounded alternative replacement and Undo/recovery remain available.

No new feature milestone starts in this publication task. Labelled extension
limits and possible default changes, planning-sequence comparisons and selective
positive local gains remain later evaluations; this publication changes no
proposal strategy, mathematics, defaults or search limits.

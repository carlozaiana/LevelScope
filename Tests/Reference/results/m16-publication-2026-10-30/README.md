# M16 publication receipt

Candidate-03 published 30 October 2026. Status: **Sequoia tested; High Sierra, video and three-hour checks pending.** The original qualification identities/status files are historical records and are not rewritten. Publication binds the exact tested manifest to a Git checkpoint; no rebuild or new runtime qualification was required for unchanged implementation.

Source manifest: `3ac2e19d8dd976a8fdc41490245ad9961af91e4bfe1a1d5afdc1ad5c03897f38`.
Measurement/build: `9b7d85d5463297b70e6edaabd5010ddabbd317ea63e2ae1eea7a8280da387bcf`.
Universal app SHA-256: `8ee1b747685a8840284dbe6992b95c1076a76fc90cfe510ac043b63c61df2286`.

The original app remains under `build/m16-qualification-2026-10-30/candidate-03/LevelScope Standalone.app`.
Main remains M13 `a9be070ab9eb5cb4f591ff9d991d8999c3a016d0`. Combined actual High Sierra/video/three-hour qualification follows M18, no later than M19, before any separately approved main integration. Physical full 7.1 listening remains separately pending. No M17 implementation is authorized.

## Published checkpoint

- Annotated tag: [`m16-proposal-preferences-sequoia-2026-10-30`](https://github.com/carlozaiana/LevelScope/tree/m16-proposal-preferences-sequoia-2026-10-30).
- Checkpoint commit: [`73a7b25e71f87299fa2d93c1b1106c118c55beec`](https://github.com/carlozaiana/LevelScope/commit/73a7b25e71f87299fa2d93c1b1106c118c55beec).
- Tag object: `3fb03c9d4ba0374f3b6be9e69b48639587e9bc82`.
- Published branch: `codex/logic-module`; main was not merged or pushed.

All 269 source/test/build files in the checkpoint match candidate-03's manifest;
all 283 frozen artifacts and 88 indexed qualification receipts remain unchanged.
Original freeze-time pending/unpublished fields identify their historical phase;
this separate receipt supplies the later Git mapping. The base commit
`63eafac9acc210d762dfd5e3dcc70defbe5bbdf3` is not the M16 implementation commit.
Original JUCE CRLF textual receipts are preserved byte-for-byte rather than
reformatted; source/document whitespace checks passed.

The initial atomic push received GitHub's `fatal error in commit_refs`; immediate
remote inspection confirmed neither ref changed. Separate normal branch and tag
pushes succeeded. No force push or tag replacement occurred.

Direct GitHub ref verification confirms the branch/tag target, unchanged main and
all **19** original annotated tag objects and peeled commit targets. Other remote
branches and the registered video worktree location/HEAD are unchanged. The GitHub
connector also retrieved the exact published commit. See
[checkpoint-verification.json](checkpoint-verification.json),
[committed-source-verification.json](committed-source-verification.json),
[github-commit-verification.json](github-commit-verification.json),
[before refs](remote-before.txt) and [after checkpoint refs](remote-after-checkpoint.txt).

The documentation-only commit containing this completed receipt updates the handover
and adds [M17's proposal](../../../../Docs/ProjectIndex/60_M17JobLraRequirementsAndPresets.md).
It follows the accepted checkpoint and leaves its tag fixed. The branch may therefore
be ahead of the checkpoint by publication documentation; no source, build, test or
frozen qualification file changes are part of that follow-up. The final remote branch
verification is retained locally in `build/m16-publication-2026-10-30/` after pushing
the receipt, avoiding a self-referential commit identity.

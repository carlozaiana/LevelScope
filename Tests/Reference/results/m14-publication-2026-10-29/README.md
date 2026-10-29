# M14 development checkpoint publication — 29 October 2026

**Sequoia tested; High Sierra, video and three-hour checks pending.**

The user approved publication of candidate-02 on `codex/logic-module`, with main
held at qualified M13 `a9be070ab9eb5cb4f591ff9d991d8999c3a016d0`.
The annotated checkpoint is [`m14-mono-rf64-sequoia-2026-10-29`](https://github.com/carlozaiana/LevelScope/tree/m14-mono-rf64-sequoia-2026-10-29).
Its target commit is `3d98284c2c0b7de92cfee70cb0007d43d1f04d1f`. Annotated tag object:
`8523392cd96b76513f10c8441d236c8f16a9e5a1`. It records the exact implementation in candidate-02 source manifest
`9dd22a598281e83431039f210202a8a24fe49517294cc6bb8d83032c4975e013`.
Publication adds documentation/evidence; it does not change the tested implementation.

Both original packages remain unchanged. The candidate-02 universal standalone
executable SHA-256 is
`188446f5752d91e40e90066d2b25133eaaa07b421ed0018da15a1be2fb60777c`,
byte-identical to the originally exercised candidate-01 app. Measurement build:
`4002f49d4d8eb3cffcf825f1624d8a00d6f48557f0755e711683cff6627f68da`.
The frozen candidate's M13 base revision is not its implementation revision;
this new checkpoint ties the manifest to Git without rebuilding the app.

[Preflight](preflight.json) records all 257 source-file matches, the historical
archive manifest hash, and all 17 prior tag objects/targets.
[Completion recheck](completion-recheck.json) validates both packages and the
complete approved Sequoia evidence. Its `qualified: false` is intentional.
Its pre-review `publication_authorized: false` is retained as tool output; the
subsequent user approval recorded here and in preflight authorizes publication only.
The [qualification archive](../m14-qualification-2026-10-29/README.md) remains
byte-for-byte historical, including review-stage status statements.

No additional build, audio-device test, Golem access or video execution is required
for unchanged implementation publication. No raw media, binaries, credentials,
connection details or unrelated Finder metadata are included. Main, prior tags,
frozen packages, installations, the separate video worktree and audio settings
must remain unchanged. Stop after publication; no next milestone is started.

Publication succeeded and was verified directly through GitHub API refs/tag data
and `git ls-remote`. All 17 prior annotated tag objects and targets, main and other
remote branches were unchanged. [Remote verification](remote-verification.json)
records the exact checkpoint, tag and main identities. A following documentation-only
commit records this receipt without moving the checkpoint tag or changing implementation.

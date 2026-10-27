# M12 development checkpoint — 27 October 2026

**Sequoia tested; High Sierra and video qualification pending.**

The user authorized publication of candidate-03 on `codex/logic-module` with a
dated annotated checkpoint. Main must remain the qualified M11 commit
`cfb1033aa7fb65a2655db276bee696161196bdeb`. This publication is not a main merge or
full runtime qualification, and does not authorize M13 implementation.

Checkpoint:
[`m12-section-trim-sequoia-2026-10-27`](https://github.com/carlozaiana/LevelScope/tree/m12-section-trim-sequoia-2026-10-27).
Its annotation records the exact commit, qualification status and artifact hashes.
The stable tag avoids a self-referential commit hash inside its own files:

```sh
git rev-parse 'm12-section-trim-sequoia-2026-10-27^{commit}'
git cat-file -p m12-section-trim-sequoia-2026-10-27
```

## Frozen candidate and evidence

| Identity | Value |
| --- | --- |
| Candidate | M12 candidate-03 |
| Source manifest SHA-256 | `d6501c0d234ccf9c02b90527fc4af6c613dd4d538a71cc5520fb6b7b21085c1f` |
| Measurement/build identity | `185566e2c3ed9835d9878643618bd6789e785dee78fd670174f1ede4eb10f5fd` |
| Original identity.json SHA-256 | `d0004ae7604b0c86158a623047c07b86a4bf8e0c85f9e462998289ae22621b6a` |
| Original universal app executable SHA-256 | `ab319348a5ca586fd9048581165f7fdefb0cb1bd7221d827c78fd3740ec2aec4` |
| Qualified three-hour RF64 artifact SHA-256 | `a57c70b0a4955556712c80706446060964a2ba8419ddb1f7cc6c6c180afcbcd6` |

[Preflight](preflight.json) records the 231 matching source inputs, 266 unchanged
frozen artifacts, 116 hashed evidence files and all fifteen original local/remote
tag objects and commit targets. Remote main and the development branch both
started at M11. No implementation integration or rebuild was needed.

The preserved tested app is:
`/Users/fuzz/Developer/LevelScope/build/m12-qualification-2026-10-27/candidate-03/LevelScope Standalone.app`.
Its original metadata describes the pre-publication build above M11; the tag maps
that source manifest to the saved commit. A later rebuild is not substituted.
Earlier packages and failed fixture/test attempts remain intact.

The [completion receipt](../m12-sequoia-phase-2026-10-27/completion-candidate-03.json)
and [handover](../../../../Docs/ProjectIndex/46_M12QualificationAndHandover.md)
record nineteen selected short runs, pinned EBU/ITU checks, the 3,304-assertion
M12 matrix, eleven silent UI cases, old-session migration and one three-hour
1,000-interval workflow with verified RF64 export. The long workflow passed 2,067
checks in 70 minutes: 87.3 MiB additional peak RSS, 15.06 ms message-loop p95 and
42.60 ms cancellation. Native measurement covers the complete output; the
independent decoder/oracle uses explicitly bounded sample regions. Fresh and
exact-binary reused receipts retain their original identities.

All actual High Sierra execution, cross-Mac exchange and video-loading/decoding/
playback tests remain deferred until after M13. The locally compiled Intel slice
targets macOS 10.13 but is not High Sierra runtime evidence. Occasional picture
skips remain an accepted development limitation; audio interruptions, drift,
sustained freezes, identity and lifecycle failures remain blockers. Physical full
7.1 listening is separately pending.

## Preservation and publication boundary

Only intended source, tests, manifests, textual evidence and documentation are
committed. Raw media, build products, tested packages, Finder metadata, credentials,
connection records and unrelated work are excluded. No historical qualification
receipt is rewritten to change its pre-publication authorization flags; this
separate record describes the later user authorization.

Push only the development branch and new tag, never main or existing tags.
Verify the remote tag object and peeled commit and all original tag targets.
Keep the video worktree, installations, plugins, Pro40 and system audio unchanged.
After publication, prepare the M13 proposal for review without implementing it.

# M8 checkpoint publication — 25 October 2026

The annotated tag [`m8-reproducible-sessions-2026-10-25`](https://github.com/carlozaiana/LevelScope/tree/m8-reproducible-sessions-2026-10-25)
identifies the accepted M8 checkpoint. Its peeled commit is also the main revision
at publication. The annotation records the exact commit and main revision, the
qualified source manifest, the original binary hashes and the evidence identity.
This stable tag reference avoids a self-referential commit hash inside the commit's
own files. Resolve it with:

```sh
git rev-parse 'm8-reproducible-sessions-2026-10-25^{commit}'
git cat-file -p m8-reproducible-sessions-2026-10-25
```

[Publication preflight](publication-preflight.json) records the original eleven
local and remote checkpoint tag objects and peeled commits, all remote branches,
the registered worktrees and unrelated-file hashes. Latest main was still the M7
base `cfbf3f981dece326d2363b911e5914ea61f1796c`; no newer main changes or integration
fixes were required. Main received the accepted checkpoint by fast-forward.

## Tested candidate → saved checkpoint

| Identity | Value |
| --- | --- |
| Frozen candidate | M8 candidate-01 |
| Source-manifest SHA-256 | `1c591122c03f40cdc7035d20837e4e39f9cf0a1c62de640149c0b6427d225229` |
| Original measurement build | `77cf1c7a5052442d4c938ed82ff1158b2bf6557d28bf99880666ec417f4c7c05` |
| Original package identity SHA-256 | `dd625b0e2241afde0e6cbd5bbe36f733d8e23bdd9afc3bba92cd62fc556e8b9f` |
| Tested universal app executable SHA-256 | `146b9a79e12387a5f1c368e1543f50e581bcd6ceea63cde2111ca4e4d9161740` |
| Checkpoint commit / main revision at publication | Peeled commit of `m8-reproducible-sessions-2026-10-25`, full SHA in tag annotation |

All 174 source/test inputs in the frozen manifest match the saved checkpoint.
All 379 original package artifacts and 116 indexed qualification evidence files
were checked before publication. The original ancestry, build and binary identities
remain unchanged. Publication adds documentation and text evidence outside the
compiled source manifest. There was no rebuild, package replacement, additional
device test or installation change.

Preserved tested app:
`/Users/fuzz/Developer/LevelScope/build/m8-qualification-2026-10-25/candidate-01/LevelScope Standalone.app`.
The separate High Sierra package and all earlier tested copies are preserved.

The [qualification summary](../m8-qualification-2026-10-25/qualification-summary.json)
and its reports retain their original pre-publication fields. They establish the
qualified state; this publication record and the annotated tag establish the later
saved state. The [handover](../../../../Docs/ProjectIndex/36_M8SessionsQualificationAndHandover.md)
describes implementation, reproduction, platform results and remaining limitations.

Both actual platforms passed session core/recovery, restored audition, native
H.264/ProRes, silent UI, bidirectional exchange and applicable M1–M7/official
regressions. Unchanged engine paths reuse the explicitly scoped historical M7
physical timing/long-measurement evidence; no new physical runs are claimed.
Occasional picture skips remain an accepted development limitation, with sustained
freezes, accumulating drift, interrupted audio, identity and lifecycle failures
still unacceptable. Physical full-eight-channel listening remains separately pending.

Only intended source, tests, manifests, text evidence and documentation are saved.
Finder metadata, raw/generated/downloaded media, build products, connection details
and unrelated work are excluded. All earlier tag objects and targets are retained;
the separate video working copy, applications, plugins, Pro40 and system audio
configuration are untouched.

No subsequent milestone is started. Separate-file stereo/5.1 import remains the
next priority; 7.1.4, stems, dialogue analysis and VST3 parity remain deferred.

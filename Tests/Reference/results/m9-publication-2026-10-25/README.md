# M9 checkpoint publication — 25 October 2026

The annotated tag [`m9-separate-file-stereo-5.1-2026-10-25`](https://github.com/carlozaiana/LevelScope/tree/m9-separate-file-stereo-5.1-2026-10-25)
identifies the accepted M9 checkpoint and main revision at publication. Its
annotation records the full commit SHA, source manifest, original binary hashes
and qualification evidence identity. As in M8, this stable tag reference avoids
a self-referential commit hash inside the commit's own files. Resolve it with:

```sh
git rev-parse 'm9-separate-file-stereo-5.1-2026-10-25^{commit}'
git cat-file -p m9-separate-file-stereo-5.1-2026-10-25
```

[Publication preflight](preflight.json) records all twelve original local/remote
checkpoint tag objects and targets, remote branches, working-copy identities,
unrelated-file hashes and qualification integrity checks. Latest remote main was
still M8 `8e984889c5e73670f370c6c568c8ea1d7902a35c`. The checkpoint was saved on
`codex/logic-module` and integrated into main by fast-forward. No integration fix
or change to the tested implementation was required.

## Frozen candidate → saved checkpoint

| Identity | Value |
| --- | --- |
| Candidate | M9 candidate-04 |
| Source-manifest SHA-256 | `ecf819093af35debedfed260a28a972183cf8484af677f6fbca749f4dbea0d34` |
| Original measurement/build identity | `6cdc517ca6618832430a7dbe678c37d3c459888674a112aa2643c2a25186f329` |
| Original identity.json SHA-256 | `374bfa3b5003d7a67651ccfea71b2a3ceceaf8edc17398573a4862ac82b27898` |
| Tested universal app executable SHA-256 | `dd15614194cb4bb40b86415ff6a93e0ddc0dd35455e5ae03adff4f4748bc6938` |
| Qualification evidence-manifest SHA-256 | `a617b6f8336e9a39676fc12c33b06809195fc61de29b60a3445d9ce09824732a` |
| Checkpoint / main at publication | Peeled checkpoint tag commit; full SHA in annotation |

All 187 source/test manifest entries match the checkpoint. All 442 candidate-04
artifacts and 244 indexed qualification files were verified before publication.
M8 and M9 candidates 01–03 also retain their original package identities and
artifact hashes. The original candidate-04 ancestry/build metadata is preserved;
publication does not relabel or substitute a later rebuild for the tested binaries.
Documentation and text evidence are outside the compiled source manifest.

Preserved tested app:
`/Users/fuzz/Developer/LevelScope/build/m9-qualification-2026-10-25/candidate-04/LevelScope Standalone.app`.
The separate High Sierra packages and earlier qualified copies are untouched.

## Qualification retained

All nine required gate groups in the [acceptance record](../m9-qualification-2026-10-25/acceptance-candidate-04.json)
pass on actual arm64 Sequoia and Intel High Sierra. They cover M1–M8 regressions,
official references, M9 grouped measurement/processing/audition and full-resolution
oracles, cross-platform comparisons, three-hour resources, bidirectional session
exchange, M8-reader rejection, silent UI and the 22 approved muted video cases.
Reuse of unchanged earlier-candidate evidence remains bounded by the recorded
equivalence checks. No new physical-device tests were needed solely for publication.

Occasional picture skips remain an accepted development limitation: five across
the Sequoia short/endurance runs and zero on High Sierra. Historical strict
zero-skip failures remain intact. Audio interruptions, accumulating drift,
sustained freezes, identity errors and memory/lifecycle failures remain blockers.
The muted two-output checks do not establish full physical 5.1/7.1 listening,
which remains a separate qualification scope.

The original qualification records, including their pre-publication status fields,
are preserved byte-for-byte. This publication record and the annotated tag describe
the subsequent saved state. The [handover](../../../../Docs/ProjectIndex/38_M9QualificationAndHandover.md)
records implementation, exact identities, reproduction, results and limitations.

Only intended source, tests, manifests, text evidence and documentation are saved.
Finder metadata, raw/generated/downloaded media, build products, connection details
and unrelated work are excluded. Every earlier tag object and target is preserved.
The video working copy, tested packages, installations, plugins, Pro40 and system
audio configuration remain unchanged.

No next milestone is started. Verified export/reporting remains future work;
stems, dialogue analysis, 7.1.4 and VST3 parity remain deferred.

# M10 checkpoint publication — 26 October 2026

The annotated checkpoint
[`m10-verified-pcm-export-2026-10-26`](https://github.com/carlozaiana/LevelScope/tree/m10-verified-pcm-export-2026-10-26)
identifies the accepted M10 commit and main revision at publication. Its annotation
records the full commit SHA, qualified source manifest, original binary hashes and
evidence identity. As in M9, the stable tag avoids a self-referential commit hash
inside its own files. Resolve the exact checkpoint/main identity with:

```sh
git rev-parse 'm10-verified-pcm-export-2026-10-26^{commit}'
git cat-file -p m10-verified-pcm-export-2026-10-26
```

[Publication preflight](preflight.json) preserves all thirteen original local and
remote tag objects/targets, branch/worktree identities, unrelated-file hashes and
qualification integrity results. Latest remote main remained the M9 checkpoint
`75843dcd017ec7990be2200454a10fc38beb44c4`. M10 is saved on `codex/logic-module`
and main receives the accepted checkpoint by fast-forward. No integration fix or
change to the tested implementation is required.

## Frozen candidate to saved checkpoint

| Identity | Value |
| --- | --- |
| Candidate | M10 candidate-03 |
| Source-manifest SHA-256 | `4ff89967e4927600d0b6defafabb28c372f4068870a0f1a9adb1ee84fe1781d8` |
| Original measurement/build identity | `f36567f5f38a47daef49bc9543252b42b3403e3bdcc2fe986b755bee1576cce9` |
| Original identity.json SHA-256 | `0b0cd42a9f4ee917d6d62f25a903edd7ea354a7d1313d6dfa64341396453cbdd` |
| Tested universal app executable SHA-256 | `92c013999bdecaf904af0888ce68dd59bcd1b3004934b9ddb36ce1a1bf511e15` |
| Tested universal export-test executable SHA-256 | `496a8952947068dbc5a9cad5c9cab7a44daa412fd610ddf557011737b4f6eaab` |
| Qualification evidence-manifest SHA-256 | `adf7a097c8c45fef584db2ef014518147e3569d67702b0e35d4fcc8de714e964` |
| Qualification acceptance SHA-256 | `077f71b3c9ae4f37c3c514ca823fa5e3e39f4738019843740578d2790b82856a` |
| Checkpoint / main at publication | Peeled checkpoint tag commit; full SHA in annotation |

All 208 source inputs and all 470 frozen candidate-03 artifacts match the qualified
manifest. All 350 indexed qualification files were verified. Candidate-01/02 source,
scripts and binaries also retain their recorded hashes. The original candidate
ancestry/build metadata is unchanged: publication does not rewrite or rebuild the
tested package to insert a later commit. Documentation and collected evidence are
outside the compiled source manifest.

Preserved tested app:
`/Users/fuzz/Developer/LevelScope/build/m10-qualification-2026-10-25/candidate-03/LevelScope Standalone.app`.
The separate High Sierra package and all earlier tested copies are untouched.

## Qualification and limits retained

All eleven aggregate gate groups in the
[acceptance receipt](../m10-qualification-2026-10-25/acceptance-candidate-03.json)
passed on actual arm64 Sequoia and Intel High Sierra. Each Mac completed 4,985 M10
core assertions, relevant M1–M9/official-reference regressions, sixteen export cases,
independent PCM24 oracle/decoder checks, ten native/injected video cases and silent
UI checks. Bidirectional exchange freshly measured sixteen opposite-Mac bundles.

Both real three-hour 5.1 RF64 runs produced the same 9,331,200,104-byte audio artifact;
metrics and per-channel peak locations agree exactly. Total times including
prerequisites and extra RF64 audition were 44.96 and 83.53 minutes. Resource,
cancellation and report-size budgets passed. Four short muted video cases and one
30-minute Current State ProRes run per Mac passed with zero skips, audio underruns
or accumulating internal-clock drift. Device configuration remained unchanged.
No rebuild or new physical-device run was needed solely for publication.

Occasional picture skips remain an accepted development limitation; earlier strict
failures remain unchanged. Audio interruptions, drift, sustained freezes, identity
errors and memory/lifecycle failures remain blockers. Muted two-output checks do
not establish full physical 5.1/7.1 listening, which remains separately pending.

Export is bounded to whole-programme PCM24 stereo/conventional 5.1 RIFF/RF64 at
unchanged rate/frame count. Dither occurs once at final quantization; zero gain is
not a bit-identical copy. Current State and the actual reopened artifact must both
pass selected programme loudness/true-peak rules before explicit publication.
This does not certify broadcaster delivery, BWF/BW64 or metadata. Other formats,
resampling, dynamics, stems/dialogue analysis, immersive expansion and VST3 parity
remain outside this checkpoint.

The original qualification records, including their pre-publication status fields,
remain byte-for-byte intact. This publication record and the annotated tag describe
the subsequent saved state. The
[handover](../../../../Docs/ProjectIndex/40_M10QualificationAndHandover.md) records
implementation, interfaces, reproduction, exact identities and limitations.

Only intended source, tests, manifests, text evidence and documentation are saved.
Finder metadata, generated/downloaded media, build products, credentials, connection
details and unrelated work are excluded. Existing checkpoint tags, tested packages,
the video worktree, installations, plugins, Pro40 and audio configuration are
preserved. No next milestone is started by this publication.

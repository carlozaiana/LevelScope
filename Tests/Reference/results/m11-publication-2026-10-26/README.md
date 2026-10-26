# M11 checkpoint publication — 26 October 2026

The annotated checkpoint
[`m11-aligned-stem-composition-2026-10-26`](https://github.com/carlozaiana/LevelScope/tree/m11-aligned-stem-composition-2026-10-26)
identifies the accepted M11 commit and main revision at publication. Its annotation
records the full commit SHA, qualified source manifest, original binary hashes and
evidence identity. As in M10, the stable tag avoids a self-referential commit hash
inside its own files. Resolve the exact checkpoint/main identity with:

```sh
git rev-parse 'm11-aligned-stem-composition-2026-10-26^{commit}'
git cat-file -p m11-aligned-stem-composition-2026-10-26
```

The user authorized publication after the coordinating chat independently reviewed
the completed qualification on 26 October. [Publication preflight](preflight.json)
preserves all fourteen original local and remote tag objects/targets, branch and
worktree identities, unrelated-file hashes and qualification integrity results.
Latest remote main remained the M10 checkpoint
`fcd4c15fe78d6352f5e96faabf3eb52f942235e1`. M11 is saved on `codex/logic-module`
and main receives the accepted checkpoint by fast-forward. No integration fix or
change to the tested implementation is required.

## Frozen candidate to saved checkpoint

| Identity | Value |
| --- | --- |
| Candidate | M11 candidate-02 |
| Source-manifest SHA-256 | `669b350652dbc3fb09d0ff724758208fec21bba40e0f6bc0188b89e6404c2c61` |
| Original measurement/build identity | `70dc29e6fabde64b0be2a0f12790528de7e088956d5817f82b6888611c723d98` |
| Original identity.json SHA-256 | `d9198c3ea9cb98a51f25f1850c4210e8351e0be52a508ba5e39e0b926e9686d6` |
| Tested universal app executable SHA-256 | `deb189e943529dd3e9e9d63a93e1d1cd2003f36752f7ab43f1e099d6c204eae5` |
| Qualification evidence-manifest SHA-256 | `75446c97004d378fdb4d4cb009d9eeb9aba735fe346e590987e490220cbf2832` |
| Qualification acceptance SHA-256 | `933b9756ee09944ee8221bbc05a553152b33934cc995c14a5f782cd034de7e6f` |
| Checkpoint / main at publication | Peeled checkpoint tag commit; full SHA in annotation |

All 222 source inputs, all 370 frozen candidate artifacts and all 146 indexed
qualification files match their recorded hashes. The 55 unchanged dependencies
supporting explicitly scoped M10 evidence reuse also match. The original candidate
ancestry/build metadata is unchanged: publication does not rewrite or rebuild the
tested package to insert a later commit. Documentation and collected text evidence
are outside the compiled source manifest. Historical candidate-01 and prior
milestone evidence are preserved, not relabelled as candidate-02 observations.

Preserved tested applications:

- Sequoia: `/Users/fuzz/Developer/LevelScope/build/m11-qualification-2026-10-26/candidate-02/LevelScope Standalone.app`.
- High Sierra: `/Users/fuzz/LevelScope-M11-qualification-2026-10-26/candidate-02/LevelScope Standalone.app`.

## Qualification and limits retained

The [acceptance receipt](../m11-qualification-2026-10-26/acceptance-candidate-02.json)
binds complete fresh evidence on actual arm64 Sequoia and Intel High Sierra.
Each Mac passed 20 explicitly selected short runs, pinned EBU/ITU checks, 24
composition cases, 72 pre-monitor captures, 24 independently decoded PCM24 exports,
six incoming session-exchange cases and eleven silent UI checks. All compared
captures and short exported artifacts are byte-identical across Macs.

Each Mac completed one consolidated three-hour composition workflow, including
measurement, persistence, resource budgets and a verified 9,331,200,104-byte RF64
export. End-to-end times were 56.54 and 105.64 minutes. The finalized RF64 hash is
identical across platforms. Independent libsndfile decoding matched 98,304 PCM
codes in selected first/middle/final spans; native artifact measurement covers the
whole programme. This is not a claim that the independent decoder read every sample.

Six bounded native/injected H.264 CFR/VFR and ProRes Source/Current cases passed
per Mac. Explicit dependency equivalence supports reuse of identified M10
candidate-03 codec/colour and physical endurance evidence. These remain original
M10 observations, not new M11 physical runs. No physical device was opened for M11.
No unaffected historical suite, long workflow or device test was repeated merely
for publication. The [handover](../../../../Docs/ProjectIndex/43_M11QualificationAndHandover.md)
records exact identities, fresh/reused evidence, commands and resource results.

Occasional picture skips remain an accepted development limitation; historical
strict failures remain unchanged. Audio interruptions, accumulating drift,
sustained freezes, identity errors and memory/lifecycle failures remain blockers.
Full physical multichannel routing/listening remains separately pending.

Composition is bounded to confirmed stereo+stereo, 5.1+mono-to-Centre and 5.1+5.1
at unity component gains and exact sample-zero alignment. There is no resampling,
automatic alignment, hidden attenuation, component DSP or dialogue analysis.
Source remains the original confirmed composition. Versioned Current plans retain
extension points for future component and whole-programme processing. Non-audible
analysis references remain a separate future identity boundary. M10's strict
measured Current/artifact export gates and single final dither stage are preserved.
Immersive expansion, intelligibility claims and VST3 parity remain deferred.

Original qualification records, including their pre-publication status fields,
remain byte-for-byte intact. This separate publication record and the annotated
tag describe the subsequent saved state. Only intended source, tests, manifests,
text evidence and documentation are saved. Finder metadata, raw/generated media,
build products, credentials, connection details and unrelated work are excluded.
Earlier tags, tested packages, the video worktree, installations, plugins, Pro40
and system audio configuration are preserved. No next milestone is started here.

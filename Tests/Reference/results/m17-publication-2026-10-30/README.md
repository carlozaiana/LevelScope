# M17 publication receipt

Candidate-04 published 30 October 2026. **Sequoia tested; High Sierra, video and
three-hour checks pending.** No rebuild or new runtime qualification was required:
the committed implementation is byte-identical to the tested manifest.

- Annotated tag: [`m17-job-lra-presets-sequoia-2026-10-30`](https://github.com/carlozaiana/LevelScope/tree/m17-job-lra-presets-sequoia-2026-10-30).
- Checkpoint commit: [`e422f85ccbf7ae11876e6a31be2409937cb8381c`](https://github.com/carlozaiana/LevelScope/commit/e422f85ccbf7ae11876e6a31be2409937cb8381c).
- Tag object: `d243e4a68503bf719fc01ab2eefb2156d982ab20`.
- Published development branch: `codex/logic-module`.
- Main remains qualified M13 `a9be070ab9eb5cb4f591ff9d991d8999c3a016d0`; no main push or integration occurred.

Source/test/build manifest SHA-256:
`853a5711c162fbbbfa9daff2dbc1acee3fee18ee94445e2ac0ae3407efdeb904`.
Measurement/build identity:
`5d1c4dc69d6df19c0fcf3dd39b53cffdd34495d23b52bed716657cb6e2c507b8`.
Original universal app executable SHA-256:
`b84db4eff7945d5e264aa7f67bab7be153aede1758b6d5818e6d2cbbab680126`.

All **277** committed source/test/build files match candidate-04. The original app
remains at `/Users/fuzz/Developer/LevelScope/build/m17-qualification-2026-10-30/candidate-04/LevelScope Standalone.app`;
all **293** artifact hashes and the package identity file are unchanged. All **128**
indexed textual qualification receipts and their index remain unchanged. Freeze-time
pending/unpublished fields are historical facts; this receipt supplies the later Git
mapping. Base HEAD `06f3f75bcf956787124079d24f8fbd60b46b2fa6` does not identify M17.
Documentation/evidence finalized outside the manifest do not substitute a later build.

Twelve explicitly selected candidate-04 automated commands passed their expected
outcomes, with separate native UI and universal/Intel 10.13 cross-build evidence.
Earlier native UI reuse is attributed to the original executed candidates with
source-equivalence rationale. The [qualification handover](../../../../Docs/ProjectIndex/62_M17QualificationAndHandover.md)
and [original evidence](../m17-qualification-2026-10-30/README.md) retain the complete
scope, historical failures and limitations. Intel cross-build is not High Sierra
runtime qualification. Combined actual High Sierra, cross-Mac, video and three-hour
qualification follows M18, no later than M19, before separately approved integration.
Physical full-7.1 listening remains separately pending.

The branch and new annotated tag were pushed atomically without force. Direct GitHub
ref checks confirm both new identities, unchanged main and all **20 pre-existing tag
objects and peeled targets**, also verified locally. The GitHub connector independently
retrieved the exact checkpoint commit. The registered video worktree location/HEAD
is unchanged. Finder metadata is preserved and excluded; no raw media, packages,
credentials or installed applications/plugins are included in the checkpoint.

See [committed source verification](committed-source-verification.json),
[checkpoint/ref verification](checkpoint-verification.json),
[GitHub commit retrieval](github-commit-verification.json),
[frozen artifact/evidence check](preserved-artifacts.json), and
[before](remote-before.txt)/[after checkpoint](remote-after-checkpoint.txt) references.

The follow-up commit containing this receipt updates publication documentation and
adds [M18's proposal](../../../../Docs/ProjectIndex/63_M18InteractiveTimelineAndProgressiveAnalysis.md)
only. It leaves the tag fixed and changes no source, test, build or original
qualification file. The branch is therefore ahead of the tested checkpoint by this
documentation follow-up. Its final remote verification is retained locally under
`build/m17-publication-2026-10-30/` to avoid a self-referential commit identity.
M18 implementation is not started and requires user approval.

# Publication of confirmed mono LRA proposals

The confirmed mono LRA checkpoint was published to GitHub and fast-forwarded into main on 24 November 2026. All four proposal methods operate on a single channel confirmed as Left, Right or Centre. The annotated tag points to the reviewed handover, including the user's functional acceptance. [Checkpoint](https://github.com/carlozaiana/LevelScope/tree/m20-mono-lra-proposals-2026-11-24) · [Qualification and method coverage](97_MonoLRAProposals.md)

## Git and package identities

| Item | Verified identity |
| --- | --- |
| Previous main | `bb862cdc207f715853702601aa51defa5daa684a` |
| Tested implementation | `81ea642b5b2870a934f3a5913c6736b0f6d41983` |
| Two-platform qualification commit | `86c69bb255ef6c3c6b9b3f0f1d4d42a22fd8146e` |
| Reviewed handover and checkpoint target | `7d529198da9f57bf0aa9e49cf2c992163e7f6055` |
| Annotated tag object | `ab5d679f76d2f29a06ee27657889e161cf00fc5c` |
| Production Source tree | `0cf444a4ffecf3712fe65955aa4eb84276bed958` |
| Frozen candidate 01 manifest SHA-256 | `d30dde12fd1b74befa017ad26c6cdaf2739fa3c03b0b337de9d01e57985d9087` |

An atomic push published the [development branch](https://github.com/carlozaiana/LevelScope/tree/codex/mono-lra-2026-11-24), the new annotated tag and the [main fast-forward](https://github.com/carlozaiana/LevelScope/commit/7d529198da9f57bf0aa9e49cf2c992163e7f6055). Direct remote readback verified the three refs and the peeled tag target. The subsequent publication-receipt commit adds only documentation and evidence. Its final GitHub readback is retained locally in `/Users/fuzz/Developer/LevelScope/build/mono-lra-2026-11-24/publication-2026-11-24/final-verification.json`. [Checkpoint verification](../../Tests/Qualification/MonoLRA/results/publication-2026-11-24/checkpoint-verification.json) · [GitHub ref readback](../../Tests/Qualification/MonoLRA/results/publication-2026-11-24/remote-checkpoint.txt)

The frozen application remains at `/Users/fuzz/Developer/LevelScope/build/mono-lra-2026-11-24/candidate-01/LevelScope Standalone.app`. Publication verified all 471 frozen package entries, the 458 files in the source snapshot against Git, all 29 previous tag objects and targets, and 6,290 protected file hashes across the 16 packages and additional evidence, manuals, sessions and media. The original Tutti WAV, finalized exports and listening notes were preserved. [Package manifest](../../Tests/Qualification/MonoLRA/results/combined-2026-11-24/candidate-01-identity.json)

## Evidence and remaining limits

Five focused groups passed **481 assertions on each actual Mac** (Sequoia 15.7.9 arm64 and High Sierra 10.13.6 Intel). Each method passed independent Current and single-channel PCM24 export verification. All 20 exact comparisons passed, including the selected stereo and conventional 5.1 controls. The separate full Tutti test ran on Sequoia with a 19.2 LU LRA ceiling. All four methods passed there, and the user reported that the workflow worked as expected. [Comparison receipts](../../Tests/Qualification/MonoLRA/results/combined-2026-11-24/comparison.json)

Publication used Git and artifact verification. The existing broad video, full-film and three-hour endurance deferrals continue with no automatic 90-minute substitute. The previously documented completion pause remains open separately.

## Agreed order after publication

The user approved this sequence in the coordinating chat, superseding the earlier proposed order in document 97. The coordinator will begin the work after a separate instruction.

1. Make Propose initialize an absent Current from a valid Source and preserve any existing Current edits. Add **Create sections from applied regions** in the same workflow slice. Creating the sections should preserve the audio and its measured Pass.
2. Add LRA diagnostic lanes.
3. Offer an explicit action to use a fully measured unsuccessful attempt as a manual starting point. Show the missed limits and the criteria used to compare candidates. A changed Current requires independent remeasurement.
4. Use those tools to guide further proposal-rule development.

The request for a quieter release UI and condensed Review proposal remains recorded in document 97.

# M18 publication and combined M14–M18 main integration — 3 November 2026

The user authorised publication and main integration after the advisory chat
reviewed the completed combined qualification. The complete M18 checkpoint and
evidence are published on `codex/ui-timeline-m18`, with the annotated tag
[`m18-interactive-timeline-2026-11-03`](https://github.com/carlozaiana/LevelScope/tree/m18-interactive-timeline-2026-11-03).
Main was fast-forwarded from qualified M13 to this checkpoint. No product change,
merge conflict, rebuild or repeated runtime qualification was needed.

## Exact identities

| Identity | Value |
| --- | --- |
| Tested product/source commit | `3665815d73fbd93c0574cdfee8d7efdd805ab7b6` |
| Qualified checkpoint and initial integrated main | `a9550d9d43ed70fde4c73f90e44081a706d30a67` |
| Annotated tag object | `d147d88fc0173a70abaa624cef686c195685901a` |
| Previous main, preserved in ancestry | `a9be070ab9eb5cb4f591ff9d991d8999c3a016d0` |
| Measurement/build identity | `79ee7b68ded355f1baf36f839c978ad692b8649cb3fe669b9fc9988060737fdb` |
| Candidate-04 identity SHA-256 | `33bc7202fc2629797bbfd3223b3849145a88bebe5c22e9157a1636f53cd69e95` |
| Frozen source-input digest | `f26a3fc790126383ea3f5d6f7fd0a9c64c0770fbdc625b3959d31366c7e5a94f` |
| Original universal app executable SHA-256 | `1b9b8db2aca5d17f7343912110d263f638216152f1bb806368220a4ac79709bf` |

The original frozen app is unchanged:

- Development Mac: `/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m18-qualification-2026-11-03/candidate-04/LevelScope Standalone.app`.
- Golem: `/Users/fuzz/LevelScope-M18-qualification-2026-11-03/candidate-04/LevelScope Standalone.app`.

Its universal executable contains arm64 and x86_64, with Intel minimum macOS
10.13. Production bounded per-channel true peak remains enabled for supported
Apple silicon stereo/5.1 work; Intel production remains serial. The corrected
Intel bounded diagnostic did not change that production policy.

All production source and CMake configuration in the checkpoint/main are identical
to the tested product commit and frozen source inputs. Later commits contain the
documented qualification-tool corrections, supplemental checks, evidence and
handover. They do not identify a later rebuilt app. The frozen manifest's
`qualification_complete: false`, and the completed qualification summary's
`publication_performed: false`, retain their historical status at recording; this
separate receipt records subsequent publication.

## Publication and preservation verification

[Preflight](preflight.json) and [remote references before publication](remote-before.txt)
showed remote main still at M13 and the published logic branch at
`fc4e4b0d54456cb32ab68672ede0a772b663d749`. Both were already ancestors of the
qualified UI branch. All 21 original annotated tag objects and peeled commits
matched locally and remotely. No newer remote work needed reconciliation.

An isolated main checkout was created at
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m18-publication-2026-11-03/integration`.
The normal fast-forward-only merge preserved every intervening M14–M18 commit.
The UI branch, new annotated tag and main were pushed atomically without force.
[Direct remote verification](remote-after-checkpoint.txt),
[checkpoint/preservation checks](checkpoint-verification.json) and an
[independent GitHub connector retrieval](github-checkpoint.json) confirm the
checkpoint, tag object/target and integrated main.

All **395** frozen artifacts and **486** indexed qualification receipts were
rechecked with zero hash mismatches. The Source tree and CMake blob match the
tested product commit. The original 21 tags and all other remote branch heads
remain unchanged. The logic and video checkouts retain their heads, branches,
indexes, status and every snapshotted tracked/non-ignored file: 4,374 and 1,839
files respectively. All 19 snapshotted `.DS_Store` files across the three original
working copies are preserved and excluded from Git. Full per-file inventories
remain in the owned local publication directory; the committed preflight records
their digests and the verification receipt records the results.

The follow-up commit containing this receipt and the project-index/handover
updates is published on both the UI branch and main. It changes only documentation
and publication metadata, leaving the fixed checkpoint tag and tested product
unchanged. Its final exact main/branch references and repeated preservation checks
are retained in `build/m18-publication-2026-11-03/final-verification.json` and
`remote-after-final.txt`, avoiding a self-referential commit hash in this receipt.
The final main identity is also returned in the publication chat.

Large generated PCM/RF64, raw reference media, frozen packages and build products
remain outside source control. Existing packages, original media, installed apps
and plugins are untouched. No app, device, Pro40 service or system audio setting
was changed during publication. No subsequent milestone was started.

## Accepted qualification and follow-ups

The [qualification report](../m18-qualification-2026-11-03/README.md),
[final summary](../m18-qualification-2026-11-03/receipts/qualification-summary.json)
and [current handover](../../../../Docs/ProjectIndex/64_M18UITimelineHandover.md)
give the complete scope. Actual arm64 Sequoia 15.7.9 and Intel High Sierra 10.13.6
passed 28 fresh short/reference groups each, selected injected/native video,
native UI review, four separately authorised muted physical cases per Mac,
bidirectional sessions, synthetic six/eight-channel display-density checks and
one actual three-hour mono/192 kHz RF64 workflow per Mac.

The three-hour workflows took 68.68 / 123.31 minutes, with an additional
9.08 / 16.53-minute complete numerical Source comparison. Both full 6.22 GB
exports are byte-identical; independent audits include samples beyond 4 GiB.
Peak positions, proposal geometry and verdicts agreed across Macs; short processed
captures matched at zero float32 ULP. Combined waveform storage was 100.70 MiB,
within 128 MiB. All eight muted physical cases had zero audio underruns,
skipped/blank frames or timing-bound failures. Original failures and test-driver
corrections remain preserved with their scope and provenance.

Two explicit follow-ups remain in the published candidate:

1. The first waveform for the three-hour mono/192 kHz input took **71.36 seconds
   on Sequoia and 120.66 seconds on High Sierra**, including payload preparation.
   These are observed-analysis timings, not native import-to-waveform timings.
   Faster initial waveform availability remains future performance work.
2. A historical-explanation heading renders its dash incorrectly. All words,
   numerical values and complete structured explanation remain accessible.
   The cosmetic issue was recorded without modifying the frozen app.

Full physical 7.1 listening remains separate. A full video codec/colour matrix
and 30-minute endurance remain conditional on relevant timing/buffering/API
changes or concrete failures. No new full-duration multichannel PCM run is
claimed. Detailed visible-passage sample inspection, device selection/reconnection,
explicit stereo monitoring of 5.1, further DSP/automation/meters and automatic
LRA/M/S-directed correction remain later work. Publication does not expand those
qualification or product claims.

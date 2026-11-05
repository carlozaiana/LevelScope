# Intel true-peak checkpoint publication — 5 November 2026

Published [annotated checkpoint](https://github.com/carlozaiana/LevelScope/tree/intel-multicore-true-peak-2026-11-05)
`intel-multicore-true-peak-2026-11-05` at
`c4eea6f1c4569ec4f31725694ddf34db8340f15e` on `codex/intel-multicore-true-peak`.
`checkpoint-verification.json` and raw `git ls-remote` receipts record the exact
remote annotated object, peeled commit, branch and all preserved tags.

The tested product is `f6fe99716618e94a3a3f5f9c066b5961a06a6075`; the checkpoint
adds qualification evidence/drivers and the next-milestone proposal without
altering product code. Source manifest:
`3733d4be90f08cb82291872057ededc4f3c794534492be05cdd85bc4a7f728fb`.
Universal app hash:
`7833642804d2b4199a31411e91a02c7067f05fb6f42168072e574d19c67f59e0`.
The preserved package is
`/Users/fuzz/Developer/LevelScope/build/intel-peaks-2026-11-05/candidate-01`, with
the unchanged Golem copy under
`/Users/fuzz/LevelScope-IntelPeaks-qualification-2026-11-05/candidate-01`.
Publication rechecked all 361 frozen artifacts, all 138 evidence hashes and the
product source manifest. No rebuild or repeated runtime qualification occurred.

Numerical and scheduling qualification passed on actual Sequoia and High Sierra.
The original absolute 200 ms native UI-gap gate failed and remains recorded.
Matched serial controls show an inherited approximately 0.8-second Golem
completion pause, accepted by the user as a development limitation for this
optimisation. A separate follow-up will profile completion work and its scaling
with programme length. See the [complete evidence](../intel-peaks-2026-11-05/README.md).

Main remains `51e47f300fcc38a386a5a92a7277029175cfab0c`. All 22 pre-existing tag
objects and targets remain unchanged. This publication includes no main merge,
installation, device operation or changes to other worktrees. The documentation
receipt commit following the tag is deliberately separate from the tested
checkpoint. M19 may continue on a fresh branch after publication verification;
its own publication/integration will require review.

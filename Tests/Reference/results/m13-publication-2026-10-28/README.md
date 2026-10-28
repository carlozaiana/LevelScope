# M13 checkpoint publication — 28 October 2026

The user authorized publication of qualified candidate-02 and integration of the
combined M12/M13 work after the completed actual-platform qualification review.
The annotated checkpoint
[`m13-bounded-section-proposals-2026-10-28`](https://github.com/carlozaiana/LevelScope/tree/m13-bounded-section-proposals-2026-10-28)
identifies both the accepted logic-branch commit and main revision at publication.
Its annotation records the full commit, source manifest, original binary and evidence
hashes. As in M11, the immutable tag avoids embedding a self-referential commit hash
inside its own files. Resolve the exact checkpoint/main identity with:

```sh
git rev-parse 'm13-bounded-section-proposals-2026-10-28^{commit}'
git cat-file -p m13-bounded-section-proposals-2026-10-28
```

## Publication and preservation checks

[Preflight](preflight.json) records all sixteen original local/remote tag objects
and peeled commits, the registered worktrees, unrelated Finder-metadata hashes and
the candidate/evidence verification. Fresh fetch showed remote main still at
M11 `cfb1033aa7fb65a2655db276bee696161196bdeb` and remote logic at
`0fc4c5a7f2364f57760228f5c5ecf2f6e0714197`. Main is an ancestor of the logic branch,
which already includes M12 `2fb0e3cfee22a527bacaa09c4ddace23129442fb`.

The complete M13 source is committed on `codex/logic-module`. Main receives that
accepted commit by fast-forward, preserving M11, M12 and all intervening history.
No implementation integration fix or new build is needed. The commit tree must
match all 250 qualified source hashes before pushing. Push only the development
branch, this new annotated tag and the fast-forwarded main; never force-push or
move an existing tag. Direct post-push verification checks both branch heads,
checkpoint tag object/peeled commit, ancestry and every original tag object/target.

## Frozen candidate to saved commit

| Identity | Value |
| --- | --- |
| Candidate | M13 candidate-02, combined M12/M13 |
| Source manifest SHA-256 | `50eba19a2d66096607ec878c8eaa05e610580c7f5e8854d7b375d846a89c5f27` |
| Original measurement/build identity | `5fd55bb1a81c80b4b556715936b6873bccae29a5d04c00b0890a9b43296a20d2` |
| Original identity.json SHA-256 | `266783d8faaabde9fbb166901d1eba8f1c0ff75c4f51a2258a9502817feed839` |
| Original universal standalone executable SHA-256 | `6f843857838c64accc92536f6ccba2c9b0f36113dacca56dc2beef55c96233c5` |
| Qualification evidence-manifest SHA-256 | `4b61b86c6ad53cea151bf2cf52512f977bca4224bdff6904e3025b19a96d1128` |
| Qualification completion receipt SHA-256 | `a314467d7e12bccf5dc1c5b136d711967739e9517ecfd8e393d82a7a1555f4fe` |
| Checkpoint and main at publication | Peeled annotated-tag commit; full SHA in the annotation |

All 250 current source inputs, 473 candidate-02 artifacts and 171 indexed qualification
files match their recorded hashes. The original completion gate was rerun unchanged
against both actual-platform evidence trees and passed identically. Candidate-01
and candidate-02 package equivalence is revalidated: only five qualification-tool
Python files differ; all production sources and native binaries are identical.
The original Sequoia evidence retains its candidate-01 identity where reused.
No qualification report, freeze-time authorization/status field or source/build
metadata is rewritten to pretend the later publication commit was compiled earlier.
Documentation and textual publication records are outside the compiled manifest.

Preserved tested applications:

- Sequoia: `/Users/fuzz/Developer/LevelScope/build/m13-qualification-2026-10-28/candidate-02/LevelScope Standalone.app`.
- High Sierra: `/Users/fuzz/LevelScope-M13-qualification-2026-10-28/candidate-02/LevelScope Standalone.app`.

Candidate-01, published M12 candidate-03 and all earlier packages remain unchanged.
No later rebuild substitutes for these tested universal binaries.

## Accepted scope and qualification

[The handover](../../../../Docs/ProjectIndex/49_M13QualificationAndHandover.md) and
[completion receipt](../m13-qualification-2026-10-28/qualification-complete.json)
record actual arm64 Sequoia 15.7.9 (24G830) and Intel High Sierra 10.13.6 (17G66).
Each Mac passed 22 explicit short/reference checks, M13's 3,671 assertions/40 cases,
M12's 3,304 assertions/32 cases, 17 silent UI checks, six foreign-session exchanges,
80 injected H.264/ProRes receipts plus timing regressions, and four coordinated
muted physical video cases. The pinned EBU/ITU reference suites passed.

Each Mac ran one consolidated three-hour composition/section-processing workflow,
including exact private assessment, explicit Apply, independent Current verification,
resource/session checks, captured audition and real RF64 render/reopen/verification.
Both passed 1,078 native assertions and independent decoder/oracle/hash checks in
117.3 / 202.1 minutes, within every resource budget. The 9,331,200,104-byte RF64
artifact is identical across platforms:
`953b2bcecbf068f5fe153461fb1ea0999939ddb40adf607052b05f4ad2146ade`.
Independent decoding covers bounded sample regions; native artifact measurement and
SHA-256 cover the full output. Retained processed captures match at zero float32 ULP;
maximum compared metric difference is `3.552713678800501e-15`. Comparison tolerances
never widen profile rules. No redundant long, reference or device run was performed
solely for publication because the implementation is unchanged.

M13 preserves already-pass no-op and global-only-first behavior, then permits one
bounded section envelope: <=3 dB additional attenuation, <=8 changed sections,
<=20% affected local duration including transitions, with the 0.5 dB reserve inside
the cap. At most three private passes precede explicit stopped-only Apply and
independent Current remeasurement. Predictions, private candidate, Current and
verified artifact authority stay distinct. M12's shared all-channel processing,
full precision, sessions, audition and M10 strict export gates remain intact.

The one-section cinema evaluation correctly produced a policy-limited failure,
with no Apply/export/playback. Its original remains unchanged; the PCM-identical
RF64 rewrap and user-reported Nuendo/Cubase provenance are documented separately.
No exact external version/settings equivalence or dialogue analysis is claimed.

Occasional picture skips retain the accepted audio/drift/identity/lifecycle conditions;
all eight fresh physical cases recorded zero skips. Historical failures are preserved.
Physical full-7.1 listening remains separately pending. Configurable preferences,
automatic segmentation, marker import, limiter/dynamics, component-relative dialogue
processing, immersive expansion and VST3 parity remain outside this checkpoint.
The inherited report caption says M12; the exact build field identifies M13.

Native report CRLF line endings and captured test-output whitespace are preserved
byte-for-byte because their hashes are part of the evidence. Edited source and
documentation pass formatting checks; raw evidence is not reformatted for Git.

Only relevant source, tests, manifests, textual evidence and documentation enter Git.
Raw media/captures, reference audio, tested packages, generated build products,
credentials, connection details, Finder metadata and unrelated work are excluded.
The video working copy, installations/plugins, Pro40 and system audio are untouched.
No next milestone is started here.

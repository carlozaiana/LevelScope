# M3 qualification — 18 October 2026

**PASS for the approved bounded M3 scope on actual Intel High Sierra and arm64
Sequoia.** Saved on `codex/logic-module` as the annotated checkpoint
[`m3-global-trim-2026-10-18`](https://github.com/carlozaiana/LevelScope/tree/m3-global-trim-2026-10-18)
and fast-forwarded into main from published M2/main
`01a34833ec06a0d237e098faadb23b398581e8f9`. The checkpoint and final main revision
at publication are the tag's peeled commit; its annotation records the full SHA.

See [qualification-summary.json](qualification-summary.json) for counts and
limitations, and [the implementation handover](../../../../Docs/ProjectIndex/22_M3GlobalTrimProposals.md)
for decision policy, identities, build commands and the silent UI checklist.

## Identity and reproduction

- [candidate-04/identity.json](candidate-04/identity.json) identifies every tested
  source file, executable, runner, pinned reference vector and video/UI fixture.
  Its SHA-256 is `14d52b1a03d068cd0875d87f458a9fa80de3b4154fb320afbdd74bcc6229581a`.
- Source digest: `5f66aafc81be73b8184a88f662f1494c9059dfc25cefe43e8e42aac7629f38fd`.
- Build fingerprint: `556349482888d4a86ebea83d1f15403045475a67d787fe422f54abaffdbb88ef`.
- Universal standalone executable SHA-256:
  `4a88b2b3c70e2887633a3202018f804754377657d0d2374b910df1239d8868f5`.

Only the named standalone and five test targets were built. JUCE 7.0.12 is pinned
to `4f43011b96eb0636104cb3e433894cda98243626`; deployment target is 10.13.
The candidate packager and platform runner are in `Tests/`. A new run must use a
fresh candidate/output folder and the manifest-pinned media. Build products and
downloaded/generated media remain excluded from Git. The source manifest, not the
base commit alone, identifies this pre-commit implementation. A later documentation
or checkpoint commit must retain this linkage rather than claiming these binaries
were rebuilt from that later commit.

## Evidence

- [Sequoia](sequoia-04/summary.json) and [High Sierra](high-sierra-04/summary.json)
  automated runs: every suite passed. Each folder includes exact OS/CPU/runtime,
  executable architecture/dependencies, integrity checks, commands, exit codes,
  full M3/M2 check results and both official reference reports.
- [Direct comparison](cross-platform-candidate-04.json): PASS within documented
  cross-platform budgets, with exact M2/M3 outcomes and proposal identities.
  EBU metrics match exactly; maximum ITU difference is 1.7763568394002505e-15 dB TP.
  These comparison budgets never widen profile acceptance rules.
- [Sequoia UI](sequoia-04/ui-smoke.json) and [High Sierra UI](high-sierra-04/ui-smoke.json):
  PASS with devices closed, routes Muted and playback stopped. Actual failure,
  cancellation, no-op, input replacement, video/monitor isolation and shutdown
  behavior are recorded. Screenshots remain local; their hashes are included.
- [Reused evidence review](reused-evidence-review.json) verifies that decoder,
  channel-layout, measurement and profile code underlying inherited M2 format
  checks is unchanged from the published M2 commit.
- [Process headers](ui-processes) independently record native ARM64 and X86-64
  app execution on the required operating systems. The High Sierra run is not
  an Intel slice under Sequoia.
- [CTest](ctest-candidate-04.log): 5/5 passed.
- [Candidate history](candidate-change-review.json), the
  [rejected candidate-02 test](ctest-candidate-02-failure.log) and
  [final source diff](candidate-03-to-04.diff) record the wording corrections and
  obsolete test assertion. Candidate-03 functional UI evidence is reused only
  for unchanged behavior; the final affected paths were rechecked on both OSes.
  Format/unsupported-layout UI coverage is explicitly inherited from M2.
- [Preservation audit](preservation-audit.json), [initial state](before.json) and
  [remote refs](remote-refs.txt): main, both worktrees and all five prior annotated
  checkpoint tags remain unchanged. Existing unrelated Finder metadata was kept.
- [Publication review](publication-review.json): all 100 source/test/manifest
  files, 245 packaged artifacts and 52 pre-publication evidence files verified;
  the retained reports pass the rerun cross-platform comparator. Fast-forward
  integration introduces no implementation changes, so the actual-platform
  qualification remains valid. Earlier preservation records refer to the state
  before publication, when main had not yet advanced to M3.

No outstanding M3 acceptance gap remains. This qualification covers static global
trim and offline measured verification only. Source audition remains original
audio. Dynamics, scenes, LRA steering, freezes, processed audition, persistence,
public export, embedded soundtrack and VST3 measurement parity remain excluded.
No new audible/device or physical picture-sync qualification is claimed. Existing
installations, preserved Golem copies, Pro40 and system audio configuration were
not changed.

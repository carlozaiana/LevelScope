# ProRes qualification evidence, 20 October 2026

See `Docs/ProjectIndex/26_ProResQualification.md` for acceptance status and
`25_ProResAcceptance.md` for the original limits and explicitly user-approved
memory-method revision. Final candidate 04 passes all required checks on both Macs,
including all six fresh 1800-second endurance runs. Earlier failures retain their
original verdicts. The accepted checkpoint is `video-prores-checkpoint-2026-10-21` on the video branch;
this task does not merge the extension into main.
All results here are text records; media, application bundles, downloaded vectors,
private desktop screenshots and credentials are excluded.
Raw JUCE JSON keeps its original CRLF endings, and captured stdout keeps its
original spacing. Do not reformat these hash-bound records to satisfy whitespace lint.

Final candidate 04 source: `6662091f1dfaff948a4a3eb3f000f47e37cfdc21`.
Production last changed at `3920d08f849d2420509e1d789c5b8208e50d79af`.
Implementation: `900edef60baaf35c4ae70fdf626a277e2406e234`.
Base main: `6d304dbecd8fa47f56b0260105aa9d07f96ba1bc`.

## Reading the results

- `identity/` binds frozen packages to source/executable hashes and records all
  checkpoint targets before work. Candidate 03 changes ProRes image backing;
  candidate 04 changes only the qualification harness and its drivers/verifier.
- `fixtures/manifest.json` identifies exact fixture bytes. Historical invalid
  cropped negatives are marked, retained, and excluded from the final required set.
- `runs/` preserves complete original outputs, including failed attempts. Never
  overwrite an old run with a repeat.
- `product-ui/` records agent-observed GUI checks separately from numeric playback.
- `memory-investigation/` preserves failed original-RSS results, diagnostic sources
  and binary identities, raw reusable-page observations and the approved revision.
  Diagnostic experiments never count as qualification.
- `candidate04-native-regression-ui-reuse.json` proves identical production source
  and eight identical executable files for reuse of Sequoia candidate-03 native,
  M1–M5 and UI results. Final playback checks are always fresh candidate-04 runs.
- `m5-cross-platform.json` is the ordinary comparator result: package identity
  mismatches are expected because candidate 02 has a different playback harness.
  `m5-cross-platform-with-binary-reuse.json` verifies exact shared binaries and
  accepts only those two identity mismatches, keeping every numeric check.
- `m5-cross-platform-candidate04-ordinary.json` preserves the final ordinary
  comparator's two package identity/source-manifest differences between Sequoia
  candidate 03 and High Sierra candidate 04. The corresponding `with-reuse` proof
  verifies all shared binaries and retains every numerical/performance check.
- `matrix-evidence-verification-candidate04.json` independently verifies all 32
  final matrix reports; `endurance-evidence-verification-candidate04.json` verifies
  all six full-duration runs, including raw memory-snapshot recomputation.
- `process-highwater-candidate04.json` additionally checks the OS process RSS peak,
  including loading, against 768 MiB in all 38 final reports.
- `preservation-after.json` records publication at checkpoint commit
  `dcc34aabba6b5b3ca4805e2965d21ffebf2e4d0a`: all eight previous tags, exact source,
  all 286 candidate artifacts on each Mac, original movie and original tracked
  checkout verified. Finder metadata differences are recorded and left as found.
  The following documentation-only commit adds this audit without retargeting the tag.
- `audio-configuration-after.json` compares the exposed device/default/rate/source
  snapshots after all test devices closed; both platforms match their baselines.
- `reproduction/` retains session scripts for fixture generation/correction,
  lossless transfer, and the explicit identical-binary regression comparison.
  They contain historical absolute paths. Adapt them only for a fresh artifact root;
  do not rerun against preserved fixtures or candidate folders.

The Sequoia native candidate-01 summary is false because three original cropped
negative fixtures were invalid. Their corrected replacements pass separately.
The candidate-01 playback matrix is false because one Proxy paint skipped a frame.
Both failed records remain visible; neither is used as an all-pass acceptance record.
Candidate 02 Sequoia Proxy endurance and candidate 03 stabilization also failed
the original RSS-growth gate; they retain their failed verdicts. Golem candidate
02 Proxy endurance passed but cannot replace any final candidate-04 endurance run.

## Preserved binary/media locations

Local root:
`/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/prores-20261020`

Golem root:
`/Users/fuzz/LevelScope-ProRes-qualification-2026-10-20`

Local candidates 01–04 and diagnostic-memory-01/02 are preserved. Golem has
candidates 01, 02 and 04. Both roots retain `fixtures` and `evidence`. Golem's actual
movie copy is `51_SyncTest.mov` at its root. The original remains
`/Users/fuzz/Developer/LevelScope/Tests/51_SyncTest.mov`.
Older qualification bundles and installations are unchanged.

Run `Tests/run_prores_native_qualification.py --help` for native checks and
`Tests/run_prores_playback_qualification.py --help` for matrix/endurance arguments.
Always use a new output directory. The playback runner requires an explicitly named
device and keeps every physical route muted. Its identity records source WAV,
movie, executable, OS and architecture. The nonzero controlled-sink isolation
check is independent of physical muting.

Build only named standalone/test targets in a fresh build directory. Do not rebuild
the frozen candidate folders, use an old relocated cache, or build/install VST3.
Internal paint/clock results are not physical optical/acoustic measurements.

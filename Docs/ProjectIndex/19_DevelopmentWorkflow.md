# Development workflow and video checkpoint locations

Recorded 17 October 2026 after publishing the video checkpoint and moving its
registered Git worktree. This is a documentation-only continuation of the tested
checkpoint; no source was rebuilt or requalified for the move.

## Working copies

| Purpose | Source root | Branch |
| --- | --- | --- |
| Original project / logic module | `/Users/fuzz/Developer/LevelScope` | `codex/logic-module` |
| Video implementation | `/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation` | `codex/video-implementation` |

Use the video root explicitly for every subsequent command in the video chat.
Each working copy keeps its normal `Source/`, `Tests/` and `Docs/` directories.
Create project worktrees outside build folders. Reserve each working copy's
`build/` for generated files, compiled applications and temporary test artifacts.
Do not switch another chat's branch, relocate the logic-module checkout, or merge
this milestone into `main`.

Git moved the registered worktree from
`/Users/fuzz/Developer/LevelScope/build/video-implementation/checkout` using
`git worktree move`. The old source directory no longer exists. Both
`git worktree list --porcelain` and the worktree administrative `gitdir` point to
the new location. No source checkout was copied and deleted.

## GitHub checkpoint

- [Video branch](https://github.com/carlozaiana/LevelScope/tree/codex/video-implementation)
- [Annotated video-checkpoint-2026-10-17 tag](https://github.com/carlozaiana/LevelScope/releases/tag/video-checkpoint-2026-10-17)
- Exact checkpoint commit: [`700b9e77e439a8ca9bb7583adf826fac01dbe924`](https://github.com/carlozaiana/LevelScope/commit/700b9e77e439a8ca9bb7583adf826fac01dbe924)
- Annotated tag object: `06f57524bedca1912c0d5c454f13561cec9ce5b4`
- Tested implementation commit: `54e2aa9d2773229dd666bb36df410ff932e9f9bf`

The checkpoint adds documentation and qualification evidence to that implementation;
its implementation source is identical. Both branch and tag were pushed together
without force and verified directly with `git ls-remote`: the branch and peeled tag
resolved to the checkpoint, and the remote tag object matched the local annotated
tag. The branch subsequently advances with this location/workflow documentation;
the checkpoint tag remains fixed. No additional outstanding video source changes
needed committing before the tag.

## Preserved tested artifacts

The historical artifact root stays at
`/Users/fuzz/Developer/LevelScope/build/video-implementation`.
Only its `checkout` source directory moved. These directories remain byte-for-byte
and mode-for-mode as inventoried before relocation:

| Directory under the historical artifact root | Retained contents |
| --- | --- |
| `build/` | Original compiled applications, test executables, generated files and CMake cache |
| `probe/` | Generated fixtures, native probes and historical binary archives |
| `evidence/` | Original raw qualification outputs, including superseded attempts |

A separate verified preservation snapshot is at
`/Users/fuzz/Developer/LevelScope/build/video-implementation/preserved/video-checkpoint-2026-10-17`:

| Directory/file under the snapshot | Contents |
| --- | --- |
| `candidate-06/` | Final Standalone and Video Qualification application bundles, plus VideoTests, MeasurementTests and ScaffoldTests executables |
| `harness-05/` | Earlier Video Qualification bundle extracted from the retained `qualification-05.tar.gz`, plus that unchanged archive |
| `committed-evidence/` | Complete committed video result records as of the checkpoint |
| `raw-evidence/` | Copy of the original raw evidence directory |
| `CMakeCache.historical.txt` | Preserved copy of the old cache, for provenance only |
| `preserved-files-sha256.json` | Snapshot file hashes, sizes, modes and directory inventory |
| `historical-artifacts-before.json` | Hash/mode inventory of original `build`, `probe` and `evidence` contents |
| `video-files-before.json` | Complete source-worktree inventory, including untracked files |
| `original-before.json`, `protected-refs-before.json` | Original checkout and protected-reference audit |
| `github-checkpoint-verification.json` | Verified remote branch, annotated tag object and peeled checkpoint |
| `verification-before.json`, `verification-after-move.json`, `verification-final.json` | Source, evidence, artifact and relocation verification results |

The current build products, rather than the older `probe/candidate-06.tar.gz`, were
copied into `candidate-06/`: that archive predates the final VideoTests expectation
fix and is not a complete final candidate package. It remains intact as a historical
artifact. The snapshot's final products match all entries in the committed
[`candidate-06-hashes.json`](../../Tests/Video/results/video-qualification-2026-10-17/identity/candidate-06-hashes.json).

Key executable SHA-256 values:

| Executable | SHA-256 |
| --- | --- |
| Final Standalone | `5d6253740a7b78266973e4db4a294d64978fb0a7f83032cd6c9809541583bed1` |
| Final Video Qualification | `291bc99b71815056ce9dcb8cf5e20974a48b79bf66ea042b56e33533e926b231` |
| Final VideoTests | `2dd4e27206b3300515d473d162caece49d5af744c6d7d53ecb7f1a9f1c736de0` |
| Final MeasurementTests | `09c89f0cb7e3cb716bd6c77fff8f91bea66bb79d7107caa0ad3780a77be585db` |
| Final ScaffoldTests | `015b069a3698037d548262fc5d0373722415a9142ebe5ce75911635d9fc871bb` |
| Endurance/listening harness 05 | `d7bef4f1f233caf2f75e992ed9ace55c962057d50c5f7dac09e3f73b1ada5c30` |

Golem's tested copies remain untouched under
`~/LevelScope-video-2026-10-17/candidate-06`, `qualification-05` and
`qualification-02`; its M1 qualification bundle remains under
`~/LevelScope-M1-qualification-2026-10-17`. Relocation does not change installed
applications/plugins, Pro40, system audio configuration or any remote copy.
Downloaded reference media, generated media, binaries and credentials stay out of Git.

## Fresh builds after relocation

The retained historical CMake cache contains the old absolute source path. It is
**not reusable after relocation**. Do not run CMake or a build in the historical
artifact directory, patch its cache, or rebuild over the tested products. Use a new
build directory inside the new video working copy, initially `build/video-dev`.
No new build directory or replacement binaries were created during this relocation.

For future configuration, first provide an external JUCE source directory at
`/Users/fuzz/Developer/LevelScope-dependencies/JUCE-7.0.12`, pinned to
`4f43011b96eb0636104cb3e433894cda98243626`. This is a future dependency location,
not a directory provisioned by this task. Verify its revision before configuration.
Keep dependency source checkouts outside `build/` as well; pass the explicit source
override so FetchContent does not create one there. The existing historical JUCE
and baseline source directories elsewhere in the original project's build area
were left untouched to preserve other work and historical artifacts.

Once that prerequisite exists, configure from the new source root:

```sh
cd /Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation
cmake -S . -B build/video-dev \
  -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_OSX_ARCHITECTURES=x86_64;arm64' \
  -DCMAKE_OSX_DEPLOYMENT_TARGET=10.13 \
  -DLEVELSCOPE_BUILD_VIDEO_TESTS=ON \
  -DLEVELSCOPE_BUILD_SCAFFOLD_TESTS=ON \
  -DLEVELSCOPE_BUILD_MEASUREMENT_TESTS=ON \
  -DFETCHCONTENT_SOURCE_DIR_JUCE=/Users/fuzz/Developer/LevelScope-dependencies/JUCE-7.0.12
cmake --build build/video-dev --target LevelScopeStandalone LevelScopeVideoTests LevelScopeVideoQualification LevelScopeMeasurementTests LevelScopeScaffoldTests -j 6
```

Build only the named standalone/test targets: the VST3 target automatically installs
its output. New binaries need their own hashes and validation; the existing tag and
qualification do not certify a rebuild. Write new test outputs to a fresh artifact
directory, leaving the retained evidence and historical recorded paths unchanged.

## Verification and qualification limits

Before relocation, 71 source/test hashes matched both implementation and checkpoint
commits; all 36 committed qualification-evidence hashes and 13 final-candidate
manifest entries matched. After the Git move, all 188 files and the complete
214-entry video inventory matched, including `.DS_Store`, `Source/.DS_Store` and
`Tests/.DS_Store`. Those untracked files remain uncommitted. The later commit changes
only the handover, this workflow document, ProjectIndex and the video test README.

The original checkout remains on `codex/logic-module` at
`b29bc7c183743f35c6271d2d7ec8fb5ba1464591`; its tracked files, index, branch and commit
are unchanged. Its untracked root `.DS_Store` changed during the operation; the audit
records before/after hashes and leaves the current Finder metadata untouched. Its
other untracked files are unchanged. `main`, the measurement checkpoint and earlier
baseline tag objects remain unchanged.

See [the qualification report](18_VideoQualificationResults.md) for the full scope.
Final candidate 06 passed actual Sequoia and High Sierra numeric/native/regression
suites and short VFR/MP4 hardware runs. The 30-minute runs and subjective headphone
observation used preceding harness 05, before the MP4 output-timestamp correction;
the final executable itself has not completed a 30-minute run. Internal timing and
the user's apparently aligned clicks/flashes do not establish an optical/acoustic
latency bound. High Sierra's product GUI observations were narrower than Sequoia's.
Embedded-soundtrack playback/selection, general codecs/HDR, persistence/export and
VST3 parity remain outside this milestone; filesystem-read cancellation has no hard
OS-read timeout. Relocation adds no new playback or qualification claims.

## Location update to paste into other chats

> Video work now uses `/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation`
> on `codex/video-implementation`. The original `/Users/fuzz/Developer/LevelScope`
> remains the logic-module checkout. Fixed tag `video-checkpoint-2026-10-17` points
> to `700b9e77e439a8ca9bb7583adf826fac01dbe924`; the branch also has location/workflow
> documentation. Historical tested artifacts stay under
> `/Users/fuzz/Developer/LevelScope/build/video-implementation`, with a verified
> snapshot in `preserved/video-checkpoint-2026-10-17`. Use fresh builds under the new
> worktree's `build/`; do not reuse the old CMake cache or switch another chat's branch.

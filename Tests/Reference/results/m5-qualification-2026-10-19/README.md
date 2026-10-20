# M5 qualification — explicit programme sections

Candidate-01 was stabilized on native arm64 Sequoia, then the identical universal
bundle was executed on actual Intel High Sierra. **M5 qualification passes on both
platforms**, including automated regressions, both official reference suites,
cross-platform comparison and silent interactive checks. Completed 20 October 2026;
the evidence directory retains its original 19 October start date.

Machine-readable sign-off: [qualification-summary.json](qualification-summary.json).
No acceptance check remains outstanding. The published annotated checkpoint is
[`m5-programme-sections-2026-10-20`](https://github.com/carlozaiana/LevelScope/tree/m5-programme-sections-2026-10-20),
saved on `codex/logic-module` and fast-forwarded into main. Its peeled commit is
both the accepted checkpoint and main revision at publication; the annotation
records the full commit SHA. No implementation changes were needed for integration.

## Tested identity

- Base commit: `9a22749901d9f6d4b9c194c9c22c0a231f2494f6` (qualified M4).
- Working copy/branch: `/Users/fuzz/Developer/LevelScope`, `codex/logic-module`.
- M5 was tested before its checkpoint commit. The base SHA alone does **not**
  identify it. [The original candidate manifest](candidate-01/identity.json)
  records every source hash and retains its historical working-candidate wording.
  All 130 source files match the checkpoint; publication only adds documentation
  and evidence outside the tested source digest.
- Source digest: `4973465a63bb59ec6c07667f96c6ad579856d89e57c213b5037067be8076d0f2`.
- Manifest SHA-256: `3f4af3ab14252ec678272d5a6fac727eafa1adafca9a4c9f2efa4459b64f87c4`.
- Production build fingerprint: `b9d810f4c9755282eecc6301d424e1413e9080b972f0b7266f74fd631231b458`.
- Universal standalone executable SHA-256:
  `c5689afab976fc5d04ab4b47208c022b74b9882ec152a3f345441182b3fc02e9`.
- Release, `x86_64;arm64`, minimum macOS 10.13, JUCE 7.0.12 commit
  `4f43011b96eb0636104cb3e433894cda98243626`.

All executable and pinned reference-vector hashes are in the manifest. Both runners
verify their artifacts before testing. Both reports identify the executing slice;
the UI process samples independently identify the actual architecture and OS.
High Sierra execution is not an Intel slice running under Sequoia.

[Publication preflight](publication-preflight.json) rechecks all 130 implementation
inputs, all 286 candidate artifacts, 26 retained successful runs and their log
hashes, both actual-platform/UI records and the cross-platform comparison. All
pass with no source changes, so the existing native qualification remains valid.
The original binaries, manifests and run records are preserved. This integrity
review is not represented as a new native execution. Raw native JSON keeps its
original CRLF line endings, and the two EBU logs retain their original trailing
spaces so their recorded hashes remain valid; source and documentation whitespace
checks are separate from those byte-preserved outputs.

## Automated results

| Suite | arm64 Sequoia | Intel High Sierra |
| --- | ---: | ---: |
| M5 sections | 2,457 passed | 2,457 passed |
| M4 diagnostics | 369 passed | 369 passed |
| M3 proposals | 302 passed | 302 passed |
| M2 compliance | 128 passed | 128 passed |
| M1 measurement | 184 passed | 184 passed |
| Scaffold | 205 passed | 205 passed |
| Video | 329,054 passed | 329,054 passed |
| EBU / M3 / M4 / M5 runner tests | 8 / 7 / 8 / 16 passed | 8 / 7 / 8 / 16 passed |
| ITU pinned vectors | 39/39 passed | 39/39 passed |
| EBU pinned vectors | 70 passed | 70 passed |

EBU includes 66 minimum-duration and four supplemental checks. Evidence resides in
[sequoia-01b](sequoia-01b/summary.json) and
[high-sierra-01](high-sierra-01/summary.json), with raw logs and structured results.
The first sandboxed Sequoia invocation could not create native temporary test
fixtures. The unchanged candidate passed outside that sandbox; see
[the invocation review](sandbox-run-review.json). This is an environment limitation,
not a changed implementation or a suppressed failing assertion.

[Cross-platform comparison](cross-platform-candidate-01.json) passes. Twenty-four
real measured section cases require exact interval, membership and availability
agreement; all section numerical differences were zero. Allowed finite min/median/max differences are at most 1e-9 LU and energy-share
differences at most 1e-12. Each platform also checks independent full-scan/sorted-value
oracles. Official-reference differences are zero except ITU TP at
1.7763568394002505e-15 dB. Inherited cross-platform budgets are 0.01 for I/M/S/TP
and 0.1 for LRA; these are never added to M2 profile allowances.

## Three-hour performance

Sequoia: macOS 15.7.9 (24G830), Apple M1 Max, native arm64.
High Sierra: macOS 10.13.6 (17G66), Intel Core i7-3770K 3.50 GHz, native x86_64.
The remote runner uses the inspected existing Python 3.8.10 runtime.

The M5 benchmark uses three hours of six-channel, 192 kHz measurement-sized data
and 1,000 sections. It exercises aggregation plus boundary edits, selection,
navigation and painting. This is distinct from PCM decoding/measurement throughput.

| M5 check | Budget | Sequoia | High Sierra |
| --- | ---: | ---: | ---: |
| Aggregation | 2,000 ms | 3.146 ms | 6.303 ms |
| Edit/navigation/paint p95 | 100 ms | 7.262 ms | 20.302 ms |
| Active worker cancellation | 1,000 ms | 0.829 ms | 1.758 ms |
| Additional section storage/workspace | 16 MiB | 1,747,195 bytes | 1,747,195 bytes |

All inherited M4 budgets also pass: analysis 14.168/28.779 ms, worker including
identity 147.793/305.087 ms, viewport p95 0.161/0.232 ms, painting p95 6.015/19.608 ms,
analysis storage 21,030,048 bytes on each, peak RSS growth 101,023,744/72,445,952 bytes
and active cancellation 1.331/1.577 ms (Sequoia/High Sierra).

## Silent UI and interpretation

The interactive checklist covers exact-frame split/move/rename/merge; invalid
boundaries; Source/Current reuse; crossing and EOF-padded windows; global-maximum
labels; short/silent populations; gain invalidation and measured rebuild; profile,
monitor and video independence; cancellation/retry; close/reopen; replacement/clear;
and three-hour navigation through the final second.

Both [Sequoia](sequoia-01b/ui-smoke.json) and
[High Sierra](high-sierra-01/ui-smoke.json) interactive checklists pass. The final
section is `[518352000,518400000)` frames, or `[10799,10800)` seconds. It contains
seven M windows and no complete S window; crossing counts are three M and ten S.
Six zoom steps plus section selection scroll the viewport to `[10631.25,10800]` s.
Numerical inspection at 10799.9 s shows M `[10799.5,10799.9)` and S
`[10796.9,10799.9)`, both -34.992474 LUFS. M remains excluded by the programme gate.
The programme duration stays 10800 s; 15 LRA EOF windows remain separately labelled.
Both isolated apps were closed after checking stopped, muted, device-closed state.

Devices stay closed, routes Muted and playback stopped. The 2,073,600,044-byte
three-hour stereo 48 kHz PCM fixture is actually imported and fully measured for
interactive inspection. Its SHA-256 is
`cf3d1fa528ce1798225b48496a91ced8f1696cbe17dbeb14f548d5356c86076a`.
This complements the six-channel 192 kHz diagnostic-data benchmark.

Sequoia's sequential Build/Cancel attempt completed before Cancel; it is not
claimed as an observed active cancellation. The final native-control tests exercise
pending completion/cancel/retry/suspend/resume on both OSes, and active-worker timing
is measured separately. High Sierra's interactive Build/Cancel visibly cancelled
section publication and an explicit retry succeeded. The underlying valid measured
M4 timeline can remain visible while section summaries are cancelled.

Unchanged media-format/layout and physical video-sync evidence is reused only for
unchanged paths, as documented in [the reuse review](reused-evidence-review.json).
All automated suites and both official sets were nevertheless rerun. No new audible
or physical synchronization claim is made.

## Scope, reproduction and preservation

[The M5 handover](../../../../Docs/ProjectIndex/24_M5ProgrammeSections.md) defines
the exact semantics, architecture, build/test commands and exclusions. Only named
standalone/test targets were built. Reference audio, UI media, build products,
connection sockets and raw desktop screenshots stay outside Git.

M5 adds explicit annotations and measured section descriptions. It does not add
section I/LRA, section compliance, local peak safety, automatic scenes, section
trims, dynamics, protection/freezes, proxy refinement, processed audition,
persistence, public export, embedded-soundtrack playback or VST3 parity. Section
edits cannot change measurements, compliance, M3 drafts/gain or audio/video transport.

Existing installations, plugins, Golem's older qualified copies, Pro40 and system
audio configuration are preserved. All seven pre-existing checkpoint tag objects
and target commits remain unchanged; the video working copy and unrelated files
are preserved. The historical [qualification-time preservation audit](preservation-after.json)
records the pre-publication branch/main base, registered worktrees, original tags
and unrelated untracked files. Publication advances only the logic branch and
main and adds the M5 annotated tag; it does not rewrite that earlier audit.
[Artifact integrity](artifact-integrity-final.json)
verifies all 286 candidate artifacts and retained run-log hashes; all 130 source
files still match the candidate manifest. The final Golem executable hash also
matches. [Capture hashes](ui-capture-manifest.json) identify the local UI evidence.

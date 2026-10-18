# M4 qualification — 18 October 2026

**PASS for the bounded M4 scope, including the three-hour extension, on actual
Intel High Sierra and arm64 Sequoia.** Published on `codex/logic-module` as the
annotated checkpoint
[`m4-measured-diagnostics-2026-10-18`](https://github.com/carlozaiana/LevelScope/tree/m4-measured-diagnostics-2026-10-18),
fast-forwarded into main from qualified M3
`62ab420cb15b3ce5e1c821c963a749a9a56ae4a7`. The checkpoint's peeled commit is the
exact main revision at publication; its annotation records the full commit SHA.
No integration change to the tested implementation was needed.

The [handover](../../../../Docs/ProjectIndex/23_M4MeasuredDiagnostics.md) defines the
read-only timeline, contribution interpretation, identity and availability rules,
build targets and deferred work. [qualification-summary.json](qualification-summary.json)
records the final results and evidence limits.

## Tested identity

- Final [candidate-03 manifest](candidate-03/identity.json): SHA-256
  `6ba208d0a4760d4c77a8e70272c7dd35e8523d6c7334ce9ab97a5c14b554187d`.
- Exact source digest:
  `59bc7a5f17cbcce63c86254e72dbba21ae280d475a19b801bc45c9089d20d9cc`.
- Production build fingerprint:
  `efd7f8afaf7817642ab2b302b2420bcec1edce4dcd481fb6ac6b85546e15bfed`.
- Universal standalone executable SHA-256:
  `aeb88e145945a305fec9be3be962ea888d127baae63f37c0cb8b2e7e3efb96d3`.

The source manifest identifies the implementation as tested before committing;
the base commit is ancestry only. Every manifested source file matches the
published checkpoint. Documentation was finalized afterward; these binaries were
not rebuilt from the publication commit. Candidate manifests and qualification-time
preservation audits retain their original historical descriptions. Raw suite
reports retain their original CRLF line endings and runner whitespace so recorded
hashes remain verifiable; source/document whitespace checks exclude raw evidence.
JUCE 7.0.12 is pinned to `4f43011b96eb0636104cb3e433894cda98243626`; Release,
universal x86_64/arm64, minimum deployment 10.13. Only the standalone and six named
test targets were built. No VST3 target or installation was invoked.

## Results

Every suite passed on both actual operating systems:

| Suite | Checks / cases per platform |
| --- | ---: |
| M4 diagnostics | 369 |
| M3 global trim | 302 |
| M2 measured compliance | 128 |
| Hardware-independent measurement | 184 |
| Standalone scaffold | 205 |
| Video timing/transport | 329,054 |
| EBU / M3 / M4 runners | 8 / 7 / 8 |
| Official ITU pinned vectors | 39 |
| Official EBU v05 | 70 (66 minimum + 4 supplementary) |

- [Sequoia results](sequoia-03b/summary.json): macOS 15.7.9 (24G830), Apple M1 Max,
  native arm64. [High Sierra results](high-sierra-03/summary.json): macOS 10.13.6
  (17G66), Intel Core i7-3770K, native x86_64. Environment records contain CPU,
  runtime, executable architecture/dependencies and artifact hashes; run records
  contain commands, exit codes and log hashes.
- [Cross-platform comparison](cross-platform-candidate-03.json): PASS. EBU metrics
  match exactly; largest ITU difference is 1.7763568394002505e-15 dB TP. Comparison
  budgets remain 0.01 for I/M/S/TP and 0.1 LU for LRA; they do not widen profile rules.
- [CTest](ctest-candidate-03.log): all six registered tests passed on Sequoia.
- [Sequoia UI](sequoia-03b/ui-smoke.json) and [High Sierra UI](high-sierra-03/ui-smoke.json)
  record actual interactive import/roles, measured Source/Current, timeline,
  ranked intervals, padding, peaks, gain invalidation/remeasurement, profile,
  M3 prediction/no-op, video/monitor isolation, silence, replacement/clear and
  close/reopen checks. Devices stayed closed, routes Muted and playback stopped.
  [Native process headers](ui-processes) independently identify both actual OSes
  and app architectures. Raw screenshots remain local, with hashes in the
  [capture manifest](ui-capture-manifest.json).

### Three-hour performance

The six-channel, 192 kHz fixture contains full-resolution measurement populations
for 10,800 seconds. The timings below cover diagnostics and the canonical snapshot
identity, not PCM measurement time. A separate 2,073,600,044-byte three-hour stereo
48 kHz WAV was fully imported/measured in each standalone app and inspected through
its final second. Its SHA-256 is
`cf3d1fa528ce1798225b48496a91ced8f1696cbe17dbeb14f548d5356c86076a`.

| Operation / resource | Limit | Sequoia | High Sierra |
| --- | ---: | ---: | ---: |
| Core analysis | 2,000 ms | 11.52 ms | 29.13 ms |
| Complete worker + identity | 2,000 ms | 147.06 ms | 306.22 ms |
| Viewport preparation p95 | 50 ms | 0.145 ms | 0.232 ms |
| Offscreen paint p95 | 100 ms | 5.79 ms | 19.90 ms |
| Additional analysis storage | 64 MiB | 20.06 MiB | 20.06 MiB |
| Additional process peak RSS | 128 MiB | 92.48 MiB | 69.59 MiB |
| Active cancellation | 1,000 ms | 1.34 ms | 1.58 ms |

Full-resolution rankings/calculations and visible extrema/gate membership are
verified independently of display reduction. These budgets qualify diagnostics
through three hours; they are not a three-hour limit imposed on the measurement API.

### Cancellation evidence boundary

Manual and automated sequential clicks on Sequoia reached Cancel after the roughly
150 ms three-hour rebuild completed; the user confirmed the same outcome. This is
**not an observed interactive active-cancellation pass**. The final native test
instantiates the actual diagnostics panel, invokes its Build/Cancel and window
close/suspend/reopen bindings at a deterministic pending-completion boundary,
and checks cleared plots, no stale publication and unchanged session state. It
passes on both actual OSes. A separate test cancels an active three-hour analysis
within the deadline. No production sleeps, artificial UI delays or test hooks were
added. This combined evidence covers cancellation behavior without mislabelling
the manual timing attempt.

## Review, fixes and preservation

- [Candidate history](candidate-change-review.json): candidate-02 fixes a false
  edge marker for an offscreen selection and misleading short-programme wording
  for long silence. Affected UI paths were rechecked on both platforms. Candidate-03
  adds native control tests only; its standalone app is byte-identical to candidate-02.
  All final automated and official-reference suites ran on both actual platforms.
- The [initial rejected fixture](rejected-diagnostic-fixture.json) used a different
  floating-point expression from the qualified loudness function. The synthetic
  fixture was corrected to use that function; no epsilon or acceptance relaxation
  was introduced. The original result remains recorded.
- M1 gate helpers were extracted with the original expression/order/comparisons;
  both official reference suites were rerun because shared measurement code changed.
- [Reused evidence review](reused-evidence-review.json): inherited format/layout
  silent UI evidence is reused for unchanged M2/M3 import/audio paths. M4 adds no
  new audible, physical device or picture-synchronization claim.
- Initial candidate-01 transfer emitted macOS tar metadata warnings; all packaged
  content hashes subsequently verified. Later overlays used portable USTAR.
  Final Sequoia runner launch initially used an incorrect working directory;
  no suite ran under that launch. Fresh `sequoia-03b` holds the successful run.
- [Preservation audit](preservation-audit.json): main, the registered logic and
  video working-copy locations, and all six pre-existing annotated tags remain
  unchanged. Unrelated Finder metadata, installed apps/plugins, prior Golem copies,
  Pro40 and system audio configuration were left untouched.
- [Final integrity check](final-integrity-check.json): all 118 manifested source
  files, candidate binaries, 22 run logs and local UI capture hashes match the
  recorded evidence. Both isolated qualification apps have exited. The retained
  reports were compared again successfully after documentation updates.

No unresolved functional M4 acceptance gap remains. Scene processing, dynamics,
LRA steering, T0–T3 placement, freezes, processed
audition, persistence, public export, embedded soundtrack and VST3 parity remain
deferred. Contributions describe overlapping measured populations; they do not
predict the causal effect of edits or establish additional delivery rules.

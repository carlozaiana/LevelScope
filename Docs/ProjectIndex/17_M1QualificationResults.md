# M1 measurement qualification — 17 October 2026

**PASS for the defined M1 scope on actual Intel High Sierra and arm64 Sequoia.**
Both official reference suites, hardware-independent regressions, runner checks and
the documented silent interactive UI checklist passed. No measurement implementation
change was needed during final qualification. This is not delivery-profile, live-meter,
physical-device, plugin or video qualification.

## Checkpoint and reproducibility

The tested implementation is commit
`b3fe9dece6b9d4cd2d8f6a3aebca69b8d0534e9c`, tree
`d9a0846a06004c8eddfe26f940da50c7d52b6a31` on `codex/logic-module`.
The dated annotated checkpoint is `measurement-milestone-2026-10-17`; it includes
the final evidence and handovers without changing that tested implementation.
Resolve its exact commit with `git rev-parse measurement-milestone-2026-10-17^{commit}`.
`main` and the two baseline/scaffold tags retain their original objects and commits.

Evidence lives in
[`Tests/Reference/results/m1-qualification-2026-10-17`](../../Tests/Reference/results/m1-qualification-2026-10-17).
Each platform has the environment report, suite JSON, regression/runner logs,
explicit-architecture scaffold log, app-process architecture header and observed UI
checklist. `identity.json` pins the committed source files, build configuration,
compiler and every delivered artifact. `evidence-sha256.json` covers the saved evidence.

Build used Xcode 16.4 (16F6), macOS SDK 15.5, CMake 4.4.3, Apple Clang 17,
Release C++17 and pinned JUCE 7.0.12
(`4f43011b96eb0636104cb3e433894cda98243626`). The universal executable deployment
minimums remain x86_64 macOS 10.13 and arm64 macOS 11.0. Measurement build fingerprint:
`e6663d640d43edcaf8c5bf53b8cc6e4b7c2111463d603fa232f435ab0a53822e`.
Engine: `bs1770-5-ebu2023-m1.1`.

| Runtime | Exact OS | CPU / model | Execution | Python |
| --- | --- | --- | --- | --- |
| Sequoia | 15.7.9 (24G830) | Apple M1 Max / Mac13,1 | arm64 | 3.9.6 |
| High Sierra | 10.13.6 (17G66) | Intel Core i7-3770K 3.50 GHz / iMac14,2 | x86_64 | existing 3.8.10 |

High Sierra dependencies were inspected before transfer/execution. Existing Python 3.8.10
and system frameworks/libc++ sufficed. No runtime installation was performed. The app
process headers verify ARM64 and X86-64 respectively; official reports record the
measurement executable's actual architecture. The scaffold was also explicitly invoked
with `arch -arm64` and `arch -x86_64` on the corresponding physical machines.

The same universal artifacts ran on both systems:

| Artifact | SHA-256 |
| --- | --- |
| App executable | `2ae5ad2702d3b1a0e832c0866440ce8fba51726b5132522f28ffc4f5c4481f10` |
| Measurement tests | `ada9a38eef62cd87a978cc7f1b37f7071c99efa170c3d5fc7dcb02b26b73d6b5` |
| Scaffold tests | `747565211d9b4c8a7ecf06712e2909060798a64e0f71cd55021ff0b2e4ea0822` |

Use the build commands in [13](13_StandaloneScaffoldHandover.md) and packaging/run
commands in [Tests/Reference/README.md](../../Tests/Reference/README.md). Build only the
named standalone/test targets. Downloaded vectors and build products remain outside Git.
On High Sierra, the test copy is under
`~/LevelScope-M1-qualification-2026-10-17/m1-qualification-2026-10-17`; the existing app
and plugin installations were preserved. The saved identity is the authority for this
run; a later rebuild produces new hashes and requires new qualification evidence.

## Results and comparisons

| Check | Sequoia arm64 | High Sierra x86_64 |
| --- | --- | --- |
| Hardware-independent measurement | 184 passed | 184 passed |
| Hardware-independent scaffold | 205 passed | 205 passed |
| EBU runner failure-detection tests | 8 passed | 8 passed |
| Official ITU BS.2217 files | 39/39 passed | 39/39 passed |
| Official EBU v05 minimum-requirements files | 66/66 passed | 66/66 passed |
| EBU minimum-requirements coverage | 72 rules / 1,052 readings | 72 rules / 1,052 readings |
| Supplementary EBU sanity files | 4/4 passed | 4/4 passed |
| Silent interactive UI checklist | passed | passed |

All vectors matched the pinned manifests. EBU loudness acceptance uses the applicable
±0.1 LU rules; LRA ±1 LU; true-peak cases 15–23 +0.2/−0.4 dB. The supplementary
alignment-file ±0.1 LU sanity limit is a project check, not an added official criterion.
Official EBU vectors are 48 kHz; generated multirate tests are separate evidence.

Cross-platform budgets are ±0.01 LU/dB for I/M/S/peaks and ±0.1 LU for LRA. Both suites
had zero comparison errors. EBU whole-file metrics matched exactly. ITU's largest
whole-file difference was `1.7763568394002505e-15 dB` in true peak; the other metrics
matched exactly. EBU curve timestamps, curve values, successive-max checkpoints and
criterion min/max readings also passed the runner's comparison. See
`cross-platform-comparison.json` for the machine-readable summary.

## Observed silent UI qualification

The agent operated the local native UI and Golem through the user's authenticated
Screen Sharing session. These are observed interactive checks, separate from headless
regression tests. Both platforms covered:

- Import and explicit file-role confirmation; measurement unavailable before confirmation.
- Source measurement, separate Current State initialization, and measured parity at 0 dB.
- −6 dB trim, immediate invalidation, revision change, cleared stale metrics, and successful
  remeasurement with the expected loudness/peak shift and unchanged LRA.
- Successful source replacement clearing Current State and requiring role confirmation.
- WAV, AIFF and FLAC import/source measurement; 12-channel import with explicit M1
  measurement refusal, including inspection of the final channel rows.
- Cancellation while a ten-minute source was actively measuring, with result discarded.
- Source clear; close while measurement was in progress; process exit and clean relaunch.
- Export readiness remaining unavailable without verification/approval.

The device stayed closed, all file routes stayed muted, playback was never started,
and monitor gain remained −30 dB. No physical audio or synchronization claim follows.
UI values are displayed to one decimal place; exact parity comes from numerical tests.
Per-platform `ui-smoke.json` records the actual fixtures and displayed observations.

## Investigated interruptions

The initial Sequoia native picker stalled during overlapping manual/automated input.
After restarting the candidate and using exclusive interaction, repeated EBU, AIFF,
FLAC, 12-channel and long-WAV imports succeeded. No reproducible application defect
remained, and the picker implementation was not changed.

An initial sandboxed regression attempt could not create a generated fixture in macOS
temporary storage. The complete portable run outside that sandbox passed; the failed
attempt is not counted as qualification. Sequoia tar extended-attribute metadata also
caused warnings from High Sierra's older tar. The transfer archive hash matched
`c22014b130fc2c86812187b393a7e00d31ecb71b2d7ba58ace10a54216d3c199` and the complete
delivered artifact inventory passed SHA-256 verification before any tests ran. No payload
was missing or modified. Future transfers can use Python's tarfile writer to omit macOS
extended attributes.

## Limits and handover

M1 is full-program offline Source/Current State measurement for mono, stereo and
conventional 5.0/5.1 at 44.1/48/96/192 kHz, with static gain only. LFE is excluded from
loudness and included in peak measurement. Short-program LRA remains marked unstable.
Reference cases 11/14 validate offline kernel maxima, not live display scheduling.

No authoritative delivery profiles/compliance verdicts, public exporter, proposals,
session persistence, processed audition, resampling, pause/seek, synchronized picture,
physical-device qualification or VST3 measurement parity are included. No plugin target
was built, updated or installed. Plugin parameters, saved-state compatibility and DSP
behavior remain unchanged. Pro40 repositories, drivers, services and system audio
configuration were untouched.

There are no outstanding checks within the documented M1 qualification gate. Continue
with [the video handover](15_VideoImplementationHandover.md); keep
[VST3 measurement parity](16_Vst3MeasurementParityMilestone.md) as a separate milestone.

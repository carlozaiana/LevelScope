# M7 completion decision — 25 October 2026

The user lifted the M7 development freeze and authorized integration of the
published diagnostic tooling, followed by the remaining qualification. Work
continues in the existing logic working copy on `codex/logic-module`. Latest
published main incorporated: `5e07bbd29930db34236d6dfcede7fea14983f90a`.
Candidate-02 and all earlier packages and evidence remain immutable. The combined
candidate is a new build, not a relabelled candidate-02 or a video playback fix.

## Revised development acceptance

Occasional skipped pictures are accepted during development provided audio stays
continuous and the picture follows the authoritative audio clock without
accumulating synchronization drift. Skipped-picture counts, event audio times,
paint intervals and frame selections remain recorded. Sustained freezes, drift,
audio interruptions, identity errors and memory/lifecycle failures remain failures.
No repeat-until-zero programme or root-cause fix is required for M7.

The ordinary harness retains its strict zero-interior-skip `passed` field and
exit status. A separate `m7-development-picture-skips-2026-10-25` evaluator records
the revised decision. Historical strict failures remain failures in their original
reports. It independently recomputes the existing memory limits, verifies media,
build, architecture, processing and eight-channel identities, and requires the
existing zero underrun, blank-picture, frame-bound and CPU-paint-deadline failures.
Coverage counts both painted and explicitly recorded skipped frames. An unexplained
raw failure, incomplete skip accounting or non-finite evidence fails closed.

For unattended runs, a maximum interval between observed paints above 250 ms is a
conservative stop/review guard against sustained stalls. This is not an expanded
audio/video synchronization allowance: the existing per-frame selection bounds
remain unchanged. First/last mean skew must remain within one frame on CFR cases.
Counts are reported for review; there is no new product guarantee for an allowable
skip rate. These are internal clock/UI observations, not optical/acoustic proof.

## Fixed qualification programme

After Sequoia stabilization, freeze one universal arm64/x86_64 candidate with
macOS 10.13 deployment support. Run the same package on actual High Sierra.
Required fresh evidence on each Mac:

- M1–M7 short native suites, sample captures and runner guards; pinned EBU/ITU
  suites; native H.264 and ProRes regressions; optional diagnostic buffer/runner,
  device-closed export and overwrite guards. Check trace symbols absent from the
  normal app/harness in both slices.
- Affected silent UI import/roles, Source/Current measurement, diagnostics,
  partial-monitoring indication, held reference picture and stopped/closed
  lifecycle checks. Earlier unchanged UI evidence is identified explicitly.
- Sixteen Current State ProRes cases (15 generated, one reference movie), four
  H.264 Source/Current CFR/VFR cases, and three 30-minute Current ProRes runs.
  Mac Studio Speakers at 48 kHz; Golem Built-in Output at 44.1 kHz; every route
  muted. The user reconfirmed availability and these inaudible device runs.

The prior real three-hour measurement/seek runs and Source ProRes matrices may be
reused only with a recorded implementation/toolchain equivalence review. Their
original source and binary identities remain attached. Fresh regressions must
confirm the new build's behavior; reused evidence does not claim old binaries
were rebuilt or rerun. No tracing runs are needed to search for a zero-skip result.

Full physical eight-channel routing/listening remains **pending**, separately
from required offline/automated qualification. Two-output all-muted runs never
establish spatial listening. M7 remains unqualified until required evidence closes.

## Roadmap boundary

After M7, defer 7.1.4 and further immersive expansion while preserving twelve-
channel interfaces. Prioritize the stereo/5.1 broadcast workflow: measurement and
diagnostics, explicit processing/remeasurement, verified audition/export, useful
reports and saved/reopened sessions. Separate-file import and large-file containers
follow demonstrated workflow needs. These priorities add no M7 implementation.

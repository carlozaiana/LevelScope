# M18 progressive per-channel PCM waveform slice

Short Sequoia/arm64 development evidence. The prior first-slice evidence remains
unchanged in `../m18-development-2026-10-31`. Main and other worktrees are unchanged.
This is waveform-slice validation, not full M18 qualification.

The new runtime overview is built from actual verified programme samples during
analysis, after the shared applied processor for Current. Every channel retains
its own minimum, maximum and earliest corresponding sample positions. The fixed
multiresolution storage coarsens its initial 64-frame resolution to fit a shared
128 MiB pool: at most 60 MiB per overview, plus an 8 MiB reserve for display buffers,
metadata and allocation overhead. Old worker/UI references remain charged until
released. Identical initial Source/Current measurements share one allocation.
Publications expose immutable completed prefixes at most every 200 ms, with an
immediate final update. Display preparation is capped at 4,096 columns. The
inspector reports the represented half-open interval, duration, extrema and exact
frames. Unread regions stay pending. M/S use their own measured windows.

The main timeline receives progressive worker data without audio callbacks or a
full PCM buffer. Runtime, input, roles, processing, job and revision checks reject
stale previews. Build waveform uses verified input and the exact applied processor
without replacing accepted measurements. No overview is persisted and no
sample-level detail read was implemented.

## Commands and evidence

`run_short.py <new-output-directory>` runs the seven explicit audio-only groups
listed in the plan impact manifest. It records exact commands, exit codes, timings,
Source/test/CMake hashes, measurement build identity, binary hashes, universal
architectures and Intel minimum versions. Every group must pass. The new waveform
suite writes partial/complete software-rendered component captures for visual QA.
JUCE uses its macOS cache directory for verified-copy fixtures, so execution needs
the existing sandbox escalation. No physical device or video is opened.

Build directory: `build/m18-dev-01`, Release universal arm64/x86_64, Intel minimum
10.13, JUCE 7.0.12 commit `4f43011b96eb0636104cb3e433894cda98243626` from the existing
local dependency checkout. Configuration adds `LEVELSCOPE_BUILD_WAVEFORM_TESTS=ON`
and `LEVELSCOPE_BUILD_SECTION_PROCESSING_TESTS=ON` to the prior first-slice options.
Only the standalone, waveform core/integration, measurement, diagnostics, short
session, section-processing smoke and export targets are built. No default/all or
plugin build/install was invoked.

`iteration-01` retains configuration/build logs, the first 977-check development
run, and standalone Clang sanitizer runs of `WaveformOverviewTests.cpp` plus
`WaveformOverview.cpp`: `-std=c++17 -O1 -g -fsanitize=thread`, and separately
`-fsanitize=address,undefined -fno-omit-frame-pointer`. Both passed all 240,560
assertions with no sanitizer findings. The waveform core was unchanged afterwards.

`candidate-01` passed all seven groups. Review then clarified one inspector label
so a missing-waveform rebuild does not suggest its valid programme measurement is
pending. That candidate is superseded by the subsequent final build and run.

`reuse.json` identifies production changes relative to the first slice. Complete
target/preset, compliance, global/section-proposal and preferences evidence is
reused from that slice because those implementations did not change. Qualified
reference evidence is reused for unchanged measurement arithmetic; fresh exact
observer on/off comparisons guard every metric, gate, population, measured window,
peak value and position through the new service boundary.

The large-allocation checks exercise three-hour-sized twelve-channel overview
storage, not three hours of PCM decoding. Full-duration RSS/timing, actual High
Sierra execution, cross-Mac exchange, physical/video tests and comprehensive native
interaction remain pending in the agreed combined qualification block. The new
captures validate component rendering, not an end-to-end native interaction test.

## Final waveform-slice development result

`candidate-02` is the final development evidence. All seven groups passed against
measurement build `d6a5d2e5bcb87f26560dadbe17deb2fdd1fca7f00de08e533df2c8ee474fb80f`:

| Group | Checks |
| --- | ---: |
| Waveform core | 240,560 |
| Waveform service/cache/UI | 1,173 |
| Measurement | 185 |
| Diagnostics | 369 |
| Short sessions | 1,071 |
| Section processing smoke | 334 |
| Export | 4,872 |

The aggregate short-fixture analysis ratio with waveform collection versus without
was 1.003558 (0.36% overhead). The largest of 100 sampled viewport
preparation plus software-paint operations was 4.479 ms at 1,200 × 650.
These are short development timings, not full-duration or native-display qualification.
All eight final release executables (standalone plus seven test targets) contain
arm64/x86_64 slices and declare Intel macOS 10.13 minimum. Actual Intel execution
remains pending. The final source, test and CMake hashes were checked again after
the run, and every report carries this build identity. Dependency warnings came
from the existing JUCE native text-input implementation.

The partial and complete component captures in `candidate-02` were visually
inspected. Final build logs remain in `iteration-01/build-final.log` and
`build-final-check.log`; the final check adds portable thread linkage for the pure
C++ test target, without changing the measurement build identity.

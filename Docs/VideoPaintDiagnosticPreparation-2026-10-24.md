# Bounded video paint diagnostic preparation — 24 October 2026

Status: candidate-03 is prepared for a fixed comparison on Sequoia and actual
Intel High Sierra. Physical playback has not started and requires fresh user
coordination. This is an instrumented diagnostic, not a video correction or an
accepted milestone. M7 remains frozen, unqualified and unpublished.

## Base and scope

The video worktree at
`/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation` was fast-forwarded
on `codex/video-implementation` to qualified published M6,
`b57a1ed23cbe7d74afe98ed16ec8fd615db9d171`. Diagnostic changes are local and
uncommitted. Main, the original `codex/logic-module` working copy, its pending M7
work and all checkpoint tags remain unchanged. The accepted workflow rename is
retained. No other checkout, branch, installed application, plugin, driver or
audio setting was changed.

The only production-source edit is a read-only test-friend declaration in
`Source/Standalone/StandaloneVideo.h`. It changes neither Player layout nor its
selection, decoding or buffering implementation. Test hooks in
`Tests/VideoPlaybackQualification.cpp` are enabled only for the new
`LevelScopeVideoPaintDiagnostic` target. The original frozen M6 executable is the
control. The ordinary qualification gates and frame/audio-clock calculations
remain intact. Diagnostic reports explicitly set `qualification_eligible=false`.

All builds used named targets in a fresh directory:
`build/paint-diagnostic-20261024-01/cmake`. Only
`LevelScopeVideoPaintDiagnostic`, `LevelScopeVideoPaintTraceTests`,
`LevelScopeVideoTests`, `LevelScopeMeasurementTests` and `LevelScopeAuditionTests`
were built. No default/all target, VST3 target or installation was run. JUCE
7.0.12 remains pinned at `4f43011b96eb0636104cb3e433894cda98243626`; its shared
checkout was not edited. Binaries contain arm64 and x86_64 slices, with the Intel
build targeting macOS 10.13.

## Exact package

- Local: `/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/paint-diagnostic-20261024-01/candidate-03`
- Golem: `/Users/fuzz/LevelScope-VideoPaintDiagnostic-2026-10-24/candidate-03`
- Source digest: `0760e080d2e5daa8c89612029236308e134237a142c879003b53c7da1d3727cf`
- `identity.json` SHA-256: `35081f6e31cd02f202162525699ecd87912f4cabbc7177a2919f4c53497cad1c`
- Diagnostic executable SHA-256: `3fd0fe58bbba684c06ac55acaf1618c318812396e3a121434e54a7dcfde1f93b`
- Frozen M6 control executable SHA-256: `c8efc1a6d0b760206c8f99105e553542a4d32fdb9722dc3e19b2172f60bb36e9`

Both locations contain the same 143 manifest artifacts. Copies of the original
M6 app, the unchanged actual `51_SyncTest.mov`, and both selected-clock WAV files
are inside the package. The originals are untouched. Source snapshots and build
metadata identify the uncommitted diagnostic independently of its base commit.
The package directories are read-only. Results go outside them. Media, binaries,
archives and build caches remain under build/test-package locations and out of Git.

The actual movie hash is
`0870d5ab93d1492f23524ed1fff9319af1a6a9cf71882543f2190929afeb8e24`.
All fixture hashes and source inputs are in the manifest. The actual movie is
1080p/25 standard ProRes 422; the proposed comparison uses it for 100 seconds.

Candidate-01 and candidate-02 are preserved preparation iterations, not
qualification runs. Candidate-01's preview found a Python serialization error
in macOS profiler metadata. The runner was corrected to serialize timestamps and
to compare device settings independently of profiler execution timestamps.
Candidate-02 was superseded before its checks. Candidate-03 includes the final
runner fixes and six passing runner tests. All three have the same diagnostic
executable hash; no playback implementation correction was made in response to
the preview failure. The failed preview remains recorded.

## Instrumentation and its limits

The diagnostic records timer entry/exit, successful or missing picture selection,
repaint requests, paint entry/full callback exit, JUCE main-thread vblank delivery,
an independent native display-link observer, and native window-state transitions.
Selection records retain exact rational source time, transport generation, view
epoch, desired/selected frame and a post-selection queue snapshot. An additional
short mutex acquisition samples Player's queue only on the UI thread. No trace
code is added to the audio callback or decoder worker.

Separate preallocated single-producer buffers hold 131072 UI events and 32768
display events, using 17039360 bytes (16.25 MiB). Recording performs no heap
allocation, file I/O, JSON serialization or memory-statistics queries. Queue
sampling remains bounded by the existing eight-frame limit. Overflow preserves
the earlier history, counts losses and invalidates the diagnostic. The extra
memory is included in the unchanged RSS ceiling; it is not subtracted from it.
Export occurs after audio disconnects, uses exclusive file creation and refuses
to overwrite an existing trace.

`ui_events` and `display_events` use the same monotonic host nanosecond clock.
The trace describes its column order, event kinds and window flags. Exact source
time is a numerator/denominator pair. Window flags identify visible, unoccluded,
active, key, minimized and live-resize states. The trace also records display ID,
nominal rate and pixel dimensions. The ordinary harness report remains alongside
the trace, with its original frame-skip and clock gates unchanged.

The independent CVDisplayLink is an additional observer, not a hook inside JUCE's
own display-link callback. It can distinguish continued refresh observations
from delayed main-thread/paint activity, but cannot by itself prove the exact
latency or loss count inside JUCE's callback-to-AsyncUpdater handoff. No shared
JUCE dependency was patched. Native invalidation/draw entry and compositor
presentation are not instrumented in this first bounded stage. Full paint-exit
timing includes the harness status text, but not the entire native drawing stack.
These observations do not measure actual scanout or physical optical/acoustic
synchronization. Traces include preroll and playback; maxima over the entire
trace must not be confused with the existing post-preroll acceptance metrics.

Instrumentation adds a display-link observer, bounded tracing and an extra UI
lock, so its overhead may perturb timing. The fixed uninstrumented M6 comparison
is necessary. A passing diagnostic run is neither a fix nor proof that the
intermittent problem is gone.

## Device-closed readiness checks

On both actual Macs, the trace-buffer checks, diagnostic startup/export,
exclusive-output guard, video timing, offline measurement, injected M6 audition,
six Python runner guard tests and preview-only preflight passed. No physical
playback ran. Local execution was arm64 on Sequoia; Golem executed the same
universal binaries as x86_64 on High Sierra 10.13.6. All 143 package artifacts
verified on both hosts, including after Golem's checks.

Local evidence:
`build/paint-diagnostic-20261024-01/local-closed-regressions.log` and
`build/paint-diagnostic-20261024-01/sequoia-preparation-03/`.
The three local regressions ran before the Python-only runner fixes; the final
package contains the same tested C++ binaries. Golem evidence:
`/Users/fuzz/LevelScope-VideoPaintDiagnostic-2026-10-24/high-sierra-preparation-03/`.
Relevant records and the manifest are copied into
`Tests/Video/results/paint-diagnostic-2026-10-24/` for handover.

The startup/export check exits before creating a window, Player or audition
object. It deliberately writes a trace with no display callbacks and
`valid=false`; this proves startup/export only and is not a playback result.
The duplicate-export attempt is expected to exit 1 while preserving the original
file. No physical audio device was opened by these readiness checks.

## Fixed comparison awaiting coordination

| Host | Explicit output | Existing rate | Selected-audio fixture | Order |
| --- | --- | ---: | --- | --- |
| Sequoia / Mac Studio | Mac Studio Speakers | 48000 Hz | selected-clock-48000.wav | M6, diagnostic, M6, diagnostic, M6, diagnostic |
| High Sierra / Golem | Built-in Output, currently headphones | 44100 Hz | selected-clock-44100.wav | M6, diagnostic, M6, diagnostic, M6, diagnostic |

Each run is 100 seconds in Source mode, with all output routes muted, monitor gain
unchanged at −30 dB and no input device. Allow roughly 11–12 minutes per Mac,
including loading. Each host is its own control/diagnostic comparison; different
host sample rates must not be treated as a test changing only the OS.

Before launch, coordinate the order and a hands-off interval with the user. Keep
the test window visible on the physical display. Preserve existing display,
power, connection and workload conditions; do not reboot, change settings or
start an unrelated stress workload within the sequence. Run one Mac's sequence
at a time so any failure or unexpected window state can be assessed before the
other starts. No audible test is included.

`Tests/VideoPaint/run_comparison.py` defaults to hash verification, device/rate
inspection and preview; it launches nothing in that mode. Physical playback
requires both `--run` and a fresh `--coordination-note`, plus a new result folder
outside the package. The runner verifies the exact device/rate before every run,
rejects changed rates instead of setting them, checks mute/source/fixture identity
in each report, records invocation and logs, and never retries. Recorded timing
failures stay in the fixed six-run sequence. Missing reports, invalid traces,
configuration mismatches or timeouts stop the sequence. A timeout terminates only
the specific child process started by the runner. Final package integrity and
audio configuration are checked again.

The same acceptance limits remain in force: zero actual painted-frame skips
after the two-second preroll, zero blanks, frame-bound violations, CPU-paint
deadline failures and underruns, along with existing attachment/queue/memory
limits. These short 100-second runs do not assess the >=120-second memory-growth
gate or replace the required 30-minute endurance qualification after any fix.
Embedded movie audio remains undecoded and unrouted. Measurement behaviour,
selected-audio authority and installed VST3 parity are outside the diagnostic's
scope of change.

After the comparison, retain every result, align skipped-frame events with the
trace, and distinguish selection readiness, main-thread dispatch, complete
harness paint duration and remaining native/compositor uncertainty. Propose a
bounded correction only when supported by that evidence. No merge, push, tag or
M7 publication is authorized by this preparation.

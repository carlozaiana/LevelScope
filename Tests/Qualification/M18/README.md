# Combined M14–M18 qualification

Authorised on 3 November 2026 through the advisory chat's delegation, following the
human's request to start the combined tests. The UI/Timeline worktree owns this
run. Main, checkpoint tags, installations, original media and Pro40 remain unchanged.

The product implementation is commit `4f02303a748b660c72471058111df8a1271758fc`,
following the tested initial-Current integration `cdac0c3`. Additional files here
are qualification drivers and orchestration, with their own frozen hashes. The
measurement build remains `93460e4031c2ec51765b4ad83dd7bc4e67551be1f94129b107ebbeea494a32ce`.
A local qualification commit and immutable candidate package are authorised.
Neither qualification preparation nor a successful build means the tests passed.

## Required fresh checks

`run_short.py` declares and checks every required short group explicitly. On each
actual host it runs measurement, compliance, global and section proposals,
preferences, M18 proposal/history UI, diagnostics/sections, injected audition,
short 7.1, session and restored playback, export, M12 processing, M14 mono/RF64,
M15 markers, M17 requirements, waveform core/service, timeline, peak parity,
video timing core, pinned EBU/ITU vectors, and native stereo/5.1 acceptance with
independent export audits. Reference workers are set to one. Other suites run
sequentially. Use actual arm64 Sequoia and x86_64 High Sierra, never Rosetta as a
replacement for Golem.

Additional required checks are:

- The universal peak-policy probes compare serial, explicitly bounded and ordinary
  app service decisions on identical stereo/5.1 PCM. Intel's ordinary service must
  retain its serial policy. Existing comparison allowances apply only between Macs
  (I/M/S/TP/SP and gains 0.01, LRA 0.1, processed float32 at most one ULP).
  Within one platform the serial/bounded oracle requires exact equality.
- `Exchange.cpp` prepares and reads global/section schema-8 sessions in both
  directions. Verified, LRA/M/S targets, all-No-Limit targets and later-processing
  states preserve exact history, requirements, preferences, roles and geometry.
  Foreign caches must be discarded. Independent receiver measurement, then local
  save/reopen, re-establishes verification only for the measured states.
- The existing M13 video mode covers signed offsets and exact picture selections
  for both audio domains across H.264 CFR/VFR and ProRes 422/LT. `VideoPressure.cpp`
  adds the real main window, native reference window and timeline during analysis
  and injected playback. During normal main-window analysis the picture is held and
  audition remains stopped. Picture skips are recorded; zero skips are not a gate.
- `Long.cpp` runs one actual three-hour mono/192 kHz PCM workflow per Mac. It
  generates full RF64 PCM, compares observation off/on, verifies initial Current
  reuse, imports 1,000 explicit section markers, proposes/applies, independently
  measures Current, inspects long diagnostics, saves/reopens history and stale
  states, checks cancellation and injected Source/Current passages, and renders,
  measures and publishes RF64. The actual timeline component is rendered offscreen
  at 5 Hz with changing viewports. This complements the short native-window checks.
- `check_long_export.py` independently parses RIFF/RF64/ds64, hashes the bundle,
  decodes selected PCM24 directly, reconstructs gain transitions and deterministic
  dither, and checks source/sample positions including the 4 GiB crossing and EOF.
- A separately coordinated muted physical subset uses an explicitly selected
  ordinary two-output device on each Mac, never a default-device shortcut or Pro40.
  Four short cases per Mac cover Original/Current with H.264 and ProRes.
- User native review covers OS file panels, preset Save/Cancel, numerical entry,
  proposal clarity and Seek → Space → Space. Driver-supplied picker results are
  not claimed as OS-dialog navigation or physical listening.

## Budgets and reuse

Keep at least 40 GiB free before the long workflow. The 2,073,600,000 mono frames
occupy 6,220,800,000 PCM bytes per complete copy. Source, ephemeral verified copy,
render staging and published artifact may coexist; retain a conservative four-copy
allowance. Dispose only of owned intermediate data after recording identities.
At most one CPU-heavy qualification workload runs on a host at a time.

The Source/Current overview allocation limit is 128 MiB combined. Extra observed
peak RSS is at most 256 MiB above the equivalent plain measurement, observation
adds at most 15% wall time, timeline preparation/paint p95 is at most 50/100 ms,
cancellation at most one second, core session save/parse at most five seconds,
complete history at most 8 MiB inside a 128 MiB session, and reports at most 16 MiB.
The single combined workflow is limited to 240 minutes on Sequoia and 300 minutes
on High Sierra. A short pilot validates the driver before the three-hour run.

Reuse the existing ThreadSanitizer and adversarial floating-environment development
receipts for unchanged product sources, alongside fresh universal core parity.
Prior full multichannel/cinema evidence is historical complementary coverage,
not a new three-hour multichannel PCM claim. New short real stereo/5.1/7.1 and
three-hour-sized synthetic multi-channel overview stress close the selected UI
coverage. Do not rerun the full LGdC cinema programme or historical long runners.

Any product fix requires a new candidate identity and explicit affected-check
selection. Preserve failed runs, previous packages and successful unchanged
receipts. Publication and main integration require subsequent approval.

## Driver development notes

The 60-second long-workflow pilot passed, with 1,021 independent decoded export
samples matching the source/gain/dither oracle. This pilot is not three-hour or
RF64-over-4-GiB evidence. An exchange preparation pilot also passed. An earlier
exchange pilot used an incorrect input folder and was retained as a setup failure.

The first native video pilot passed its interaction checks but used an inappropriate
256 MiB increment gate for the entire GUI/video process. The driver was corrected
to use the existing 768 MiB total video RSS budget, and to drain its Objective-C
autorelease pool around each manually driven snapshot iteration. This test launcher
has no normal application event loop around those calls. Production event-loop code
is unchanged. The native check explicitly selects Current inspection so actual
progressive PCM is painted during Current remeasurement. All failed pilot receipts
remain in the owned qualification build directory.

Candidate 01's first six short groups passed on both actual Macs. Its multichannel
suite then failed at the same old helper which set `job.state` directly and never
began an analysis job. The current session correctly rejects that missing ticket/
stage. Only `Tests/MultichannelTests.cpp` changed to use the public begin/finish
lifecycle. The corrected focused suite passed 43,202 checks. Candidate 02 retains
the identical product source and app executable. Its runner may reuse only the
first six self-contained passing receipts after verifying both candidate manifests,
identical product-source inventories, report hashes and exact test executable hashes.
All remaining groups execute in the new candidate's result tree. Failed candidate-01
receipts remain untouched.

The user separately approved four 10-second muted physical-output cases per Mac,
using Mac Studio Speakers and Golem's Built-in Output at 48 kHz. No further
permission is needed for that exact subset. System-default selections and Pro40
remain unchanged.

Candidate 02 passed all groups through timeline on Sequoia, then found a marker
migration assertion still expecting the text “M17 reader”. The production warning
correctly names schema 8/M18. The test now checks the current schema number, Save As
and preservation of the original, and its focused rerun passed 487 checks in ten
groups. Candidate 03 changes no product code. Reuse now covers any passing native
executable group with unchanged executable bytes and a verified report, copying
its retained fixture tree for later checks and recording the reuse chain. Changed
executables run again. Reference scripts and plain-text groups run fresh.

The native short driver's worker-job wait is 120 seconds for Intel's bounded
three-pass serial proposal. Visual/interaction waits remain 30 seconds and actual
elapsed times are still recorded. This changes only the diagnostic timeout, not
application behaviour or the qualification resource/response budgets.

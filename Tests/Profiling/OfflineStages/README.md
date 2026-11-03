# Offline job profiling

This diagnostic exercises the development app's serial measurement services on the
development Mac. Use the configured M18 Release build. The harness links its
existing Release libraries and JUCE objects. A second executable overlays four
translation units compiled from the same source with coarse timing probes. The
normal app and its libraries are not rebuilt. The control links the unmodified
service objects. The instruments time block boundaries, whole private passes
and verification stages. There are no per-sample timing calls.

From the repository root, choose a fresh scratch directory, binary directory and
evidence directory. The fixture generator requires a new audio directory and
the runner requires a new runs directory.

```sh
scratch_dir=$(mktemp -d /tmp/levelscope-profile.XXXXXX)
probe_dir=build/m18-offline-probe
evidence_dir=build/offline-profile-evidence
python3 Tests/Profiling/OfflineStages/fixtures.py "$scratch_dir/audio"
python3 Tests/Profiling/OfflineStages/build_probe.py build/m18-refresh-fix "$probe_dir"
python3 Tests/Profiling/OfflineStages/run_trials.py "$probe_dir" "$scratch_dir/audio" "$evidence_dir/runs"
python3 Tests/Profiling/OfflineStages/build_probe.py build/m18-refresh-fix "$probe_dir" --sampling-only
python3 Tests/Profiling/OfflineStages/sample_workers.py "$probe_dir" "$scratch_dir/audio" "$scratch_dir" "$evidence_dir/samples"
python3 Tests/Profiling/OfflineStages/summarize.py "$evidence_dir" "$probe_dir"
```

The programmes have 60 seconds of 48 kHz interleaved 24-bit PCM and 20 equal
sections. A single crest in section 11 produces the local correction. The 5.1
fixture uses the LFE channel for this crest. The regular job has the EBU
programme targets of −23 LUFS and −1 dBTP. The global-only control uses the
same measured programme with a deliberately attainable IL target. Fixture
creation, section setup and disk report writing are outside the timed jobs.

Each repetition validates the original audio, runs Source analysis with both
progressive previews, accepts initial Current from identical Source, generates
a global-only proposal, generates a fresh section proposal, and repeats it with
the exact peak attribution available. It then applies the proposal and runs an
independent Current measurement against the processed samples. The harness
requires the expected pass counts, a passing independently measured Current,
and the same applied plan on both section proposals. It fixes test identities
so the runner can check the exact Source, Current, candidate and proposal
digests and plans across repetitions and builds.

Jobs execute sequentially. Case order alternates between stereo and 5.1 and
between the control and instrumented executable. Each process has a fresh
session and makes new private copies. Input files are warm in the OS cache.
Use control timings for user-visible work and instrumented timings for stage
attribution. Source/Current times measure the service job. Proposal reports
also include synchronous global preflight and acceptance. Dispatch/polling
latency of the real window is outside those totals. Process CPU time includes
the rendering harness. Worker CPU time is obtained from the worker's Mach thread.

During background jobs a 100 ms JUCE timer refreshes the real proposal controls
and gain lane and paints an offscreen 1160 × 460 timeline. Source and Current
publish their actual waveform and loudness progress. The timer's gaps measure
message-loop responsiveness under this load. Real device playback and full
application interaction use different workloads.

Sampling runs use the unmodified service objects with a startup delay in a
separate launcher so the sampler attaches before the first job. The sampling script
retains aggregate inclusive stack counts for the named job workers. For
example, a stack inside the section true-peak observer can also be inside
addPeak and push. Keep these categories separate when interpreting percentages.

## Compare serial and bounded peaks

Configure a fresh Release build with the proposal tests and measurement tests
enabled and `LEVELSCOPE_ENABLE_BOUNDED_TRUE_PEAK=OFF`. Build the named
`LevelScopeSectionProposalTests` and `LevelScopePeakChannelParity` targets. For
the comparison, `--compare-peaks` links the existing service objects and driver
with only the execution choice in the measurement service changed. Both
executables have the same measurement build binding. `app-mode` links the
normal service objects with the build's configured policy.

```sh
python3 Tests/Profiling/OfflineStages/build_probe.py "$peak_build" "$probe_dir" --compare-peaks
python3 Tests/Profiling/OfflineStages/compare_peaks.py decisions "$probe_dir" "$scratch_dir/audio" "$evidence_dir"
python3 Tests/Profiling/OfflineStages/compare_peaks.py timings "$probe_dir" "$scratch_dir/audio" "$evidence_dir"
```

Use a fresh evidence directory for each candidate. Run the parity target and
review performance gates before changing the app build option. An enabled build
gets a distinct measurement build binding. Rebuild the probe against that build
and use the `enabled` stage to compare its normal service with its own serial oracle.

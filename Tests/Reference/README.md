# Offline measurement reference evidence

`itu-bs2217-manifest.json` pins the official ITU download and every extracted WAV
by SHA-256. Expected integrated results and ±0.1 LU tolerances come from the
**file-based** column of ITU-R BS.2217-1 (2012). LFE-only results use `null` to
represent absence of finite loudness. No supplied audio is committed or shipped.

The fetched archive contains 39 mono/stereo/5.1 files, not the newer extended-layout
vectors listed separately on ITU's site. Those extended layouts are outside M1.

After obtaining the archive from the manifest's official URL, verify its hash before
extracting to `build/measurement-reference/itu`. Run from the repository root:

```sh
python3 Tests/run_reference_measurements.py --binary build/scaffold/LevelScopeMeasurementTests_artefacts/Release/LevelScopeMeasurementTests --vectors build/measurement-reference/itu --report build/measurement-reference/results.json
```

The runner refuses missing/modified inputs. Its output includes all measurements,
platform, architecture, manifest hash and binary hash. Add `--compare previous.json`
to compare numeric results across platforms. `--arch x86_64` on macOS requests that
slice of a universal binary; an emulated Sequoia run is **not** High Sierra evidence.
Only integrated loudness has official expected values in this manifest. Other metrics
in the report support reproducibility comparisons, not additional ITU pass claims.
The 39/39 passing baseline for the current test binary is saved as
[itu-bs2217-sequoia-arm64-2026-10-15.json](results/itu-bs2217-sequoia-arm64-2026-10-15.json).

`OfflineMeasurementTests.cpp` independently generates filter/level/peak/LRA signals
from published definitions and exercises the service lifecycle. It does not embed
or copy EBU audio. Generated Tech 3342 cases 1–4 and Tech 3341 true-peak cases 15–23
are distinct from the official EBU v5 archive run described below.

## Official EBU v05 qualification

The user supplied the official v05 archive after the automated download returned
HTTP 403. The 70 extracted WAVs and README were verified byte-for-byte against that
archive, whose ZIP CRC check also passes. `ebu-v05-manifest.json` pins the archive,
README and individual WAV hashes, decoded metadata, explicit file-channel roles,
expected readings and tolerances. These hashes identify the user-supplied copy;
they are not claimed to be independently published EBU checksums. Keep the original
filenames, including the archive's doubled `.wav.wav` suffixes.

```sh
python3 Tests/run_ebu_measurements.py --binary build/scaffold/LevelScopeMeasurementTests_artefacts/Release/LevelScopeMeasurementTests --vectors Tests/Ebu_Loudness_Testset_V05/ebu-loudness-test-setv05 --report build/measurement-reference/ebu-results.json
python3 -B Tests/EbuReferenceRunnerTests.py
```

The local audio directory is ignored by Git; no supplied audio is shipped. The
runner checks hashes before and after measurement, refuses missing/extra WAVs and
changed metadata, and fails missing/non-finite readings or incomplete time-series
coverage. Results include the binary/runner/manifest hashes, kernel/build identifiers,
actual executable architecture, every criterion and the checked time-series readings.
Use `--compare prior-report.json` for cross-platform numeric comparisons and optional
`--arch x86_64` to select that binary slice on macOS.

The saved result is [ebu-v05-sequoia-arm64-2026-10-15.json](results/ebu-v05-sequoia-arm64-2026-10-15.json):

- **66/66 minimum-requirements files pass**, covering all Tech 3341 cases 1–23 and
  Tech 3342 cases 1–6, with 72 metric rules and 1,052 checked readings.
- **4/4 supplementary files pass** nominal loudness sanity checks at a project-chosen
  ±0.1 LU tolerance. These are three alignment sines and the reference listening noise;
  they are not four additional EBU minimum-requirements cases.
- Cases 1/2/9/12 check every 10 Hz curve point in the specified steady interval.
  Cases 11/14 check 20 successive cumulative maxima at 6 s / 0.8 s boundaries. The
  harness finishes a copy of the engine at each boundary, retaining uninterrupted
  filter state in the original. This verifies the offline kernel on those signals;
  it does not qualify a live display or scheduler.
- All supplied WAVs are 48 kHz. Generated 44.1/96/192 kHz checks remain separate evidence.
- Eight runner tests exercise asymmetric limits, non-finite readings, incomplete
  time series and failed/changed comparison evidence. No DSP correction was needed.

Normative values come from [Tech 3341 v4, Table 1](https://tech.ebu.ch/docs/tech/tech3341.pdf)
and [Tech 3342 v4, Table 1](https://tech.ebu.ch/docs/tech/tech3342.pdf). This result
qualifies these numeric tests on arm64 Sequoia, not the whole product or delivery profiles.

High Sierra qualification: transfer/build the universal test binary and run the
hardware-independent suite there, then run both reference runners against the same
pinned vectors and compare their reports with Sequoia. Record exact OS/toolchain/CPU. These
tests do not require an audio device or change system audio configuration.

Silent UI smoke checklist: import a supported file, confirm roles, Measure Source,
initialize Current State, re-measure at 0 dB, change trim to -6 dB (old values must
clear), re-measure, cancel an in-flight job, replace/clear the source, and close while
measuring. Confirm all routes remain muted, the device stays closed, and export stays
disabled. Repeat with AIFF/FLAC and verify a 12-channel source imports but reports
unsupported measurement. No playback is needed for this checklist.

## Portable final-qualification run

Completed 17 October 2026 on actual High Sierra 10.13.6 x86_64 and Sequoia 15.7.9
arm64. Both official suites, 184 measurement checks, 205 scaffold checks, eight runner
tests and both silent UI checklists passed. Exact evidence is in
[`results/m1-qualification-2026-10-17`](results/m1-qualification-2026-10-17);
see [the qualification report](../../Docs/ProjectIndex/17_M1QualificationResults.md).

`LevelScopeMeasurementTests --write-ui-fixtures /absolute/new-directory` generates
WAV/AIFF/FLAC stereo tones, a ten-minute WAV for cancellation/close tests, and an
explicitly tagged 12-channel WAV outside M1's measurement layouts. It refuses an
existing directory. The generated audio stays in ignored build storage.

After committing the tested source and reports, prepare an isolated folder:

```sh
python3 -B Tests/package_m1_qualification.py --output build/m1-qualification-2026-10-17 --fixtures build/m1-ui-fixtures-2026-10-17-v2 --ebu-baseline Tests/Reference/results/ebu-v05-sequoia-arm64-2026-10-17.json --itu-baseline Tests/Reference/results/itu-bs2217-sequoia-arm64-2026-10-17.json
```

This copies only the standalone app, two test executables, runners/manifests, pinned
reference WAVs, generated UI fixtures and numeric baselines. `identity.json` records
the committed source tree/files, compiler/build settings and every delivered artifact's
SHA-256. The directory contains external audio/build products and must never be committed.
Inspect the older Mac's OS, CPU and available Python/runtime dependencies before transfer
or execution; the current runners need Python 3.7 or newer. Preserve installed software.

From the transferred folder, using the inspected Python interpreter:

```sh
sh Tests/run_m1_bundle.sh /absolute/path/to/python3 high-sierra results/high-sierra-2026-10-17
```

Use `sequoia` for the arm64 local run. The platform guard requires actual 10.13.x/x86_64
or 15.x/arm64 respectively; it rejects an Intel Sequoia substitute. The runner verifies
the complete delivered artifact inventory and runs both executable suites, the runner
tests and both official suites with comparisons to the pinned Sequoia reports.
The folder receives logs plus OS build, CPU, architecture, dependency and source/build
identity evidence. The silent UI checklist remains a separate observed/manual check.

For transfers to older macOS, Python's archive writer avoids Sequoia-specific extended
attributes that High Sierra's tar warns about. From `build`, before adding run results:

```sh
python3 -m tarfile -c m1-qualification-2026-10-17.tar.gz m1-qualification-2026-10-17
```

Verify the archive hash at the destination and require the bundle runner's complete
artifact-hash check to pass before accepting any test evidence. Do not transfer SSH
credentials or install the bundle over an existing working app/plugin.

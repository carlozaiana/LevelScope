# M18 progressive loudness and timeline navigation — development evidence

The previous waveform slice remains in `../m18-waveform-2026-10-31`. This slice runs
only the standalone and selected audio-only tests, using Release universal arm64 /
x86_64, Intel deployment target 10.13 and the existing JUCE 7.0.12 checkout.

`run_focused.py <new-output-directory>` records hashes of sources, tests, build
configuration and binaries. It checks the build identity on every report and runs
five named tests: TimelineProgress, Measurement, Diagnostics, Waveform and
JobRequirements. It also runs the four local pinned EBU v05 cases with measured
M/S curve acceptance rules through the existing reference checker. Fixtures are
new owned directories in `/tmp`; JUCE verified-copy cache creation requires the
existing macOS sandbox escalation. Only the test double for an audio output is
opened. No physical device, video or network access is required.

The new interaction suite examines a real analysis worker before EOF, cancels it,
and verifies that both previews are retired without a completed measurement. It
exercises the existing stopped seek service with no open device, fractional numeric
input, the exact last ruler pixel/EOF, preparation/playback/job/loading/stale guards,
independent Current inspection and Original audition, drag vs click, navigation,
Target availability and the two divider bounds. A synthetic three-hour-sized M/S
window set checks bounded storage and preparation time without a three-hour decode.
Waveform integration additionally compares every preview window end and energy
against independently completed observer-off results at all supported rates and
layouts, including Current processing and transitions, assembly and composition.

The first development pass (`iteration-01`) is superseded by the final candidate.
It caught no transport/measurement failures. Subsequent refinements keep the old
default combined comparison/info height, put both preset controls into the unused
header cell only at maximum compression, preserve exact inspected frames behind
the decimal display, and retain a consistent right LUFS axis without a waveform.

Full native interaction, actual Intel High Sierra and cross-Mac comparison,
physical audio/video and the combined full-duration/RSS workflow remain deferred.
The sample-level visible-passage reader, gain/section inspector integration and
saved applied explanations are outside this slice. No main integration or plugin
installation is part of this development run.

### Final development result

[Final summary, reports and hashes](final-summary.json): **2,692 checks passed**
across timeline/progress (343), measurement (185), waveform integration (1,305),
diagnostics (370) and Target/presets (489). All four pinned EBU v05 M/S curve
files passed their original cadence, value and tolerance checks. The Source and
processed Current preview end frames and energies matched the final windows
exactly across the supported rates/layouts.

The synthetic three-hour preview used 3,455,488 bytes for 215,968 real-cadence
test windows, with p95 preparation of 2.273 ms. This tests viewport preparation
and bounded storage, not a three-hour decode. On the short integration fixtures,
the observed/plain analysis time ratio was 1.000 to three decimals, with a
largest sampled viewport preparation plus software paint time of 4.695 ms.

The [incomplete](candidate-01/partial.png), [completed](candidate-01/navigation.png)
and [compressed comparison](candidate-01/compact.png) component captures were
visually inspected. The first candidate's diagnostics assertion compared the
entire component, including the newly enabled **Fit selection** button. The
[diagnostics rerun](diagnostics-rerun/diagnostics.json) checks equality of the
actual plot for an offscreen selection and a visible-selection positive control.
Only that test changed. The app, all production sources and the other test
binaries were checked against the original hashes and did not change.

Build identity: `aa06a460f743ba647f966f886be4cb76e9abbddb1588013f67eca522ee8c2cbe`.
Development app: `/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m18-dev-01/LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app`.

### Manual interaction checklist

1. Open and confirm a mix long enough to watch analysis. On **Inspect Original**,
   the 400 ms M and 3 s S curves appear as their measured windows arrive. Unread
   time is pending and the M gate is pending until acceptance. On the first
   successful import, Original and untouched Current then show Ready. Short files
   can finish between screen updates.
2. Choose **Inspect Current** and **Audition Original**. Click a channel or enter
   seconds (e.g. `1/7`) and press Return. Check the inspected frame, bin and M/S
   windows; the white audition playhead should stay put. While stopped, click the
   upper **Timecode** ruler or press **Seek**: the playhead moves without playback
   or a domain switch. This also works with the device closed. During playback,
   a ruler click must leave audition running from its existing position.
3. Use **Fit selection**, **Zoom +**, drag to pan, or use the scrollbar; **Whole
   programme** restores the full clip. Try **Follow** during playback, then pan
   manually to release it. Check the timecode above the plot and seconds below it.
4. Compress the comparison rows and drag the divider below the scrollbar. The
   metric font and Target help line should stay readable; the lower inspector
   should gain the height released by the graph. **Session / tools → Reset layout**
   restores both. Play and Stop, or finish a trim/proposal job: the Target dropdown
   and **Save as preset** should unlock without discarding a draft or changing
   programme applicability.

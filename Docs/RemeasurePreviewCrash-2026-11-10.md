# Re-measure crash after a gain edit — 10 November 2026

The direct-editing candidate-02 is superseded for use by the corrected package below.
All earlier packages remain preserved; nothing is published or merged into main.

## Cause and reproduction

The user's report matches candidate-02's arm64 UUID
`1C55E59F-16B6-350E-85E7-24A310783EC9`, built from
`082d45ff7d88768f44bd9b557f16da561f2b7d93`.

Starting Re-measure cleared the timeline's waveform pointer while leaving its gain
preview binding alive. A paint before the next timer refresh then read `binFrames`
through the null pointer while formatting the preview heading. The measurement
worker was still preparing its verified source; the crash was on the UI paint thread.

A synthetic programme with injected output reproduced the report: set global gain
to +5 dB with Inspect Current showing its gain preview, click Re-measure Current
State and paint immediately. Before the fix the process exited with SIGSEGV (139),
fault address `0x10`, and the same `MeasurementTimeline::paint` offsets +776 and +524
as the user's report. The saved diagnostic summary omits the report's machine and
peripheral identifiers.

Earlier tests covered gain edits and successful remeasurement but did not force a
paint between clearing the Current overview and the next timer refresh.

## Correction

Retiring waveform data now also retires its gain binding and cached gain key.
Changing programme identity clears those bindings and accepted-overview state.
Rendering and provisional sample-peak authority require actual waveform data as
well as a gain binding. A staged binding with no data paints safely as unavailable.
The sections window binds the new programme before supplying its new gain preview.

Only three product files change: MeasurementTimeline.cpp/.h and the sections-view
binding order in StandaloneMainComponent.cpp. Measurement, processing, proposal,
export, routing and device behaviour are unchanged.

## Focused verification and build

Preserved corrected app:
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/ui-remeasure-crash-2026-11-10/candidate-01/LevelScope Standalone.app`

Exact tested commit: `5ac6ac57530004ab28167e65e63fda7475b09b68`.
Measurement build: `cbe021a9ba51b864d8234bd1a66d5a04d459e5738da7a26825eeeabfdb10e271`.
The package is universal arm64/x86_64 with Intel minimum macOS 10.13. Its Golem copy
is `/Users/fuzz/LevelScope-UI-remeasure-crash-2026-11-10/candidate-01/LevelScope Standalone.app`.

| Focused group | Sequoia 15.7.9 arm64 | High Sierra 10.13.6 x86_64 |
|---|---:|---:|
| Timeline, gain-preview lifecycle and existing editing/navigation | 2,001 pass | 2,001 pass |
| Main UI, immediate repaint, remeasurement and cancellation | 153 pass | 153 pass |
| Total | **2,154 pass** | **2,154 pass** |

Both actual hosts ran the same immutable package with synthetic audio, isolated
settings and injected output. No physical audio interface was opened. Screenshots
confirm that the formerly crashing transition renders preparation/cancellation and
unavailable peak data without presenting an incomplete measurement as verified.

The regression paints synchronously after +5 dB and Re-measure, then verifies that
measurement completes for the exact applied +5 dB plan. It also checks cancellation
immediately after another gain edit. Timeline checks cover cleared data, a staged
binding awaiting data, programme replacement, unavailable peak data and rebinding
the same numeric gain key. Existing unaffected target, core measurement/DSP,
monitoring/device, video, long-duration and listening qualification is reused.

[Acceptance receipt](../Tests/Reference/results/ui-remeasure-crash-2026-11-10/acceptance.json)
contains package identity, actual-host receipts and preservation checks. All 428
artifacts in the previous candidate-02 remain byte-identical. `main` and `origin/main`
remain at `22961f03de7919e1de631f8b13eb4b852108ac0b`; publication is not performed.

## Manual retry

1. Open the corrected app and the same session; select **Inspect Current**.
2. Set global gain to **5 dB** and click **Re-measure Current State**. Preparation
   should remain responsive and measurement should complete.
3. Change the gain again, start Re-measure and cancel during preparation. Current
   should remain **Needs re-measurement**, with playback stopped.

The waveform can be temporarily unavailable while preparation has no overview data;
progress and cancellation remain available. This patch does not add a new audio
read or change the combined 128 MiB waveform budget.

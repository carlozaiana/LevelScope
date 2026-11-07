# Timeline gate membership follow-up — 7 November 2026

Continuation on `codex/m20-ui-integration`. This corrects the interpretation of
request 11 in the [earlier inspection handover](UIInspectionConvenience-2026-11-07.md)
and replaces its paired-switch shortcut. Earlier frozen packages and listening
sessions remain untouched. Publication and main integration remain pending.

## Behaviour

The relative IL line remains visible in both Original and Current. Only Current
adds red, double-width segments, and only where the corresponding measured
momentary windows change IL gate membership: included to excluded, or excluded
to included. Each domain uses its own completed measurement's gate, with exact
energy comparisons and the existing strict IL population rule. An amplitude
change with unchanged membership is not red. An unchanged momentary value can
be red when the programme gate moves past it.

Comparison remains at the actual 100 ms measurement cadence. A coloured cell
represents the cadence interval ending at that measured window; it does not
claim an interpolated sub-cadence crossing time. Cells remain visible when zoom
clips their window end beyond the viewport. Missing/incomplete or unverified
measurements cannot produce red comparisons. The Original relative line and
waveform/loudness blend are unchanged.

Plain **A** switches inspection and audition together while stopped and idle,
in the main and reusable Sections timeline. Text editors and modal dialogs keep
keyboard input; active playback/preparation/jobs keep their existing guards.
The previous Cmd+Shift+A chord no longer invokes this action. Play remains explicit.

Short-term colours are unchanged by this follow-up. Their classification still
comes from the measured LRA population at full-resolution window ends. The
short-term series itself remains ungated. The focused comparison audit records
finite displayed S windows, exclusions, minimum S and the LRA relative gate for
Original and Current separately; EOF-only drain windows are not displayed.

No DSP, measurement engine, proposal search/decisions, history schema, export,
monitoring policy or waveform storage budget change. Source-file hashing gives
the new app a new measurement build identity. Fresh comparison copies are
therefore independently remeasured, checked against the previous six scalars
and sample/true-peak positions, and saved separately with the original complete
historical explanations and exact processing plans.

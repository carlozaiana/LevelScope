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

## Measured comparison findings

Each of the three distinct states has **2,727 finite displayed short-term
windows and zero exclusions**. The Original result is identical in both sessions.

| State | Lowest displayed S (LUFS) | Relative LRA gate (LUFS) |
| --- | ---: | ---: |
| Original | -41.22650001310821 | -41.95005352865540 |
| M19 Current | -42.19079245304594 | -44.24610304240639 |
| M20 Current | -42.86442383015300 | -44.30113281687761 |

Thus the absence of blue-grey in these examples is expected, independently of
colour perception. The existing synthetic regression also verifies an actually
excluded S window remains in the excluded display category.

M19 contains 88 changed-membership cadence cells; M20 contains 59. Both include
the cells ending at 131.6, 261.7, 261.9 and 262.0 seconds, changing from excluded
Original to included Current. These are actual sampled-window classifications,
not a claim that every interpolated point in the wider passage crosses a gate.
M19's IL relative gate is -31.95981413624387 LUFS for Original and
-34.25312986301650 LUFS for Current. At the 131.6-second window end, M is
-32.31645274536555 LUFS (Original, excluded) and -33.28074519078556 LUFS
(Current, included).

The first follow-up package (candidate 03, commit
`bb1bd06a2b27a88d592901c26016fe9d1732ec3d`) passed 919 timeline and 105 native UI
checks on each actual Mac. Its independent comparison refresh passed 185 checks.
Visual inspection then found the legend still described the previous
amplitude-difference rule. The final package corrects only that legend string;
the gate/shortcut implementation is unchanged. The final timeline rendering is
checked again, while unchanged native workflow evidence is reused.

## Final preserved build and focused qualification

- Final package: `build/ui-inspection-2026-11-07/candidate-04`.
- App: `/Users/fuzz/.codex/worktrees/a411/LevelScope/build/ui-inspection-2026-11-07/candidate-04/LevelScope Standalone.app`.
- Exact source commit: `6120b65caf9d88219ad5ded4914a6d070776f149`.
- Measurement build: `af2af77617baf191ba57d1bc9393ef39c818e566dbf9a33e36a030c4ac1f3bea`.
- App executable SHA-256: `275c909a324b02bc8dde301f76fd9760cb02b19976dd5447819579dbc96384d6`.
- Universal arm64/x86_64; Intel minimum macOS 10.13.

Final timeline checks: **919 pass on Sequoia 15.7.9 and 919 pass on actual High
Sierra 10.13.6**. The preceding candidate's **105 native UI/shortcut checks on each
Mac** are reused: the preservation audit verifies that the only product/test
source difference is the corrected legend string. Both final rendered gate
screenshots were visually inspected. Independent refresh of both listening copies
passes **185 checks**, reproducing the six measured scalars and all sample/true-peak
positions exactly, preserving gains/requirements/explanations, and reopening with
verified Current. The final gate reports exactly match the first follow-up audit.

The existing cold-setup/large-history heartbeat limitation and all prior
qualification deferrals remain unchanged. No physical devices were opened. All
four inspection packages and both original listening copies have verified hashes;
main/origin remain `22961f03de7919e1de631f8b13eb4b852108ac0b`. No publication or
main integration occurred. Listening quality remains for the user's review.

New session copies, already applied and independently verified:

- M19: `/Users/fuzz/.codex/worktrees/a411/LevelScope/build/ui-inspection-2026-11-07/listening-gate-04/m19/inspection.levelscope`.
- M20: `/Users/fuzz/.codex/worktrees/a411/LevelScope/build/ui-inspection-2026-11-07/listening-gate-04/m20/inspection.levelscope`.

Manual check: open the new app and one new session. While stopped, press A to
switch Original/Current. Zoom near 131.6 s and 261.7–262.0 s: the relative line
remains, with red membership-change cells only in Current. Type A into a text
field to check that text entry remains independent. Blue-grey S exclusions are
absent in these particular sessions; no Apply or proposal regeneration is needed.

[Exact receipts and preservation audit](../Tests/Reference/results/ui-gate-followup-2026-11-07/acceptance.json).

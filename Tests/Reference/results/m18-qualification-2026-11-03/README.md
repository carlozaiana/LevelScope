# Combined M14–M18 qualification, 3 November 2026

**Passed for the agreed combined M14–M18 scope.** Candidate 04 passed the short,
reference, selected video, muted physical, native-review and bidirectional session
checks on both actual Macs. Both also passed their actual three-hour workflows,
independent RF64 audits, complete numerical comparisons and supplemental synthetic
display-density checks. The final cross-Mac comparison passed. Main integration
and publication have not been performed.

The user authorised the combined qualification through the advisory chat and
separately approved four ten-second muted physical cases per Mac. The UI/Timeline
worktree owns this run. [Qualification selection and driver notes](../../../Qualification/M18/README.md)
describe the commands, budgets, earlier failures and corrections.

## Frozen application

| Identity | Value |
| --- | --- |
| Source/product commit | `3665815d73fbd93c0574cdfee8d7efdd805ab7b6` |
| Measurement build | `79ee7b68ded355f1baf36f839c978ad692b8649cb3fe669b9fc9988060737fdb` |
| Candidate identity SHA-256 | `33bc7202fc2629797bbfd3223b3849145a88bebe5c22e9157a1636f53cd69e95` |
| Universal app executable SHA-256 | `1b9b8db2aca5d17f7343912110d263f638216152f1bb806368220a4ac79709bf` |
| Development host | M1 Max, 10 cores, 32 GiB, arm64 Sequoia 15.7.9 |
| Golem | Intel iMac14,2, 8 logical CPUs, 32 GiB, x86_64 High Sierra 10.13.6 |
| Intel deployment minimum | macOS 10.13 |
| Main, unchanged | `a9be070ab9eb5cb4f591ff9d991d8999c3a016d0` |

Exact development app:
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m18-qualification-2026-11-03/candidate-04/LevelScope Standalone.app`.
The same frozen app on Golem is
`/Users/fuzz/LevelScope-M18-qualification-2026-11-03/candidate-04/LevelScope Standalone.app`.

The frozen manifest's `qualification_complete: false` records its status at freeze;
it is not rewritten after testing. Subsequent qualification receipts establish
the tested outcome. Test-only checker/probe changes have separate provenance.

## Completed fresh checks

All 28 declared short/reference groups passed on each actual host. They cover
measurement, compliance, global/section proposals, diagnostics, sections,
audition, short real multichannel, sessions, restored playback, export,
section processing, preferences, proposal/history UI, mono, requirements,
waveform service/core, timeline, markers, peak parity, video core, EBU/ITU,
native stereo/5.1 workflows and independent export audits. There was no reuse of
older candidate groups in this candidate's short run.

| Result | Sequoia | High Sierra |
| --- | ---: | ---: |
| 28 short/reference groups | 399.64 s | 1,075.65 s |
| Pinned EBU cases | 70 passed | 70 passed |
| Pinned ITU cases | 39 passed | 39 passed |
| Injected video matrix | 80 cases passed | 80 cases passed |
| Native-window video/analysis pressure | 8 cases passed | 8 cases passed |
| Native video peak RSS | 560.94 MiB | 225.47 MiB |
| Muted physical video | 4 cases passed | 4 cases passed |
| Receiving session variants | 8 passed | 8 passed |
| Explanation inspection heartbeat p95 | 19.98 ms | 89.27 ms |

The complete cross-Mac comparison, including the actual three-hour job, passed
with exact Source/Current peak positions, proposal geometry and verdicts.
Captured short processed float32 samples differed by
0 ULP; the largest compared metric difference was `3.552713678800501e-15`.
Existing cross-Mac tolerances remain unchanged. Within-host serial/bounded
comparisons require exact equality.

The original Intel service probe's `bounded` label was incorrect because its
private build retained the production arm64 restriction. The corrected diagnostic
exercises bounded workers on Intel and matches serial/app measurements, peak
positions and all seven proposal outcomes exactly for both layouts. The ordinary
Intel application remains serial. Corrected probe timings were 13.01 s versus
21.65 s serial for stereo, and 22.65 s versus 63.49 s for 5.1. These are diagnostic
probe totals, not the separate native workflow timings below.

## Fresh native workflow timing

Each job is 60 seconds at 48 kHz. The stereo case uses a global proposal; the 5.1
case uses a section proposal with three private candidate passes. Native controls
are driver-assisted; output timing uses injected callbacks, not physical latency.

| Interaction | Sequoia stereo | Sequoia 5.1 | High Sierra stereo | High Sierra 5.1 |
| --- | ---: | ---: | ---: | ---: |
| Import to first partial waveform | 0.639 s | 1.222 s | 2.279 s | 3.219 s |
| Confirm roles to first waveform | 0.263 s | 0.662 s | 0.692 s | 1.316 s |
| Import to measurements ready | 1.467 s | 2.150 s | 5.235 s | 10.727 s |
| Play to first non-zero output sample | 0.304 s | 0.704 s | 0.491 s | 1.187 s |
| Import to first non-zero output sample | 1.770 s | 2.854 s | 5.726 s | 11.914 s |
| Propose | 0.019 s | 4.102 s | 0.088 s | 23.891 s |
| Independent Current remeasurement | 0.993 s | 1.502 s | 3.046 s | 8.357 s |

## Actual three-hour mono/192 kHz workflow

Both Macs passed all 87 workflow checks: 68.68 minutes on Sequoia and 123.31
minutes on High Sierra. The programme contains
2,073,600,000 actual PCM frames; each full RF64 file is 6,220,800,104 bytes. It used
1,000 explicitly imported sections, one local correction covering 10.9 seconds
including transitions (0.100926%), and a separate +1.970132700246 dB global gain.
The measured private candidate passed at -23 LUFS and -1.498709 dBTP. Apply left
Current unverified until the subsequent independent measurement passed.

| Stage or resource | Sequoia | High Sierra |
| --- | ---: | ---: |
| Plain Source measurement | 504.34 s | 950.30 s |
| Source with waveform/M/S observation | 509.82 s | 948.65 s |
| Observed versus plain wall time | +1.09% | -0.17% (within timing noise) |
| Three private proposal passes | 1,696.35 s | 2,943.88 s |
| Independent Current measurement | 519.12 s | 963.54 s |
| Core session save / parse | 0.110 / 0.483 s | 0.235 / 1.026 s |
| Session reopen including source verification | 37.16 s | 59.31 s |
| Render and independently measure artifact | 642.33 s | 1,185.79 s |
| Publish verified artifact | 100.61 s | 152.36 s |
| Combined workflow | 4,120.69 s | 7,398.31 s |
| Additional full-result equality check | 545.03 s | 991.96 s |
| Workflow plus numerical closure | 77.76 min | 139.84 min |
| Cancellation acknowledgement | 10.50 ms | 36.35 ms |
| UI heartbeat p95 | 15.03 ms | 21.82 ms |
| Timeline preparation / paint p95 | 1.40 / 6.91 ms | 3.84 / 25.69 ms |
| Source/Current overview allocation | 100.70 MiB | 100.70 MiB |
| Additional observed Source peak RSS | 85.62 MiB | 66.72 MiB |
| Session / JSON+HTML report bytes | 10,725,495 / 872,091 | 10,725,481 / 872,175 |

The Source observation comparison used peak RSS of 72,777,728 bytes plain versus
162,562,048 bytes observed on Sequoia, and 43,945,984 versus 113,909,760 on High
Sierra. The separately reported whole-workflow RSS growth was 303,611,904 /
222,720,000 bytes; it includes proposal, session and export work and is not the
incremental-preview comparison. The actual timeline rendered offscreen at 5 Hz,
with 20,305 / 36,287 updates and 4,411 / 8,273 partial-preview updates respectively.

The first partial waveform arrived 71.36 seconds (Sequoia) / 120.66 seconds
(High Sierra) after the observed analysis started, including full-payload preparation. This is not an import-to-waveform
native-window timing and does not replace the short-job timings above. Faster
initial waveform availability remains a separate follow-up.

Each independent audit streamed all bundle hashes and reconstructed 1,030 selected
PCM24 samples, including gain transitions, EOF and 314 samples beyond 4 GiB. All
matched exactly. The two full exported RF64 files have the identical SHA-256
`055ae11446acfe0e9c405102f0a78f403d017075bb561c69435c27d908888449`. The separate Source/Current injected playback check matched all
eight selected 4,096-frame passages exactly, reused the prepared copy and recorded
zero underruns. No physical device was opened by this long workflow.

The original long driver's observation assertion compared input identity/content
and peak positions; it did not compare every numerical result. The separately
frozen `ExactLongSource.cpp` closes that evidence gap using the existing full
short-test oracle. On each Mac, one additional plain Source pass matched every
metric/status, gate, population, percentile, peak value/position and all 323,954
M/S/LRA windows exactly against the verified observed Source cache. This added
545.03 / 991.96 seconds including session verification (507.43 / 931.88 seconds
measuring). It reused each existing file and did not repeat proposal, Current,
playback or export. The strengthened final cross-Mac checker required and passed
these receipts on both hosts.

The supplemental synthetic display check passed 238 assertions on each Mac. Both
domains together occupied 77.53 MiB for six channels and 100.70 MiB for eight.
Preparation/paint p95 was 1.74/24.29 ms and 1.90/26.26 ms respectively on
Sequoia, and 2.83/79.00 ms and 3.03/88.81 ms on High Sierra. These are
three-hour-sized synthetic extrema and loudness data, not three hours of
multichannel PCM decoding. Its separately frozen universal executable is recorded
in `receipts/density-01-provenance.json`.

## Session and native review

Both directions tested global/stereo and section/5.1 sessions in verified,
LRA/M/S-target, all-No-Limit and later-processing-unverified states. Complete
encoded historical explanations, requirements, preferences, roles and geometry
survived exactly. Foreign measurement caches were discarded; receiver-side
independent measurement and a subsequent save/reopen established valid local
verification. History status remained separate from Current verification.

Actual native Open/Save panels were used on both Macs. Preset Cancel closed;
preset Save closed and selected the complete definition, including No Limit and
M/S values. Existing presets were preserved, and the disposable test presets were
archived. Numerical entry changed the dropdown to Unsaved preset. Reopening the
owned 5.1 session showed Ready/PASS, all three explanation pages were accessible,
Propose reported Already meets requirements, and Proposal limits opened in front
and cancelled edits correctly. Save As created distinct files. Apps were closed
after the reviews. This is computer-driven native review, not human aesthetic
approval or physical listening. Seek/Space/Space key dispatch is covered by the
injected transport/native-component tests.

A cosmetic issue remains: the historical explanation heading renders its dash
incorrectly. The words, numerical results and complete structured explanation
remain accessible. It was recorded without changing the frozen candidate.

## Muted physical and video scope

Four ten-second cases per host used explicitly selected Mac Studio Speakers or
Built-in Output at 48 kHz, with every route muted: H.264 CFR and ProRes 422, each
with Source and Current. All eight had zero audio underruns, skipped/blank frames
and timing-bound failures. System-default device selections stayed unchanged;
Pro40 was not opened. Golem's chosen output ran at the authorised 48 kHz.

The separate injected/native-window matrix includes H.264 CFR/VFR and ProRes
422/LT, signed offsets, seek, stop/resume/restart and EOF. Analysis paints actual
partial waveforms while picture is held; audition is stopped during ordinary
analysis. Native picture/output checks then exercise the playback path. The
native pressure cases recorded late polls but zero skipped frames or underruns.
The agreed allowance for occasional picture skips was not needed in these runs.

## Evidence and remaining scope

`receipts/` preserves raw JSON, logs, selected captures and candidate manifests;
`receipt-manifest.json` records their sizes and SHA-256 hashes. Full owned PCM,
RF64 exports, sessions and frozen executables remain in each host's qualification
directory. Failed runs and superseded packages are retained, not overwritten.
The pagination probe records why candidate 04 reduced explanation pages from
12,000 to 4,000 bytes while preserving complete text. Test lifecycle, migration
assertion, cross-Mac checker and Intel diagnostic corrections are documented in
the driver notes and separate provenance receipts.

All 395 frozen artifact hashes were rechecked successfully on each Mac. Golem's
Python audit runs created two ordinary `__pycache__` files outside that frozen
inventory; they are recorded in its integrity receipt. No hashed source or
executable changed. [Machine-readable final summary](receipts/qualification-summary.json)
and [complete cross-Mac comparison](receipts/cross-mac-complete-04.json) record the
final gates. Main and all 21 existing tag objects/targets remain unchanged.

Full physical 7.1 listening remains separate. A complete codec/colour matrix and
30-minute video endurance remain conditional on timing/buffering/API changes or
concrete failures. No new full-duration multichannel PCM run is claimed. The
visible-passage sample reader, device reconnection, stereo monitoring fold-down,
future DSP/automation/meters and automatic LRA/M/S correction remain later work.

For a brief manual check of the candidate:

1. Import stereo then 5.1, inspect the proposed roles and their provenance, and
   confirm analysis. Inspect the initial Current measurement origin.
2. Save and cancel a target preset; check that each dialog closes and direct
   edits select Unsaved preset.
3. Drag a section, adjust its boundaries and numerical gain, and inspect the
   draft before explicit Apply or Discard.
4. Seek while stopped using the timecode ruler, then use Space twice for Play/Stop
   in the main and Sections windows, with an explicitly selected output.
5. Review and Apply a proposal, confirm Current becomes unverified, then remeasure.
   Save As and reopen to inspect the last-applied explanation.

Audio-device choice and any audible playback remain explicit user actions.

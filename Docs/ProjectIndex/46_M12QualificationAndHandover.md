# M12 section processing — implementation and Sequoia-phase handover

M13 qualification update, 28 October 2026: bounded candidate-02 is implemented
and qualified on actual Sequoia/High Sierra; [ProjectIndex 47](47_M13BoundedSectionTrimProposals.md) is authoritative for this scope. The broader
no-change rule covers all applicable explicitly selected, method-defined objectives;
M13 still evaluates only supported programme rules. The intended cinema comparison
is user-reported Nuendo statistics supplied 28 October, while the historical screenshot
is Cubase 10 supplementary evidence. No exact Nuendo version is claimed. The user authorized checkpoint publication
and combined M12/M13 main integration on 28 October. The preserved M12 qualification
status remains unchanged; [49](49_M13QualificationAndHandover.md) records the completed combined
M12/M13 qualification. The [published M13 checkpoint](https://github.com/carlozaiana/LevelScope/tree/m13-bounded-section-proposals-2026-10-28)
identifies the combined main revision at publication.

Status: **“Sequoia tested; High Sierra and video qualification pending.”**
The approved Sequoia phase is complete on candidate-03, 27 October 2026.
This is not full two-platform/video qualification. The user subsequently authorized
publication of this development checkpoint, with main held at M11.

At the M12-only publication, main remained the published M11 commit
`cfb1033aa7fb65a2655db276bee696161196bdeb`. Candidate-03 is saved on
`codex/logic-module` as annotated checkpoint
[`m12-section-trim-sequoia-2026-10-27`](https://github.com/carlozaiana/LevelScope/tree/m12-section-trim-sequoia-2026-10-27).
The checkpoint commit is `2fb0e3cfee22a527bacaa09c4ddace23129442fb`; direct GitHub
verification confirmed all fifteen earlier tags unchanged and main still at M11.
The [publication record](../../Tests/Reference/results/m12-publication-2026-10-27/README.md)
binds the tested manifest and original binaries to that checkpoint. At that
publication, main integration and M13 implementation were not authorized. M13 was
subsequently approved on 28 October; its separate candidate/evidence record is
[49](49_M13QualificationAndHandover.md). The user subsequently approved combined M12/M13 publication and main integration on 28 October; see the M13 publication record.

Approved contract: [44](44_M12SectionTrimProcessing.md). The pre-implementation
change-impact and resource plan is [45](45_M12ImpactQualificationPlan.md).

## Frozen candidate

Candidate-03 is preserved at
`build/m12-qualification-2026-10-27/candidate-03`. Its source manifest is
`d6501c0d234ccf9c02b90527fc4af6c613dd4d538a71cc5520fb6b7b21085c1f`;
measurement build identity is
`185566e2c3ed9835d9878643618bd6789e785dee78fd670174f1ede4eb10f5fd`.
The standalone executable SHA-256 is
`ab319348a5ca586fd9048581165f7fdefb0cb1bd7221d827c78fd3740ec2aec4`.
This identifies the tested source above M11, saved by the checkpoint tag; it is not
source already present in M11. The original package retains its pre-publication
build metadata and bytes. Documentation and evidence are outside the compiled
source manifest, and publication does not substitute a later rebuild.

The package contains arm64 and x86_64 slices. The Intel executable declares macOS
10.13; the established JUCE 7.0.12 dependency is unchanged. This is local build
evidence, not actual High Sierra execution. No plugin target was built or installed.

## Implemented behavior

The manual Programme sections editor now creates a separate gain draft. Gains
are finite doubles from -12 to +12 dB. Green shows applied settings; orange shows
the draft. The curves are processing settings on audio-source time, not measured
loudness. Set draft gain does not change audio. Apply requires stopped playback
and no conflicting job. Different applied processing discards the previous Current
measurement and displays UNVERIFIED; remeasurement is a separate explicit action.
Cancellation or failure keeps the applied gain and the unverified state. Completed
measurement displays its actual result. No rollback, retry or hidden correction.

The sample order is original programme/composition, section gain, then global
trim. Every programme channel, including LFE, receives the same coefficient.
Composition inputs remain independently accessible and unmodified. M&E/dialogue
balance does not change. The implementation adds no lookahead, duration change,
resampling, playback buffering or presentation-clock change.

`section-gain-v1` is a typed Core plan, used by one standalone worker processor
for measurement, audition and export. It uses exact absolute frames and the
versioned centered 100 ms dB ramp specified in document 44. Each stage has its own
float32 rounding boundary. The original global-gain arithmetic is retained.
Equal adjacent gains canonicalize; names and redundant equal-gain boundaries are
outside audio identity. An all-zero envelope uses the legacy identity. Applying
an identical canonical plan preserves revisions and valid measurements.

Annotations remain independent. Geometry changes stale a dependent draft without
changing the applied envelope. Applied intervals retain their original annotation
provenance. After editing annotations, a new draft begins at zero unless its
geometry still matches the applied provenance; the UI explains this. Labels alone
do not stale a draft. Input/role/composition changes cannot rebind an old envelope.

The M3 continuation is versioned `global-trim-m12.1`. It requires valid measured
Source and Current on matching inputs, binds the fixed section-plan identity,
and proposes only the final global gain. A measured pass is a no-op. Predictions
remain separate from measurement and cannot satisfy export gates.

## Persistence, reports and product boundaries

Session schema 4 reads schemas 1–3, retains exact section/global doubles and
frame coordinates, and requires an M12 reader when saved. Save As preserves an
older-readable original. Unapplied drafts must be Applied or Discarded before
Save; leaving a draft offers Discard/Cancel. A draft created during a background
save remains unsaved, preventing a pending close/open action from dropping it.
Damaged/incompatible caches remain replaceable while valid core settings recover.
Reopening is stopped, device-closed, routes-muted, with a new runtime job identity.

Report schema `levelscope-measurement-report-v3` records the processing engine,
stage order, exact gains, applied frame intervals and transition intervals.
The retained scalar `processing_gain_db`/`processing_gain_hex` fields describe the
final global trim; `current_processing_plan` records the complete section-plus-
global plan. There is no single combined gain for a varying envelope.
Measurement-only reports still distinguish unavailable/failing results and invent
no artifact hash or verdict. Export still requires passing measured Current and
passing finalized, reopened output with no additional processing. One deterministic
TPDF stage occurs at final PCM24 quantization, never in Current processing.

New section processing is stereo/conventional 5.1 at 44.1/48/96/192 kHz, at most
three hours and 10,000 intervals. Mono/5.0/7.1 static-gain behavior remains supported
as before. Mono section processing and verified export are a near-term extension;
they are not M12 claims. Shared processing retains channel-driven twelve-channel
capacity. Physical full-7.1 listening remains separately pending.

Future component processing belongs before composition and must enter Current
identity; whole-programme processing remains a separate stage. Analysis-only
references stay non-audible with independent analysis invalidation. No component
DSP, dialogue analysis, automatic section proposals or intelligibility claims are
implemented here.

## Evidence and deferred work

Candidate-01 stopped at an old M5 test assumption that every editor control must
fit inside an 820-pixel panel. M12 deliberately uses scrollable content. The
corrected regression verifies control bounds and that the last controls can be
reached in the minimum window. Candidate-02 changed that test only, preserving
production source and app bytes.

Candidate-02 then stopped at the M8 migration check: the qualification packager's
audio-only derivative had lost the session's canonical JSON spacing. The reader
correctly rejected it. Candidate-03 changes only packaging/runner code to remove
the attachment without reserializing the rest of the core, and tests that rule.
All native executables match candidate-02. Both historical failures remain failures;
neither required changing production code or loading video.

Passing native checks are reused only after exact production-source, measurement-
build, executable-hash and actual-host equivalence checks. Receipts retain their
original candidate identities and reuse chain. Silent UI was observed on candidate-01
and reused on candidate-03 with the identical application. The new layout test and
unfinished or affected packaging/migration checks are run separately. Original
candidate packages and results remain unchanged.

Development evidence under `build/m12-qualification-2026-10-27` includes the
3,304-assertion final development audio matrix and an independent PCM24 decoder/quantizer check.
These are development results; formal frozen-candidate receipts are recorded separately.
The first sandboxed smoke attempt failed because JUCE could not create its cache
directory. The same audio-only test with normal cache access passed; the failed
attempt remains recorded. No physical device or video test was involved.

The required fresh selection covers M12, affected M1–M11 audio callers, pinned EBU
and ITU vectors, one consolidated three-hour composition workflow, independent
RIFF/RF64 verification and the silent UI checklist. It invokes explicit executable
modes rather than historical qualification runners. M8 migration fixtures with an
attached picture are copied to separately identified audio-only test documents;
the original is unchanged, and rational offsets remain serialized. This does not
qualify video migration.

All actual High Sierra work, cross-platform exchange and video-loading/decoding/
playback checks are deferred until after M13. Local Intel compilation with the
10.13 deployment target is build evidence only. M11 results retain their original
candidate identity and do not qualify M12 features. Occasional picture skips
remain an accepted development limitation; audio interruption, sustained freezes,
drift, identity or lifecycle failures remain blockers when video testing resumes.

The required Sequoia short, long and UI receipts passed the strict completion
validator. The phase is **“Sequoia tested; High Sierra and video qualification pending.”**
`Tests/complete_m12_sequoia.py` rejects missing or mismatched evidence and records
full qualification as false. Subsequent authorization permits only the development
checkpoint publication and an M13 proposal. Combined post-M13 qualification must
precede a request to integrate into main.

## Completed Sequoia results

The [evidence index](../../Tests/Reference/results/m12-sequoia-phase-2026-10-27/README.md)
and [completion receipt](../../Tests/Reference/results/m12-sequoia-phase-2026-10-27/completion-candidate-03.json)
identify every required run, original/reused receipt, build and artifact hash.
Actual host: Apple M1 Max, arm64, Sequoia 15.7.9 (24G830).

- All nineteen required short runs passed. M12 covers 3,304 assertions and 32
  rate/layout/input cases, including all-channel offline/audition/export parity,
  exact transitions, invalidation, failure recovery and session migration.
- The selected M1–M11 audio regressions passed, including mono and existing 5.1/7.1
  static-gain behavior. Pinned EBU v05 passed all 70 files; ITU BS.2217 passed all
  39 results and its documented baseline comparison. EBU and M12 runner checks pass.
- Two supplemental schema-2 stereo/5.1 migrations passed with the frozen reader,
  completing explicit coverage of old schemas 1–3. Exact settings and alignment
  confirmation survived; original bytes stayed unchanged. The first supplemental
  oracle incorrectly required pre-archive millisecond mtimes to remain identical.
  Its failure is retained; the corrected check permits verified filesystem metadata
  refresh without relaxing content identity or settings preservation.
- Eleven silent UI entries pass, with the exact-app reuse provenance described
  above. Independent verification checked all 384,000 PCM24 codes in the UI export.
  Devices stayed closed and routes muted; no physical output or video test ran.
- One consolidated three-hour, 48 kHz, six-channel composition with 1,000 applied
  intervals passed 2,067 checks in **70.0 minutes**, without repeating an earlier
  long workflow. Eight injected Source/Current captures matched the independent
  sample oracle at start, transitions, middle and EOF, with zero underruns and
  verified-source-copy reuse.
- Final Current measured -22.99999999067813 LUFS / -21.061581588646064 dBTP. The
  finalized RF64 artifact, measured at zero additional gain, measured
  -22.99999999081686 LUFS / -21.06157025973893 dBTP. Both **pass selected programme
  loudness and true-peak rules**; Source's failing loudness was not an export blocker.

| Sustained resource check | Observed | Limit |
| --- | --- | --- |
| Additional peak working RSS | 87.3 MiB | 256 MiB |
| Message-loop heartbeat p95 | 15.06 ms | 100 ms |
| Owned-copy cancellation | 42.60 ms | 1,000 ms |
| Save / parse | 95.59 / 463.47 ms | 5,000 ms each |
| Session size | 10,589,230 bytes | 128 MiB |
| JSON + HTML reports | 1,005,198 bytes | 16 MiB |
| Peak disk growth (informative) | 30,049,669,120 bytes | Preflight required at least 65 GiB free |
| Consolidated runtime | 4,199.545 seconds | 180 minutes |

The real RF64 file is 9,331,200,104 bytes, 518,400,000 frames, six channels, 48 kHz.
SHA-256: `a57c70b0a4955556712c80706446060964a2ba8419ddb1f7cc6c6c180afcbcd6`.
The native measurement scans the entire artifact. Independent libsndfile 1.2.2
checks decoded 16,384 selected frames; the additional known-input section/global/
single-TPDF oracle checked 188,304 PCM24 codes, including three transition regions.
These are bounded independent sample checks, not a claim that the independent
decoder compared every sample in the three-hour file. Complete file hashes and
JSON/HTML/manifest consistency also pass.

No required Sequoia check is outstanding. Actual High Sierra runtime, cross-Mac
exchange and all selected video checks remain pending until after M13. No physical
listening qualification is inferred from injected captures. Main, existing tags,
earlier packages, installations and the video worktree are preserved.

## Reproduction and review

The frozen package contains `identity.json`, the exact `source/` snapshot, the
standalone app, named test executables and Python runners. Verify every artifact
against `identity.json` before reuse. A rebuild is a new candidate, not a replacement
for the tested package. `build/logic-m12/CMakeCache.txt` and the identity record
retain the universal Release / macOS 10.13 / JUCE 7.0.12 build configuration.

From the frozen package directory, `Tests/run_m12_qualification.py --phase short`
requires a fresh `--output` directory and `--decoder-library` path. It enumerates
all nineteen required runs. The optional `--reuse-short` and `--reuse-package`
arguments require matching production source, build, executable hashes and host;
they preserve the original result identities. No historical milestone runner is
called recursively.

`--phase long` additionally requires the explicit immutable M11
`candidate-02-results/sequoia/long/long-files/composition.levelscope` input. Run it
only for an affected qualification need, not merely to repeat a passing sustained
workflow. Allow 65 GiB free headroom and up to 180 minutes. It uses injected devices
and audio-only captures; it never launches the historical composition/video runner.
The independent decoder remains the preserved M10 libsndfile 1.2.2 dependency,
identified by hash in the environment receipts.

`Tests/complete_m12_sequoia.py` takes the package, short and long result roots and
the observed silent UI receipt. Missing, failed or mismatched evidence prevents
phase completion. The script cannot declare High Sierra/video qualification or
authorize publication or main integration.
The final coverage audit additionally requires the two schema-2 migrations:
`Tests/Reference/results/m12-sequoia-phase-2026-10-27/auxiliary-tools/complete_with_migration.py`
calls that frozen validator and checks the supplemental reports, binary identity,
original hashes and exact settings before producing the final completion receipt.

For operator review, open the candidate app explicitly by its package path. Import
or reopen supported audio, confirm roles and initialize Current. Open Measured
diagnostics → Programme sections, edit annotations, create a section-gain draft,
set each desired gain and Apply while stopped. Re-measure Current explicitly.
If needed, use the existing final-global-trim proposal and explicitly apply and
remeasure again. Applied audio remains unchanged by later annotation edits. Save
requires resolving the draft; Render and verify / Publish retain the M10 gates.

## Follow-up order

Design revision, subsequently approved for implementation on 28 October 2026:
[47 — M13 proposal](47_M13BoundedSectionTrimProposals.md)
is the first bounded **section-led** selection method. The revised
[shared specification](../ProposalEngineSpec_v2.md) separates dynamics-led,
section-led and mixed strategies from their algorithms and processing preferences.
M13 still uses no limiter/dynamics and retains its conservative defaults; the
0.5 dB reserve is included inside the 3 dB cap, with partial reserve disclosed.
The bounded cinema evaluation and current M13 qualification status are recorded in
[49](49_M13QualificationAndHandover.md). No preserved M12 evidence or tested package
has changed; combined M13 results retain their own candidate identity.

After combined M12/M13 qualification, separately propose mono section/export support
and basic optional marker import; qualify standalone limiter, downward and upward
processing in small steps; then dynamics-led proposals, enhanced section-led
proposals and the mixed strategy. Keep future milestone numbers open. Marker/EDL
import uses explicit time origin/rate, exact audio frames and preview/confirmation;
imports are annotations, not gain changes, and useful proposals need no marker track.
Component-relative dialogue/M&E processing and validated intelligibility objectives
follow later. 7.1.4/Atmos and VST3 parity remain deferred. Main stays at M11 until
combined actual High Sierra/selected video qualification and subsequent integration
approval; the published M12 checkpoint and all historical records remain preserved.

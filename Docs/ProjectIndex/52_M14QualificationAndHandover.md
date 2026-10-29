# M14 implementation and qualification handover

Status, 29 October 2026: **Sequoia tested; High Sierra, video and three-hour checks pending.**
Candidate-02 is accepted for publication on `codex/logic-module` as the annotated
development checkpoint [m14-mono-rf64-sequoia-2026-10-29](https://github.com/carlozaiana/LevelScope/tree/m14-mono-rf64-sequoia-2026-10-29). No main integration is authorized. Main remains qualified M13
`a9be070ab9eb5cb4f591ff9d991d8999c3a016d0`.

Contract: [50](50_M14MonoProcessingAndVerifiedExport.md). Pre-implementation
impact selection: [51](51_M14ImpactQualificationPlan.md).

## Candidate and fresh evidence

Candidate-02 is the accepted development-checkpoint candidate under
`build/m14-qualification-2026-10-29/candidate-02`. The approved reduced Sequoia
phase is complete. This is not full M14 qualification.

- Base: `a9be070ab9eb5cb4f591ff9d991d8999c3a016d0` (M13 base; candidate implementation is identified by the manifest below).
- Candidate-02 source manifest: `9dd22a598281e83431039f210202a8a24fe49517294cc6bb8d83032c4975e013`.
- Original tested candidate-01 manifest: `d584a7a1b95a33b7f9befe75129d59fa7ccea04ec1f4b5624498616af4edacf9`.
- Measurement build: `4002f49d4d8eb3cffcf825f1624d8a00d6f48557f0755e711683cff6627f68da`.
- Universal app executable SHA-256: `188446f5752d91e40e90066d2b25133eaaa07b421ed0018da15a1be2fb60777c`.
- Candidate-02 identity-file SHA-256: `c8b0184b700be7ada4963b0b80bcfca1e577936855291e1280e5f6aa0ccc913d`.
- Candidate-01 identity-file SHA-256: `50b81b6be16000ed70304b3f48b8abab18a06df6015e1df75b11a217f1749655`.
- Actual execution: Apple M1 Max, arm64, macOS Sequoia 15.7.9 (24G830).
- Universal arm64/x86_64 app and named tests; Intel load command minimum 10.13,
  arm64 minimum 11.0, SDK 15.5, existing JUCE 7.0.12. Intel is cross-built only.

The manifest contains source/tests/build configuration, not qualification prose.
Frozen packages are never rewritten to add later reports. Results live separately
in `sequoia-candidate-01/{short,cinema,long}` and the UI checklist. Each phase
validates frozen artifacts and records commands, log/report hashes and environment.
`complete_m14_qualification.py` requires every phase and all silent UI observations;
even success explicitly leaves full qualification false.

Candidate-02 repairs only three qualification Python files: the runner, completion
gate and runner tests. Every production/C++ test source and every standalone,
test and legacy binary is byte-identical to candidate-01. The explicit equivalence
record permits retaining candidate-01 short/cinema/UI and the single sustained
workflow with their original identities. Fresh candidate-02 evidence is in the
`legacy` phase and `candidate-02-runner-tests.log`.

The first legacy-helper invocation was **inconclusive**: its grouped-input
assertion could also reject an ordinary interleaved session. That historical receipt
is preserved, not used as proof of parser rejection. The corrected test invokes
the old M13 parser and requires both exact messages: successful core read for the
mono baseline, and `Applied section plan does not match this programme` for the
applied mono sections. The old schema-1 helper intentionally exits nonzero in both
cases; its parser messages distinguish readable core from rejected processing.
Completion requires this correction; no production rebuild or repeated long run
is justified by a qualification-tool-only change.

Fresh short coverage passed 22 selected run groups: 1,666 M14 assertions over
12 rate/role workflows; 48 raw odd/even PCM cases; focused M12/M13 and composition
tests; measurement/scaffold/compliance/proposal/diagnostics/section/audition/7.1/
session/grouped-input/export regressions; 70 pinned EBU and 39 ITU files; and
qualification-runner checks. The original ambiguous legacy row is superseded by
the separately required corrected phase. Independent libsndfile 1.2.2/Python checks
cover every sample of all short exports (11,403,012 mono PCM codes), plus 2,472 raw
quantizer codes. Its retained library SHA-256 is
`e9fd4d76bc457dcfb956617e8d467df85087ce069747fccd78e1e702d8629c14`.

All 11 silent Sequoia UI cases passed with the candidate-01 app: import/roles,
measurement/profile, manual draft/Apply, independent actual failure, proposal/
Apply/verified pass, already-pass preservation, reader warning/Save As/reopen,
explicit verified export, RF64 notices, invalid-input preservation, and stopped/
muted/device-closed lifecycle. The exported 480,001-frame mono PCM file has
1,440,072 bytes including its single non-audio pad byte and SHA-256
`99c09351c6f047615e697f9955f9cb585dff2d0c09c997c976ea4ac0cc112690`.
Independent decoding and JSON/HTML/manifest checks pass. The test app was closed
without opening any physical audio device; no video was loaded.

## Sustained results and resources

The single 90-minute mono/192 kHz workflow passed 1,065 assertions. It covered
1,036,800,000 frames, 1,000 annotations, three private proposal passes, explicit
Apply and independent Current measurement, diagnostics, session save/reopen,
cancellation, eight injected Source/Current capture/seek cases and verified PCM24
export. There was no second long run after the tool correction.

| Observation | Result | Approved limit |
| --- | --- | --- |
| Complete native workflow | 2,109.37 s / 35.16 min | 120 min |
| Additional peak process RSS | 117.5 MiB | 256 MiB |
| Message-loop heartbeat p95 | 15.03 ms | 100 ms |
| Owned-copy cancellation response | 11.30 ms | 1 s |
| Session save / core parse | 53.87 / 234.69 ms | 5 s each |
| Complete media re-verification on reopen | 18.51 s | Recorded; distinct from core parse |
| Session size | 5,287,294 bytes | 128 MiB |
| JSON + HTML size | 357,921 bytes | 16 MiB |
| Observed filesystem free-space reduction | 10,580,905,984 bytes / 9.85 GiB | >=20 GiB initial headroom |

The filesystem delta includes concurrently written candidate-02/qualification
artifacts; it is a conservative observed free-space change, not isolated audio
allocation. Source, Current and private measurement storage is included in the RSS
growth. The separate 10,000-section retained-plan/attribution bounds also passed
in the selected M12/M13 regression checks.

Actual Current and finalized artifact passed the selected rules; the private
candidate's pass was never substituted for either. The final output remains one
Centre channel at 192 kHz, ordinary RIFF, 3,110,400,068 bytes. SHA-256:
`94070c5af19fd502f417dcb9bb7bd69dcda07df8a03f87d1ff4c1c5ff79be240`.
The qualified worker measured every output frame at zero additional gain.
Independent libsndfile checks covered 16,384 decoded frames; the independent
section/global/single-TPDF oracle checked 60,480 codes at bounded locations,
including transitions and the unique peak. Complete file/hash/report checks passed.
These independent long checks are bounded, not a second full decode or a
three-hour/RF64 end-to-end claim. All 48 legacy stereo/5.1 quantizer cases also
passed the changed independent parser, with 196,800 exact codes.

## Evidence, reproduction and preservation

The [text-only evidence archive](../../Tests/Reference/results/m14-qualification-2026-10-29/README.md)
includes frozen identities, exact commands/receipts, equivalence, build inspection,
all completed reports and the completion gate. Raw media, binaries, captures and
connection details are excluded. Preserved app:
`build/m14-qualification-2026-10-29/candidate-02/LevelScope Standalone.app`;
the actually exercised candidate-01 app is byte-identical and remains preserved.

`Tests/package_m14_qualification.py` freezes only the standalone and 15 named test
targets from the universal Release build, pinned reference files and source manifest.
`Tests/run_m14_qualification.py` enumerates short/cinema/long/legacy phases and their
required tests; it does not call a historical milestone runner. Reproduction requires
the local pinned vectors, the explicit original/copy cinema paths and the retained
independent decoder library. Exact invocations are in each phase's `runs.json`.
Use the candidate-02 completion tool with `--candidate candidate-02`,
`--tested-candidate candidate-01` and `--results sequoia-candidate-01`; it validates
the approved tool-only reuse and all original hashes. Its result is
`sequoia_phase_complete: true`, `qualified: false` with explicit deferrals.

At qualification completion, main and the logic branch still pointed to M13 and
all **17 existing annotated tag objects and target commits** matched the before
record. The subsequent publication approval advances only `codex/logic-module`
and adds the new annotated development tag. The [publication record](../../Tests/Reference/results/m14-publication-2026-10-29/README.md)
connects this commit/tag to the frozen manifest and tested binaries. Publication
rechecks hashes and completion receipts without rebuilding or rerunning devices.
Original qualification receipts retain their pre-approval status fields; those
fields are historical, not the current publication authorization. No required
Sequoia check remains outstanding. No main integration, video working-copy change,
installation/plugin change, Golem access, physical-device opening or system-audio/
Pro40 change is part of publication.

## Implemented behavior

Confirmed one-channel Left/Right/Centre uses the unchanged section/global sample
evaluator, M13 policy and measured Apply/verify lifecycle. Mono output has one
explicit route/mapping and preserves its speaker mask. No duplication, compensation,
new target rule or new arithmetic is introduced. Odd PCM24 payloads have one pad
byte outside the audio frames. Export retains one final TPDF quantization stage,
Current/artifact measured gates and strict canonical artifact validation.

RF64 input has a separate bounded compatibility policy. Zero-filled extended ds64,
zero sample count and structurally valid trailing metadata are handled as specified
in document 50. Compatibility notices are visible and retained in provenance.
The strict export inspector does not accept those input exceptions. Mono section
sessions retain schema 4 and warn that an M14-capable reader is required; Save As
preserves the earlier file. The old-reader test uses the frozen M13 binary.

Direct import of the original cinema file passed on candidate-01: exact original
SHA-256 `1fb7f93a33d600bf0a62d63e642e5703f03aec8e6920f8a6e510cf65f61fdeaf`,
269,875,200 frames at 48 kHz (5,622.4 seconds), PCM24, six ordered L/R/C/LFE/Ls/Rs
channels. Six 4,096-frame decoded comparisons against the verified M13 copy agree
exactly, including crossing and beyond 4 GiB and the final span. Original size and
modification time remain unchanged. This is bounded decoded equivalence with full
content identity, not a new complete cinema proposal/export workflow.

## Reused evidence

M13 provenance and original cinema payload-identical-copy records are available;
reuse is conditional on documented unchanged dependencies and fresh original/copy
identity/parity checks. M13 runtime evidence is not new-feature qualification.

## Pending coverage and review boundary

- Actual Intel High Sierra execution/UI and cross-Mac session exchange.
- All video load/decode/playback; rational timing data may be checked without media.
- New three-hour mono/RF64 full-duration workflow. M14's fresh sustained workflow
  is limited to 90 minutes at 192 kHz and ordinary RIFF export; sparse/container
  checks beyond 4 GiB are not full-duration end-to-end evidence.

Plan combined qualification after about five post-M13 milestones (through M18),
no later than six (M19), with one frozen combined candidate and impact-selected
manifest. Future milestones need their own approval. Explain concrete dependency,
arithmetic, clock/buffering/latency or failure triggers before earlier deferred
testing. Accepted occasional picture skips retain existing audio/drift/lifecycle
requirements; full physical 7.1 listening remains separately pending.

Readiness review is complete and the user authorized development-checkpoint
publication on 29 October 2026. Stop after publication; main integration and the
next milestone remain unauthorized. Pending qualification is unchanged.

# M16 qualification and handover

30 October 2026. **Candidate-03: Sequoia tested; High Sierra, video and three-hour checks pending.** The approved Sequoia phase is complete. Candidate-03 is published on `codex/logic-module` at commit `73a7b25e71f87299fa2d93c1b1106c118c55beec`, with annotated tag [`m16-proposal-preferences-sequoia-2026-10-30`](https://github.com/carlozaiana/LevelScope/tree/m16-proposal-preferences-sequoia-2026-10-30). Main remains qualified M13 `a9be070ab9eb5cb4f591ff9d991d8999c3a016d0`; no main integration is authorized. The separate [publication receipt](../../Tests/Reference/results/m16-publication-2026-10-30/README.md) records the committed-source mapping and remote verification.

## Implemented boundary

ProjectIndex 56 records the revised approval; 57 records the impact-selected plan, and 58 records the LRA-focused roadmap/UI priorities without implementing them.

`ProposalPolicy` is immutable per request, separate from measurement profiles and audio processing. Defaults: no region-count restriction, newly affected duration 20%, additional attenuation 3 dB, cumulative local attenuation 12 dB, reserve 0.5 dB inside both gain caps. Optional explicit count is 1–10,000 within the independently enforced interval capacity. Duration supports 0–100%, additional 0–6 dB, cumulative 0–12 dB and reserve 0–1 dB (no larger than the additional cap), all exact tenths in the document. These are preservation preferences, not EBU or universal cinema rules. No silent clamp or subset selection.

The method compares canonical envelope plateau/ramp segments, joins touching changed support, includes transition endpoints conservatively, and clips to the existing programme. Equivalent canonical processing ignores redundant marker splits. Different boundaries can change attribution/processing legitimately. It reports changed annotation rows separately from connected processing regions, full correction details, canonical complexity, additional/cumulative reduction, cumulative support and neighbour/slope descriptors. A zero-headroom correction remains explained without being counted as a changed annotation. No coordinated neighbour adjustment or automatic segmentation is added.

Already-pass remains a genuine no-op and global correction remains first. There is one proposed envelope and at most three private passes. All-channel processing, full gain precision and independent Current remeasurement are unchanged. No delivery threshold is relaxed; 20% still rejects the historical cinema-like 89.53% support. Preferences alone preserve valid measurements, evaluation and audio/render identity, while binding both global and section proposals to the new policy digest. Old completions/Apply are rejected. Report approval includes metadata identity.

The explicit limits panel edits a draft and confirms it only while stopped/no file or measurement/proposal job. Cancel and identical confirmation preserve active state. Complete proposal inspection uses bounded pages with scrolling, without the old 32-row cutoff; final outcome leads the summary. Reports expose complete transient proposal inspection without treating private measurement/planning as Current authority. Reports preserve active policy and separate historical applied policy/binding. Sessions retain exact applied envelope and historical policy/binding; transient private results/explanation are not restored as authority. Save a report before discarding/closing if the complete original proposal explanation is needed later.

Session schema 6 stores policy. Pre-M16 schemas 1–5 automatically adopt defaults in memory while retaining media, exact processing and annotations; no legacy executable mode/adoption prompt. Opening never rewrites the original; explicit Save/Save As requires an M16 reader. Invalid/unknown independently identifiable policy recovers the otherwise valid core with proposals unavailable until explicit replacement. Structural core corruption still fails. Conservative build/platform caches remain replaceable. Audio remeasurement may be needed across builds.

## Exact tested identity

- Frozen package: `build/m16-qualification-2026-10-30/candidate-03/LevelScope Standalone.app`.
- Source/test/build manifest SHA-256: `3ac2e19d8dd976a8fdc41490245ad9961af91e4bfe1a1d5afdc1ad5c03897f38`.
- Measurement/build identity: `9b7d85d5463297b70e6edaabd5010ddabbd317ea63e2ae1eea7a8280da387bcf`.
- Universal app executable SHA-256: `8ee1b747685a8840284dbe6992b95c1076a76fc90cfe510ac043b63c61df2286`.
- Published checkpoint commit: `73a7b25e71f87299fa2d93c1b1106c118c55beec`. All 269 committed source/test/build files match the frozen manifest. The original 283 package artifacts and 88 indexed qualification receipts were hash-verified unchanged; publication required no rebuild or repeated runtime checks. The publication documentation follow-up does not change tested code.
- Development base HEAD: `63eafac9acc210d762dfd5e3dcc70defbe5bbdf3` on `codex/logic-module`, following the published M15 checkpoint. This is a base commit, **not** a commit identifying the M16 implementation. The frozen manifest is the tested source authority. Qualification prose and textual receipts are recorded after freezing and are outside that source digest.
- Actual execution: arm64 Apple M1 Max, macOS Sequoia 15.7.9 / build 24G830. Release universal arm64+x86_64; every required binary contains the Intel slice with minimum macOS 10.13. Established Apple clang/JUCE 7.0.12 dependency identities and individual binary hashes are in the receipts. Intel cross-build success does not qualify High Sierra runtime.

The preserved `identity.json` freeze-time pending status is not rewritten. Final status is recorded separately in [candidate-03/status.json](../../Tests/Reference/results/m16-qualification-2026-10-30/candidate-03/status.json).

## Fresh qualification

All **ten explicitly selected commands** passed their expected outcomes in **108.11 seconds**. The M15 reader command must fail for the valid schema 6 input; it is an expected negative test. No recursive historical runner, long media workflow, video execution, Golem access or physical-device opening occurred.

| Required check | Result |
| --- | --- |
| M16 preferences, support/accounting, capacity/full inspection, UI responsiveness, identity/session/report, large measured proposal, mono roles | PASS: all seven required groups, 256 assertions |
| Global proposals | PASS: 302 checks |
| Short session recovery/coherence/compatibility | PASS: 1,006 checks |
| Programme sections | PASS: 2,468 checks; synthetic duration coordinates only |
| M12 short processing/offline/pre-monitor/export/timing regression | PASS: 334 checks |
| M13 short proposal/apply/independent verification regression | PASS: 366 checks |
| Marker mapping, import/session migration and resource checks | PASS: 487 checks, using approved minimal Nuendo/Pro Tools fixtures |
| Short export/report/gate regression | PASS: 4,872 checks |
| Frozen M15 reader rejects schema 6 | PASS: expected failure; exact old build/message plus existing input/core checksum/schema/applied-plan/cache proof |
| Qualification runner self-tests | PASS: six tests; missing/duplicate/failed/mismatched acceptance evidence blocks completion |

The large measured synthetic case has **134 changed annotation rows**, **one connected processing region**, and **66,720 of 480,000 source frames**: **1.390000 seconds / 13.900000%**. Its exact private candidate passes in three private passes; Current stays unchanged until Apply. Separate tests prove explicit 65-region acceptance, rejection under a 64-region user limit, the 10,000-interval boundary, density invariance for identical canonical processing, and the cinema-like 89.53% support rejection under the default 20%. No full cinema rerun or file-specific tuning was performed. Confirmed mono Left/Right/Centre each pass explicit application and independent Current remeasurement without duplication.

Silent native UI checks were repeated with **candidate-03 itself**: safe stopped/device-closed/muted reopen; measured proposal generation without application; all four inspection pages including Section 134 and final private verdict; accepted support in frames/seconds/percentage; cancellation; invalid reserve rejection; explicit 65/100% settings; stale Apply disabled while measurements and 0 dB applied gain remain current/unchanged; SaveAs/New/reopen with exact settings; explicit unrestricted-default restoration. The superseded empty global-draft notice is absent after section proposal completion. The test app was closed after checks. The saved test session/media remain outside Git.

### Resources and boundaries

| Short/structural observation | Observed | Approved bound |
| --- | --- | --- |
| 10,000-interval support selection |43.17 ms;8.97 MiB added peak RSS |500 ms;32 MiB |
|Retained attribution/plan/full text |8.18 MiB |16 MiB |
|Inspection workload added peak RSS |54.75 MiB |256 MiB |
|Complete 10,000-row text |22.86 ms |1 second |
|Complete JSON+HTML |5.38 MiB;47.47 ms |16 MiB;5 seconds |
|Inspection heartbeat |p95 48.65 ms; maximum 91.49 ms |p95≤100 ms |
|Cancelled request acknowledgment |0.71 ms |≤1 second; active private-job cancellation also passes M13 regression |
|Native UI 1,000-annotation session SaveAs/reopen |464 ms / 578 ms including UI observation;56,644 bytes |5 seconds;128 MiB file cap |
|Disposable candidate-03 test data |203.68 MiB |2 GiB |

These are short synthetic/structural observations. Maximum-cache session storage and full three-hour audio were **not** rebenchmarked; unchanged core/cache resource guards remain, with targeted session regressions fresh. Synthetic 10,800-second diagnostic coordinates do not constitute decoding/processing three hours of audio. Full-duration and actual High Sierra resources remain pending.

### Historical results and justified reuse

Candidate-01's runner stopped on an obsolete wording assertion (`above 20%` vs`above 20.0%`); its failure remains recorded. Candidate-02 corrected that test and passed all selected automated/UI checks. Final display review then added accepted support seconds/percentage and removed an irrelevant empty global-draft notice. Candidate-03 also strengthens the reader-evidence gate; the full short manifest/UI was rerun. Both earlier frozen packages and original receipts remain unchanged. Development failures and fixture corrections remain in the evidence index; no prior failure was relabelled as a pass.

The audit confirms measurement arithmetic, sample processing, programme decoders/input/composition, PCM24/RF64 implementation, audition, clocks and video code remain byte-identical to the appropriate earlier candidates. All 24 Core/Measurement and Core/Processing files also match qualified M13. The only CMake change adds the new preferences panel to the existing UI library; dependency/optimization/deployment configuration is unchanged. Shared session/proposal/report callers were freshly exercised despite lower-level reuse.

- M13 candidate-02 manifest `50eba19a2d66096607ec878c8eaa05e610580c7f5e8854d7b375d846a89c5f27`: original pinned EBU/ITU and audio/video baseline evidence remains attributed to M13, not M16.
- M14 candidate-02 manifest `9dd22a598281e83431039f210202a8a24fe49517294cc6bb8d83032c4975e013`: unchanged mono/RF64 arithmetic/decoder and prior bounded cinema/90-minute evidence, retaining its candidate-01 binary-equivalence provenance.
- M15 candidate-05 manifest `38676a45d47ddba6eb54f7ab299f52f622f6be26554256781b8960749c609564`: actual DAW export-setting provenance and unchanged marker interpretation; focused import/migration tests are fresh.

These reused records do not qualify M14–M16 on High Sierra or qualify new features against video/full-duration requirements. No dependency or observed failure triggered earlier deferred testing.

## Reproduction and preservation

Use `Tests/package_m16_qualification.py` to freeze a **new** directory after building only the eight named app/test targets listed in its manifest. Never overwrite a frozen candidate. To rerun the accepted selection into a new disposable directory:

```sh
python3 -B build/m16-qualification-2026-10-30/candidate-03/Tests/run_m16_qualification.py --candidate build/m16-qualification-2026-10-30/candidate-03 --output build/m16-qualification-2026-10-30/new-sequoia-run
```

This runner intentionally leaves full qualification false until the separate silent UI and integrity/phase records are combined. It verifies candidate hashes before execution and requires actual arm64 macOS. Native tests require access to their owned temporary verified-source copies. No device approval is implied by this command.

Full textual evidence: [M16 evidence index](../../Tests/Reference/results/m16-qualification-2026-10-30/README.md), including candidate-03 identity, command receipts, UI observations, resources, reuse audit and preservation record. At qualification, all 19 existing annotated tag objects and targets, the then-current logic HEAD, main M13 and registered video worktree HEAD/location were unchanged; unrelated Finder metadata was preserved. M13/M14/M15 frozen artifacts were hash-verified. No installs, plugin builds, system audio changes or Pro40 access occurred.

## Deferred and remaining scope

Actual Intel High Sierra/cross-Mac exchange, selected video and three-hour workflows remain deferred to a frozen combined candidate after M18, no later than M19. A universal Intel 10.13 cross-build is not runtime qualification. Occasional picture skips remain accepted when those tests resume; audio interruptions, drift, sustained freezes and identity/lifecycle failures remain blockers. Physical full-7.1 listening remains separately pending.

Selective LRA requirements/presets, progressive automatic analysis, interactive waveform/measurement timeline, region suggestions, selective raising, protection semantics, coordinated neighbours, limiting/compression and alternative strategies require separately reviewed milestones. M16 adds no LRA acceptance rule or dialogue/intelligibility claim. Main and all existing tags/packages/installations remain preserved. Checkpoint publication is complete; main integration remains deferred and requires later review/approval. All 19 pre-existing tag objects and targets and main/video refs were verified unchanged directly on GitHub. See the separate publication receipt for exact identities. [60 — M17 job LRA requirements and presets](60_M17JobLraRequirementsAndPresets.md) is a proposal only, awaiting approval; no M17 implementation is included.

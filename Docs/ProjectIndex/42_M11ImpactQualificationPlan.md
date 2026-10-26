# M11 pre-implementation change-impact and qualification plan

Recorded 26 October 2026 **before implementation**. Base/source publication is M10
`fcd4c15fe78d6352f5e96faabf3eb52f942235e1`; approved scope is [41](41_M11AlignedStemComposition.md).
Historical status when recorded: plan recorded; no M11 candidate or qualification
pass yet. Current results are maintained separately in [43](43_M11QualificationAndHandover.md).

## Expected changes and required fresh evidence

| Affected boundary and indirect consumers | Required frozen-candidate checks on both actual Macs |
| --- | --- |
| Programme definition/reader and verified member copies; called by measurement, audition and export | M11 three-presets × four-rates × supported component representations; independent sum oracle, chunk/seek parity, non-finite/over-full-scale input, all-or-nothing failures, disk/cancel/reuse; pinned EBU/ITU through current adapter |
| ServiceContext and immutable Source/Current identities; compliance, M3, diagnostics and approvals depend on them | Existing focused compliance/proposal unit binaries plus M11 mismatch/invalidation/no-op/failed-verification cases; focused diagnostics/sections equivalence to composed oracle |
| Shared audition preparation/reader, without planned callback/video-scheduler change | Composed pre-monitor capture at all rates/gains/channels; existing overload/restart/cancel/lifecycle tests; representative H.264 and ProRes Source/Current seek/offset/restart/EOF integration using injected devices |
| Session schema/cache/context serialization and UI; old interleaved/group callers | Schema 1/2 migration and schema-3 round trips; frozen-M10 reader rejection; missing/identical/different member relink; corruption/recoverable cache, save consistency, stale jobs; bidirectional Mac exchange |
| Export member checks, artifact context, report provenance/approval | Short actual composed PCM24 exports with independent decoder/seeded oracle, exact artifact/report hashes and strict gates; existing short export numeric/transaction/failure tests invoked directly, no historical runner |
| Standalone UI and CMake source/build identity | Silent UI workflows on both Macs; explicitly named standalone/test targets and final source/binary manifest; existing interleaved stereo/5.1/System I 7.1 and M9 groups compatibility |

Development uses short individual tests by area. Frozen selection is a literal
required-case inventory in the M11 runner/aggregate validator: nonzero exit, missing
case, wrong candidate/platform, skipped required phase or missing UI/cross-Mac
receipt prevents qualification. Do not call `run_m8_qualification.py`,
`run_m9_qualification.py` or `run_m10_qualification.py` recursively. Existing small
unit binaries and runner/parser tests can be invoked directly where affected.
Record exact executed commands, identities, runtimes and results. Refine selection
as actual edits reveal dependencies, explaining additions/removals before freeze.

## Reused evidence, with limits

Reuse M10 candidate-03 source digest
`4ff89967e4927600d0b6defafabb28c372f4068870a0f1a9adb1ee84fe1781d8`, app executable
`92c013999bdecaf904af0888ce68dd59bcd1b3004934b9ddb36ce1a1bf511e15`, and its indexed
[results](../../Tests/Reference/results/m10-qualification-2026-10-25/README.md)
only for unchanged RF64/container/quantizer specification, prior broad codec/colour
coverage and historical physical endurance. They remain M10 evidence, not M11 runs.
Independent M11 output decoding and fresh focused numerical/transaction regressions
are still required because shared callers/identities change.

Preserve earlier measurement/profile algorithm qualification and original package
identities; fresh pinned EBU/ITU checks establish the changed adapter boundary.
Reuse full video codec/colour and physical endurance evidence only if the final diff
confirms unchanged decoder/scheduler/RT buffering/timing behaviour. Changed buffering
or timing paths, integration failures or uncertain dependencies trigger appropriate
additional tests and recorded rationale, with physical-device coordination first.
Ordinary worker-side composition is new sustained preparation work: injected capture
and the consolidated long workflow test it directly. No automatic device opening.

## One consolidated long run per Mac

One real three-hour 48 kHz 5.1+5.1 composition, with 1,000 sections: import/verification,
Source/Current measurement, first/middle/final injected captures, session save/reopen
and resources, PCM24 RF64 render/reopen/measure/publish and independent decoding.
Do not run earlier M4/M8/M9/M10 three-hour workflows separately. Use algebraic/oracle
samples without allocating another full-programme floating reference file.

Maintain approved budgets: <=180 min automated workflow, <=256 MiB additional
composition/export working RSS above existing result storage, UI p95 <=100 ms,
cancel acknowledgement <=1 s, save/core-cache parse <=5 s, session <=128 MiB,
JSON+HTML <=16 MiB. Report actual stages and disk high-water marks. Tests of
three-hour/192 kHz size arithmetic do not claim that sustained workload was run.

## Estimates and available resources

Planning estimates, not measured M11 results:

- Short development loops: seconds to a few minutes, only affected targets.
- Frozen targeted offline/regression/reference block: approximately 20–45 min
  per Mac; UI/exchange an additional 20–40 min, depending on manual interaction.
- Consolidated long workflow: approximately 60–100 min Sequoia and 95–160 min
  High Sierra, based on M10's 45/84 min workflow plus second-component I/O/summing.
  The acceptance budget remains 180 min; measure rather than assume compliance.
- Two PCM24 5.1 long inputs: ~18.7 GB; private verified copies ~18.7 GB; one RF64
  output ~9.3 GB. Allow about 65 GB working headroom including short fixtures,
  package/evidence and reserve, with no duplicate full-length composition render.
  Inspect actual filesystem free space before starting; share unchanged pinned
  media by validated references rather than duplicate historical packages.

Read-only preflight: Sequoia has ~126.9 GB available; Golem internal has ~193.0 GB
and RazzFazz ~1.24 TB. Golem is actual Intel High Sierra 10.13.6 (17G66), reachable
over the existing authenticated SSH connection. Inspect explicit Python/decoder
paths before qualification because noninteractive PATH does not find python3.
The relocated ProRes folder remains accessed through its verified original-path
symlink; do not modify that package or its M5 dependencies.

## Extension boundary and completion

Keep original components accessible. Version component processing, composition and
whole-programme processing independently; M11 accepts only identity component plans
and the existing whole-programme static gain. Future component DSP belongs to
Current, never original Source. Analysis-only identities stay non-audible and
separate. Future dialogue targets need explicit methods/criteria and measured
verification of both dialogue improvement and complete-programme consequences.

Stabilize Sequoia, freeze the universal package, qualify the same candidate on actual
High Sierra. Preserve every frozen package and historical failure. After changes
rerun affected checks and broaden only where justified. Keep occasional picture
skips recorded/accepted; audio interruption, drift, sustained freezes and identity/
lifecycle/memory failures block completion. Physical full-multichannel listening
remains separate. Final report must distinguish fresh, reused and outstanding
evidence. Stop before checkpoint publication or merge for user readiness review.

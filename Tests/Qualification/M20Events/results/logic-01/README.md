# M20 complete-event Logic interface evidence

19 November 2026. **26 required groups / 3,919 assertions passed on native arm64
Sequoia. Not combined M20 qualification.** No physical devices or video were used.
Actual High Sierra, new UI controls/fallback/native responsiveness, broader detector
controls, independent holdouts and coordinated listening remain pending.

See [implementation and outcomes](../../../../../Docs/ProjectIndex/81_M20EventImplementationAndImpact.md),
[UI interface handover](../../../../../Docs/ProjectIndex/82_M20CompleteEventLogicHandover.md)
and the [local checkpoint receipt](checkpoint.md).

- Final package: `/Users/fuzz/Developer/LevelScope/build/m20-events-2026-11-19/logic-handover-03`.
- Source/test manifest: `131ddde3204d54da6cfe63645004332a358ab9ce26e0561958ad2be76e6c8596`.
- Production build: `12d44f32d59488dcc03c9d88edea96d4c03a98b439fd9cd37d8d15078276c6a7`.
- App executable: `b20070886b34336ed971db2f404cbe19d362c417aef63e808f4a29c651fd7fd3`.
- Source base: `dcdff7ddbbe91abb1a64048f4c3486d88807537d`, containing accepted UI
  continuation `607a3f54a24474ebcb6bf7cb5daa59a4e64e4bf7`.
- Host: Apple M1 Max / 32 GiB, macOS 15.7.9 build 24G830, arm64 execution.
  Universal binaries; Intel minimum OS 10.13 verified in the load command.

`summary.json` is the explicit acceptance-checker output. `receipts/` contains each
selected receipt without rewriting its original build identity or measured values.
`execution-record.json` maps each receipt to its actual frozen package and records
RSS and byte-equivalence. `packages/` preserves all three source/artifact manifests
and the initial core receipt. Initial manifest/runs, supplemental runs and final core
command are separate records, with `/usr/bin/time -l` logs. `reports/` contains the
full textual proposal evidence, exact gains, measured/forecast differences, coverage
reasons, limits and omissions. Raw media, sessions, exports and binaries are excluded.

Package 02 added only fallback/population tests and required manifest entries.
Package 03 added five core assertions for old/shortened-ramp coverage. All production
sources, build identities and complete app contents are byte-identical across 01–03.
The affected core group was rerun on package 03. The other receipts are deliberately
reused for unchanged tested behavior, with their original identity preserved; no
later rebuild replaces the original apps or earlier evidence.

The real music cases in event mode return **Review required**, not an automatic
complete-event recommendation. The gate and cinema examples return measured failures
within the bounded search. A passing test receipt verifies the expected safe behavior;
it does not mean that content meets its selected requirements or sounds acceptable.
The source event detector has not reliably identified the user's first musical hit.

`preservation-before.json` and `integration-base.json` retain the starting refs and
unrelated-file hashes. Final preservation verification and the implementation commit
relationship are in the separately committed checkpoint receipt. No push, tag or main
merge is part of this Logic handover.

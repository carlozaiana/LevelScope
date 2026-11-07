# Timeline inspection evidence — 7 November 2026

[Complete handover](../../../../Docs/UIInspectionConvenience-2026-11-07.md).

`acceptance.json` records 1,505 focused checks per actual Mac, the exact candidate
identities, app equivalence across the test-only recovery revision, 65 music
reproducibility checks, and preservation of the existing session files/main.

- `sequoia-01`, `high-sierra-01`: unmodified test receipts/logs and rendered PNGs.
- `identities`: immutable package manifests; candidate 02's app is byte-identical
  to candidate 01 tested on both platforms. Only two test files changed.
- `listening-recovery`: complete new fixed-protocol music outcomes, plans and reasons.
  The stored original comparison evidence remains in the prior M20 evidence folder.
- `development`: clean expected preflight rejection of the changed comparison file;
  the earlier refresh harness's empty crash log is retained. `acceptance.json`
  records the null-history test failure and its cause. It was not an app crash.
- `evidence-sha256.json`: byte-exact inventory. Generated audio, binary packages and
  `.levelscope` listening copies remain in the owned build tree, outside Git.

No physical output was opened. No M20 publication or main integration occurred.
Listening remains pending; numerical equivalence is not an audio-quality judgment.

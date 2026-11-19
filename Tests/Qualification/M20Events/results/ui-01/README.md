# Complete-event UI integration evidence — 20 November 2026

Local tested implementation `8b28348dee77179b395253728d71818404bd5022` on
`codex/m20-complete-event-ui-2026-11-20`, integrating Logic tip
`8e50899d5487cef0466a8bd7ba30559047f972be` and accepted UI tip `607a3f5`.
Main remains `22961f03de7919e1de631f8b13eb4b852108ac0b`. No publication/tag/merge.

The exact app, manual checklist, bounds and remaining gaps are in
[ProjectIndex 83](../../../../../Docs/ProjectIndex/83_M20CompleteEventUIIntegration.md).
The final candidate is `build/m20-complete-event-ui-2026-11-20/candidate-01` in the
existing UI worktree. Files in this directory are receipt copies; original
candidate, gate, full synthetic fixture runs and earlier packages remain preserved.

- `candidate-identity.json`: full source/artifact manifest, exact commit/build,
  universal architectures and Intel minimum 10.13. This is not actual High Sierra
  execution of the new method.
- `summary.json`, six logs and their JSON receipts: finite final Sequoia programme,
  **2,767 passing assertions**. Native event 144, legacy UI 163, timeline 2,299,
  event stereo 55, event 5.1 52, explicit fallback 54. Native images are under
  `event-ui-files`; complete service reports are under `reports`.
- `authority-gate-identity.json` / `authority-gate-result.json`: preserved
  pre-enable product-off build with 141 passing native assertions. A test-only
  friend override exercised the real main-window callbacks before enabling the
  product action. The final candidate adds immediate status-refresh, redundant
  reselection and enabled-default assertions; all earlier authority checks remain.
- `initial-selector-failure.json`: initial failed check exposed JUCE changing the
  automatic item label without refreshing its selected text/ID. Fixed by restoring
  that ID with no change notification; the successful gate then preceded enablement.
- `acceptance.json`: `check_ui.py` verified every required receipt, exact build,
  hashes, all per-check results and the preserved gate. It deliberately reports
  `combined_m20_qualification_complete: false`.
- `impact-reuse.json`: 73 unchanged Core/service/test source hashes from frozen
  Logic package 03 and the rationale for reusing unchanged qualification. Current
  UI callers and representative combined-runtime service paths were tested fresh.
- `preservation.json`: main/origin-main, all 25 existing tag objects/targets,
  other worktree heads, accepted design sources and original Logic packages 01–03
  plus the accepted UI package verified unchanged after testing.
- `execution-record.json`: exact local paths, host OS/architecture and scope.
- `text-copy-hashes.json`: text log/report copies use LF and omit trailing
  whitespace; their original run hashes are recorded. Raw runs are unchanged.
  All JSON receipts, package manifests and images are byte-identical copies.

Final event native measurements: review switching p95 13.25 ms, heartbeat p95
40.55 ms / maximum 424.07 ms, cancellation 19.69 ms (Propose) and 63.89 ms
(remeasurement). The heartbeat includes test-driven report creation/snapshots.
These establish the requested p95/short-cancellation checks, not a worst-case or
long-duration responsiveness guarantee.

Only named standalone/test targets were built. Synthetic fixtures and injected
outputs only; no physical device/listening, plugin installation, original-media
change, installed-app replacement, new worktree or publication. The prior
keyboard focus/Undo, sparse selection, transition bands and boundary batches passed.

Actual High Sierra/new-method cross-platform execution, broader detector controls,
independent holdouts and coordinated listening remain pending. Known music
partial/uncertain coverage and cinema/gate measured failures were not retuned away.
Existing video qualification and the accepted occasional picture skips are unchanged;
audio/drift/lifecycle blockers remain blockers.

To repeat only the finite frozen-candidate programme into a **new** receipt directory:

```sh
python3 Tests/Qualification/M20/run_ui.py \
  build/m20-complete-event-ui-2026-11-20/candidate-01 NEW_RECEIPTS \
  --case event-ui --case ui --case timeline \
  --case event-stereo --case event-surround --case event-fallback
python3 Tests/Qualification/M20Events/check_ui.py \
  build/m20-complete-event-ui-2026-11-20/candidate-01 NEW_RECEIPTS \
  build/m20-complete-event-ui-2026-11-20/authority-gate-01 NEW_ACCEPTANCE.json
```

GUI and the executable-owned verified-source cache must be accessible. Neither
command opens a physical device or reruns historical endurance/video/holdout suites.

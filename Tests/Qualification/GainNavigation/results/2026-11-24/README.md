# Gain entry and boundary navigation — focused qualification

The same candidate-01 passed 2,475 assertions on each actual host: Sequoia 15.7.9 arm64 and High Sierra 10.13.6 x86_64. The accompanying acceptance.json records exact implementation/package identities and per-group receipts. User acceptance is pending. Tests use owned fixtures and injected output; no physical device was opened. The user's running unsaved Tutti session was untouched.

Local preserved candidate:
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/timeline-input-boundary-fixes-2026-11-24/candidate-01`

Run with a fresh destination using the frozen candidate's runner (Python 3.8+):

```sh
python3 CANDIDATE/source/Tests/Qualification/M20/run_ui.py CANDIDATE FRESH_RESULTS --case gain-input-ui --case gain-entry --case graph-boundaries --case current-sections-ui
```

The runner verifies the frozen package's 473 artifact hashes before/after and requires matching build identity, successful exit and every assertion PASS. Each result JSON and log is retained with the host summary. Owned native windows need a logged-in desktop; the new click test deliberately positions its fixture in front so macOS hit testing dispatches the event. Do not substitute an open user session.

The new native group has 221 assertions. Existing focused gain-entry/manual checks contribute 1,079, graph/section/batch regressions 943, and preceding Current/section-copy native checks 232. No unrelated full historical suites were repeated. Documents 91–99 retain unaffected DSP/planner/model, source-integrity, transport and export qualification.

Development failure receipts are retained under development: the initial Mixed-field assertion was corrected to compare the exact prior display; two obscured-window native mouse injection attempts were corrected by positioning the owned window. The final development and frozen runs passed without relaxing the required behavior.

Broad video, full-film and three-hour endurance remain deferred, with no 90-minute substitute. The inherited completion pause is separate. Existing limited dropped-frame acceptance does not waive audio interruptions, accumulating drift or sustained freezes. No new audible, independent-programme correction or physical-latency claim follows from these UI checks. Main/tags/manuals/prior builds remain unchanged; no publication is authorized by this receipt.

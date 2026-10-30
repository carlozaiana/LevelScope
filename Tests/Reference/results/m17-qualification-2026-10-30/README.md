# M17 textual qualification evidence

30 October 2026. **Candidate-04: Sequoia tested; High Sierra, video and three-hour checks pending.** Stop for review before checkpoint publication/main integration. Main remains M13.

- [Final phase status](final-status.json): exact source/build/app identity, 12 selected commands, native UI/reuse and resources.
- [Candidate-04 immutable identity](candidate-04/identity.json): source digest `853a5711c162fbbbfa9daff2dbc1acee3fee18ee94445e2ac0ae3407efdeb904`; freeze-time pending wording intentionally preserved.
- [Fresh candidate-04 commands/results](sequoia-candidate-04/summary.json): all required outcomes in 153.69 seconds; the M16-reader exit 1 is the checked negative test.
- [Native candidate-04 UI observations](ui-candidate-04/observations.json) and [independent saved-session audit](ui-candidate-04/session-audit.json).
- [Universal/Intel deployment inspection](cross-build-candidate-04.json): not High Sierra runtime evidence.
- [Reuse audit](reuse-audit.json), [preservation](preservation.json), [preflight](preflight.json) and [direct GitHub refs](remote-preservation.txt).
- [Handover](../../../../Docs/ProjectIndex/62_M17QualificationAndHandover.md) and [pre-implementation impact plan](../../../../Docs/ProjectIndex/61_M17ImpactQualificationPlan.md).

`evidence-index.json` contains SHA-256 and byte size for 128 copied textual receipts. Raw audio, sessions with caches, app/test binaries, downloaded media, credentials and Finder metadata are not included. Frozen packages and synthetic working media remain in the local build tree. No original cinema/DAW project listing was used or included.

## Historical evidence remains historical

Candidate-01 stopped on an obsolete marker-test migration wording assertion. Candidate-02 fixes that test, reuses eight byte-identical passing commands and runs four remaining commands; see its equivalence record and `complete_candidate02.py`. Both original summaries remain unchanged.

Candidate-03 corrects the warning's text encoding and passes all twelve commands. Native review then finds an obsolete draft notice after identical confirmation. Candidate-04 refreshes confirmed/identically loaded definitions in the panel and adds a regression. All twelve commands are fresh on candidate-04; its affected native paths are also fresh.

Candidate-01 native cache/invalid-requirements recovery and cancellation, and candidate-03 native manual/proposal/artifact/report workflows are reused only for source-equivalent paths. Their executed source/build/app identities and defects remain in the original observations. Candidate-03 reports were independently hash-checked and visually reviewed in native Safari. This does not claim those UI actions ran on candidate-04.

Development failures include the initial sandbox restriction on owned temporary source copies, a compile-time test typo, and the base-versus-combined evaluation digest regression; later corrected runs are separate. No failed historical record was relabelled as a pass.

Actual High Sierra/Golem, cross-Mac exchange, video load/decode/playback and three-hour workflows remain deferred to the combined post-M18 candidate, no later than M19. No physical device was opened. No publication or integration occurred. `collect_evidence.py` copies only selected text, with equality checks before preserving existing receipts.

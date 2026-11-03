# M18 main proposal workflow and saved applied explanation

Development work in the managed UI/Timeline worktree. This slice is not the
combined M18 qualification or approval to integrate main.

Development app:
`/Users/fuzz/.codex/worktrees/a411/LevelScope/build/m18-dev-01/LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app`

## Implementation

- Main Propose, Review, Apply, Discard, limits and saved-explanation controls.
  The existing global-first numerical proposal algorithm and pass budget remain.
- Outcome-first progress, completion and rejection reasons. Duration diagnostics
  name the percentage, allowance and longest selected sections. Unsupported LRA,
  M/S or absent IL/TP requirements retain the existing capability limits. Baseline,
  prediction and private measurement have separate labels.
- Automatic local preview in the existing gain lane, with separately disclosed
  global gain, exact combined inspection and original per-section reasons. Stale
  and rejected previews cannot Apply. Manual gain drafts remain separate.
- A versioned, independently checksummed historical attachment in schema 8.
  Both changing Apply paths capture before invalidation. Binary64 gains/scalars
  remain exact. Original annotation names/intervals, full requirements/preferences,
  baseline statuses/provenance and available private candidate results survive
  reopening. No-op/Discard preserve history, later manual edits mark its association,
  and replacing Source clears it. History never authorizes verification or export.
- An 8 MiB history budget inside the unchanged 128 MiB session budget. The 4 MiB
  core and 16 MiB report bounds remain. Invalid or missing history is recoverable
  independently. No rows are silently omitted.

## Manual interaction checklist

1. Open and analyse a stereo or 5.1 mix, confirm applicability and set Target values.
   Press **Propose**. Check the outcome and progress below the waveform, and use
   **Review proposal** to read the complete baseline, prediction/candidate and reasons.
2. Click a proposed section in the gain lane, or use its arrows to visit changed
   sections. Inspect local gain, proposed global gain and their combined value.
   Open Proposal limits or change a Target field to confirm that the old preview
   is labelled stale and Apply becomes unavailable.
3. Generate a fresh ready proposal and choose **Apply proposal**. Confirm that
   Current becomes unverified. Use **Re-measure Current State** separately.
   **Discard proposal** should leave the applied sound and measurements intact.
4. Save and reopen the `.levelscope` session. Open **Last applied explanation**.
   Change global gain or a requirement and inspect the historical label. Restoring
   the old gain must not restore verification without independent remeasurement.
5. With a changing correction required, try an explicit **18 LU LRA** maximum.
   The outcome should explain the unsupported correction and say no candidate was
   measured. A duration rejection should show its affected percentage and selected
   allowance, with every proposed section still accessible through Review.

## Validation selection and limits

The runner builds no targets. It executes only the named proposal workflow,
section smoke, proposal preferences, global proposal, short session, timeline and
requirement suites. Synthetic audio, private fixtures and injected output devices
are used. No physical device is opened and no reference video is loaded.

The preliminary logs are retained in `build/`. Pilot 01 had an incorrect test
expectation about the successful save acknowledgement for a coherent older
revision; the newer edit correctly remained dirty. Pilot 02's injected peak map
put its synthetic peaks in a channel removed by the test's stereo resize; the
fixture was corrected. Pilot 03 passed the new workflow checks. Final results
and binary identities are recorded separately and supersede these pilots.

Actual High Sierra, physical device/video, cross-Mac and full-duration/RSS testing
remain deferred to the combined qualification block after M18, no later than M19.
The Intel deployment target alone does not establish actual High Sierra support.

Candidate 01 caught a real dense-display regression: repeated construction of the
10,000-section identity used growing UTF-8 string concatenation. The geometry key
now uses a linear byte stream shared by the model. A separate integer-format
oracle verifies byte-for-byte identity compatibility; Apply still recomputes and
checks the complete binding. Candidate 02 reruns the affected focused selection.

Candidate 02 passed the new workflow and section smoke checks but stopped at a
legacy test expecting the old uppercase first line. The measured candidate and
its extent remained exact; the assertion now expects the approved outcome-first
heading. Candidate 03 passed proposal workflow, section smoke, preferences,
global proposals and short sessions. Its evidence runner then encountered the
session suite's integer check count (other suites return arrays). The runner now
supports both report shapes. These intermediate logs are retained.

## Final development results

**Candidate 04 is final.** All 3,737 checks passed against source/build identity
`875fe8b80d1813513fc09c6ee0a4e5bd8e246bbd56d83128ac97a2b7794dbd8a`. The app and all five test executables are universal arm64/x86_64 and
retain Intel minimum macOS 10.13. Source/test hashes, exact commands, binary hashes
and architecture checks are in [the final manifest](candidate-04/summary.json).

| Focused group | Checks |
| --- | ---: |
| proposal-workflow | 363 |
| section-smoke | 366 |
| preferences | 256 |
| global | 303 |
| session-short | 1,071 |
| timeline | 889 |
| requirements | 489 |

The 10,000-change automatic proposal display completed controls + lane paint at
14.042 ms p95. A complete 10,000-row historical
record with 128-character four-byte Unicode names occupied
6,807,676 bytes, within 8 MiB. Encoding, reading,
formatting and complete paginated-access verification together took
428.879 ms on the development Mac. These are
short structural fixtures, not the deferred full-duration/RSS qualification.

The native main window was visually inspected with transport stopped and device
closed. Component captures cover [ready proposal](candidate-04/ready-proposal.png),
[automatic gains](candidate-04/proposed-gain-lane.png),
[duration rejection](candidate-04/duration-rejection.png),
[unsupported LRA](candidate-04/unsupported-lra.png) and
[measured rejection](candidate-04/measured-candidate-rejection.png).
No plugin target, main integration or external publication was performed.

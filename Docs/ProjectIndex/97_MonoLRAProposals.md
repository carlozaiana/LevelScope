# Confirmed mono in LRA proposals

The standalone app can now propose, apply and independently verify an LRA correction on a single channel confirmed as Left, Right or Centre. Measurement, processing and PCM24 export use that one channel with its assigned role. Playback routing remains a separate monitor choice.

## Method coverage

| Proposal settings | Existing method identifier | Tested mono operation |
| --- | --- | --- |
| Standard proposals (default) with an LRA ceiling | `lra-selective-trim-v1` | Whole programme, local correction and global gain |
| Measured-passage refinement | `region-gain-search-v1` | Automatic, selected section, temporary time selection and replacement by a measured passing alternative |
| Prefer complete-event treatment | `region-gain-search-v2` | Automatic, selected section and temporary time selection |
| Competing passage boundaries | `region-gain-search-v3` | Automatic, selected section, temporary time selection and manual refinement of an applied cut |

The three proposal services now use the same confirmed-layout check as section gain processing. The region prediction collector accepts one input channel and builds its energy bins from that channel alone. The selectable methods retain their search order, gains, margins and budgets. A changing LRA proposal requires confirmed roles, selected IL/true-peak/LRA requirements, stable and available LRA, and available values for the other selected metrics. The numerical qualification covers the unavailable and unstable cases, bounded no-pass, unconfirmed mono and excluded roles/layouts. Stereo and conventional 5.1 continue to use the same channel processing.

The existing Propose workflow verifies its starting Current as needed, privately measures candidates, applies the selected absolute plan and independently measures the resulting Current. Alternative replacement uses the original baseline. Cancellation before Apply preserves that baseline. Cancellation of post-Apply measurement leaves the plan available for remeasurement. A manual edit expires incompatible alternatives and requires a fresh Current measurement. Saved explanations survive reopening. Finalized export has its own measurement and file hash.

## Qualification on the frozen universal application

| Group | Sequoia 15.7.9 arm64 | High Sierra 10.13.6 x86_64 |
| --- | ---: | ---: |
| Mono service, independent sample oracle and verified export | 134 | 134 |
| Metric availability, stability and layout exclusions | 77 | 77 |
| Native mono workflow, all four methods, recovery and persistence | 153 | 153 |
| Established stereo control | 60 | 60 |
| Established conventional 5.1 control | 57 | 57 |
| **Total passed** | **481** | **481** |

On the 65-second, 48 kHz single-channel fixture, Source measured 5.910634 LU LRA. The job requested −23 LUFS ±0.2 LU, a −1 dBTP ceiling and at most 4.9 LU LRA. In method order, the independent Current results were **4.290034, 4.400000, 4.290034 and 4.290034 LU**. The measured-passage value is for the explicitly chosen replacement. Each method passed independent verification of the finalized one-channel WAV. The native checks exercised Centre, Left, Right and Centre respectively, with a single channel and the corresponding WAV speaker mask. All three region methods also respected both one-second extension limits for a selected section and a temporary selection.

Twenty exact comparisons passed. For each mono method the two machines produced identical settings, ordered roles, processing, input identity, complete Source and Current numerical payloads (including windows and peak positions) and finalized WAV bytes. Both multichannel controls matched the earlier qualified processing and payloads on their respective hosts. The preserved per-host qualification reports and exact fingerprints are in [the combined results](../../Tests/Qualification/MonoLRA/results/combined-2026-11-24/comparison.json).

A development run stopped after 17 checks because its fixture requested a measurement of an uninitialized Current after confirming a role. The corrected fixture initializes Current through the model before measuring it. The [first attempt](../../Tests/Qualification/MonoLRA/results/combined-2026-11-24/development-fixture-first-attempt.json) is retained. All five frozen groups passed on their first run on each host.

## Build and session identity

- Development branch `codex/mono-lra-2026-11-24` starts at published main `bb862cdc207f715853702601aa51defa5daa684a`.
- Tested implementation `81ea642b5b2870a934f3a5913c6736b0f6d41983`, Git Source tree `0cf444a4ffecf3712fe65955aa4eb84276bed958`.
- Frozen package manifest SHA-256 `d30dde12fd1b74befa017ad26c6cdaf2739fa3c03b0b337de9d01e57985d9087`. Each host checked all 471 artifacts before and after execution. [Candidate manifest](../../Tests/Qualification/MonoLRA/results/combined-2026-11-24/candidate-01-identity.json)
- Measurement build `17ab6500434d86940ac75133bfa614574a2464369bb842a0b2f99713c574ef9e`. The source-derived build identity advances. Existing measurement method, settings digests and session formats apply. Reopening on this build was tested with document schema 10 for the first three selections and schema 11 for competing boundaries. Earlier builds retain their original proposal eligibility and stored measurements remain subject to their build and identity checks.
- Sequoia package and output root `/Users/fuzz/Developer/LevelScope/build/mono-lra-2026-11-24`. Golem output root `/Users/fuzz/LevelScope-Mono-LRA-20261124`. Complete sessions and WAVs stay alongside the host receipts. The editable listening session is a byte-identical copy at handover. [App, session and source paths](../../Tests/Qualification/MonoLRA/results/combined-2026-11-24/user-test-session.json)
- Local main and tracking main stayed at the published baseline. The audit verified all 29 existing tag objects/targets and 5,666 retained file hashes across 15 older packages. [Preservation receipt](../../Tests/Qualification/MonoLRA/results/combined-2026-11-24/preservation.json)

## User check

1. Launch `candidate-01/LevelScope Standalone.app`. Open the editable `user-test/Mono competing boundaries.levelscope`. The generated single-channel source is referenced in place. Current should show **Passed**, with LRA about **4.29 LU**.
2. On the editable copy, choose **Initialize Current From Source**, retain **Competing passage boundaries** and the **4.9 LU** ceiling, then click **Propose**. Allow independent Current measurement to finish. The starting Source is about **5.91 LU**.
3. For listening, import a mono programme with an audible loud/quiet passage, confirm its actual Left, Right or Centre role and wait for Source analysis. With stable LRA, choose the job ceiling and click **Propose**. Compare Source and Current across the proposed cuts. An unsuccessful search should keep the starting processing. For a passing case, verify the exported WAV has one channel.

The two-host evidence covers device-free numerical and native UI operation. At the time of the frozen qualification, listening on a real mono mix was pending. The existing broad video, full-film and three-hour endurance deferrals continue without a 90-minute substitute. The previously recorded completion pause remains open separately. Publication and main integration require separate authorization.

## Tutti mono follow-up and user feedback

On Sequoia, all four methods passed independent Current and finalized one-channel PCM24 export verification on the full 16:56.85 `Tutti_MonoMix.wav`. The comparison interpreted the unlabelled channel as Centre and used a 19.2 LU ceiling against a measured Source LRA of 20.2990 LU. Competing passage boundaries selected the third of three passing attempts and produced 18.5393 LU LRA. The complete results and listening cues are preserved locally in [Tutti results.md](</Users/fuzz/Developer/LevelScope/build/mono-lra-2026-11-24/Tutti-test/Tutti results.md>).

After testing the workflow, the user reported **“I tested and all is working as expected.”** This records functional acceptance. The user gave no separate ranking of the four methods or listening verdict for individual boundaries.

## Proposed development sequence

1. **Keep Propose actionable after Clear Current State.** The disabled state was reproduced on an isolated copy of the saved Tutti session. Setting the LRA ceiling to 18.0 kept Propose available. Clearing Current kept Original measured and disabled Propose. Initializing Current from Source enabled it again. For a valid, confirmed Source with Current absent, an explicit Propose click should initialize a clean Current and continue through the usual measurements. An existing Current, including manual edits, must remain the starting point. When another prerequisite needs user action, show the reason beside the button.
2. **Create sections from applied regions.** Let the operator add the starts and ends of applied local gain regions to the editorial section map. Use the actual applied cuts, including any manually moved ones. Retain existing section boundaries and names, make repeated use idempotent and support Undo as one action. Creating sections alone should preserve audio and its verified Current result. The user can then edit a section’s level or refine a cut and independently remeasure the changed processing.
3. **Use a measured near miss as a manual starting point.** If a bounded search has no pass, identify the closest fully measured attempt under the stated comparison criteria and display which targets it missed and by how much. A deliberate manual edit command could load that plan from its original baseline and allow section-level refinement. Preserve recovery, stale-result checks and independent Current and export verification. Forecasts alone cannot supply that evidence.
4. **Simplify the release UI after the workflow settles.** Give the main view the result, next action and any immediate blocker. Condense Review proposal around audible locations, gain changes, measured outcomes and boundary uncertainty. Keep long identifiers, frame indices and full-precision values available in expandable technical detail for development and support.

These are follow-on proposals recorded from the user’s feedback. The qualified mono package remains the reference for this handover.

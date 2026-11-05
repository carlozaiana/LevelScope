# M19 — implementation and qualification handover

**Publication update, 5 November 2026:** the user subsequently approved publication.
Checkpoint [`m19-lra-passages-2026-11-05`](https://github.com/carlozaiana/LevelScope/tree/m19-lra-passages-2026-11-05)
targets `a14c88e8a94412eb2e2536860aa7cbd79de0dd39` and is included in main, together
with the Intel optimisation. The [documentation-only publication receipt](../../Tests/Reference/results/m19-publication-2026-11-05/README.md)
records exact refs and preservation checks. The original qualification-time status
and candidate identities below are retained as historical evidence, not current
publication restrictions. No rebuild or product change occurred during publication.

## Original qualification-time handover

5 November 2026. **Candidate-02 meets the approved impact-selected M19 qualification
scope on actual Sequoia and High Sierra; the coordinated boundary listening and
user file-inspection review is complete. Ready for checkpoint review.** No M19 tag, push or main merge
has occurred. Branch: `codex/m19-lra-passages`, in the established logic checkout.
Main stays at qualified M18 `51e47f300fcc38a386a5a92a7277029175cfab0c`.

The preceding Intel optimisation was published separately as
[`intel-multicore-true-peak-2026-11-05`](https://github.com/carlozaiana/LevelScope/tree/intel-multicore-true-peak-2026-11-05),
checkpoint `c4eea6f1c4569ec4f31725694ddf34db8340f15e`, followed by publication
receipt `09c8b8f6345c7c889f21f5399a38284ce1d353b2`. Its inherited completion pause
remains a separately accepted limitation, not a claimed fix in M19.

## Candidate identity

Candidate-02 source commit: `cf275afb133118dd2a06e7bca59f82f98898b3cd`.

- Source manifest: `0cc78f242abc6b1b799c04c2154fa25fb58ca57be295af2e047866c7d0a431e7`.
- Measurement build: `993f26a9295201a87fedea54a60bc8ca57877e93bb8bebc11729ab5bf5d2af1a`.
- Universal app executable SHA-256: `4bab04604bce5e7137dfb103ef880dd127a6cb402844281473307775d1c42536`.
- Local app: `/Users/fuzz/Developer/LevelScope/build/m19-qualification-2026-11-05/candidate-02/LevelScope Standalone.app`.
- Golem app: `/Users/fuzz/LevelScope-M19-qualification-2026-11-05/candidate-02/LevelScope Standalone.app`.
- Apple clang 17.0.0 (`clang-1700.0.13.5`), JUCE 7.0.12 commit
  `4f43011b96eb0636104cb3e433894cda98243626`; arm64/x86_64 universal, Intel minimum 10.13.

Candidate-01 remains unchanged: implementation commit
`2da3ac4ccb03a75976403c03bbc3461d00639752`, source manifest
`4cfde8ef7436d066e3c7d0e070db9158778d7dceed20afdfefbd94b641523b2d`, build
`0d1334c201daf85587e8d00dc261ed43d375a52e9f1d0a4136e4ac599479e18e`, app hash
`d9429446212fc299ee0efeafe4992676d635b78acd2bc73e632000971d292901`.
Both packages are preserved beside each other on both Macs. All 377/380 recorded
artifacts respectively were verified. Later documentation/evidence commits do
not replace these frozen identities.

## Implemented behaviour

The approved deterministic contract is in [66](66_M19LraPassagesAndMeasuredCandidates.md).
`lra-passages-v1` groups the existing validated LRA-window population by actual
three-second source-time support. Quiet/high population rows identify counts,
padding, score, percentile context, exact frame bounds and included-population
share. They are passages to investigate, not additive LRA contributions or
editorial scene judgments. No job LRA ceiling means descriptive inspection only.

`lra-selective-trim-v1` uses measured Current, matching Source inputs, selected
IL/TP/LRA requirements and active versioned proposal preferences. All selected
rules passing returns a genuine no-op. Otherwise it tries the defined global-only
case first when LRA already meets its ceiling, or a deterministic attenuation
envelope over high-population regions. It uses at most three complete private
assessments and one verified copy, retains the better actually measured attempt
if refinement worsens the result, and explicitly assesses any final global gain.
No policy or profile bounds are relaxed. Policy restrictions and method/search
limits are explained separately; failure is not a claim that no solution exists.

Apply requires an exact complete measured all-rule pass and stopped transport.
It invalidates Current, which must be measured independently. Cancelled verification
retains applied gains unverified. Predictions/private measurements never verify
Current or exported audio. Export still renders PCM24 with one final dither stage
and independently measures the finalized reopened artifact.

Source/Current inspection does not change audition selection, seek, play, open a
device or change audio. Propose with LRA opens the passage/estimate window;
**Assess candidate** is a separate explicit action. Tools also exposes descriptive
passage inspection. Selecting a row highlights its source interval. Complete
measured reasons remain in Review and the existing historical explanation model.

Editorial markers are unchanged. Dense/redundant markers do not alter the proposed
audio. Applied LRA processing geometry is recorded separately from annotation
geometry. The existing manual draft editor still initializes a clearly labelled
zero-gain draft when its annotation geometry differs from applied processing;
only an explicit later Apply can replace audio. M19 does not add a general editor
for analysis regions or silently create editorial sections.

## Fresh, reused and historical evidence

The explicit selection in [67](67_M19ImpactQualificationPlan.md) ran on:

- Sequoia 15.7.9, build 24G830, Apple M1 Max, native arm64.
- Actual High Sierra 10.13.6, build 17G66, Intel Core i7-3770K, native x86_64.

Candidate-01 passed all eleven direct groups on each Mac: core populations/limits;
stereo; 5.1; 44.1/96/192 kHz short workflows alongside 48 kHz; cancellation/failure/
staleness; separate-file/composition equivalence; the ten-minute gating case;
the eight-minute cinema excerpt; offline measurement; section-proposal smoke;
and job requirements. No recursive historical runner was invoked.

These checks include independent Current verification, finalized PCM24 verification,
all-channel/LFE pre-monitor sample parity at seeks/transitions/EOF, monitor-mute
identity preservation, history/session recovery, positive final global gain,
exact additional/cumulative/duration limits, marker-density independence, and
source replacement during assessment. A deeper-attenuation case crosses a gate,
worsens measured LRA, retains the earlier better attempt and still rejects Apply
when the retained candidate fails.

Candidate-02 changes only `LraPassagesPanel.cpp/.h` in product code: a custom-drawn
legend rendered missing characters and is now a standard label. Source comparison
is recorded in `local-integrity-and-reuse.json`. All measurement, proposal,
processing, requirements, session, export, audition and timing code is identical.
Candidate-01 numerical receipts are reused with their original build/binary hashes.
Fresh candidate-02 native tests on **both Macs** cover actual Source completion,
one-pass Current initialization, requirement reuse, descriptive Source inspection,
Current estimates, highlighting without seeking, explicit assessment/cancel,
measured Apply, independent Current verification and historical explanation.
The corrected screenshots were inspected on both platforms.

The final coverage audit added three short 5.1 workflows at 44.1/96/192 kHz,
linked to candidate-02's frozen objects, on both Macs. Together with the original
stereo rate group and separate 48 kHz workflows, these cover all four supported
rates for both layouts. They repeat measured assessment, independent Current,
session/history, finalized export and injected all-channel sample checks. The
selected envelopes/global gains and reported before/after LRA values agree
exactly across hosts. Supplemental helpers have their own source/executable
hashes in `supplemental-provenance.json`; they do not replace frozen binaries.

Cross-platform results selected exactly identical envelopes and global gains at
all tested rates/layouts. The compared bounded-cinema I/LRA/TP/P10/P95 metrics and
quiet/high passage arrays were identical. The reviewer retains established outer
comparison tolerances, which never widen profile acceptance. Within each platform,
candidate versus independent Current and pre-monitor samples versus the independent
processing oracle are checked exactly.

Unchanged M18/Intel measurement arithmetic, pinned EBU/ITU, input formats, large-file
handling and timing/video evidence is reused with its original identity. No video
was loaded or decoded for M19. This is not fresh video or three-hour audio-runtime
qualification; no measurement kernel, buffering, clock or video code changed.

Development failures remain in the evidence: a sandbox temp-directory failure;
an overly narrow expected error code for source replacement; missing text-change
notification in the first UI driver; and an initial synthetic gating fixture that
crossed the gate one attempt earlier than intended. Test corrections are identified,
not relabelled historical passes. The original unreadable legend screenshot is
retained. A resource wrapper also hit a sandbox `time`/`sysctl` denial after its
probe passed; the direct probe was recorded separately without that wrapper.

## Bounded resources and cinema result

| Check | Sequoia | High Sierra |
| --- | ---: | ---: |
| Eleven automated groups, elapsed | 304.26 s | 575.08 s |
| Supplemental 5.1 44.1/96/192 kHz complete workflows | 119.27 s | 233.93 s |
| Three-hour-sized population: grouping + envelope selection | 0.625 s | 1.485 s |
| Resource probe resident memory | 27.66 MiB | 23.17 MiB |
| Native inspection heartbeat p95, 20 ms timer | 24.12 ms | 29.19 ms |
| Native inspection maximum gap | 49.83 ms | 138.02 ms |
| Native active assessment cancellation | 81.74 ms | 69.56 ms |
| Eight-minute cinema workflow, peak RSS | 56.19 MiB | 39.68 MiB |

The population probe contains 107,985 windows and 1,080 high regions; 540 fit the
default envelope-duration budget. It includes planning, not just grouping, and
meets the 2 s/5 s and 64 MiB budgets. It decodes no three-hour audio. The old
section-proposal regression peaks at 500.47/421.85 MiB in its own fixture checks;
that is distinct from new passage-analysis storage. Local development, build,
packages and selected media total approximately 4.08 GiB at final review,
including the supplemental rate fixtures and requested listening reference.

Original cinema SHA-256:
`1fb7f93a33d600bf0a62d63e642e5703f03aec8e6920f8a6e510cf65f61fdeaf`.
Excerpt: frames 77,760,000–100,800,000, 27:00–35:00 at 48 kHz, L/R/C/LFE/Ls/Rs,
PCM24 payload copied byte-for-byte. Excerpt SHA-256:
`cfd57bea9a108db850b4e371de0f299f3396fbc1730e58e357380042a3064dd5`.
The original was not modified; original and excerpts remain outside Git.

The excerpt measured −20.909512 LUFS, 23.435333 LU LRA and −1.925859 dBTP.
With the explicitly chosen **15.0 LU test ceiling**, the retained complete attempt
measured 22.385113 LU after two private passes. The selected region and transitions
covered 19.0625% of the excerpt; its 3 dB cap left uncovered estimated reduction,
and other regions exceeded the remaining duration budget. Apply remained unavailable;
Current stayed unchanged. Preferences were not tuned to force a pass. These are
excerpt results, not full-film results or a Nuendo comparison.

## Listening and readiness boundary

Eight-second baseline/candidate float captures surround one onset and return for
stereo and the bounded cinema attempt. Each has exact source-frame/plan/hash
provenance and no dither or hidden normalization. They are diagnostic captures,
not published exports; the cinema candidate remains a measured failure.
Optional cinema Left/Right-only files contain those two original channels, with
no downmix. They cannot qualify Centre/LFE/surround or full 5.1 transparency.

With explicit coordination, the user selected **Saffire Pro 40 on Sequoia,
48 kHz, outputs 1/2**, initially at −5 dB monitor volume. The four synthetic stereo
clips played twice. The user described the audible sine-wave transitions and
subtle click-like effects at transitions/endings as normal for that signal. This
is not recorded as inaudible or click-free. The four cinema Left/Right clips then
played at −5 dB; the user heard no clicks and requested more context at 0 dB.

The requested repeat played a 99.4-second unchanged Left/Right reference, then
the four eight-second clips at 0 dB, with explicit pauses. Every play reached its
exact EOF with zero recorded underruns; the device closed after each sequence.
These counters establish application transport results, not independently captured
analogue output. The reference is excerpt frames 8,592,000–13,363,200; its first
and last eight seconds exactly match the two unchanged boundary clips.

The user initially reported no clicks but ambiguous/missing-sounding clip halves.
Read-only sample inspection found complete captures, exact original PCM agreement,
extremely quiet first-half Left/Right audio in the onset pair, and a quiet second
half in the return pair. Stronger onset first-half content is in Centre, excluded
from the explicitly partial audition. This and unannounced quiet gaps made the
sequence difficult to follow; no truncation was found. The original observation
and analysis remain in `listening-cinema-content-inspection.json`.

After receiving both the stereo and complete six-channel files, the user reported:
**“I checked all 4 stereo and 5.1 files and they are flawless.”** This completes the
approved bounded clip review. The inspection method and hardware were not specified;
this is user-reported file inspection, not a new physical 5.1/7.1 routing certification
or a claim of transparency for every programme. No further playback was required.

`supplemental/listen` uses candidate-01's unchanged audition objects with the exact
`--allow-audio-saffire-minus5` opt-in. `supplemental/listen-reference-zero` uses
candidate-02 objects and the exact `--allow-audio-saffire-zero-reference` opt-in.
Both require explicit device/input/fresh-receipt arguments, stop between captures
and close the device. Their individual hashes and original receipts are preserved.
Any future physical listening still requires coordination; 100 ms transitions
alone do not establish perceptual transparency.

All required M19 checks are complete; `completion-review.json` records the final
readiness checks separately from the original numerical-only reviewer. Publication/tag/main merge
remain separately authorized actions. Full physical 7.1 listening, new LRA strategies,
quiet-specific boost, neighbouring-region coordination, dialogue analysis,
dynamics, 7.1.4 and VST3 parity remain outside M19. The accepted occasional picture
skip limitation and blockers for audio interruption, accumulating drift, sustained
freeze and identity/lifecycle failures remain in force.

## Entry points and reproduction

Core grouping/planning: `Source/Core/Proposal/LraPassages.*`; private full-programme
assessment/ranking: `Source/Standalone/StandaloneLraProposal.*`; dispatch and Apply:
`StandaloneSectionProposal.*` and `StandaloneSessionModel.h`; native inspection:
`LraPassagesPanel.*`; overlay/history: `TimelineGainLane.*`,
`AppliedProposalExplanation.cpp`. No measurement or processing kernel changed.

Configure the existing standalone build with
`LEVELSCOPE_BUILD_LRA_PROPOSAL_TESTS=ON`, `LEVELSCOPE_BUILD_MEASUREMENT_TESTS=ON`,
`LEVELSCOPE_BUILD_SECTION_PROPOSAL_TESTS=ON`, `LEVELSCOPE_BUILD_JOB_REQUIREMENTS_TESTS=ON`,
universal `arm64;x86_64`, deployment target 10.13 and the established JUCE 7.0.12.
Build only `LevelScopeStandalone`, `LevelScopeLraProposalTests`,
`LevelScopeMeasurementTests`, `LevelScopeSectionProposalTests`,
`LevelScopeJobRequirementsTests`. The frozen package contains the exact CMake cache,
flags, link commands and sources. No plugin/install target is needed.

`Tests/Qualification/M19/run.py` explicitly selects the eleven groups;
`Native.cpp` uses real product UI/components with controlled picker completions;
`review.py` checks both platforms and exact/toleranced comparisons. Missing groups
cannot mark numerical completion, and that reviewer deliberately leaves overall
`qualified=false` until package integrity, resources and coordinated listening
have been reviewed. Its original receipt is retained; the separate final
`completion-review.json` and `current-status.json` record completed qualification.
See the textual evidence in
`Tests/Reference/results/m19-qualification-2026-11-05` and the complete preserved
working receipts under `build/m19-qualification-2026-11-05`.

# M20 complete-event preference — Logic interface handover

19 November 2026. **Logic interface implemented and focused-tested on Sequoia.
New UI integration, actual High Sierra and combined musical qualification pending.**
The local checkpoint identity is recorded separately in the publication-free
[checkpoint receipt](../../Tests/Qualification/M20Events/results/logic-01/checkpoint.md).

The accepted contract is [80](80_M20CompleteEventFinalReview.md), implemented from
the combined base recorded in [81](81_M20EventImplementationAndImpact.md). Logic owns
Core/service/session work in `/Users/fuzz/Developer/LevelScope` on
`codex/m20-region-planning`. The existing UI working copy and all historical packages
remain separate and untouched. Main remains `22961f03de7919e1de631f8b13eb4b852108ac0b`.
No publication/main merge or new physical-device use is authorized by this checkpoint.

## Exact tested package

Implementation/evidence commit: `4cbcca7f6b63b792b9a8ccbe73dbd5702ee311c5` on
`codex/m20-region-planning`. The checkpoint receipt is a separate documentation-only
follow-up; no tested implementation changed after that commit.

`/Users/fuzz/Developer/LevelScope/build/m20-events-2026-11-19/logic-handover-03/LevelScope Standalone.app`

Source/test manifest `131ddde3204d54da6cfe63645004332a358ab9ce26e0561958ad2be76e6c8596`;
production build `12d44f32d59488dcc03c9d88edea96d4c03a98b439fd9cd37d8d15078276c6a7`;
app executable `b20070886b34336ed971db2f404cbe19d362c417aef63e808f4a29c651fd7fd3`.
The implementation commit records these exact files; later documentation-only receipt
commits do not replace the tested binaries. Package 01/02/03 app contents are identical;
test-only evolution and receipt reuse are listed in the [evidence](../../Tests/Qualification/M20Events/results/logic-01/README.md).

The UI owner should deliberately integrate the local Logic continuation into the
existing UI branch, preserving later accepted UI changes. No checkout has been moved
and no UI dispatch is performed by this handover. Begin with the schema-2 activation
control and explicit measured-alternative selection; do not enable fallback until
the integrated authority tests pass. Existing schema-1 sessions continue their old
behaviour.

## Integration contract

`RegionPlanningSettings::defaults()` remains schema 1 with its original serialized
method/digests, inward-only automatic boundaries and legacy search. Explicitly call
`RegionPlanningSettings::eventDefaults(mode, inputIdentity)` for schema 2. This selects
complete-event preference and comparative assessment; it does not enable outward
margins. In automatic mode, set `allowAutomaticMargins` explicitly before positive
left/right margins. Reading context never grants processing permission. The approved
1 s experiment and 60% policy are explicit test settings; the product's current 20%
policy and existing gains/reserve are unchanged.

Relevant entry points:

- `Source/Core/Proposal/EventFeatureDefinition.*`: full versioned heuristic definition.
- `Source/Core/Proposal/EventPlanning.h`, `EventFeatures.cpp`, `EventCoverage.cpp`:
  bounded streaming Source evidence, frozen need association, exact additional
  envelope coverage, deterministic complete-first selection and assessment scheduling.
- `Source/Standalone/StandaloneEventRegionProposal.cpp`: one auxiliary source read,
  bounded prediction/search and at most three exact private assessments.
- `Source/Standalone/StandaloneRegionProposal.h`, `RegionAssessmentReceipt.cpp`:
  immutable attempt review and privately constructed runtime measured authority.
- `StandaloneSessionModel::selectMeasuredRegionAlternative(resultId, attemptId)`:
  explicit selection only; existing Apply remains a separate stopped-only action.
- `SessionDocument`, `AppliedProposalExplanation`, `AppliedProposalSerialization`:
  versioned settings/history, with private runtime authority never restored.

Correction need comes from frozen Current loudness populations. Source-event evidence
comes from verified PCM after composition/assembly and before whole-programme gain.
The matching Current auxiliary envelope distinguishes processing-associated rises.
The 100 ms prediction sidecar, qualified meter, profile and processor remain separate.
No Source reconstruction by division, new DSP, hidden gain or new compliance rule.

An onset/decay bracket is conservative evidence, not a musical timing guarantee.
Complete means the chosen intended additional correction is reached by the earliest
plausible onset and maintained through the upper supported decay bracket. The actual
old/new processors determine coverage, including existing ramps and internal changes.
A late ramp cannot gain credit merely by avoiding a detected peak. LFE cannot create
loudness correction need; applied gain and all-channel peak constraints still include it.

## Required UI state distinction

| Service state | Display / allowed action |
| --- | --- |
| `noChange` | Every selected measured rule passes; no processing or remeasurement |
| `ready`, recommended complete attempt | Explain complete treatment evidence and actual measured pass; ordinary explicit Apply |
| `ready`, measured global-only | Separate label; no invented local-event coverage |
| `reviewRequired` | Passing alternatives exist but complete coverage is not established; inspect only until explicit selection |
| `ready`, explicitly selected alternative | Keep its supported-complete / known-partial / uncertain / mixed label; do not relabel as an automatic complete recommendation |
| `noSolution` / unavailable / cancelled / stale | No new Apply authority; show exact limits and observed results; no impossibility claim |

Keep `recommendedAttempt`, `selectedAttempt`, and the UI's inspected attempt separate.
Review-menu navigation never calls the selection service. Add **Use this measured
alternative** only for an owned actual passing attempt, displaying its exact gains,
coverage reasons and measured rules. Pass the immutable review's `resultIdentity`
and attempt's `attemptIdentity`, not a mutable index/metrics bundle. Selection changes
no samples, measurement or revision. The existing guard against applying while a
different diagnostic plan is displayed must survive integration.

The receipt has no public constructor or decoder. It retains exact scalar statuses,
per-channel peaks/locations, gates/percentiles and context provenance after the private
full window arrays are released. Selection and Apply revalidate the same requirements,
input/roles/range, processing and settings binding. A saved historical narrative is
not a receipt. Source/processing/requirements/settings/manual edits and New/Open stale
a proposal; monitor controls and video offset do not alter audio identity.

The immutable review's fresh UUID is its proposal-invocation generation. Keep the
existing async completion ticket/session-revision guard: the UI captures the request
before `job.begin`, so comparing the captured request ticket directly with the later
active job ID would incorrectly reject ordinary completions. New/Open changes the
runtime session identity. Inspecting a row must never mint a new selection token.

Apply constructs bounded history before mutation, changes processing once and marks
Current UNVERIFIED. Cancellation/failure of independent Current measurement retains
the applied gain and honest status. Finalized-export verification still measures the
actual file without reapplying gain. No candidate result becomes Current or export
verification.

## History and compatibility

New writer: document schema 10, programme representation 9. Explicitly accepted old
schema pairings remain unchanged; schema 9 cannot claim settings-v2. Older sessions
open without opting into event preference or rewriting the original. Saving requires
an event-preference reader; Save As retains an older-readable original.

Applied history attachment 3 records schema-2 feature/permission definitions, exact
settings, historical coverage/selection reasons and baseline/result provenance.
Legacy attachments 1/2 retain their definitions and digests. The accepted large-
narrative fix remains within the existing 8 MiB cap. Private selection tokens and
candidate window arrays are not serialized or recreated as current verification.
Core/session caps remain 4/128 MiB; compatible measurement caches remain replaceable.

## Focused evidence and next stage

The required test manifest is `Tests/Qualification/M20Events/check_logic.py`.
It requires every named fresh receipt, matching source/build identities and no
physical devices/video; absent or failed checks cannot produce Logic readiness.
No recursive historical runner is invoked. Build only explicitly named standalone
and test targets, never the default plugin/all target.

The new pure tests cover complete/late/early/internal/uncertain treatment, ranking,
permitted margins and fixed ownership. Streaming checks compare 5/10/20 ms cadence,
selected sample rates, chunking, antiphase and channel-role behaviour. Service tests
cover immutable receipt ownership, stale rejection, explicit selection, independent
Current peak/location equality, history/session round trips and actual PCM export.
Legacy region, LRA, manual-edit/timeline and native integration tests are selected by
changed callers. The existing bounded music, gate and cinema cases are development
material, not independent holdouts or mandatory target successes.

The frozen combined UI candidate must still complete new native UI/fallback checks,
actual High Sierra/cross-platform checks, and coordinated independent listening.
Intel 10.13 cross-build success is not High Sierra runtime evidence. Preserve the
accepted occasional picture-skip limitation and the existing audio-continuity,
authoritative-clock/drift and lifecycle blockers. No new video implementation or
physical test is part of this Logic-stage dispatch.

## Known effectiveness limit and required presentation

The supplied music/dialogue material produced passing measured alternatives but no
supported-complete recommendation in any of the four new-mode configurations.
The tentative early activity bracket at 10.05–10.16 s remains associated with later
activity; this is not detection of the user's approximate 12.283 s first musical hit.
The 20/40/60% comparison did not resolve that evidence limitation. Display the actual
**known partial + uncertain** result, with the omitted entrance/ending reasons, rather
than describing a passing LRA number as musical success. The cinema and gate examples
remain measured failures. Detailed numeric/forecast reports are preserved with the
test evidence; no new listening judgment has been obtained.

The detector constants remain an initial versioned hypothesis. Broader gradual,
quiet, noisy and overlapping-event controls and independent holdouts remain required
before a general musical-quality claim. Do not retune one example or change shipped
defaults during UI integration. A substantive method change should return to Strategy
with the concrete failure and proposed rule revision.

## Reproduce the selected build and tests

From the Logic checkout, with the existing JUCE 7.0.12 dependency:

```sh
cmake -S . -B build/m20-events-review \
  -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_OSX_ARCHITECTURES=arm64;x86_64' \
  -DCMAKE_OSX_DEPLOYMENT_TARGET=10.13 \
  -DFETCHCONTENT_SOURCE_DIR_JUCE=/Users/fuzz/Developer/LevelScope-dependencies/JUCE-7.0.12 \
  -DFETCHCONTENT_FULLY_DISCONNECTED=ON \
  -DLEVELSCOPE_BUILD_JOB_REQUIREMENTS_TESTS=ON \
  -DLEVELSCOPE_BUILD_LRA_PROPOSAL_TESTS=ON \
  -DLEVELSCOPE_BUILD_REGION_PLANNING_TESTS=ON \
  -DLEVELSCOPE_BUILD_REGION_UI_TESTS=ON \
  -DLEVELSCOPE_BUILD_TIMELINE_PROGRESS_TESTS=ON
cmake --build build/m20-events-review --target \
  LevelScopeStandalone LevelScopeEventPlanningTests LevelScopeRegionPlanningTests \
  LevelScopeLraProposalTests LevelScopeRegionUITests LevelScopeTimelineProgressTests -j 6
```

Use a fresh owned test root. Event test groups are `core`, `features`, `stereo`,
`surround`, `heterogeneous`, `fallback`, `population`, and `external` with an explicit
input/experiment. `run_selected.py` spells out the finite 26-group programme;
`check_logic.py` rejects missing/failed/wrong-build receipts and changed source files.
No command targets the plugin or installs an app. JUCE uses its executable-owned
cache for verified-source temporaries; a sandbox blocking that cache fails preparation
and requires authorized cache access, not a fabricated successful result.

The frozen configuration and link/flags records remain in package 03. Fresh future
builds are not substitutes for that tested package. Short test fixtures consumed
about 0.73 GiB; the compile tree about 0.55 GiB and each package about 0.12 GiB.
New feature reservation is bounded by 10 MiB within 64 MiB planning; actual measured
figures and process RSS are separated in 81. At most three complete assessments and
64 proxy plans remain the shared caps. Native p95 ≤100 ms and actual High Sierra
planning/cancellation targets must be measured after UI integration.

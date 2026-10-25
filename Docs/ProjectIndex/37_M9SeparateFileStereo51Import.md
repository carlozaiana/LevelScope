# M9 proposal — separate-file stereo/5.1 programme import

Status: **approved 25 October 2026; candidate-04 qualified, published and integrated into main.**
See [38](38_M9QualificationAndHandover.md) for the final scope and evidence.
Publication and integration require a separate readiness review.
Prepared 25 October 2026 from published M8 commit
`8e984889c5e73670f370c6c568c8ea1d7902a35c`, tag
`m8-reproducible-sessions-2026-10-25`. Read with [35](35_PostM7BroadcastWorkflowProposal.md),
[36](36_M8SessionsQualificationAndHandover.md), the qualified measurement contract
and `Docs/ProposalEngineSpec_v2.md`.

## Objective and completion boundary

Import two mono files as one stereo programme, or six mono files as one conventional
5.1 programme. The assembled programme must work through the existing Source →
measured Current State → global-trim proposal/application/verification → audition
and diagnostics workflow, and survive session saving, reopening and relinking.
No physical audio device is required for import, measurement, evaluation or processing.

This is channel assembly of one existing mix. Every input supplies one channel;
there is no summing, stem composition, dialogue analysis or processing expansion.
The EBU profile definition, digest, applicability and narrow verdict wording remain
unchanged. Assembly capability/identity changes must be versioned separately from
the numerical loudness method and delivery profile.

Completion requires the offline, automated routing/sample-capture, session and
silent UI gates on both actual Macs, plus the coordinated bounded muted video
regression block below. Full physical 5.1 listening is separate evidence and is not
a prerequisite. Physical full-7.1 listening remains pending independently.

## Recommended input and alignment contract

| Programme | Required role slots | Assembly order |
| --- | --- | --- |
| Stereo | Left, Right | L, R |
| Conventional 5.1 | Left, Right, Centre, LFE, Left surround, Right surround | L, R, C, LFE, Ls, Rs |

The table defines programme-channel order, not an inference about file-selection
order. Keep each filename visibly associated with its assigned role. Each role
appears exactly once and every selected file is used exactly once. Two/six channels
or matching filenames alone do not confirm layout. Mono decoder labels, often
Centre, do not establish that file's role in the assembled programme. Existing
qualified conventional-5.1 surround aliases retain their interpretation and 1.41
loudness weight; System I rear roles must not replace them.

Recommended defaults:

- Exactly two or six distinct mono files. Reject duplicate selections of the same
  underlying file, including path/symlink/hard-link aliases. Distinct files with
  identical content, including silent channels, are legitimate and remain distinct
  assignments; do not deduplicate them by hash.
- Use the already qualified mono encodings in WAV, AIFF and native FLAC at
  44.1/48/96/192 kHz. Matching sample rate and decoded frame count are mandatory.
  Containers and supported bit depths may differ; retain each file's provenance.
  No new codec, rate, floating-point precision or container-size capability is implied.
- File sample zero maps to programme sample zero for every channel. Require explicit
  confirmation that the files represent the same aligned mix/version. Equal rates
  and lengths establish technical compatibility, not semantic or acoustic alignment.
- No automatic BWF/timecode alignment, arbitrary offsets, resampling, padding,
  truncation, normalization, polarity changes or silence substitution. State that
  metadata offsets are not applied. Unsupported format/rate/count, a non-mono input
  or even a one-sample length mismatch fails with a file-specific explanation.
- New import, role changes and relinking are stopped-only. Stage validation and
  confirmation before replacing the live programme. Cancel/failure keeps the prior
  valid session stopped. Successful new-group import follows explicit source
  replacement: old audio results/proposals and old programme-bound sections are
  cleared, and Current State must be initialized for the new programme.

## Architecture and dependencies

M8's serialized core already separates programme, files and assignments, but its
runtime `SourceDocument`, `ServiceContext`, `VerifiedSource` and measurement/audition
readers still carry one `sourceFile`. M9 must remove that assumption at the source
boundary, rather than manufacture a temporary interleaved delivery file.

Recommended internal boundaries (names are provisional):

1. **Programme source descriptor:** immutable, versioned assembly mode, per-file
   IDs/content hashes/format/locators and ordered channel assignments. An interleaved
   adapter retains the existing file-channel order; a mono-files adapter maps each
   file's sole channel into its explicit role slot. Retain twelve-channel capacity
   and existing interleaved 7.1 behavior without implementing eight/twelve-file import.
2. **Verified programme source:** an all-or-nothing set of privately owned, verified
   file copies. Verify every original identity and decoded format; expose the group
   only after all members succeed. Reuse that immutable group across seeks and
   stopped Source/Current or gain changes. Preflight aggregate disk requirements,
   check every write/flush and remove only owned partials on failure/cancellation.
3. **Shared worker-only programme reader:** read a bounded block at one programme
   frame across all channels. Measurement and audition consume the same ordered
   samples and the existing full-precision static-gain processor. A short read or
   error in any channel fails the group; it never becomes a zero-filled channel.
   Decoder work stays off the real-time audio callback. No duration-sized PCM buffer
   or additional full-length interleaved copy is needed.
4. **Existing consumers:** adapt Source/Current measurement, M2 binding, M3 matching
   inputs, M4/M5 locations/sections and M6 audition to the descriptor. Keep numerical
   DSP and profile rules unchanged wherever possible. The internal render test
   adapter may consume this reader; it does not expose public export.

Implementation order after approval: descriptor/identity and interleaved parity;
verified group/shared reader and capture parity; import/session UI and migration;
Sequoia stabilization, candidate freeze and actual High Sierra qualification.
Keep the established logic working copy and preserve the separate video worktree.

## Identity, freshness and safe reopening

Introduce a versioned programme-content identity covering assembly mode, every
member's content/decoded format, ordered channel mapping/roles, exact alignment,
rate, frame count and programme range. Pathnames, display names, random document
IDs and monitor routing are not audio-content identity. An equivalent interleaved
file may yield identical samples while retaining different provenance and identity;
numerical parity must not silently authorize cache reuse across representations.

| Change | Required effect |
| --- | --- |
| Any member's bytes, assignment, roles or alignment change | Reject stale measurements/evaluations/proposals and prepared audition; no partial group accepted |
| Processing gain change | Source remains valid; Current measurement/evaluation becomes unavailable and audition UNVERIFIED |
| Same-content verified relink | Update locator/document revision; compatible measurements can be retained after the complete group is verified |
| Profile/applicability change | Reevaluate compatible measurements without changing programme audio |
| Monitor volume/routes or video attachment/offset | Preserve programme measurement identity; retain existing stopped-only routing and live monitor-volume rules |
| Future analysis-only reference change | Invalidate only dependent future analysis, not unchanged audible programme measurements; references remain unimplemented |

Alignment confirmation is persisted against the exact versioned assembly definition
and member content identities. Unchanged reopen and identical verified relinks retain
it; any different content, mapping or alignment requires renewed confirmation.
Names, matching rates and matching lengths never establish alignment.

Source and Current snapshots must describe the same complete ordered group, range
and roles; Current must additionally match the active processing revision. Bind all
jobs to the runtime opening/ticket. Late verification or measurement results cannot
publish into a newly opened/replaced session. Predictions remain non-authoritative.

Extend the session core and assembly schemas explicitly for multiple files; keep
the existing container/checksum boundary if it needs no change. Recommended policy:
read and migrate M8 schema-1 sessions, write a version-2 core/assembly, and clearly
identify the format upgrade. Preserve the original file until an explicit save;
offer Save As for users retaining an M8-readable copy. Older M8 builds must reject
the new core safely. Keep conservative exact build/platform cache compatibility;
M8 caches will normally require remeasurement in an M9 build.

Save one coherent revision of the whole group, settings, sections and video offset.
Reopening verifies each file and reconstructs evaluations from eligible measurements.
Missing/changed files leave the entire programme explicitly incomplete: no partial
audition, processing or inherited compliance pass. Settings and sections remain
recoverable. Relink by stable file slot with content verification; no silent filename
matching. Saving must reject destinations aliasing any member or the reference video.
Retain independently replaceable caches and M8's save/open cancellation guarantees.

## Minimal UI

Add an explicit **Import channel files…** action alongside interleaved import. Its
staging view selects Stereo or 5.1 and shows a row for each programme channel:
role, filename, format, rate, frame count/duration and validation/availability.
Keep role assignment explicit and provide one alignment confirmation. Confirm only
when the group is complete and valid. Report aggregate and per-file progress and
allow cancellation.

The main view identifies an assembled programme, with filenames beside roles in
mapping/routing and per-channel inspection. Interleaved inputs retain “Channel 1…”
labels. Session reopening lists missing/changed members and supports explicit
per-slot relinking. Unavailable groups and UNVERIFIED Current State remain obvious;
audition never implies measured compliance. No folder scanning or automatic relink
heuristics are required in M9.

## Acceptance tests

1. **Deterministic mapping and validation.** Exercise stereo permutations and all
   six-file selection permutations with distinct impulses/tones and explicit role
   assignment. Confirm exact output order independently of selection order. Reject
   duplicate/missing roles, reused file aliases, wrong file counts, non-mono inputs,
   unsupported encodings/rates and one-sample mismatches. Preserve valid identical-
   content files. Test a deliberately shifted equal-length fixture: the importer
   must neither correct it nor claim content-based alignment verification.
2. **Interleaved numerical oracle.** Generate matched mono groups and independently
   interleaved fixtures from identical PCM values. Compare every assembled sample,
   all full-resolution loudness/gate/LRA results and per-channel peak values/positions
   through Source and Current measurement. Preserve unavailable/silent/short-file
   statuses. Cover the existing rate matrix, supported containers/bit depths and
   mixed-container/bit-depth groups with known expected decoded samples.
3. **All-channel processing and monitoring.** Capture before monitor controls at
   zero gain, positive/negative non-round gains and existing gain limits. Require
   exact same-platform offline/audition parity on every channel, including LFE.
   Use unique highest peaks on a surround and on LFE, mute that channel only for
   monitoring, and verify unchanged measurements, compliance and M3 constraints.
   Retain overload whole-callback silence, explicit restart and monitor-ramp recovery.
4. **Identity and lifecycle.** Mutate any member between import, copy, validation,
   measurement and adoption; reject same-name changed content, stale tickets and
   individually valid but mismatched Source/Current groups. Test profile-only reuse,
   processing invalidation, routing/video invariance and failed/cancelled measured
   verification retaining the applied gain with an unverified/actual-failure status.
5. **Group I/O failures.** Inject space/create/read/write/flush failures at the first,
   middle and last file, as well as short reads, non-finite samples and cancellation.
   Preserve the previous live session and saved document; clean up only owned partials.
   Test missing members without substituting silence, including disappearance during
   reopening. No successful preparation may expose a mixed old/new file group.
6. **Session migration and exchange.** Open preserved M8 interleaved stereo/5.1/7.1
   sessions; migrate settings/sections/profile/rational alignment without rounding.
   Round-trip new groups, gain bits, role order, file identities, caches and dirty-state
   behavior. Transfer both directions between Macs, relink identical moved files,
   reject changed files, recover valid cores with unusable caches, and reject the
   new schema in the preserved M8 reader. All opens start stopped/device-closed/muted.
7. **Timing and long media.** Exercise sample-zero, arbitrary and near-EOF seeks,
   stop/resume/restart, cancellation and EOF across all files with a reused verified
   group. Use a real three-hour six-mono 48 kHz programme, complete Source/Current
   measurements, save/reopen and random/final-sample captures. Preserve exact duration
   and programme/audio-source coordinates, including diagnostic rankings and sections.
8. **Regressions and UI.** Run relevant M1–M8 suites on both platforms, including
   official pinned EBU/ITU suites because the shared source/measurement adapter is
   affected. Retain interleaved 5.1 interpretation, System I 7.1 and optional video
   diagnostic isolation. Silent UI covers staging, confirmation, cancel/failure,
   full measured workflow, missing-file relink, reopen and close during owned jobs.

Same-platform parity is exact where identical decoded float32 input and arithmetic
are used. Across Macs, retain zero-gain bit equality, at most one float32 ULP for
nonzero processed captures, 0.01 LU/dB metric differences and 0.1 LU for LRA, with
availability/status, role order, frame counts and peak locations exact. These test
budgets never widen profile acceptance allowances.

## Resource and video qualification budgets

Retain M8's three-hour/1,000-section document budgets: save and parse each ≤5 s,
additional persistence memory ≤256 MiB, UI p95 ≤100 ms, cancellation acknowledgement
≤1 s and document ≤128 MiB. Separately measure import/verification/preparation time
against total bytes; hashing entire media is not part of the five-second parse budget.

For the new assembly reader, require additional working memory ≤256 MiB above the
existing measurement/diagnostic result storage, bounded independently of duration.
Owned snapshot payload must be no greater than the sum of the selected file bytes;
free-space checks include explicit overhead/reserve and still check actual writes.
On the prepared three-hour group, seek/prefill must complete within two seconds and
cooperative preparation/open cancellation within one second on both Macs. Record
actual values and cleanup time rather than hide them in an overall pass.

Stabilize on Sequoia first, freeze one universal candidate, then run that exact
candidate on actual Intel High Sierra. Record source manifest, original build and
executable hashes and all reused evidence explicitly. No Intel-on-Sequoia substitute.

Run native/injected Source/Current H.264 CFR/VFR and ProRes Proxy/LT/422 checks using
the assembled audio. Verify exact rational offsets/seeks, authoritative audio-clock
and EOF behavior, and zero embedded-video-audio decoder outputs. Coordinate a
bounded physical block before opening devices: those ten short Source/Current
video cases and one 30-minute six-file Current/ProRes endurance run per Mac, all
routes muted. Do not rerun unaffected video format/color matrices merely to obtain
zero skips. Occasional skips retain the approved development acceptance and counts;
sustained freezes, accumulating drift, audio interruptions, identity and memory/
lifecycle failures remain blockers. A video-engine change requires revisiting the
affected regression scope, not silently reusing obsolete evidence.

Injected captures verify all channels even with no multichannel hardware. Partial
physical monitoring is explicit routing only, with remaining routes visibly muted;
there is no implicit downmix, channel dropping in processing, or reassignment.
Stereo listening or selected-channel playback cannot establish full 5.1/7.1 listening.

## Exclusions and following work

No 5.0/arbitrary-count/eight-file/7.1.4 importer, mixed mono/interleaved assembly,
offset editor, BWF/timecode auto-alignment, RF64/BW64 extension, media consolidation,
stem summing, analysis references, speech/dialogue analysis, section processing,
dynamics, new profile, automatic correction, public export/report, autosave, embedded
soundtrack or VST3 parity. Existing qualified interleaved capabilities remain intact.

Keep the model extensible for later audible components and separate analysis-only
references as described in [35](35_PostM7BroadcastWorkflowProposal.md), without
implementing them. After M9, propose verified PCM export and artifact-bound reporting;
resolve container-size needs for actual delivery formats before broad export claims.
Preserve tags, tested packages, installations, the video working copy, Pro40 and
system audio configuration. Build only named standalone/test targets after approval.

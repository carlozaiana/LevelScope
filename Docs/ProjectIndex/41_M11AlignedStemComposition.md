# M11 proposal — explicit aligned M&E and dialogue composition

Status: **approved for implementation and qualification, 26 October 2026**, with
the user's change-impact qualification revision and future-processing clarification
below. Checkpoint publication/main integration require a later readiness review.

Base: published M10 `fcd4c15fe78d6352f5e96faabf3eb52f942235e1`, checkpoint
`m10-verified-pcm-export-2026-10-26`. Remote main and `codex/logic-module` were
checked and both identify that commit. Continue in the existing logic working copy;
preserve all checkpoints, tested packages, installations and the video worktree.
Read with [35](35_PostM7BroadcastWorkflowProposal.md),
[40](40_M10QualificationAndHandover.md), [14](14_OfflineMeasurementContract.md)
and `Docs/ProposalEngineSpec_v2.md`.

## Recommendation and reason

Create one explicitly composed programme from an M&E component and one selected
dialogue component, then use the qualified measurement, global-trim, diagnostics,
audition, session and verified-export workflow on that programme.

M8–M10 now cover reproducible sessions, separate-file channel assembly and verified
PCM delivery/reporting. Stem composition is the next previously identified workflow
extension. It is distinct from M9: channel assembly places files into channels of
one mix; composition adds explicitly mapped component samples. It is also distinct
from an analysis-only dialogue reference, which must never be added to playback.
M11 does not infer whether a file actually contains M&E or dialogue, detect speech,
or decide a creative balance. Explicit confirmation covers complementary,
correctly balanced stems intended to be combined, and their exact alignment.

## Bounded scope and defaults proposed for approval

Exactly two audible components: one M&E and one dialogue version. No arbitrary
mixer, nested compositions, alternative-version library or per-component processing.

| Output programme | M&E component | Dialogue component | Explicit mapping |
| --- | --- | --- | --- |
| Stereo | Stereo | Stereo | Matching Left/Right roles are summed |
| Conventional 5.1 | Conventional 5.1 | Mono | Mono contributes to Centre only; other channels retain M&E |
| Conventional 5.1 | Conventional 5.1 | Conventional 5.1 | Matching roles are summed, including LFE |

Stereo/5.1 components can be existing qualified interleaved sources or complete M9
two/six-mono groups. Mono dialogue is one qualified mono file with explicit Centre
assignment. Retain existing format boundaries, including only the M10-qualified
stereo/5.1 PCM24 RF64 subset; this does not add mono or arbitrary RF64 import.
Require explicit ordered roles and a displayed mapping, irrespective of picker
order, filenames or channel count. Conventional surround aliases keep their
qualified interpretation and 1.41 loudness weight. Existing 7.1 single-programme
behaviour is preserved; 7.1 composition/export is outside M11.

Both components must have the same supported rate (44.1/48/96/192 kHz), exact frame
count and complete programme range. Sample zero is shared. The user confirms exact
alignment and intended composition against the versioned definition and every
member's verified identity. Names, matching duration and matching rates are not
proof. No resampling, offsetting, automatic alignment, padding, truncation, panning,
downmixing or inferred LFE feed. Mono-to-stereo dialogue is deferred rather than
inventing a centre-pan law. An unmapped contribution in the mono-to-Centre preset
is explicit; a missing required file is always an incomplete programme.

Both components have fixed unity gain. The existing full-precision, all-channel
global trim follows composition. Component gain, mute/solo, fades, dialogue ducking
and automated correction are deferred. Monitor routing and volume remain separate.
Language and version are user-entered provenance, not detected content attributes.

Unity component gains are only M11's implementation boundary. Preserve accessible,
separately identifiable original components and typed extension points for component
processing, composition and whole-programme processing. Source is always the
original confirmed composition; Current is a versioned processing plan, not a
permanent global-trim-only model. Future component DSP settings must bind the
appropriate Current audio identity, measurements, proposals and export approvals.
They must never rewrite Source or silently become audible through an analysis-only
reference. Unsupported future plans fail explicitly in M11.

The visible states are **Source — composed programme** and **Current State —
composed programme plus global trim**. Component measurements or an earlier final
mix cannot stand in for either. Source need not pass the selected rules; the M10
Current-and-artifact publication gates remain intact. No automatic normalization of
either component or of the sum is allowed.

## Shared sample path and numerical contract

Generalize the existing programme-reader boundary; do not add separate summing
implementations in measurement, audition and export. Use bounded worker buffers,
per-file verified private copies and one immutable composition request. The realtime
callback continues to consume prepared audio and perform monitor-only operations.
Do not render an extra full-length intermediate merely to support composition.

Freeze a versioned `aligned-two-component-sum-v1` operation. Decode each component
through the existing float32 reader. At mapped channels compute
`float(double(meSample) + double(dialogueSample))` in fixed M&E/dialogue order.
For mono-to-Centre, copy non-Centre M&E samples unchanged. Then apply the existing
`static-gain-v1` operation to every output channel, including LFE. Preserve this
rounding boundary; do not fuse the sum with trim or use block-dependent reductions.

Same-build composed samples must match an independently generated oracle exactly,
irrespective of block size, seek history or offline/audition/export consumer.
Use binary-exact fixtures for cross-platform exact summation tests; retain the
existing at-most-one-float32-ULP allowance for nonzero global trim. This is a test
comparison allowance, never an extra compliance tolerance.

Finite sums above full scale remain unchanged in the floating processing stream.
Reject non-finite inputs, results and overflow. Never clamp, limit, attenuate or
normalize secretly. M6's monitor-overload stop/restart policy applies to audition;
M10's full-scale, Current pass, encoded-artifact measurement and publication gates
apply to export. Explicit user global trim can reduce level, followed by actual
remeasurement. Dither still occurs exactly once at final PCM24 quantization, never
in composition, Current measurement or audition. Different composition provenance
may produce a different export seed than an equivalent interleaved oracle: compare
pre-quantization samples and each artifact against its own recorded seed.

Programme duration, frame zero and coordinates do not change. Reference video
continues to follow presented audio, with existing rational offsets, seeks,
stop/resume/restart and EOF behaviour. Embedded video audio stays undecoded/unrouted.

## Architecture, identities and invalidation

Current implementation has useful boundaries but no stem mixer:

- `ProgrammeAssembly` and `ProgrammeReader` describe one interleaved file or M9
  mono-channel assembly; they deliberately do not sum channels.
- `VerifiedSource` verifies and owns source/member copies; reuse this lifecycle
  across seeks and stopped Source/Current or gain changes.
- `ServiceContext`, `SourceDocument`, session serialization, export input checks
  and report provenance currently branch between single-file and M9 assembly.
  They need a shared programme definition, not scattered two-stem special cases.
- `StaticGainProcessor`, the measurement kernel, M2 rules and M10 quantizer should
  retain their numerical semantics. Adapt input identity at their boundaries.

Introduce a bounded, versioned component/recipe representation: stable component
ID, purpose, language/version labels, interleaved/group member definition, ordered
channel assignments, exact zero alignment and explicit output-role mapping.
Keep legacy definitions readable without pretending they are compositions.
Interfaces remain driven by channel count/ordered roles and preserve twelve-channel
capacity; M11's qualified output layouts are still just stereo and conventional 5.1.

Distinguish file SHA-256s, the composition-definition digest and any encoded output
SHA-256. A composition digest is not a claim that a single source file exists.
Audible identity includes every participating member's verified content, format,
roles, range, mapping, alignment and summing version. Current identity adds exact
global trim. Profile/applicability bind evaluations separately. Runtime/job identity
remains separate from persistent document identity.

| Change | Required consequence |
| --- | --- |
| Audible member content, role mapping, alignment or recipe | Reconfirm composition; invalidate Source/Current measurements, evaluations, proposals, prepared audition and export approval |
| Global processing gain | Preserve Source; invalidate affected Current results/approval; retain the edited gain on cancellation or failed verification |
| Compatible profile/applicability | Reevaluate valid measured audio; no unnecessary audio remeasurement |
| Identical verified relink | Retain composition/alignment identity; apply existing cache compatibility rules |
| Display label, language/version text correction | Dirty document/report context; no audio change; refresh report/approval |
| Video attachment/offset or section annotation | Preserve audio identities; update dependent display/report context |
| Monitor volume/routing | Preserve all processing and measurement identities |

Import/edit a composition in a transactional draft. Failure or cancellation keeps
the previous valid session stopped. Explicit commit is required to replace Source;
never combine a new component with old valid results. New programmes use the normal
Current initialization at 0 dB. An edit of an existing composition retains its
explicit global gain but marks Current unverified. Programme-bound section maps
must not be silently rebound to different input identity; preserve the prior saved
session and make replacement of annotations clear before commit.

Structural edits and processing changes are stopped-only. Monitor volume remains
smoothly adjustable during playback. Full missing/read-failure handling applies
to both components and all grouped members; no partial programme or substitute
silence may reach measurement, playback or export. Reject late completions after
New/Open/replacement, including reopening the same saved document.

Future analysis references must have a distinct non-audible purpose and dependency
identity. M11 does not expose or implement them. Changing such a future reference
must invalidate dependent dialogue analysis only, not unchanged programme audio.
Unknown component purposes/recipe versions fail explicitly; they never become
audible by default.

## Sessions, reports and standards boundary

Recommend session schema 3 for compositions, with readers for M8 schema 1 and M9/M10
schema 2. Preserve old documents until explicit Save, explain that new-format saves
require M11, and offer Save As to retain the old readable document. Conservative
build/platform cache compatibility and independently recoverable cache records
remain unchanged. Reopen stopped, device-closed and routes-muted. Incomplete media
groups retain recoverable settings and references but cannot supply an inherited
pass, partial playback or export.

Save one coherent immutable revision, including recipe, per-member identity,
alignment confirmation, exact trim, annotations, profile/applicability and rational
video alignment. Concurrent edits remain unsaved; obsolete jobs cannot modify a
new runtime session. Fresh measurements remain necessary after incompatible builds.

Version the JSON/HTML report schema for component provenance and explicit composed
Source wording. Include every member hash, mapping, summing identity and user
language/version labels. Keep measurement-only reports available for failing or
unavailable states, with no invented artifact identity. Export remains a separate
Render and verify / Publish export transaction, with fresh approval for changed
report context. Historical published bundles remain unchanged.

The existing EBU profile definition/digest stays unchanged. Evaluate the completed
whole programme, not each component as if it were a delivery master. Applicability
still requires the user's confirmation. [EBU R128-2023](https://tech.ebu.ch/docs/r/r128.pdf)
specifies whole-programme measurement without privileging speech/music/effects;
composition does not introduce a dialogue-gated method. Retain the -23 LUFS target,
explicit +/-0.2 LU allowance and strict -1 dBTP evaluation ceiling, with informative
LRA. Keep [BS.1770-5](https://www.itu.int/rec/R-REC-BS.1770-5-202311-I/en) measurement
and existing peak coverage on every channel. This does not certify stem balance,
dialogue intelligibility, metadata or complete broadcaster delivery.

## Acceptance tests and resource budgets

1. **Definition and confirmation.** All three presets; shuffled file/channel order;
   conventional WAV/FLAC 5.1 aliases; interleaved and grouped components. Reject
   mismatched rate/frame count, unknown/duplicate/missing roles, unsupported maps,
   duplicate physical files/aliases across required inputs, nonzero offsets,
   unsupported formats, missing/unreadable members and unconfirmed alignment.
   Matching filenames/lengths do not confirm anything. Identical verified relinks
   retain confirmation; changed member bytes or mapping require renewed consent.
2. **Independent sample oracle.** Per-channel identifiable impulses/tones, correlated
   addition, phase cancellation, silence, non-finite inputs, overflow and above-full-
   scale sums. Cover all four rates, all presets and zero/positive/negative exact
   global gains. Compare independent precomposed float32 files with shared-reader
   output, offline processing, pre-monitor captures and pre-quantization export.
   Test chunk boundaries, seeks and first/middle/final frames. LFE and muted
   monitoring channels must still affect samples, peaks and M3 gain constraints.
3. **Measured workflow.** Source and Current must describe the same exact composition;
   reject individually valid mismatched snapshots. Run M2 full-precision boundaries,
   M3 explicit Apply/remeasure, already-passes no-op and unsuccessful-verification
   cases. Actual failure stays visible; predictions never establish a pass.
   Compare M4/M5 diagnostics, peak locations and time coordinates to the oracle.
   Keep section statistics descriptive and separate from new compliance rules.
4. **Audition/lifecycle.** Injected all-channel captures, Source/Current selection,
   stopped-only edits, smooth monitor volume, overload silencing and explicit
   restart after volume reduction, queue generations and cancellation. No physical
   device is needed for composition, measurement or processing qualification.
5. **Sessions/reports/export.** Both-direction cross-Mac exchange, old-schema read,
   explicit save migration, exact gain/recipe/alignment round trips, member/group
   relink, incomplete groups and damaged-cache recovery. Test coherent save and
   late jobs. JSON/HTML agree on provenance and domains. Reopen and measure actual
   PCM24 output with zero extra gain; use M10's independent decoder/quantizer checks,
   exact hashes, tampering tests and strict publication gates. Test failing Current
   measurement-only reports separately. Preserve existing stereo/5.1 export and
   interleaved stereo/5.1/7.1 import/measurement/audition behaviour.
6. **Failure safety.** Inject disk-full/read/copy/write/flush/finalize/report/rename
   failures, source mutation, cancellation at each stage and app close. Clean only
   owned partials; preserve previous session and sources; never publish a partial
   or stale export. Reuse verified copies only when every participating identity
   matches. No devices open automatically on import, remeasurement or reopen.
7. **Actual long programme.** On each Mac run a three-hour 48 kHz 5.1+5.1 composition
   through Source/Current measurement, bounded injected seek/capture checks, RF64
   render/reopen/measurement and explicit publication, with 1,000 annotations.
   Target total automated workflow time <=180 minutes, excluding human review and
   setup of independent oracle media. Additional composition/export working RSS
   <=256 MiB above existing measurement-result storage; interaction p95 <=100 ms;
   cancellation acknowledgement <=1 second. Preserve session <=128 MiB,
   save/core-cache parse <=5 seconds and combined JSON/HTML <=16 MiB budgets.
   Time unavoidable member hashing/copying separately, including it in total
   workflow time. Do not extrapolate 48 kHz long-run performance to 192 kHz.

Temporary disk budget is the sum of required unique verified member copies, one
export artifact and bounded report/metadata overhead/reserve on each affected
filesystem. Reuse unchanged copies across consumers; no extra full-size composed
render. Independently test arithmetic up to three hours at 192 kHz and deliberate
space failures. Record actual peak temporary bytes and cleanup latency. Plan fixture
reuse and free space before qualification; do not duplicate large historical video
packages or remove earlier evidence to make a test fit.

## Qualification sequence and completion criteria

First implement and stabilize on Sequoia. Freeze a new universal candidate with
source manifest, binary hashes, versions and exact test media identities. Then run
that same candidate on actual Intel High Sierra; no Intel-on-Sequoia substitute.
Any subsequent implementation fix requires a new candidate and affected reruns on
both Macs. Build only explicitly named standalone/test targets; no VST3/install.

Use the pre-implementation [impact plan](42_M11ImpactQualificationPlan.md), recording
affected behaviour, shared callers/dependencies/build changes, selected fresh tests,
reused evidence and estimated time/disk needs. Use short targeted development tests.
On the frozen candidate run required M11 and affected regressions explicitly; do not
invoke recursive historical qualification runners. Required-test absence must fail
the aggregate gate. Because the shared programme input boundary changes, retain
pinned official EBU/ITU checks through the affected adapter, while preserving the
original measurement algorithm and profile rules.
Keep existing cross-platform measurement budgets (0.01 LU/dB, 0.1 LU LRA), exact
statuses/roles/ranges and identical-oracle peak locations. Tests do not widen
acceptance thresholds. Reuse unrelated codec/colour evidence only with documented
component equivalence and the original tested identities.

Use representative H.264/ProRes integration checks with composed audio and relevant
seek/offset/EOF/lifecycle cases. Full codec/colour matrices and physical 30-minute
endurance runs are conditional on affected timing/buffering/video behaviour or
observed problems, not automatic M11 gates. Coordinate any physical-device use
beforehand. Record skipped-picture counts/timing under the existing accepted
development policy; audio interruptions, accumulating drift, sustained freezes,
identity/lifecycle and memory failures remain blockers. No zero-skip retry campaign.
Run one consolidated three-hour composition workflow per Mac covering sustained
measurement, session/resources and verified RF64 export; do not repeat equivalent
earlier long tests. After fixes rerun affected checks and broaden only when failure
or uncertain dependencies justify it. Preserve original identities of reused evidence.

Complete silent interactive composition/mapping/confirmation, edit invalidation,
measure/trim/remeasure, audition readiness, save/reopen/relink and render/report/
publication checks on both actual systems. Record offline qualification, injected
all-channel playback/routing and any physical listening evidence separately. Stereo
device checks cannot establish physical 5.1/7.1 listening; missing multichannel
hardware does not block the approved offline/injected scope.

M11 is ready for checkpoint review only when all required automated, long-resource,
cross-platform and silent UI gates pass for the same candidate and any coordinated
physical block is complete. Publish neither a qualified claim nor main integration
while required checks remain outstanding. Prepare exact identities, limitations
and a handover, then stop for checkpoint publication approval.

## Exclusions and subsequent order

No speech detection, isolated-stem measurement product, speech-gated mix measurement,
dialogue/background diagnostics or intelligibility score. No creative balancing,
per-stem gain/solo/automation, dynamics, section processing, new delivery profiles,
arbitrary offsets, automatic alignment, mono-to-stereo panning, nested/multi-stem
mixers, alternative-language batch exports, extra container/metadata formats,
network publication, 7.1.4/Atmos or VST3 parity. Preserve twelve-channel architectural
capacity and the separate pending physical full-7.1 listening qualification.

After M11, propose analysis-only dialogue references and clearly identified
isolated-reference measurements. The longer-term goal includes proposals for
dialogue levelling/DSP relative to M&E after the broader whole-programme
proposal/apply/verify workflow is established. Future targets must name their
analysis method and acceptance criteria; explicit application must be followed by
verification of both the intended dialogue improvement and the complete programme.
Speech-gated mix measurements, dialogue/background diagnostics and intelligibility
remain distinct methods needing separate validation. No component DSP or dialogue
analysis is implemented by M11. Further processing and report/delivery additions
remain bounded proposals, driven by actual workflow needs.

Approved decisions: the three mappings above, fixed unity component gains,
sample-zero-only alignment, schema-3 migration and impact-selected qualification.
Mono dialogue contributes only to Centre in the explicit 5.1 preset.

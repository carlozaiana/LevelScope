# M14 proposal — complete mono processing and verified PCM export

Status, 29 October 2026: **Sequoia tested; High Sierra, video and three-hour checks pending.** Candidate-02 is accepted for development-checkpoint publication on `codex/logic-module`; main remains M13. See [52](52_M14QualificationAndHandover.md) for the publication record.

The user-approved revision below supersedes the original two-platform/three-hour
M14 proposal. Main stays at qualified M13. See [51](51_M14ImpactQualificationPlan.md)
for the pre-implementation selection and [52](52_M14QualificationAndHandover.md)
for fresh, reused and pending evidence.

Start from published M13 candidate-02, checkpoint
[`m13-bounded-section-proposals-2026-10-28`](https://github.com/carlozaiana/LevelScope/tree/m13-bounded-section-proposals-2026-10-28),
commit `a9be070ab9eb5cb4f591ff9d991d8999c3a016d0`. Remote `main` and
`codex/logic-module` both identified this commit during proposal review. Continue
in `/Users/fuzz/Developer/LevelScope` on the established logic branch. Preserve
the separate video working copy, every tag, frozen package and historical receipt.

Read [49 — M13 handover](49_M13QualificationAndHandover.md),
[44 — section processing](44_M12SectionTrimProcessing.md),
[47 — bounded proposals](47_M13BoundedSectionTrimProposals.md),
[39 — verified export](39_M10VerifiedPcmExportAndReporting.md), and the
[shared proposal specification](../ProposalEngineSpec_v2.md).

## Recommended outcome and reason for this boundary

Complete the existing **single-channel programme** workflow: import and confirm,
measure, annotate, manually trim or generate an existing bounded proposal,
explicitly Apply, independently remeasure, audition, save/reopen and publish a
verified mono PCM24 export with truthful reports.

This closes the explicit mono gap left by M12/M13 before adding another processing
method. It is a layout/capability extension of qualified behavior, not a new
automatic algorithm. Basic optional marker import remains the next separately
reviewed workflow improvement, followed by bounded EDL import and the separately
qualified DSP stages described in the roadmap. Marker import is not a prerequisite
for mono support or useful proposals.

## Repository findings and dependencies

| Existing implementation | M14 dependency / intended change |
| --- | --- |
| `Core/Measurement/OfflineMeasurement.cpp` already accepts one Left, Right or Centre channel, with weight 1.0 | Preserve its arithmetic and existing mono measurements; qualify the new downstream paths. |
| `Standalone/StandaloneProcessing.h` has a channel-count-driven sample evaluator, but section-plan validation allows only stereo/5.1 | Add the explicit supported mono-role capability without changing the M12 envelope, operation order or float32 rounding. |
| `Standalone/StandaloneSectionProposal.cpp` and the shared section method already use the context's channel count | Extend fixtures, guards and acceptance coverage to mono; retain M13 policy and private-candidate separation. |
| `Standalone/StandaloneExport.cpp` assumes two or six output roles | Add an explicit one-channel mapping that preserves the confirmed role. |
| `Standalone/Pcm24Wave.cpp` accepts stereo/5.1 masks and assumes an even payload | Qualify mono masks, odd-payload padding, size calculations, writer finalization and reopened inspection. |
| `Standalone/StandaloneMedia.cpp` restricts RF64 validation to the existing PCM24 stereo/5.1 subset | Extend to mono and add a separate bounded input parser for the cinema layout; keep finalized-export inspection strict. |
| `Standalone/SessionDocument.cpp` already stores ordered roles and section plans and rejects unsupported applied layouts | Preserve schema-4 representation if sufficient; qualify old-reader rejection and display the minimum-reader requirement for mono section plans. |
| Audition, diagnostics and reports consume shared programme/measurement identities | Exercise mono end to end and remove misleading layout-specific wording; do not introduce alternate processing. |

No plugin DSP, new dependency, new programme-assembly mode or measurement algorithm
is expected. Inspect shared callers and build dependencies before implementing;
record an explicit impact-selected test manifest rather than launching historical
milestone runners recursively.

## Exact supported scope

- One file containing **one programme channel**, using the existing qualified
  WAV/AIFF/FLAC input paths at 44.1, 48, 96 or 192 kHz, up to three hours. Extend
  the narrow PCM24 RF64 path to the same mono case. Existing input bit-depth rules
  remain in force; this does not add arbitrary RF64/BW64/BWF variants.
- Preserve the existing supported mono interpretations: explicitly confirmed
  **Left, Right or Centre**. These are distinct identities, not interchangeable
  labels. Centre is the recommended interpretation for an ordinary unlabelled
  mono programme, subject to user confirmation; do not overwrite supplied roles.
  An isolated LFE, surround, height or unassigned channel is not newly qualified
  as a mono programme. Preserve broader existing Source-audition behavior.
- Keep a mono programme distinct from an isolated dialogue stem or analysis
  reference. Selecting the programme QC profile still requires the user's
  applicability confirmation; one channel does not establish programme purpose.
- Enable M12 manual section gains, -12 to +12 dB, existing versioned 100 ms
  transitions and final global trim. Annotation edits remain separate from applied
  audio. Use the same absolute-frame processing evaluator for measurement,
  private candidates, audition and export.
- Extend the **existing M13 method** to mono without changing policy: already-pass
  is a genuine no-op; global-only correction first; at most 3 dB additional local
  attenuation, eight changed sections, 20% affected duration including transitions,
  and a 0.5 dB reserve inside the cap. One envelope and at most three private passes.
  Retain partial-reserve/uncovered-deficit explanations and honest policy failure.
- Explicit stopped-only Apply, Current marked UNVERIFIED, then independent Current
  remeasurement. Cancelled/failed verification retains the applied plan; completed
  failure shows actual measured values. No rollback, automatic correction or retry.
- Retain measured diagnostics, peak locations, annotations, exact gain/frame values,
  source/current identity, session saving/reopening and measurement-only reports.
  Show one channel consistently; unavailable values and LRA stability stay truthful.

## Measurement, evaluation and identity

Use the existing mono BS.1770 path with **one** weighted channel, no duplication or
dual-mono compensation. Left, Right and Centre each have coefficient 1.0 in
[ITU-R BS.1770-5, Table 3](https://www.itu.int/dms_pubrec/itu-r/rec/bs/R-REC-BS.1770-5-202311-I%21%21PDF-E.pdf).
Sample and true peak refer to that same single channel.

Keep the authoritative EBU profile definition and digest unchanged: intended
-23 LUFS target, explicit +/-0.2 LU QC allowance and strict -1 dBTP ceiling,
informative LRA, existing pre-produced whole-programme applicability and narrow
verdict wording. Do not add meter/test tolerances or a mono-specific target.
The rule provenance remains [EBU R128-2023](https://tech.ebu.ch/docs/r/r128.pdf).

Bind Source, Current and proposals to matching content, ordered roles, range and
processing revision. A role change invalidates dependent results even when its
numerical loudness weight is equal. Profile changes reuse compatible measurements;
monitor routing/volume and video offset preserve audio identity. Predictions and
private candidates never satisfy Current or artifact publication gates.

Keep `section-gain-v1`, the measurement method, M13 policy and dither arithmetic
unchanged when their behavior is unchanged. Identify expanded product capability
and exact build separately. Preserve legacy stereo/5.1 processing samples and
deterministic export seed/PCM behavior; do not gratuitously change their dither by
folding a display/capability version into the seed. Any unavoidable semantic change
must be versioned explicitly and added to the affected qualification set.

## Mono export, monitoring and session contracts

**Export:** one-channel PCM24 WAVE_FORMAT_EXTENSIBLE, unchanged sample rate/frame
count. The explicit output map contains the sole input channel. Preserve its
confirmed speaker role: Left mask `0x1`, Right `0x2`, Centre `0x4`; do not silently
relabel it or turn it into stereo. These mask definitions come from
[Microsoft's WAVEFORMATEXTENSIBLE documentation](https://learn.microsoft.com/en-us/windows/win32/api/mmreg/ns-mmreg-waveformatextensible).

Use ordinary RIFF when the complete container fits and the existing limited RF64
form otherwise. Mono PCM24 payload length is `3 * frames`: an odd frame count needs
one trailing pad byte, which is **not** another audio sample and is excluded from
the data-chunk length. File/RIFF sizes include it, including the RIFF/RF64 boundary
decision. This follows the [RIFF chunk-padding rule](https://learn.microsoft.com/en-us/windows/win32/xaudio2/resource-interchange-file-format--riff-).
Audit disk estimates, parser bounds, finalization, hashes and report byte counts
together; merely allowing `channels == 1` is insufficient.

Retain one final deterministic unshaped TPDF stage at PCM24 quantization. Zero-gain
dithered export is not a bit-identical copy. Render and verify remains separate
from Publish export: measured Current **and** the finalized, reopened artifact
at zero extra gain must pass selected rules. Source need not pass. Preserve
overflow/non-finite rejection and no hidden clipping/attenuation. Failing or
unavailable Current still permits an accurately labelled measurement-only report,
with no invented artifact identity or verdict. No broad delivery/BWF claim.

**Monitoring:** reuse the existing explicit route from the sole source channel to
one selected physical output, or Muted. No automatic two-speaker duplication,
panning law, downmix or loudness compensation. A stereo device can monitor the
chosen output, but the signal and export remain mono. All offline operations work
device-closed. Processing/selection/role/routing changes are stopped-only; monitor
volume stays smoothly adjustable live. Preserve overload callback silencing and
explicit safe restart without changing processing.

**Sessions:** prefer existing schema 4 because ordered roles and section plans are
already represented. Warn before saving a mono section plan that it requires an
M14-capable reader; offer Save As to preserve the original file. Prove the frozen
M13 reader rejects that unsupported plan safely rather than dropping its processing.
An unchanged mono static-gain session should retain existing compatibility.
If inspection finds a structural extension necessary, version it and test migration
explicitly; never rely on an ignored optional field for safety. Preserve coherent
save revisions, dirty edits, late-job rejection, original provenance and recoverable
cache disposal. Reopen stopped, device-closed and routes-muted. Cross-Mac execution is deferred.

## Acceptance tests

| Required area | Passing evidence |
| --- | --- |
| Mono interpretation / levels | All four rates and confirmed L/R/C; ordinary mono WAV/AIFF/FLAC plus the bounded RF64 subset; invalid/missing roles and masks rejected. No extra channel or 3 dB adjustment. An independently calculated tone test and mono-versus-explicit-dual-mono fixture expose accidental duplication. |
| Exact processing samples | Independent section/global oracle across zero, positive, negative and full-precision gains, every transition edge, short spans, odd/even frame counts, varying chunks and random seeks. Offline, candidate and pre-monitor/pre-quantization captures agree exactly within one build. Cross-platform processed captures stay within the inherited one-float32-ULP bound. |
| Existing proposal policy | Already-pass preserves revision/results and performs no private work; global-only precedes local correction; success, no-solution, cap/count/duration boundaries, partial reserve and three-pass ceiling; transitions/EOF peak attribution, matching inputs, stale drafts and unsupported objectives. No new strategy or relaxed rule. |
| Truth and invalidation | Manual/proposed Apply makes Current unverified; independent measurement alone restores authority. Cancellation/failure retains the plan; measured failure remains visible. Role/content/plan edits invalidate; monitor/video/profile-only operations follow their existing rules. Muting the sole channel must not remove it from measurement or processing. |
| Containers / quantization | Complete independent decoding and PCM-code oracle for short mono exports; all rates/roles, one-frame and other odd/even lengths, exact padding and masks, malformed/truncated input, arithmetic just below/at/above the RIFF limit. One TPDF stage, deterministic seed, no extra frames, valid hashes and correct file lengths. |
| Publication / reports | Unverified/failing Current, failed finalized artifact, stale approval, missing media, disk-full, short write, flush/close failure, cancellation and destination collision cannot publish audio. A passing artifact requires explicit Publish. Source/Current/artifact remain distinct in JSON/HTML; measurement-only missing/unavailable values remain honest. |
| Session / UI | Old mono static-gain sessions, applied mono sections, exact doubles and frames, Save As warning, M13-reader rejection, changed/missing source, identical relink, damaged cache recovery, concurrent edits and late completions. Silent UI on Sequoia; actual cross-Mac exchange/remeasurement is deferred. No incidental device opening. |
| Audition / timing | Audio-only injected Source/Current captures, one-route/muted behavior, live monitor gain, overload/restart, cancellation/disk failure, device lifecycle, seek/resume/restart/EOF and continuous positions. Preserved audio-source clock, duration and video offset identities. |
| Existing workflows | Focused stereo/5.1 section/proposal/export parity, interleaved/grouped input and all three composition mappings; mono/5.0/7.1 static-gain regression, 5.1 surround weights and 7.1 muted peak constraints. No expansion of 5.0/7.1 section/export claims. |
| References and completeness | Pinned EBU/ITU suites through the current worker on Sequoia; required test list and receipts validated against frozen source/binary identities. Missing required checks cannot yield qualified status. |

Preserve exact verdict comparisons. Cross-platform checks retain M13's 0.01 LU/dB
I/M/S/TP/SP/gain and 0.1 LU LRA comparison bounds; they are not acceptance allowances
for programme rules. Explain any discrepancy and rerun affected checks after fixes.

## Bounded RF64 input compatibility

Container-field provenance is [EBU Tech 3306 (2009), section 3.4 and Annex A.2](https://tech.ebu.ch/docs/tech/tech3306v1_1.pdf):
ds64 supplies 64-bit outer/data sizes and sample count, with an optional size
table. M14's zero-count derivation, reserved-extension handling and stale outer
size recovery are explicitly bounded input compatibility policies, not claims
that those producer values satisfy the specification.

The original cinema file is `/Volumes/Tucan/GuerraDiCesare/LGdC_Masters/LGdC__51_Master_V2.wav`.
Read-only inspection confirms a 92-byte ds64 body (28 defined bytes plus 64 zero
bytes), zero sample count, table length zero, extensible six-channel PCM24, and a
complete 4,857,753,600-byte PCM payload at offset 16,384. The trailing `regn`,
`umid` and `DGDA` chunks account for the remaining physical bytes. Its declared
outer RIFF extent ends inside the first trailing metadata chunk, not at physical
EOF. This is an explicitly reported compatibility case, not a conformance claim.

Use separate input and exported-artifact parsing policies. Input support remains
PCM24 extensible mono (confirmed L/R/C), stereo or conventional 5.1, the existing
four rates and <=3 hours. Accept ds64 bodies 28..4096 bytes only with table length
zero and zero-filled reserved extension. Reject nonzero size tables, unknown
extension data, duplicate essential chunks, unsupported formats and bad masks.
For zero sample count, derive frames from validated data bytes / block alignment;
a nonzero count must agree exactly. Validate PCM subtype, bits, rate, alignment,
byte rate, divisibility, physical bounds, odd padding and all checked arithmetic.

Walk the complete file as bounded RIFF chunks, skipping non-audio payloads without
allocating them. Limit metadata to 128 MiB total and 4096 chunks. Accept an outer
size smaller than physical EOF only when it includes the complete padded audio
payload and all remaining bytes are structurally valid non-audio chunks. The
outer extent may end inside such metadata, as in the supplied file. Report this
size discrepancy and derived-count behavior explicitly. Reject outer sizes beyond
physical EOF or before the audio end, partial headers/padding, chunk overruns,
multiple data chunks, extra essential chunks and unresolved sentinel sizes. Metadata
never contributes samples. No filename, hash or producer-name exception is allowed.

Strict artifact inspection continues to require the canonical header, exact count,
physical/declared size agreement and finalized PCM structure. An imported file that
is compatible is not automatically a valid LevelScope export. Original content and
metadata are never rewritten; no metadata-preservation certification is implied.

Compare original duration, confirmed roles and bounded decoded samples with the
preserved M13 payload-identical copy, including first/last samples and offsets
beyond 4 GiB. Use the archived original/payload hashes as provenance and verify
original identity. Do not repeat the complete cinema proposal/export workflow.
Small synthetic fixtures cover the compatibility rules and malformed cases.

## Approved Sequoia qualification and resources

Before coding, record the impact table in document 51. Use short audio-only checks
for development, then freeze a new universal candidate. Run required M14 checks,
affected regressions, pinned references, independent PCM/report verification and
silent UI on actual arm64 Sequoia. Build the x86_64 slice at deployment target
10.13 using existing dependencies; this is not High Sierra runtime evidence.
No historical milestone runner is invoked recursively.

Use **one consolidated 90-minute mono/192 kHz workflow on Sequoia** with 1,000
annotations, proposal/Apply/independent verification, seek captures, diagnostics,
session/resource checks and verified ordinary PCM24 WAVE export. Payload size is
3,110,400,000 bytes, below the RIFF limit. Short cases cover other rates, L/R/C,
odd/even frames, 10,000 sections and strict failure handling. Separate constructed
or sparse files exercise RIFF/RF64 boundary arithmetic and reads beyond 4 GiB;
they are bounded container checks, not full-duration export qualification.

Retain intended three-hour product bounds but label **new three-hour mono/RF64
end-to-end qualification pending**. Defer all Golem/actual High Sierra execution,
cross-Mac exchange and all tests loading, decoding or playing video. Audio-only
clock/seek/EOF/continuity/routing tests remain required. Rational alignment can be
saved/compared without media decoding. Do not access Golem during this phase.

Plan one frozen combined qualification after approximately five milestones after
M13 (M14-M18), no later than six (M19). Each future milestone still requires its
own approval. Deferred manifest: actual High Sierra tests/UI/exchange, selected
video checks and new full-duration mono/RF64 workflow, selected by combined impact.
Earlier checks may be justified by changed macOS APIs/dependencies, measurement
arithmetic, buffering, clocks/latency or concrete failures; explain the specific
trigger before revising the deferral. Do not automatically restart historical suites.

Budgets: proposal/attribution/plan data <=16 MiB at 10,000 sections; additional
working RSS <=256 MiB above existing result storage; UI p95 <=100 ms; cancellation
<=1 s; save/parse <=5 s each; session <=128 MiB; JSON+HTML <=16 MiB. No full-duration
float cache. Proposed 90-minute workflow runtime limit <=120 minutes, with an
estimated 30-70 minutes execution excluding build/operator time. Short/reference
and UI work estimated 45-90 minutes. Preflight >=20 GiB free and calculated
source/copy/staging needs; sparse logical lengths are not physical disk budgets.
Record actual peak disk/RSS and clean only owned temporaries. Preserve packages,
symlinks and original media. No physical device opening is part of this phase.

Accepted occasional video skips remain a later qualification limitation, with
counts/timing preserved. Audio interruption, accumulating drift, sustained freeze,
identity/lifecycle/memory failures remain blockers. Full physical 7.1 listening
stays separate. Missing required Sequoia evidence blocks Sequoia-phase completion;
deferred tests cannot be labelled passed or silently omitted from the handover.

## Exclusions and review decisions

No marker/EDL import, automatic segmentation, new composition or separate-file
assembly preset, stereo duplication mode, pan-law setting, per-component DSP,
dialogue analysis, limiter/compressor, configurable proposal controls, LRA target,
multistrategy comparison, new export bit depth/dither mode, general RF64/BWF/BW64
support, additional video feature, immersive expansion or VST3 parity.

The user approved complete mono plus the bounded direct RF64 input compatibility
above, with the reduced Sequoia phase. Marker import remains a separate future
proposal. Keep main at M13 and use `codex/logic-module`. On Sequoia-phase completion,
report candidate identity, fresh/reused evidence, actual resource measurements and
pending High Sierra/video/three-hour coverage; stop for review before checkpoint
publication, main integration or the next milestone. No device opening, publication
or integration is authorized by implementation approval.

Subsequent review decision, 29 October: candidate-02 development-checkpoint
publication is authorized with the explicit pending-qualification status above.
This does not authorize main integration or the next milestone.

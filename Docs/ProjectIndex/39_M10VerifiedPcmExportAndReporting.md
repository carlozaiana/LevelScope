# M10 proposal — verified PCM export and reporting

Status: **approved 25 October 2026; candidate-03 qualified on both actual platforms 26 October 2026.**
This document retains the approved scope and acceptance plan; implemented results
and exact tested identities are in [40](40_M10QualificationAndHandover.md).
Checkpoint publication and main integration were authorized on 26 October 2026;
the published checkpoint and evidence mapping are recorded in [40](40_M10QualificationAndHandover.md).
Base: published M9 `75843dcd017ec7990be2200454a10fc38beb44c4`, checkpoint
`m9-separate-file-stereo-5.1-2026-10-25`. Read with [14](14_OfflineMeasurementContract.md),
[35](35_PostM7BroadcastWorkflowProposal.md), [38](38_M9QualificationAndHandover.md)
and `Docs/ProposalEngineSpec_v2.md`. Implementation stays in the
established logic working copy; preserve every checkpoint, package and video worktree.

## Objective and proposed completion boundary

Export a complete stereo or conventional 5.1 Current State, reopen and measure the
actual encoded PCM, and publish it with reports identifying exactly those bytes.
Source, measured Current State and measured exported artifact remain distinct.
Export never consumes monitor gain, physical routing, video audio or predictions.

Recommended scope includes limited large-file support: PCM24 RIFF/WAVE when it
fits, RF64 otherwise, plus validation/reimport of that same stereo/5.1 PCM24 subset.
This is an explicit extension beyond M9, requiring its own container tests. It is
needed for the existing three-hour programme workflow: at 48 kHz, PCM24 stereo is
3,110,400,000 payload bytes, whereas six channels are 9,331,200,000 bytes (8.69 GiB).
Headers add to these sizes. Silently truncating, splitting or rejecting all long
5.1 programmes would not meet this recommended scope.

Completion requires offline export, artifact verification, reporting, failure
handling, silent UI and relevant regression gates on actual Sequoia and High Sierra.
It does not require a physical six/eight-output interface or establish multichannel
listening qualification. Publication of M10 itself remains a later checkpoint action.

## Recommended export contract

| Property | M10 proposal |
| --- | --- |
| Audio state | Explicit Current State only; zero processing gain is allowed |
| Inputs | Qualified interleaved stereo/5.1 and M9 two/six-mono programmes |
| Duration/range | Complete programme, up to the existing three-hour supported envelope; identical frame count and sample zero |
| Rate | Preserve 44.1, 48, 96 or 192 kHz; no resampling or assumed 48 kHz conversion |
| Samples | Packed little-endian signed PCM24, 24 valid bits |
| Container | WAVEFORMATEXTENSIBLE in RIFF/WAVE or RF64; choice and estimated size shown before rendering |
| Output order | Stereo L/R; conventional 5.1 L/R/C/LFE/Ls/Rs |
| Speaker mask | Stereo 0x3; conventional 5.1 0x3f, with the existing back/surround alias interpretation |
| Quantization | Explicit round-to-nearest, ties-to-even, with deterministic unshaped TPDF dither; no noise shaping |
| Output set | One new folder containing programme.wav, report.json, report.html and a hash manifest |
| Publication policy | Current State and the written artifact must both pass the selected authoritative programme rules |

The export mapping is an explicit role-based permutation, displayed before work
starts and included in export identity. It does not change stored source order or
role identity. Conventional side/back surround aliases map to the established 5.1
surround pair and retain loudness weight 1.41; this must not introduce System I
rear weighting. Every channel, including muted channels and LFE, is exported.
No downmix, summing, channel omission or filename-based assignment is permitted.
Measurement-only EOF/filter drains and video offsets never add samples or padding
to the exported programme.

Microsoft specifies WAVEFORMATEXTENSIBLE's valid bits, mask and ascending speaker-bit
channel order in its [primary format documentation](https://learn.microsoft.com/en-us/windows-hardware/drivers/ddi/ksmedia/ns-ksmedia-waveformatextensible).
The RF64 compatibility target is the PCM/ds64 subset of
[EBU Tech 3306 v1.1](https://tech.ebu.ch/files/live/sites/tech/files/shared/tech/tech3306v1_1.pdf).
This is a deliberate interoperability choice, not a claim to implement the current
international exchange recommendation: [EBU Tech 3306 v2](https://tech.ebu.ch/docs/tech/tech3306.pdf)
redirects standardization to [ITU-R BS.2088](https://www.itu.int/rec/R-REC-BS.2088/en)
(currently BS.2088-2, November 2025). BW64/ADM remains outside M10.
Do not label plain PCM WAV/RF64 as complete BWF or broadcaster delivery compliance.
No copying or authoring of bext/iXML/ADM, timecode origins, loudness metadata or
section markers in the audio file is proposed. Relevant provenance goes in the report.

### Numerical contract

The pre-quantization stream must use the existing shared operation:
`float(double(sample) * pow(10.0, gainDb / 20.0))`, with exact stored gain and all
channels present. Compare this stream with offline measurement and pre-monitor
audition captures. Quantization is a separate versioned output operation, not a
change to Current State processing or its measurement identity.
Dither is applied exactly once, at the final PCM24 quantization boundary. Even at
zero processing gain a dithered export can change sample values; it is not a
bit-identical copy. Source/Current measurement and audition receive no export dither.

Use a documented TPDF generator with independent per-channel streams and one-LSB
peak amplitude (two-LSB peak-to-peak) at PCM24 precision. Derive and record a stable
seed from versioned programme/processing/format identity, excluding pathnames,
monitoring and profile choice. Freeze the generator, seed derivation, rounding
and known-answer vectors before qualification. Output must be independent of
worker block size, scheduling and repeated runs on the same qualified build/platform.
No hidden random seeding, silence gating or automatic dither changes. Dither can
make digital-zero input nonzero; report actual resulting metrics/statuses.

Reject non-finite samples, processed samples outside [-1, 1), and any rounded code
outside [-8388608, 8388607]. Never clamp, wrap, limit, lower gain or retry with an
undisclosed correction. For accepted samples, quantization plus dither error is
bounded by 1.5 PCM24 LSB. Silence and edge cases are exercised by the internal
render tests even when unavailable programme loudness prevents public export.

Same-platform pre-quantization parity is exact. Independently unpacked output PCM
must match the specified quantizer's expected integer codes exactly. Across Macs,
retain the qualified pre-quantization tolerance (at most one float32 ULP for
nonzero gain); allow at most one PCM24 code difference with the same dither stream.
Do not promise identical cross-platform file hashes where permitted floating-point
differences affect rounding. Each artifact has its own exact SHA-256 and measurement.
Equivalent grouped/interleaved inputs retain different provenance and may have
different derived dither seeds; compare their audio before quantization and each
encoded output against its own independent seeded oracle.

## Measured readiness and explicit publication

1. **Prerequisites:** complete verified programme, confirmed ordered roles and group
   alignment, matching valid Source/Current measurements, supported Current plan,
   authoritative profile and confirmed applicability. Current State must pass its
   selected rules. Source need not pass. No requirement to generate an M3 proposal
   if manually adjusted or already acceptable Current State passes.
2. **Render and verify:** user reviews format, channel map, dither, predicted disk
   use and destination, then starts one cancellable offline job. Playback stays
   stopped and the device closed. Hold a coherent immutable request; block session
   and processing edits until the job completes/cancels, while inspection remains
   available. Reuse the validated private programme copy where identity permits.
3. **Actual output:** finalize and close the encoded file, independently check
   format/mask/order/frame count and hash, reopen it, and measure decoded samples
   with the qualified kernel. Apply zero additional gain: the trim is already in
   the file. Keep original processing provenance separately from output decoding.
4. **Review:** display measured Source / Current State / Exported artifact and a
   separate artifact verdict. A passing Current result does not override a failing
   encoded result. Quantization-induced boundary failures stop here with actual
   measurements shown; retain the processing gain and do not automatically correct.
5. **Publish export:** a separate explicit action approves that exact artifact hash,
   request, profile digest/applicability and report revision. Recheck source-member
   content identities, staged output bytes and request freshness before publication.
   Publish the entire new folder only after all checks and report writes succeed.

The EBU profile definition/digest and full-precision limits remain unchanged:
-23 LUFS normalization target, explicit +/-0.2 LU QC allowance, strict -1 dBTP ceiling,
and no universal LRA limit. The standard's separate true-peak meter-tolerance
provision is documentary, not extra output headroom. See
[EBU R128-2023](https://tech.ebu.ch/docs/r/r128.pdf).
Missing/non-finite required metrics, unsupported applicability, placeholder profiles,
stale identities and predictions cannot authorize publication. Unavailable LRA
does not prevent a pass when there is no LRA rule.

Suggested successful wording: **"Exported artifact: passes selected programme
loudness and true-peak rules."** State that this describes these PCM bytes only;
it does not verify broadcaster delivery, metadata or a future lossy/video encode.
Failure/cancellation preserves the previous session and valid Source/Current
results: export itself has not changed their audio. No final audio set is published.
A failed-verification report may be saved explicitly, clearly marked **NOT EXPORTED**.
An I/O failure without complete output measurement must not acquire a measured verdict.

An independent **Save measurement report** action is available even when audio
export is unavailable, including a failing Current State before any render. It
captures only valid Source/Current measurements (or their explicit unavailable/
stale reasons), profile evaluation and session context. Its JSON and HTML are
labelled **MEASUREMENT-ONLY — NO EXPORTED AUDIO**. There is no artifact hash,
artifact measurement or artifact verdict unless an actual verified artifact is
included in a distinct artifact report. This report action cannot authorize or
bypass either audio-publication gate. Test no source, unsupported export layout,
missing/stale Current, failed profile, and report-write failure independently.

## Architecture and identity changes required

Inspection of M9 finds useful foundations, but not a ready public exporter:

- `ProgrammeReader`, `VerifiedSource` and `StaticGainProcessor` already provide the
  all-channel path. Retain them rather than introducing a parallel gain engine.
- `verifyTemporaryRender` is an M1 internal PCM24 harness. It uses a channel-count
  writer and is not an export format/metadata/transaction contract. Preserve its
  regression role while adding an explicitly specified production encoder boundary.
- `RenderService` and `ReportService` are unavailable interfaces. Add a versioned
  export specification, immutable render request, artifact-measurement receipt,
  report model and publication coordinator; workers receive no mutable UI references.
- The core evaluator has a `renderedOutput` enum but currently rejects it; the
  adapter also assumes the measured hash is the original programme hash, and the
  display labels non-Source results as Current State. Add explicit artifact-domain
  binding and wording, with an evaluation-version change. Preserve M2 numerical
  rules/profile digest and all existing Source/Current semantics. Do not fake a
  rendered artifact as Current State to bypass this boundary.
- `ExportApproval::matches` is only an identity scaffold. Extend/rework it to require
  full profile digest, applicability, export format/quantizer/channel-map identity,
  exact artifact hash, verified measurement identity, current runtime opening/job
  and coherent revision. Existing booleans alone cannot authorize publication.
- Pinned JUCE 7.0.12 has RF64 reading/writing, but capability is not qualification.
  Its writer selects RF64 from data size and finalizes headers through a void
  helper. Explicitly check the total RIFF size including headers, all 64-bit
  arithmetic, write/seek/flush/finalization errors and the near-4-GiB boundary.
  Prefer a bounded standalone encoder/container adapter over changing the shared
  JUCE dependency or plugin. Validate the output with an independent test parser.

Export identity covers programme/member hashes, assembly/alignment, input roles,
range, full-precision plan, output role permutation, format/rate/bit depth,
quantizer/dither version and seed, and encoder/build identity. Artifact measurement
adds output bytes/hash, decoded layout and measurement provenance. Evaluation
separately adds profile definition and applicability. Report identity adds its
captured session/document revision and report schema/version.

Audio or format edits invalidate pending export approval. A compatible profile-only
change reuses valid audio measurements and an unchanged encoded artifact, but
reevaluates and regenerates the report/approval. Video alignment and annotation
changes preserve audio identities; if included in the report, refresh that context
and approval without rerendering. Monitor gain/routes never affect exported samples,
measurements or reports. After publication a receipt describes historical bytes;
later edits cannot make it a pass for the new Current State.

Do not persist active jobs or approvals in `.levelscope`. Keep schema-2 session
compatibility unless a separately justified change becomes necessary; no export
history subsystem is needed. An unsaved session can export an explicitly identified
captured revision, and the report must say it was unsaved. Reopening remains
stopped/device-closed/routes-muted and cannot restart or authorize an old export.

## Reports and filesystem transaction

One immutable report model produces versioned UTF-8 JSON and self-contained static
HTML. Include session/programme revision, original file basenames and SHA-256s,
assembly/alignment and role mapping, exact gain/processing identity, format/container,
dither identity, artifact SHA-256, exact frame count/rate/duration, app/build/OS/CPU
and measurement provenance, profile edition/digest/applicability and per-rule reasons.

Compare measured Source, Current and artifact I, sample/true peak, LRA and its
stability, maximum M/S, plus per-channel peaks/positions. Distinguish absent values
with explicit statuses/reasons; serialize finite values with round-trip precision,
and gain additionally with its exact representation. Rounded HTML numbers never
determine verdicts. Include section labels and exact frame intervals as annotations,
and rational video offset as context only, not proof of synchronization. Do not
invent section/dialogue measurements or include proposal predictions as results.
Keep filenames/labels HTML-escaped and default to basenames, without external assets,
scripts or network requests. PDF generation and interactive report plots are deferred.

Use a uniquely owned hidden staging directory under the selected destination parent.
Preflight output bytes, any new verified-source copies, report overhead and reserve
per affected filesystem; still check every actual operation. Finalize/flush and
validate all audio/report files before an exclusive same-filesystem directory rename
to a new destination. Never overwrite existing export sets, sources, group members,
session documents or reference video, including symlink/hard-link aliases. The manifest
hashes audio and both reports; reports bind the audio hash without circular self-hashes.

Cancellation, ENOSPC, short writes, failed finalization/hash/measurement/report/rename
and stale completions remove only owned partials and publish no successful set.
Crash tests must distinguish an incomplete staging directory from a completed bundle;
recovery must not auto-publish or delete unowned files. No resumable export is proposed.
Scope atomic publication to a local filesystem; network/cloud destinations are deferred.

## Acceptance tests

1. **Readiness and verdict boundaries.** Exact and immediately adjacent I/TP limits;
   passing Current but failing encoded output; no target/placeholder/unsupported
   applicability; missing/non-finite metrics, optional LRA, stale snapshots and
   predictions. Verify domain-specific labels and no automatic correction/retry.
2. **Samples and mapping.** Stereo/six-channel impulses, polarity and unique tones,
   arbitrary confirmed input order, conventional 5.1 back/side metadata aliases,
   grouped/interleaved equivalents, all supported rates, zero/positive/negative
   non-round gains and plan limits. Compare pre-quantization measurement/audition
   captures exactly and independently unpack every output code. Mute unique highest
   surround/LFE peaks for monitoring and prove measurements, M3 constraints and
   exported samples remain unchanged. Reject unsupported export layouts explicitly.
3. **Quantizer and container.** Known-answer packed PCM24, tie rounding, dither
   distribution/bias/channel independence, repeated-run/chunk-size determinism,
   silence/sub-LSB/near-full-scale and non-finite/overflow cases. Independently parse
   format, mask, GUID, data length/padding, RF64 ds64 and frame counts. Exercise total
   RIFF size just below/at/above the limit, malformed/overflow/truncated RF64 and
   incorrect masks. Use real multi-gigabyte output as well as virtual-size tests.
4. **Actual-file verification.** Reopen decoded output, compare it to the independent
   PCM oracle, measure with the qualified kernel and check exact output hash binding.
   Test deliberate sample/header/tail tampering and swaps before publication; no
   double gain or reuse of a pre-encode pass. Reimport supported produced WAV/RF64
   and verify frame/rate/role/measurement/audition parity. Obtain an independent
   decoder interoperability check, recording its version; the production reader
   and writer must not be the sole oracle.
5. **Lifecycle and filesystems.** Inject create/space/read/write/seek/flush/finalize/
   report/rename failures and cancellation at each stage, including source-member
   changes and report-only failure after successful audio verification. Exercise
   destination races, aliases, app close and process termination. Verify coherent
   revisions, no obsolete completion after New/Open, preserved input/session bytes,
   no published partial set and bounded cleanup. Frozen artifact/report approval
   cannot survive an incompatible edit or reopen.
6. **Reports and reuse.** JSON/HTML agree with full-precision measured records,
   domain and profile identity; null/silent/short-LRA states remain distinct. Verify
   escaping, exact rational/frame locations, hashes, unsaved-state annotation and
   report failure statuses. Profile-only reevaluation avoids remeasurement; monitor
   and video changes cannot alter audio identity. Exchange exported bundles between
   both Macs and verify hashes/decode/measurements independently; do not trust saved
   report verdicts as session measurement authority.
7. **Long-media/resources.** Render, reopen, fully measure and publish a real three-hour
   six-channel 48 kHz PCM24 RF64 programme on each Mac; exercise an ordinary RIFF
   stereo case, first/random/final-sample captures and report generation with 1,000
   section annotations. Additional export working memory <=256 MiB above existing
   measurement-result storage, UI interaction p95 <=100 ms, cancellation acknowledgement
   <=1 s, and combined JSON/HTML <=16 MiB. Target full preparation/render/verification
   wall time <=programme duration for that long fixture on each Mac, excluding human
   review. Measure each stage, peak RSS, actual temporary bytes and cleanup latency.
   Disk use is one output payload plus required verified-source copies and bounded
   overhead; no additional full-size interleaved render copy. Test rate/layout size
   arithmetic through three hours at 192 kHz without claiming untested long-run speed.
8. **Regressions and silent UI.** Relevant M1–M9 suites on both actual platforms,
   pinned EBU/ITU vectors for the affected reader/artifact adapter, M2/M3 identity
   and boundary tests, M4/M5 diagnostics/sections, M6 parity/overload, M8/M9 reopening
   and grouped relinks, and existing interleaved stereo/5.1/System I 7.1 behavior.
   Silent UI covers readiness, format/map/dither/size preview, progress/cancel/close,
   actual failure, publish approval, existing-folder rejection and readable report.
   Verify export does not start playback, open devices or change transport position.

Metric cross-platform budgets remain 0.01 LU/dB and 0.1 LU LRA, with exact statuses,
role order, frame counts and locations for the same decoded oracle. Differences
between separately quantized files must be identified, not hidden as measurement
engine error. Profile comparisons always use each artifact's full-precision values;
test tolerances never widen its acceptance limits.

## Qualification sequence and video coverage

Approved sequence: specify/test the encoder and RF64 boundary first; add immutable
rendering and output-domain evaluation; then transaction/report/UI integration.
Stabilize on Sequoia, freeze one universal candidate and source/artifact manifest,
then qualify that exact candidate on actual Intel High Sierra. Keep test-only
independent parsers/decoders separate from product dependencies and inspect available
runtimes first. Build named standalone/test targets only; no VST3 or installation.

Run the ten native/injected Source/Current H.264 CFR/VFR and ProRes Proxy/LT/422
cases on both Macs, including reimported exported audio, seek/offset/restart/EOF,
authoritative audio presentation and zero embedded-audio outputs. Retain diagnostic
trace isolation. Reuse unaffected codec/color evidence with a written equivalence.

Proposed physical block, requiring coordination before device opening: four short
muted Source/Current cases (H.264 VFR and ProRes 422), plus one 30-minute Current/
ProRes run per Mac with reimported exported PCM. Use the previously preferred
built-in outputs only after confirming availability/permission for this milestone.
No audible tests are needed for export qualification. Record picture skips/timing
under the accepted development policy; do not chase zero-skip runs. Audio interruptions,
accumulating drift, sustained freezes, identity errors and lifecycle/memory failures
remain blockers. An unexpected shared playback change requires affected coverage
to be expanded deliberately, not wholesale reuse of old evidence.

## Approved defaults and exclusions

Approved defaults: **PCM24 + deterministic unshaped TPDF;
limited RF64 for long files; canonical stereo/5.1 output mapping; strict measured
Current and artifact pass gates; explicit final publication; new-folder JSON/HTML
report bundle.** RF64 is the main scope addition relative to M9. If it is deferred,
M10 must instead explicitly reject outputs exceeding RIFF capacity and cannot claim
complete three-hour 5.1 export. Required recipient-specific metadata or BW64 delivery
should be scoped separately rather than implied by the selected EBU level rules.

Outside M10: separate mono-file export, mono/5.0/7.1 export qualification, 16-bit or
float output, resampling, lossy formats, video mux/export, BWF metadata authoring,
BW64/ADM/Atmos, stem composition, dialogue analysis, 7.1.4, section processing,
dynamics/limiting, automatic correction, new delivery profiles, general report
design/PDF, noncompliant-audio export bypass, batch processing, network publication
and VST3 parity. Existing broader
import/measurement/audition capabilities and twelve-channel architecture remain intact.

After M10, decide further delivery formats/metadata and report improvements from
actual recipient requirements. The original approval covered M10 implementation
and qualification; the later publication authorization and saved checkpoint are
recorded in [40](40_M10QualificationAndHandover.md). A next milestone still requires
separate user approval.

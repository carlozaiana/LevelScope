# M7 proposal — channel-based PCM 7.1 support

Approved 24 October 2026. **Combined candidate-03 qualified on both actual Macs on 25 October
under [the revised development acceptance](34_M7CompletionAcceptance.md), published
as `m7-pcm-7.1-2026-10-25` and integrated into main.** See [32](32_M7QualificationAndHandover.md)
for completed evidence and [33](33_M7VideoPaintTimingHandover.md) for the historical
user-requested video handover. This document makes stage one of
[29](29_MultichannelPcmMilestone.md) concrete and carries forward
[M6's qualified contract](30_M6QualificationAndHandover.md).

## Baseline and recommendation

Fetched `origin/main` and the established `codex/logic-module` working copy both
resolve to `b57a1ed23cbe7d74afe98ed16ec8fd615db9d171`, annotated checkpoint
`m6-verified-current-state-audition-2026-10-24`. Continue in
`/Users/fuzz/Developer/LevelScope`; preserve the separate VideoImplementation
working copy, all ten existing checkpoint tags, installed apps/plugins and frozen
test packages. Fetch again after approval and review any newer main changes.

Qualify one named eight-channel speaker arrangement through the existing offline
Source/Current, measured evaluation, global trim, diagnostics and audition paths.
Keep three independently reported evidence levels: offline, injected playback,
and physical eight-channel listening. Missing eight-output hardware must not block
offline qualification or be disguised as a listening pass.

Implementation findings at the proposal baseline, before M7 changes:

- `StandaloneChannelLayout.h` already provides twelve-channel capacity, distinct
  side/rear/height role choices, original-file-indexed roles and explicit output
  routes. Duplicate outputs are rejected; unmapped routes start muted.
- `OfflineMeasurement.h/.cpp` only accepts the M1 mono/stereo/5.0/5.1 sets and has
  one surround pair. `StandaloneMeasurementService.cpp` aliases generic and side
  surrounds for M1, but rejects rear roles. This mapping cannot describe M7.
- `StandaloneMedia.cpp` uses the reader's `AudioChannelSet` for suggestions. Pinned
  JUCE 7.0.12 supplies count-based defaults: an unlabelled eight-channel WAV can
  suggest SDDS, whereas its generic eight-channel set uses side/rear surrounds.
  WAVE back-channel mask bits also use JUCE's generic surround enum values.
  These suggestions are not a reliable ordered 7.1 interpretation by themselves.
- `Core/Proposal/AnalysisEngine.cpp` explicitly rejects more than six peak entries.
  Replace this diagnostic storage restriction deliberately; do not simply open
  measurement and leave eight-channel diagnostics broken.
- `StaticGainProcessor`, verified-source preparation, FIFO and captured samples
  already use channel count. Current audition delegates capability validation to
  offline measurement. M2/M3 already bind measured identities and full precision.

## Exact layout and file interpretation

Recommended capability ID: `pcm-7.1-bs2051-system-i-v1`. It denotes BS.2051-3
Sound System I, 0+7+0 plus one LFE: three front, two side and two rear speakers,
with no height, wide or extra centre channels. Nominal positions identify the
measurement interpretation; M7 does not measure the user's room geometry.

| Canonical role | Nominal azimuth / elevation | Accepted explicit labels | Loudness energy coefficient |
| --- | --- | --- | ---: |
| L | +30° / 0° | L, FL, Left, Front Left, M+030 | 1.00 |
| R | −30° / 0° | R, FR, Right, Front Right, M−030 | 1.00 |
| C | 0° / 0° | C, FC, Centre, Center, Front Centre/Center, M+000 | 1.00 |
| LFE | not position-weighted | LFE, LFE1, Low Frequency Effects | 0 |
| Lss | +90° / 0° | Lss, SL, Side Left, Left Side Surround, M+090 | 1.41 |
| Rss | −90° / 0° | Rss, SR, Side Right, Right Side Surround, M−090 | 1.41 |
| Lrs | +135° / 0° | Lrs, BL, Back Left, Rear Left, Left Rear Surround, M+135 | 1.00 |
| Rrs | −135° / 0° | Rrs, BR, Back Right, Rear Right, Right Rear Surround, M−135 | 1.00 |

These are a finite interpretation table, not fuzzy guesses from filenames. Match
aliases case-insensitively after trimming surrounding whitespace and collapsing
internal whitespace. ASCII minus in ITU labels is equivalent to the printed minus;
other punctuation is not removed or guessed. JUCE `left/rightSurroundSide` map to Lss/Rss;
`left/rightSurroundRear` map to Lrs/Rrs. Generic Ls/Rs, “Surround”, and generic JUCE
surround roles are ambiguous for eight-channel manual assignments and require
explicit side/rear resolution. They retain their qualified meaning in M1 5.0/5.1.
The WAVE BACK bits are a format-specific exception with explicit provenance:
they mean Lrs/Rrs in the complete supported 7.1 mask, despite JUCE's enum names.

The arrangement and position names follow [BS.2051-3, Table 11](https://www.itu.int/dms_pubrec/itu-r/rec/bs/R-REC-BS.2051-3-202205-I!!PDF-E.pdf).
Coefficients follow [BS.1770-5, Annex 3, Tables 4–5, System I](https://www.itu.int/dms_pubrec/itu-r/rec/bs/R-REC-BS.1770-5-202311-I!!PDF-E.pdf).
In particular, **rear channels are 1.00, not the side pair's 1.41**. Use the
standard's named-layout table; do not infer a different coefficient from a room
placement tolerance or apply these coefficients as playback gains.

There is no mandatory internal PCM reorder. Every buffer, peak location, identity,
route and UI row continues to refer to the original decoded file-channel index.
Accept any explicitly confirmed permutation of exactly the eight canonical roles.
Provide two clearly named order presets, which fill roles and still require Confirm:

| Interpretation | File channels 1–8 |
| --- | --- |
| WAVE mask / standard FLAC order | L, R, C, LFE, Lrs, Rrs, Lss, Rss |
| Explicit side-before-rear order | L, R, C, LFE, Lss, Rss, Lrs, Rrs |

For WAVEFORMATEXTENSIBLE, accept the complete `0x0000063f` speaker mask and its
ascending-bit order. Do not mistake the SDDS/front-centre variant for this layout.
[Microsoft's WAVEFORMATEXTENSIBLE definition](https://learn.microsoft.com/en-us/windows-hardware/drivers/ddi/ksmedia/ns-ksmedia-waveformatextensible)
defines the relevant mask bits and ordering. A validated native FLAC default
channel assignment uses the first order; inspect any override channel mask before
making that suggestion. Only the equivalent complete mask is supported for this
M7 layout. See [RFC 9639 §§8.6.2 and 9.1.3](https://www.rfc-editor.org/rfc/rfc9639.html#section-9.1.3).
This relies on the format-defined assignment plus user confirmation, not a generic
channel-count guess.

For unlabelled WAV and AIFF, leave eight-channel positions unconfirmed and require
an explicit preset or individual assignment. Do not advertise arbitrary AIFF CHAN
metadata interpretation; unsupported metadata must be identified as such. A
reported metadata/manual-assignment conflict requires an explicit, recorded manual
interpretation; malformed or channel-count-inconsistent masks cannot silently
fall back to a qualified layout. An explicit unsupported layout (SDDS, wide,
7.0, 5.1.2, 7.1.4, two LFE, duplicate or missing roles) never acquires M7 support
just because it contains eight samples per frame. Broader existing Source import
and audition remain separate capabilities.

Bounded format/rate matrix: existing interleaved WAV PCM16/24/32 and float32,
AIFF PCM16/24/32, native FLAC PCM16/24 with the existing integrity requirements;
44100, 48000, 96000 and 192000 Hz. Preserve other existing Source behavior without
new qualification claims. No new decoder or sample-rate conversion. Use a native
FLAC fixture for the three-hour test rather than making a new large-WAV/AIFF
container-support claim as part of M7.

## Measurement and profile semantics

Extend the standalone worker's existing BS.1770-5 algorithm through Annex 3's
channel summation. Keep K filtering, gate rules, window cadence, EOF conventions
and true-peak implementation unchanged. Measure unweighted sample/true peaks and
their locations independently for all eight channels, including LFE; the programme
peak is their maximum. An LFE-only programme still has peaks but no usable
programme loudness. Zero-valued/missing loudness must not be invented.

Retain M/S and LRA behavior from [Tech 3341-2023](https://tech.ebu.ch/docs/tech/tech3341.pdf)
and [Tech 3342-2023](https://tech.ebu.ch/docs/tech/tech3342.pdf), including short-programme
interpretability and silence statuses. Loudness coefficients never alter the PCM.
Static gain remains the same full-precision operation on every channel, including
LFE, with no automatic bass management, +10 dB LFE gain, limiting or downmix.

**Recommendation: the existing EBU profile can evaluate the selected scalar rules
for this qualified layout**, when the user confirms its existing pre-produced,
whole-programme linear-audio applicability. This is the implementation conclusion
from R128's BS.1770 reference and the Annex 3 channel capability, not an assertion
that any broadcaster accepts a particular eight-channel deliverable.
[R128-2023 clauses h–n](https://tech.ebu.ch/docs/r/r128.pdf) support the unchanged
−23 LUFS target, ±0.2 LU QC allowance and strict −1 dBTP ceiling. Its separate
±0.3 dB meter provision for 20 kHz-band-limited signals does not relax that ceiling.
No universal LRA limit is added.

Keep `ebu-r128-2023-programme-linear-qc` definition `1.0.0` and its digest unchanged
because its rules/applicability are unchanged. Add qualified layout/method support
in the measurement capability layer, not an eight-output prerequisite or a new
delivery rule. Bump the measurement implementation/mapping identity and require
compatible fresh snapshots. If implementation actually changes any profile
definition field, it must receive a new definition version instead of changing
the existing version in place.

Use the current narrow result wording: “Source/Current State evaluation” and
“Passes selected programme loudness and true-peak rules.” Metadata, full delivery,
monitor setup and encoded-export verification remain separate. A valid eight-channel
measurement can exist with no profile selected. Inapplicable/invalid profiles,
missing/non-finite metrics, unsupported layouts and stale snapshots cannot pass.

## Architecture and required changes

1. Add a small, explicit, versioned layout interpretation/capability descriptor.
   Resolve format suggestions and explicit aliases to ordered speaker positions;
   use separate legacy M1 and System I tables. Reuse it for Source/Current requests,
   audition validation and displayed capability. Keep generic arrays and buffers
   sized by count within the existing twelve-channel limit. Do not import a general
   arbitrary-speaker editor or enable twelve-channel measurement.
2. Extend the measurement role representation with distinct side/rear semantics,
   support the exact eight-role set at the four rates, and keep all eight peak
   vectors and locations. Version the changed mapping/measurement contract.
   Preserve M1 numerical behavior; no plugin detector/filter changes.
3. Bind source hash, ordered roles, layout-definition version, range, decoder/build,
   method and exact processing plan/revision into the existing snapshot/audition
   identities. Role edits invalidate both Source and Current interpretation/results;
   gain edits invalidate Current only. Reject late completions and mismatched
   Source/Current pairs even when individually valid. Profile changes only
   reevaluate compatible measurements; monitor/routing/video changes do not
   change measurement identity. Audition routing identity remains separately bound.
4. Extend M2/M3 integration tests and capability diagnostics. Retain deterministic
   trim policy, exact boundaries, prediction-only drafts, stopped explicit Apply
   and real post-Apply measurement. A peak on LFE or either rear channel can constrain
   the proposal. Already-pass is a no-op; failed/cancelled verification retains
   applied gain and shows unverified/actual measured failure, without retry/rollback.
5. Extend M6 audition to the supported layout, keeping the shared processor and
   verified-copy reuse. Source/Current, processing, roles and routing remain
   stopped-only; monitor volume remains smoothly adjustable. Keep whole-callback
   overload silence, no presentation advance for the rejected block, latched stop
   and explicit restart with reset buffers/ramp. No physical device is needed for
   preparing a source or computing Current State.
6. Remove diagnostics' six-channel storage assumption using the validated count
   and role binding. Retain all eight channel maxima and source-time locations in
   M4 timeline/peak views and M5 section summaries. Keep full-resolution calculations,
   extrema-preserving plots, separate LRA distribution, explicit padding and
   read-only navigation. Add no section processing or per-channel loudness steering.
7. Minimal UI: explicit supported/unsupported layout, per-row canonical role and
   interpretation origin, `Channel 1 … Channel 8` labels, clear Source/Current
   gain/freshness, and monitoring coverage. Fix the deferred interleaved “File N”
   wording here; separate-file import remains excluded.

## Three evidence levels and monitoring behavior

| Level / capability label after passing | Required evidence | Needs eight physical outputs? |
| --- | --- | --- |
| **7.1 offline measurement and processing — qualified** | Actual two-platform import, role confirmation, Source/Current measurement, M2/M3, diagnostics and lifecycle with device closed | No |
| **7.1 audition/routing — automated verification passed** | Actual two-platform injected-device runs; all eight pre-monitor streams and routed callback outputs checked independently | No |
| **7.1 physical routing and listening — qualified on [named setup]** | Coordinated, simultaneous eight-output routing to an identified 7.1 monitoring chain, channel identification and listening on each claimed host/setup | Yes |

These labels must not be combined into an unconditional “7.1 fully qualified”.
The release report lists each level's pass/pending/fail separately; product UI
shows current measurement/monitoring capability and coverage, not a stored claim
that an arbitrary newly connected device was listening-qualified.

With fewer outputs, permit explicit channel selection using the existing one-to-one
routes and explicit Muted rows. Show, for example, **“Partial monitoring: L, R;
6 of 8 channels muted. No downmix.”** This is the isolated L/R content, not a stereo
mix of the programme. The same applies to sequential checks of centre, sides,
rears or LFE through selected outputs. Do not automatically create those routes,
reassign invalid routes, fold channels together or call partial monitoring 7.1
listening. A device change invalidates or requires review of the output map; it
does not invalidate offline measurements. Wrong-rate devices fail explicitly.

All eight channels continue through processing and pre-monitor capture even when
some are deliberately muted at the output stage. Full electrical coverage requires
eight distinct valid routes; physical speaker correspondence still needs level 3.
No LFE redirection/crossover/calibration is supplied by this app. All-muted
instrumented timing runs are labelled as such; preserve the product's existing
all-muted Play restriction rather than silently changing it for qualification.

Pending physical checks must identify host/OS, interface, exact output-to-speaker
map, actual sample rate, monitoring chain and results. Eight advertised interface
outputs or a successful stereo check alone are insufficient. Electrical capture of
eight outputs can qualify that routing test, but cannot substitute for an identified
seven-speaker-plus-LFE listening setup. No existing device or audible authorization
is assumed for M7; coordinate opening and sound separately before each agreed group.

## Acceptance tests

All proposed limits below are engineering acceptance criteria, separate from
profile allowances. Preserve per-case official reference limits unchanged.

| Area | Required acceptance |
| --- | --- |
| Layout/order | Unique channel-ID fixtures verify raw decoded order for every supported format. Verify complete WAVE mask, FLAC default/override and manual orders; all 40,320 role permutations at the mapping layer plus representative PCM permutations. Reject duplicate, unassigned, mixed side/rear aliases, wide/SDDS/height/extra-LFE/count-inconsistent cases. No count-only qualification. |
| Independent loudness oracle | Isolate each role at every rate; compare front/rear energy 1, side 1.41 and LFE exclusion using a separate oracle, not the production weight helper. Steady side-minus-rear result equals `10*log10(1.41)` within 1e-9 LU. An equal-energy seven-main-channel sum uses coefficient sum 7.82. Exercise uneven channels, silence, gate crossings, short programmes, chunking and repeated runs. |
| All-channel peaks | Place the unique maximum on each of eight channels, including LFE and rears; verify programme maximum, original index, sample-frame and true-peak fractional location. Include impulses at both boundaries, intersample peaks, ties, missing/silent values and non-finite rejection. Independently generated accuracy cases use existing M1 limits. |
| Device independence | Run the full import → Source measurement → Current gain → remeasurement → M2/M3 → diagnostics flow with no device selected and with an injected backend that fails/counts every open attempt: zero opens. Offline work succeeds even if the available device has zero/two outputs or cannot run the source rate. |
| Compliance/proposals | Exact thresholds and adjacent representable values; side/rear/LFE peak-limited and impossible trim; no optional-LRA blockage; stale/mismatched content/order/range/plan rejection; no-change preservation; cancelled/failed verification behavior. Same profile definition/digest; no test budget widens its limits. |
| Sample parity | Source/Current zero and positive/negative/non-rounded gains across all rates/formats; LFE gain included. Same-platform offline versus actual pre-monitor audition samples must be bit-exact on all eight channels. Independent long-double reference within one float32 ULP. Cross-platform zero gain bit-exact, finite nonzero at most one ULP; report actual maxima. |
| Injected routing | Per-channel signatures into 8- and 12-output devices, arbitrary/non-contiguous maps, unassigned output silence, duplicates/out-of-range rejection. On two outputs explicitly route a selected pair and prove all other channels remain in pre-monitor capture and intentionally muted post-route. Zero-output/rate mismatch/loss/hot replacement never changes roles or silently falls back. |
| Audition lifecycle | Seeks at non-block boundaries/final sample, stopped switching and Apply, content mutation, stale workers, EOF/restart, device loss, cancellation and close during work. Overload on each mapped channel silences the entire callback; reducing monitor and explicit restart resumes correctly without stale ramp/FIFO data. No hidden clipping, attenuation or resampler. |
| Identity/UI | Role edits invalidate affected measured results; profiles reuse compatible measurements; gain affects Current only. Monitor, explicit partial routing, section annotations and video offset preserve measured identities. Show UNVERIFIED before remeasurement and actual FAIL distinctly. Channel 7/8 cannot disappear in scroll/peak/section views. |
| Diagnostics | All eight peak owners/locations, exact audio-source interval labels, EOF padding, unavailable/LRA stability, rankings and full-resolution statistics. Mixed gate inclusion/extrema survive plot reduction. Navigation does not play or process audio. |

Retain cross-platform metric budgets of 0.01 LU/dB and 0.1 LU LRA, and existing
exact status/frame identity rules. PCM/role permutations should preserve aggregate
metrics within 1e-9 LU/dB for deliberately non-boundary fixtures; exact reordering
of per-channel peaks is required. Numerical accumulation order is documented;
near-gate/boundary behavior is tested separately, never hidden by that tolerance.

### Long programmes and resource limits

On each actual Mac, use one real three-hour eight-channel 48 kHz native FLAC
fixture through validation, verified preparation, Source and nonzero Current
measurement, final-second inspection and random seek/capture. The fixture must
contain distinguishable channels, changing levels and late peaks; record encoded
size, decoded frame count, hashes and temporary disk demand. Shorter fixtures
exercise all four rates. A three-hour/192 kHz measurement-sized snapshot with 1,000
sections exercises downstream time/index/storage limits without claiming that it
is a second complete PCM throughput run.

Retain M4/M5 budgets for the eight-channel snapshot: diagnostic core/worker <=2 s,
viewport p95 <=50 ms, paint/edit/refresh p95 <=100 ms, diagnostic storage <=64 MiB,
diagnostic extra peak RSS <=128 MiB, section work <=2 s and <=16 MiB additional
workspace, active diagnostic/section cancellation <=1 s. Audition/preparation
additional owned buffer/workspace remains <=16 MiB, independent of duration;
cooperative storage cancellation <=2 s on responsive storage. Record full import,
measurement and preparation elapsed times, with no new faster-than-real-time claim.
Inject insufficient-space/create/write/flush failures, clean owned partials and
retain the prior stopped valid session. Do not fill a real disk. Abrupt-termination
temporary-file limitations remain documented.

## Two-platform and video qualification sequence

1. After approval, stabilize in a new M7 build directory on native arm64 Sequoia,
   retaining pinned JUCE 7.0.12, universal arm64/x86_64, macOS minimum 10.13 and the
   established floating-point configuration. Build only named Standalone,
   Measurement, Scaffold, Compliance, Proposal, Diagnostics, Sections, Audition,
   VideoTests, VideoQualification and a named M7 test target. No default/all,
   VST3, installer or installed-app replacement.
2. Run all M1–M6/runner regressions and both pinned official EBU V05 and ITU BS.2217
   suites because shared measurement/import interpretation changes. Compare M1
   behavior against its qualified evidence. New eight-channel fixtures are
   independently generated extension tests; embedding official vectors into eight
   channels is supplemental evidence, not a new official EBU/ITU 7.1 test suite.
3. Freeze candidate source/build/binary/vector manifests. Inspect Golem's available
   dependencies and run that same candidate natively on actual Intel High Sierra
   in a separate folder. Intel-on-Sequoia is not a substitute. Rerun affected tests
   on both Macs after fixes; preserve failures and older candidates.
4. Complete silent UI workflow on both Macs, device closed. Verify offline support,
   eight distinct roles/peaks, M2/M3, unverified Current, partial-monitor labels,
   source replacement/clear, cancellation, source-time diagnostics and video offset
   neutrality. Automated injected-device checks supply level 2 evidence separately.
5. Run existing native H.264 CFR/VFR and ProRes profile/metadata/lifecycle regressions,
   now also with eight-channel Source and nonzero Current captures. Preserve rational
   signed offsets, presentation clock, seeks/generations, stop/resume/restart and
   EOF/duration. Embedded movie audio remains undecoded and unrouted.
6. Coordinate physical, all-muted timing runs on available built-in outputs; an
   eight-output interface is unnecessary. Exercise the existing short ProRes
   profile/rate matrix and H.264 CFR/VFR in both domains, then three 30-minute
   ProRes Current runs (Proxy/LT/422) per Mac with eight channels processed upstream
   of monitoring. Use agreed source/device rates, e.g. 48 kHz Sequoia and 44.1 kHz
   High Sierra if available. Preserve document 25's paint/clock/underrun gates,
   queue <=8, peak RSS <=768 MiB and both growth limits <=32 MiB with the accepted
   memory method. These qualify timing under the recorded device configuration;
   they do not qualify spatial listening or physical optical/acoustic sync.
7. Schedule physical eight-channel mapping/listening only when the user has an
   appropriate setup and authorizes it. Verify each channel reaches its explicitly
   named speaker/LFE path, then Source versus nonzero Current, clean transport and
   monitor-volume behavior. Record each actual Mac/setup/rate separately. Stereo,
   sequential subset or all-muted runs never close this evidence level.

## Completion, dependencies and exclusions

M7's recommended mandatory release scope is **offline 7.1 plus automated audition/
routing and the relevant video/UI regressions on both actual Macs**. Those must all
pass for that scope to be called qualified. If physical timing is unavailable,
report it separately and withhold the combined audition/video qualification;
completed offline evidence remains valid. Missing eight-channel listening alone
does not block the mandatory scope: publish “Physical 7.1 listening: pending” with
the exact untested setup/rates. A failed mandatory test remains a failure, not a
hardware-pending exception.

Dependencies are the accepted M6 baseline, pinned JUCE/toolchain, primary mapping
tables above, deterministic extension fixtures, actual High Sierra access and
enough temporary storage. Eight-output monitoring hardware is a dependency only
of the separately claimed physical evidence. No hardware purchase is assumed.

Preserve twelve-channel/height-role architecture for a later 7.1.4 milestone, but
do not enable it. Its position mapping, weights, profile applicability, all-channel
peaks, diagnostics, routing and two-platform qualification remain separate. FLAC's
eight-channel limit also means future twelve-channel PCM requires an appropriate
container rather than extending FLAC by assumption.

RF64 and other large-file container extensions, Dolby Atmos objects/rendering/ADM import, separate-file import, section processing,
dynamics, automatic correction, export, persistence, new video features and VST3
measurement parity remain excluded. No installed plugin behavior/state changes,
system audio configuration changes, Pro40 repository/driver/service changes, or
alteration of preserved video work are part of M7.

## Approved decisions

- Select exactly BS.2051 System I side-plus-rear 7.1; other eight-channel layouts
  stay unsupported for M7 measurement/Current State.
- Require explicit order/role confirmation, with the two named presets and bounded
  format interpretations above; include the small “Channel N” wording correction.
- Keep EBU profile 1.0.0 unchanged, while versioning the expanded measurement/layout
  capability and preserving narrow evaluation language.
- Accept mandatory offline/automated/video qualification independently of full
  physical 7.1 listening. Record physical listening as pending until an actual
  eight-output, correctly connected monitoring setup is available and coordinated.

No decision about purchasing or configuring monitoring hardware is needed to
approve implementation of the mandatory scope. The user approved implementation
and qualification, then separately authorized checkpoint publication and integration;
the completed publication is recorded in document 32.

## Additional approved regression requirements

- Concrete WAV and native FLAC 5.1 fixtures with conventional back/surround metadata
  or decoder labels retain their established ordered roles and 1.41 surround weight.
  The new 1.00 rear coefficient is restricted to explicitly identified System I 7.1.
  Compare their measured results with the unchanged M1 role interpretation.
- Separate supported 7.1 fixtures place the unique highest true peak on a rear
  channel and on LFE. Explicitly mute that channel for partial monitoring: programme
  measurement, M2 verdict and M3 peak-constrained proposal must remain identical.
  All eight channels must be present in processing and pre-monitor capture.
- Muted physical video/device checks are authorized on Mac Studio Speakers at
  48 kHz and Golem Built-in Output at 44.1 kHz. All routes remain muted. This does
  not authorize or establish full physical 7.1 listening qualification.

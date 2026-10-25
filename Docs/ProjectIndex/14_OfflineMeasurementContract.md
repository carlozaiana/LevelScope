# First offline measurement milestone

Decision date: 15 October 2026. Development: `codex/logic-module`, from
`standalone-scaffold-2026-10-15`. The plugin and both repository checkpoints are preserved.

## Scope and standards

M1 adds a worker-only, reusable measurement kernel and source/current-state analysis.
Current State supports identity and explicit static gain, applied to all file channels
including LFE. This trim is distinct from dynamics detector/apply masks and monitor gain.
An internal render/verification harness exercises encoded output; public delivery export,
proposals, dialogue detection and full dynamics integration remain subsequent work.

Pinned measurement references: ITU-R BS.1770-5 (November 2023), EBU Tech 3341 v4
(November 2023), EBU Tech 3342 v4 (November 2023). Delivery profiles are a separate
versioned layer. Placeholder targets never authorize compliance or export.
ATSC A/85:2026-07 requires separate dialogue measurement for long-form content; ordinary
level-gated full-mix loudness cannot be relabelled as that measurement.

Initial measurement layouts: mono C/L/R, stereo L/R, 5.0/5.1 L/R/C/Ls/Rs with optional
LFE. Surround and surround-side pairs are explicit accepted aliases for conventional
5.1; those legacy decoder labels keep their 1.41 coefficient, including WAVE back
bits and native FLAC's conventional 5.1 assignment. File
order is retained. L/R/C energy weights are 1; Ls/Rs are 1.41; LFE is 0 for loudness.
All channels enter sample/true peak. Rates: 44100, 48000, 96000, 192000 Hz; reject
other rates rather than substitute a fallback. Import support remains broader.

M7 implementation under qualification adds only explicitly identified BS.2051
System I 7.1, capability `pcm-7.1-bs2051-system-i-v1`, role mapping
`ordered-pcm-roles-m7.1` and measurement method `bs1770-5-ebu2023-m7.1`.
L/R/C and the distinct M±135 rear pair use 1.00, M±090 sides use 1.41, LFE uses
0 for loudness. All eight channels enter peaks and the unchanged full-precision
static-gain operation, regardless of monitor routes. The original file order is
retained; channel count alone never establishes System I. See the exact metadata,
alias and confirmation rules in [31](31_M7ChannelBasedPcm71Proposal.md).
Height layouts and RF64/BW64 capability extensions remain outside this scope.
The EBU profile's definition and digest are unchanged; evaluation requires the
compatible versioned measurement. Qualification status belongs in the M7 report,
separately for offline processing, injected routing and physical 7.1 listening.

## Numerical and boundary contract

- New standalone kernel; do not change the shared plugin's existing filter or stats.
- Double-precision K filtering/energy accumulation with fixed sample traversal. The
  48 kHz filter uses the published coefficients. Other rates use a response-preserving
  bilinear design; independent frequency fixtures qualify those rates.
- Source frame zero begins a new measurement with zero filter state. No host transport
  guards, start ramps, output routing or monitor gain enter the analysis.
- Integrated: complete 400 ms blocks every 100 ms, strict BS.1770 gate comparisons,
  -70 LUFS absolute gate, relative gate 10 LU below absolute-gated block energy mean.
  Incomplete final blocks are discarded; display histogram bins never determine gates.
- M/S: rectangular 400 ms / 3 s windows. Publish curves every 100 ms; calculate maxima
  at every sample so curve decimation cannot miss offset-sensitive peaks. Only complete
  source windows enter maxima. Curves carry source-frame end positions.
- LRA: complete 3 s windows every 100 ms, including a separate 1.5 s zero-input drain
  after source EOF following Tech 3342's file-based reference note. This drain does not
  change source duration, integrated loudness, source M/S maxima or sample peak.
  LRA uses its own absolute-gated short-term energy mean and -20 LU relative gate;
  percentile indices follow Tech 3342's reference rounding. Values for programs below
  60 s are marked unstable. No eligible population means unavailable, never numeric zero.
- True peak: independent interpolation on unweighted samples, all channels. Zero
  extension at both file boundaries; drain FIR state at EOF. Report method/version.
  Sample peak remains separate. Digital silence has zero linear peak and no finite dB
  value; expose a silence reason. Reject non-finite samples and accumulator overflow.
- Full-program only in M1. Segment requests fail explicitly until context rules exist.

## Identity and lifecycle

Copy source bytes to a uniquely owned temporary snapshot with bounded buffers and
cooperative cancellation. Verify the snapshot SHA-256 against imported identity before
decoding; validate its metadata and payload. Workers decode only that private snapshot.
Original path/stamp are diagnostic, not content identity. Normal/error/cancel teardown
deletes owned temporary files. Abrupt process termination may leave temporary files.

Requests/results bind source hash, ordered roles, range, typed processing plan digest,
measurement version, decoder/build/platform identifiers and session/job ticket. Invalid
or stale results cannot replace current truth. Processing edits invalidate measured
Current State. Profiles do not alter base measurements but invalidate verdict/approval.
Predictions have no route into measured state. No on-disk session schema is frozen.

Render verification measures a reopened encoded artifact. Its identity includes output
SHA-256 and format. Any future approval must match session revision, measurement/plan
identity, profile edition and output hash. Output paths must be distinct from input,
including links/aliases; M1 internal rendering only uses owned temporary artifacts.

## Acceptance and evidence

Official vectors: EBU test set v5 and ITU BS.2217. Keep downloaded audio outside Git
and outside shipped artifacts; record URL, archive/file hashes and expected values.
Synthetic tests are independently generated from published signal definitions and
must be labelled separately from downloaded official vectors. Missing external vectors
are reported as outstanding evidence, not silently counted as passed.

Published tolerances: applicable loudness tests ±0.1 LU; Tech 3342 LRA ±1 LU;
Tech 3341 true-peak cases 15–23 +0.2/-0.4 dB. Also test exact boundaries, silence,
short files, LFE-only, channel permutations, chunk independence, repeated runs,
identity/static-gain parity, encoded-output remeasurement, mutation/cancellation and
stale-ticket rejection. Cross-platform target: 0.01 LU/dB, 0.1 LU LRA; qualification
requires actual suite runs on both arm64 Sequoia and Intel High Sierra.

Qualification update (15 October 2026): the supplied EBU v05 archive now passes all
66 minimum-requirements files on arm64 Sequoia (72 metric rules, 1,052 readings), plus
four separately labelled supplementary sanity checks. Official files are 48 kHz;
generated multirate checks remain separate. Cases 11/14 validate successive kernel
maxima, without claiming live display/scheduler behavior. See the manifest, runner and
saved report under `Tests/Reference` and the updated handover for exact evidence.
Final qualification update (17 October 2026): actual Intel High Sierra 10.13.6 and
arm64 Sequoia 15.7.9 both pass the 184 measurement checks, 205 scaffold checks,
eight runner tests, both pinned official suites and the silent interactive UI checklist.
Cross-platform comparisons have zero errors within the budgets above. See
[17_M1QualificationResults.md](17_M1QualificationResults.md) for exact identities,
results and scope. The measurement checkpoint is `measurement-milestone-2026-10-17`.

## Video and transport contract (design only)

Use a persistent session UUID, separate audio-content identity and reference-video
identity (hash/path/stream ID), and a monotonically increasing transport generation.
Audio timeline positions are signed 64-bit sample frames at the source rate. Store
video PTS/timebase and nominal frame rate as reduced integer rationals; preserve VFR
timestamps rather than infer them from frame count. Timecode metadata includes start
frame/timecode origin and an explicit drop-frame flag; drop-frame affects labels only.

Define offset as `video presentation time = audio timeline time - videoStartOnAudio`:
a positive offset places video start later on the audio timeline. Perform conversions
with checked rational arithmetic and an explicit rounding rule at display boundaries.
Selected audio is authoritative; embedded video audio is muted. Playing is driven by
presented audio frames plus measured device latency, not worker decode position.
Stop retains position; restart seeks to zero; seek flushes queued audio/video, bumps
generation, prerolls and resumes only if previously playing. Audio EOF stops transport;
outside video coverage show no frame. Decoder work/display run outside the RT callback.
Reference-video attachment/offset does not change audio measurement identity.

Video implementation follows M1 and precedes polished UI. Test fractional rates,
drop-frame minute/ten-minute boundaries, signed offsets, overflow, seek generations,
VFR and known sync points before claiming frame accuracy.

## Primary references

- https://www.itu.int/rec/R-REC-BS.1770-5-202311-I/en
- https://tech.ebu.ch/docs/tech/tech3341.pdf
- https://tech.ebu.ch/docs/tech/tech3342.pdf
- https://tech.ebu.ch/docs/r/r128.pdf
- https://www.atsc.org/wp-content/uploads/2026/07/A85-2026-07.pdf
- https://tech.ebu.ch/publications/ebu_loudness_test_set
- https://www.itu.int/oth/R1102000001/en

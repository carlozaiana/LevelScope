# M9 separate-file programme import — implementation and qualification

Status: **candidate-04 qualified on actual arm64 Sequoia and Intel High Sierra
on 25 October 2026; published and integrated into main.**
The annotated checkpoint is
[`m9-separate-file-stereo-5.1-2026-10-25`](https://github.com/carlozaiana/LevelScope/tree/m9-separate-file-stereo-5.1-2026-10-25).
Its peeled commit is the accepted M9 checkpoint and main revision at publication.
The annotation records the full commit SHA, qualified source manifest, original
binary hashes and evidence identity. The [publication record](../../Tests/Reference/results/m9-publication-2026-10-25/README.md)
documents that mapping and the preservation checks.
The approved contract is [37](37_M9SeparateFileStereo51Import.md). This work starts
from published M8 `8e984889c5e73670f370c6c568c8ea1d7902a35c` on the established
`codex/logic-module` working copy. The M8 tag, tested universal package and earlier
evidence are preserved.

## Implemented boundaries

`ProgrammeAssembly` represents an ordered set of two mono stereo members or six
conventional 5.1 members. The programme order is L/R or L/R/C/LFE/Ls/Rs. Picker
order is not an assignment. The import staging screen requires explicit roles and
confirmation that these files represent the same aligned mix. Every member must
be mono, at an existing qualified rate, with exactly the same decoded frame count.
Sample zero in every member is programme sample zero; metadata/timecode offsets
are not applied. No resampling, summing, padding, truncation or alignment estimation
is performed. Distinct identical-content files are permitted; repeated filesystem
identities, including symlinks and hardlinks, are rejected.

The version-2 assembly digest covers member content hashes and decoded format,
ordered roles/assignments and zero-frame offsets. It excludes locators, display
names and persistent slot IDs. Alignment confirmation stores that exact digest.
Changing content, assignment or alignment invalidates confirmation. An unchanged
reopen or complete same-content relink retains it. Equal names, lengths or rates
never establish alignment. A different-content relink is rejected; use a new staged
import and renewed alignment confirmation to replace programme content.

`VerifiedSource` owns an all-or-nothing group of private copies. The aggregate
free-space preflight includes reserve space and every write/flush is checked.
Failure/cancellation destroys only owned partial copies. `ProgrammeReader` reads
bounded blocks from those files in programme order. Measurement, the internal
render verification adapter and audition use this same worker-only reader and the
existing full-precision, all-channel static-gain processor, including LFE. There
is no extra full-duration interleaved temporary. Prepared groups are reused across
seeks and stopped Source/Current/gain changes. The real-time callback does no file
access or decoding. The existing twelve-channel interfaces and interleaved 7.1
support remain; the separate-file importer is deliberately two/six channels only.

Measurement identities explicitly include the assembly capability. Numerically
identical interleaved audio keeps different provenance and cannot silently lend
its caches to a group. Group member changes during preparation, measurement or
adoption reject obsolete results. Monitoring and video remain outside programme
measurement identity. Predictions still cannot satisfy measured compliance rules.
The EBU target definition/digest and numerical measurement method are unchanged.

The import window stages files before replacing the live programme. Failure or
cancellation retains the previous session stopped. Successful replacement clears
old measurements, Current State and actionable proposals and changes programme
identity. Mapping/routing and per-channel diagnostics display filenames alongside
roles. Existing interleaved inputs keep their Channel 1/Channel 2 labels.

## Sessions and recovery

The existing checksum container is retained. M9 reads M8 core/assembly schema 1
and writes schema 2. Opening an M8 document migrates it in memory only. Save warns
that the new format requires M9, and offers Save As with an M9-suffixed suggestion
so the original M8-readable file can remain intact. Older M8 readers reject the
new core. Conservative build/platform cache compatibility remains in force; M8
measurements require remeasurement in the changed build.

Every group member is verified before an opened programme becomes usable. An
unresolved member makes the whole group incomplete: no partial playback, silent
replacement channel or inherited pass. Settings/sections remain recoverable.
Relinking identifies a stable file slot and verifies exact content/format. Partial
relinks can be saved while the programme remains incomplete. Dormant replaceable
group caches cannot be evaluated while unresolved; only a complete verified group
and compatible build/platform can restore them. Corrupt independent caches are
discarded without losing a valid core. Session writes reject aliases of **any**
member or reference picture. Persistent document identity remains separate from
runtime opening/job identity, and coherent captured revisions preserve M8's dirty
state and late-completion rules.

## Qualification gates

The frozen universal package is
`build/m9-qualification-2026-10-25/candidate-04`. Its recorded base is the M8 commit
above; qualification ran before the M9 commit was created. All 187 source-manifest
entries match the saved checkpoint. The original package retains its historical
base/build identity; it has not been rebuilt or rewritten to insert the later
publication commit. The source manifest identifies the implementation under test.

| Identity | SHA-256 |
| --- | --- |
| Source manifest | `ecf819093af35debedfed260a28a972183cf8484af677f6fbca749f4dbea0d34` |
| Measurement/build identity | `6cdc517ca6618832430a7dbe678c37d3c459888674a112aa2643c2a25186f329` |
| Package identity.json | `374bfa3b5003d7a67651ccfea71b2a3ceceaf8edc17398573a4862ac82b27898` |
| Universal standalone executable | `dd15614194cb4bb40b86415ff6a93e0ddc0dd35455e5ae03adff4f4748bc6938` |

Candidate-01 remains preserved. Its UI checks found that a rejected replacement's
incomplete-group notice did not clearly distinguish the still-valid previous
session. Candidate-02 changes only that notice in `StandaloneSessionUI.cpp`:
“Relink rejected; previous session retained.” It does not change the reader, DSP,
measurement, session schema, identity algorithm or video behavior. The build
identity changes, so candidate-01 caches are conservatively discarded. Both
platforms rerun the automated suites against candidate-02; unchanged candidate-01
UI observations may be reused only with the recorded one-file equivalence, and
the changed notice and cache/reopening behavior are explicitly retested.

Candidate-02 is also preserved. Its physical runner failed **before device
opening** because the path-only hash helper received `__file__` as a string.
Candidate-03 changes only `Tests/run_m9_muted_video.py` to pass `Path(__file__)`.
All application/test binaries, compiled source inputs and fixtures are byte-identical
to candidate-02, with the same measurement/build identity; no rebuild substitutes
for the tested binaries. Candidate-02's source manifest is
`367f667b23e193ff4b233ca9c5a3db966304511f03bdb5caa584e5b4a17ffbf3` and its package
identity hash is `16ac81d6b4a740a32486bc5610841f02300ecbfa05e5c27baf21640229d00f2e`.
The explicit `candidate-03-equivalence.json` record authorizes reuse of unchanged
offline/UI evidence without rewriting original report identities. Physical runs
use candidate-03 and the corrected runner; the prior failure log remains intact.

Candidate-04 preserves candidate-03 and adds the full-resolution stable-LRA/silence
oracle matrix described below, its comparator and fail-closed runner coverage.
Only `LevelScopeProgrammeTests` is rebuilt. All application binaries, other test
executables, production source and fixtures remain byte-identical to candidate-03;
the measurement/build identity is unchanged. `candidate-04-equivalence.json`
lists every changed input/artifact. Candidate-02 regression/long-media/UI evidence
and candidate-03 physical evidence retain their original identities and apply to
the identical production implementation. The new test binary must pass the complete
M9 short suite and supplemental oracles on each actual platform. This is additional
coverage, not an application or playback fix.

Raw reports and work files remain under ignored
`build/m9-qualification-2026-10-25`. The final acceptance and evidence index are
`acceptance-candidate-04.json`, `status.json` and `evidence-manifest.json` in the
results directory below. Earlier sub-suite summaries retain their original
incomplete overall status; the aggregate acceptance record closes the separate
gates without rewriting historical results.

The user also confirmed the exact roles and sample-zero alignment of both Tutti
Uccidono stereo and 5.1 mono sets. Supplementary silent development checks passed
for both 48 kHz/24-bit, 48,808,914-frame programmes. The originals were hashed
before and after and remained unchanged. These checks supplement the generated
numerical oracles; they do not replace frozen-candidate qualification or establish
physical multichannel listening.

Required qualification gates (individual results below):

- Sequoia stabilization and frozen universal candidate; identical binaries on
  actual Intel High Sierra, with exact OS/CPU/build/artifact identities.
- M1–M8 regressions and pinned EBU/ITU references, since the reader adapter changed.
- All two/six-member picker permutations; complete numerical interleaved oracles;
  all-channel offline/audition parity, including LFE and muted highest peaks.
- Exact alignment binding, group failure/cleanup, stale identity rejection,
  explicit global trim and actual measured verification, diagnostics and sessions.
- M8 migration, new-format rejection by preserved M8 reader, cache recovery,
  bidirectional grouped exchange and identical-content per-slot relinking.
- Real three-hour six-mono 48 kHz Source/Current measurement, coherent save/reopen,
  1,000 sections, bounded memory, responsive persistence, active cancellation and
  prepared arbitrary/final-sample seeks.
- Silent UI qualification on both Macs.
- Ten short Source/Current H.264 CFR/VFR and ProRes Proxy/LT/422 cases, and one
  30-minute Current/ProRes run per Mac, using explicitly approved built-in devices
  and all routes muted. Announce each physical block before opening a device.

Completed candidate-02 evidence:

- M1–M8 regression blocks pass on actual arm64 Sequoia 15.7.9 and Intel High Sierra
  10.13.6. This includes pinned official EBU/ITU references and isolated optional
  paint diagnostics. Cross-platform M1–M5, 124-case M6, 218-case M7 and 48-case M8
  comparisons pass; the compared processed captures are bit-identical.
- Each platform's M9 suite passes 5,261 assertions, 56 Source/Current sample cases,
  ten native reference-picture cases and three preserved M8 migration cases. The
  56-case cross-platform comparison passes with bit-identical captures (zero ULP).
- Both directions of stereo/six-file session exchange pass: exact settings and
  alignment confirmation survive, incompatible caches are rejected, and both
  measured states can be rebuilt. Preserved M8 readers reject schema 2 on both
  actual platforms; these deliberately failing positive-exchange assertions remain
  raw historical results with separate expected-rejection assessments.
- Sequoia's silent UI gate passes, using the explicitly documented candidate-01
  equivalence for unchanged interactions plus candidate-02 reopening, migration,
  missing-member recovery, cancellation and corrected rejection notice checks.
- Sequoia's actual three-hour, six-file Source/Current run passes 1,041 assertions:
  verification 74.874 s; measurement 722.536/769.950 s; save/parse 94.45/455.36 ms;
  UI callback p95 7.58 ms; peak RSS growth 133,758,976 bytes; session 10,433,268
  bytes; cancellation 0.287 ms and preparation cancellation 0.942 ms; prepared
  seek/prefill maximum 6.52 ms. Owned payload equals the original files' summed
  1,014,119,716 bytes, with no extra full-duration interleaved copy.
- High Sierra's same actual three-hour programme also passes 1,041 assertions:
  verification 91.516 s; Source/Current measurement 1,299.615/1,343.717 s;
  save/parse 215.47/1,024.89 ms; UI callback p95 7.34 ms; RSS growth 82,939,904
  bytes; session 10,432,888 bytes; cancellation 0.100 ms and preparation
  cancellation 5.215 ms; prepared seek/prefill maximum 8.60 ms. Owned payload
  equals the original summed payload. All stated three-hour resource budgets pass.

Candidate-03's ten short physical Sequoia cases are accepted: Source/Current for
H.264 CFR, H.264 VFR and ProRes Proxy/LT/422. Source CFR, VFR and ProRes 422 each
recorded one actual painted-frame skip; the other seven cases recorded zero.
The strict zero-skip harness failures remain unchanged, alongside the separate
accepted-policy assessment. All routes were muted and audio configuration was
unchanged.

Sequoia's 30-minute Current/ProRes endurance run is accepted with two recorded
picture skips, a maximum painted interval of 50.378 ms and zero selected-frame
PTS skew. The strict harness failure is retained separately from the accepted
policy result. Audio continuity, identity, lifecycle and memory requirements pass;
system audio configuration is unchanged.

High Sierra's ten short cases and 30-minute Current/ProRes endurance run are also
accepted, with zero recorded skips, zero underrun blocks, zero selected-frame PTS
skew and passing memory/lifecycle gates. The endurance maximum painted interval
is 34.685 ms. Both Macs' audio configurations are unchanged. The 22-case summary
links to preserved strict reports and separate policy assessments; no case was
repeated to obtain zero skips. The physical block used the approved built-in
outputs with every route muted and is now complete.

The final coverage review identified that the short interleaved oracle matrix
exercises unavailable LRA but does not directly compare a stable nonzero LRA.
Supplemental tests add 64-second dynamic programmes with excluded gate windows,
and silent programmes, at all four rates for stereo/5.1 Source/Current. They compare
every M/S/LRA window, gates, percentiles, metric status and per-channel peak location
against independently written interleaved audio. Only test tooling changes; no
application binary is rebuilt. Candidate-04 and its both-platform results are
recorded separately, preserving the earlier evidence and identities.

Candidate-04's actual Sequoia supplemental gate passes: 370 assertions across
32 full-resolution dynamic/silence cases, then the original 5,261 assertions,
56 captures, ten native reference-picture cases, three M8 migrations and nine
runner guard tests. Every original capture hash, metric result, availability and
peak location matches candidate-02 exactly. The dynamic fixtures have stable LRA
about 10.46 LU; silent metrics remain unavailable. All frozen artifacts were verified
before execution.

Candidate-04's same supplemental gate also passes on actual Intel High Sierra:
370 extended assertions/32 cases, 5,261 baseline assertions/56 captures, ten native
video cases, three M8 migrations and nine runner guard tests. Its baseline reports
and capture hashes exactly match candidate-02. The extended cross-platform comparison
passes across 31,168 full-resolution windows and 128 per-channel peak pairs.
Programme loudness/LRA/peak metrics, window energies, gate counts and peak locations
compare exactly; the maximum per-channel true-peak difference is
`3.857309866213148e-15 dB`, within the existing 0.01 dB comparison budget. This test
budget does not widen the profile's strict true-peak acceptance ceiling.

**High Sierra's interactive UI gate is complete.** Twenty-one recorded observations
cover stereo and mixed-format six-file staging, duplicate-role rejection, mapping
invalidation of alignment confirmation, atomic adoption, measured Source/Current
and explicit M3 application/verification, device-closed preparation, all-channel
filename/role diagnostics, save/reopen, missing-member recovery, identical-content
relink and mismatched-content rejection. Both groups restored eligible caches and
their exact gain without remeasurement. The six-file measured trim was
`-2.596466531665 dB`, verified at `-23.000000 LUFS / -23.361284 dBTP`.

Opening an M8 document recovered its settings and rejected incompatible caches.
The explicit format-upgrade dialog offered Save As; the saved schema-2 copy
preserved document identity, processing, profile, sections, audition selection and
video settings, while the original schema-1 SHA-256 remained unchanged. Cancelling
active three-hour group opening and import retained the previous session. Close
clicks during active opening left the app running; the transient close-notice text
was not captured. Clear Source removed the entire group and its measurement
authority. The final unsaved-close guard worked and the test app closed normally.
No device was opened for this UI block.

Native multi-file selection and one stereo save/reopen step used coordinated user
assistance because Screen Sharing dropped input. Later actions were observed
directly through Screen Sharing. The executed candidate-02 app is byte-identical
to candidate-04, verified again after the checks. See
`high-sierra-ui-candidate-04-observations.json` and
`high-sierra-ui-final-artifacts.json`; the earlier pending record and unscored
navigation attempt remain historical. **All required M9 gates now pass within the
approved scope.** Physical multichannel listening is a separate evidence level.

Occasional picture skips retain the accepted development policy. Record counts
and timestamps; preserve the strict harness result separately. No repetition just
to obtain zero skips. Audio interruptions, accumulating drift, sustained freezes,
identity errors and lifecycle/memory failures remain blockers. Optional paint tracing
remains confined to its separate diagnostic executable and is not a playback fix.

Full physical multichannel listening is separate and is not established by these
muted two-output checks. Full physical 7.1 listening remains pending independently.

## Reproduction and evidence

The preserved package contains the exact source manifest, universal executables,
test runners and pinned generated/reference fixtures. Run its automated runner
from the package root, writing to a **new sibling directory**, never into the
frozen package:

```sh
python3 -B Tests/run_m9_qualification.py --platform sequoia --output /absolute/new-results
# On actual High Sierra, use its available Python 3.8 and --platform high-sierra.
```

The runner checks OS, executable architecture, source/build identity and artifact
hashes. It runs hardware-independent tests only. Bidirectional exchange, interactive
UI and the explicitly coordinated physical block have separate records and gates.
`compare_m9_reports.py` compares the two result folders; the M1–M8 comparators retain
their existing metric/sample budgets. The optional paint diagnostic remains isolated
from ordinary standalone and qualification executables.

The build uses JUCE 7.0.12 (`4f43011b96eb0636104cb3e433894cda98243626`), Apple clang
17, Release `-O3 -DNDEBUG`, universal `arm64;x86_64`, and deployment minimum 10.13.
`build/logic-m9` has the named measurement/scaffold/compliance/proposal/diagnostics/
sections/audition/multichannel/session/programme/video test options enabled. Build
only explicit standalone/test targets; never use the default/all or VST3 target:

```sh
cmake --build build/logic-m9 --target LevelScopeStandalone LevelScopeProgrammeTests \
  LevelScopeMeasurementTests LevelScopeScaffoldTests LevelScopeComplianceTests \
  LevelScopeProposalTests LevelScopeDiagnosticsTests LevelScopeSectionsTests \
  LevelScopeAuditionTests LevelScopeMultichannelTests LevelScopeSessionTests \
  LevelScopeSessionPlaybackTests LevelScopeVideoTests LevelScopeVideoQualification \
  LevelScopeVideoPaintDiagnostic LevelScopeVideoPaintTraceTests -j 8
```

A rebuild is a new artifact and must not replace the frozen tested package.
Current evidence is collected under
[`Tests/Reference/results/m9-qualification-2026-10-25`](../../Tests/Reference/results/m9-qualification-2026-10-25).
The final integrity audit matches all 187 current source-manifest entries and all
442 candidate-04 artifacts. All twelve pre-existing checkpoint tag objects and
targets are unchanged, as are the registered working-copy locations and the video
working copy. Raw media, sample captures and build products remain in the ignored
local test folders. The latest remote main was still M8
`8e984889c5e73670f370c6c568c8ea1d7902a35c`; main received the accepted M9 checkpoint
by fast-forward, without integration changes to the tested implementation.

The preserved tested app is
`/Users/fuzz/Developer/LevelScope/build/m9-qualification-2026-10-25/candidate-04/LevelScope Standalone.app`.
It must not be replaced by a later rebuild. The publication annotation links the
checkpoint commit to this source manifest, original binary hashes and evidence.
Publication documentation is outside the tested source manifest. No additional
device runs or installation changes were needed for this unchanged implementation.
The original qualification records retain their pre-publication status fields;
the publication record and annotated tag establish the subsequent saved state.
Any later implementation change must be assessed before reusing qualification.

## Deferred work

No stems, dialogue analysis/references, arbitrary file-count import, offsets,
RF64/BW64 extension, automatic correction, section processing, dynamics, export,
reporting or VST3 parity. 7.1.4/immersive expansion remains deferred. The next
bounded proposal after M9 is verified PCM export and artifact-bound reporting,
with large-container requirements resolved for actual broadcast delivery needs.

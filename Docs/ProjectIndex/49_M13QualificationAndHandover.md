# M13 bounded section-led proposals — implementation and qualification record

Status, 28 October 2026: **candidate-02 qualified on actual arm64 Sequoia and Intel High Sierra; published combined M12/M13 checkpoint**.
The frozen package and its source manifest are preserved under
`build/m13-qualification-2026-10-28/candidate-02`; candidate-01 remains preserved.
The user authorized commits, an annotated tag, main integration and pushes after
review. The checkpoint is
[`m13-bounded-section-proposals-2026-10-28`](https://github.com/carlozaiana/LevelScope/tree/m13-bounded-section-proposals-2026-10-28).
Its peeled commit is both the saved logic checkpoint and main revision at publication;
the annotation records that full commit and the exact tested identities below.
The [publication record](../../Tests/Reference/results/m13-publication-2026-10-28/README.md)
records preservation, preflight and the unchanged-source fast-forward from M11
`cfb1033aa7fb65a2655db276bee696161196bdeb`. Published M12 candidate-03 and its
`m12-section-trim-sequoia-2026-10-27` tag retain their original objects and historical status.

Read [47](47_M13BoundedSectionTrimProposals.md),
[48](48_M13ImpactQualificationPlan.md), and the [shared specification](../ProposalEngineSpec_v2.md).

Frozen candidate-02 source manifest: `50eba19a2d66096607ec878c8eaa05e610580c7f5e8854d7b375d846a89c5f27`.
Measurement/build identity: `5fd55bb1a81c80b4b556715936b6873bccae29a5d04c00b0890a9b43296a20d2`.
Universal standalone executable: `6f843857838c64accc92536f6ccba2c9b0f36113dacca56dc2beef55c96233c5`.
Base documentation HEAD at freeze: `0fc4c5a7f2364f57760228f5c5ecf2f6e0714197`;
the manifest identifies the then-uncommitted implementation. The checkpoint maps
that exact source to a saved commit without replacing the tested package. All named
test/standalone targets were rebuilt with the same production identity, JUCE 7.0.12,
arm64+x86_64 and deployment target 10.13. No plugin or installed app was built/replaced.
Freeze-time `qualification_status` fields in the preserved identity files are
historical and must not be edited after testing. Current qualification is recorded
separately by the phase receipts and final completion gate; a later checkpoint
commit must bind to this exact manifest and these original binaries.

## Implemented development path

The shared measurement kernel has an optional continuous peak observer. It sees the
existing sample/interpolation calculations without resetting filter history or changing
ordinary loudness/peak arithmetic. `SectionPeakCollector` retains all-channel maxima,
source-frame/eighth-frame locations and padding attribution. Sample attribution uses
half-open section intervals. True-peak neighborhoods conservatively extend 50 ms plus
64 frames on each side; only the bucket lookup clamps outside-programme coordinates.
The original location/padding flag remains visible, and programme duration never grows.
A bucket/range-maximum representation bounds storage even for 10,000 short sections.

`section-trim-proposal-v1` preserves no-change and M3 global-only outcomes first.
A local correction requires current matching Source/Current inputs, compatible explicit
manual geometry and no outstanding manual draft. It selects the whole affected set,
requests deficit plus 0.5 dB, caps additional attenuation at 3 dB, and enforces eight
changed sections, 20% local envelope/transition support and M12 cumulative gain limits.
A representable stored plateau cannot exceed the cap through rounding. Capped reserve
and uncovered deficit are disclosed. Over-budget sets are rejected without selecting
an arbitrary subset. Large failure explanations have a labelled bounded text view;
full-resolution correction data is retained. Neighbour differences and transition
slopes are descriptive, not new compliance rules.

The private service reuses one verified programme copy and at most three full passes:
optional baseline attribution, exact local plan, and exact final section/global plan.
Compatible profile-independent attribution can be reused. An identical final plan
skips the redundant pass. Candidate snapshots have a separate `privateCandidate`
domain; they cannot enter the Current cache, audition selection or artifact approval.
The measured evaluator version records this domain extension; the EBU profile and its
digest are unchanged. A candidate must strictly pass the selected programme rules
before its new section envelope is offered for Apply. Global-only drafts retain M3's
prediction/apply/verify semantics.

Explicit stopped-only Apply changes section and global settings atomically, invalidates
Current truth and retains Source. Independent Current remeasurement remains required.
Cancellation/failure after Apply retains the applied plan as UNVERIFIED; completed
failing measurement displays the actual failure. No rollback, retry or new correction
is automatic. Already-pass performs no copy/measurement, changes no audio revision or
valid measurement, and does not create an unsaved actionable draft or block saving.

The existing section timeline displays the proposed cyan envelope beside applied green
and manual-draft orange curves. Selecting an annotation highlights its source interval
and shows the associated proposed correction. Geometry changes make a proposal stale;
cosmetic names and video/monitor edits preserve audio identity. Transient actionable
proposals must be discarded before saving; only explicitly applied M12 schema-4 plans
persist. Reopening keeps the established stopped/device-closed/muted behavior and
conservative cache validation. Current and actual exported-audio gates, and one final
PCM24 dither stage, remain unchanged.

## Development evidence, not frozen qualification

Records remain under `build/m13-qualification-2026-10-28`:

- First 273-check end-to-end run: `dev-short-03.json`, build identity
  `aed13b8e09654bab09d8f58ca02d8cbfc6e9d09629c50c6eaaa189fbae23fa70`.
- Forty-case rate/representation matrix: `dev-matrix-01.json`, 2,719 checks, build
  `c9b23ba34a7a48868c603666020929a4d4f254512ad2c186a24e00f88dea4ff0`.
  Four rates, stereo/5.1, interleaved/grouped sources and all three composition mappings.
- Extended injected-audition/export smoke: `dev-smoke-04.json`, 201 checks, same
  development build. Independent libsndfile/section/global/TPDF check matched every
  one of 2,880,000 5.1 PCM24 codes in `dev-smoke-04-independent-pcm.json`.
- The sandboxed first workflow attempt could not access JUCE's temporary-copy
  directory. Normal-cache audio-only execution passed; both attempts remain recorded.
- Subsequent safeguards and test additions require fresh frozen-candidate execution.
  These development reports are not substituted for that evidence.
- Extended full matrix `dev-matrix-02.json`: 3,464 checks / 40 cases passed, build
  `5fd55bb1a81c80b4b556715936b6873bccae29a5d04c00b0890a9b43296a20d2`.
  Final test-only additions pin synthetic PCM across platforms, retain bounded
  pre-monitor captures and exercise two-pass/no-solution assessment branches.
- Silent development UI exposed an obsolete job-status line after Apply. The model
  now replaces it with explicit applied/UNVERIFIED wording; a regression checks it.
  No processing arithmetic changed. Frozen UI qualification must verify this fix.

## Completed frozen-candidate qualification

[Archived machine-readable evidence](../../Tests/Reference/results/m13-qualification-2026-10-28/README.md)
retains phase identities, command receipts, hashes and UI observations. Media,
executables and connection details stay outside this archive.

| Required phase | Sequoia arm64 | Actual High Sierra Intel |
| --- | --- | --- |
| Explicit short/regression/reference set | PASS, 22 required checks (candidate-01 binary-equivalent evidence) | PASS, 22 required checks on candidate-02 |
| Silent combined M12/M13 UI | PASS, 17 cases | PASS, 17 cases |
| Injected H.264/ProRes timing and section audition | PASS, 80 receipts / 382 assertions plus timing suite | PASS, 80 receipts / 382 assertions plus timing suite |
| Single consolidated three-hour workflow | PASS, 117.3 min including independent RF64 checks | PASS, 202.1 min including independent RF64 checks |
| Cross-platform session exchange and numerical comparison | Six actual High Sierra sessions received/remeasured PASS; final comparison PASS | Six actual Sequoia sessions received/remeasured PASS; final comparison PASS |
| Four coordinated muted physical video cases | PASS on candidate-02, 0 skips, configuration unchanged | PASS on candidate-02, 0 skips, configuration unchanged |
| Single cinema evaluation | Completed safely; proposal rejected by policy | Sequoia only by design |


Candidate-02 changes only five qualification-tool Python files. The first physical
runner stopped before opening a device because a composition report uses
`measured_programme_definition_digest`, not the single-file content-hash field.
The corrected runner validates the composition definition; nine runner regressions
pass. The failed preflight is preserved as a failure under `historical`, not rewritten.
All production sources, C++ tests, dependencies and executable bytes are identical
to candidate-01 (manifest `2398d0106975e57f8d1f7484b4141ff5becb418649cc65b5f42e8d981c6f6a30`).
`candidate-02-equivalence.json` enumerates the five tool-only differences. Completion
and comparison tools require exact equivalence and reject any changed product/test
binary or non-tool source; fresh tool tests are required on both Macs. Completed
candidate-01 Sequoia phases, UI and completed long workflow retain their original
identities and are reused explicitly. No redundant cinema or long run is required.

The Sequoia short phase passed all 22 explicitly required checks on actual arm64
macOS 15.7.9 (24G830), Apple M1 Max. M13 reports 3,671 assertions across 40 cases. Fresh checks include the
M12 processing suite, selected M1–M11 regressions, pinned EBU/ITU vectors, runner
validation, older-reader rejection, and independent final PCM comparison. The phase
records command receipts, hashes and environment in
`candidate-01-results/sequoia/short`; its passing summary does not establish complete
M13 qualification.

Silent frozen-candidate UI checks have verified private candidate separation,
explicit Apply/UNVERIFIED/independent measurement, already-pass preservation,
global-first correction, stale rejection, exact manual-gain persistence, annotation
independence, recoverable Current-cache damage and strict artifact publication gates.
All 17 required Sequoia silent UI items now pass, including cancellation on a
200-second disposable source and honest 100%-duration rejection. Exact observations
are recorded in `sequoia/ui-checklist.json`. The device remained closed throughout.
The frozen injected H.264 CFR/VFR and ProRes exact/fractional-rate phase passed
382 checks with 80 Source/Current/offset/position receipts, plus the timing suite.
These are not physical presentation tests. The single cinema evaluation completed with an explained policy-limited failure
and no Apply, export or playback; details follow below.

High Sierra candidate-02 completed all 22 short checks in 1,762.03 seconds on
10.13.6 (17G66), x86_64 i7-3770K, using existing Python 3.8.10 and libsndfile 1.2.2.
Its M13 suite passed 3,671 assertions/40 cases and M12 passed 3,304 assertions/32
cases. All 17 silent UI cases passed through the actual High Sierra app, including
the combined M12 manual-processing checks. Screen Sharing image updates lagged some
input; each outcome was observed after updating, with no product changes needed.
The saved manual gain retained exact double bits `bfbf9add3746f62e`; a deliberately
damaged Current cache was discarded with settings and valid Source preserved.
The UI export required explicit publication after separate Current/artifact passes;
all three PCM/JSON/HTML hashes matched the manifest. The app was closed afterward.

Both six-case foreign-session exchanges passed: caches were conservatively discarded,
exact settings retained, and the receiving Mac independently remeasured the result.
High Sierra injected video passed in 48.08 seconds. Both coordinated physical blocks
completed four eight-second H.264 VFR/ProRes 422 Source/Current cases, every route
muted, with zero recorded skips and unchanged audio configuration. These are bounded
physical presentation tests, not full multichannel listening qualification.

Sequoia's single consolidated three-hour workflow passed 1,078 assertions plus the
independent RF64 checks in 117.3 minutes. The 1,000-section 5.1 composition used three
private passes, explicit Apply and independent Current verification. Source/Current
captures at eight positions matched the all-channel oracle without underruns.
The actual 9,331,200,104-byte RF64 output passed reopened measurement and independent
libsndfile checks, with full artifact hash
`953b2bcecbf068f5fe153461fb1ea0999939ddb40adf607052b05f4ad2146ade`.
Additional peak working RSS was 83.0 MiB; UI heartbeat p95 was 15.03 ms and cancellation
33.64 ms. Save/parse took 98.53/475.26 ms, the session was 10,494,401 bytes and combined
JSON/HTML 388,938 bytes. All planned Sequoia budgets passed. Independent decoding
uses bounded sample locations and complete-file hashing, not an assertion that
every output sample was independently decoded.

High Sierra's same consolidated workload passed 1,078 assertions and the independent
RF64 checks in 202.1 minutes. Its three private passes produced the same one-section
correction and final global gain; independent Current verification then passed.
Both actual-platform RF64 artifacts have the identical full-file SHA-256 above.
Each independent decoder run checked 98,304 PCM codes at bounded decode locations;
the separate section/global/single-dither oracle checked 182,880 codes. No historical
long workflow was recursively launched or repeated.

| Long-workflow resource / outcome | Sequoia | High Sierra | Required bound |
| --- | ---: | ---: | --- |
| Workflow plus independent decoder | 117.3 min | 202.1 min | C++ workflow <=240 / <=300 min |
| Additional peak working RSS | 83.0 MiB | 49.6 MiB | <=256 MiB |
| UI heartbeat p95 | 15.03 ms | 15.04 ms | <=100 ms |
| Cancellation acknowledgement | 33.64 ms | 65.54 ms | <=1,000 ms |
| Save / parse | 98.53 / 475.26 ms | 263.47 / 982.04 ms | <=5,000 ms each |
| Session bytes | 10,494,401 | 10,493,973 | <=128 MiB |
| JSON + HTML bytes | 388,938 | 389,010 | <=16 MiB |
| Observed disk growth | 26,651,189,248 bytes | 31,721,431,040 bytes | 65 GiB free preflight passed |

The final [cross-platform comparison](../../Tests/Reference/results/m13-qualification-2026-10-28/cross-platform-comparison.json)
passed all 40 cases, with **zero float32 ULP difference** in the retained processed
captures and maximum measured-metric difference `3.552713678800501e-15`. Input,
ordered-role, plan, verdict and LRA-stability checks passed. These comparison
tolerances did not widen the profile's evaluation bounds.

The [required-evidence completion gate](../../Tests/Reference/results/m13-qualification-2026-10-28/qualification-complete.json)
passed after both long phases. **No required M13 qualification check remains
outstanding.** Actual High Sierra M12 processing/UI and the selected deferred video
checks are now qualified as part of this combined M13 candidate. This does not
rewrite the preserved M12 candidate-03's historical Sequoia-only status.

Post-test integrity checks verified all 250 working source files against candidate-02,
both preserved packages on both Macs, the original standalone executable bytes,
all 16 pre-existing tag objects/targets, branch heads and registered worktree heads.
At qualification completion, main remained M11 and no publication had been performed.
The user subsequently approved publication and integration. Fresh publication
preflight revalidated all source, candidate/evidence hashes and the complete gate;
remote main still identified M11, an ancestor of the logic branch. The combined
checkpoint therefore integrates by fast-forward with no implementation change,
rebuild or repeated device/long tests. Original qualification receipts retain their
pre-publication authorization fields; the separate publication record documents
the later approval. No next milestone, installation or audio-configuration change
is included.

## Cinema provenance and container boundary

User-reported Nuendo values dated 28 October 2026 are recorded in document 47;
exact version/settings/range equivalence remain unconfirmed. The historical screenshot
is **Cubase 10**, kept separately as supplementary evidence. Dialogue-gated observations
remain future context and do not become M13 objectives or reference-vector tolerances.

Original: `/Volumes/Tucan/GuerraDiCesare/LGdC_Masters/LGdC__51_Master_V2.wav`.
SHA-256: `1fb7f93a33d600bf0a62d63e642e5703f03aec8e6920f8a6e510cf65f61fdeaf`.
The original has an extended 92-byte ds64 chunk, zero ds64 sample count and trailing
non-payload metadata outside the limited product-qualified RF64 subset. Do not broaden
that importer as part of M13 or call this an audio-engine failure.

A qualification-only rewrap preserves the declared PCM payload exactly and leaves the
original unchanged. Payload SHA-256:
`1a4860f5a3b8d59bcc7401e5f16d57686674475658b36a90f2f6baca4c19300e`.
Canonical RF64 copy SHA-256:
`2b00ba5ae88ef06b5d949adf7745144abf8b38d4cb97f78e8abe94d5916b9af6`.
The declared payload is 269,875,200 frames at 48 kHz/24-bit, L R C LFE Ls Rs,
range `[0,269875200)` (5,622.4 seconds). No sample conversion, gain or dither occurred.
`cinema-input/provenance.json` records the copy, original/container and payload hashes.
This is not a product RF64-support extension, BWF-compliance claim or proof of exact
Nuendo range/settings equivalence. Media and screenshots remain outside GitHub.

The user selected **one whole-programme section** for the single bounded evaluation.
The single frozen-candidate evaluation completed in 1,412.12 seconds (23.54 minutes).
Source and unchanged 0 dB Current produced identical measured values:

| Metric | LevelScope candidate-01 | User-reported Nuendo observation |
| --- | ---: | ---: |
| Integrated loudness | −25.978754 LUFS | −26.0 LUFS |
| LRA | 20.005237 LU | 19.9 LU |
| Maximum true peak | −0.826621 dBTP | −0.85 dBTP |
| Maximum momentary | −5.198827 LUFS | −5.4 LUFS |
| Maximum short-term | −9.781314 LUFS | −9.8 LUFS |

LevelScope per-channel TP, L/R/C/LFE/Ls/Rs: −0.826621, −1.473165,
−1.037365, −1.980466, −1.961937, −1.994212 dBTP. The report preserves
full-precision values, statuses, identities and locations. These differences do not
establish an engine error: external version/settings/range equivalence is still
unconfirmed. Dialogue-gated measures were not implemented or evaluated.

The proposal outcome is **No solution found within this policy**. Global-only
loudness/peak constraints do not intersect. The sole local section covers 100% of
the programme, exceeding 20%. At the intended normalization target its requested
reduction including reserve is 3.652133407380 dB; the 3 dB cap leaves no reserve and
0.152133407380 dB uncovered deficit. The whole envelope is rejected. One private
baseline-attribution pass was used; no processed candidate was measured, applied,
exported or played. The qualification check passed because it completed this
honest bounded assessment safely; the programme itself did not pass the selected
rules and no corrective plan was accepted. Do not infer scenes or retry to force
a pass. Listening remains separately coordinated.

## Remaining limitations and next proposals

### Reproduction and evidence integrity

Use the preserved candidate package for exact-binary reproduction; rebuilding does
not replace its tested identity. `identity.json` records each source/artifact hash,
the compiler and JUCE commit. The local build directory is `build/logic-m13`:
Release, universal `arm64;x86_64`, deployment target 10.13, with the named
`LevelScopeStandalone`, `LevelScopeVideoQualification` and test targets enumerated
by `Tests/package_m13_qualification.py`. Never build the default all/plugin target
or install a rebuild as part of reproducing these checks.

Each phase is explicit and creates a new output folder. From the frozen package,
use `python3 -B Tests/run_m13_qualification.py --platform sequoia --phase short
--output <new-folder> --decoder-library <recorded-libsndfile>`; use
`--platform high-sierra` and Golem's existing `/usr/local/bin/python3` for actual
Intel execution. Other modes are `long`, `video`, `exchange` and `cinema`; their
input paths and exact arguments are retained in each phase's `runs.json`.
Do not invoke historical milestone runners recursively or rerun an accepted long
phase merely to regenerate a report. UI and coordinated physical receipts remain
separate from these command-line phases.

From the repository, the final integrity/comparison gates are:

```sh
python3 -B build/m13-qualification-2026-10-28/candidate-02/Tests/compare_m13_reports.py \
  --candidate build/m13-qualification-2026-10-28/candidate-02 \
  --equivalent-candidate build/m13-qualification-2026-10-28/candidate-01 \
  --sequoia build/m13-qualification-2026-10-28/candidate-02-results/sequoia \
  --high-sierra build/m13-qualification-2026-10-28/candidate-02-results/high-sierra \
  --output build/m13-qualification-2026-10-28/candidate-02-results/cross-platform-comparison.json
python3 -B build/m13-qualification-2026-10-28/candidate-02/Tests/complete_m13_qualification.py \
  --candidate build/m13-qualification-2026-10-28/candidate-02 \
  --equivalent-candidate build/m13-qualification-2026-10-28/candidate-01 \
  --results build/m13-qualification-2026-10-28/candidate-02-results \
  --output build/m13-qualification-2026-10-28/candidate-02-results/qualification-complete.json
```

The comparison tolerances are cross-platform checks only: 0.01 LU/dB for
I/M/S/TP/SP and gain, 0.1 LU for LRA, and at most one float32 ULP for processed
captures. They do not relax the profile's evaluation bounds. Completion requires
the exact required run sets, report/log hashes, actual OS/architecture, 17 UI cases,
four physical cases per Mac, fresh repaired-tool tests, session exchange and the
bounded cinema result. A missing or changed required receipt fails the gate.

### Product limitations

The report app caption retains the inherited `LevelScope Standalone M12` text;
its exact `build` field identifies this M13 candidate. Processing remains
`section-gain-v1`, and this caption is not a claim that an M12 binary ran the test.

M13 has one fixed section-led algorithm; no configurable controls, segmentation,
marker/EDL import, coordinated neighbour adjustment, limiter, compressor, multi-strategy
comparison or component-relative dialogue processing. Three additional private passes
are a work limit, not a proof of optimality. The 20% local-support budget excludes a
uniform global gain, which changes the full programme. LRA stays informative.

Later separately reviewed work follows the practical stereo/5.1 roadmap: complete mono
section/export support and basic markers with exact programme-time mapping; qualified
standalone limiting and downward/upward processing; dynamics-led proposals; enhanced
section-led and then mixed strategies. Dialogue/M&E objectives require explicit methods,
selected objectives and acceptance criteria, plus measured verification of the complete
programme. Basic programme-rule pass must not suppress a future explicitly selected,
validated objective. 7.1.4/immersive expansion and VST3 parity remain separate/deferred.

Occasional picture skips remain an accepted development limitation only when audio
continues, authoritative-clock synchronization does not accumulate drift, and identity,
memory and lifecycle requirements pass. Preserve skipped-count/timing records and
historical strict failures. Full physical 7.1 listening remains separately pending.

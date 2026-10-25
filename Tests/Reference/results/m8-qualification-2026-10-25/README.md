# M8 candidate-01 qualification evidence

Status: **candidate-01 qualified on both actual platforms; ready for checkpoint
publication.** No M8 commit, tag or merge has been made.

The [combined qualification summary](qualification-summary.json) records every
required gate. See the [implementation and handover](../../../../Docs/ProjectIndex/36_M8SessionsQualificationAndHandover.md)
for scope, exact identities, reproduction instructions and resource results.

The exact candidate is identified by `candidate-identity.json`. Its ancestry is
M7 commit `cfbf3f981dece326d2363b911e5914ea61f1796c`; M8 is currently uncommitted.
The source manifest, rather than the ancestry commit, binds these results to the
working implementation. Package identity SHA-256:
`dd625b0e2241afde0e6cbd5bbe36f733d8e23bdd9afc3bba92cd62fc556e8b9f`.

Sequoia ran natively as arm64 on macOS 15.7.9 (24G830), Apple M1 Max. Actual High
Sierra ran natively as x86_64 on macOS 10.13.6 (17G66), Intel i7-3770K. Both
environment reports record CPU, executable architecture, build, source and artifact hashes. No
physical audio device was opened. Captured audition uses injected test devices.

Fresh results on each Mac:

- M8 persistence: 2,073 assertions, including real three-hour source verification.
- M8 restored playback: 2,431 assertions, 48 full-channel sample-capture cases,
  five H.264/ProRes restored reference-picture cases, no embedded audio outputs.
- M1–M7 numerical, model, diagnostics, sections, proposal, compliance, audition,
  channel-layout and lifecycle regressions; official ITU (39) and EBU (70) vectors.
- Runner guards and optional paint-diagnostic isolation/exclusive-output tests.
- Silent session UI: see [Sequoia](sequoia-ui.json) and [High Sierra](high-sierra-ui.json) for exact observed cases and runtime
  executable paths. The development and packaged UI executables have identical
  SHA-256. Test fixture/audio files remain outside Git.

- Bidirectional session exchange/relink passed; exact settings/sections retained
  while foreign caches were discarded and Current remained UNVERIFIED. Native
  reports: [Sequoia → High Sierra](exchange-from-sequoia.json) and
  [High Sierra → Sequoia](exchange-from-high-sierra.json). Each exchange run passed
  2,042 assertions; the combined summary records its inputs and invocation.
- Cross-platform [M8 restored samples](comparison-m8.json), [M1–M5](comparison-m1-m5.json),
  [M6 captures](comparison-m6.json) and [M7 captures/native video](comparison-m7.json)
  passed. All compared captured streams were bit-identical (maximum float32 ULP 0).
- Final [local preservation audit](preservation-local.json) and
  [High Sierra package audit](preservation-high-sierra.json) passed: 379 package
  artifacts unchanged on each Mac, 174 working-source entries match the manifest,
  all eleven earlier tag objects/targets and unrelated work preserved.

Both raw `summary.json` files deliberately leave `qualification_complete` false
because that runner alone does not establish cross-Mac exchange or interactive completion.
The combined summary adds those gates without altering the original reports.
The nested M7 runner's generic physical/UI gate does not imply that M8 repeated
physical endurance: `m7-reuse-review.json` bounds reuse of unchanged engine paths.
The accepted M7 timing/long-measurement records retain their historical identities.
Occasional picture skips remain an accepted development limitation; sustained
freezes, drift, interrupted audio, identity and lifecycle failures are not accepted.
Full physical 7.1 listening remains separately pending.

Retained locally: the immutable universal package, downloaded/generated reference
media, raw captured `.f32` streams, UI session fixtures and full execution folders.
Excluded here: binaries, raw audio/video, credentials, connection setup and Finder
metadata. These text reports and manifests are prepared for a later authorized
checkpoint publication, not an M8 tag or merge.

`evidence-sha256.json` inventories the retained text evidence (excluding itself).

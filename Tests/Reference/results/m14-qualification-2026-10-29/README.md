# M14 Sequoia phase evidence

Candidate-02 review identity:
`9dd22a598281e83431039f210202a8a24fe49517294cc6bb8d83032c4975e013`.
**Sequoia phase complete; High Sierra, video and full-duration mono/RF64 qualification pending.**
Ready for development-checkpoint review; no publication or integration.
Main remains qualified M13 `a9be070ab9eb5cb4f591ff9d991d8999c3a016d0`.

See [the handover](../../../../Docs/ProjectIndex/52_M14QualificationAndHandover.md),
[approved contract](../../../../Docs/ProjectIndex/50_M14MonoProcessingAndVerifiedExport.md)
and [pre-implementation impact selection](../../../../Docs/ProjectIndex/51_M14ImpactQualificationPlan.md).

Candidate-01 supplied the tested production/C++ binaries. Candidate-02 changes
only three qualification Python files; [equivalence](candidate-equivalence.json)
records exact source and binary comparisons. Both frozen identities are retained.
The original generic grouped-reader failure is inconclusive and remains historical;
the required corrected `sequoia/legacy` phase checks explicit parser messages from
the unchanged M13 binary, including a readable mono baseline control. Both helpers
intentionally exit 1; the old applied-mono parser failure is distinguished from the
helper's schema-1-only assertion. No audio processing was dropped or rewritten.

Actual test platform: Apple M1 Max, arm64, Sequoia 15.7.9 (24G830). All 16 named
app/test executables contain arm64 and x86_64; Intel minimum OS is 10.13. This is
cross-build evidence, not actual High Sierra runtime qualification. No physical
device was opened and no video was loaded or decoded. Every automated phase records
arguments, environment, exit status, logs and report hashes. The completion tool
requires short, cinema, sustained, corrected legacy and all 11 UI cases, while
explicitly refusing to claim full platform/video/duration qualification.

Fresh bounded cinema evidence validates original SHA-256, exact format/duration,
ordered roles and six decoded spans against the preserved M13 copy, including
offsets beyond 4 GiB. Archived M13 payload identity and prior cinema measurements
may inform comparison; they are not new Current measurements of the imported file.
The complete cinema proposal/export workflow was not repeated. Small RF64 fixtures
and six sparse boundary files exercise parsing, metadata exclusion, malformed input,
padding and bounded seeks. These are not full-duration mono/RF64 evidence.

The single sustained 90-minute mono/192 kHz workflow passed 1,065 assertions and
independent bounded PCM/hash/report verification, with ordinary RIFF output. It took
35.16 minutes, used 117.5 MiB additional peak RSS, recorded 15.03 ms heartbeat p95
and 11.30 ms cancellation response. Session save/parse took 53.87/234.69 ms.
See [completion](sequoia-phase-complete.json) and the handover for exact budgets
and the distinction between native full measurement and bounded independent reads.
New three-hour mono/RF64 end-to-end testing, actual High Sierra/UI, cross-Mac exchange
and all video execution remain pending for the combined window after M18, no later
than M19. Future milestones require approval. Occasional picture skips remain
accepted only under the existing audio-continuity/drift/identity/lifecycle limits.
Full physical 7.1 listening remains separately pending.

The development sandbox failure under `historical-development` was inability to
write an owned JUCE temporary copy; the authorized rerun passed. It is not hidden
or described as a product pass. Other development records are supplemental, not
substitutes for frozen-candidate evidence.

This is a text-only Git archive: no audio/video, binaries, raw captures, connection
details or credentials. Preserved local packages and results under
`build/m14-qualification-2026-10-29` remain authoritative for reproduction. Reports
and manifests retain their original identities; the archive does not rewrite them
to suggest execution of a different binary. Nothing is committed, tagged or merged.

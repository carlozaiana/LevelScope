# M7 checkpoint publication — 25 October 2026

The annotated tag [`m7-pcm-7.1-2026-10-25`](https://github.com/carlozaiana/LevelScope/tree/m7-pcm-7.1-2026-10-25)
identifies the accepted M7 checkpoint. Its peeled commit is also the main revision
at this publication. The tag annotation records that exact commit and the frozen
candidate identities; it does not change the binaries' original build identity.

[Publication preflight](publication-preflight.json) records verification against
latest main `5e07bbd29930db34236d6dfcede7fea14983f90a`, already incorporated by
candidate-03. No source or integration fix was needed. All 165 source inputs match
the frozen manifest, all 365 package artifacts match, and all 317 completion
evidence files retain their recorded hashes. The 46 fresh physical reports were
reevaluated under the approved policy and matched the original accepted audit.
Publication did not rebuild the candidate or reopen audio devices.

[Historical evidence hashes](historical-evidence-sha256.json) preserve all 510
earlier M7 evidence files, including the original strict video-skip failures.
The ten earlier local and remote tag objects and targets are recorded in the
preflight for direct comparison after publication. The existing video working
copy, unrelated files, tested packages and installations are preserved.

The [qualification report](../m7-completion-2026-10-25/README.md) and its summary
retain their pre-publication status, including `published: false`; those immutable
records describe when qualification completed. This publication record and the
annotated tag establish the later saved state. See the
[handover](../../../../Docs/ProjectIndex/32_M7QualificationAndHandover.md) for scope,
identity mapping, evidence reuse and remaining limitations.

Offline measurement/processing, automated audition/routing and the required
UI/video scope passed on actual arm64 Sequoia and Intel High Sierra. Full physical
eight-channel listening remains separately pending. Occasional picture skips are
an accepted development limitation subject to the retained audio, timing,
identity, lifecycle and memory requirements; this checkpoint claims no paint fix.

Only source, tests, manifests, documentation and text evidence are committed.
Reference media, raw sample captures, generated build products, credentials,
connection details and Finder metadata are excluded.

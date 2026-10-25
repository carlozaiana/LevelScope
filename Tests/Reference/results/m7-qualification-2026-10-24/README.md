# M7 qualification evidence

**M7 is frozen and not qualified: the required Sequoia video paint-skip gate
failed in a controlled repeat. This is not an accepted checkpoint.**

Candidate-01 failed the final-sample audition seek in the real three-hour run on
both actual Macs. Its original reports remain preserved. Candidate-02 contains
the bounded FLAC seek correction. Its complete automated suites, official vectors,
real three-hour runs, cross-platform comparisons and silent UI checks passed on
both actual Macs. Physical video timing is failed/incomplete; later device groups
are deferred. Physical eight-channel listening remains separately pending.

`candidate-02/identity.json` records the replacement package. Its source digest is
`c3df888f567a69adb83cbb6c7d2c7a1dcdf1da9498d43d688ffe890786c0e43e`.
`candidate-02/source-differences.json` records the four changed inputs and the
basis for reusing unaffected UI observations. Candidate-02's affected silent UI
observations and completed short comparisons are recorded separately.
Its full M5/M6/M7 cross-platform reports and final platform reports are now also
collected under `candidate-02`. The first and controlled failed Sequoia Source
matrices are preserved independently. A diagnostic comparison reproduced missed
paints using preserved M6 stereo as well as M7 stereo; those new observations do
not overwrite earlier M6 qualification evidence.

The user requested a video handover and a frozen M7 candidate. See
[ProjectIndex 33](../../../../Docs/ProjectIndex/33_M7VideoPaintTimingHandover.md).
The final [candidate-02 summary](candidate-02/qualification-summary.json) retains
an overall unqualified result. High Sierra's physical Source matrix passed all
sixteen cases; later physical groups on both Macs remain deferred. Final source,
package, shutdown and local/remote tag-preservation audits passed. No M7
checkpoint or merge was created.

The files described below at the top level belong to preserved candidate-01.
`candidate-01-final-high-sierra` adds its completed failed long-run report without
replacing the earlier partial record. Raw screenshots containing unrelated desktop
content are retained locally; the repository contains observations and hashes.

See [the M7 qualification and handover](../../../../Docs/ProjectIndex/32_M7QualificationAndHandover.md)
for scope, remaining checks, reproduction instructions and limitations.

`candidate-identity.json` identifies the exact working source manifest, ancestry,
compiler, build configuration, pinned dependencies and package artifact hashes.
Its ancestry commit does not contain the uncommitted M7 implementation. The
source-manifest digest is
`dd9f3392c3744c246951e33ca426545eb9d1aa8b97cafdf1ec0a2d21e73cad45`.

The `sequoia` and `high-sierra` folders contain actual-platform environment records,
native test reports, official-vector reports and runner logs. Sample-capture rows
record the hashes of raw samples retained in ignored qualification folders. Raw
audio, video, captures, downloaded vectors, build products and connection details
are excluded from this evidence directory.

Evidence levels remain separate:

1. Offline import, measurement, evaluation and processing.
2. Injected audition/routing and sample capture; silent UI and muted device/video
   timing have their own records.
3. Physical eight-channel routing and listening: **pending**, not established by
   stereo, subset or all-muted testing.

An incomplete required UI, long-programme or video check prevents the combined
M7 scope from being called qualified even when the short automated suites pass.

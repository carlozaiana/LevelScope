# M10 evidence index

The acceptance decision is in `acceptance-candidate-03.json`; until it says
`qualified: true`, M10 is incomplete. `status.json` is the human-readable status.
The complete narrative, limitations and tested binary identities are in
`Docs/ProjectIndex/40_M10QualificationAndHandover.md`.

## Frozen candidate and provenance

- `candidate-03-identity.json`: base Git ancestry, source manifest, build identity,
  universal executable hashes and all frozen package inputs.
- `candidate-03-transition.json`: exact changes from candidate-02, including the
  live applicability verdict fix and qualification-tool corrections.
- `local-integrity-candidate-03.json`, `high-sierra-integrity-candidate-03.json`:
  package integrity on both Macs and current workspace source equivalence.
- `repository-preservation-candidate-03.json`: unchanged checkpoint tags, branch,
  worktrees and unrelated files. `m9-package-preservation-candidate-03.json` records
  preservation of the previous tested source and binaries.
- `m9-component-equivalence-candidate-03.json`: precise boundaries for reuse of
  unchanged component evidence; it does not relabel an old run as a new run.

Candidate-01/02 identities and historical failures remain separate. In particular,
candidate-02 failed an interactive display check, and its interrupted High Sierra
long run has no completed result. Neither candidate is qualified.

## Current-candidate gates

1. `sequoia-candidate-03/` and `high-sierra-candidate-03/`: actual-platform core,
   relevant M1–M9 regressions, official vectors, export matrix, independent PCM24
   oracle/decoder and native/injected video evidence.
2. `*-candidate-03-long/`: real three-hour 5.1 RF64, all-channel reimport/audition,
   independent decoding and measured resource budgets.
3. `*-cross-platform-candidate-03.json`, `m7-m9-core-captures-candidate-03.json`
   and `m10-long-detail-comparison-candidate-03.json`: numeric/status/identity
   comparisons. Comparison tolerances never enlarge profile acceptance limits.
4. `exchange-*-candidate-03*.json`: both directions of fresh import/measurement
   and independent recipient decoding. Saved report verdicts are not authority.
5. `ui-*-candidate-03/observations.json`: silent interactive tests, identifying
   observed behavior separately from automated worker tests and reused HTML
   readability evidence.
6. `*-physical-*-candidate-03/` and `physical-video-summary-candidate-03.json`:
   the approved muted physical-device video cases. Raw strict verdicts and skip
   events remain intact. The accepted development policy is evaluated separately.

Each raw automated sub-suite intentionally leaves `qualification_complete: false`:
it does not independently cover UI, exchange or physical gates. The aggregate
acceptance receipt combines those separate records without modifying raw results.

`auxiliary-tools/` contains post-freeze collection/comparison scripts, not product
changes. Their output is separate from the frozen package's own tools. The evidence
manifest hashes the text records; it does not assert that missing gates passed.

No reference audio, generated PCM, private connection data, build products or
installed applications belong in this directory. Original packages, binaries,
sample captures and rendered audio remain in the local qualification folders.
Full physical multichannel listening is outside this muted-device evidence.

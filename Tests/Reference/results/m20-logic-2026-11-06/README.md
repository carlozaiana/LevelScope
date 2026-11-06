# M20 Logic integration evidence — 6 November 2026

**Focused Logic checks pass; combined M20 qualification is incomplete.** No native
M20 UI, actual High Sierra, physical audio or video execution was performed here.
No publication/tag/main integration is authorized.

Accepted local package: `/Users/fuzz/Developer/LevelScope/build/m20-logic-2026-11-06/logic-handover-03`.
Source manifest: `fd200ebdc2d1e1f43715bafcba4d2f8e8a166d28283176fcc7c408f990ca6249`.
Product build: `8c2397c14e757cd9f3c8580cd89c0a2bb065bc2e3de5a9eb97e30f03eb8f5254`.
Universal app executable SHA-256: `85775cfc34bfe77ca8b55d0b92addc9fa56e46cfae5cda4079c1f179a9a421c9`.
Accepted RegionPlanningTests executable SHA-256: `abe16f04bab5c1b364eec09d851e68bf45d4125eafc6753c297eb508fb80e0a9`.
The final local source commit is recorded in `checkpoint.json`; the package's base
commit predates the implementation and is not itself the tested implementation.

- `logic-check-summary.json`: 17 explicitly required groups, 2,144 checks, all pass.
  `m20_combined_qualified` remains false. `missing-receipts-guard.json` demonstrates
  that omitting required evidence fails the guard.
- `receipt-provenance.json`: original package identity and reuse rationale per group.
  `accepted/` assembles those receipts without rewriting their build/results.
- `receipts-01/`: original complete run, including the failed gate assertion.
  `receipts-02/gating.json`: corrected direct model counterexample.
  `receipts-03/core.json`: fresh added prior-envelope/short-radius checks.
- `accepted-package-equivalence.json`: production source/app byte-identical in all
  three packages. Only the region test source/executable changed after package 01.
  Other test bodies and executables were unchanged; no needless qualification rerun.
- `development/`: earlier failed geometry discovery and sandbox-temp-directory
  results. These remain failures, not retroactively relabelled passes.
- `explanations/`: actual measured candidate results, forecasts, geometry, policy
  limitations and split-off comparison. No generated audio or private full-window
  candidate cache is included.
- `host-and-binary.json`, package identity manifests and build logs: Sequoia arm64,
  universal binary architectures and Intel 10.13 deployment target. Cross-build
  does not substitute for actual High Sierra execution.
- `preservation.json` and `checkpoint-tags-preserved.txt`: unchanged main and all
  25 pre-existing tag objects/targets. Frozen packages/media/installations preserved.

The eight-minute cinema excerpt is unchanged and its SHA-256 matches the preserved
M19 provenance. The M20 attempt fails the selected 15 LU ceiling; a 6.7 LU gating
case likewise cannot Apply a result near 7 LU. Those expected measured product
failures are passing safety/qualification checks, not compliance passes.

See [ProjectIndex 72](../../../../Docs/ProjectIndex/72_M20LogicHandover.md) for scope,
API/identity/session contracts, numerical results, resources and UI-stage work.
Combined qualification after UI integration must still cover native interaction,
actual High Sierra/cross-platform, selected control variants and independent real
holdouts with coordinated listening. Earlier M19/UI evidence is dependency evidence,
not fresh M20 feature qualification.

Local implementation commit: `6ce361e2e0cefd6e7fe9e9710a2bf09354e8c4bb`. All 386 manifest files
were verified against its committed Git blobs. Raw native receipts retain their
original CRLF/trailing-space output so their hashes and historical identities remain
exact; source/document whitespace checks pass separately.

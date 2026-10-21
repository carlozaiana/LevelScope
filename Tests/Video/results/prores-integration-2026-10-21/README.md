# ProRes main integration — 21 October 2026

Fresh, inaudible integration regressions on actual Sequoia arm64 and High Sierra
x86_64 use the unchanged, hash-verified candidate 04. `comparison.json` enforces
the ordinary M1–M5 cross-platform numerical/performance gates with identical package
identities; no identity exceptions are needed. Both native runs contain 52 cases,
including actual-file decoding/seeking and independent nonzero audio isolation.
The unchanged large-file loading, UI, matrix and six endurance runs remain in
`../prores-2026-10-20/`; this directory does not replace or regrade those records.

`build-summary.json` and CTest output describe a separate fresh universal build of
nine named standalone/test targets. That rebuild is a CI/integration regression
check; it is not the frozen, two-platform endurance-qualified package. Its build
log and new binaries remain in the recorded fresh build directory. No VST3 target
was built or packaged. Every hardware-independent test keeps audio devices closed.

`before.json` records all nine tag objects/targets and original working-copy hashes.
`preservation-pre-integration.json` verifies unchanged tested source, local candidate
and unrelated worktree contents. Golem's `integration.json` verifies all 286 package
artifacts before and after its runs. Historical raw reports retain their bytes,
paths and tested identities. Media, compiled products and credentials are excluded.

See `Docs/ProjectIndex/27_ProResMainIntegration.md` for integration disposition and
remaining limits. The handover follow-up changes documentation/evidence only.

Publication audit: `preservation-published-strict.json` retains the strict inventory
mismatch for the original root `.DS_Store`. `publication-audit.json` records the
before/after Finder metadata hashes and distinguishes that unrelated change from
the unchanged source, candidate, other files, branches and nine checkpoint tags.
The current metadata was not reverted.

`github-ci.json`, `github-ci-artifacts.json` and `github-ci-tests.log` record the
successful GitHub run for integration publication `83bc6adc64fc495f5a2ada3e2974e6fe508d583f`.
The full build log is retained locally with its SHA-256; CI uploads only the
standalone artifact. The final follow-up changes documentation/evidence only.
